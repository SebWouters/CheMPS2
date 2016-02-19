#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libfock/jk.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.hpp>
#include <libmints/writer_file_prefix.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "chemps2/Irreps.h"
#include "chemps2/Problem.h"
#include "chemps2/CASSCF.h"
#include "chemps2/Initialize.h"
#include "chemps2/EdmistonRuedenberg.h"
#include "chemps2/CASPT2.h"
#include "chemps2/Lapack.h"

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints->DPD_ID(x)

namespace psi{ namespace dmrg{

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "DMRG"|| options.read_globals()) {

        /*- The DMRG wavefunction multiplicity in the form (2S+1) -*/
        options.add_int("WFN_MULTP", -1);

        /*- The DMRG wavefunction irrep uses the same conventions as PSI4. How convenient :-).
            Just to avoid confusion, it's copied here. It can also be found on
            http://sebwouters.github.io/CheMPS2/classCheMPS2_1_1Irreps.html .

            Symmetry Conventions        Irrep Number & Name
            Group Number & Name         0 	1 	2 	3 	4 	5 	6 	7
            0: c1                       A 							
            1: ci                       Ag 	Au 						
            2: c2                       A 	B 						
            3: cs                       A' 	A'' 						
            4: d2                       A 	B1 	B2 	B3 				
            5: c2v                      A1 	A2 	B1 	B2 				
            6: c2h                      Ag 	Bg 	Au 	Bu 				
            7: d2h                      Ag 	B1g 	B2g 	B3g 	Au 	B1u 	B2u 	B3u    
        -*/
        options.add_int("WFN_IRREP", -1);
        
        /*- The number of reduced renormalized basis states to be
            retained during successive DMRG instructions -*/
        options.add_array("DMRG_STATES");
        
        /*- The energy convergence to stop an instruction
            during successive DMRG instructions -*/
        options.add_array("DMRG_E_CONVERGENCE");
        
        /*- The maximum number of sweeps to stop an instruction
            during successive DMRG instructions -*/
        options.add_array("DMRG_MAXSWEEPS");
        
        /*- The noiseprefactors for successive DMRG instructions -*/
        options.add_array("DMRG_NOISEPREFACTORS");
        
        /*- Whether or not to print the correlation functions after the DMRG calculation -*/
        options.add_bool("DMRG_PRINT_CORR", false);
        
        /*- Whether or not to create intermediary MPS checkpoints -*/
        options.add_bool("DMRG_CHKPT", false);

        /*- Doubly occupied frozen orbitals for DMRG, per irrep. Same
            conventions as for other MR methods -*/
        options.add_array("FROZEN_DOCC");

        /*- Active space orbitals for DMRG, per irrep. Same conventions as for other MR methods. -*/
        options.add_array("ACTIVE");
        
        /*- Convergence threshold for the gradient norm. -*/
        options.add_double("D_CONVERGENCE", 1e-6);
        
        /*- Whether or not to store the unitary on disk (convenient for restarting). -*/
        options.add_bool("DMRG_STORE_UNIT", true);
        
        /*- Whether or not to use DIIS for DMRG. -*/
        options.add_bool("DMRG_DO_DIIS", false);
        
        /*- When the update norm is smaller than this value DIIS starts. -*/
        options.add_double("DMRG_DIIS_BRANCH", 1e-2);
        
        /*- Whether or not to store the DIIS checkpoint on disk (convenient for restarting). -*/
        options.add_bool("DMRG_STORE_DIIS", true);
        
        /*- Maximum number of DMRG iterations -*/
        options.add_int("DMRG_MAX_ITER", 100);
        
        /*- Which root is targeted: 1 means ground state, 2 first excited state, etc. -*/
        options.add_int("DMRG_WHICH_ROOT", 1);
        
        /*- Whether or not to use state-averaging for roots >=2 with DMRG-SCF. -*/
        options.add_bool("DMRG_STATE_AVG", true);
        
        /*- Which active space to use for DMRG calculations:
               --> input with SCF rotations (INPUT);
               --> natural orbitals (NO);
               --> localized and ordered orbitals (LOC) -*/
        options.add_str("DMRG_ACTIVE_SPACE", "INPUT", "INPUT NO LOC");
        
        /*- Whether to start the active space localization process from a random unitary or the unit matrix. -*/
        options.add_bool("DMRG_LOC_RANDOM", true);
        
        /*- Whether to calculate the DMRG-CASPT2 energy after the DMRGSCF calculations are done. -*/
        options.add_bool("DMRG_CASPT2", false);
        
        /*- CASPT2 IPEA shift -*/
        options.add_double("DMRG_IPEA", 0.0);
        
        /*- CASPT2 Imaginary shift -*/
        options.add_double("DMRG_IMAG_SHIFT", 0.0);

    }

    return true;
}


int chemps2_groupnumber(const string SymmLabel){

    int SyGroup = 0;
    bool stopFindGN = false;
    const int magic_number_max_groups_chemps2 = 8;
    do {
        if ( SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup)) == 0 ){ stopFindGN = true; }
        else { SyGroup++; }
    } while (( !stopFindGN ) && ( SyGroup < magic_number_max_groups_chemps2 ));

    (*outfile) << "Psi4 symmetry group was found to be <" << SymmLabel.c_str() << ">." << endl;
    if ( SyGroup >= magic_number_max_groups_chemps2 ){
        (*outfile) << "CheMPS2 did not recognize this symmetry group name. CheMPS2 only knows:" << endl;
        for (int cnt=0; cnt<magic_number_max_groups_chemps2; cnt++){
            (*outfile) << "   <" << (CheMPS2::Irreps::getGroupName(cnt)).c_str() << ">" << endl;
        }
        throw PSIEXCEPTION("CheMPS2 did not recognize the symmetry group name!");
    }
    return SyGroup;

}


void buildJK(SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, boost::shared_ptr<JK> myJK, boost::shared_ptr<Wavefunction> wfn){

    const int nso    = wfn->nso();
    int * nsopi      = wfn->nsopi();
    const int nmo    = wfn->nmo();
    int * nmopi      = wfn->nmopi();
    const int nirrep = wfn->nirrep();

    // nso can be different from nmo
    SharedMatrix SO_RDM;     SO_RDM = SharedMatrix( new Matrix( "SO RDM",   nirrep, nsopi, nsopi ) );
    SharedMatrix Identity; Identity = SharedMatrix( new Matrix( "Identity", nirrep, nsopi, nsopi ) );
    SharedMatrix SO_JK;       SO_JK = SharedMatrix( new Matrix( "SO JK",    nirrep, nsopi, nsopi ) );
    SharedMatrix work;         work = SharedMatrix( new Matrix( "work",     nirrep, nsopi, nmopi ) );
    
    work->gemm(false, false, 1.0, Cmat, MO_RDM, 0.0);
    SO_RDM->gemm(false, true, 1.0, work, Cmat, 0.0);
    
    std::vector<SharedMatrix> & CL = myJK->C_left();
    CL.clear();
    CL.push_back( SO_RDM );
    
    std::vector<SharedMatrix> & CR = myJK->C_right();
    CR.clear();
    Identity->identity();
    CR.push_back( Identity );
    
    myJK->set_do_J(true);
    myJK->set_do_K(true);
    myJK->set_do_wK(false);
    myJK->compute();

    SO_JK->copy( myJK->K()[0] );
    SO_JK->scale( -0.5 );
    SO_JK->add( myJK->J()[0] );
    
    work->gemm(false, false, 1.0, SO_JK, Cmat, 0.0);
    MO_JK->gemm(true, false, 1.0, Cmat, work,  0.0);

}


void copyPSIMXtoCHEMPS2MX( SharedMatrix source, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::DMRGSCFmatrix * target ){

    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set(irrep, orb1, orb2, source->get(irrep, orb1, orb2));
            }
        }
    }
    
}


void copyCHEMPS2MXtoPSIMX( CheMPS2::DMRGSCFmatrix * source, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix target ){

    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set(irrep, orb1, orb2, source->get(irrep, orb1, orb2));
            }
        }
    }
    
}


void buildQmatOCC( CheMPS2::DMRGSCFmatrix * theQmatOCC, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, boost::shared_ptr<JK> myJK, boost::shared_ptr<Wavefunction> wfn ){

    MO_RDM->zero();
    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb = 0; orb < iHandler->getNOCC(irrep); orb++){
            MO_RDM->set(irrep, orb, orb, 2.0);
        }
    }
    buildJK( MO_RDM, MO_JK, Cmat, myJK, wfn );
    copyPSIMXtoCHEMPS2MX( MO_JK, iHandler, theQmatOCC );

}


void buildQmatACT( CheMPS2::DMRGSCFmatrix * theQmatACT, CheMPS2::DMRGSCFindices * iHandler, double * DMRG1DM, SharedMatrix MO_RDM, SharedMatrix MO_JK, SharedMatrix Cmat, boost::shared_ptr<JK> myJK, boost::shared_ptr<Wavefunction> wfn ){

    MO_RDM->zero();
    const int nOrbDMRG = iHandler->getDMRGcumulative(iHandler->getNirreps());
    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        const int NOCC = iHandler->getNOCC(irrep);
        const int shift = iHandler->getDMRGcumulative(irrep);
        for (int orb1 = 0; orb1 < iHandler->getNDMRG(irrep); orb1++){
            for (int orb2 = orb1; orb2 < iHandler->getNDMRG(irrep); orb2++){
                const double value = DMRG1DM[ shift + orb1 + nOrbDMRG * ( shift + orb2 ) ];
                MO_RDM->set(irrep, NOCC+orb1, NOCC+orb2, value);
                MO_RDM->set(irrep, NOCC+orb2, NOCC+orb1, value);
            }
        }
    }
    buildJK( MO_RDM, MO_JK, Cmat, myJK, wfn );
    copyPSIMXtoCHEMPS2MX( MO_JK, iHandler, theQmatACT );

}


void buildHamDMRG( boost::shared_ptr<IntegralTransform> ints, boost::shared_ptr<MOSpace> Aorbs_ptr, CheMPS2::DMRGSCFmatrix * theTmatrix, CheMPS2::DMRGSCFmatrix * theQmatOCC, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::Hamiltonian * HamDMRG, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> wfn ){

    ints->update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints->transform_tei( Aorbs_ptr, Aorbs_ptr, Aorbs_ptr, Aorbs_ptr );
    dpd_set_default(ints->get_dpd_id());
    const int nirrep = wfn->nirrep();
    
    // Econstant and one-electron integrals
    double Econstant = wfn->molecule()->nuclear_repulsion_energy();
    for (int h = 0; h < iHandler->getNirreps(); h++){
        const int NOCC = iHandler->getNOCC(h);
        for (int froz = 0; froz < NOCC; froz++){
            Econstant += 2 * theTmatrix->get(h, froz, froz) + theQmatOCC->get(h, froz, froz);
        }
        const int shift = iHandler->getDMRGcumulative(h);
        const int NDMRG = iHandler->getNDMRG(h);
        for (int orb1 = 0; orb1 < NDMRG; orb1++){
            for (int orb2 = orb1; orb2 < NDMRG; orb2++){
                HamDMRG->setTmat( shift+orb1, shift+orb2, theTmatrix->get(h, NOCC+orb1, NOCC+orb2) + theQmatOCC->get(h, NOCC+orb1, NOCC+orb2) );
            }
        }
    }
    HamDMRG->setEconst( Econstant );
    
    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[S,S]"), ID("[S,S]"), ID("[S>=S]+"), ID("[S>=S]+"), 0, "MO Ints (SS|SS)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                HamDMRG->setVmat( p, r, q, s, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}

void buildTmatrix( CheMPS2::DMRGSCFmatrix * theTmatrix, CheMPS2::DMRGSCFindices * iHandler, boost::shared_ptr<PSIO> psio, SharedMatrix Cmat, boost::shared_ptr<Wavefunction> wfn ){

    const int nirrep = wfn->nirrep();
    const int nmo    = wfn->nmo();
    const int nTriMo = nmo * (nmo + 1) / 2;
    const int nso    = wfn->nso();
    const int nTriSo = nso * (nso + 1) / 2;
    int * mopi       = wfn->nmopi();
    int * sopi       = wfn->nsopi();
    double * work1   = new double[ nTriSo ];
    double * work2   = new double[ nTriSo ];
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_T, work1, nTriSo, 0, 0, "outfile");
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_SO_V, work2, nTriSo, 0, 0, "outfile");
    for (int n = 0; n < nTriSo; n++){ work1[n] += work2[n]; }
    delete [] work2;

    SharedMatrix soOei; soOei = SharedMatrix( new Matrix("SO OEI", nirrep, sopi, sopi) );
    SharedMatrix half;   half = SharedMatrix( new Matrix(  "Half", nirrep, mopi, sopi) );
    SharedMatrix moOei; moOei = SharedMatrix( new Matrix("MO OEI", nirrep, mopi, mopi) );

    soOei->set( work1 );
    half->gemm(true, false, 1.0, Cmat, soOei, 0.0);
    moOei->gemm(false, false, 1.0, half, Cmat, 0.0);
    delete [] work1;

    copyPSIMXtoCHEMPS2MX( moOei, iHandler, theTmatrix );

}


void fillRotatedTEI_coulomb( boost::shared_ptr<IntegralTransform> ints, boost::shared_ptr<MOSpace> OAorbs_ptr, CheMPS2::DMRGSCFintegrals * theRotatedTEI, CheMPS2::DMRGSCFindices * iHandler, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> wfn ){

    ints->update_orbitals();
    // Since we don't regenerate the SO ints, we don't call sort_so_tei, and the OEI are not updated !!!!!
    ints->transform_tei( OAorbs_ptr, OAorbs_ptr, MOSpace::all, MOSpace::all );
    dpd_set_default(ints->get_dpd_id());
    const int nirrep = wfn->nirrep();

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[Q,Q]"), ID("[A,A]"), ID("[Q>=Q]+"), ID("[A>=A]+"), 0, "MO Ints (QQ|AA)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym];
            const int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym];
                const int srel = s - K.params->soff[ssym];
                theRotatedTEI->set_coulomb( psym, qsym, rsym, ssym, prel, qrel, rrel, srel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


void fillRotatedTEI_exchange( boost::shared_ptr<IntegralTransform> ints, boost::shared_ptr<MOSpace> OAorbs_ptr, boost::shared_ptr<MOSpace> Vorbs_ptr, CheMPS2::DMRGSCFintegrals * theRotatedTEI, CheMPS2::DMRGSCFindices * iHandler, boost::shared_ptr<PSIO> psio ){

    ints->update_orbitals();
    ints->transform_tei( Vorbs_ptr, OAorbs_ptr, Vorbs_ptr, OAorbs_ptr );
    dpd_set_default(ints->get_dpd_id());

    // Two-electron integrals
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    //int buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum, int anti, const char *label);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), ID("[T,Q]"), 0, "MO Ints (TQ|TQ)");
    for(int h = 0; h < iHandler->getNirreps(); ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym] + iHandler->getNOCC(psym) + iHandler->getNDMRG(psym);
            const int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym] + iHandler->getNOCC(rsym) + iHandler->getNDMRG(rsym);
                const int srel = s - K.params->soff[ssym];
                theRotatedTEI->set_exchange( qsym, ssym, psym, rsym, qrel, srel, prel, rrel, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

}


void copyUNITARYtoPSIMX( CheMPS2::DMRGSCFunitary * unitary, CheMPS2::DMRGSCFindices * iHandler, SharedMatrix target ){

    for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
        for (int orb1 = 0; orb1 < iHandler->getNORB(irrep); orb1++){
            for (int orb2 = 0; orb2 < iHandler->getNORB(irrep); orb2++){
                target->set( irrep, orb1, orb2, unitary->getBlock(irrep)[ orb1 + iHandler->getNORB(irrep) * orb2 ] );
            }
        }
    }
    
}


void update_WFNco( CheMPS2::DMRGSCFmatrix * Coeff_orig, CheMPS2::DMRGSCFindices * iHandler, CheMPS2::DMRGSCFunitary * unitary, boost::shared_ptr<Wavefunction> wfn, SharedMatrix work1, SharedMatrix work2 ){

    copyCHEMPS2MXtoPSIMX( Coeff_orig, iHandler, work1 );
    copyUNITARYtoPSIMX( unitary, iHandler, work2 );
    wfn->Ca()->gemm(false, true, 1.0, work1, work2, 0.0);
    wfn->Cb()->copy(wfn->Ca());

}


extern "C"
SharedWavefunction dmrg(SharedWavefunction wfn, Options& options)
{

    /* This plugin is able to perform a DMRG calculation in a molecular orbital active space. */
    
    /*******************************
     *   Environment information   *
     *******************************/
    boost::shared_ptr<PSIO> psio(_default_psio_lib_); // Grab the global (default) PSIO object, for file I/O
    if (!wfn){ throw PSIEXCEPTION("SCF has not been run yet!"); }

    /*************************
     *   Fetch the options   *
     *************************/
     
    const int wfn_irrep               = options.get_int("WFN_IRREP");
    const int wfn_multp               = options.get_int("WFN_MULTP");
    int * dmrg_states                 = options.get_int_array("DMRG_STATES");
    const int ndmrg_states            = options["DMRG_STATES"].size();
    double * dmrg_econv               = options.get_double_array("DMRG_E_CONVERGENCE");
    const int ndmrg_econv             = options["DMRG_E_CONVERGENCE"].size();
    int * dmrg_maxsweeps              = options.get_int_array("DMRG_MAXSWEEPS");
    const int ndmrg_maxsweeps         = options["DMRG_MAXSWEEPS"].size();
    double * dmrg_noiseprefactors     = options.get_double_array("DMRG_NOISEPREFACTORS");
    const int ndmrg_noiseprefactors   = options["DMRG_NOISEPREFACTORS"].size();
    const bool dmrg_print_corr        = options.get_bool("DMRG_PRINT_CORR");
    const bool mps_chkpt              = options.get_bool("DMRG_CHKPT");
    int * frozen_docc                 = options.get_int_array("FROZEN_DOCC");
    int * active                      = options.get_int_array("ACTIVE");
    const double d_convergence        = options.get_double("D_CONVERGENCE");
    const bool dmrg_store_unit        = options.get_bool("DMRG_STORE_UNIT");
    const bool dmrg_do_diis           = options.get_bool("DMRG_DO_DIIS");
    const double dmrg_diis_branch     = options.get_double("DMRG_DIIS_BRANCH");
    const bool dmrg_store_diis        = options.get_bool("DMRG_STORE_DIIS");
    const int dmrg_max_iter           = options.get_int("DMRG_MAX_ITER");
    const int dmrg_which_root         = options.get_int("DMRG_WHICH_ROOT");
    const bool dmrg_state_avg         = options.get_bool("DMRG_STATE_AVG");
    const string dmrg_active_space    = options.get_str("DMRG_ACTIVE_SPACE");
    const bool dmrg_loc_random        = options.get_bool("DMRG_LOC_RANDOM");
    const bool dmrg_caspt2            = options.get_bool("DMRG_CASPT2");
    const double dmrg_ipea            = options.get_double("DMRG_IPEA");
    const double dmrg_imag_shift      = options.get_double("DMRG_IMAG_SHIFT");
    const int dmrg_num_vec_diis       = CheMPS2::DMRGSCF_numDIISvecs;
    const std::string unitaryname     = psi::get_writer_file_prefix( wfn->molecule()->name() ) + ".unitary.h5";
    const std::string diisname        = psi::get_writer_file_prefix( wfn->molecule()->name() ) + ".DIIS.h5";
    
    /****************************************
     *   Check if the input is consistent   *
     ****************************************/
    
    const int SyGroup= chemps2_groupnumber( wfn->molecule()->sym_label() );
    const int nmo    = wfn->nmo();
    const int nirrep = wfn->nirrep();
    int * orbspi     = wfn->nmopi();
    int * docc       = wfn->doccpi();
    int * socc       = wfn->soccpi();
    if ( wfn_irrep<0 )                            { throw PSIEXCEPTION("Option WFN_IRREP (integer) may not be smaller than zero!"); }
    if ( wfn_multp<1 )                            { throw PSIEXCEPTION("Option WFN_MULTP (integer) should be larger or equal to one: WFN_MULTP = (2S+1) >= 1 !"); }
    if ( ndmrg_states==0 )                        { throw PSIEXCEPTION("Option DMRG_STATES (integer array) should be set!"); }
    if ( ndmrg_econv==0 )                         { throw PSIEXCEPTION("Option DMRG_E_CONVERGENCE (double array) should be set!"); }
    if ( ndmrg_maxsweeps==0 )                     { throw PSIEXCEPTION("Option DMRG_MAXSWEEPS (integer array) should be set!"); }
    if ( ndmrg_noiseprefactors==0 )               { throw PSIEXCEPTION("Option DMRG_NOISEPREFACTORS (double array) should be set!"); }
    if ( ndmrg_states!=ndmrg_econv )              { throw PSIEXCEPTION("Options DMRG_STATES (integer array) and DMRG_ECONV (double array) should contain the same number of elements!"); }
    if ( ndmrg_states!=ndmrg_maxsweeps )          { throw PSIEXCEPTION("Options DMRG_STATES (integer array) and DMRG_MAXSWEEPS (integer array) should contain the same number of elements!"); }
    if ( ndmrg_states!=ndmrg_noiseprefactors )    { throw PSIEXCEPTION("Options DMRG_STATES (integer array) and DMRG_NOISEPREFACTORS (double array) should contain the same number of elements!"); }
    if ( options["FROZEN_DOCC"].size() != nirrep ){ throw PSIEXCEPTION("Option FROZEN_DOCC (integer array) should contain as many elements as there are irreps!"); }
    if ( options["ACTIVE"].size()      != nirrep ){ throw PSIEXCEPTION("Option ACTIVE (integer array) should contain as many elements as there are irreps!"); }
    for ( int cnt=0; cnt<ndmrg_states; cnt++ ){
       if ( dmrg_states[cnt] < 2 ){
          throw PSIEXCEPTION("Entries in DMRG_STATES (integer array) should be larger than 1!");
       }
    }
    if ( d_convergence<=0.0 )                     { throw PSIEXCEPTION("Option D_CONVERGENCE (double) must be larger than zero!"); }
    if ( dmrg_diis_branch<=0.0 )                  { throw PSIEXCEPTION("Option DMRG_DIIS_BRANCH (double) must be larger than zero!"); }
    if ( dmrg_max_iter<1 )                        { throw PSIEXCEPTION("Option DMRG_MAX_ITER (integer) must be larger than zero!"); }
    if ( dmrg_which_root<1 )                      { throw PSIEXCEPTION("Option DMRG_WHICH_ROOT (integer) must be larger than zero!"); }
    if (( dmrg_caspt2 ) && ( dmrg_ipea < 0.0 ))   { throw PSIEXCEPTION("Option DMRG_IPEA (double) must be larger than zero!"); }
    
    /*******************************************
     *   Create a CheMPS2::ConvergenceScheme   *
     *******************************************/

    CheMPS2::Initialize::Init();
    CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme( ndmrg_states );
    for (int cnt=0; cnt<ndmrg_states; cnt++){
       OptScheme->setInstruction( cnt, dmrg_states[cnt], dmrg_econv[cnt], dmrg_maxsweeps[cnt], dmrg_noiseprefactors[cnt] );
    }

    /******************************************************************************
     *   Print orbital information; check consistency of frozen_docc and active   *
     ******************************************************************************/

    int * nvirtual = new int[nirrep];
    bool virtualsOK = true;
    for (int cnt=0; cnt<nirrep; cnt++){
       nvirtual[cnt] = orbspi[cnt] - frozen_docc[cnt] - active[cnt];
       if ( nvirtual[cnt] < 0 ){ virtualsOK = false; }
    }
    (*outfile) << "wfn_irrep   = " << wfn_irrep << endl;
    (*outfile) << "wfn_multp   = " << wfn_multp << endl;
    (*outfile) << "numOrbitals = [ " << orbspi[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << orbspi[cnt];      } (*outfile) << " ]" << endl;
    (*outfile) << "R(O)HF DOCC = [ " << docc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << docc[cnt];        } (*outfile) << " ]" << endl;
    (*outfile) << "R(O)HF SOCC = [ " << socc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << socc[cnt];        } (*outfile) << " ]" << endl;
    (*outfile) << "frozen_docc = [ " << frozen_docc[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << frozen_docc[cnt]; } (*outfile) << " ]" << endl;
    (*outfile) << "active      = [ " << active[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << active[cnt];      } (*outfile) << " ]" << endl;
    (*outfile) << "virtual     = [ " << nvirtual[0];
    for (int cnt=1; cnt<nirrep; cnt++){ (*outfile) << " , " << nvirtual[cnt];    } (*outfile) << " ]" << endl;
    if ( !virtualsOK ){ throw PSIEXCEPTION("For at least one irrep: frozen_docc[ irrep ] + active[ irrep ] > numOrbitals[ irrep ]!"); }

    /*******************************************
     *   Create another bit of DMRG preamble   *
     *******************************************/
    CheMPS2::DMRGSCFindices * iHandler = new CheMPS2::DMRGSCFindices(nmo, SyGroup, frozen_docc, active, nvirtual);
    CheMPS2::DMRGSCFunitary * unitary = new CheMPS2::DMRGSCFunitary(iHandler);
    CheMPS2::DIIS * theDIIS = NULL;
    CheMPS2::DMRGSCFintegrals * theRotatedTEI = new CheMPS2::DMRGSCFintegrals( iHandler );
    const int nOrbDMRG = iHandler->getDMRGcumulative(nirrep);
    double * DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
    double * DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];
    CheMPS2::DMRGSCFmatrix * theFmatrix = new CheMPS2::DMRGSCFmatrix( iHandler ); theFmatrix->clear();
    CheMPS2::DMRGSCFmatrix * theQmatOCC = new CheMPS2::DMRGSCFmatrix( iHandler ); theQmatOCC->clear();
    CheMPS2::DMRGSCFmatrix * theQmatACT = new CheMPS2::DMRGSCFmatrix( iHandler ); theQmatACT->clear();
    CheMPS2::DMRGSCFmatrix * theTmatrix = new CheMPS2::DMRGSCFmatrix( iHandler ); theTmatrix->clear();
    CheMPS2::DMRGSCFwtilde * wmattilde  = new CheMPS2::DMRGSCFwtilde( iHandler );
    delete [] nvirtual;

    /***************************************************
     *   Create the active space Hamiltonian storage   *
     ***************************************************/

    int nElectrons = 0;
    for (int cnt=0; cnt<nirrep; cnt++){ nElectrons += 2 * docc[cnt] + socc[cnt]; }
    (*outfile) << "nElectrons  = " << nElectrons << endl;

    // Number of electrons in the active space
    int nDMRGelectrons = nElectrons;
    for (int cnt=0; cnt<nirrep; cnt++){ nDMRGelectrons -= 2 * frozen_docc[cnt]; }
    (*outfile) << "nEl. active = " << nDMRGelectrons << endl;

    // Create the CheMPS2::Hamiltonian --> fill later
    int * orbitalIrreps = new int[ nOrbDMRG ];
    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; h++){
       for (int cnt=0; cnt<active[h]; cnt++){ //Only the active space is treated with DMRG-SCF!
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }
    CheMPS2::Hamiltonian * HamDMRG = new CheMPS2::Hamiltonian(nOrbDMRG, SyGroup, orbitalIrreps);
    delete [] orbitalIrreps;

    /* Create the CheMPS2::Problem
       You can fill Ham later, as Problem only keeps a pointer to the Hamiltonian object.
       Since only doubly occupied frozen orbitals are allowed, wfn_multp and wfn_irrep do not change. */
    CheMPS2::Problem * Prob = new CheMPS2::Problem( HamDMRG , wfn_multp-1 , nDMRGelectrons , wfn_irrep );
    if ( !(Prob->checkConsistency()) ){ throw PSIEXCEPTION("CheMPS2::Problem : No Hilbert state vector compatible with all symmetry sectors!"); }
    Prob->SetupReorderD2h(); // Does nothing if group not d2h

    /**************************************
     *   Input is parsed and consistent   *
     *   Start with DMRG                  *
     **************************************/

    SharedMatrix work1; work1 = SharedMatrix( new Matrix("work1", nirrep, orbspi, orbspi) );
    SharedMatrix work2; work2 = SharedMatrix( new Matrix("work2", nirrep, orbspi, orbspi) );
    boost::shared_ptr<JK> myJK; myJK = boost::shared_ptr<JK>(new DiskJK(wfn->basisset(), options));
    myJK->set_cutoff(0.0);
    myJK->initialize();
    CheMPS2::DMRGSCFmatrix * Coeff_orig  = new CheMPS2::DMRGSCFmatrix( iHandler );
    copyPSIMXtoCHEMPS2MX(wfn->Ca(), iHandler, Coeff_orig);

    std::vector<int> OAorbs; // Occupied + active
    std::vector<int> Aorbs;  // Only active
    std::vector<int> Vorbs;  // Virtual
    std::vector<int> empty;
    for (int h = 0; h < iHandler->getNirreps(); h++){
       for (int orb = 0; orb < iHandler->getNOCC(h) + iHandler->getNDMRG(h); orb++){
          OAorbs.push_back( iHandler->getOrigNOCCstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler->getNDMRG(h); orb++){
          Aorbs.push_back( iHandler->getOrigNDMRGstart(h) + orb );
       }
       for (int orb = 0; orb < iHandler->getNVIRT(h); orb++){
          Vorbs.push_back( iHandler->getOrigNVIRTstart(h) + orb );
       }
    }
    boost::shared_ptr<MOSpace> OAorbs_ptr; OAorbs_ptr = boost::shared_ptr<MOSpace>( new MOSpace( 'Q', OAorbs, empty ) );
    boost::shared_ptr<MOSpace>  Aorbs_ptr;  Aorbs_ptr = boost::shared_ptr<MOSpace>( new MOSpace( 'S',  Aorbs, empty ) );
    boost::shared_ptr<MOSpace>  Vorbs_ptr;  Vorbs_ptr = boost::shared_ptr<MOSpace>( new MOSpace( 'T',  Vorbs, empty ) );
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back( OAorbs_ptr   );
    spaces.push_back(  Aorbs_ptr   );
    spaces.push_back(  Vorbs_ptr   );
    spaces.push_back( MOSpace::all );
    // CheMPS2 requires RHF or ROHF orbitals.
    boost::shared_ptr<IntegralTransform> ints;
    ints = boost::shared_ptr<IntegralTransform>( new IntegralTransform( wfn, spaces, IntegralTransform::Restricted ) );
    ints->set_keep_iwl_so_ints( true );
    ints->set_keep_dpd_so_ints( true );
    //ints->set_print(6);

    (*outfile) << "###########################################################" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###                       DMRG-SCF                      ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###            CheMPS2 by Sebastian Wouters             ###" << endl;
    (*outfile) << "###        https://github.com/SebWouters/CheMPS2        ###" << endl;
    (*outfile) << "###   Comput. Phys. Commun. 185 (6), 1501-1514 (2014)   ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###########################################################" << endl;
    (*outfile) << endl;
    (*outfile) << "Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;

    //Convergence variables
    double gradNorm = 1.0;
    double updateNorm = 1.0;
    double * theupdate = new double[ unitary->getNumVariablesX() ];
    for (int cnt=0; cnt<unitary->getNumVariablesX(); cnt++){ theupdate[cnt] = 0.0; }
    double * theDIISparameterVector = NULL;
    double Energy = 1e8;

    int theDIISvectorParamSize = 0;
    int maxlinsize = 0;
    for (int irrep=0; irrep<nirrep; irrep++){
        const int linsize_irrep = iHandler->getNORB(irrep);
        theDIISvectorParamSize += linsize_irrep*(linsize_irrep-1)/2;
        if (linsize_irrep>maxlinsize){ maxlinsize = linsize_irrep; }
    }

    const int nOrbDMRG_pow4    = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
    const int unitary_worksize = 4 * maxlinsize * maxlinsize;
    const int sizeWorkMem      = ( nOrbDMRG_pow4 > unitary_worksize ) ? nOrbDMRG_pow4 : unitary_worksize;
    double * mem1 = new double[sizeWorkMem];
    double * mem2 = new double[sizeWorkMem];

    CheMPS2::EdmistonRuedenberg * theLocalizer = NULL;
    if ( dmrg_active_space.compare("LOC")==0 ){ theLocalizer = new CheMPS2::EdmistonRuedenberg( HamDMRG->getVmat(), iHandler->getGroupNumber() ); }

    //Load unitary from disk
    if ( dmrg_store_unit ){
        struct stat stFileInfo;
        int intStat = stat( unitaryname.c_str(), &stFileInfo );
        if (intStat==0){ unitary->loadU( unitaryname ); }
    }

    //Load DIIS from disk
    if (( dmrg_do_diis ) && ( dmrg_store_diis )){
        struct stat stFileInfo;
        int intStat = stat( diisname.c_str(), &stFileInfo );
        if (intStat==0){
            if (theDIIS == NULL){
                theDIIS = new CheMPS2::DIIS( theDIISvectorParamSize, unitary->getNumVariablesX(), dmrg_num_vec_diis );
                theDIISparameterVector = new double[ theDIISvectorParamSize ];
            }
            theDIIS->loadDIIS( diisname );
        }
    }

    int nIterations = 0;

    /*****************************
     ***   Actual DMRG loops   ***
     *****************************/
    while ((gradNorm > d_convergence) && (nIterations < dmrg_max_iter)){

        nIterations++;

        //Update the unitary transformation
        if (unitary->getNumVariablesX() > 0){

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            unitary->updateUnitary(mem1, mem2, theupdate, true, true); //multiply = compact = true
            if (( dmrg_do_diis ) && ( updateNorm <= dmrg_diis_branch )){
                if ( dmrg_active_space.compare("NO")==0 ){
                    cout << "DIIS has started. Active space not rotated to NOs anymore!" << endl;
                }
                if ( dmrg_active_space.compare("LOC")==0 ){
                    cout << "DIIS has started. Active space not rotated to localized orbitals anymore!" << endl;
                }
                if (theDIIS == NULL){
                    theDIIS = new CheMPS2::DIIS( theDIISvectorParamSize, unitary->getNumVariablesX(), dmrg_num_vec_diis );
                    theDIISparameterVector = new double[ theDIISvectorParamSize ];
                    unitary->makeSureAllBlocksDetOne(mem1, mem2);
                }
                unitary->getLog(theDIISparameterVector, mem1, mem2);
                theDIIS->appendNew(theupdate, theDIISparameterVector);
                theDIIS->calculateParam(theDIISparameterVector);
                unitary->updateUnitary(mem1, mem2, theDIISparameterVector, false, false); //multiply = compact = false
            }

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

        }
        if (( dmrg_store_unit ) && (gradNorm!=1.0)){ unitary->saveU( unitaryname ); }
        if (( dmrg_store_diis ) && (updateNorm!=1.0) && (theDIIS!=NULL)){ theDIIS->saveDIIS( diisname ); }

        //Fill HamDMRG
        update_WFNco( Coeff_orig, iHandler, unitary, wfn, work1, work2 );
        buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
        buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
        buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );

        //Localize the active space and reorder the orbitals within each irrep based on the exchange matrix
        if (( dmrg_active_space.compare("LOC")==0 ) && (theDIIS==NULL)){ //When the DIIS has started: stop

            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            theLocalizer->Optimize( mem1, mem2, dmrg_loc_random );
            theLocalizer->FiedlerExchange(maxlinsize, mem1, mem2);
            CheMPS2::CASSCF::fillLocalizedOrbitalRotations(theLocalizer->getUnitary(), iHandler, mem1);
            unitary->rotateActiveSpaceVectors(mem1, mem2);

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());

            update_WFNco( Coeff_orig, iHandler, unitary, wfn, work1, work2 );
            buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
            buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
            buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );
            (*outfile) << "Rotated the active space to localized orbitals, sorted according to the exchange matrix." << endl;

        }

        //Do the DMRG sweeps, and calculate the 2DM
        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, mps_chkpt, psi4TMPpath);
            for (int state = 0; state < dmrg_which_root; state++){
                if (state > 0){ theDMRG->newExcitation( fabs( Energy ) ); }
                Energy = theDMRG->Solve();
                if ( dmrg_state_avg ){ // When SA-DMRGSCF: 2DM += current 2DM
                    theDMRG->calc2DMandCorrelations();
                    CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
                }
                if ((state == 0) && (dmrg_which_root > 1)){ theDMRG->activateExcitations( dmrg_which_root-1 ); }
            }
            if ( !(dmrg_state_avg) ){ // When SS-DMRGSCF: 2DM += last 2DM
                theDMRG->calc2DMandCorrelations();
                CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
            }
            if ( dmrg_print_corr ){ theDMRG->getCorrelations()->Print(); }
            if ( CheMPS2::DMRG_storeRenormOptrOnDisk ){ theDMRG->deleteStoredOperators(); }
            delete theDMRG;
            if ((dmrg_state_avg) && (dmrg_which_root > 1)){
                const double averagingfactor = 1.0 / dmrg_which_root;
                for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] *= averagingfactor; }
            }
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }

        if (( dmrg_active_space.compare("NO")==0 ) && (theDIIS==NULL)){ //When the DIIS has started: stop
            CheMPS2::CASSCF::copy_active( DMRG1DM, theFmatrix, iHandler, true );
            CheMPS2::CASSCF::block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, true, DMRG2DM ); // Unitary is updated and DMRG2DM rotated
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );     
            update_WFNco( Coeff_orig, iHandler, unitary, wfn, work1, work2 );
            buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
            buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
            (*outfile) << "Rotated the active space to natural orbitals, sorted according to the NOON." << endl;
        }

        if (dmrg_max_iter == nIterations){
            if ( dmrg_store_unit ){ unitary->saveU( unitaryname ); }
            break;
        }

        buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
        fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio, wfn );
        fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            CheMPS2::CASSCF::buildFmat( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);
            CheMPS2::CASSCF::buildWtilde(wmattilde, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);
            CheMPS2::CASSCF::augmentedHessianNR(theFmatrix, wmattilde, iHandler, unitary, theupdate, &updateNorm, &gradNorm);

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
        }
    }

    outfile->Printf("The DMRG-SCF energy = %3.10f \n", Energy);
    Process::environment.globals["CURRENT ENERGY"] = Energy;
    Process::environment.globals["DMRGSCF ENERGY"] = Energy;

    if (( dmrg_caspt2 ) && ( nIterations > 0 )){

        (*outfile) << "###########################" << endl;
        (*outfile) << "###                     ###" << endl;
        (*outfile) << "###     DMRG-CASPT2     ###" << endl;
        (*outfile) << "###                     ###" << endl;
        (*outfile) << "###########################" << endl;

        // Tmatrix, Qocc, and Qact are built in the while loop above
        outfile->Printf("Rotating to pseudocanonical orbitals");
        CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
        CheMPS2::CASSCF::block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL );
        CheMPS2::CASSCF::block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, DMRG2DM );
        CheMPS2::CASSCF::block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL );

        update_WFNco( Coeff_orig, iHandler, unitary, wfn, work1, work2 );
        buildTmatrix( theTmatrix, iHandler, psio, wfn->Ca(), wfn );
        buildQmatOCC( theQmatOCC, iHandler, work1, work2, wfn->Ca(), myJK, wfn );
        buildHamDMRG( ints, Aorbs_ptr, theTmatrix, theQmatOCC, iHandler, HamDMRG, psio, wfn );

        const int tot_dmrg_power6 = nOrbDMRG_pow4 * nOrbDMRG * nOrbDMRG;
        double * contract = new double[ tot_dmrg_power6 ];
        double * three_dm = new double[ tot_dmrg_power6 ];

        {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );

            for (int cnt = 0; cnt < nOrbDMRG_pow4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
            const string psi4TMPpath = PSIOManager::shared_object()->get_default_path();
            CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, false, psi4TMPpath); // Rotated orbital space --> do not use checkpoint
            for (int state = 0; state < dmrg_which_root; state++){
                if (state > 0){ theDMRG->newExcitation( fabs( Energy ) ); }
                const double E_CASSCF = theDMRG->Solve();
                if ((state == 0) && (dmrg_which_root > 1)){ theDMRG->activateExcitations( dmrg_which_root-1 ); }
            }
            theDMRG->calc_rdms_and_correlations( true );
            CheMPS2::CASSCF::copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM  );        // 2-RDM
            CheMPS2::CASSCF::setDMRG1DM( nDMRGelectrons, nOrbDMRG, DMRG1DM, DMRG2DM );    // 1-RDM
            buildQmatACT( theQmatACT, iHandler, DMRG1DM, work1, work2, wfn->Ca(), myJK, wfn );
            CheMPS2::CASSCF::construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
            CheMPS2::CASSCF::copy_active( theFmatrix, mem2, iHandler );                   // Fock
            //CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem2, contract );
            for ( int cnt = 0; cnt < tot_dmrg_power6; cnt++ ){ contract[ cnt ] = 0.0; }
            for ( int ham_orbz = 0; ham_orbz < nOrbDMRG; ham_orbz++ ){
                theDMRG->Diag4RDM( three_dm, ham_orbz, false );
                int size = tot_dmrg_power6;
                double f_zz = mem2[ ham_orbz + nOrbDMRG * ham_orbz ];
                int inc1 = 1;
                daxpy_( &size, &f_zz, three_dm, &inc1, contract, &inc1 ); // trace( Fock * 4-RDM )
            }
            CheMPS2::CASSCF::copy3DMover( theDMRG->get3DM(), nOrbDMRG, three_dm );        // 3-RDM --> three_dm was used as work space for the constracted 4-RDM
            if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
            if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
            delete theDMRG;

            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
       }

       fillRotatedTEI_coulomb(  ints, OAorbs_ptr, theRotatedTEI, iHandler, psio, wfn );
       fillRotatedTEI_exchange( ints, OAorbs_ptr, Vorbs_ptr,  theRotatedTEI, iHandler, psio );

       (*outfile) << "CASPT2 : Norm F - F_pseudocan = " << CheMPS2::CASSCF::deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
       double E_CASPT2 = 0.0;
       {
            std::ofstream capturing;
            std::streambuf * cout_buffer;
            string chemps2filename = outfile_name + ".chemps2";
            (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
            capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
            cout_buffer = cout.rdbuf( capturing.rdbuf() );
     
            CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract, dmrg_ipea );
            E_CASPT2 = myCASPT2->solve( dmrg_imag_shift );
            delete myCASPT2;
            
            cout.rdbuf(cout_buffer);
            capturing.close();
            std::ifstream copying;
            copying.open( chemps2filename , ios::in ); // read only
            if (copying.is_open()){
                string line;
                while( getline( copying, line ) ){ (*outfile) << line << endl; }
                copying.close();
            }
            system(("rm " + chemps2filename).c_str());
            
       }
       
       delete [] three_dm;
       delete [] contract;

       outfile->Printf("The DMRG-CASPT2 variational correction energy = %3.10f \n", E_CASPT2);
       outfile->Printf("The DMRG-CASPT2 energy = %3.10f \n", Energy + E_CASPT2);
       Process::environment.globals["CURRENT ENERGY"]    = Energy + E_CASPT2;
       Process::environment.globals["DMRGCASPT2 ENERGY"] = Energy + E_CASPT2;

    }

    delete [] mem1;
    delete [] mem2;
    delete [] theupdate;
    if (theDIISparameterVector!=NULL){ delete [] theDIISparameterVector; }
    if (theLocalizer!=NULL){ delete theLocalizer; }
    if (theDIIS!=NULL){ delete theDIIS; }
    delete Coeff_orig;

    delete wmattilde;
    delete theTmatrix;
    delete theQmatOCC;
    delete theQmatACT;
    delete theFmatrix;
    delete [] DMRG1DM;
    delete [] DMRG2DM;
    delete theRotatedTEI;
    delete unitary;
    delete iHandler;

    delete OptScheme;
    delete Prob;
    delete HamDMRG;

    return wfn;
}

}} // End Namespaces
