#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libmints/typedefs.h>
//Header above this comment contains typedef boost::shared_ptr<psi::Matrix> SharedMatrix;
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <libfock/jk.h>

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "chemps2/Irreps.h"
#include "chemps2/Hamiltonian.h"
#include "chemps2/Problem.h"
#include "chemps2/ConvergenceScheme.h"
#include "chemps2/DMRG.h"
#include "chemps2/Initialize.h"

using namespace std;

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

INIT_PLUGIN

namespace psi{ namespace dmrgci{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "DMRGCI"|| options.read_globals()) {

        /*- The DMRGCI wavefunction multiplicity in the form (2S+1) -*/
        options.add_int("WFN_MULTP", -1);

        /*- The DMRGCI wavefunction irrep uses the same conventions as PSI4. How convenient :-).
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
        options.add_array("DMRG_ECONV");
        
        /*- The maximum number of sweeps to stop an instruction
            during successive DMRG instructions -*/
        options.add_array("DMRG_MAXSWEEPS");
        
        /*- The noiseprefactors for successive DMRG instructions -*/
        options.add_array("DMRG_NOISEPREFACTORS");
        
        /*- Whether or not to print the correlation functions after the DMRG calculation -*/
        options.add_bool("DMRG_PRINT_CORR", true);
        
        /*- Whether or not to create intermediary MPS checkpoints -*/
        options.add_bool("MPS_CHKPT", false);

        /*- Doubly occupied frozen orbitals for DMRGCI, per irrep. Same
            conventions as for other MR methods -*/
        options.add_array("FROZEN_DOCC");

        /*- Active space orbitals for DMRGCI, per irrep. Same conventions as for other MR methods. -*/
        options.add_array("ACTIVE");
        
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


extern "C" PsiReturnType
dmrgci(Options &options)
{

    /*
     * This plugin is able to perform a DMRGCI calculation in a molecular orbital active space.
     * TODO: Provide a DMRGCI option to select other single-particle basis states as well (MO, SO, MP2-NO).
     */

    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if (!wfn){ throw PSIEXCEPTION("SCF has not been run yet!"); }

    /****************************************************************************
     *   Parse input and see if everything is consistent.                       *
     *   Only afterwards rotate the integrals and start with expensive stuff.   *
     ****************************************************************************/
    
    // Fetching the options
    const int wfn_irrep             = options.get_int("WFN_IRREP");
    const int wfn_multp             = options.get_int("WFN_MULTP");
    
    int * dmrg_states               = options.get_int_array("DMRG_STATES");
    const int ndmrg_states          = options["DMRG_STATES"].size();
    double * dmrg_econv             = options.get_double_array("DMRG_ECONV");
    const int ndmrg_econv           = options["DMRG_ECONV"].size();
    int * dmrg_maxsweeps            = options.get_int_array("DMRG_MAXSWEEPS");
    const int ndmrg_maxsweeps       = options["DMRG_MAXSWEEPS"].size();
    double * dmrg_noiseprefactors   = options.get_double_array("DMRG_NOISEPREFACTORS");
    const int ndmrg_noiseprefactors = options["DMRG_NOISEPREFACTORS"].size();
    
    const bool dmrg_print_corr      = options.get_bool("DMRG_PRINT_CORR");
    const bool mps_chkpt            = options.get_bool("MPS_CHKPT");
    
    int * frozen_docc               = options.get_int_array("FROZEN_DOCC");
    int * active                    = options.get_int_array("ACTIVE");
    
    // Fetching basic irrep information
    const int SyGroup = chemps2_groupnumber(Process::environment.molecule()->sym_label());
    int nmo       = wfn->nmo();
    int nirrep    = wfn->nirrep();
    int * orbspi  = wfn->nmopi();
    int * docc    = wfn->doccpi();
    int * socc    = wfn->soccpi();

    // See if the other stuff has been well initialized
    if ( wfn_irrep<0 )                            { throw PSIEXCEPTION("Option WFN_IRREP (integer) may not be smaller than zero!"); }
    if ( wfn_multp<1 )                            { throw PSIEXCEPTION("Option WFN_MULTP (integer) should be larger or equal to one: WFN_MULTP = (2S+1) >= 1 !"); }
    if ( ndmrg_states==0 )                        { throw PSIEXCEPTION("Option DMRG_STATES (integer array) should be set!"); }
    if ( ndmrg_econv==0 )                         { throw PSIEXCEPTION("Option DMRG_ECONV (double array) should be set!"); }
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

    // Check consistency of frozen_docc and active with orbspi, but wait with throwing the exception
    int * nvirtual = new int[nirrep];
    bool virtualsOK = true;
    for (int cnt=0; cnt<nirrep; cnt++){
       nvirtual[cnt] = orbspi[cnt] - frozen_docc[cnt] - active[cnt];
       if ( nvirtual[cnt] < 0 ){ virtualsOK = false; }
    }
  
    // Print basic active space information
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
    delete [] nvirtual;
    
    // Now that users can actually see what they did wrong with frozen_docc and active, throw exception
    if ( !virtualsOK ){ throw PSIEXCEPTION("For at least one irrep: frozen_docc[ irrep ] + active[ irrep ] > numOrbitals[ irrep ]!"); }
    
    // Total number of electrons
    int nElectrons = 0;
    for (int cnt=0; cnt<nirrep; cnt++){ nElectrons += 2 * docc[cnt] + socc[cnt]; }
    (*outfile) << "nElectrons  = " << nElectrons << endl;
    
    // Number of electrons in the active space
    int nDMRGelectrons = nElectrons;
    for (int cnt=0; cnt<nirrep; cnt++){ nDMRGelectrons -= 2 * frozen_docc[cnt]; }
    (*outfile) << "nEl. active = " << nDMRGelectrons << endl;

    // Create the CheMPS2 Hamiltonian
    int nDMRGorbitals = 0;
    for (int cnt=0; cnt<nirrep; cnt++){ nDMRGorbitals += active[cnt]; }
    int * orbitalIrreps = new int[ nDMRGorbitals ];
    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; h++){
       for (int cnt=0; cnt<active[h]; cnt++){ //Only the active space is treated with DMRG-CI!
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }
    CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(nDMRGorbitals, SyGroup, orbitalIrreps);
    delete [] orbitalIrreps;
    
    // Create the CheMPS2 Problem object. You can fill Ham later, as Problem only keeps a pointer to the Hamiltonian object.
    // Important remark: since only doubly occupied frozen orbitals are allowed, wfn_multp and wfn_irrep do not change.
    //                   (only Abelian symmetry groups with real-valued character tables are allowed...)
    CheMPS2::Problem * Prob = new CheMPS2::Problem( Ham , wfn_multp-1 , nDMRGelectrons , wfn_irrep );
    if ( !(Prob->checkConsistency()) ){ throw PSIEXCEPTION("CheMPS2::Problem : No Hilbert state vector compatible with all symmetry sectors!"); }
    Prob->SetupReorderD2h(); // Does nothing if group not d2h
    
    /******************************************
     *   Input is parsed and consistent.      *
     *   Rotate the integrals and fill Ham.   *
     ******************************************/
     
    (*outfile) << "Start filling the active space Hamiltonian." << endl;
    
    SharedMatrix MO_RDM; MO_RDM = SharedMatrix( new Matrix("MO frozen RDM", nirrep, orbspi, orbspi) );
    SharedMatrix MO_JK;   MO_JK = SharedMatrix( new Matrix("MO frozen JK",  nirrep, orbspi, orbspi) );
    MO_RDM->zero();
    for (int irrep = 0; irrep < nirrep; irrep++){
        for (int orb = 0; orb < frozen_docc[irrep]; orb++){
            MO_RDM->set(irrep, orb, orb, 2.0);
        }
    }
    boost::shared_ptr<JK> myJK; myJK = boost::shared_ptr<JK>( new DiskJK( wfn->basisset() ) );
    myJK->set_cutoff( 0.0 );
    myJK->initialize();
    buildJK( MO_RDM, MO_JK, wfn->Ca(), myJK, wfn );

    // CheMPS2 requires RHF or ROHF orbitals.
    // Generate only the two-electron integrals for the active space.
    std::vector<int> onlyActive;
    std::vector<int> empty;
    int jump = 0;
    for (int h=0; h<nirrep; h++){ // Tell the two-electron rotator which integrals we'd like
       for (int orb=frozen_docc[h]; orb<frozen_docc[h]+active[h]; orb++){ onlyActive.push_back( jump + orb ); }
       jump += orbspi[h];
    }
    boost::shared_ptr<MOSpace> onlyActive_ptr;
    onlyActive_ptr = boost::shared_ptr<MOSpace>( new MOSpace( 'S', onlyActive, empty ) );
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back( onlyActive_ptr );
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei( onlyActive_ptr, onlyActive_ptr, onlyActive_ptr, onlyActive_ptr );
    dpd_set_default(ints.get_dpd_id());
    
    // Constant part of the energy: due to nuclear repulsion and doubly occupied orbitals
    double Econstant = Process::environment.molecule()->nuclear_repulsion_energy();
    
    // Tmat: due to active space hopping and Coulomb/Exchange interaction with the doubly occupied orbitals
    double ** Tmat = new double * [ nirrep ];
    for (int h=0; h<nirrep; h++){ Tmat[ h ] = new double[ active[h] * active[h] ]; }
    
    // Loop over the one electron integrals
    int nTriMo = nmo * (nmo + 1) / 2;
    double * temp = new double[nTriMo];
    Matrix moOei("MO OEI", nirrep, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, "outfile");
    moOei.set(temp);
    for (int h=0; h<nirrep; h++){
       // A part contributes to the constant part of the energy
       for (int froz=0; froz<frozen_docc[h]; froz++){
          Econstant += 2 * moOei[h][froz][froz] + MO_JK->get(h, froz, froz);
       }
       // And a part can just be directly copied
       for (int row=0; row<active[h]; row++){
          for (int col=row; col<active[h]; col++){
             Tmat[h][row + active[h]*col] = moOei[h][frozen_docc[h]+row][frozen_docc[h]+col]
                                          + MO_JK->get(h, frozen_docc[h]+row, frozen_docc[h]+col);
             Tmat[h][col + active[h]*row] = Tmat[h][row + active[h]*col];
          }
       }
    }
    delete [] temp;

    /*
     * Now, loop over the DPD buffer
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    //global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[S,S]"), ID("[S,S]"), ID("[S>=S]+"), ID("[S>=S]+"), 0, "MO Ints (SS|SS)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            const int p = K.params->roworb[h][pq][0];
            const int q = K.params->roworb[h][pq][1];
            /*const int psym = K.params->psym[p]; const int prel = p - K.params->poff[psym]; and idem for qrs */
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                Ham->setVmat( p, r, q, s, K.matrix[h][pq][rs] );
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Vmat has been entirely set. Now do Tmat and Econstant
    Ham->setEconst( Econstant );
    int nJump = 0;
    for (int h=0; h<nirrep; h++){
       for (int row=0; row<active[h]; row++){
          for (int col=row; col<active[h]; col++){
             Ham->setTmat( nJump + row , nJump + col , Tmat[h][ row + active[h] * col ] );
          }
       }
       nJump += active[h];
    }
    for (int h=0; h<nirrep; h++){ delete [] Tmat[h]; }
    delete [] Tmat;
    
    (*outfile) << "Finished filling the active space Hamiltonian." << endl;
    (*outfile) << "###########################################################" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###                       DMRG-CI                       ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###            CheMPS2 by Sebastian Wouters             ###" << endl;
    (*outfile) << "###        https://github.com/SebWouters/CheMPS2        ###" << endl;
    (*outfile) << "###   Comput. Phys. Commun. 185 (6), 1501-1514 (2014)   ###" << endl;
    (*outfile) << "###                                                     ###" << endl;
    (*outfile) << "###########################################################" << endl;
    (*outfile) << endl;
    
    CheMPS2::Initialize::Init();

    std::ofstream capturing;
    std::streambuf * cout_buffer;
    string chemps2filename = outfile_name + ".chemps2";
    (*outfile) << "CheMPS2 output is temporarily written to the file " << chemps2filename << " and will be copied here." << endl;
    capturing.open( chemps2filename.c_str() , ios::trunc ); // truncate
    cout_buffer = cout.rdbuf( capturing.rdbuf() );

    // The convergence scheme
    CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme( ndmrg_states );
    for (int cnt=0; cnt<ndmrg_states; cnt++){
       OptScheme->setInstruction( cnt, dmrg_states[cnt], dmrg_econv[cnt], dmrg_maxsweeps[cnt], dmrg_noiseprefactors[cnt] );
    }

    // Run the DMRGCI calculation
    CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme, mps_chkpt);
    const double EnergyDMRG = theDMRG->Solve();
    if ( dmrg_print_corr ){
       theDMRG->calc2DMandCorrelations();
       theDMRG->getCorrelations()->Print();
       
       if (true){ // Example usage of CheMPS2's 2-RDM object
          CheMPS2::TwoDM * the2RDM = theDMRG->get2DM();
          double Energy2RDM = Ham->getEconst();
          for (int orbi=0; orbi<nDMRGorbitals; orbi++){
             for (int orbj=0; orbj<nDMRGorbitals; orbj++){
                double OneRDM_ij = 0.0;
                /* GammaA(i,j,k,l) = sum_s1,s2 < a^+_i,s1  a^+_j,s2  a_l,s2  a_k,s1 >
                   GammaB(i,j,k,l) = sum_s1,s2 (-1)^{s1-s2} < a^+_i,s1  a^+_j,s2  a_l,s2  a_k,s1 > */
                for (int orbk=0; orbk<nDMRGorbitals; orbk++){ OneRDM_ij += the2RDM->getTwoDMA_HAM(orbi, orbk, orbj, orbk); }
                OneRDM_ij = OneRDM_ij / ( nDMRGelectrons - 1 );
                Energy2RDM += Ham->getTmat(orbi, orbj) * OneRDM_ij;
                for (int orbk=0; orbk<nDMRGorbitals; orbk++){
                   for (int orbl=0; orbl<nDMRGorbitals; orbl++){
                      /* CheMPS2 uses physics notation: on the following line
                         ERI = int_dr1,dr2 orbi(r_1) orbj(r_2) orbk(r_1) orbl(r_2) / | r_1 - r_2 | */
                      Energy2RDM += 0.5 * Ham->getVmat(orbi, orbj, orbk, orbl) * the2RDM->getTwoDMA_HAM(orbi, orbj, orbk, orbl);
                   }
                }
             }
          }
          cout << "Energy( DMRG-CI )     = " << EnergyDMRG << endl;
          cout << "Energy( 2-RDM * Ham ) = " << Energy2RDM << endl;
       }
       
    }
    
    //Clean up
    if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
    delete theDMRG;
    delete OptScheme;
    delete Prob;
    delete Ham;
    
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

    outfile->Printf("The DMRG-CI energy = %3.10f \n", EnergyDMRG);
    Process::environment.globals["CURRENT ENERGY"]    = EnergyDMRG;
    Process::environment.globals["DMRG TOTAL ENERGY"] = EnergyDMRG;
    
    return Success;
}

}} // End Namespaces
