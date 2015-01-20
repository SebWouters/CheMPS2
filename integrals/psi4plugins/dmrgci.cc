#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>

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
        options.add_int("DMRG_PRINT_CORR", 1);

        /*- Doubly occupied frozen orbitals for DMRGCI, per irrep. Same
            conventions as for other MR methods -*/
        options.add_array("FROZEN_DOCC");

        /*- Active space orbitals for DMRGCI, per irrep. Same conventions as for other MR methods. -*/
        options.add_array("ACTIVE");
        
    }

    return true;
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

    /* TODO: Sparse rotations for DMRG-SCF require to update the wfn coefficients (skew orbital basis)
    boost::shared_ptr<psi::Matrix> CoeffAlpha = wfn->Ca();
    boost::shared_ptr<psi::Matrix> CoeffBeta  = wfn->Cb();
    CoeffAlpha->print();
    CoeffAlpha->identity();
    wfn->Ca()->print();*/

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
    
    const int dmrg_print_corr       = options.get_int("DMRG_PRINT_CORR");
    
    int * frozen_docc               = options.get_int_array("FROZEN_DOCC");
    int * active                    = options.get_int_array("ACTIVE");
    
    // Fetching basic irrep information
    std::string SymmLabel = Process::environment.molecule()->sym_label();
    int nmo       = wfn->nmo();
    int nirrep    = wfn->nirrep();
    int * orbspi  = wfn->nmopi();
    int * docc    = wfn->doccpi();
    int * socc    = wfn->soccpi();
    
    // See if CheMPS2 finds the symmetry group
    int SyGroup = 0;
    {
        bool stopFindGN = false;
        const int magic_number_max_groups_chemps2 = 8;
        do {
            if ( SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup)) == 0 ){ stopFindGN = true; }
            else { SyGroup++; }
        } while (( !stopFindGN ) && ( SyGroup < magic_number_max_groups_chemps2 ));

        fprintf(outfile, "Psi4 symmetry group was found to be <"); fprintf(outfile, SymmLabel.c_str()); fprintf(outfile, ">.\n");
        if ( SyGroup >= magic_number_max_groups_chemps2 ){
           fprintf(outfile, "CheMPS2 did not recognize this symmetry group name. CheMPS2 only knows:\n");
           for (int cnt=0; cnt<magic_number_max_groups_chemps2; cnt++){
               fprintf(outfile, "   <"); fprintf(outfile, (CheMPS2::Irreps::getGroupName(cnt)).c_str()); fprintf(outfile, ">\n");
           }
           throw PSIEXCEPTION("CheMPS2 did not recognize the symmetry group name!");
        }
    }

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
    fprintf(outfile, "wfn_irrep   = %1d \n",wfn_irrep);
    fprintf(outfile, "wfn_multp   = %1d \n",wfn_multp);
    fprintf(outfile, "numOrbitals = [ %1d ", orbspi[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", orbspi[cnt]     ); } fprintf(outfile, "]\n");
    fprintf(outfile, "R(O)HF DOCC = [ %1d ", docc[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", docc[cnt]       ); } fprintf(outfile, "]\n");
    fprintf(outfile, "R(O)HF SOCC = [ %1d ", socc[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", socc[cnt]       ); } fprintf(outfile, "]\n");
    fprintf(outfile, "frozen_docc = [ %1d ", frozen_docc[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", frozen_docc[cnt]); } fprintf(outfile, "]\n");
    fprintf(outfile, "active      = [ %1d ", active[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", active[cnt]     ); } fprintf(outfile, "]\n");
    fprintf(outfile, "virtual     = [ %1d ", nvirtual[0]);
    for (int cnt=1; cnt<nirrep; cnt++){ fprintf(outfile, ", %1d ", nvirtual[cnt]   ); } fprintf(outfile, "]\n");
    delete [] nvirtual;
    
    // Now that users can actually see what they did wrong with frozen_docc and active, throw exception
    if ( !virtualsOK ){ throw PSIEXCEPTION("For at least one irrep: frozen_docc[ irrep ] + active[ irrep ] > numOrbitals[ irrep ]!"); }
    
    // Total number of electrons
    int nElectrons = 0;
    for (int cnt=0; cnt<nirrep; cnt++){ nElectrons += 2 * docc[cnt] + socc[cnt]; }
    fprintf(outfile, "nElectrons  = %1d \n", nElectrons);
    
    // Number of electrons in the active space
    int nDMRGelectrons = nElectrons;
    for (int cnt=0; cnt<nirrep; cnt++){ nDMRGelectrons -= 2 * frozen_docc[cnt]; }
    fprintf(outfile, "nEl. active = %1d \n", nDMRGelectrons);

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
     
    fprintf(outfile, "Start filling the active space Hamiltonian.\n");

    // CheMPS2 requires RHF or ROHF orbitals.
    // Generate only the two-electron integrals for the frozen and active spaces.
    std::vector<int> frozenActive;
    std::vector<int> empty;
    int jump = 0;
    for (int h=0; h<nirrep; h++){ // Tell the two-electron rotator which integrals we'd like
       for (int orb=0; orb<frozen_docc[h]+active[h]; orb++){ frozenActive.push_back( jump + orb ); }
       jump += orbspi[h];
    }
    boost::shared_ptr<MOSpace> frozenActive_ptr;
    frozenActive_ptr = boost::shared_ptr<MOSpace>( new MOSpace( 'S', frozenActive, empty ) );
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back( frozenActive_ptr );
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei( frozenActive_ptr, frozenActive_ptr, frozenActive_ptr, frozenActive_ptr );
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
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, outfile);
    moOei.set(temp);
    for (int h=0; h<nirrep; h++){
       // A part contributes to the constant part of the energy
       for (int froz=0; froz<frozen_docc[h]; froz++){
          Econstant += 2 * moOei[h][froz][froz];
       }
       // And a part can just be directly copied
       for (int row=0; row<active[h]; row++){
          for (int col=row; col<active[h]; col++){
             Tmat[h][row + active[h]*col] = moOei[h][frozen_docc[h]+row][frozen_docc[h]+col];
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
            const int psym = K.params->psym[p];
            const int qsym = K.params->qsym[q];
            const int prel = p - K.params->poff[psym];
            const int qrel = q - K.params->qoff[qsym];
            const int p_active = prel - frozen_docc[psym];
            const int q_active = qrel - frozen_docc[qsym];
            int p_glob = p_active; for (int prev_irr=0; prev_irr<psym; prev_irr++){ p_glob += active[prev_irr]; }
            int q_glob = q_active; for (int prev_irr=0; prev_irr<qsym; prev_irr++){ q_glob += active[prev_irr]; }
            const bool p_DMRG_AS = ((p_active>=0) && (p_active < active[psym]));
            const bool q_DMRG_AS = ((q_active>=0) && (q_active < active[qsym]));
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                const int r = K.params->colorb[h][rs][0];
                const int s = K.params->colorb[h][rs][1];
                const int rsym = K.params->rsym[r];
                const int ssym = K.params->ssym[s];
                const int rrel = r - K.params->roff[rsym];
                const int srel = s - K.params->soff[ssym];
                const int r_active = rrel - frozen_docc[rsym];
                const int s_active = srel - frozen_docc[ssym];
                const bool r_DMRG_AS = ((r_active>=0) && (r_active < active[rsym]));
                const bool s_DMRG_AS = ((s_active>=0) && (s_active < active[ssym]));

                // From the plugin to dump all integrals : Ham->setVmat(p,r,q,s,K.matrix[h][pq][rs]);
                
                if (( p_active<0 ) && ( p==q )){ // Orbital p is frozen and equals q : Coulomb
                
                   // Orbitals r and s are frozen as well : Econstant
                   if (( r_active<0 ) && (r==s)){ Econstant += 2*K.matrix[h][pq][rs]; }
                   
                   // Since p==q, this automatically should imply ssym==rsym; Orbitals r and s are active
                   if (( r_DMRG_AS ) && ( s_DMRG_AS )){ Tmat[ rsym ][ r_active + active[rsym] * s_active ] += 2*K.matrix[h][pq][rs]; }
                }
                
                if (( p_active<0 ) && ( p==s )){ // Orbital p is frozen and equals s : Exchange
                
                   // Orbitals r and q are frozen as well : Econstant
                   if (( r_active<0 ) && ( r==q )){ Econstant -= K.matrix[h][pq][rs]; }
                   
                   // Since p==s, this automatically should imply qsym==rsym; Orbitals r and q are active
                   if (( r_DMRG_AS ) && ( q_DMRG_AS )){ Tmat[ rsym ][ r_active + active[rsym] * q_active ] -= K.matrix[h][pq][rs]; }
                }
                
                if (( p_DMRG_AS ) && ( q_DMRG_AS ) && ( r_DMRG_AS ) && ( s_DMRG_AS )){ //Everything is active
                   int r_glob = r_active; for (int prev_irr=0; prev_irr<rsym; prev_irr++){ r_glob += active[prev_irr]; }
                   int s_glob = s_active; for (int prev_irr=0; prev_irr<ssym; prev_irr++){ s_glob += active[prev_irr]; }
                   Ham->setVmat( p_glob , r_glob , q_glob , s_glob , K.matrix[h][pq][rs] );
                }
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
    
    fprintf(outfile, "Finished filling the active space Hamiltonian.\n");
    fprintf(outfile, "###########################################################\n");
    fprintf(outfile, "###                                                     ###\n");
    fprintf(outfile, "###                       DMRG-CI                       ###\n");
    fprintf(outfile, "###                                                     ###\n");
    fprintf(outfile, "###            CheMPS2 by Sebastian Wouters             ###\n");
    fprintf(outfile, "###        https://github.com/SebWouters/CheMPS2        ###\n");
    fprintf(outfile, "###   Comput. Phys. Commun. 185 (6), 1501-1514 (2014)   ###\n");
    fprintf(outfile, "###                                                     ###\n");
    fprintf(outfile, "###########################################################\n");
    fprintf(outfile, "\n");
    
    CheMPS2::Initialize::Init();

    std::ofstream psi4outfile;
    std::streambuf * cout_buffer;
    if ( outfile_name != "stdout" ){
        fclose(outfile);
        outfile = NULL;
        psi4outfile.open( outfile_name.c_str() , ios::app ); // append
        cout_buffer = cout.rdbuf( psi4outfile.rdbuf() );
    }

    // The convergence scheme
    CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme( ndmrg_states );
    for (int cnt=0; cnt<ndmrg_states; cnt++){
       OptScheme->setInstruction( cnt, dmrg_states[cnt], dmrg_econv[cnt], dmrg_maxsweeps[cnt], dmrg_noiseprefactors[cnt] );
    }

    // Run the DMRGCI calculation
    CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme);
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
    if (CheMPS2::DMRG_storeMpsOnDisk){ theDMRG->deleteStoredMPS(); }
    if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
    delete theDMRG;
    delete OptScheme;
    delete Prob;
    delete Ham;
    
    if ( outfile_name != "stdout" ){
        cout.rdbuf(cout_buffer);
        psi4outfile.close();
        outfile = fopen(outfile_name.c_str(), "a");
        if (outfile == NULL){
            throw PSIEXCEPTION("PSI4: Unable to reopen output file.");
        }
    }

    fprintf(outfile, "The DMRG-CI energy = %3.10f \n", EnergyDMRG);
    Process::environment.globals["CURRENT ENERGY"]    = EnergyDMRG;
    Process::environment.globals["DMRG TOTAL ENERGY"] = EnergyDMRG;
    
    return Success;
}

}} // End Namespaces
