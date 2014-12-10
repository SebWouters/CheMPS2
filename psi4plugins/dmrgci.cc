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
        
        /*- Whether or not noise is on during a particular DMRG instruction.
            If ( dmrg_noise[instruction] ) noise is added. This array should
            of course contain as many elements as DMRG_STATES -*/
        options.add_array("DMRG_NOISE");
        
        /*- The energy convergence for sweeps -*/
        options.add_double("DMRG_E_CONVERGENCE", 1e-6);
        
        /*- The noise prefactor during the initial sweeps -*/
        options.add_double("DMRG_NOISE_FACTOR", 0.03);
        
        /*- The number of sweeps when the noise is still on -*/
        options.add_int("DMRG_MAXITER_NOISE", 5);
        
        /*- The number of sweeps when the noise is turned off -*/
        options.add_int("DMRG_MAXITER_SILENT", 20);
        
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
     * TODO: Is it possible to rotate only the relevent parts of the MO space for the DMRG-CI calculation (e.g. no virtuals)?
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
    int * dmrg_noise                = options.get_int_array("DMRG_NOISE");
    const int ndmrg_noise           = options["DMRG_NOISE"].size();
    const double dmrg_e_convergence = options.get_double("DMRG_E_CONVERGENCE");
    const double dmrg_noise_factor  = options.get_double("DMRG_NOISE_FACTOR");
    const int dmrg_maxiter_noise    = options.get_int("DMRG_MAXITER_NOISE");
    const int dmrg_maxiter_silent   = options.get_int("DMRG_MAXITER_SILENT");
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
    if ( ndmrg_states!=ndmrg_noise )              { throw PSIEXCEPTION("Options DMRG_NOISE and DMRG_STATES (integer arrays) should contain the same number of elements!"); }
    if ( dmrg_e_convergence<=0.0 )                { throw PSIEXCEPTION("Option DMRG_E_CONVERGENCE (double/float) should be positive!"); }
    if ( dmrg_maxiter_noise<=0 )                  { throw PSIEXCEPTION("Option DMRG_MAXITER_NOISE (integer) should be positive!"); }
    if ( dmrg_maxiter_silent<=0 )                 { throw PSIEXCEPTION("Option DMRG_MAXITER_SILENT (integer) should be positive!"); }
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
       for (int cnt=0; cnt<active[h]; cnt++){ //Remember that only the active space is considered!
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
    
    /******************************************
     *   Input is parsed and consistent.      *
     *   Rotate the integrals and fill Ham.   *
     ******************************************/
     
    fprintf(outfile, "Start filling the active space Hamiltonian.\n");
    
    // TODO: Only rotate the active space !!!
    // TODO: Option to work in other orthonormal single particle bases than MO !
    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the LibTrans object, check out the plugin_mp2 example.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
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
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
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
       if ( dmrg_noise[cnt] ){ OptScheme->setInstruction( cnt , dmrg_states[cnt] , dmrg_e_convergence , dmrg_maxiter_noise  , dmrg_noise_factor ); }
       else                  { OptScheme->setInstruction( cnt , dmrg_states[cnt] , dmrg_e_convergence , dmrg_maxiter_silent , 0.0               ); }
    }

    // Run the DMRGCI calculation
    CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob, OptScheme);
    const double EnergyDMRG = theDMRG->Solve();
    if ( dmrg_print_corr ){
       theDMRG->calc2DMandCorrelations();
       theDMRG->getCorrelations()->Print();
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

    fprintf(outfile, "The DMRG-CI energy = %3.10f", EnergyDMRG);
    Process::environment.globals["CURRENT ENERGY"]    = EnergyDMRG;
    Process::environment.globals["DMRG TOTAL ENERGY"] = EnergyDMRG;
    
    return Success;
}

}} // End Namespaces
