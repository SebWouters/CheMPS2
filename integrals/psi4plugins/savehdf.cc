#include <psi4/psi4-dec.h>
#include <psi4/libparallel/parallel.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libdpd/dpd.h>
#include <psi4/libtrans/integraltransform.h>
#include <psi4/libiwl/iwl.hpp>
#include <psi4/libciomr/libciomr.h>

#include <stdlib.h>
#include <iostream>

#include "chemps2/Irreps.h"
#include "chemps2/Hamiltonian.h"
#include "chemps2/Problem.h"

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

using namespace std;

namespace psi{ namespace savehdf{

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "SAVEHDF"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}


extern "C"
SharedWavefunction savehdf(SharedWavefunction wfn, Options& options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");

    // Grab the global (default) PSIO object, for file I/O
    std::shared_ptr<PSIO> psio(_default_psio_lib_);

    /*MoldenWriter mollie(wfn);
    mollie.write("infoForMolden.txt");*/ //Please call the MoldenWriter from the psi4 input file from now on.

    // Quickly check that there are no open shell orbitals here...
    int nirrep  = wfn->nirrep();

    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example.
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());
    
    //Readin the MO OEI in moOei & print everything
    int nmo       = wfn->nmo();
    int nIrreps   = wfn->nirrep();
    int * orbspi  = init_int_array(nirrep);
    int * docc    = init_int_array(nirrep);
    int * socc    = init_int_array(nirrep);
    for ( int h = 0; h < nirrep; ++h ){
        orbspi[h] = wfn->nmopi()[h];
        docc[h] = wfn->doccpi()[h];
        socc[h] = wfn->soccpi()[h];
    }
    
    int nTriMo = nmo * (nmo + 1) / 2;
    double *temp = new double[nTriMo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, "outfile");
    moOei.set(temp);
    
    int * orbitalIrreps = new int[nmo];

    int SyGroup = 0;
    bool stopFindGN = false;
    std::string SymmLabel = wfn->molecule()->sym_label();
    do {
        if (SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup))==0) stopFindGN = true;
        else SyGroup += 1;
    } while ((!stopFindGN) && (SyGroup<42));
    outfile->Printf( "If anything went wrong: Is ");
    outfile->Printf( SymmLabel.c_str());
    outfile->Printf( " equal to ");
    outfile->Printf( (CheMPS2::Irreps::getGroupName(SyGroup)).c_str());
    outfile->Printf( " ?\n");

    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }

    CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(nmo,SyGroup,orbitalIrreps);
    delete [] orbitalIrreps;
    
    double NuclRepulsion = wfn->molecule()->nuclear_repulsion_energy();
    Ham->setEconst(NuclRepulsion);
    double EnergyHF = NuclRepulsion;
    
    int nTot = 0;
    for(int h = 0; h < nirrep; ++h){
       for (int cnt = 0; cnt < moOei.rowspi(h); cnt++){
          for (int cnt2 = cnt; cnt2 < moOei.colspi(h); cnt2++){
             Ham->setTmat(nTot+cnt, nTot+cnt2, moOei[h][cnt][cnt2]);
          }
          if (cnt <docc[h])                            EnergyHF += 2*moOei[h][cnt][cnt];
          if ((cnt>=docc[h]) && (cnt<docc[h]+socc[h])) EnergyHF +=   moOei[h][cnt][cnt];
       }
       nTot += moOei.rowspi(h);
    }
    
    outfile->Printf( "DOCC = [ ");
    for (int cnt=0; cnt<nirrep; cnt++){
        outfile->Printf( "%2d  ", docc[cnt]);
    }
    outfile->Printf( "] \nSOCC = [ ");
    for (int cnt=0; cnt<nirrep; cnt++){
        outfile->Printf( "%2d  ", socc[cnt]);
    }
    outfile->Printf( "] \n");
 
    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                /*outfile->Printf( "%1d %1d %1d %1d %16.48f \n",
                                 p, q, r, s, K.matrix[h][pq][rs]);*/
                Ham->setVmat(p,r,q,s,K.matrix[h][pq][rs]);
                if ((p==q) && (r==s)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF += 2*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += 0.5*K.matrix[h][pq][rs];
                }
                if ((p==s) && (r==q)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF -= K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    outfile->Printf( "****  HF Energy = %16.48f \n", EnergyHF);

    outfile->Printf( "****  The debug check of Hamiltonian ****\n");
    outfile->Printf( "Econst = %16.24f \n", Ham->getEconst());
   
    double test = 0.0;
    double test2 = 0.0;
    for (int i=0; i<Ham->getL(); i++){
       for (int j=0; j<Ham->getL(); j++){
          test += Ham->getTmat(i,j);
          if (i<=j) test2 += Ham->getTmat(i,j);
       }
    }
    outfile->Printf( "1-electron integrals: Sum over all elements  : %16.24f \n", test);
    outfile->Printf( "1-electron integrals: Sum over Tij with i<=j : %16.24f \n", test2);
      
    test = 0.0;
    test2 = 0.0;
    for (int i=0; i<Ham->getL(); i++){
       for (int j=0; j<Ham->getL(); j++){
          for (int k=0; k<Ham->getL(); k++){
             for (int l=0; l<Ham->getL(); l++){
                test += Ham->getVmat(i,j,k,l);
                if ((i<=j) && (j<=k) && (k<=l)) test2 += Ham->getVmat(i,j,k,l);
             }
          }
       }
    }
    outfile->Printf( "2-electron integrals: Sum over all elements          : %16.24f \n", test);
    outfile->Printf( "2-electron integrals: Sum over Vijkl with i<=j<=k<=l : %16.24f \n", test2);
 
    cout.precision(15);
    Ham->save();

    delete Ham;
    delete [] temp;

    return wfn;
}

}} // End Namespaces
