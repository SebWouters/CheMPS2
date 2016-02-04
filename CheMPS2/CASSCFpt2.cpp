/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <assert.h>

#include "CASSCF.h"
#include "DMRG.h"
#include "FCI.h"
#include "DMRGSCFVmatRotations.h"
#include "Lapack.h"
#include "Cumulant.h"
#include "CASPT2.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

double CheMPS2::CASSCF::caspt2( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * theDMRGSCFoptions, const bool g1correction ){

   // Determine the maximum NORB(irrep); and the maximum NORB(irrep) which is OK according to the cutoff.
   int maxlinsize   = 0;
   int maxlinsizeOK = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int linsize_irrep = iHandler->getNORB(irrep);
      if  (linsize_irrep > maxlinsize  )                                                         {  maxlinsize  = linsize_irrep; }
      if ((linsize_irrep > maxlinsizeOK) && (linsize_irrep <= CheMPS2::DMRGSCF_maxlinsizeCutoff)){ maxlinsizeOK = linsize_irrep; }
   }
   const bool doBlockWise = (maxlinsize <= CheMPS2::DMRGSCF_maxlinsizeCutoff) ? false : true; //Only if bigger, do we want to work blockwise
   
   // Determine the blocksize for the 2-body transformation
   int maxBlockSize = maxlinsize;
   if (doBlockWise){
      int factor   = (int) (ceil( (1.0 * maxlinsize) / CheMPS2::DMRGSCF_maxlinsizeCutoff ) + 0.01);
      maxBlockSize = max( (int) (ceil( (1.0 * maxlinsize) / factor ) + 0.01) , maxlinsizeOK ); //If a particular index can be rotated at once....
   }
   
   // Allocate 2-body rotation memory: One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
   const int maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!
   const int nOrbDMRGpower4 = nOrbDMRG*nOrbDMRG*nOrbDMRG*nOrbDMRG;
   const int sizeWorkmem = max( max( maxBSpower4 , maxlinsize*maxlinsize*4 ) , nOrbDMRGpower4 ); //For (2-body tfo, updateUnitary, calcNOON, rotate2DM, rotateUnitaryNOeigenvecs)
   double * mem1 = new double[sizeWorkmem];
   double * mem2 = new double[sizeWorkmem];
   double * mem3 = NULL;
   if (doBlockWise){ mem3 = new double[maxBSpower4]; }
   
   // The two-body rotator
   DMRGSCFVmatRotations theRotator(HamOrig, iHandler);
   
   // Load unitary from disk
   if (theDMRGSCFoptions->getStoreUnitary()){
      struct stat stFileInfo;
      int intStat = stat((theDMRGSCFoptions->getUnitaryStorageName()).c_str(),&stFileInfo);
      if (intStat==0){ unitary->loadU(theDMRGSCFoptions->getUnitaryStorageName()); }
   }

   // Fill active space Hamiltonian
   Hamiltonian * HamAS = new Hamiltonian(nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital());
   int N = Nelectrons;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ N -= 2*iHandler->getNOCC(irrep); }
   Problem * Prob = new Problem(HamAS, TwoS, N, Irrep);
   Prob->SetupReorderD2h(); //Doesn't matter if the group isn't D2h, Prob checks it.
   buildQmatOCC();
   buildTmatrix();
   fillConstAndTmatDMRG(HamAS);
   if (doBlockWise){ theRotator.fillVmatDMRGBlockWise(HamAS, unitary, mem1, mem2, mem3, maxBlockSize); }
   else {            theRotator.fillVmatDMRG(HamAS, unitary, mem1, mem2); }
   double Energy = 0.0;
   double * three_dm = new double[ nOrbDMRGpower4 * nOrbDMRG * nOrbDMRG ];
   double * contract = new double[ nOrbDMRGpower4 * nOrbDMRG * nOrbDMRG ];
   DMRGSCFmatrix * theKmatAndersson = NULL;
   if ( g1correction ){ theKmatAndersson = new DMRGSCFmatrix( iHandler ); }
   
   // Solve the active space problem
   if (( OptScheme == NULL ) && ( rootNum == 1 )){ // Do FCI

      const int nalpha = ( N + TwoS ) / 2;
      const int nbeta  = ( N - TwoS ) / 2;
      const double workmem = 1000.0; // 1GB
      const int verbose = 2;
      CheMPS2::FCI * theFCI = new CheMPS2::FCI( HamAS, nalpha, nbeta, Irrep, workmem, verbose );
      double * inoutput = new double[ theFCI->getVecLength(0) ];
      theFCI->ClearVector( theFCI->getVecLength(0), inoutput );
      inoutput[ theFCI->LowestEnergyDeterminant() ] = 1.0;
      Energy = theFCI->GSDavidson( inoutput );
      theFCI->Fill2RDM( inoutput, DMRG2DM );                                  // 2-RDM
      theFCI->Fill3RDM( inoutput, three_dm );                                 // 3-RDM
      setDMRG1DM( N, nOrbDMRG, DMRG1DM, DMRG2DM );                            // 1-RDM
      buildQmatACT(); // <--- after 1-RDM calculation
      pseudocanonical_occupied( theTmatrix, theQmatOCC, theQmatACT, unitary, mem1, mem2, iHandler );
      pseudocanonical_virtual(  theTmatrix, theQmatOCC, theQmatACT, unitary, mem1, mem2, iHandler );
      buildQmatOCC();
      buildTmatrix();
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      if ( g1correction ){
         buildKmatAndersson( theKmatAndersson );
         add_g1_to_fock( theFmatrix, theKmatAndersson, mem1, mem2, DMRG1DM, iHandler );
      }
      copy_active( theFmatrix, iHandler, mem1 );                              // Fock
      theFCI->Fock4RDM( inoutput, three_dm, mem1, contract );                 // trace( Fock * 4-RDM )
      delete theFCI;
      delete [] inoutput;

   } else { // Do the DMRG sweeps

      for (int cnt = 0; cnt < nOrbDMRGpower4; cnt++){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM
      CheMPS2::DMRG * theDMRG = new DMRG(Prob, OptScheme);
      for (int state = 0; state < rootNum; state++){
         if (state > 0){ theDMRG->newExcitation( fabs( Energy ) ); }
         Energy = theDMRG->Solve();
         if ((state == 0) && (rootNum > 1)){ theDMRG->activateExcitations( rootNum-1 ); }
      }
      const double calculate_3rdm = true;
      theDMRG->calc_rdms_and_correlations( calculate_3rdm );
      copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM  );                                                      // 2-RDM
      copy3DMover( theDMRG->get3DM(), nOrbDMRG, three_dm );                                                      // 3-RDM
      setDMRG1DM( N, nOrbDMRG, DMRG1DM, DMRG2DM );                                                               // 1-RDM
      buildQmatACT(); // <--- after 1-RDM calculation
      pseudocanonical_occupied( theTmatrix, theQmatOCC, theQmatACT, unitary, mem1, mem2, iHandler );
      pseudocanonical_virtual(  theTmatrix, theQmatOCC, theQmatACT, unitary, mem1, mem2, iHandler );
      buildQmatOCC();
      buildTmatrix();
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      if ( g1correction ){
         buildKmatAndersson( theKmatAndersson );
         add_g1_to_fock( theFmatrix, theKmatAndersson, mem1, mem2, DMRG1DM, iHandler );
      }
      copy_active( theFmatrix, iHandler, mem1 );                                                                 // Fock
      CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem1, contract ); // trace( Fock * 4-RDM )
      if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
      if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
      delete theDMRG;

   }
   
   const double E_CASSCF = Energy;
   delete Prob;
   delete HamAS;
   
   //Calculate the matrix elements needed to calculate the gradient and hessian
   if (doBlockWise){ theRotator.fillRotatedTEIBlockWise(theRotatedTEI, unitary, mem1, mem2, mem3, maxBlockSize); }
   else {            theRotator.fillRotatedTEI( theRotatedTEI, unitary, mem1, mem2 ); }

   delete [] mem1;
   delete [] mem2;
   if (doBlockWise){ delete [] mem3; }
   if ( g1correction ){ delete theKmatAndersson; }
   
   CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract );
   
   delete myCASPT2;
   delete [] three_dm;
   delete [] contract;
   
   return Energy;

}

void CheMPS2::CASSCF::construct_fock( DMRGSCFmatrix * Fock, const DMRGSCFmatrix * Tmat, const DMRGSCFmatrix * Qocc, const DMRGSCFmatrix * Qact, const DMRGSCFindices * idx ){

   const int n_irreps = idx->getNirreps();
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){
      const int NORB = idx->getNORB( irrep );
      for (int row = 0; row < NORB; row++){
         for (int col = 0; col < NORB; col++){
            Fock->set( irrep, row, col, Tmat->get( irrep, row, col )
                                      + Qocc->get( irrep, row, col )
                                      + Qact->get( irrep, row, col ) );
         }
      }
   }

}

void CheMPS2::CASSCF::add_g1_to_fock( DMRGSCFmatrix * Fock, DMRGSCFmatrix * Kmat, double * work1, double * work2, double * one_rdm, const DMRGSCFindices * idx ){

   const int n_irreps = idx->getNirreps();

   // g1 = - 1/4 [ Dkd + dKD ]
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NORB = idx->getNORB( irrep );
      if ( NORB > 0 ){

         // work1 = D
         const int NOCC = idx->getNOCC( irrep );
         const int NACT = idx->getNDMRG( irrep );
         const int NVIR = idx->getNVIRT( irrep );
         const int N_OA = NOCC + NACT;
         const int JUMP = idx->getDMRGcumulative( irrep );
         const int LAS  = idx->getDMRGcumulative( n_irreps );
         for ( int cnt = 0; cnt < NORB * NORB; cnt++ ){ work1[ cnt ] = 0.0; }
         for ( int occ = 0; occ < NOCC; occ++ ){ work1[ occ + NORB * occ ] = 2.0; }
         for ( int row = 0; row < NACT; row++ ){
            for ( int col = 0; col < NACT; col++ ){
               work1[ NOCC + row + NORB * ( NOCC + col ) ] = one_rdm[ JUMP + row + LAS * ( JUMP + col ) ];
            }
         }

         // work2 = DK
         int linsize = NORB;
         char notrans = 'N';
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &alpha, work1, &linsize, Kmat->getBlock( irrep ), &linsize, &beta, work2, &linsize );

         // work3 = DK
         int copysize = NORB * NORB;
         int inc1 = 1;
         double * work3 = work2 + copysize;
         dcopy_( &copysize, work2, &inc1, work3, &inc1 );

         // work3 = DK(2-D) = 2 * work3 - work2 * work1
         alpha = -1.0;
         beta = 2.0;
         dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &alpha, work2, &linsize, work1, &linsize, &beta, work3, &linsize );

         // Not an approximation: the occupied-occupied and virtual-virtual block of work3 should be 0
         for ( int row = 0; row < NOCC; row++ ){
            for ( int col = 0; col < NOCC; col++ ){
               work3[ row + NORB * col ] = 0.0;
            }
         }
         for ( int row = 0; row < NVIR; row++ ){
            for ( int col = 0; col < NVIR; col++ ){
               work3[ N_OA + row + NORB * ( N_OA + col ) ] = 0.0;
            }
         }

         // Fock = Fock - 0.25 ( DK(2-D) + (2-D)KD )
         double * Fblock = Fock->getBlock( irrep );
         for ( int row = 0; row < NORB; row++ ){
            for ( int col = 0; col < NORB; col++ ){
               Fblock[ row + NORB * col ] -= 0.25 * ( work3[ row + NORB * col ] + work3[ col + NORB * row ] );
            }
         }
      }
   }

}

void CheMPS2::CASSCF::copy_active( const DMRGSCFmatrix * mat, const DMRGSCFindices * idx, double * result ){

   const int n_irreps = idx->getNirreps();
   const int LAS      = idx->getDMRGcumulative( n_irreps );
   for ( int cnt = 0; cnt < LAS * LAS; cnt++ ){ result[ cnt ] = 0.0; }

   int passed = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){
      const int NOCC = idx->getNOCC(  irrep );
      const int NACT = idx->getNDMRG( irrep );
      for ( int row = 0; row < NACT; row++ ){
         for ( int col = 0; col < NACT; col++ ){
            result[ passed + row + LAS * ( passed + col ) ] = mat->get( irrep, NOCC + row, NOCC + col );
         }
      }
      passed += NACT;
   }

   assert( passed == LAS );

}



