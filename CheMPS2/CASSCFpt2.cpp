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

double CheMPS2::CASSCF::caspt2( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * theDMRGSCFoptions ){

   const int num_elec = Nelectrons - 2 * iHandler->getNOCCsum();
   assert( num_elec >= 0 );

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

   // Rotate to pseudocanonical orbitals
   buildTmatrix();
   buildQmatOCC();
   buildQmatACT(); // DMRG1RDM needs to be set by CASSCF::solve
   construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
   block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL );
   block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, DMRG2DM );
   block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL );

   // Fill active space Hamiltonian
   Hamiltonian * HamAS = new Hamiltonian(nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital());
   Problem * Prob = new Problem(HamAS, TwoS, num_elec, Irrep);
   Prob->SetupReorderD2h(); //Doesn't matter if the group isn't D2h, Prob checks it.
   buildQmatOCC();
   buildTmatrix();
   fillConstAndTmatDMRG(HamAS);
   DMRGSCFVmatRotations theRotator(HamOrig, iHandler);
   if (doBlockWise){ theRotator.fillVmatDMRGBlockWise(HamAS, unitary, mem1, mem2, mem3, maxBlockSize); }
   else {            theRotator.fillVmatDMRG(HamAS, unitary, mem1, mem2); }
   double Energy = 0.0;
   double * three_dm = new double[ nOrbDMRGpower4 * nOrbDMRG * nOrbDMRG ];
   double * contract = new double[ nOrbDMRGpower4 * nOrbDMRG * nOrbDMRG ];

   // Solve the active space problem
   if (( OptScheme == NULL ) && ( rootNum == 1 )){ // Do FCI

      const int nalpha = ( num_elec + TwoS ) / 2;
      const int nbeta  = ( num_elec - TwoS ) / 2;
      const double workmem = 1000.0; // 1GB
      const int verbose = 2;
      CheMPS2::FCI * theFCI = new CheMPS2::FCI( HamAS, nalpha, nbeta, Irrep, workmem, verbose );
      double * inoutput = new double[ theFCI->getVecLength(0) ];
      theFCI->ClearVector( theFCI->getVecLength(0), inoutput );
      inoutput[ theFCI->LowestEnergyDeterminant() ] = 1.0;
      Energy = theFCI->GSDavidson( inoutput );
      theFCI->Fill2RDM( inoutput, DMRG2DM );                                  // 2-RDM
      theFCI->Fill3RDM( inoutput, three_dm );                                 // 3-RDM
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );                     // 1-RDM
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      cout << "CASPT2 : Fock operator RMS deviation from block-diagonal = " << deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
      copy_active( theFmatrix, mem1, iHandler );                              // Fock
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
      theDMRG->calc_rdms_and_correlations( true );
      copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM  );                                                      // 2-RDM
      copy3DMover( theDMRG->get3DM(), nOrbDMRG, three_dm );                                                      // 3-RDM
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );                                                        // 1-RDM
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      cout << "CASPT2 : Fock operator RMS deviation from block-diagonal = " << deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
      copy_active( theFmatrix, mem1, iHandler );                                                                 // Fock
      CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem1, contract );
      if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
      if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
      delete theDMRG;

   }

   delete Prob;
   delete HamAS;
   
   //Calculate the matrix elements needed to calculate the gradient and hessian
   if (doBlockWise){ theRotator.fillRotatedTEIBlockWise(theRotatedTEI, unitary, mem1, mem2, mem3, maxBlockSize); }
   else {            theRotator.fillRotatedTEI( theRotatedTEI, unitary, mem1, mem2 ); }

   delete [] mem1;
   delete [] mem2;
   if (doBlockWise){ delete [] mem3; }
   
   CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract );
   const double E_CASPT2 = myCASPT2->solve();
   
   delete myCASPT2;
   delete [] three_dm;
   delete [] contract;
   
   return E_CASPT2;

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
   
double CheMPS2::CASSCF::deviation_from_blockdiag( DMRGSCFmatrix * matrix, const DMRGSCFindices * idx ){
   
   double rms_deviation = 0.0;
   const int n_irreps = idx->getNirreps();
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){
      const int NOCC = idx->getNOCC( irrep );
      for ( int row = 0; row < NOCC; row++ ){
         for ( int col = row + 1; col < NOCC; col++ ){
            rms_deviation += matrix->get( irrep, row, col ) * matrix->get( irrep, row, col );
            rms_deviation += matrix->get( irrep, col, row ) * matrix->get( irrep, col, row );
         }
      }
      const int NACT = idx->getNDMRG( irrep );
      for ( int row = 0; row < NACT; row++ ){
         for ( int col = row + 1; col < NACT; col++ ){
            rms_deviation += matrix->get( irrep, NOCC + row, NOCC + col ) * matrix->get( irrep, NOCC + row, NOCC + col );
            rms_deviation += matrix->get( irrep, NOCC + col, NOCC + row ) * matrix->get( irrep, NOCC + col, NOCC + row );
         }
      }
      const int NVIR = idx->getNVIRT( irrep );
      const int N_OA = NOCC + NACT;
      for ( int row = 0; row < NVIR; row++ ){
         for ( int col = row + 1; col < NVIR; col++ ){
            rms_deviation += matrix->get( irrep, N_OA + row, N_OA + col ) * matrix->get( irrep, N_OA + row, N_OA + col );
            rms_deviation += matrix->get( irrep, N_OA + col, N_OA + row ) * matrix->get( irrep, N_OA + col, N_OA + row );
         }
      }
   }
   rms_deviation = sqrt( rms_deviation );
   return rms_deviation;

}


