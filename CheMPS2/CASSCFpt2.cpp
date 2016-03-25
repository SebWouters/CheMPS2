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
#include "DMRGSCFrotations.h"
#include "Lapack.h"
#include "Cumulant.h"
#include "CASPT2.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

double CheMPS2::CASSCF::caspt2( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * theDMRGSCFoptions, const double IPEA, const double IMAG ){

   const int num_elec = Nelectrons - 2 * iHandler->getNOCCsum();
   assert( num_elec >= 0 );

   //Determine the maximum NORB(irrep) and the max_block_size for the ERI orbital rotation
   const int maxlinsize      = iHandler->getNORBmax();
   const long long fullsize  = ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize );
   const string tmp_filename = CheMPS2::defaultTMPpath + "/" + CheMPS2::DMRGSCF_eri_storage_name;
   const int dmrgsize_power4 = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
   //For (ERI rotation, update unitary, block diagonalize, orbital localization)
   const int temp_work_size = (( fullsize > CheMPS2::DMRGSCF_max_mem_eri_tfo ) ? CheMPS2::DMRGSCF_max_mem_eri_tfo : fullsize );
   const int work_mem_size  = max( max( temp_work_size , maxlinsize * maxlinsize * 4 ) , dmrgsize_power4 );
   double * mem1 = new double[ work_mem_size ];
   double * mem2 = new double[ work_mem_size ];
   const int tot_dmrg_power6 = dmrgsize_power4 * nOrbDMRG * nOrbDMRG;

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
   fillConstAndTmatDMRG( HamAS );
   DMRGSCFrotations::rotate( VMAT_ORIG, HamAS->getVmat(), NULL, 'A', 'A', 'A', 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
   double E_CASSCF = 0.0;
   double * three_dm = new double[ tot_dmrg_power6 ];
   double * contract = new double[ tot_dmrg_power6 ];

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
      E_CASSCF = theFCI->GSDavidson( inoutput );
      theFCI->Fill2RDM( inoutput, DMRG2DM );                     // 2-RDM
      theFCI->Fill3RDM( inoutput, three_dm );                    // 3-RDM
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );        // 1-RDM
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      copy_active( theFmatrix, mem2, iHandler );                 // Fock
      theFCI->Fock4RDM( inoutput, three_dm, mem2, contract );    // trace( Fock * 4-RDM )
      delete theFCI;
      delete [] inoutput;

   } else { // Do the DMRG sweeps

      for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] = 0.0; } // Clear the 2-RDM
      CheMPS2::DMRG * theDMRG = new DMRG(Prob, OptScheme);
      for (int state = 0; state < rootNum; state++){
         if (state > 0){ theDMRG->newExcitation( fabs( E_CASSCF ) ); }
         E_CASSCF = theDMRG->Solve();
         if ((state == 0) && (rootNum > 1)){ theDMRG->activateExcitations( rootNum-1 ); }
      }
      theDMRG->calc_rdms_and_correlations( true );
      copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM  );        // 2-RDM
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );          // 1-RDM
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      copy_active( theFmatrix, mem2, iHandler );                   // Fock
      //CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem2, contract );
      for ( int cnt = 0; cnt < tot_dmrg_power6; cnt++ ){ contract[ cnt ] = 0.0; }
      for ( int ham_orbz = 0; ham_orbz < nOrbDMRG; ham_orbz++ ){
         theDMRG->Diag4RDM( three_dm, ham_orbz, false );
         int size = tot_dmrg_power6;
         double f_zz = mem2[ ham_orbz + nOrbDMRG * ham_orbz ];
         int inc1 = 1;
         daxpy_( &size, &f_zz, three_dm, &inc1, contract, &inc1 ); // trace( Fock * 4-RDM )
      }
      copy3DMover( theDMRG->get3DM(), nOrbDMRG, three_dm );        // 3-RDM --> three_dm was used as work space for the constracted 4-RDM
      if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
      if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
      delete theDMRG;

   }

   delete Prob;
   delete HamAS;

   //Calculate the matrix elements needed to calculate the CASPT2 V-vector
   DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'C', 'F', 'F', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
   DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'V', 'C', 'V', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
   delete_file( tmp_filename );

   delete [] mem1;
   delete [] mem2;

   cout << "CASPT2 : Norm F - F_pseudocan = " << deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
   CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract, IPEA );
   const double E_CASPT2 = myCASPT2->solve( IMAG );

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


