/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

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
#include "Lapack.h"
#include "DMRGSCFrotations.h"
#include "EdmistonRuedenberg.h"
#include "Davidson.h"
#include "FCI.h"
#include "MPIchemps2.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

void CheMPS2::CASSCF::delete_file( const string filename ){

   struct stat file_info;
   const int thestat = stat( filename.c_str(), &file_info );
   if ( thestat == 0 ){
      const string temp = "rm " + filename;
      int info = system( temp.c_str() );
      cout << "Info on system( " << temp << " ) = " << info << endl;
   } else {
      cout << "No file " << filename << " found." << endl;
   }

}

double CheMPS2::CASSCF::solve( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * scf_options ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   const int num_elec = Nelectrons - 2 * iHandler->getNOCCsum();
   assert( num_elec >= 0 );
   assert(( OptScheme != NULL ) || ( rootNum == 1 ));

   // Convergence variables
   double gradNorm = 1.0;
   double updateNorm = 1.0;
   double * gradient = NULL;
   if ( am_i_master ){
      gradient = new double[ unitary->getNumVariablesX() ];
      for ( int cnt = 0; cnt < unitary->getNumVariablesX(); cnt++ ){ gradient[ cnt ] = 0.0; }
   }
   double * diis_vec = NULL;
   double Energy = 1e8;

   // The CheMPS2::Problem for the inner DMRG calculation
   Hamiltonian * HamDMRG = new Hamiltonian( nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital() );
   Problem * Prob = new Problem( HamDMRG, TwoS, num_elec, Irrep );
   Prob->SetupReorderD2h(); // Doesn't matter if the group isn't D2h, Prob checks it.

   // Determine the maximum NORB( irrep ) and the max_block_size for the ERI orbital rotation
   const int maxlinsize      = iHandler->getNORBmax();
   const long long fullsize  = ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize );
   const string tmp_filename = tmp_folder + "/" + CheMPS2::DMRGSCF_eri_storage_name;
   const int dmrgsize_power4 = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
   // For ( ERI rotation, update unitary, block diagonalize, orbital localization )
   DMRGSCFintegrals * theRotatedTEI = new DMRGSCFintegrals( iHandler );
   DMRGSCFwtilde * wmattilde = new DMRGSCFwtilde( iHandler );
   const int temp_work_size = (( fullsize > CheMPS2::DMRGSCF_max_mem_eri_tfo ) ? CheMPS2::DMRGSCF_max_mem_eri_tfo : fullsize );
   const int work_mem_size = max( max( temp_work_size , maxlinsize * maxlinsize * 4 ) , dmrgsize_power4 );
   double * mem1 = new double[ work_mem_size ];
   double * mem2 = new double[ work_mem_size ];

   // Only MASTER process has Edmiston-Ruedenberg active space localizer
   EdmistonRuedenberg * theLocalizer = NULL;
   if (( scf_options->getWhichActiveSpace() >= 2 ) && ( am_i_master )){ theLocalizer = new EdmistonRuedenberg( HamDMRG->getVmat(), iHandler->getGroupNumber() ); }

   // Load unitary from disk
   if ( scf_options->getStoreUnitary() ){
      if ( am_i_master ){
         struct stat file_info;
         int master_stat = stat( (scf_options->getUnitaryStorageName()).c_str(), &file_info );
         if ( master_stat == 0 ){ unitary->loadU( scf_options->getUnitaryStorageName() ); }
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      unitary->broadcast( MPI_CHEMPS2_MASTER ); // Unitary was set to the identity at construction
      #endif
   }

   // Only master has DIIS; load DIIS from disk
   DIIS * diis = NULL;
   if (( scf_options->getDoDIIS() ) && ( scf_options->getStoreDIIS() ) && ( am_i_master )){
      struct stat file_info;
      int master_stat = stat( (scf_options->getDIISStorageName()).c_str(), &file_info );
      if ( master_stat == 0 ){
         const int diis_vec_size = iHandler->getROTparamsize();
         diis = new DIIS( diis_vec_size, unitary->getNumVariablesX(), scf_options->getNumDIISVecs() );
         diis->loadDIIS( scf_options->getDIISStorageName() );
         diis_vec = new double[ diis_vec_size ];
      }
   }

   int nIterations = 0;

   /*******************************
   ***   Actual DMRGSCF loops   ***
   *******************************/
   while (( gradNorm > scf_options->getGradientThreshold() ) && ( nIterations < scf_options->getMaxIterations() )){

      nIterations++;

      // Update the unitary transformation
      if (( unitary->getNumVariablesX() > 0 ) && ( am_i_master )){

         unitary->updateUnitary( mem1, mem2, gradient, true, true ); //multiply = compact = true

         if (( scf_options->getDoDIIS() ) && ( updateNorm <= scf_options->getDIISGradientBranch() )){
            if ( scf_options->getWhichActiveSpace() == 1 ){
               cout << "DMRGSCF::solve : DIIS has started. Active space not rotated to NOs anymore!" << endl;
            }
            if ( scf_options->getWhichActiveSpace() == 2 ){
               cout << "DMRGSCF::solve : DIIS has started. Active space not rotated to localized orbitals anymore!" << endl;
            }
            if ( scf_options->getWhichActiveSpace() == 3 ){
               cout << "DMRGSCF::solve : DIIS has started. Active space not ordered according to the Fiedler vector of the exchange matrix anymore!" << endl;
            }
            if ( diis == NULL ){
               const int diis_vec_size = iHandler->getROTparamsize();
               diis = new DIIS( diis_vec_size, unitary->getNumVariablesX(), scf_options->getNumDIISVecs() );
               diis_vec = new double[ diis_vec_size ];
               unitary->makeSureAllBlocksDetOne( mem1, mem2 );
            }
            unitary->getLog( diis_vec, mem1, mem2 );
            diis->appendNew( gradient, diis_vec );
            diis->calculateParam( diis_vec );
            unitary->updateUnitary( mem1, mem2, diis_vec, false, false ); //multiply = compact = false
         }
      }
      if (( am_i_master ) && ( scf_options->getStoreUnitary() ) && ( gradNorm != 1.0 )){ unitary->saveU( scf_options->getUnitaryStorageName() ); }
      if (( am_i_master ) && ( scf_options->getStoreDIIS() ) && ( updateNorm != 1.0 ) && ( diis != NULL )){ diis->saveDIIS( scf_options->getDIISStorageName() ); }
      int master_diis = (( diis != NULL ) ? 1 : 0 );
      #ifdef CHEMPS2_MPI_COMPILATION
      MPIchemps2::broadcast_array_int( &master_diis, 1, MPI_CHEMPS2_MASTER );
      unitary->broadcast( MPI_CHEMPS2_MASTER );
      #endif

      // Fill HamDMRG
      buildTmatrix();
      buildQmatOCC();
      fillConstAndTmatDMRG( HamDMRG );
      if ( am_i_master ){
         DMRGSCFrotations::rotate( VMAT_ORIG, HamDMRG->getVmat(), NULL, 'A', 'A', 'A', 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      HamDMRG->getVmat()->broadcast( MPI_CHEMPS2_MASTER );
      #endif

      // Localize the active space and reorder the orbitals within each irrep based on the exchange matrix
      if (( scf_options->getWhichActiveSpace() == 2 ) && ( master_diis == 0 )){ // When the DIIS has started: stop
         if ( am_i_master ){
            theLocalizer->Optimize(mem1, mem2, scf_options->getStartLocRandom());
            theLocalizer->FiedlerExchange(maxlinsize, mem1, mem2);
            fillLocalizedOrbitalRotations(theLocalizer->getUnitary(), iHandler, mem1);
            unitary->rotateActiveSpaceVectors(mem1, mem2);
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         unitary->broadcast( MPI_CHEMPS2_MASTER );
         #endif
         buildTmatrix(); //With an updated unitary, the Qocc, Tmat, and HamDMRG objects need to be updated as well.
         buildQmatOCC();
         fillConstAndTmatDMRG( HamDMRG );
         if ( am_i_master ){
            DMRGSCFrotations::rotate( VMAT_ORIG, HamDMRG->getVmat(), NULL, 'A', 'A', 'A', 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
            cout << "DMRGSCF::solve : Rotated the active space to localized orbitals, sorted according to the exchange matrix." << endl;
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         HamDMRG->getVmat()->broadcast( MPI_CHEMPS2_MASTER );
         #endif
      }

      // Reorder the orbitals based on the Fiedler vector of the exchange matrix
      if (( scf_options->getWhichActiveSpace() == 3 ) && ( master_diis == 0 )){ // When the DIIS has started: stop
         int * dmrg2ham = new int[ nOrbDMRG ];
         if ( am_i_master ){
            theLocalizer->FiedlerGlobal( dmrg2ham );
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         MPIchemps2::broadcast_array_int( dmrg2ham, nOrbDMRG, MPI_CHEMPS2_MASTER );
         #endif
         Prob->setup_reorder_custom( dmrg2ham );
         delete [] dmrg2ham;
      }

      if (( OptScheme == NULL ) && ( rootNum == 1 )){ // Do FCI, and calculate the 2DM

         if ( am_i_master ){
            const int nalpha = ( num_elec + TwoS ) / 2;
            const int nbeta  = ( num_elec - TwoS ) / 2;
            const double workmem = 1000.0; // 1GB
            const int verbose = 2;
            CheMPS2::FCI * theFCI = new CheMPS2::FCI( HamDMRG, nalpha, nbeta, Irrep, workmem, verbose );
            double * inoutput = new double[ theFCI->getVecLength(0) ];
            theFCI->ClearVector( theFCI->getVecLength(0), inoutput );
            inoutput[ theFCI->LowestEnergyDeterminant() ] = 1.0;
            Energy = theFCI->GSDavidson( inoutput );
            theFCI->Fill2RDM( inoutput, DMRG2DM );
            delete theFCI;
            delete [] inoutput;
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         MPIchemps2::broadcast_array_double( &Energy, 1, MPI_CHEMPS2_MASTER );
         MPIchemps2::broadcast_array_double( DMRG2DM, dmrgsize_power4, MPI_CHEMPS2_MASTER );
         #endif

      } else { // Do the DMRG sweeps, and calculate the 2DM

         assert( OptScheme != NULL );
         for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] = 0.0; } // Clear the 2-RDM ( to allow for state-averaged calculations )
         DMRG * theDMRG = new DMRG( Prob, OptScheme, CheMPS2::DMRG_storeMpsOnDisk, tmp_folder );
         for ( int state = 0; state < rootNum; state++ ){
            if ( state > 0 ){ theDMRG->newExcitation( fabs( Energy ) ); }
            Energy = theDMRG->Solve();
            if ( scf_options->getStateAveraging() ){ // When SA-DMRGSCF: 2DM += current 2DM
               theDMRG->calc2DMandCorrelations();
               copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
            }
            if (( state == 0 ) && ( rootNum > 1 )){ theDMRG->activateExcitations( rootNum - 1 ); }
         }
         if ( !( scf_options->getStateAveraging() )){ // When SS-DMRGSCF: 2DM += last 2DM
            theDMRG->calc2DMandCorrelations();
            copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
         }
         if (( scf_options->getDumpCorrelations() ) && ( am_i_master )){ theDMRG->getCorrelations()->Print(); }
         if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
         if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
         delete theDMRG;
         if (( scf_options->getStateAveraging() ) && ( rootNum > 1 )){
            const double averagingfactor = 1.0 / rootNum;
            for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] *= averagingfactor; }
         }

      }
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );

      // Possibly rotate the active space to the natural orbitals
      if (( scf_options->getWhichActiveSpace() == 1 ) && ( master_diis == 0 )){ // When the DIIS has started: stop
         copy_active( DMRG1DM, theQmatWORK, iHandler, true );
         block_diagonalize( 'A', theQmatWORK, unitary, mem1, mem2, iHandler, true, DMRG2DM, NULL, NULL ); // Unitary is updated and DMRG2DM rotated
         setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );
         buildTmatrix(); // With an updated unitary, the Qocc and Tmat matrices need to be updated as well.
         buildQmatOCC();
         if ( am_i_master ){ cout << "DMRGSCF::solve : Rotated the active space to natural orbitals, sorted according to the NOON." << endl; }
      }

      if (( nIterations == scf_options->getMaxIterations() ) && ( am_i_master ) && ( scf_options->getStoreUnitary() )){
         unitary->saveU( scf_options->getUnitaryStorageName() );
      }

      // Calculate the matrix elements needed to calculate the gradient and hessian
      buildQmatACT();
      if ( am_i_master ){
         DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'C', 'F', 'F', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
         DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'V', 'C', 'V', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
         buildFmat(  theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM );
         buildWtilde( wmattilde, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM );
         augmentedHessianNR( theFmatrix, wmattilde, iHandler, unitary, gradient, &updateNorm, &gradNorm ); // On return the gradient contains the update
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      MPIchemps2::broadcast_array_double( &updateNorm, 1, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_double( &gradNorm,   1, MPI_CHEMPS2_MASTER );
      #endif

   }

   delete [] mem1;
   delete [] mem2;
   delete theRotatedTEI;
   delete wmattilde;
   if ( am_i_master ){ delete_file( tmp_filename ); }

   delete Prob;
   delete HamDMRG;
   if ( gradient != NULL ){ delete [] gradient; }
   if ( diis_vec != NULL ){ delete [] diis_vec; }
   if ( diis != NULL ){ delete diis; }
   if ( theLocalizer != NULL ){ delete theLocalizer; }

   successful_solve = true;
   return Energy;

}

void CheMPS2::CASSCF::augmentedHessianNR( DMRGSCFmatrix * localFmat, DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * theupdate, double * updateNorm, double * gradNorm ){

   /* A good read to understand
            (1) how the augmented Hessian arises from a rational function optimization
            (2) where the parameter lambda in Eq. (22) of Yanai, IJQC 109, 2178-2190 (2009) comes from
            (3) why the smallest algebraic eigenvalue + corresponding eigenvector should be retained for minimizations
      Banerjee, Adams, Simons, Shepard, "Search for stationary points on surfaces",
      J. Phys. Chem. 1985, volume 89, pages 52-57, doi:10.1021/j100247a015  */

   //Calculate the gradient
   const int x_linearlength = localUmat->getNumVariablesX();
   gradNorm[ 0 ] = construct_gradient( localFmat, localIdx, theupdate );

   //Find the lowest eigenvalue and corresponding eigenvector of the augmented hessian
   {
      Davidson deBoskabouter( x_linearlength + 1, CheMPS2::DAVIDSON_NUM_VEC,
                                                  CheMPS2::DAVIDSON_NUM_VEC_KEEP,
                                                  CheMPS2::DAVIDSON_FCI_RTOL,
                                                  CheMPS2::DAVIDSON_PRECOND_CUTOFF, false ); // No debug printing
      double ** whichpointers = new double*[ 2 ];
      char instruction = deBoskabouter.FetchInstruction( whichpointers );
      assert( instruction == 'A' );
      diag_hessian( localFmat, localwtilde, localIdx, whichpointers[ 1 ] );
      whichpointers[ 1 ][ x_linearlength ] = 0.0;
      for ( int cnt = 0; cnt < x_linearlength; cnt++ ){ // Initial guess = [ -gradient / diag(hessian) , 1 ]
         const double denom = ( whichpointers[ 1 ][ cnt ] > CheMPS2::DAVIDSON_PRECOND_CUTOFF ) ? whichpointers[ 1 ][ cnt ] : CheMPS2::DAVIDSON_PRECOND_CUTOFF;
         whichpointers[ 0 ][ cnt ] = - theupdate[ cnt ] / denom;
      }
      whichpointers[ 0 ][ x_linearlength ] = 1.0;
      instruction = deBoskabouter.FetchInstruction( whichpointers );
      while ( instruction == 'B' ){
         augmented_hessian( localFmat, localwtilde, localIdx, whichpointers[ 0 ], whichpointers[ 1 ], theupdate, x_linearlength );
         instruction = deBoskabouter.FetchInstruction( whichpointers );
      }
      assert( instruction == 'C' );
      double scalar = 1.0 / whichpointers[ 0 ][ x_linearlength ];
      cout << "DMRGSCF::augmentedHessianNR : Augmented Hessian update found with " << deBoskabouter.GetNumMultiplications() << " Davidson iterations." << endl;
      if ( CheMPS2::DMRGSCF_debugPrint ){
         cout << "DMRGSCF::augmentedHessianNR : Lowest eigenvalue = " << whichpointers[ 1 ][ 0 ] << endl;
         cout << "DMRGSCF::augmentedHessianNR : The last number of the eigenvector (which will be rescaled to one) = " << scalar << endl;
      }
      for ( int cnt = 0; cnt < x_linearlength; cnt++ ){ theupdate[ cnt ] = scalar * whichpointers[ 0 ][ cnt ]; }
      delete [] whichpointers;
   }

   //Calculate the update norm
   updateNorm[ 0 ] = 0.0;
   for ( int cnt = 0; cnt < x_linearlength; cnt++ ){ updateNorm[ 0 ] += theupdate[ cnt ] * theupdate[ cnt ]; }
   updateNorm[ 0 ] = sqrt( updateNorm[ 0 ] );
   cout << "DMRGSCF::augmentedHessianNR : Norm of the update = " << updateNorm[ 0 ] << endl;

}

double CheMPS2::CASSCF::construct_gradient( DMRGSCFmatrix * Fmatrix, const DMRGSCFindices * idx, double * gradient ){

   const int n_irreps = idx->getNirreps();

   int jump = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NORB = idx->getNORB( irrep );
      const int NOCC = idx->getNOCC( irrep );
      const int NACT = idx->getNDMRG( irrep );
      const int NVIR = idx->getNVIRT( irrep );
      const int N_OA = NOCC + NACT;
      double * FMAT  = Fmatrix->getBlock( irrep );

      // Type 0: act-occ
      for ( int occ = 0; occ < NOCC; occ++ ){
         for ( int act = 0; act < NACT; act++ ){
            gradient[ jump + act + NACT * occ ] = 2 * ( FMAT[ NOCC + act + NORB * occ ] - FMAT[ occ + NORB * ( NOCC + act ) ] );
         }
      }
      jump += NOCC * NACT;

      // Type 1: vir-act
      for ( int act = 0; act < NACT; act++ ){
         for ( int vir = 0; vir < NVIR; vir++ ){
            gradient[ jump + vir + NVIR * act ] = 2 * ( FMAT[ N_OA + vir + NORB * ( NOCC + act ) ] - FMAT[ NOCC + act + NORB * ( N_OA + vir ) ] );
         }
      }
      jump += NACT * NVIR;

      // Type 2: vir-occ
      for ( int occ = 0; occ < NOCC; occ++ ){
         for ( int vir = 0; vir < NVIR; vir++ ){
            gradient[ jump + vir + NVIR * occ ] = 2 * ( FMAT[ N_OA + vir + NORB * occ ] - FMAT[ occ + NORB * ( N_OA + vir ) ] );
         }
      }
      jump += NOCC * NVIR;
   }

   double gradient_norm = 0.0;
   for ( int cnt = 0; cnt < jump; cnt++ ){ gradient_norm += gradient[ cnt ] * gradient[ cnt ]; }
   gradient_norm = sqrt( gradient_norm );
   cout << "DMRGSCF::construct_gradient : Norm of the gradient = " << gradient_norm << endl;
   return gradient_norm;

}

void CheMPS2::CASSCF::augmented_hessian( DMRGSCFmatrix * Fmatrix, DMRGSCFwtilde * Wtilde, const DMRGSCFindices * idx, double * origin, double * target, double * gradient, const int linsize ){

   for ( int cnt = 0; cnt < linsize; cnt++ ){
      target[ cnt ] = gradient[ cnt ] * origin[ linsize ];
   }
   add_hessian( Fmatrix, Wtilde, idx, origin, target );
   target[ linsize ] = 0.0;
   for ( int cnt = 0; cnt < linsize; cnt++ ){
      target[ linsize ] += gradient[ cnt ] * origin[ cnt ];
   }

}

void CheMPS2::CASSCF::diag_hessian( DMRGSCFmatrix * Fmatrix, const DMRGSCFwtilde * Wtilde, const DMRGSCFindices * idx, double * diagonal ){

   const int n_irreps = idx->getNirreps();

   int jump = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NORB = idx->getNORB( irrep );
      const int NOCC = idx->getNOCC( irrep );
      const int NACT = idx->getNDMRG( irrep );
      const int NVIR = idx->getNVIRT( irrep );
      const int N_OA = NOCC + NACT;
      double * FMAT = Fmatrix->getBlock( irrep );

      for ( int occ = 0; occ < NOCC; occ++ ){
         const double F_occ = FMAT[ occ * ( NORB + 1 ) ];
         for ( int act = 0; act < NACT; act++ ){
            const double F_act = FMAT[ ( NOCC + act ) * ( NORB + 1 ) ];
            diagonal[ jump + act + NACT * occ ] = - 2 * ( F_occ + F_act ) + ( Wtilde->get( irrep, irrep, NOCC + act, occ, NOCC + act, occ )
                                                                            - Wtilde->get( irrep, irrep, NOCC + act, occ, occ, NOCC + act )
                                                                            - Wtilde->get( irrep, irrep, occ, NOCC + act, NOCC + act, occ ) 
                                                                            + Wtilde->get( irrep, irrep, occ, NOCC + act, occ, NOCC + act ) );
         }
      }
      jump += NOCC * NACT;

      for ( int act = 0; act < NACT; act++ ){
         const double F_act = FMAT[ ( NOCC + act ) * ( NORB + 1 ) ];
         for ( int vir = 0; vir < NVIR; vir++ ){
            const double F_vir = FMAT[ ( N_OA + vir ) * ( NORB + 1 ) ];
            diagonal[ jump + vir + NVIR * act ] = - 2 * ( F_act + F_vir ) + Wtilde->get( irrep, irrep, NOCC + act, N_OA + vir, NOCC + act, N_OA + vir );
         }
      }
      jump += NACT * NVIR;

      for ( int occ = 0; occ < NOCC; occ++ ){
         const double F_occ = FMAT[ occ * ( NORB + 1 ) ];
         for ( int vir = 0; vir < NVIR; vir++ ){
            const double F_vir = FMAT[ ( N_OA + vir ) * ( NORB + 1 ) ];
            diagonal[ jump + vir + NVIR * occ ] = - 2 * ( F_occ + F_vir ) + Wtilde->get( irrep, irrep, occ, N_OA + vir, occ, N_OA + vir );
         }
      }
      jump += NOCC * NVIR;
   }

}

void CheMPS2::CASSCF::add_hessian( DMRGSCFmatrix * Fmatrix, DMRGSCFwtilde * Wtilde, const DMRGSCFindices * idx, double * origin, double * target ){

   const int n_irreps = idx->getNirreps();

   int jump = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NORB = idx->getNORB( irrep );
      const int NOCC = idx->getNOCC( irrep );
      const int NACT = idx->getNDMRG( irrep );
      const int NVIR = idx->getNVIRT( irrep );
      const int N_OA = NOCC + NACT;
      double * FMAT = Fmatrix->getBlock( irrep );
      const int ptr_AO = jump;
      const int ptr_VA = jump + NACT * NOCC;
      const int ptr_VO = jump + NACT * NOCC + NVIR * NACT;

      // OpenMP parallelism via dgemm_

      DGEMM_WRAP( +1.0, 'T', 'N', origin + ptr_VA,           FMAT + N_OA,        target + ptr_AO, NACT, NOCC, NVIR, NVIR, NORB, NACT );
      DGEMM_WRAP( +1.0, 'T', 'T', origin + ptr_VA,           FMAT + NORB * N_OA, target + ptr_AO, NACT, NOCC, NVIR, NVIR, NORB, NACT );
      DGEMM_WRAP( -1.0, 'N', 'N', origin + ptr_AO,           FMAT,               target + ptr_AO, NACT, NOCC, NOCC, NACT, NORB, NACT );
      DGEMM_WRAP( -1.0, 'N', 'T', origin + ptr_AO,           FMAT,               target + ptr_AO, NACT, NOCC, NOCC, NACT, NORB, NACT );
      DGEMM_WRAP( -1.0, 'N', 'N', FMAT + NOCC + NORB * NOCC, origin + ptr_AO,    target + ptr_AO, NACT, NOCC, NACT, NORB, NACT, NACT );
      DGEMM_WRAP( -1.0, 'T', 'N', FMAT + NOCC + NORB * NOCC, origin + ptr_AO,    target + ptr_AO, NACT, NOCC, NACT, NORB, NACT, NACT );
      DGEMM_WRAP( -1.0, 'N', 'N', FMAT + NOCC + NORB * N_OA, origin + ptr_VO,    target + ptr_AO, NACT, NOCC, NVIR, NORB, NVIR, NACT );
      DGEMM_WRAP( -1.0, 'T', 'N', FMAT + N_OA + NORB * NOCC, origin + ptr_VO,    target + ptr_AO, NACT, NOCC, NVIR, NORB, NVIR, NACT );

      DGEMM_WRAP( +1.0, 'N', 'T', FMAT + N_OA,               origin + ptr_AO,           target + ptr_VA, NVIR, NACT, NOCC, NORB, NACT, NVIR );
      DGEMM_WRAP( +1.0, 'T', 'T', FMAT + NORB * N_OA,        origin + ptr_AO,           target + ptr_VA, NVIR, NACT, NOCC, NORB, NACT, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', origin + ptr_VO,           FMAT + NORB * NOCC,        target + ptr_VA, NVIR, NACT, NOCC, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'T', origin + ptr_VO,           FMAT + NOCC,               target + ptr_VA, NVIR, NACT, NOCC, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', origin + ptr_VA,           FMAT + NOCC + NORB * NOCC, target + ptr_VA, NVIR, NACT, NACT, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'T', origin + ptr_VA,           FMAT + NOCC + NORB * NOCC, target + ptr_VA, NVIR, NACT, NACT, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', FMAT + N_OA + NORB * N_OA, origin + ptr_VA,           target + ptr_VA, NVIR, NACT, NVIR, NORB, NVIR, NVIR );
      DGEMM_WRAP( -1.0, 'T', 'N', FMAT + N_OA + NORB * N_OA, origin + ptr_VA,           target + ptr_VA, NVIR, NACT, NVIR, NORB, NVIR, NVIR );

      DGEMM_WRAP( -1.0, 'N', 'N', origin + ptr_VO,           FMAT,               target + ptr_VO, NVIR, NOCC, NOCC, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'T', origin + ptr_VO,           FMAT,               target + ptr_VO, NVIR, NOCC, NOCC, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', origin + ptr_VA,           FMAT + NOCC,        target + ptr_VO, NVIR, NOCC, NACT, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'T', origin + ptr_VA,           FMAT + NORB * NOCC, target + ptr_VO, NVIR, NOCC, NACT, NVIR, NORB, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', FMAT + N_OA + NORB * N_OA, origin + ptr_VO,    target + ptr_VO, NVIR, NOCC, NVIR, NORB, NVIR, NVIR );
      DGEMM_WRAP( -1.0, 'T', 'N', FMAT + N_OA + NORB * N_OA, origin + ptr_VO,    target + ptr_VO, NVIR, NOCC, NVIR, NORB, NVIR, NVIR );
      DGEMM_WRAP( -1.0, 'N', 'N', FMAT + N_OA + NORB * NOCC, origin + ptr_AO,    target + ptr_VO, NVIR, NOCC, NACT, NORB, NACT, NVIR );
      DGEMM_WRAP( -1.0, 'T', 'N', FMAT + NOCC + NORB * N_OA, origin + ptr_AO,    target + ptr_VO, NVIR, NOCC, NACT, NORB, NACT, NVIR );

      jump += NACT * NOCC + NVIR * NACT + NVIR * NOCC;
   }

   #pragma omp parallel
   {

      int jump_row = 0;
      for ( int irrep_row = 0; irrep_row < n_irreps; irrep_row++ ){

         const int NORB_row = idx->getNORB( irrep_row );
         const int NOCC_row = idx->getNOCC( irrep_row );
         const int NACT_row = idx->getNDMRG( irrep_row );
         const int NVIR_row = idx->getNVIRT( irrep_row );
         const int N_OA_row = NOCC_row + NACT_row;
         double * resAO = target + jump_row;
         double * resVA = resAO  + NACT_row * NOCC_row;
         double * resVO = resVA  + NVIR_row * NACT_row;

         int jump_col = 0;
         for ( int irrep_col = 0; irrep_col < n_irreps; irrep_col++ ){

            const int NOCC_col = idx->getNOCC( irrep_col );
            const int NACT_col = idx->getNDMRG( irrep_col );
            const int NVIR_col = idx->getNVIRT( irrep_col );
            const int N_OA_col = NOCC_col + NACT_col;
            double * vecAO = origin + jump_col;
            double * vecVA = vecAO  + NACT_col * NOCC_col;
            double * vecVO = vecVA  + NVIR_col * NACT_col;

            #pragma omp for schedule(static)
            for ( int act_row = 0; act_row < NACT_row; act_row++ ){
               for ( int occ_col = 0; occ_col < NOCC_col; occ_col++ ){
                  double * mat = Wtilde->getBlock( irrep_row, irrep_col, NOCC_row + act_row, occ_col );
                  DGEMV_WRAP( -1.0, mat +            NORB_row * NOCC_col, resAO + act_row,            vecAO + NACT_col * occ_col, NOCC_row, NACT_col, NORB_row, NACT_row, 1 );
                  DGEMV_WRAP( -1.0, mat +            NORB_row * N_OA_col, resAO + act_row,            vecVO + NVIR_col * occ_col, NOCC_row, NVIR_col, NORB_row, NACT_row, 1 );
                  DGEMV_WRAP(  1.0, mat + N_OA_row + NORB_row * NOCC_col, resVA + NVIR_row * act_row, vecAO + NACT_col * occ_col, NVIR_row, NACT_col, NORB_row, 1,        1 );
                  DGEMV_WRAP(  1.0, mat + N_OA_row + NORB_row * N_OA_col, resVA + NVIR_row * act_row, vecVO + NVIR_col * occ_col, NVIR_row, NVIR_col, NORB_row, 1,        1 );
               }
               for ( int act_col = 0; act_col < NACT_col; act_col++ ){
                  double * mat = Wtilde->getBlock( irrep_row, irrep_col, NOCC_row + act_row, NOCC_col + act_col );
                  DGEMV_WRAP(  1.0, mat,                                  resAO + act_row,            vecAO + act_col,            NOCC_row, NOCC_col, NORB_row, NACT_row, NACT_col );
                  DGEMV_WRAP( -1.0, mat            + NORB_row * N_OA_col, resAO + act_row,            vecVA + NVIR_col * act_col, NOCC_row, NVIR_col, NORB_row, NACT_row, 1        );
                  DGEMV_WRAP( -1.0, mat + N_OA_row,                       resVA + NVIR_row * act_row, vecAO + act_col,            NVIR_row, NOCC_col, NORB_row, 1,        NACT_col );
                  DGEMV_WRAP(  1.0, mat + N_OA_row + NORB_row * N_OA_col, resVA + NVIR_row * act_row, vecVA + NVIR_col * act_col, NVIR_row, NVIR_col, NORB_row, 1,        1        );
               }
            }

            #pragma omp for schedule(static)
            for ( int occ_row = 0; occ_row < NOCC_row; occ_row++ ){
               for ( int occ_col = 0; occ_col < NOCC_col; occ_col++ ){
                  double * mat = Wtilde->getBlock( irrep_row, irrep_col, occ_row, occ_col );
                  DGEMV_WRAP( 1.0, mat + NOCC_row + NORB_row * NOCC_col, resAO + NACT_row * occ_row, vecAO + NACT_col * occ_col, NACT_row, NACT_col, NORB_row, 1, 1 );
                  DGEMV_WRAP( 1.0, mat + NOCC_row + NORB_row * N_OA_col, resAO + NACT_row * occ_row, vecVO + NVIR_col * occ_col, NACT_row, NVIR_col, NORB_row, 1, 1 );
                  DGEMV_WRAP( 1.0, mat + N_OA_row + NORB_row * NOCC_col, resVO + NVIR_row * occ_row, vecAO + NACT_col * occ_col, NVIR_row, NACT_col, NORB_row, 1, 1 );
                  DGEMV_WRAP( 1.0, mat + N_OA_row + NORB_row * N_OA_col, resVO + NVIR_row * occ_row, vecVO + NVIR_col * occ_col, NVIR_row, NVIR_col, NORB_row, 1, 1 );
               }
               for ( int act_col = 0; act_col < NACT_col; act_col++ ){
                  double * mat = Wtilde->getBlock( irrep_row, irrep_col, occ_row, NOCC_col + act_col );
                  DGEMV_WRAP( -1.0, mat + NOCC_row,                       resAO + NACT_row * occ_row, vecAO + act_col,            NACT_row, NOCC_col, NORB_row, 1, NACT_col );
                  DGEMV_WRAP(  1.0, mat + NOCC_row + NORB_row * N_OA_col, resAO + NACT_row * occ_row, vecVA + NVIR_col * act_col, NACT_row, NVIR_col, NORB_row, 1, 1        );
                  DGEMV_WRAP( -1.0, mat + N_OA_row,                       resVO + NVIR_row * occ_row, vecAO + act_col,            NVIR_row, NOCC_col, NORB_row, 1, NACT_col );
                  DGEMV_WRAP(  1.0, mat + N_OA_row + NORB_row * N_OA_col, resVO + NVIR_row * occ_row, vecVA + NVIR_col * act_col, NVIR_row, NVIR_col, NORB_row, 1, 1        );
               }
            }
            jump_col += NACT_col * NOCC_col + NVIR_col * NACT_col + NVIR_col * NOCC_col;
         }
         jump_row += NACT_row * NOCC_row + NVIR_row * NACT_row + NVIR_row * NOCC_row;
      }
   }

}

void CheMPS2::CASSCF::DGEMM_WRAP( double prefactor, char transA, char transB, double * A, double * B, double * C, int m, int n, int k, int lda, int ldb, int ldc ){

   if ( m * k * n == 0 ){ return; }
   double add = 1.0;
   dgemm_( &transA, &transB, &m, &n, &k, &prefactor, A, &lda, B, &ldb, &add, C, &ldc );

}

void CheMPS2::CASSCF::DGEMV_WRAP( double prefactor, double * matrix, double * result, double * vector, int rowdim, int coldim, int ldmat, int incres, int incvec ){

   if ( rowdim * coldim == 0 ){ return; }
   char notrans = 'N';
   double add = 1.0;
   dgemv_( &notrans, &rowdim, &coldim, &prefactor, matrix, &ldmat, vector, &incvec, &add, result, &incres );

}

void CheMPS2::CASSCF::buildWtilde(DMRGSCFwtilde * localwtilde, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM){

   localwtilde->clear();
   const int numIrreps  = localIdx->getNirreps();
   const int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   for (int irrep_pq = 0; irrep_pq < numIrreps; irrep_pq++){
   
      const int NumOCCpq  = localIdx->getNOCC(  irrep_pq );
      const int NumDMRGpq = localIdx->getNDMRG( irrep_pq );
      const int NumORBpq  = localIdx->getNORB(  irrep_pq );
      const int NumCOREpq = NumOCCpq + NumDMRGpq;
      
      //If irrep_pq == irrep_rs and P == R occupied --> QS only active or virtual
      #pragma omp parallel for schedule(static)
      for (int relindexP = 0; relindexP < NumOCCpq; relindexP++){
         double * subblock = localwtilde->getBlock( irrep_pq, irrep_pq, relindexP, relindexP);
         for (int relindexS = NumOCCpq; relindexS < NumORBpq; relindexS++){
            for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 4 * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                   + localJKocc->get(irrep_pq, relindexQ, relindexS)
                                                                   + localJKact->get(irrep_pq, relindexQ, relindexS) );
            }
         }
      }
      
      //If irrep_pq == irrep_rs and P,R active --> QS only occupied or virtual
      #pragma omp parallel for schedule(static)
      for (int combined = 0; combined < NumDMRGpq*NumDMRGpq; combined++){
         
         const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
         const int relindexR  = NumOCCpq + ( combined / NumDMRGpq );
         double * subblock    = localwtilde->getBlock( irrep_pq, irrep_pq, relindexP, relindexR );
         const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
         const int DMRGindexR = relindexR - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
         const double OneDMvalue = local1DM[ DMRGindexP + totOrbDMRG * DMRGindexR ];
         
         for (int relindexS = 0; relindexS < NumOCCpq; relindexS++){
            for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
            for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
         }
         for (int relindexS = NumCOREpq; relindexS < NumORBpq; relindexS++){
            for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
            for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
         }
      }
   
      for (int irrep_rs = 0; irrep_rs < numIrreps; irrep_rs++){
      
         const int NumOCCrs  = localIdx->getNOCC(  irrep_rs );
         const int NumDMRGrs = localIdx->getNDMRG( irrep_rs );
         const int NumORBrs  = localIdx->getNORB(  irrep_rs );
         const int NumCORErs = NumOCCrs + NumDMRGrs;
         const int productirrep = Irreps::directProd( irrep_pq, irrep_rs );
      
         // P and R occupied --> QS only active or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumOCCpq*NumOCCrs; combined++){
         
            const int relindexP = combined % NumOCCpq;
            const int relindexR = combined / NumOCCpq;
            double * subblock   = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR);
            
            for (int relindexS = NumOCCrs; relindexS < NumORBrs; relindexS++){
               for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
                  subblock[ relindexQ + NumORBpq * relindexS ] +=  
                      4 * ( 4 * theInts->FourIndexAPI(irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relindexR)
                              - theInts->FourIndexAPI(irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relindexR)
                              - theInts->FourIndexAPI(irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relindexP) );
               }
            }
         } // End combined P and R occupied
         
         // P and R active --> QS only occupied or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumDMRGpq*NumDMRGrs; combined++){
         
            const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
            const int relindexR  = NumOCCrs + ( combined / NumDMRGpq );
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
            const int DMRGindexR = relindexR - NumOCCrs + localIdx->getDMRGcumulative( irrep_rs );
            
            for (int irrep_alpha = 0; irrep_alpha < numIrreps; irrep_alpha++){
               
               const int irrep_beta   = Irreps::directProd( irrep_alpha, productirrep );
               const int NumDMRGalpha = localIdx->getNDMRG( irrep_alpha );
               const int NumDMRGbeta  = localIdx->getNDMRG( irrep_beta  );
               
               for (int alpha = 0; alpha < NumDMRGalpha; alpha++){
               
                  const int DMRGalpha = localIdx->getDMRGcumulative( irrep_alpha ) + alpha;
                  const int relalpha  = localIdx->getNOCC( irrep_alpha ) + alpha;
                  
                  for (int beta = 0; beta < NumDMRGbeta; beta++){
                  
                     const int DMRGbeta = localIdx->getDMRGcumulative( irrep_beta ) + beta;
                     const int relbeta  = localIdx->getNOCC( irrep_beta ) + beta;
                     
                     const double TwoDMvalue1  = local2DM[ DMRGindexR + totOrbDMRG * ( DMRGalpha + totOrbDMRG * ( DMRGindexP + totOrbDMRG * DMRGbeta ) ) ];
                     const double TwoDMvalue23 = local2DM[ DMRGindexR + totOrbDMRG * ( DMRGalpha + totOrbDMRG * ( DMRGbeta + totOrbDMRG * DMRGindexP ) ) ]
                                               + local2DM[ DMRGindexR + totOrbDMRG * ( DMRGindexP + totOrbDMRG * ( DMRGbeta + totOrbDMRG * DMRGalpha ) ) ];
                     
                     for (int relindexS = 0; relindexS < NumOCCrs; relindexS++){
                        for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                        for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                     }
                     for (int relindexS = NumCORErs; relindexS < NumORBrs; relindexS++){
                        for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                        for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                     }
                  }
               }
            }
         } // End combined P and R active
         
         // P active and R occupied  -->  Q occupied or virtual  //  S active or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumDMRGpq*NumOCCrs; combined++){
            
            const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
            const int relindexR  = combined / NumDMRGpq;
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
            
            for (int alpha = 0; alpha < NumDMRGpq; alpha++){
            
               const int DMRGalpha     = localIdx->getDMRGcumulative( irrep_pq ) + alpha;
               const int relalpha      = localIdx->getNOCC( irrep_pq ) + alpha;
               const double OneDMvalue = local1DM[ DMRGalpha + totOrbDMRG * DMRGindexP ];
               
               for (int relindexS = NumOCCrs; relindexS < NumORBrs; relindexS++){
                  for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relalpha, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relalpha, relindexS, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relalpha) );
                  }
                  for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relalpha, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relalpha, relindexS, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relalpha) );
                  }
               }
            }
         } // End combined P active and R occupied
         
         // P occupied and R active  -->  Q active or virtual  //  S occupied or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumOCCpq*NumDMRGrs; combined++){
            
            const int relindexP  = combined % NumOCCpq;
            const int relindexR  = NumOCCrs + ( combined / NumOCCpq );
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexR = relindexR - NumOCCrs + localIdx->getDMRGcumulative( irrep_rs );
            
            for (int beta = 0; beta < NumDMRGrs; beta++){
            
               const int DMRGbeta      = localIdx->getDMRGcumulative( irrep_rs ) + beta;
               const int relbeta       = localIdx->getNOCC( irrep_rs ) + beta;
               const double OneDMvalue = local1DM[ DMRGindexR + totOrbDMRG * DMRGbeta ];
               
               for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
                  for (int relindexS = 0; relindexS < NumOCCrs; relindexS++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relbeta, relindexP) );
                  }
                  for (int relindexS = NumCORErs; relindexS < NumORBrs; relindexS++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relbeta, relindexP) );
                  }
               }
            }
         } // End combined P occupied and R active
         
      }
   }

}

void CheMPS2::CASSCF::buildFmat(DMRGSCFmatrix * localFmat, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM){
   
   localFmat->clear();
   const int numIrreps  = localIdx->getNirreps();
   const int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   for (int irrep_pq = 0; irrep_pq < numIrreps; irrep_pq++){
   
      const int NumORB  = localIdx->getNORB(  irrep_pq );
      const int NumOCC  = localIdx->getNOCC(  irrep_pq );
      const int NumDMRG = localIdx->getNDMRG( irrep_pq );
      const int NumOCCDMRG = NumOCC + NumDMRG;
      
      #pragma omp parallel for schedule(static)
      for (int p = 0; p < NumOCC; p++){
         for (int q = 0; q < NumORB; q++){
            localFmat->set( irrep_pq, p, q, 2 * ( localTmat->get(  irrep_pq, q, p )
                                                + localJKocc->get( irrep_pq, q, p )
                                                + localJKact->get( irrep_pq, q, p ) ) );
         }
      }
      
      #pragma omp parallel for schedule(static)
      for (int p = NumOCC; p < NumOCCDMRG; p++){
         const int DMRGindex_p = p - NumOCC + localIdx->getDMRGcumulative( irrep_pq );
         
         //One-body terms --> matrix multiplication?
         for (int r = NumOCC; r < NumOCCDMRG; r++){
            const double OneDMvalue = local1DM[ DMRGindex_p + totOrbDMRG * ( DMRGindex_p + r - p ) ];
            for (int q = 0; q < NumORB; q++){
               localFmat->getBlock(irrep_pq)[ p + NumORB * q ] += OneDMvalue * ( localTmat->get( irrep_pq, q, r ) + localJKocc->get( irrep_pq, q, r) );
            }
         }
         
         //Two-body terms --> matrix multiplication possible?
         for (int irrep_r = 0; irrep_r < numIrreps; irrep_r++){
            const int irrep_product = Irreps::directProd(irrep_pq, irrep_r);
            for (int irrep_s = 0; irrep_s < numIrreps; irrep_s++){
               const int irrep_t = Irreps::directProd(irrep_product, irrep_s);
               for (int r = localIdx->getNOCC(irrep_r); r < localIdx->getNOCC(irrep_r) + localIdx->getNDMRG(irrep_r); r++){
                  const int DMRGindex_r = r - localIdx->getNOCC( irrep_r ) + localIdx->getDMRGcumulative( irrep_r );
                  for (int s = localIdx->getNOCC(irrep_s); s < localIdx->getNOCC(irrep_s) + localIdx->getNDMRG(irrep_s); s++){
                     const int DMRGindex_s = s - localIdx->getNOCC( irrep_s ) + localIdx->getDMRGcumulative( irrep_s );
                     for (int t = localIdx->getNOCC(irrep_t); t < localIdx->getNOCC(irrep_t) + localIdx->getNDMRG(irrep_t); t++){
                        const int DMRGindex_t = t - localIdx->getNOCC( irrep_t ) + localIdx->getDMRGcumulative( irrep_t );
                        const double TwoDMvalue = local2DM[ DMRGindex_p + totOrbDMRG * ( DMRGindex_r + totOrbDMRG * ( DMRGindex_s + totOrbDMRG * DMRGindex_t ) ) ];
                        for (int q = 0; q < NumORB; q++){
                           localFmat->getBlock(irrep_pq)[ p + NumORB * q ] += TwoDMvalue * theInts->FourIndexAPI(irrep_pq, irrep_r, irrep_s, irrep_t, q, r, s, t);
                        }
                     }
                  }
               }
            }
         }
      }
   }

}


