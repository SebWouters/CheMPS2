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
#include "DMRG.h"
#include "FCI.h"
#include "DMRGSCFrotations.h"
#include "Lapack.h"
#include "Cumulant.h"
#include "CASPT2.h"
#include "MPIchemps2.h"
#include "EdmistonRuedenberg.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

void CheMPS2::CASSCF::write_f4rdm_checkpoint( const string f4rdm_file, int * hamorb1, int * hamorb2, const int tot_dmrg_power6, double * contract ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){

      hid_t file_id  = H5Fcreate( f4rdm_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
      hid_t group_id = H5Gcreate( file_id, "/F4RDM", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

      hsize_t dimarray1   = 1;
      hid_t dataspace1_id = H5Screate_simple( 1, &dimarray1, NULL );
      hid_t dataset1_id   = H5Dcreate( group_id, "hamorb1", H5T_NATIVE_INT, dataspace1_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      H5Dwrite( dataset1_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hamorb1 );
      H5Dclose( dataset1_id );
      H5Sclose( dataspace1_id );

      hsize_t dimarray2   = 1;
      hid_t dataspace2_id = H5Screate_simple( 1, &dimarray2, NULL );
      hid_t dataset2_id   = H5Dcreate( group_id, "hamorb2", H5T_NATIVE_INT, dataspace2_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      H5Dwrite( dataset2_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hamorb2 );
      H5Dclose( dataset2_id );
      H5Sclose( dataspace2_id );

      hsize_t dimarray3   = tot_dmrg_power6;
      hid_t dataspace3_id = H5Screate_simple( 1, &dimarray3, NULL );
      hid_t dataset3_id   = H5Dcreate( group_id, "contract", H5T_NATIVE_DOUBLE, dataspace3_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      H5Dwrite( dataset3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, contract );
      H5Dclose( dataset3_id );
      H5Sclose( dataspace3_id );

      H5Gclose( group_id );
      H5Fclose( file_id );

      cout << "Created F.4-RDM checkpoint file " << f4rdm_file << " at next orbitals ( " << hamorb1[ 0 ] << " , " << hamorb2[ 0 ] << " )." << endl;

   }

}

bool CheMPS2::CASSCF::read_f4rdm_checkpoint( const string f4rdm_file, int * hamorb1, int * hamorb2, const int tot_dmrg_power6, double * contract ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   // Check whether the file exists
   int exists = 0;
   if ( am_i_master ){
      struct stat file_info;
      const int file_stat = stat( f4rdm_file.c_str(), &file_info );
      if ( file_stat == 0 ){ exists = 1; }
   }
   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_array_int( &exists, 1, MPI_CHEMPS2_MASTER );
   #endif
   if ( exists == 0 ){
      return false; // Not loaded
   }

   if ( am_i_master ){

      hid_t file_id  = H5Fopen( f4rdm_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
      hid_t group_id = H5Gopen( file_id, "/F4RDM", H5P_DEFAULT );

      hid_t dataset1_id = H5Dopen( group_id, "hamorb1", H5P_DEFAULT );
      H5Dread( dataset1_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hamorb1 );
      H5Dclose( dataset1_id );

      hid_t dataset2_id = H5Dopen( group_id, "hamorb2", H5P_DEFAULT );
      H5Dread( dataset2_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, hamorb2 );
      H5Dclose( dataset2_id );

      hid_t dataset3_id = H5Dopen( group_id, "contract", H5P_DEFAULT );
      H5Dread( dataset3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, contract );
      H5Dclose( dataset3_id );

      H5Gclose( group_id );
      H5Fclose( file_id );

   }
   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_array_int( hamorb1, 1, MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_array_int( hamorb2, 1, MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_array_double( contract, tot_dmrg_power6, MPI_CHEMPS2_MASTER );
   #endif

   return true; // Loaded

}

void CheMPS2::CASSCF::fock_dot_4rdm( double * fockmx, CheMPS2::DMRG * dmrgsolver, CheMPS2::Hamiltonian * ham, int next_orb1, int next_orb2, double * work, double * result, const bool CHECKPOINT, const bool PSEUDOCANONICAL ){

   const int LAS = ham->getL();
   int size      = LAS * LAS * LAS * LAS * LAS * LAS;
   int inc1      = 1;

   for ( int diag = 0; diag < LAS; diag++ ){
      if (( next_orb1 == diag ) && ( next_orb2 == diag )){
         double prefactor = 0.5 * fockmx[ diag + LAS * diag ];
         if ( fabs( prefactor ) > 0.0 ){
            dmrgsolver->Symm4RDM( work, diag, diag, false );
            daxpy_( &size, &prefactor, work, &inc1, result, &inc1 );
         }
         if ( diag == LAS - 1 ){
            next_orb1 = 0;
            next_orb2 = 1;
         } else {
            next_orb1 = diag + 1;
            next_orb2 = diag + 1;
         }
         if ( CHECKPOINT ){ write_f4rdm_checkpoint( CheMPS2::DMRGSCF_f4rdm_name, &next_orb1, &next_orb2, size, result ); }
      }
   }

   if ( PSEUDOCANONICAL == false ){
      for ( int orb1 = 0; orb1 < LAS; orb1++ ){
         for ( int orb2 = orb1 + 1; orb2 < LAS; orb2++ ){
            if (( next_orb1 == orb1 ) && ( next_orb2 == orb2 )){
               double prefactor = 0.5 * ( fockmx[ orb1 + LAS * orb2 ] + fockmx[ orb2 + LAS * orb1 ] );
               if (( ham->getOrbitalIrrep( orb1 ) == ham->getOrbitalIrrep( orb2 ) ) && ( fabs( prefactor ) > 0.0 )){
                  dmrgsolver->Symm4RDM( work, orb1, orb2, false );
                  daxpy_( &size, &prefactor, work, &inc1, result, &inc1 );
               }
               if ( orb2 == LAS - 1 ){
                  next_orb1 = next_orb1 + 1;
                  next_orb2 = next_orb1 + 1;
               } else {
                  next_orb2 = next_orb2 + 1;
               }
               if (( ham->getOrbitalIrrep( orb1 ) == ham->getOrbitalIrrep( orb2 ) ) && ( CHECKPOINT )){
                  write_f4rdm_checkpoint( CheMPS2::DMRGSCF_f4rdm_name, &next_orb1, &next_orb2, size, result );
               }
            }
         }
      }
   }

}

double CheMPS2::CASSCF::caspt2( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * scf_options, const double IPEA, const double IMAG, const bool PSEUDOCANONICAL, const bool CHECKPOINT, const bool CUMULANT ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   const int num_elec = Nelectrons - 2 * iHandler->getNOCCsum();
   assert( num_elec >= 0 );

   if ( CASPT2::vector_length( iHandler ) == 0 ){
      if ( am_i_master ){
         cout << "CheMPS2::CASSCF::caspt2 : There are no CASPT2 excitations between the CORE, ACTIVE, and VIRTUAL orbital spaces." << endl;
      }
      return 0.0;
   }

   //Determine the maximum NORB(irrep) and the max_block_size for the ERI orbital rotation
   const int maxlinsize      = iHandler->getNORBmax();
   const long long fullsize  = ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize );
   const string tmp_filename = tmp_folder + "/" + CheMPS2::DMRGSCF_eri_storage_name;
   const int dmrgsize_power4 = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
   //For (ERI rotation, update unitary, block diagonalize, orbital localization)
   DMRGSCFintegrals * theRotatedTEI = new DMRGSCFintegrals( iHandler );
   const int temp_work_size = (( fullsize > CheMPS2::DMRGSCF_max_mem_eri_tfo ) ? CheMPS2::DMRGSCF_max_mem_eri_tfo : fullsize );
   const int work_mem_size  = max( max( temp_work_size , maxlinsize * maxlinsize * 4 ) , dmrgsize_power4 );
   const int tot_dmrg_power6 = dmrgsize_power4 * nOrbDMRG * nOrbDMRG;
   double * mem1 = new double[ work_mem_size ];
   double * mem2 = new double[ ( PSEUDOCANONICAL ) ? work_mem_size : max( work_mem_size, tot_dmrg_power6 ) ];

   // If you did not run CheMPS2::CASSCF::solve, you NEED to load the unitary from disk
   if ( successful_solve == false ){
      assert( scf_options->getStoreUnitary() );
      if ( am_i_master ){
         struct stat file_info;
         const int file_stat = stat( (scf_options->getUnitaryStorageName()).c_str(), &file_info );
         assert( file_stat == 0 );
         unitary->loadU( scf_options->getUnitaryStorageName() );
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      unitary->broadcast( MPI_CHEMPS2_MASTER );
      #endif
   } else {
      // Rotate to pseudocanonical orbitals
      if ( PSEUDOCANONICAL ){ // successful_solve needs to be true, because the 1-RDM is needed for buildQmatACT();
         buildTmatrix();
         buildQmatOCC();
         buildQmatACT();
         construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
         block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
         block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
         block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
         if (( am_i_master ) && ( scf_options->getStoreUnitary() )){ unitary->saveU( scf_options->getUnitaryStorageName() ); }
      }
   }

   // Fill active space Hamiltonian
   Hamiltonian * HamAS = new Hamiltonian( nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital() );
   Problem * Prob = new Problem( HamAS, TwoS, num_elec, Irrep );
   Prob->SetupReorderD2h(); // Doesn't matter if the group isn't D2h, Prob checks it.
   buildTmatrix();
   buildQmatOCC();
   fillConstAndTmatDMRG( HamAS );
   if ( am_i_master ){
      DMRGSCFrotations::rotate( VMAT_ORIG, HamAS->getVmat(), NULL, 'A', 'A', 'A', 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
   }
   #ifdef CHEMPS2_MPI_COMPILATION
   HamAS->getVmat()->broadcast( MPI_CHEMPS2_MASTER );
   #endif

   // Reorder the orbitals based on the Fiedler vector of the exchange matrix
   if ( scf_options->getWhichActiveSpace() == 3 ){
      int * dmrg2ham = new int[ nOrbDMRG ];
      if ( am_i_master ){
         EdmistonRuedenberg * theLocalizer = new EdmistonRuedenberg( HamAS->getVmat(), iHandler->getGroupNumber() );
         theLocalizer->FiedlerGlobal( dmrg2ham );
         delete theLocalizer;
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      MPIchemps2::broadcast_array_int( dmrg2ham, nOrbDMRG, MPI_CHEMPS2_MASTER );
      #endif
      Prob->setup_reorder_custom( dmrg2ham );
      delete [] dmrg2ham;
   }

   double E_CASSCF = 0.0;
   double * three_dm = new double[ tot_dmrg_power6 ];
   double * contract = new double[ tot_dmrg_power6 ];
   for ( int cnt = 0; cnt < tot_dmrg_power6; cnt++ ){ contract[ cnt ] = 0.0; }

   int next_hamorb1 = 0;
   int next_hamorb2 = 0;
   const bool make_checkpt = (( CUMULANT == false ) && ( CHECKPOINT ));
   bool checkpt_loaded = false;
   if ( make_checkpt ){
      assert(( OptScheme != NULL ) || ( rootNum > 1 ));
      checkpt_loaded = read_f4rdm_checkpoint( CheMPS2::DMRGSCF_f4rdm_name, &next_hamorb1, &next_hamorb2, tot_dmrg_power6, contract );
   }

   // Solve the active space problem
   if (( OptScheme == NULL ) && ( rootNum == 1 )){ // Do FCI

      if ( am_i_master ){
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
      }
      #ifdef CHEMPS2_MPI_COMPILATION
      MPIchemps2::broadcast_array_double( &E_CASSCF, 1, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_double(  DMRG2DM, dmrgsize_power4, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_double( three_dm, tot_dmrg_power6, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_double( contract, tot_dmrg_power6, MPI_CHEMPS2_MASTER );
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );
      #endif

   } else { // Do the DMRG sweeps

      assert( OptScheme != NULL );
      for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] = 0.0; } // Clear the 2-RDM
      CheMPS2::DMRG * theDMRG = new DMRG( Prob, OptScheme, make_checkpt, tmp_folder );
      for ( int state = 0; state < rootNum; state++ ){
         if ( state > 0 ){ theDMRG->newExcitation( fabs( E_CASSCF ) ); }
         if ( checkpt_loaded == false ){ E_CASSCF = theDMRG->Solve(); }
         if (( state == 0 ) && ( rootNum > 1 )){ theDMRG->activateExcitations( rootNum - 1 ); }
      }
      theDMRG->calc_rdms_and_correlations( true );
      copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM ); // 2-RDM
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM ); // 1-RDM
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler );
      copy_active( theFmatrix, mem2, iHandler ); // Fock
      if ( CUMULANT ){
         CheMPS2::Cumulant::gamma4_fock_contract_ham( Prob, theDMRG->get3DM(), theDMRG->get2DM(), mem2, contract );
      } else {
         fock_dot_4rdm( mem2, theDMRG, HamAS, next_hamorb1, next_hamorb2, three_dm, contract, make_checkpt, PSEUDOCANONICAL );
      }
      theDMRG->get3DM()->fill_ham_index( 1.0, false, three_dm, 0, nOrbDMRG );
      if (( CheMPS2::DMRG_storeMpsOnDisk ) && ( make_checkpt == false )){ theDMRG->deleteStoredMPS(); }
      if ( CheMPS2::DMRG_storeRenormOptrOnDisk ){ theDMRG->deleteStoredOperators(); }
      delete theDMRG;

   }

   delete Prob;
   delete HamAS;

   if ( PSEUDOCANONICAL == false ){
      if ( am_i_master ){ cout << "CASPT2 : Deviation from pseudocanonical = " << deviation_from_blockdiag( theFmatrix, iHandler ) << endl; }
      block_diagonalize( 'O', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
      block_diagonalize( 'A', theFmatrix, unitary, mem1, mem2, iHandler, false, DMRG2DM, three_dm, contract ); // 2-RDM, 3-RDM, and trace( Fock * cu(4)-4-RDM )
      block_diagonalize( 'V', theFmatrix, unitary, mem1, mem2, iHandler, false, NULL, NULL, NULL );
      setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM ); // 1-RDM
      buildTmatrix();
      buildQmatOCC();
      buildQmatACT();
      construct_fock( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler ); // Fock
   }

   // Calculate the matrix elements needed to calculate the CASPT2 V-vector
   if ( am_i_master ){
      DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'C', 'F', 'F', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
      DMRGSCFrotations::rotate( VMAT_ORIG, NULL, theRotatedTEI, 'C', 'V', 'C', 'V', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename );
      delete_file( tmp_filename );
   }

   delete [] mem1;
   delete [] mem2;

   double E_CASPT2 = 0.0;
   if ( am_i_master ){
      cout << "CASPT2 : Deviation from pseudocanonical = " << deviation_from_blockdiag( theFmatrix, iHandler ) << endl;
      CheMPS2::CASPT2 * myCASPT2 = new CheMPS2::CASPT2( iHandler, theRotatedTEI, theTmatrix, theFmatrix, DMRG1DM, DMRG2DM, three_dm, contract, IPEA );
      delete theRotatedTEI;
      delete [] three_dm;
      delete [] contract;
      E_CASPT2 = myCASPT2->solve( IMAG );
      delete myCASPT2;
   } else {
      delete theRotatedTEI;
      delete [] three_dm;
      delete [] contract;
   }
   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_array_double( &E_CASPT2, 1, MPI_CHEMPS2_MASTER );
   #endif

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


