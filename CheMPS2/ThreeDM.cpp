/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <sstream>

#include <hdf5.h>

#include "ThreeDM.h"
#include "Lapack.h"
#include "Options.h"
#include "MPIchemps2.h"
#include "Wigner.h"
#include "Special.h"

using std::max;

CheMPS2::ThreeDM::ThreeDM( const SyBookkeeper * book_in, const Problem * prob_in, const bool disk_in ){

   book = book_in;
   prob = prob_in;
   disk = disk_in;

   L = book->gL();
   {
      const long long linsize = ( long long ) L;
      const long long size    = (( disk ) ? linsize * linsize * linsize * linsize * linsize
                                          : linsize * linsize * linsize * linsize * linsize * linsize );
      assert( INT_MAX >= size );
      array_size = size;
   }

   elements = new double[ array_size ];
   #pragma omp simd
   for ( int cnt = 0; cnt < array_size; cnt++ ){ elements[ cnt ] = 0.0; }

   if ( disk ){
      temp_disk_orbs = new int[ 6 * array_size ];
      temp_disk_vals = new double [ array_size ];
      create_file();
   } else {
      temp_disk_orbs = NULL;
      temp_disk_vals = NULL;
   }

}

CheMPS2::ThreeDM::~ThreeDM(){

   delete [] elements;
   if ( disk ){ delete [] temp_disk_orbs;
                delete [] temp_disk_vals; }

}

#ifdef CHEMPS2_MPI_COMPILATION
void CheMPS2::ThreeDM::mpi_allreduce(){

   if ( disk ){

      for ( int orb = 0; orb < L; orb++ ){
         read_file( orb );
         MPIchemps2::allreduce_array_double( elements, temp_disk_vals, array_size );
         #pragma omp simd
         for ( int cnt = 0; cnt < array_size; cnt++ ){ elements[ cnt ] = temp_disk_vals[ cnt ]; }
         write_file( orb );
      }

   } else {

      double * temp = new double[ array_size ];
      MPIchemps2::allreduce_array_double( elements, temp, array_size );
      #pragma omp simd
      for ( int cnt = 0; cnt < array_size; cnt++ ){ elements[ cnt ] = temp[ cnt ]; }
      delete [] temp;

   }

}
#endif

void CheMPS2::ThreeDM::set_dmrg_index( const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6, const double value ){

   // prob assumes you use DMRG orbs, while elements is stored in hamiltonian orbitals
   // irrep sanity checks are performed in ThreeDM::fill_site

   const int orb1 = (( prob->gReorder() ) ? prob->gf2( cnt1 ) : cnt1 );
   const int orb2 = (( prob->gReorder() ) ? prob->gf2( cnt2 ) : cnt2 );
   const int orb3 = (( prob->gReorder() ) ? prob->gf2( cnt3 ) : cnt3 );
   const int orb4 = (( prob->gReorder() ) ? prob->gf2( cnt4 ) : cnt4 );
   const int orb5 = (( prob->gReorder() ) ? prob->gf2( cnt5 ) : cnt5 );
   const int orb6 = (( prob->gReorder() ) ? prob->gf2( cnt6 ) : cnt6 );

   if ( disk ){
      int private_counter = -1;
      #pragma omp critical
      {
         private_counter = temp_disk_counter;
         temp_disk_counter++;
      }
      assert( private_counter < array_size );
      temp_disk_orbs[ 6 * private_counter + 0 ] = orb1;
      temp_disk_orbs[ 6 * private_counter + 1 ] = orb2;
      temp_disk_orbs[ 6 * private_counter + 2 ] = orb3;
      temp_disk_orbs[ 6 * private_counter + 3 ] = orb4;
      temp_disk_orbs[ 6 * private_counter + 4 ] = orb5;
      temp_disk_orbs[ 6 * private_counter + 5 ] = orb6;
      temp_disk_vals[ private_counter ] = value;
      return;
   }

   elements[ orb1 + L * ( orb2 + L * ( orb3 + L * ( orb4 + L * ( orb5 + L * orb6 )))) ] = value;
   elements[ orb2 + L * ( orb3 + L * ( orb1 + L * ( orb5 + L * ( orb6 + L * orb4 )))) ] = value;
   elements[ orb3 + L * ( orb1 + L * ( orb2 + L * ( orb6 + L * ( orb4 + L * orb5 )))) ] = value;
   elements[ orb2 + L * ( orb1 + L * ( orb3 + L * ( orb5 + L * ( orb4 + L * orb6 )))) ] = value;
   elements[ orb3 + L * ( orb2 + L * ( orb1 + L * ( orb6 + L * ( orb5 + L * orb4 )))) ] = value;
   elements[ orb1 + L * ( orb3 + L * ( orb2 + L * ( orb4 + L * ( orb6 + L * orb5 )))) ] = value;

   elements[ orb4 + L * ( orb5 + L * ( orb6 + L * ( orb1 + L * ( orb2 + L * orb3 )))) ] = value;
   elements[ orb5 + L * ( orb6 + L * ( orb4 + L * ( orb2 + L * ( orb3 + L * orb1 )))) ] = value;
   elements[ orb6 + L * ( orb4 + L * ( orb5 + L * ( orb3 + L * ( orb1 + L * orb2 )))) ] = value;
   elements[ orb5 + L * ( orb4 + L * ( orb6 + L * ( orb2 + L * ( orb1 + L * orb3 )))) ] = value;
   elements[ orb6 + L * ( orb5 + L * ( orb4 + L * ( orb3 + L * ( orb2 + L * orb1 )))) ] = value;
   elements[ orb4 + L * ( orb6 + L * ( orb5 + L * ( orb1 + L * ( orb3 + L * orb2 )))) ] = value;

}

double CheMPS2::ThreeDM::get_ham_index( const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6 ) const{

   assert( disk == false );
   return elements[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * ( cnt4 + L * ( cnt5 + L * cnt6 )))) ];

}

void CheMPS2::ThreeDM::fill_ham_index( const double alpha, const bool add, double * storage, const int last_orb_start, const int last_orb_num ){

   assert( last_orb_start >= 0 );
   assert( last_orb_num   >= 1 );
   assert( last_orb_start + last_orb_num <= L );

   if ( disk ){

      for ( int ham_orb = last_orb_start; ham_orb < ( last_orb_start + last_orb_num ); ham_orb++ ){
         read_file( ham_orb );
         const int shift = ( ham_orb - last_orb_start ) * array_size;
         if ( add == false ){
            #pragma omp simd
            for ( int cnt = 0; cnt < array_size; cnt++ ){ storage[ shift + cnt ]  = alpha * elements[ cnt ]; }
         } else {
            #pragma omp simd
            for ( int cnt = 0; cnt < array_size; cnt++ ){ storage[ shift + cnt ] += alpha * elements[ cnt ]; }
         }
      }

   } else {

      const int shift = last_orb_start * L * L * L * L * L;
      const int size  = last_orb_num   * L * L * L * L * L;
      if ( add == false ){
         #pragma omp simd
         for ( int cnt = 0; cnt < size; cnt++ ){ storage[ cnt ]  = alpha * elements[ shift + cnt ]; }
      } else {
         #pragma omp simd
         for ( int cnt = 0; cnt < size; cnt++ ){ storage[ cnt ] += alpha * elements[ shift + cnt ]; }
      }

   }

}

double CheMPS2::ThreeDM::trace(){

   double value = 0.0;

   if ( disk ){

      for ( int cnt3 = 0; cnt3 < L; cnt3++ ){
         read_file( cnt3 );
         for ( int cnt2 = 0; cnt2 < L; cnt2++ ){
            for ( int cnt1 = 0; cnt1 < L; cnt1++ ){
               value += elements[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * ( cnt1 + L * cnt2 ))) ];
            }
         }
      }

   } else {

      for ( int cnt3 = 0; cnt3 < L; cnt3++ ){
         for ( int cnt2 = 0; cnt2 < L; cnt2++ ){
            for ( int cnt1 = 0; cnt1 < L; cnt1++ ){
               value += get_ham_index( cnt1, cnt2, cnt3, cnt1, cnt2, cnt3 );
            }
         }
      }

   }

   return value;

}

void CheMPS2::ThreeDM::create_file() const{

   #ifdef CHEMPS2_MPI_COMPILATION
      const int mpi_rank = MPIchemps2::mpi_rank();
   #else
      const int mpi_rank = 0;
   #endif

   assert( disk == true );

   std::stringstream filename;
   filename << CheMPS2::THREE_RDM_storage_prefix << mpi_rank << ".h5";

   hid_t file_id  = H5Fcreate( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   hid_t group_id = H5Gcreate( file_id, "three_rdm", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   for ( int orb = 0; orb < L; orb++ ){

      std::stringstream storagename;
      storagename << "elements_" << orb;

      hsize_t dimarray   = array_size;
      hid_t dataspace_id = H5Screate_simple( 1, &dimarray, NULL );
      hid_t dataset_id   = H5Dcreate( group_id, storagename.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, elements );

      H5Dclose( dataset_id );
      H5Sclose( dataspace_id );

   }

   H5Gclose( group_id );
   H5Fclose( file_id );

}

void CheMPS2::ThreeDM::write_file( const int last_ham_orb ) const{

   #ifdef CHEMPS2_MPI_COMPILATION
      const int mpi_rank = MPIchemps2::mpi_rank();
   #else
      const int mpi_rank = 0;
   #endif

   assert( disk == true );

   std::stringstream filename;
   filename << CheMPS2::THREE_RDM_storage_prefix << mpi_rank << ".h5";

   hid_t file_id  = H5Fopen( filename.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
   hid_t group_id = H5Gopen( file_id, "three_rdm", H5P_DEFAULT );

      std::stringstream storagename;
      storagename << "elements_" << last_ham_orb;

      hid_t dataset_id = H5Dopen( group_id, storagename.str().c_str(), H5P_DEFAULT );
      H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, elements );

      H5Dclose( dataset_id );

   H5Gclose( group_id );
   H5Fclose( file_id );

}

void CheMPS2::ThreeDM::read_file( const int last_ham_orb ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const int mpi_rank = MPIchemps2::mpi_rank();
   #else
      const int mpi_rank = 0;
   #endif

   assert( disk == true );

   std::stringstream filename;
   filename << CheMPS2::THREE_RDM_storage_prefix << mpi_rank << ".h5";

   hid_t file_id  = H5Fopen( filename.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
   hid_t group_id = H5Gopen( file_id, "three_rdm", H5P_DEFAULT );

      std::stringstream storagename;
      storagename << "elements_" << last_ham_orb;

      hid_t dataset_id = H5Dopen( group_id, storagename.str().c_str(), H5P_DEFAULT );
      H5Dread( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, elements );

      H5Dclose( dataset_id );

   H5Gclose( group_id );
   H5Fclose( file_id );

}

void CheMPS2::ThreeDM::correct_higher_multiplicities(){

   if ( prob->gTwoS() != 0 ){
      double alpha = 1.0 / ( prob->gTwoS() + 1.0 );
      int inc1 = 1;
      if ( disk ){
         for ( int ham_orb = 0; ham_orb < L; ham_orb++ ){
            read_file( ham_orb );
            dscal_( &array_size, &alpha, elements, &inc1 );
            write_file( ham_orb );
         }
      } else {
         dscal_( &array_size, &alpha, elements, &inc1 );
      }
   }

}

void CheMPS2::ThreeDM::flush_disk(){

   assert( disk == true );
   for ( int ham_orb = 0; ham_orb < L; ham_orb++ ){

      read_file( ham_orb );

      for ( int counter = 0; counter < temp_disk_counter; counter++ ){

         const int orb1 = temp_disk_orbs[ 6 * counter + 0 ];
         const int orb2 = temp_disk_orbs[ 6 * counter + 1 ];
         const int orb3 = temp_disk_orbs[ 6 * counter + 2 ];
         const int orb4 = temp_disk_orbs[ 6 * counter + 3 ];
         const int orb5 = temp_disk_orbs[ 6 * counter + 4 ];
         const int orb6 = temp_disk_orbs[ 6 * counter + 5 ];
         const double value = temp_disk_vals[ counter ];

         if ( orb1 == ham_orb ){ elements[ orb5 + L * ( orb6 + L * ( orb4 + L * ( orb2 + L * orb3 ))) ] = value;
                                 elements[ orb6 + L * ( orb5 + L * ( orb4 + L * ( orb3 + L * orb2 ))) ] = value; }
         if ( orb2 == ham_orb ){ elements[ orb6 + L * ( orb4 + L * ( orb5 + L * ( orb3 + L * orb1 ))) ] = value;
                                 elements[ orb4 + L * ( orb6 + L * ( orb5 + L * ( orb1 + L * orb3 ))) ] = value; }
         if ( orb3 == ham_orb ){ elements[ orb4 + L * ( orb5 + L * ( orb6 + L * ( orb1 + L * orb2 ))) ] = value;
                                 elements[ orb5 + L * ( orb4 + L * ( orb6 + L * ( orb2 + L * orb1 ))) ] = value; }
         if ( orb4 == ham_orb ){ elements[ orb2 + L * ( orb3 + L * ( orb1 + L * ( orb5 + L * orb6 ))) ] = value;
                                 elements[ orb3 + L * ( orb2 + L * ( orb1 + L * ( orb6 + L * orb5 ))) ] = value; }
         if ( orb5 == ham_orb ){ elements[ orb3 + L * ( orb1 + L * ( orb2 + L * ( orb6 + L * orb4 ))) ] = value;
                                 elements[ orb1 + L * ( orb3 + L * ( orb2 + L * ( orb4 + L * orb6 ))) ] = value; }
         if ( orb6 == ham_orb ){ elements[ orb1 + L * ( orb2 + L * ( orb3 + L * ( orb4 + L * orb5 ))) ] = value;
                                 elements[ orb2 + L * ( orb1 + L * ( orb3 + L * ( orb5 + L * orb4 ))) ] = value; }
      }

      write_file( ham_orb );

   }

}

void CheMPS2::ThreeDM::save_HAM( const string filename ) const{

   assert( disk == false );
   save_HAM_generic( filename, L, "3-RDM", elements );

}

void CheMPS2::ThreeDM::save_HAM_generic( const string filename, const int LAS, const string tag, double * array ){

   hid_t   file_id      = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   long long linsize    = ( long long ) LAS;
   hsize_t dimarray     = ( linsize * linsize * linsize * linsize * linsize * linsize );
   hid_t   group_id     = H5Gcreate( file_id, tag.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   hid_t   dataspace_id = H5Screate_simple( 1, &dimarray, NULL );
   hid_t   dataset_id   = H5Dcreate( group_id, "elements", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array );

   H5Dclose( dataset_id );
   H5Sclose( dataspace_id );
   H5Gclose( group_id );
   H5Fclose( file_id );

   std::cout << "Saved the " << tag << " to the file " << filename << std::endl;

}

void CheMPS2::ThreeDM::fill_site( TensorT * denT, TensorL *** Ltensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors,
                                  Tensor3RDM **** dm3_a_J0_doublet, Tensor3RDM **** dm3_a_J1_doublet, Tensor3RDM **** dm3_a_J1_quartet,
                                  Tensor3RDM **** dm3_b_J0_doublet, Tensor3RDM **** dm3_b_J1_doublet, Tensor3RDM **** dm3_b_J1_quartet,
                                  Tensor3RDM **** dm3_c_J0_doublet, Tensor3RDM **** dm3_c_J1_doublet, Tensor3RDM **** dm3_c_J1_quartet,
                                  Tensor3RDM **** dm3_d_J0_doublet, Tensor3RDM **** dm3_d_J1_doublet, Tensor3RDM **** dm3_d_J1_quartet ){

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   temp_disk_counter = 0;
   const int orb_i = denT->gIndex();
   const int DIM = max(book->gMaxDimAtBound( orb_i ), book->gMaxDimAtBound( orb_i+1 ));
   const double sq3 = sqrt( 3.0 );

   #pragma omp parallel
   {

      double * workmem  = new double[ DIM * DIM ];
      double * workmem2 = new double[ DIM * DIM ];

      const int upperbound1 = ( orb_i * ( orb_i + 1 )) / 2;
      int jkl[] = { 0, 0, 0 };
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic) nowait
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for ( int global = 0; global < upperbound1; global++ ){
         Special::invert_triangle_two( global, jkl );
         const int orb_j = jkl[ 0 ];
         const int orb_k = jkl[ 1 ];
         if ( book->gIrrep( orb_j ) == book->gIrrep( orb_k )){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               const double d1 = diagram1( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], workmem );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_k, orb_i, orb_i,  4 * d1 );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_k, orb_i, -2 * d1 );
            }
         }
      }

      const int triangle3   = L - orb_i - 1;
      const int upperbound3 = ( triangle3 * ( triangle3 + 1 ) ) / 2;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic) nowait
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for ( int global = 0; global < upperbound3; global++ ){
         Special::invert_triangle_two( global, jkl );
         const int orb_j = L - 1 - jkl[ 1 ];
         const int orb_k = orb_j + jkl[ 0 ];
         if ( book->gIrrep( orb_j ) == book->gIrrep( orb_k )){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               const double d3 = diagram3( denT, F0tensors[orb_i][orb_k-orb_j][orb_j-1-orb_i], workmem );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_k, orb_i, orb_i,  4 * d3 );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_k, orb_i, -2 * d3 );
            }
         }
      }
      
      const int upperbound4 = ( orb_i * ( orb_i + 1 ) * ( orb_i + 2 ) ) / 6;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic) nowait
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for ( int global = 0; global < upperbound4; global++ ){
         Special::invert_triangle_three( global, jkl );
         const int orb_j = jkl[ 0 ];
         const int orb_k = jkl[ 1 ];
         const int orb_l = jkl[ 2 ];
         const int recalculate_global = orb_j + ( orb_k * ( orb_k + 1 ) ) / 2 + ( orb_l * ( orb_l + 1 ) * ( orb_l + 2 ) ) / 6;
         assert( global == recalculate_global );
         if ( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ) == Irreps::directProd( book->gIrrep( orb_l ), book->gIrrep( orb_i ) ) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
            #endif
            {
               const int cnt1 = orb_k - orb_j;
               const int cnt2 = orb_l - orb_k;
               const int cnt3 = orb_i - 1 - orb_l;
               if ( cnt2 > 0 ){
                  const double d4 =                        diagram4_5_6_7_8_9( denT, dm3_b_J0_doublet[cnt1][cnt2][cnt3], workmem, 'B' );
                  const double d5 = (( cnt1 == 0 ) ? 0.0 : diagram4_5_6_7_8_9( denT, dm3_b_J1_doublet[cnt1][cnt2][cnt3], workmem, 'B' ));
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_i, d4 + d5 );
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_i, d4 - d5 );
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_i, orb_l, -2 * d4 );
               }
               if ( cnt1 + cnt2 > 0 ){
                  const double d6 = diagram4_5_6_7_8_9( denT, dm3_c_J0_doublet[cnt1][cnt2][cnt3], workmem, 'C' );
                  const double d7 = diagram4_5_6_7_8_9( denT, dm3_c_J1_doublet[cnt1][cnt2][cnt3], workmem, 'C' );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_i, -2 * d6 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_i, d6 - d7 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_i, orb_k, d6 + d7 );
               }
               {
                  const double d8 = diagram4_5_6_7_8_9( denT, dm3_d_J0_doublet[cnt1][cnt2][cnt3], workmem, 'D' );
                  const double d9 = diagram4_5_6_7_8_9( denT, dm3_d_J1_doublet[cnt1][cnt2][cnt3], workmem, 'D' );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_j, orb_i, orb_i, -2 * d8 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_j, orb_i, d8 + d9 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_i, orb_j, d8 - d9 );
               }
            }
         }
      }
      
      const int upperbound10 = upperbound1 * ( L - orb_i - 1 );
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic) nowait
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int combined = 0; combined < upperbound10; combined++){
         const int global = combined % upperbound1;
         Special::invert_triangle_two( global, jkl );
         const int orb_l = orb_i + 1 + ( combined / upperbound1 );
         const int orb_j = jkl[ 0 ];
         const int orb_k = jkl[ 1 ];
         if ( Irreps::directProd(book->gIrrep( orb_j ), book->gIrrep( orb_k )) == Irreps::directProd(book->gIrrep( orb_l ), book->gIrrep( orb_i )) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
            #endif
            {
               const double d10 = diagram10( denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 );
               const double d11 = (( orb_j == orb_k ) ? 0.0 :
                                  diagram11( denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ));
               set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_i, d10 + d11 );
               set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_i, d10 - d11 );
               set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_i, orb_l, - 2 * d10 );
            }
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               {
                  const double d12 = diagram12( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 );
                  const double d13 = diagram13( denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_i, - 2 * d12 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_i, d12 - d13 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_i, orb_k, d12 + d13 );
               }
               if ( orb_j < orb_k ){
                  const double d14 = diagram14( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 );
                  const double d15 = diagram15( denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_j, orb_i, orb_i, - 2 * d14 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_j, orb_i, d14 - d15 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_i, orb_j, d14 + d15 );
               }
            }
         }
      }
      
      const int upperbound16 = upperbound3 * orb_i;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic) nowait
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for ( int combined = 0; combined < upperbound16; combined++ ){
         const int orb_j  = combined % orb_i;
         const int global = combined / orb_i;
         Special::invert_triangle_two( global, jkl );
         const int orb_m = L - 1 - jkl[ 1 ];
         const int orb_n = orb_m + jkl[ 0 ];
         if ( Irreps::directProd(book->gIrrep( orb_j ), book->gIrrep( orb_i )) == Irreps::directProd(book->gIrrep( orb_m ), book->gIrrep( orb_n )) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( orb_m, orb_n ) )
            #endif
            {
               const double d16 = diagram16( denT, Ltensors[orb_i-1][orb_i-1-orb_j], S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
               const double d17 = (( orb_m == orb_n ) ? 0.0 :
                                  diagram17( denT, Ltensors[orb_i-1][orb_i-1-orb_j], S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ));
               set_dmrg_index( orb_m, orb_n, orb_i, orb_j, orb_i, orb_i, d16 - d17 );
               set_dmrg_index( orb_m, orb_n, orb_i, orb_i, orb_j, orb_i, d16 + d17 );
               set_dmrg_index( orb_m, orb_n, orb_i, orb_i, orb_i, orb_j, - 2 * d16 );
            }
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_m, orb_n ) )
            #endif
            {
               {
                  const double d18 = diagram18( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                  const double d19 = diagram19( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_i, orb_i, d18 + d19 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_n, orb_i, - 2 * d18 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_i, orb_n, d18 - d19 );
               }
               if ( orb_m < orb_n ){
                  const double d20 = diagram20( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                  const double d21 = diagram21( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                  set_dmrg_index( orb_j, orb_n, orb_i, orb_m, orb_i, orb_i, d20 + d21 );
                  set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_m, orb_i, - 2 * d20 );
                  set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_i, orb_m, d20 - d21 );
               }
            }
         }
      }
      
      for ( int orb_m = orb_i+1; orb_m < L; orb_m++ ){
         for ( int orb_n = orb_m; orb_n < L; orb_n++ ){
         
            /*
             * Strategy for MPI of partitioning 2-2-2 and 3-1-2
             *   - owner of left renormalized operator is going to calculate
             *   - outer loop is (m,n)
             *   - in (m,n) loop make a temporary duplicate of S_mn & F_mn
            */
            
            const int irrep_mn  = Irreps::directProd( book->gIrrep( orb_m ), book->gIrrep( orb_n ) );
            const int irrep_imn = Irreps::directProd( book->gIrrep( orb_i ), irrep_mn );
            
            /************************************************************
             *  Make sure every process has the relevant S_mn and F_mn  *
             ************************************************************/
            #ifdef CHEMPS2_MPI_COMPILATION
            #pragma omp barrier // Everyone needs to be done before tensors are created and communicated
            #pragma omp single
            if ( orb_m > orb_i + 1 ){ // All processes own Fx/Sx[ index ][ n - m ][ m - i - 1 == 0 ]
               const int own_S_mn = MPIchemps2::owner_absigma( orb_m, orb_n );
               const int own_F_mn = MPIchemps2::owner_cdf(  L, orb_m, orb_n );
               if ( MPIRANK != own_F_mn ){ F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = new TensorF0( orb_i+1, irrep_mn, false, book );
                                           F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = new TensorF1( orb_i+1, irrep_mn, false, book ); }
               if ( MPIRANK != own_S_mn ){ S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = new TensorS0( orb_i+1, irrep_mn, false, book );
                    if ( orb_m != orb_n ){ S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = new TensorS1( orb_i+1, irrep_mn, false, book ); }}
                                      MPIchemps2::broadcast_tensor( F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], own_F_mn );
                                      MPIchemps2::broadcast_tensor( F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], own_F_mn );
                                      MPIchemps2::broadcast_tensor( S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], own_S_mn );
               if ( orb_m != orb_n ){ MPIchemps2::broadcast_tensor( S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], own_S_mn ); }
            }
            #endif
            
            /**************************************************************
             *  2-2-2 : Find out how many contributions each process has  *
             **************************************************************/
            int counter_Fjk = 0;
            int counter_Sjk = 0;
            for (int global = 0; global < upperbound1; global++){
               Special::invert_triangle_two( global, jkl );
               const int orb_j = jkl[ 0 ];
               const int orb_k = jkl[ 1 ];
               const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
               if ( irrep_jk == irrep_mn ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) ){ counter_Sjk++; }
                  if ( MPIRANK == MPIchemps2::owner_cdf(  L, orb_j, orb_k ) ){ counter_Fjk++; }
                  #else
                  counter_Fjk++;
                  counter_Sjk++;
                  #endif
               }
            }
            
            /*******************************
             *  2-2-2 : contributions F-F  *
             *******************************/
            if ( counter_Fjk > 0 ){
            
               TensorF0 * tens_29_33 = new TensorF0( orb_i, irrep_mn, true, book );
               TensorF1 * tens_30_32 = new TensorF1( orb_i, irrep_mn, true, book );
               TensorF1 * tens_34_40 = new TensorF1( orb_i, irrep_mn, true, book );
               TensorF0 * tens_35_41 = new TensorF0( orb_i, irrep_mn, true, book );
               TensorF1 * tens_36_42 = new TensorF1( orb_i, irrep_mn, true, book );
               TensorF1 * tens_37_44 = new TensorF1( orb_i, irrep_mn, true, book );
               TensorF1 * tens_38_43 = new TensorF1( orb_i, irrep_mn, true, book );
               fill_tens_29_33( denT, tens_29_33, F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               fill_tens_30_32( denT, tens_30_32, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               fill_tens_36_42( denT, tens_36_42, F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               fill_tens_34_35_37_38( denT, tens_34_40, tens_35_41, tens_37_44, tens_38_43, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
         
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for ( int global = 0; global < upperbound1; global++ ){
                  Special::invert_triangle_two( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
                  if ( irrep_jk == irrep_mn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
                     #endif
                     {
                        if ( orb_m < orb_n ){
                           const double d31 = tens_29_33->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d32 = tens_30_32->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d42 = tens_36_42->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d40 = tens_34_40->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d41 = tens_35_41->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d44 = tens_37_44->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           const double d43 = tens_38_43->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_k, orb_m, orb_i,   4*d31 );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_m, orb_k, orb_i, -2*(d31 + d32 + d40) );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_k, orb_i, orb_m, -2*(d31 + d41) );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_m, orb_i, orb_k,     d31 + d32 + d41 + d42 + d43 );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_k, orb_m,     d31 + d32 + d41 + d42 + d44 );
                           set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_m, orb_k, -2*(d31 + d42) );
                        }
                        const double d29 = tens_29_33->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d30 = tens_30_32->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d36 = tens_36_42->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d34 = tens_34_40->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d35 = tens_35_41->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d37 = tens_37_44->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        const double d38 = tens_38_43->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_n, orb_i,   4*d29 );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_k, orb_i, -2*(d29 + d30 + d34) );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_i, orb_n, -2*(d29 + d35) );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_i, orb_k,     d29 + d30 + d35 + d36 + d37 );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_k, orb_n,     d29 + d30 + d35 + d36 + d38 );
                        set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_n, orb_k, -2*(d29 + d36) );
                     }
                  }
               }
               
               delete tens_29_33;
               delete tens_30_32;
               delete tens_34_40;
               delete tens_35_41;
               delete tens_36_42;
               delete tens_37_44;
               delete tens_38_43;
               
            }
            
            /***********************************
             *  2-2-2 : contributions F-Sigma  *
             ***********************************/
            if ( counter_Fjk > 0 ){
            
               TensorF0 * tens_49_51 =                              new TensorF0( orb_i, irrep_mn, true, book );
               TensorF1 * tens_50_52 = (( orb_m == orb_n ) ? NULL : new TensorF1( orb_i, irrep_mn, true, book ));
                                      fill_tens_49_51( denT, tens_49_51, S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               if ( orb_m != orb_n ){ fill_tens_50_52( denT, tens_50_52, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem ); }
            
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for ( int global = 0; global < upperbound1; global++ ){
                  Special::invert_triangle_two( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
                  if ( irrep_jk == irrep_mn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
                     #endif
                     {
                        if ( orb_j < orb_k ){
                           const double d51 =                       tens_49_51->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' );
                           const double d52 = ((orb_m==orb_n)?0.0 : tens_50_52->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'N' ));
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_i, orb_i, - 2 * d51 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_j, orb_i, d51 + d52 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_i, orb_j, d51 - d52 );
                        }
                        const double d49 =                       tens_49_51->inproduct( F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' );
                        const double d50 = ((orb_m==orb_n)?0.0 : tens_50_52->inproduct( F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], 'T' ));
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_i, orb_i, - 2 * d49 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_k, orb_i, d49 + d50 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_i, orb_k, d49 - d50 );
                     }
                  }
               }
               
                                      delete tens_49_51;
               if ( orb_m != orb_n ){ delete tens_50_52; }
            
            }
            
            /***************************************
             *  2-2-2 : contributions Sigma-Sigma  *
             ***************************************/
            if ( counter_Sjk > 0 ){
            
               TensorS0 * tens_22 =                              new TensorS0( orb_i, irrep_mn, true, book );
               TensorS1 * tens_23 = (( orb_m == orb_n ) ? NULL : new TensorS1( orb_i, irrep_mn, true, book ));
               TensorS1 * tens_25 = (( orb_m == orb_n ) ? NULL : new TensorS1( orb_i, irrep_mn, true, book ));
               TensorS1 * tens_26 = (( orb_m == orb_n ) ? NULL : new TensorS1( orb_i, irrep_mn, true, book ));
               TensorS0 * tens_27 = (( orb_m == orb_n ) ? NULL : new TensorS0( orb_i, irrep_mn, true, book ));
               TensorS1 * tens_28 =                              new TensorS1( orb_i, irrep_mn, true, book );
               fill_tens_22_24( denT, tens_22, S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
                  fill_tens_28( denT, tens_28, S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               if ( orb_m != orb_n ){ fill_tens_23( denT, tens_23,                   S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
                                fill_tens_25_26_27( denT, tens_25, tens_26, tens_27, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ); }
            
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound1; global++){
                  Special::invert_triangle_two( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
                  if ( irrep_jk == irrep_mn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
                     #endif
                     {
                        const int cnt1   = orb_k - orb_j;
                        const double d22 =                                   tens_22->inproduct(S0tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N');
                        const double d23 = (((cnt1==0)||(orb_m==orb_n))?0.0: tens_23->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        const double d25 = (((cnt1==0)||(orb_m==orb_n))?0.0: tens_25->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        const double d26 = (((cnt1==0)||(orb_m==orb_n))?0.0: tens_26->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        const double d27 = (            (orb_m==orb_n) ?0.0: tens_27->inproduct(S0tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        const double d28 = ( (cnt1==0)                 ?0.0: tens_28->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_n, orb_i,  2*d22 + 2*d23 + d25             );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_n, orb_m, orb_i,  2*d22 - 2*d23 - d25             );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_i, orb_n,   -d22 -   d23 + d26 + d27 + d28 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_n, orb_i, orb_m,   -d22 +   d23 - d26 - d27 + d28 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_m, orb_n,   -d22 +   d23 - d26 + d27 - d28 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_n, orb_m,   -d22 -   d23 + d26 - d27 - d28 );
                     }
                  }
               }
               
                                      delete tens_22;
               if ( orb_m != orb_n ){ delete tens_23;
                                      delete tens_25;
                                      delete tens_26;
                                      delete tens_27; }
                                      delete tens_28;
            
            }
            
            /***********************************
             *  2-2-2 : contributions Sigma-F  *
             ***********************************/
            if ( counter_Sjk > 0 ){
            
               TensorS0 * tens_45 =                              new TensorS0( orb_i, irrep_mn, true, book );
               TensorS1 * tens_46 =                              new TensorS1( orb_i, irrep_mn, true, book );
               TensorS0 * tens_47 = (( orb_m == orb_n ) ? NULL : new TensorS0( orb_i, irrep_mn, true, book ));
               TensorS1 * tens_48 = (( orb_m == orb_n ) ? NULL : new TensorS1( orb_i, irrep_mn, true, book ));
                                      fill_tens_45_47( denT, tens_45, F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, true  );
                                      fill_tens_46_48( denT, tens_46, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, true  );
               if ( orb_m != orb_n ){ fill_tens_45_47( denT, tens_47, F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, false );
                                      fill_tens_46_48( denT, tens_48, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, false ); }
               
            
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound1; global++){
                  Special::invert_triangle_two( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
                  if ( irrep_jk == irrep_mn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
                     #endif
                     {
                        const int cnt1 = orb_k - orb_j;
                        if ( orb_m < orb_n ){
                           const double d47 =                        tens_47->inproduct(S0tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N');
                           const double d48 = (( cnt1 == 0 ) ? 0.0 : tens_48->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_m, orb_i, orb_i, d47 + d48 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_m, orb_i, d47 - d48 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_i, orb_m,  -2 * d47 );
                        }
                        const double d45 =                       tens_45->inproduct(S0tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N');
                        const double d46 = (( cnt1 == 0 ) ? 0.0: tens_46->inproduct(S1tensors[orb_i-1][cnt1][orb_i-1-orb_k],'N'));
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_n, orb_i, orb_i, d45 + d46 );
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_n, orb_i, d45 - d46 );
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_i, orb_n,  -2 * d45 );
                     }
                  }
               }
               
                                      delete tens_45;
                                      delete tens_46;
               if ( orb_m != orb_n ){ delete tens_47;
                                      delete tens_48; }
            
            }
            
            /**************************************************************
             *  3-1-2 : Find out how many contributions each process has  *
             **************************************************************/
            int counter_jkl = 0;
            for (int global = 0; global < upperbound4; global++){
               Special::invert_triangle_three( global, jkl );
               const int orb_j = jkl[ 0 ];
               const int orb_k = jkl[ 1 ];
               const int orb_l = jkl[ 2 ];
               const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
               if ( irrep_jkl == irrep_imn ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK ){ counter_jkl++; }
                  #else
                  counter_jkl++;
                  #endif
               }
            }
            
            /***********************************
             *  3-1-2 : contributions A-Sigma  *
             ***********************************/
            if ( counter_jkl > 0 ){
            
               Tensor3RDM * a_S0_doublet =                              new Tensor3RDM( orb_i, -2, 1, 3, irrep_imn, true, book );
               Tensor3RDM * a_S1_doublet = (( orb_m == orb_n ) ? NULL : new Tensor3RDM( orb_i, -2, 1, 3, irrep_imn, true, book ));
               Tensor3RDM * a_S1_quartet = (( orb_m == orb_n ) ? NULL : new Tensor3RDM( orb_i, -2, 3, 3, irrep_imn, true, book ));
                                      fill_a_S0( denT, a_S0_doublet,               S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               if ( orb_m != orb_n ){ fill_a_S1( denT, a_S1_doublet, a_S1_quartet, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ); }
               
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound4; global++){
                  Special::invert_triangle_three( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int orb_l = jkl[ 2 ];
                  const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
                  if ( irrep_jkl == irrep_imn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                     #endif
                     {
                        const int cnt1 = orb_k - orb_j;
                        const int cnt2 = orb_l - orb_k;
                        const int cnt3 = orb_i - 1 - orb_l;
                        if ( cnt1 + cnt2 > 0 ){
                           const double d90_95 =                                        a_S0_doublet->contract(dm3_a_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d93_98 = (                 (cnt1==0) ?0.0:  sq3*a_S0_doublet->contract(dm3_a_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d91_96 = (((orb_n==orb_m)||(cnt1==0))?0.0:      a_S1_doublet->contract(dm3_a_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d94_99 = ( (orb_n==orb_m)            ?0.0: -sq3*a_S1_doublet->contract(dm3_a_J0_doublet[cnt1][cnt2][cnt3]));
                           const double d92_97 = (((orb_n==orb_m)||(cnt1==0))?0.0:      a_S1_quartet->contract(dm3_a_J1_quartet[cnt1][cnt2][cnt3]));
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_n, orb_i, -2*d90_95 + 2*d91_96 + d92_97 );
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_n, orb_m, orb_i, -2*d90_95 - 2*d91_96 - d92_97 );
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_i, orb_n,    d90_95 +   d91_96 - d92_97 + d93_98 + d94_99 );
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_n, orb_i, orb_m,    d90_95 -   d91_96 + d92_97 + d93_98 - d94_99 );
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_m, orb_n,    d90_95 -   d91_96 + d92_97 - d93_98 + d94_99 );
                           set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_n, orb_m,    d90_95 +   d91_96 - d92_97 - d93_98 - d94_99 );
                        }
                     }
                  }
               }
               
                                      delete a_S0_doublet;
               if ( orb_m != orb_n ){ delete a_S1_doublet;
                                      delete a_S1_quartet; }
            
            }
            
            /*************************************
             *  3-1-2 : contributions BCD-Sigma  *
             *************************************/
            if ( counter_jkl > 0 ){
            
               Tensor3RDM * bcd_S0_doublet =                              new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               Tensor3RDM * bcd_S1_doublet = (( orb_m == orb_n ) ? NULL : new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book ));
               Tensor3RDM * bcd_S1_quartet = (( orb_m == orb_n ) ? NULL : new Tensor3RDM( orb_i, -2, 3, 1, irrep_imn, true, book ));
                                      fill_bcd_S0( denT, bcd_S0_doublet,                 S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               if ( orb_m != orb_n ){ fill_bcd_S1( denT, bcd_S1_doublet, bcd_S1_quartet, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ); }
               
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound4; global++){
                  Special::invert_triangle_three( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int orb_l = jkl[ 2 ];
                  const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
                  if ( irrep_jkl == irrep_imn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                     #endif
                     {
                        const int cnt1 = orb_k - orb_j;
                        const int cnt2 = orb_l - orb_k;
                        const int cnt3 = orb_i - 1 - orb_l;
                        if ( cnt2 > 0 ){ // dm3_b if NOT k==l
                           const double d120_125 =                                        bcd_S0_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d124_129 = (                 (cnt1==0) ?0.0:  sq3*bcd_S0_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d123_128 = ( (orb_n==orb_m)            ?0.0:  sq3*bcd_S1_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]));
                           const double d121_126 = (((orb_n==orb_m)||(cnt1==0))?0.0:      bcd_S1_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d122_127 = (((orb_n==orb_m)||(cnt1==0))?0.0:      bcd_S1_quartet->contract(dm3_b_J1_quartet[cnt1][cnt2][cnt3]));
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_j, orb_k, orb_i,  -d120_125 + 3*d121_126 +   d123_128 - d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_k, orb_j, orb_i,  -d120_125 - 3*d121_126 +   d123_128 + d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_j, orb_i, orb_k,  -d120_125 - 3*d121_126 -   d123_128 - d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_k, orb_i, orb_j,  -d120_125 + 3*d121_126 -   d123_128 + d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_i, orb_j, orb_k, 2*d120_125 + 2*d121_126 + 2*d122_127 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_i, orb_k, orb_j, 2*d120_125 - 2*d121_126 - 2*d122_127 );
                        }
                        if ( cnt1 + cnt2 > 0 ){ // dm3_c if NOT j==k==l
                           const double d110_115 =                           bcd_S0_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d112_117 =                      -sq3*bcd_S0_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]);
                           const double d111_116 = ((orb_n==orb_m)?0.0: -sq3*bcd_S1_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]));
                           const double d113_118 = ((orb_n==orb_m)?0.0:      bcd_S1_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d114_119 = ((orb_n==orb_m)?0.0:   -2*bcd_S1_quartet->contract(dm3_c_J1_quartet[cnt1][cnt2][cnt3]));
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_l, orb_i, 2*d110_115 + 2*d111_116 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_l, orb_j, orb_i,  -d110_115 -   d111_116 - d112_117 - 3*d113_118 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_i, orb_l, 2*d110_115 - 2*d111_116 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_l, orb_i, orb_j,  -d110_115 +   d111_116 - d112_117 + 3*d113_118 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_j, orb_l,  -d110_115 +   d111_116 + d112_117 +   d113_118 + d114_119 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_l, orb_j,  -d110_115 -   d111_116 + d112_117 -   d113_118 - d114_119 );
                        }
                        // dm3_d for all j <= k <= l
                        const double d100_105 =                                       -bcd_S0_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]);
                        const double d102_107 =                                   -sq3*bcd_S0_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]);
                        const double d101_106 = ( (orb_n==orb_m)            ?0.0:  sq3*bcd_S1_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]));
                        const double d103_108 = ( (orb_n==orb_m)            ?0.0:      bcd_S1_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]));
                        const double d104_109 = (((orb_n==orb_m)||(cnt2==0))?0.0:    2*bcd_S1_quartet->contract(dm3_d_J1_quartet[cnt1][cnt2][cnt3]));
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_l, orb_i, 2*d100_105 + 2*d101_106 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_l, orb_k, orb_i,  -d100_105 -   d101_106 - d102_107 - 3*d103_108 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_i, orb_l, 2*d100_105 - 2*d101_106 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_l, orb_i, orb_k,  -d100_105 +   d101_106 - d102_107 + 3*d103_108 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_k, orb_l,  -d100_105 +   d101_106 + d102_107 +   d103_108 + d104_109 );
                        set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_l, orb_k,  -d100_105 -   d101_106 + d102_107 -   d103_108 - d104_109 );
                     }
                  }
               }
               
                                      delete bcd_S0_doublet;
               if ( orb_m != orb_n ){ delete bcd_S1_doublet;
                                      delete bcd_S1_quartet; }
            
            }
            
            /***************************************
             *  3-1-2 : contributions F transpose  *
             ***************************************/
            if ( counter_jkl > 0 ){
            
               Tensor3RDM * F0_T_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               Tensor3RDM * F1_T_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               Tensor3RDM * F1_T_quartet = new Tensor3RDM( orb_i, -2, 3, 1, irrep_imn, true, book );
               fill_F0_T( denT, F0_T_doublet,               F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               fill_F1_T( denT, F1_T_doublet, F1_T_quartet, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
               
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound4; global++){
                  Special::invert_triangle_three( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int orb_l = jkl[ 2 ];
                  const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
                  if ( irrep_jkl == irrep_imn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                     #endif
                     {
                        const int cnt1 = orb_k - orb_j;
                        const int cnt2 = orb_l - orb_k;
                        const int cnt3 = orb_i - 1 - orb_l;
                        if ( cnt2 > 0 ){ // dm3_b if NOT k==l
                           const double d130_135 =                      F0_T_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d131_136 = ((cnt1==0)?0.0:  sq3*F0_T_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d132_137 =                  sq3*F1_T_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d133_138 = ((cnt1==0)?0.0:      F1_T_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d134_139 = ((cnt1==0)?0.0:      F1_T_quartet->contract(dm3_b_J1_quartet[cnt1][cnt2][cnt3]));
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_l, orb_n, orb_i,    d130_135 +   d131_136 + d132_137 + 3*d133_138 );
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_n, orb_l, orb_i,    d130_135 -   d131_136 + d132_137 - 3*d133_138 );
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_l, orb_i, orb_n, -2*d130_135 - 2*d131_136 );
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_n, orb_i, orb_l,    d130_135 +   d131_136 - d132_137 +   d133_138 + 2*d134_139 );
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_l, orb_n, -2*d130_135 + 2*d131_136 );
                           set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_n, orb_l,    d130_135 -   d131_136 - d132_137 -   d133_138 - 2*d134_139 );
                        }
                        if ( cnt1 + cnt2 > 0 ){ // dm3_c if NOT j==k==l
                           const double d150_155 =      F0_T_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d151_156 = -sq3*F0_T_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]);
                           const double d152_157 =  sq3*F1_T_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d153_158 =      F1_T_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]);
                           const double d154_159 =     -F1_T_quartet->contract(dm3_c_J1_quartet[cnt1][cnt2][cnt3]);
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_k, orb_n, orb_i, -2*d150_155 - 2*d152_157 );
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_n, orb_k, orb_i,    d150_155 + d151_156 + d152_157 - 3*d153_158 );
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_k, orb_i, orb_n,  4*d150_155 );
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_n, orb_i, orb_k, -2*d150_155 + 2*d153_158 + 2*d154_159 );
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_k, orb_n, -2*d150_155 - 2*d151_156 );
                           set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_n, orb_k,    d150_155 + d151_156 + d152_157 + d153_158 - 2*d154_159 );
                        }
                        // dm3_d for all j <= k <= l
                        const double d170_175 =                -F0_T_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]);
                        const double d171_176 =            -sq3*F0_T_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]);
                        const double d172_177 =            -sq3*F1_T_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]);
                        const double d173_178 =                 F1_T_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]);
                        const double d174_179 = ((cnt2==0)?0.0: F1_T_quartet->contract(dm3_d_J1_quartet[cnt1][cnt2][cnt3]));
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_j, orb_n, orb_i, -2*d170_175 - 2*d172_177 );
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_n, orb_j, orb_i,    d170_175 + d171_176 + d172_177 - 3*d173_178 );
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_j, orb_i, orb_n,  4*d170_175 );
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_n, orb_i, orb_j, -2*d170_175 + 2*d173_178 + 2*d174_179 );
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_j, orb_n, -2*d170_175 - 2*d171_176 );
                        set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_n, orb_j,    d170_175 + d171_176 + d172_177 + d173_178 - 2*d174_179 );
                     }
                  }
               }
               
               delete F0_T_doublet;
               delete F1_T_doublet;
               delete F1_T_quartet;
            
            }
            
            /*************************************
             *  3-1-2 : contributions F regular  *
             *************************************/
            if (( orb_m != orb_n ) && ( counter_jkl > 0 )){
            
               Tensor3RDM * F0_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               Tensor3RDM * F1_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               Tensor3RDM * F1_quartet = new Tensor3RDM( orb_i, -2, 3, 1, irrep_imn, true, book );
               fill_F0( denT, F0_doublet,             F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               fill_F1( denT, F1_doublet, F1_quartet, F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
               
               #ifdef CHEMPS2_MPI_COMPILATION
                  #pragma omp for schedule(dynamic) nowait
               #else
                  #pragma omp for schedule(static) nowait
               #endif
               for (int global = 0; global < upperbound4; global++){
                  Special::invert_triangle_three( global, jkl );
                  const int orb_j = jkl[ 0 ];
                  const int orb_k = jkl[ 1 ];
                  const int orb_l = jkl[ 2 ];
                  const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
                  if ( irrep_jkl == irrep_imn ){
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                     #endif
                     {
                        const int cnt1 = orb_k - orb_j;
                        const int cnt2 = orb_l - orb_k;
                        const int cnt3 = orb_i - 1 - orb_l;
                        if ( cnt2 > 0 ){ // dm3_b if NOT k==l
                           const double d140_145 =                      F0_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d141_146 = ((cnt1==0)?0.0:  sq3*F0_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d142_147 =                  sq3*F1_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d143_148 = ((cnt1==0)?0.0:      F1_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3]));
                           const double d144_149 = ((cnt1==0)?0.0:      F1_quartet->contract(dm3_b_J1_quartet[cnt1][cnt2][cnt3]));
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_l, orb_m, orb_i,    d140_145 +   d141_146 + d142_147 + 3*d143_148 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_m, orb_l, orb_i,    d140_145 -   d141_146 + d142_147 - 3*d143_148 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_l, orb_i, orb_m, -2*d140_145 - 2*d141_146 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_m, orb_i, orb_l,    d140_145 +   d141_146 - d142_147 +   d143_148 + 2*d144_149 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_l, orb_m, -2*d140_145 + 2*d141_146 );
                           set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_m, orb_l,    d140_145 -   d141_146 - d142_147 -   d143_148 - 2*d144_149 );
                        }
                        if ( cnt1 + cnt2 > 0 ){ // dm3_c if NOT j==k==l
                           const double d160_165 =      F0_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d161_166 = -sq3*F0_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]);
                           const double d162_167 =  sq3*F1_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3]);
                           const double d163_168 =      F1_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3]);
                           const double d164_169 =     -F1_quartet->contract(dm3_c_J1_quartet[cnt1][cnt2][cnt3]);
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_k, orb_m, orb_i, -2*d160_165 - 2*d162_167 );
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_m, orb_k, orb_i,    d160_165 + d161_166 + d162_167 - 3*d163_168 );
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_k, orb_i, orb_m,  4*d160_165 );
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_m, orb_i, orb_k, -2*d160_165 + 2*d163_168 + 2*d164_169 );
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_i, orb_k, orb_m, -2*d160_165 - 2*d161_166 );
                           set_dmrg_index( orb_j, orb_l, orb_n, orb_i, orb_m, orb_k,    d160_165 + d161_166 + d162_167 + d163_168 - 2*d164_169 );
                        }
                        // dm3_d for all j <= k <= l
                        const double d180_185 =                 -F0_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]);
                        const double d181_186 =             -sq3*F0_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]);
                        const double d182_187 =             -sq3*F1_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3]);
                        const double d183_188 =                  F1_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3]);
                        const double d184_189 = ((cnt2==0)?0.0:  F1_quartet->contract(dm3_d_J1_quartet[cnt1][cnt2][cnt3]));
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_j, orb_m, orb_i, -2*d180_185 - 2*d182_187 );
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_m, orb_j, orb_i,    d180_185 + d181_186 + d182_187 - 3*d183_188 );
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_j, orb_i, orb_m,  4*d180_185 );
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_m, orb_i, orb_j, -2*d180_185 + 2*d183_188 + 2*d184_189 );
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_i, orb_j, orb_m, -2*d180_185 - 2*d181_186 );
                        set_dmrg_index( orb_k, orb_l, orb_n, orb_i, orb_m, orb_j,    d180_185 + d181_186 + d182_187 + d183_188 - 2*d184_189 );
                     }
                  }
               }
               
               delete F0_doublet;
               delete F1_doublet;
               delete F1_quartet;
               
            }
            
            /***************************************
             *  Clean up the double S_mn and F_mn  *
             ***************************************/
            #ifdef CHEMPS2_MPI_COMPILATION
            #pragma omp barrier // Everyone needs to be done before tensors are deleted
            #pragma omp single
            if ( orb_m > orb_i + 1 ){ // All processes own Fx/Sx[ index ][ n - m ][ m - i - 1 == 0 ]
               const int own_S_mn = MPIchemps2::owner_absigma( orb_m, orb_n );
               const int own_F_mn = MPIchemps2::owner_cdf(  L, orb_m, orb_n );
               if ( MPIRANK != own_F_mn ){ delete F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i]; F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = NULL;
                                           delete F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i]; F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = NULL; }
               if ( MPIRANK != own_S_mn ){ delete S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i]; S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = NULL;
                    if ( orb_m != orb_n ){ delete S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i]; S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i] = NULL; } }
            }
            #endif
            
         }
      }
      
      for ( int orb_m = orb_i+1; orb_m < L; orb_m++ ){
      
         const int irrep_m = book->gIrrep( orb_m );
      
         /**************************************************************
          *  3-2-1 : Find out how many contributions each process has  *
          **************************************************************/
         int counter_jkl = 0;
         for (int global = 0; global < upperbound4; global++){
            Special::invert_triangle_three( global, jkl );
            const int orb_j = jkl[ 0 ];
            const int orb_k = jkl[ 1 ];
            const int orb_l = jkl[ 2 ];
            const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
            if ( irrep_jkl == irrep_m ){
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK ){ counter_jkl++; }
               #else
               counter_jkl++;
               #endif
            }
         }
         
         /********************
          *  3-2-1 : part 1  *
          ********************/
         if ( counter_jkl > 0 ){
         
            Tensor3RDM *   A_T_doublet = new Tensor3RDM( orb_i, -2, 1, 3, irrep_m, true, book );
            Tensor3RDM * BCD_T_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
               fill_53_54( denT,   A_T_doublet, Ltensors[orb_i][orb_m-1-orb_i], workmem );
            fill_55_to_60( denT, BCD_T_doublet, Ltensors[orb_i][orb_m-1-orb_i], workmem );
         
            #ifdef CHEMPS2_MPI_COMPILATION
               #pragma omp for schedule(dynamic) nowait
            #else
               #pragma omp for schedule(static) nowait
            #endif
            for ( int global = 0; global < upperbound4; global++ ){
               Special::invert_triangle_three( global, jkl );
               const int orb_j = jkl[ 0 ];
               const int orb_k = jkl[ 1 ];
               const int orb_l = jkl[ 2 ];
               const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
               if ( irrep_jkl == irrep_m ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                  #endif
                  {
                     const int cnt1 = orb_k - orb_j;
                     const int cnt2 = orb_l - orb_k;
                     const int cnt3 = orb_i - 1 - orb_l;
                     if ( cnt1 + cnt2 > 0 ){
                        const double d53 =     A_T_doublet->contract( dm3_a_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d54 = sq3*A_T_doublet->contract( dm3_a_J1_doublet[cnt1][cnt2][cnt3] );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_i, orb_i, d53 - d54 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_m, orb_i, d53 + d54 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_i, orb_m, - 2 * d53 );
                     }
                     if ( cnt2 > 0 ){
                        const double d55 =     BCD_T_doublet->contract( dm3_b_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d56 = sq3*BCD_T_doublet->contract( dm3_b_J1_doublet[cnt1][cnt2][cnt3] );
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_l, orb_i, orb_i, d55 + d56 );
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_l, orb_i, d55 - d56 );
                        set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_i, orb_l, - 2 * d55 );
                     }
                     if ( cnt1 + cnt2 > 0 ){
                        const double d57 =     BCD_T_doublet->contract( dm3_c_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d58 = sq3*BCD_T_doublet->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        set_dmrg_index( orb_j, orb_l, orb_m, orb_k, orb_i, orb_i, - 2 * d57 );
                        set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_k, orb_i, d57 - d58 );
                        set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_i, orb_k, d57 + d58 );
                     }
                     const double d59 =     -BCD_T_doublet->contract( dm3_d_J0_doublet[cnt1][cnt2][cnt3] );
                     const double d60 = -sq3*BCD_T_doublet->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     set_dmrg_index( orb_k, orb_l, orb_m, orb_j, orb_i, orb_i, - 2 * d59 );
                     set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_j, orb_i, d59 + d60 );
                     set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_i, orb_j, d59 - d60 );
                  }
               }
            }
            
            delete   A_T_doublet;
            delete BCD_T_doublet;
         
         }
         
         /**************************************************************
          *  1-4-1 : Find out how many contributions each process has  *
          **************************************************************/
         int counter_j = 0;
         for ( int orb_j = 0; orb_j < orb_i; orb_j++ ){
            const int irrep_j = book->gIrrep( orb_j );
            if ( irrep_j == irrep_m ){
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIRANK == MPIchemps2::owner_q( L, orb_j ) ){ counter_j++; }
               #else
               counter_j++;
               #endif
            }
         }
         
         /***********************
          *  3-2-1 : part 2     *
          *  1-4-1 : L^T - L^T  *
          ***********************/
         if (( counter_jkl > 0 ) || ( counter_j > 0 )){
         
            Tensor3RDM * tens_61 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_63 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_65 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_67 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_68 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_76 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_77 = new Tensor3RDM( orb_i, -2, 1, 1, irrep_m, true, book );
            Tensor3RDM * tens_69 = new Tensor3RDM( orb_i, -2, 3, 1, irrep_m, true, book );
            Tensor3RDM * tens_78 = new Tensor3RDM( orb_i, -2, 3, 1, irrep_m, true, book );
            Tensor3RDM * tens_79 = new Tensor3RDM( orb_i, -2, 3, 1, irrep_m, true, book );
            fill_61( denT, tens_61, Ltensors[orb_i][orb_m-1-orb_i], workmem );
            fill_63_65( denT, tens_63, tens_65, tens_67, tens_68, tens_76, tens_77, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2 );
            fill_69_78_79( denT, tens_69, tens_78, tens_79, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2 );
            
            #ifdef CHEMPS2_MPI_COMPILATION
               #pragma omp for schedule(dynamic) nowait
            #else
               #pragma omp for schedule(static) nowait
            #endif
            for ( int orb_j = 0; orb_j < orb_i; orb_j++ ){
               const int irrep_j = book->gIrrep( orb_j );
               if ( irrep_j == irrep_m ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_q( L, orb_j ) ) //Everyone owns the L-tensors --> task division based on Q-tensor ownership
                  #endif
                  {
                     const double d2 = sqrt( 2.0 ) * tens_61->inproduct( Ltensors[orb_i-1][orb_i-1-orb_j], 'N' );
                     set_dmrg_index( orb_j, orb_i, orb_i, orb_m, orb_i, orb_i, 2 * d2 );
                     set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_m, orb_i,   - d2 );
                  }
               }
            }
      
            #ifdef CHEMPS2_MPI_COMPILATION
               #pragma omp for schedule(dynamic) nowait
            #else
               #pragma omp for schedule(static) nowait
            #endif
            for ( int global = 0; global < upperbound4; global++ ){
               Special::invert_triangle_three( global, jkl );
               const int orb_j = jkl[ 0 ];
               const int orb_k = jkl[ 1 ];
               const int orb_l = jkl[ 2 ];
               const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
               if ( irrep_jkl == irrep_m ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                  #endif
                  {
                     const int cnt1 = orb_k - orb_j;
                     const int cnt2 = orb_l - orb_k;
                     const int cnt3 = orb_i - 1 - orb_l;
                     if ( cnt2 > 0 ){
                        const double d61 =     tens_61->contract( dm3_b_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d62 = sq3*tens_61->contract( dm3_b_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d63 =     tens_63->contract( dm3_b_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d64 = sq3*tens_63->contract( dm3_b_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d65 =     tens_65->contract( dm3_b_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d66 = sq3*tens_65->contract( dm3_b_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d67 =     tens_67->contract( dm3_b_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d68 =     tens_68->contract( dm3_b_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d69 =     tens_69->contract( dm3_b_J1_quartet[cnt1][cnt2][cnt3] );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_m, orb_i, -2*d61 - 2*d62 + d63 + d64 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_l, orb_i, -2*d61 + 2*d62 + d63 - d64 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_m, d61 + d62 + d65 + d66 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_i, orb_l, d61 - d62 + d67 + d68 + d69 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_m, orb_l, d61 + d62 + d67 - d68 - d69 );
                        set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_m, d61 - d62 + d65 - d66 );
                     }
                     if ( cnt1 + cnt2 > 0 ){
                        const double d70 =      tens_61->contract( dm3_c_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d71 =  sq3*tens_61->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d72 =     -tens_63->contract( dm3_c_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d73 = -sq3*tens_63->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d74 =     -tens_65->contract( dm3_c_J0_doublet[cnt1][cnt2][cnt3] );
                        const double d75 = -sq3*tens_65->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d76 =      tens_76->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d77 =      tens_77->contract( dm3_c_J1_doublet[cnt1][cnt2][cnt3] );
                        const double d78 =      tens_78->contract( dm3_c_J1_quartet[cnt1][cnt2][cnt3] );
                        const double d79 =      tens_79->contract( dm3_c_J1_quartet[cnt1][cnt2][cnt3] );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_m, orb_i, 4*d70 + 2*d72 );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_m, orb_k, orb_i, -2*d70 + 2*d71 - d72 + d73 );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_m, -2*d70 + 2*d74 );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_m, orb_i, orb_k, d70 - d71 - d74 + d76 + d78 );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_m, d70 - d71 - d74 + d75 );
                        set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_m, orb_k, -2*d70 - d72 + d77 + d79 );
                     }
                     const double d80 =     -tens_61->contract( dm3_d_J0_doublet[cnt1][cnt2][cnt3] );
                     const double d81 = -sq3*tens_61->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     const double d82 =      tens_63->contract( dm3_d_J0_doublet[cnt1][cnt2][cnt3] );
                     const double d83 =  sq3*tens_63->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     const double d84 =      tens_65->contract( dm3_d_J0_doublet[cnt1][cnt2][cnt3] );
                     const double d85 =  sq3*tens_65->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     const double d86 =      tens_76->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     const double d87 =      tens_77->contract( dm3_d_J1_doublet[cnt1][cnt2][cnt3] );
                     const double d88 =     -tens_78->contract( dm3_d_J1_quartet[cnt1][cnt2][cnt3] );
                     const double d89 =     -tens_79->contract( dm3_d_J1_quartet[cnt1][cnt2][cnt3] );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_l, orb_i, 4*d80 + 2*d82 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_l, orb_k, orb_i, -2*d80 - 2*d81 - d82 - d83 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_i, orb_l, -2*d80 + 2*d84 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_l, orb_i, orb_k, d80 + d81 - d84 - d85 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_k, orb_l, d80 + d81 - d84 + d86 + d88 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_l, orb_k, -2*d80 - d82 + d87 + d89 );
                  }
               }
            }
            
            delete tens_61;
            delete tens_63;
            delete tens_65;
            delete tens_67;
            delete tens_68;
            delete tens_69;
            delete tens_76;
            delete tens_77;
            delete tens_78;
            delete tens_79;
         }
      
      }
      
      delete [] workmem;
      delete [] workmem2;

   }

   if (( disk ) && ( temp_disk_counter > 0 )){ flush_disk(); }

}

double CheMPS2::ThreeDM::diagram1(TensorT * denT, TensorF0 * denF0, double * workmem) const{

   assert( denF0->get_irrep() == 0 );
   const int orb_i = denT->gIndex();
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
               
            int dimL = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimR = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimL > 0 ) && ( dimR > 0 )){
               
               double * Tblock =  denT->gStorage( NL, TwoSL, IL, NL+2, TwoSL, IL );
               double * Fblock = denF0->gStorage( NL, TwoSL, IL, NL,   TwoSL, IL );
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta  = 0.0; //set
               dgemm_( &notrans, &notrans, &dimL, &dimR, &dimL, &alpha, Fblock, &dimL, Tblock, &dimL, &beta, workmem, &dimL );
               
               int length = dimL * dimR;
               int inc = 1;
               total += ( TwoSL + 1 ) * ddot_( &length, workmem, &inc, Tblock, &inc );

            }
         }
      }
   }
   
   total *= sqrt( 0.5 );
   return total;

}

double CheMPS2::ThreeDM::diagram3(TensorT * denT, TensorF0 * denF0, double * workmem) const{

   assert( denF0->get_irrep() == 0 );
   const int orb_i = denT->gIndex();
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
               
            int dimL = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimR = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimL > 0 ) && ( dimR > 0 )){
               
               double * Tblock =  denT->gStorage( NL,   TwoSL, IL, NL+2, TwoSL, IL );
               double * Fblock = denF0->gStorage( NL+2, TwoSL, IL, NL+2, TwoSL, IL );
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta  = 0.0; //set
               dgemm_( &notrans, &notrans, &dimL, &dimR, &dimR, &alpha, Tblock, &dimL, Fblock, &dimR, &beta, workmem, &dimL );
               
               int length = dimL * dimR;
               int inc = 1;
               total += ( TwoSL + 1 ) * ddot_( &length, workmem, &inc, Tblock, &inc );

            }
         }
      }
   }
   
   total *= sqrt( 0.5 );
   return total;

}

double CheMPS2::ThreeDM::diagram4_5_6_7_8_9(TensorT * denT, Tensor3RDM * d3tens, double * workmem, const char type) const{

   const int orb_i = denT->gIndex();
   assert( d3tens->get_irrep()  == book->gIrrep( orb_i ) );
   assert( d3tens->get_two_j2() == 1 );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
               
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               const int ILdown = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
            
                  int dimR     = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILdown );
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILdown );
                  
                  if (( dimLdown > 0 ) && ( dimR > 0 )){
                     
                     double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,     NL+1, TwoSLprime, ILdown );
                     double * Tdown =   denT->gStorage( NL-1, TwoSLprime, ILdown, NL+1, TwoSLprime, ILdown );
                     double * block = d3tens->gStorage( NL-1, TwoSLprime, ILdown, NL,   TwoSL,      IL     );
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimR, &dimLup, &alpha, block, &dimLdown, Tup, &dimLup, &beta, workmem, &dimLdown );
                     
                     int length = dimLdown * dimR;
                     int inc = 1;
                     const double factor = ((type =='D') ? ( sqrt( 0.5 * ( d3tens->get_two_j1() + 1 ) ) * ( TwoSLprime + 1 ) )
                                                         : ( Special::phase( TwoSL + 1 - TwoSLprime )
                                                           * sqrt( 0.5 * ( d3tens->get_two_j1() + 1 ) * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) ));
                     total += factor * ddot_( &length, workmem, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram10(TensorT * denT, TensorS0 * denS0, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denS0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denS0->get_irrep()    );
            
            int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSL, ILxIixIl );
            int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSL, ILxIixIl );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
            
               double * Tdown  =  denT->gStorage( NL-2, TwoSL, ILxIixIl, NL, TwoSL, ILxIixIl );
               double * Sblock = denS0->gStorage( NL-2, TwoSL, ILxIixIl, NL, TwoSL, IL       );
         
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                  
                  if ( dimRup > 0 ){
                     
                     double * Tup    = denT->gStorage( NL, TwoSL, IL,       NL+1, TwoSR, ILxIi );
                     double * Lblock = denL->gStorage( NL, TwoSL, ILxIixIl, NL+1, TwoSR, ILxIi );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Sblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRdown, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total -= ( TwoSR + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   total *= sqrt( 0.5 );
   return total;

}

double CheMPS2::ThreeDM::diagram11(TensorT * denT, TensorS1 * denS1, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denS1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denS1->get_irrep()    );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
                  int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSLprime, ILxIixIl );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSLprime, ILxIixIl );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  =  denT->gStorage( NL-2, TwoSLprime, ILxIixIl, NL, TwoSLprime, ILxIixIl );
                     double * Sblock = denS1->gStorage( NL-2, TwoSLprime, ILxIixIl, NL, TwoSL,      IL       );
               
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                     
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if (( dimRup > 0 ) && ( abs( TwoSLprime - TwoSR ) == 1 )){
                           
                           double * Tup    = denT->gStorage( NL, TwoSL,      IL,       NL+1, TwoSR, ILxIi );
                           double * Lblock = denL->gStorage( NL, TwoSLprime, ILxIixIl, NL+1, TwoSR, ILxIi );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Sblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                           dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRdown, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                  * Special::phase( TwoSR + 3 + TwoSLprime )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );

                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram12(TensorT * denT, TensorF0 * denF0, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denF0->get_irrep()    );
            
            int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSL, ILxIixIl );
            int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSL, ILxIixIl );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
            
               double * Tdown  =  denT->gStorage( NL, TwoSL, ILxIixIl, NL+2, TwoSL, ILxIixIl );
               double * Fblock = denF0->gStorage( NL, TwoSL, ILxIixIl, NL,   TwoSL, IL       );
         
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                  
                  if ( dimRup > 0 ){
                     
                     double * Tup    = denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR, ILxIi    );
                     double * Lblock = denL->gStorage( NL+1, TwoSR, ILxIi, NL+2, TwoSL, ILxIixIl );
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Fblock,  &dimLdown, Tup,    &dimLup, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRup, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total += sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSR + 1 ) )
                            * Special::phase( TwoSL + 3 - TwoSR )
                            * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram13(TensorT * denT, TensorF1 * denF1, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denF1->get_irrep()    );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
                  int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLprime, ILxIixIl );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSLprime, ILxIixIl );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  =  denT->gStorage( NL, TwoSLprime, ILxIixIl, NL+2, TwoSLprime, ILxIixIl );
                     double * Fblock = denF1->gStorage( NL, TwoSLprime, ILxIixIl, NL,   TwoSL,      IL       );
               
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                     
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if (( dimRup > 0 ) && ( abs( TwoSLprime - TwoSR ) == 1 )){
                           
                           double * Tup    = denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR,      ILxIi    );
                           double * Lblock = denL->gStorage( NL+1, TwoSR, ILxIi, NL+2, TwoSLprime, ILxIixIl );
                           
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Fblock,  &dimLdown, Tup,    &dimLup, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRup, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSL + 1 ) * ( TwoSR + 1 ) * ( TwoSLprime + 1 ) )
                                  * Special::phase( 2 * TwoSR + 2 )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );

                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram14(TensorT * denT, TensorF0 * denF0, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denF0->get_irrep()    );
            
            int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSL, ILxIixIl );
            int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSL, ILxIixIl );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
            
               double * Tdown  =  denT->gStorage( NL, TwoSL, ILxIixIl, NL+2, TwoSL, ILxIixIl );
               double * Fblock = denF0->gStorage( NL, TwoSL, IL,       NL,   TwoSL, ILxIixIl );
         
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                  
                  if ( dimRup > 0 ){
                     
                     double * Tup    = denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR, ILxIi    );
                     double * Lblock = denL->gStorage( NL+1, TwoSR, ILxIi, NL+2, TwoSL, ILxIixIl );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_(   &trans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Fblock,  &dimLup,   Tup,    &dimLup, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRup, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total += sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSR + 1 ) )
                            * Special::phase( TwoSL + 3 - TwoSR )
                            * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram15(TensorT * denT, TensorF1 * denF1, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIixIl = Irreps::directProd( IL, denF1->get_irrep()    );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
                  int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLprime, ILxIixIl );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSLprime, ILxIixIl );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  =  denT->gStorage( NL, TwoSLprime, ILxIixIl, NL+2, TwoSLprime, ILxIixIl );
                     double * Fblock = denF1->gStorage( NL, TwoSL,      IL,       NL,   TwoSLprime, ILxIixIl );
               
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                     
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if (( dimRup > 0 ) && ( abs( TwoSLprime - TwoSR ) == 1 )){
                           
                           double * Tup    = denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR,      ILxIi    );
                           double * Lblock = denL->gStorage( NL+1, TwoSR, ILxIi, NL+2, TwoSLprime, ILxIixIl );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_(   &trans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Fblock,  &dimLup,   Tup,    &dimLup, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Lblock, &dimRup, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSR + 1 ) ) * ( TwoSLprime + 1 )
                                  * Special::phase( TwoSL + TwoSLprime )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );

                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram16(TensorT * denT, TensorL * denL, TensorS0 * denS0, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denS0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL+1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+3, TwoSLprime, ILxIj );
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi );
                  
                  if (( dimRup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup    =  denT->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIi  );
                     double * Tdown  =  denT->gStorage( NL+1, TwoSLprime, ILxIj, NL+3, TwoSLprime, ILxIj  );
                     double * Sblock = denS0->gStorage( NL+1, TwoSLprime, ILxIi, NL+3, TwoSLprime, ILxIj  );
                     double * Lblock =  denL->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIj  );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_(   &trans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLup,   Tup,    &dimLup, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Sblock, &dimRup, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total -= ( TwoSLprime + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   total *= sqrt( 0.5 );
   return total;

}

double CheMPS2::ThreeDM::diagram17(TensorT * denT, TensorL * denL, TensorS1 * denS1, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denS1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL+1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+3, TwoSLprime, ILxIj );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  = denT->gStorage( NL+1, TwoSLprime, ILxIj, NL+3, TwoSLprime, ILxIj );
                     double * Lblock = denL->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIj );
                  
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if ( dimRup > 0 ){
                        
                           double * Tup    =  denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR,      ILxIi );
                           double * Sblock = denS1->gStorage( NL+1, TwoSR, ILxIi, NL+3, TwoSLprime, ILxIj );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_(   &trans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLup,   Tup,    &dimLup, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Sblock, &dimRup, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSR + 1 ) ) * ( TwoSLprime + 1 )
                                  * Special::phase( TwoSL + TwoSLprime + 3 )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );

                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram18(TensorT * denT, TensorL * denL, TensorF0 * denF0, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIj );
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi );
                  
                  if (( dimRup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup    =  denT->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIi  );
                     double * Tdown  =  denT->gStorage( NL-1, TwoSLprime, ILxIj, NL+1, TwoSLprime, ILxIj  );
                     double * Fblock = denF0->gStorage( NL+1, TwoSLprime, ILxIj, NL+1, TwoSLprime, ILxIi  );
                     double * Lblock =  denL->gStorage( NL-1, TwoSLprime, ILxIj, NL,   TwoSL,      IL     );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Fblock, &dimRdown, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total += sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) )
                            * Special::phase( TwoSL + 1 - TwoSLprime )
                            * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram19(TensorT * denT, TensorL * denL, TensorF1 * denF1, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIj );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  = denT->gStorage( NL-1, TwoSLprime, ILxIj, NL+1, TwoSLprime, ILxIj );
                     double * Lblock = denL->gStorage( NL-1, TwoSLprime, ILxIj, NL,   TwoSL,      IL    );
                  
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if ( dimRup > 0 ){
                     
                           double * Tup    =  denT->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSR, ILxIi );
                           double * Fblock = denF1->gStorage( NL+1, TwoSLprime, ILxIj, NL+1, TwoSR, ILxIi );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                           dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Fblock, &dimRdown, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) * ( TwoSR + 1 ) )
                                  * Special::phase( 2 * TwoSR )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );
                                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram20(TensorT * denT, TensorL * denL, TensorF0 * denF0, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF0->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIj );
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi );
                  
                  if (( dimRup > 0 ) && ( dimLdown > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup    =  denT->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIi  );
                     double * Tdown  =  denT->gStorage( NL-1, TwoSLprime, ILxIj, NL+1, TwoSLprime, ILxIj  );
                     double * Fblock = denF0->gStorage( NL+1, TwoSLprime, ILxIi, NL+1, TwoSLprime, ILxIj  );
                     double * Lblock =  denL->gStorage( NL-1, TwoSLprime, ILxIj, NL,   TwoSL,      IL     );
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Fblock, &dimRup, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total += sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) )
                            * Special::phase( TwoSL + 1 - TwoSLprime )
                            * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::ThreeDM::diagram21(TensorT * denT, TensorL * denL, TensorF1 * denF1, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   assert( denF1->get_irrep() == Irreps::directProd( book->gIrrep( orb_i ), denL->get_irrep() ) );
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIi = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxIj = Irreps::directProd( IL, denL->get_irrep()     );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
         
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIj );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIj );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tdown  = denT->gStorage( NL-1, TwoSLprime, ILxIj, NL+1, TwoSLprime, ILxIj );
                     double * Lblock = denL->gStorage( NL-1, TwoSLprime, ILxIj, NL,   TwoSL,      IL    );
                  
                     for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                        int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
                        
                        if ( dimRup > 0 ){
                     
                           double * Tup    =  denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR,      ILxIi );
                           double * Fblock = denF1->gStorage( NL+1, TwoSR, ILxIi, NL+1, TwoSLprime, ILxIj );
                           
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Fblock, &dimRup, &beta, workmem2, &dimLdown );
                           
                           int length = dimLdown * dimRdown;
                           int inc = 1;
                           total += sqrt( 3.0 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                  * Special::phase( TwoSR + TwoSLprime )
                                  * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );
                                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;

}

void CheMPS2::ThreeDM::fill_a_S0( TensorT * denT, Tensor3RDM * tofill, TensorS0 * denS0, double * workmem ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS0->get_irrep(), book->gIrrep( orb_i ) );
   assert( tofill->get_irrep()  == ImxInxIi );
   assert( tofill->get_nelec()  == 3 );
   assert( tofill->get_two_j2() == 1 );
   
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denS0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-3, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL,   TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-2, TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL,   TwoSL,      IL       );
                     double * Tdown   =   denT->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL-2, TwoSL,      ILxImxIn );
                     double * Sblock  =  denS0->gStorage( NL-2, TwoSL,      ILxImxIn,    NL,   TwoSL,      IL       );
                     double * Wblock  = tofill->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL,   TwoSL,      IL       );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = -0.5 * ( TwoSL + 1 );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Sblock, &dimRdown, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               
                  dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi       );
                  dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSLprime, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxIi       );
                     double * Tdown   =   denT->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL-1, TwoSLprime, ILxImxInxIi );
                     double * Sblock  =  denS0->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxIi       );
                     double * Wblock  = tofill->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL,   TwoSL,      IL          );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSL + 1 - TwoSLprime );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Sblock, &dimRdown, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_bcd_S0( TensorT * denT, Tensor3RDM * tofill, TensorS0 * denS0, double * workmem ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS0->get_irrep(), book->gIrrep( orb_i ) );
   assert( tofill->get_irrep()  == ImxInxIi );
   assert( tofill->get_nelec()  == 1 );
   assert( tofill->get_two_j2() == 1 );
   
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denS0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL+1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL,   TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){ // Type 100, 102, 110, 112, 120, 124
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL,   TwoSL,      IL          );
                     double * Tdown   =   denT->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+2, TwoSL,      ILxImxIn    );
                     double * Sblock  =  denS0->gStorage( NL,   TwoSL,      IL,          NL+2, TwoSL,      ILxImxIn    );
                     double * Wblock  = tofill->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxImxInxIi );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSL + 1 - TwoSLprime );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimRup, &alpha, Tup, &dimLup, Sblock, &dimRup, &beta, workmem, &dimLup );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLup, &dimLdown, &dimRdown, &alpha, workmem, &dimLup, Tdown, &dimLdown, &beta, Wblock, &dimLup );
                  
                  }
               
                  dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi       );
                  dimRdown = book->gCurrentDim( orb_i+1, NL+3, TwoSLprime, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){ // Type 105, 107, 115, 117, 125, 129
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxIi       );
                     double * Tdown   =   denT->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+3, TwoSLprime, ILxImxInxIi );
                     double * Sblock  =  denS0->gStorage( NL+1, TwoSLprime, ILxIi ,      NL+3, TwoSLprime, ILxImxInxIi );
                     double * Wblock  = tofill->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxImxInxIi );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = - 0.5 * ( TwoSLprime + 1 );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimRup, &alpha, Tup, &dimLup, Sblock, &dimRup, &beta, workmem, &dimLup );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLup, &dimLdown, &dimRdown, &alpha, workmem, &dimLup, Tdown, &dimLdown, &beta, Wblock, &dimLup );
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_a_S1( TensorT * denT, Tensor3RDM * doublet, Tensor3RDM * quartet, TensorS1 * denS1, double * workmem, double * workmem2 ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS1->get_irrep(), book->gIrrep( orb_i ) );
   assert( doublet->get_irrep() == ImxInxIi ); assert( doublet->get_nelec() == 3 ); assert( doublet->get_two_j2() == 1 );
   assert( quartet->get_irrep() == ImxInxIi ); assert( quartet->get_nelec() == 3 ); assert( quartet->get_two_j2() == 3 );
   
   doublet->clear();
   quartet->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denS1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-3; TwoSLprime <= TwoSL+3; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-3, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
               
                     int dimRup   = book->gCurrentDim( orb_i+1, NL,   TwoSL,      IL       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL-2, TwoSRprime, ILxImxIn );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSL - TwoSRprime ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL,   TwoSL,      IL       );
                        double * Tdown   =   denT->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL-2, TwoSRprime, ILxImxIn );
                        double * Sblock  =  denS1->gStorage( NL-2, TwoSRprime, ILxImxIn,    NL,   TwoSL,      IL       );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Sblock, &dimRdown, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSRprime + 1 ) ) * ( TwoSL + 1 )
                                            * Special::phase( TwoSL + TwoSLprime + 1 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = 2 * sqrt( TwoSRprime + 1.0 ) * ( TwoSL + 1 )
                                            * Special::phase( TwoSL + TwoSLprime + 3 )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSLprime, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                        double * Tdown   =   denT->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL-1, TwoSLprime, ILxImxInxIi );
                        double * Sblock  =  denS1->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSR,      ILxIi       );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Sblock, &dimRdown, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                            * Special::phase( TwoSR + TwoSLprime )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-3, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = 2 * sqrt( TwoSL + 1.0 ) * ( TwoSR + 1 )
                                            * Special::phase( TwoSR + TwoSLprime )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_bcd_S1( TensorT * denT, Tensor3RDM * doublet, Tensor3RDM * quartet, TensorS1 * denS1, double * workmem, double * workmem2 ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS1->get_irrep(), book->gIrrep( orb_i ) );
   assert( doublet->get_irrep() == ImxInxIi ); assert( doublet->get_nelec() == 1 ); assert( doublet->get_two_j2() == 1 );
   assert( quartet->get_irrep() == ImxInxIi ); assert( quartet->get_nelec() == 1 ); assert( quartet->get_two_j2() == 3 );
   
   doublet->clear();
   quartet->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denS1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-3; TwoSLprime <= TwoSL+3; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL+1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
               
                     int dimRup   = book->gCurrentDim( orb_i+1, NL,   TwoSL,      IL       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSRprime, ILxImxIn );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSL - TwoSRprime ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL,   TwoSL,      IL       );
                        double * Tdown   =   denT->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+2, TwoSRprime, ILxImxIn );
                        double * Sblock  =  denS1->gStorage( NL,   TwoSL,      IL,          NL+2, TwoSRprime, ILxImxIn );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimRup,   &alpha, Tup,     &dimLup, Sblock, &dimRup,   &beta, workmem,  &dimLup );
                        dgemm_( &notrans, &trans,   &dimLup, &dimLdown, &dimRdown, &alpha, workmem, &dimLup, Tdown,  &dimLdown, &beta, workmem2, &dimLup );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( 0.5 * ( TwoSLprime + 1 ) ) * ( TwoSRprime + 1 )
                                            * Special::phase( TwoSL + TwoSRprime )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( TwoSLprime + 1.0 ) * ( TwoSRprime + 1 )
                                            * Special::phase( TwoSL + TwoSRprime )
                                            * Wigner::wigner6j( 1, 2, 3, TwoSL, TwoSLprime, TwoSRprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL+3, TwoSLprime, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                        double * Tdown   =   denT->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+3, TwoSLprime, ILxImxInxIi );
                        double * Sblock  =  denS1->gStorage( NL+1, TwoSR,      ILxIi ,      NL+3, TwoSLprime, ILxImxInxIi );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimRup,   &alpha, Tup,     &dimLup, Sblock, &dimRup,   &beta, workmem,  &dimLup );
                        dgemm_( &notrans, &trans,   &dimLup, &dimLdown, &dimRdown, &alpha, workmem, &dimLup, Tdown,  &dimLdown, &beta, workmem2, &dimLup );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( 0.5 * ( TwoSR + 1 ) ) * ( TwoSLprime + 1 )
                                            * Special::phase( TwoSL + TwoSLprime + 3 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( TwoSR + 1.0 ) * ( TwoSLprime + 1 )
                                            * Special::phase( TwoSL + TwoSLprime + 1 )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_F0_T( TensorT * denT, Tensor3RDM * tofill, TensorF0 * denF0, double * workmem ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denF0->get_irrep(), book->gIrrep( orb_i ) );
   assert( tofill->get_irrep()  == ImxInxIi );
   assert( tofill->get_nelec()  == 1 );
   assert( tofill->get_two_j2() == 1 );
   
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denF0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL, TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL, TwoSL,      IL       );
                     double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL,      ILxImxIn );
                     double * Fblock  =  denF0->gStorage( NL,   TwoSL,      ILxImxIn,    NL, TwoSL,      IL       );
                     double * Wblock  = tofill->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL,      IL       );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * ( TwoSL + 1 );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Fblock, &dimRdown, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               
                  dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi       );
                  dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxIi       );
                     double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxImxInxIi );
                     double * Fblock  =  denF0->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxIi       );
                     double * Wblock  = tofill->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL,   TwoSL,      IL          );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSLprime + 1 - TwoSL );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Fblock, &dimRdown, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_F1_T( TensorT * denT, Tensor3RDM * doublet, Tensor3RDM * quartet, TensorF1 * denF1, double * workmem, double * workmem2 ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denF1->get_irrep(), book->gIrrep( orb_i ) );
   assert( doublet->get_irrep() == ImxInxIi ); assert( doublet->get_nelec() == 1 ); assert( doublet->get_two_j2() == 1 );
   assert( quartet->get_irrep() == ImxInxIi ); assert( quartet->get_nelec() == 1 ); assert( quartet->get_two_j2() == 3 );
   
   doublet->clear();
   quartet->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denF1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-3; TwoSLprime <= TwoSL+3; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
               
                     int dimRup   = book->gCurrentDim( orb_i+1, NL, TwoSL,      IL       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSRprime, ILxImxIn );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSL - TwoSRprime ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL, TwoSL,      IL       );
                        double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSRprime, ILxImxIn );
                        double * Fblock  =  denF1->gStorage( NL,   TwoSRprime, ILxImxIn,    NL, TwoSL,      IL       );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Fblock, &dimRdown, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSL + 1 ) ) * ( TwoSRprime + 1 )
                                            * Special::phase( TwoSLprime + TwoSRprime + 3 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( TwoSL + 1.0 ) * ( TwoSRprime + 1 )
                                            * Special::phase( TwoSLprime + TwoSRprime + 3 )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                        double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxImxInxIi );
                        double * Fblock  =  denF1->gStorage( NL+1, TwoSLprime, ILxImxInxIi, NL+1, TwoSR,      ILxIi       );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Fblock, &dimRdown, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) * ( TwoSR + 1 ) )
                                            * Special::phase( 2*TwoSLprime + 2 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) * ( TwoSR + 1 ) )
                                            * Special::phase( 2*TwoSLprime )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_F0( TensorT * denT, Tensor3RDM * tofill, TensorF0 * denF0, double * workmem ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denF0->get_irrep(), book->gIrrep( orb_i ) );
   assert( tofill->get_irrep()  == ImxInxIi );
   assert( tofill->get_nelec()  == 1 );
   assert( tofill->get_two_j2() == 1 );
   
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denF0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL, TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL, TwoSL,      IL       );
                     double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL,      ILxImxIn );
                     double * Fblock  =  denF0->gStorage( NL,   TwoSL,      IL,          NL, TwoSL,      ILxImxIn );
                     double * Wblock  = tofill->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL,      IL       );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * ( TwoSL + 1 );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Fblock, &dimRup, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               
                  dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIi       );
                  dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSLprime, ILxIi       );
                     double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxImxInxIi );
                     double * Fblock  =  denF0->gStorage( NL+1, TwoSLprime, ILxIi,       NL+1, TwoSLprime, ILxImxInxIi );
                     double * Wblock  = tofill->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL,   TwoSL,      IL          );
                  
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 0.5 * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSLprime + 1 - TwoSL );
                     double beta  = 0.0; //SET
                     dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown, &dimLdown, Fblock, &dimRup, &beta, workmem, &dimLdown );
                     alpha        = 1.0;
                     beta         = 1.0; //ADD
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup, &alpha, workmem, &dimLdown, Tup, &dimLup, &beta, Wblock, &dimLdown );
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_F1( TensorT * denT, Tensor3RDM * doublet, Tensor3RDM * quartet, TensorF1 * denF1, double * workmem, double * workmem2 ) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denF1->get_irrep(), book->gIrrep( orb_i ) );
   assert( doublet->get_irrep() == ImxInxIi ); assert( doublet->get_nelec() == 1 ); assert( doublet->get_two_j2() == 1 );
   assert( quartet->get_irrep() == ImxInxIi ); assert( quartet->get_nelec() == 1 ); assert( quartet->get_two_j2() == 3 );
   
   doublet->clear();
   quartet->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxInxIi = Irreps::directProd( IL, ImxInxIi              );
            const int ILxImxIn    = Irreps::directProd( IL, denF1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            
            for ( int TwoSLprime = TwoSL-3; TwoSLprime <= TwoSL+3; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL          );
               int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLprime, ILxImxInxIi );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
               
                     int dimRup   = book->gCurrentDim( orb_i+1, NL, TwoSL,      IL       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSRprime, ILxImxIn );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSL - TwoSRprime ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL, TwoSL,      IL       );
                        double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSRprime, ILxImxIn );
                        double * Fblock  =  denF1->gStorage( NL,   TwoSL,      IL,          NL, TwoSRprime, ILxImxIn );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Fblock, &dimRup, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup, &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSRprime + 1 ) ) * ( TwoSL + 1 )
                                            * Special::phase( TwoSLprime + TwoSL + 3 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( TwoSRprime + 1.0 ) * ( TwoSL + 1 )
                                            * Special::phase( TwoSLprime + TwoSL + 3 )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) <= 2 )){
                     
                        double * Tup     =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                        double * Tdown   =   denT->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL+1, TwoSLprime, ILxImxInxIi );
                        double * Fblock  =  denF1->gStorage( NL+1, TwoSR,      ILxIi,       NL+1, TwoSLprime, ILxImxInxIi );
                     
                        char notrans = 'N';
                        char trans   = 'T';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET
                        dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Fblock, &dimRup, &beta, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup, &beta, workmem2, &dimLdown );
                        
                        if ( abs( TwoSL - TwoSLprime ) == 1 ){
                        
                           double * Wblock  = doublet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( 0.5 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                            * Special::phase( TwoSLprime + TwoSR + 2 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL-1, TwoSLprime, ILxImxInxIi, NL, TwoSL, IL );
                           double prefactor = sqrt( TwoSL + 1.0 ) * ( TwoSR + 1 )
                                            * Special::phase( TwoSLprime + TwoSR )
                                            * Wigner::wigner6j( 1, 3, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                     
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_29_33(TensorT * denT, TensorF0 * tofill, TensorF0 * denF0, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denF0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denF0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denF0->get_irrep() );
            
            int dimLup   = book->gCurrentDim( orb_i, NL, TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i, NL, TwoSL, ILxImxIn );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 )){
            
               {
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL, IL,       NL+2, TwoSL, IL       );
                     double * Tdown =   denT->gStorage( NL,   TwoSL, ILxImxIn, NL+2, TwoSL, ILxImxIn );
                     double * right =  denF0->gStorage( NL+2, TwoSL, ILxImxIn, NL+2, TwoSL, IL       );
                     double * left  = tofill->gStorage( NL,   TwoSL, ILxImxIn, NL,   TwoSL, IL       );
                  
                     double factor = TwoSL + 1.0;
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
                     
                  }
               }
            
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL, IL,          NL+1, TwoSR, ILxIi       );
                     double * Tdown =   denT->gStorage( NL,   TwoSL, ILxImxIn,    NL+1, TwoSR, ILxImxInxIi );
                     double * right =  denF0->gStorage( NL+1, TwoSR, ILxImxInxIi, NL+1, TwoSR, ILxIi       );
                     double * left  = tofill->gStorage( NL,   TwoSL, ILxImxIn,    NL,   TwoSL, IL          );
                  
                     double factor = 0.5 * ( TwoSR + 1 );
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );

                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_30_32(TensorT * denT, TensorF1 * tofill, TensorF1 * denF1, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denF1->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denF1->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLprime, ILxImxIn );
               int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL,      IL       );
               int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                  double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,       NL+2, TwoSL,      IL       );
                  double * Tdown =   denT->gStorage( NL,   TwoSLprime, ILxImxIn, NL+2, TwoSLprime, ILxImxIn );
                  double * right =  denF1->gStorage( NL+2, TwoSLprime, ILxImxIn, NL+2, TwoSL,      IL       );
                  double * left  = tofill->gStorage( NL,   TwoSLprime, ILxImxIn, NL,   TwoSL,      IL       );
               
                  double factor = sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSL - TwoSLprime );
                  char notrans  = 'N';
                  char trans    = 'T';
                  double zero   = 0.0;
                  double one    = 1.0;
                  dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                  dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
               
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_36_42(TensorT * denT, TensorF1 * tofill, TensorF0 * denF0, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denF0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denF0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denF0->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL, TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i, NL, TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) == 1 )){
                     
                        double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR, ILxIi       );
                        double * Tdown =   denT->gStorage( NL,   TwoSLprime, ILxImxIn,    NL+1, TwoSR, ILxImxInxIi );
                        double * right =  denF0->gStorage( NL+1, TwoSR,      ILxImxInxIi, NL+1, TwoSR, ILxIi       );
                        double * left  = tofill->gStorage( NL,   TwoSLprime, ILxImxIn,    NL,   TwoSL, IL          );
                     
                        double factor = 0.5 * sqrt( 6.0 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                      * Special::phase( TwoSR + TwoSLprime + 1 )
                                      * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR );
                        char notrans  = 'N';
                        char trans    = 'T';
                        double zero   = 0.0;
                        double one    = 1.0;
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );

                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_34_35_37_38(TensorT * denT, TensorF1 * fill34, TensorF0 * fill35, TensorF1 * fill37, TensorF1 * fill38, TensorF1 * denF1, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   fill34->clear();
   fill35->clear();
   fill37->clear();
   fill38->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denF1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denF1->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL, TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i, NL, TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                     for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
                  
                        int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                        int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSRprime, ILxImxInxIi );
                        
                        if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSR - TwoSRprime ) <= 2 )){
                        
                           double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                           double * Tdown =   denT->gStorage( NL,   TwoSLprime, ILxImxIn,    NL+1, TwoSRprime, ILxImxInxIi );
                           double * right =  denF1->gStorage( NL+1, TwoSRprime, ILxImxInxIi, NL+1, TwoSR,      ILxIi       );
                           
                           char notrans  = 'N';
                           char trans    = 'T';
                           double zero   = 0.0;
                           double one    = 1.0;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &one, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one, workmem, &dimLdown, Tup,   &dimLup,   &zero, workmem2, &dimLdown );
                           
                           {
                              double * left = fill34->gStorage( NL, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = 0.5 * ( TwoSRprime + 1 ) * sqrt( 1.0 * ( TwoSR + 1 ) * ( TwoSL + 1 ) )
                                            * Special::phase( TwoSL + TwoSR + 3 )
                                            * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRprime, TwoSLprime, 2 );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           if ( TwoSL == TwoSLprime ){
                              double * left = fill35->gStorage( NL, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = 0.5 * ( TwoSRprime + 1 ) * sqrt( 6.0 * ( TwoSR + 1 ) )
                                            * Special::phase( TwoSL + TwoSRprime + 3 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSRprime, TwoSL );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double * left = fill37->gStorage( NL, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = 3 * ( TwoSRprime + 1 ) * sqrt( 1.0 * ( TwoSR + 1 ) * ( TwoSL + 1 ) )
                                            * Special::phase( 2*TwoSL + TwoSR + TwoSRprime )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR      )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSRprime, TwoSLprime );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double * left = fill38->gStorage( NL, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = 3 * ( TwoSRprime + 1 ) * sqrt( 1.0 * ( TwoSR + 1 ) * ( TwoSL + 1 ) )
                                            * Special::phase( 2*TwoSR + TwoSL + TwoSLprime )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSRprime )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSRprime, TwoSL      );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_49_51(TensorT * denT, TensorF0 * tofill, TensorS0 * denS0, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denS0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denS0->get_irrep() );
            
            int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSL, ILxImxIn );
            int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL       );
            int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSL, ILxImxIn );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRup > 0 ) && ( dimRdown > 0 )){
            
               double * Tup   =   denT->gStorage( NL, TwoSL, IL,       NL+2, TwoSL, IL       );
               double * Tdown =   denT->gStorage( NL, TwoSL, ILxImxIn, NL,   TwoSL, ILxImxIn );
               double * right =  denS0->gStorage( NL, TwoSL, ILxImxIn, NL+2, TwoSL, IL       );
               double * left  = tofill->gStorage( NL, TwoSL, ILxImxIn, NL,   TwoSL, IL       );
            
               double factor = - ( TwoSL + 1.0 );
               char notrans  = 'N';
               char trans    = 'T';
               double zero   = 0.0;
               double one    = 1.0;
               dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
               dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
               
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_50_52(TensorT * denT, TensorF1 * tofill, TensorS1 * denS1, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denS1->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denS1->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLprime, ILxImxIn );
               int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL,      IL       );
               int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRup > 0 ) && ( dimRdown > 0 )){
               
                  double * Tup   =   denT->gStorage( NL, TwoSL,      IL,       NL+2, TwoSL,      IL       );
                  double * Tdown =   denT->gStorage( NL, TwoSLprime, ILxImxIn, NL,   TwoSLprime, ILxImxIn );
                  double * right =  denS1->gStorage( NL, TwoSLprime, ILxImxIn, NL+2, TwoSL,      IL       );
                  double * left  = tofill->gStorage( NL, TwoSLprime, ILxImxIn, NL,   TwoSL,      IL       );
               
                  double factor = - ( TwoSL + 1.0 );
                  char notrans  = 'N';
                  char trans    = 'T';
                  double zero   = 0.0;
                  double one    = 1.0;
                  dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                  dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
                  
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_22_24(TensorT * denT, TensorS0 * tofill, TensorS0 * denS0, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denS0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denS0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denS0->get_irrep() );
            
            int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i, NL-2, TwoSL, ILxImxIn );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 )){
            
               {
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSL, ILxImxIn );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL, IL,       NL+2, TwoSL, IL       );
                     double * Tdown =   denT->gStorage( NL-2, TwoSL, ILxImxIn, NL,   TwoSL, ILxImxIn );
                     double * right =  denS0->gStorage( NL,   TwoSL, ILxImxIn, NL+2, TwoSL, IL       );
                     double * left  = tofill->gStorage( NL-2, TwoSL, ILxImxIn, NL,   TwoSL, IL       );
                  
                     double factor = TwoSL + 1.0;
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
                     
                  }
               }
            
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
               
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi       );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSR, ILxImxInxIi );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL, IL,          NL+1, TwoSR, ILxIi       );
                     double * Tdown =   denT->gStorage( NL-2, TwoSL, ILxImxIn,    NL-1, TwoSR, ILxImxInxIi );
                     double * right =  denS0->gStorage( NL-1, TwoSR, ILxImxInxIi, NL+1, TwoSR, ILxIi       );
                     double * left  = tofill->gStorage( NL-2, TwoSL, ILxImxIn,    NL,   TwoSL, IL          );
                  
                     double factor = 0.5 * ( TwoSR + 1 );
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );

                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_28(TensorT * denT, TensorS1 * tofill, TensorS0 * denS0, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denS0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denS0->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denS0->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i, NL-2, TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                     int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi       );
                     int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSR, ILxImxInxIi );
                     
                     if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSLprime - TwoSR ) == 1 )){
                     
                        double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR, ILxIi       );
                        double * Tdown =   denT->gStorage( NL-2, TwoSLprime, ILxImxIn,    NL-1, TwoSR, ILxImxInxIi );
                        double * right =  denS0->gStorage( NL-1, TwoSR,      ILxImxInxIi, NL+1, TwoSR, ILxIi       );
                        double * left  = tofill->gStorage( NL-2, TwoSLprime, ILxImxIn,    NL,   TwoSL, IL          );
                     
                        double factor = sqrt( 1.5 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                      * Special::phase( TwoSLprime + TwoSR + 1 )
                                      * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR );
                        char notrans  = 'N';
                        char trans    = 'T';
                        double zero   = 0.0;
                        double one    = 1.0;
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );

                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_23(TensorT * denT, TensorS1 * tofill, TensorS1 * denS1, double * workmem) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denS1->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denS1->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSLprime, ILxImxIn );
               int dimRup   = book->gCurrentDim( orb_i+1, NL+2, TwoSL,      IL       );
               int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRup > 0 ) && ( dimRdown > 0 )){
               
                  double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,       NL+2, TwoSL,      IL       );
                  double * Tdown =   denT->gStorage( NL-2, TwoSLprime, ILxImxIn, NL,   TwoSLprime, ILxImxIn );
                  double * right =  denS1->gStorage( NL,   TwoSLprime, ILxImxIn, NL+2, TwoSL,      IL       );
                  double * left  = tofill->gStorage( NL-2, TwoSLprime, ILxImxIn, NL,   TwoSL,      IL       );
               
                  double factor = TwoSL + 1.0;
                  char notrans  = 'N';
                  char trans    = 'T';
                  double zero   = 0.0;
                  double one    = 1.0;
                  dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                  dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,     &dimLdown );
                        
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_25_26_27(TensorT * denT, TensorS1 * fill25, TensorS1 * fill26, TensorS0 * fill27, TensorS1 * denS1, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   fill25->clear();
   fill26->clear();
   fill27->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn    = Irreps::directProd( IL, denS1->get_irrep()    );
            const int ILxIi       = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxInxIi = Irreps::directProd( ILxIi, denS1->get_irrep() );
            
            for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
            
               int dimLup   = book->gCurrentDim( orb_i, NL,   TwoSL,      IL       );
               int dimLdown = book->gCurrentDim( orb_i, NL-2, TwoSLprime, ILxImxIn );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
                  for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                     for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
                  
                        int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR,      ILxIi       );
                        int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSRprime, ILxImxInxIi );
                        
                        if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSRprime - TwoSR ) <= 2 )){
                        
                           double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,          NL+1, TwoSR,      ILxIi       );
                           double * Tdown =   denT->gStorage( NL-2, TwoSLprime, ILxImxIn,    NL-1, TwoSRprime, ILxImxInxIi );
                           double * right =  denS1->gStorage( NL-1, TwoSRprime, ILxImxInxIi, NL+1, TwoSR,      ILxIi       );
                           
                           char notrans  = 'N';
                           char trans    = 'T';
                           double zero   = 0.0;
                           double one    = 1.0;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &one, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one, workmem, &dimLdown, Tup,   &dimLup,   &zero, workmem2, &dimLdown );
                           
                           {
                              double * left = fill25->gStorage( NL-2, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = ( TwoSR + 1 ) * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) )
                                            * Special::phase( TwoSL + TwoSRprime + 3 )
                                            * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRprime, TwoSLprime, 2 );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double * left = fill26->gStorage( NL-2, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = 3 * ( TwoSR + 1 ) * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) )
                                            * Special::phase( TwoSR + TwoSRprime + TwoSL + TwoSLprime + 2 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSRprime, TwoSL      )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSRprime );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           if ( TwoSLprime == TwoSL ){
                              double * left = fill27->gStorage( NL-2, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                              double factor = ( TwoSR + 1 ) * sqrt( 1.5 * ( TwoSRprime + 1 ) )
                                            * Special::phase( TwoSL + TwoSR + 3 )
                                            * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSRprime, TwoSL );
                              int length = dimLup * dimLdown;
                              int inc    = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_45_47(TensorT * denT, TensorS0 * tofill, TensorF0 * denF0, double * workmem, const bool first) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denF0->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denF0->get_irrep() );
            
            int dimLup   = book->gCurrentDim( orb_i,   NL,   TwoSL, IL       );
            int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSL, ILxImxIn );
            int dimRup   = book->gCurrentDim( orb_i+1, NL,   TwoSL, IL       );
            int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSL, ILxImxIn );
                  
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimRup > 0 ) && ( dimRdown > 0 )){
            
               double * Tup   =   denT->gStorage( NL,   TwoSL, IL,       NL,   TwoSL, IL       );
               double * Tdown =   denT->gStorage( NL-2, TwoSL, ILxImxIn, NL,   TwoSL, ILxImxIn );
               double * left  = tofill->gStorage( NL-2, TwoSL, ILxImxIn, NL,   TwoSL, IL       );
            
               double factor = -( TwoSL + 1.0 );
               char notrans  = 'N';
               char trans    = 'T';
               double zero   = 0.0;
               double one    = 1.0;
               if ( first ){
                  double * right =  denF0->gStorage( NL, TwoSL, ILxImxIn, NL, TwoSL, IL );
                  dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem, &dimLdown );
                  dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,    &dimLdown );
               } else {
                  double * right =  denF0->gStorage( NL, TwoSL, IL, NL, TwoSL, ILxImxIn );
                  dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRup, &zero, workmem, &dimLdown );
                  dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup, &one,  left,    &dimLdown );
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_tens_46_48(TensorT * denT, TensorS1 * tofill, TensorF1 * denF1, double * workmem, const bool first) const{

   const int orb_i = denT->gIndex();
   assert( tofill->get_irrep() == denF1->get_irrep() );
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxImxIn = Irreps::directProd( IL, denF1->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               for ( int TwoSLprime = TwoSL-2; TwoSLprime <= TwoSL+2; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSLprime, ILxImxIn );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSLprime, ILxImxIn );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,       NL,   TwoSL,      IL       );
                     double * Tdown =   denT->gStorage( NL-2, TwoSLprime, ILxImxIn, NL,   TwoSLprime, ILxImxIn );
                     double * left  = tofill->gStorage( NL-2, TwoSLprime, ILxImxIn, NL,   TwoSL,      IL       );
                  
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     if ( first ){
                        double factor  = sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSL - TwoSLprime + 2 );
                        double * right =  denF1->gStorage( NL, TwoSLprime, ILxImxIn, NL, TwoSL, IL );
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem, &dimLdown );
                        dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,    &dimLdown );
                     } else {
                        double factor  = - ( TwoSL + 1.0 );
                        double * right =  denF1->gStorage( NL, TwoSL, IL, NL, TwoSLprime, ILxImxIn );
                        dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRup, &zero, workmem, &dimLdown );
                        dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup, &one,  left,    &dimLdown );
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_53_54(TensorT * denT, Tensor3RDM * tofill, TensorL * denL, double * workmem) const{

   const int orb_i = denT->gIndex();
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIm = Irreps::directProd( IL, denL->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-3, TwoSLprime, ILxIm );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSLprime, ILxIm );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,    NL,   TwoSL,      IL    );
                     double * Tdown =   denT->gStorage( NL-3, TwoSLprime, ILxIm, NL-1, TwoSLprime, ILxIm );
                     double * left  = tofill->gStorage( NL-3, TwoSLprime, ILxIm, NL,   TwoSL,      IL    );
                     double * right =   denL->gStorage( NL-1, TwoSLprime, ILxIm, NL,   TwoSL,      IL    );
                  
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     double factor = - sqrt( 0.5 ) * ( TwoSL + 1 );
                     
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem, &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,    &dimLdown );

                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_55_to_60(TensorT * denT, Tensor3RDM * tofill, TensorL * denL, double * workmem) const{

   const int orb_i = denT->gIndex();
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIm = Irreps::directProd( IL, denL->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIm );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIm );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,    NL,   TwoSL,      IL    );
                     double * Tdown =   denT->gStorage( NL-1, TwoSLprime, ILxIm, NL+1, TwoSLprime, ILxIm );
                     double * left  = tofill->gStorage( NL-1, TwoSLprime, ILxIm, NL,   TwoSL,      IL    );
                     double * right =   denL->gStorage( NL,   TwoSL,      IL,    NL+1, TwoSLprime, ILxIm );
                  
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     double factor = sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * Special::phase( TwoSL + 1 - TwoSLprime );
                     
                     dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRup, &zero, workmem, &dimLdown );
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup, &one,  left,    &dimLdown );

                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_61(TensorT * denT, Tensor3RDM * tofill, TensorL * denL, double * workmem) const{

   const int orb_i = denT->gIndex();
   tofill->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIm = Irreps::directProd( IL, denL->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
               
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIm );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLprime, ILxIm );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup   =   denT->gStorage( NL,   TwoSL,      IL,    NL+2, TwoSL,      IL    );
                     double * Tdown =   denT->gStorage( NL-1, TwoSLprime, ILxIm, NL+1, TwoSLprime, ILxIm );
                     double * left  = tofill->gStorage( NL-1, TwoSLprime, ILxIm, NL,   TwoSL,      IL    );
                     double * right =   denL->gStorage( NL+1, TwoSLprime, ILxIm, NL+2, TwoSL,      IL    );
                  
                     char notrans  = 'N';
                     char trans    = 'T';
                     double zero   = 0.0;
                     double one    = 1.0;
                     double factor = sqrt( 0.5 ) * ( TwoSL + 1 );
                     
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &factor, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem, &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one,    workmem, &dimLdown, Tup,   &dimLup,   &one,  left,    &dimLdown );

                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_63_65(TensorT * denT, Tensor3RDM * fill63, Tensor3RDM * fill65, Tensor3RDM * fill67, Tensor3RDM * fill68, Tensor3RDM * fill76, Tensor3RDM * fill77, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   fill63->clear();
   fill65->clear();
   fill67->clear();
   fill68->clear();
   fill76->clear();
   fill77->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIm    = Irreps::directProd( IL, denL->get_irrep()     );
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxIi = Irreps::directProd( ILxIi, denL->get_irrep()  );
            
            for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
            
               int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL    );
               int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
               
               if (( dimLup > 0 ) && ( dimRup > 0 )){
               
                  for ( int TwoSLprime = TwoSL-1; TwoSLprime <= TwoSL+1; TwoSLprime+=2 ){
                     for ( int TwoSRprime = TwoSLprime-1; TwoSRprime <= TwoSLprime+1; TwoSRprime+=2 ){
                  
                        int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIm    );
                        int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSRprime, ILxImxIi );
                        
                        if (( dimLdown > 0 ) && ( dimRdown > 0 ) && ( abs( TwoSR - TwoSRprime ) == 1 )){
                        
                           double * Tup   = denT->gStorage( NL,   TwoSL,      IL,       NL+1, TwoSR,      ILxIi    );
                           double * Tdown = denT->gStorage( NL-1, TwoSLprime, ILxIm,    NL,   TwoSRprime, ILxImxIi );
                           double * right = denL->gStorage( NL,   TwoSRprime, ILxImxIi, NL+1, TwoSR,      ILxIi    );
                        
                           char notrans  = 'N';
                           char trans    = 'T';
                           double zero   = 0.0;
                           double one    = 1.0;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &one, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one, workmem, &dimLdown, Tup,   &dimLup,   &zero, workmem2, &dimLdown );
                           
                           {
                              double factor  = sqrt( 0.5 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) ) * ( TwoSR + 1 )
                                             * Special::phase( TwoSL + TwoSRprime + 2 )
                                             * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRprime, TwoSLprime, 1 );
                              double * left  = fill63->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           if ( TwoSRprime == TwoSL ){
                              double factor  = - sqrt( 0.5 ) * ( TwoSR + 1 );
                              double * left  = fill65->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           if ( TwoSR == TwoSLprime ){
                              double factor  = sqrt( 0.5 * ( TwoSRprime + 1 ) * ( TwoSL + 1 ) ) * Special::phase( TwoSL - TwoSRprime );
                              double * left  = fill67->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double factor  = sqrt( 6.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) ) * ( TwoSR + 1 )
                                             * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSLprime, TwoSRprime )
                                             * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSLprime, TwoSL      );
                              double * left  = fill68->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double factor  = sqrt( 6.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) ) * ( TwoSR + 1 )
                                             * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime )
                                             * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSR      )
                                             * Special::phase( TwoSLprime + TwoSRprime + 2 - TwoSL - TwoSR );
                              double * left  = fill76->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double factor  = - sqrt( 6.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) ) * ( TwoSR + 1 )
                                             * Wigner::wigner9j( 1, 1, 2, TwoSLprime, TwoSL, 1, TwoSRprime, TwoSR, 1 );
                              double * left  = fill77->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }

                        }
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::ThreeDM::fill_69_78_79(TensorT * denT, Tensor3RDM * fill69, Tensor3RDM * fill78, Tensor3RDM * fill79, TensorL * denL, double * workmem, double * workmem2) const{

   const int orb_i = denT->gIndex();
   fill69->clear();
   fill78->clear();
   fill79->clear();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILxIm    = Irreps::directProd( IL, denL->get_irrep()     );
            const int ILxIi    = Irreps::directProd( IL, book->gIrrep( orb_i ) );
            const int ILxImxIi = Irreps::directProd( ILxIi, denL->get_irrep()  );
            
            for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
            
               int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL    );
               int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSR, ILxIi );
               
               if (( dimLup > 0 ) && ( dimRup > 0 )){
               
                  for ( int TwoSRprime = TwoSR-1; TwoSRprime <= TwoSR+1; TwoSRprime+=2 ){
                     for ( int TwoSLprime = TwoSRprime-1; TwoSLprime <= TwoSRprime+1; TwoSLprime+=2 ){
                  
                        int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLprime, ILxIm    );
                        int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSRprime, ILxImxIi );
                        
                        if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                        
                           double * Tup   = denT->gStorage( NL,   TwoSL,      IL,       NL+1, TwoSR,      ILxIi    );
                           double * Tdown = denT->gStorage( NL-1, TwoSLprime, ILxIm,    NL,   TwoSRprime, ILxImxIi );
                           double * right = denL->gStorage( NL,   TwoSRprime, ILxImxIi, NL+1, TwoSR,      ILxIi    );
                        
                           char notrans  = 'N';
                           char trans    = 'T';
                           double zero   = 0.0;
                           double one    = 1.0;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &one, Tdown,   &dimLdown, right, &dimRdown, &zero, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &one, workmem, &dimLdown, Tup,   &dimLup,   &zero, workmem2, &dimLdown );
                           
                           const double prefactor = 2 * sqrt( 3.0 * ( TwoSL + 1 ) * ( TwoSRprime + 1 ) ) * ( TwoSR + 1 );
                           {
                              double factor  = prefactor * Wigner::wigner6j( 1, 1, 2, TwoSR, TwoSLprime, TwoSRprime )
                                                         * Wigner::wigner6j( 1, 2, 3, TwoSLprime, TwoSL, TwoSR );
                              double * left  = fill69->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double factor  = prefactor * Wigner::wigner6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSR )
                                                         * Wigner::wigner6j( 1, 3, 2, TwoSL, TwoSRprime, TwoSLprime )
                                                         * Special::phase( TwoSR + TwoSRprime + TwoSL + TwoSLprime + 2 );
                              double * left  = fill78->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }
                           {
                              double factor  = - prefactor * Wigner::wigner9j( 1, 1, 2, TwoSLprime, TwoSL, 3, TwoSRprime, TwoSR, 1 );
                              double * left  = fill79->gStorage( NL-1, TwoSLprime, ILxIm, NL, TwoSL, IL );
                              int length     = dimLup * dimLdown;
                              int inc        = 1;
                              daxpy_( &length, &factor, workmem2, &inc, left, &inc );
                           }

                        }
                     }
                  }
               }
            }
         }
      }
   }

}



