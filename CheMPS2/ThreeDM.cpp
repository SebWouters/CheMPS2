/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2015 Sebastian Wouters

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
#include <iostream>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <gsl/gsl_sf_coupling.h>

#include "ThreeDM.h"
#include "Lapack.h"
#include "MyHDF5.h"
#include "Options.h"
#include "MPIchemps2.h"

using std::max;
using std::cout;
using std::endl;

CheMPS2::ThreeDM::ThreeDM(const SyBookkeeper * book_in, const Problem * prob_in){

   book = book_in;
   prob = prob_in;
   L = book->gL();
   const int max_integer = INT_MAX;
   const long long size  = L*L*L*L*L*L;
   assert( max_integer >= size );
   elements = new double[ size ];
   for (int cnt = 0; cnt < size; cnt++){ elements[ cnt ] = 0.0; }

}

CheMPS2::ThreeDM::~ThreeDM(){

   delete [] elements;

}

#ifdef CHEMPS2_MPI_COMPILATION
void CheMPS2::ThreeDM::mpi_allreduce(){

   const int size = L*L*L*L*L*L; // L^6 fits in integer (tested in constructor)
   double * temp = new double[ size ];
   MPIchemps2::allreduce_array_double( elements, temp, size );
   for (int cnt = 0; cnt < size; cnt++){ elements[ cnt ] = temp[ cnt ]; }
   delete [] temp;

}
#endif

void CheMPS2::ThreeDM::set_dmrg_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6, const double value){

   //Prob assumes you use DMRG orbs...
   //Irrep sanity checks are performed in ThreeDM::fill_site
   elements[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * ( cnt4 + L * ( cnt5 + L * cnt6 ) ) ) ) ] = value;
   elements[ cnt2 + L * ( cnt3 + L * ( cnt1 + L * ( cnt5 + L * ( cnt6 + L * cnt4 ) ) ) ) ] = value;
   elements[ cnt3 + L * ( cnt1 + L * ( cnt2 + L * ( cnt6 + L * ( cnt4 + L * cnt5 ) ) ) ) ] = value;
   elements[ cnt2 + L * ( cnt1 + L * ( cnt3 + L * ( cnt5 + L * ( cnt4 + L * cnt6 ) ) ) ) ] = value;
   elements[ cnt3 + L * ( cnt2 + L * ( cnt1 + L * ( cnt6 + L * ( cnt5 + L * cnt4 ) ) ) ) ] = value;
   elements[ cnt1 + L * ( cnt3 + L * ( cnt2 + L * ( cnt4 + L * ( cnt6 + L * cnt5 ) ) ) ) ] = value;
   
   elements[ cnt4 + L * ( cnt5 + L * ( cnt6 + L * ( cnt1 + L * ( cnt2 + L * cnt3 ) ) ) ) ] = value;
   elements[ cnt5 + L * ( cnt6 + L * ( cnt4 + L * ( cnt2 + L * ( cnt3 + L * cnt1 ) ) ) ) ] = value;
   elements[ cnt6 + L * ( cnt4 + L * ( cnt5 + L * ( cnt3 + L * ( cnt1 + L * cnt2 ) ) ) ) ] = value;
   elements[ cnt5 + L * ( cnt4 + L * ( cnt6 + L * ( cnt2 + L * ( cnt1 + L * cnt3 ) ) ) ) ] = value;
   elements[ cnt6 + L * ( cnt5 + L * ( cnt4 + L * ( cnt3 + L * ( cnt2 + L * cnt1 ) ) ) ) ] = value;
   elements[ cnt4 + L * ( cnt6 + L * ( cnt5 + L * ( cnt1 + L * ( cnt3 + L * cnt2 ) ) ) ) ] = value;

}

double CheMPS2::ThreeDM::get_dmrg_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = prob->gIrrep(cnt1);
   const int irrep2 = prob->gIrrep(cnt2);
   const int irrep3 = prob->gIrrep(cnt3);
   const int irrep4 = prob->gIrrep(cnt4);
   const int irrep5 = prob->gIrrep(cnt5);
   const int irrep6 = prob->gIrrep(cnt6);
   if ( Irreps::directProd(Irreps::directProd(irrep1, irrep2), irrep3) == Irreps::directProd(Irreps::directProd(irrep4, irrep5), irrep6) ){
      return elements[ cnt1 + L * ( cnt2 + L * ( cnt3 + L * ( cnt4 + L * ( cnt5 + L * cnt6 ) ) ) ) ];
   }
   
   return 0.0;

}

double CheMPS2::ThreeDM::get_ham_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( prob->gReorderD2h() ){
      return get_dmrg_index( prob->gf1(cnt1), prob->gf1(cnt2), prob->gf1(cnt3), prob->gf1(cnt4), prob->gf1(cnt5), prob->gf1(cnt6) );
   }
   return get_dmrg_index( cnt1, cnt2, cnt3, cnt4, cnt5, cnt6 );

}

double CheMPS2::ThreeDM::trace() const{

   double value = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         for (int cnt3=0; cnt3<L; cnt3++){
            value += get_dmrg_index( cnt1, cnt2, cnt3, cnt1, cnt2, cnt3 );
         }
      }
   }
   return value;

}

void CheMPS2::ThreeDM::save() const{

   hid_t file_id = H5Fcreate(CheMPS2::THREE_RDM_storagename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      hid_t group_id = H5Gcreate(file_id, "three_rdm", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

         hsize_t dimarray       = L*L*L*L*L*L;
         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "elements", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, elements);

         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::ThreeDM::read(){

   hid_t file_id = H5Fopen(CheMPS2::THREE_RDM_storagename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t group_id = H5Gopen(file_id, "three_rdm", H5P_DEFAULT);

         hid_t dataset_id = H5Dopen(group_id, "elements", H5P_DEFAULT);
         H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, elements);
         H5Dclose(dataset_id);

      H5Gclose(group_id);
   H5Fclose(file_id);
   
   std::cout << "ThreeDM::read : Everything loaded!" << std::endl;

}

int CheMPS2::ThreeDM::trianglefunction(const int k, const int glob){

   int cnt2tilde = 1;
   while(cnt2tilde*(cnt2tilde+1)/2 <= glob){ cnt2tilde++; }
   return k - cnt2tilde;
   
}

void CheMPS2::ThreeDM::tripletrianglefunction(const int global, int * jkl){

   int cnt3 = 0;
   while ( ((cnt3+1)*(cnt3+2)*(cnt3+3))/6 <= global ){ cnt3++; }
   const int globalmin = global - (cnt3*(cnt3+1)*(cnt3+2))/6;
   int cnt2 = 0;
   while ( ((cnt2+1)*(cnt2+2))/2 <= globalmin ){ cnt2++; }
   int cnt1 = globalmin - (cnt2*(cnt2+1))/2;
   jkl[0] = cnt1;
   jkl[1] = cnt2;
   jkl[2] = cnt3;
   
}

void CheMPS2::ThreeDM::fill_site( TensorT * denT, TensorL *** Ltensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors,
                                  Tensor3RDM **** dm3_a_J0_doublet, Tensor3RDM **** dm3_a_J1_doublet, Tensor3RDM **** dm3_a_J1_quartet,
                                  Tensor3RDM **** dm3_b_J0_doublet, Tensor3RDM **** dm3_b_J1_doublet, Tensor3RDM **** dm3_b_J1_quartet,
                                  Tensor3RDM **** dm3_c_J0_doublet, Tensor3RDM **** dm3_c_J1_doublet, Tensor3RDM **** dm3_c_J1_quartet,
                                  Tensor3RDM **** dm3_d_J0_doublet, Tensor3RDM **** dm3_d_J1_doublet, Tensor3RDM **** dm3_d_J1_quartet ){

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   const int orb_i = denT->gIndex();
   const int DIM = max(book->gMaxDimAtBound( orb_i ), book->gMaxDimAtBound( orb_i+1 ));
   const double prefactor_spin = 1.0/(prob->gTwoS() + 1.0);
   
   #pragma omp parallel
   {
   
      double * workmem  = new double[DIM*DIM];
      double * workmem2 = new double[DIM*DIM];
      
      /*for (int orb_j = 0; orb_j < orb_i; orb_j++){
         for (int orb_k = orb_j; orb_k < orb_i; orb_k++){*/
      const int upperbound1 = ( orb_i * ( orb_i+1 )) / 2;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int global = 0; global < upperbound1; global++){
         const int orb_k = orb_i - trianglefunction( orb_i, global ) - 1; // cnt2tilde - 1
         const int orb_j = global - (orb_k*(orb_k+1))/2;
         //cout << orb_j << " <= " << orb_k << " < " << orb_i << endl;
         if ( book->gIrrep( orb_j ) == book->gIrrep( orb_k )){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               const double d1 = diagram1( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], workmem ) * prefactor_spin;
               set_dmrg_index( orb_j, orb_i, orb_i, orb_k, orb_i, orb_i,  4 * d1 );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_k, orb_i, -2 * d1 );
            }
         }
      }
      
      const int upperbound2 = orb_i * ( L - 1 - orb_i );
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int global = 0; global < upperbound2; global++){
         const int orb_j = global % orb_i;
         const int orb_k = orb_i + 1 + ( global / orb_i );
         if ( book->gIrrep( orb_j ) == book->gIrrep( orb_k )){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) ) //Everyone owns the L-tensors --> task division based on ABSigma-tensor ownership
            #endif
            {
               const double d2 = diagram2_22_23( denT, Ltensors[orb_i-1][orb_i-1-orb_j], Ltensors[orb_i][orb_k-1-orb_i], workmem, workmem2 ) * prefactor_spin;
               set_dmrg_index( orb_j, orb_i, orb_i, orb_k, orb_i, orb_i, 2 * d2 );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_k, orb_i,   - d2 );
            }
         }
      }
      
      /*for (int j_index=theindex+1; j_index<L; j_index++){
         for (int k_index=j_index; k_index<L; k_index++){*/
      const int triangle3   = L - orb_i - 1;
      const int upperbound3 = (triangle3*(triangle3+1))/2;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int global = 0; global < upperbound3; global++){
         const int row = trianglefunction( triangle3, global );
         const int col = global - ((triangle3-row)*(triangle3-1-row))/2;
         const int orb_j = orb_i + 1 + row;
         const int orb_k = orb_j + col;
         //cout << orb_i << " <  " << orb_j << " <= " << orb_k << endl;
         if ( book->gIrrep( orb_j ) == book->gIrrep( orb_k )){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               const double d3 = diagram3( denT, F0tensors[orb_i][orb_k-orb_j][orb_j-1-orb_i], workmem ) * prefactor_spin;
               set_dmrg_index( orb_j, orb_i, orb_i, orb_k, orb_i, orb_i,  4 * d3 );
               set_dmrg_index( orb_j, orb_i, orb_i, orb_i, orb_k, orb_i, -2 * d3 );
            }
         }
      }
      
      const int upperbound4 = ( orb_i * ( orb_i+1 ) * ( orb_i + 2 ) ) / 6;
      int jkl[] = { 0, 0, 0 };
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int global = 0; global < upperbound4; global++){
         tripletrianglefunction( global, jkl );
         const int orb_j = jkl[0];
         const int orb_k = jkl[1];
         const int orb_l = jkl[2];
         const int recalculate_global = orb_j + (orb_k*(orb_k+1))/2 + (orb_l*(orb_l+1)*(orb_l+2))/6;
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
                  const double d4 =                        diagram4_5_6_7_8_9( denT, dm3_b_J0_doublet[cnt1][cnt2][cnt3], workmem, 'B' ) * prefactor_spin;
                  const double d5 = (( cnt1 == 0 ) ? 0.0 : diagram4_5_6_7_8_9( denT, dm3_b_J1_doublet[cnt1][cnt2][cnt3], workmem, 'B' ) * prefactor_spin);
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_i, d4 + d5 );
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_i, d4 - d5 );
                  set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_i, orb_l, -2 * d4 );
               }
               if ( cnt1 + cnt2 > 0 ){
                  const double d6 = diagram4_5_6_7_8_9( denT, dm3_c_J0_doublet[cnt1][cnt2][cnt3], workmem, 'C' ) * prefactor_spin;
                  const double d7 = diagram4_5_6_7_8_9( denT, dm3_c_J1_doublet[cnt1][cnt2][cnt3], workmem, 'C' ) * prefactor_spin;
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_i, -2 * d6 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_i, d6 - d7 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_i, orb_k, d6 + d7 );
               }
               {
                  const double d8 = diagram4_5_6_7_8_9( denT, dm3_d_J0_doublet[cnt1][cnt2][cnt3], workmem, 'D' ) * prefactor_spin;
                  const double d9 = diagram4_5_6_7_8_9( denT, dm3_d_J1_doublet[cnt1][cnt2][cnt3], workmem, 'D' ) * prefactor_spin;
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_j, orb_i, orb_i, -2 * d8 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_j, orb_i, d8 + d9 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_i, orb_j, d8 - d9 );
               }
            }
         }
      }
      
      const int upperbound10 = upperbound1 * ( L - orb_i - 1 );
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int combined = 0; combined < upperbound10; combined++){
         const int global = combined % upperbound1;
         const int orb_l  = orb_i + 1 + ( combined / upperbound1 );
         const int orb_k  = orb_i - trianglefunction( orb_i, global ) - 1;
         const int orb_j  = global - (orb_k*(orb_k+1))/2;
         if ( Irreps::directProd(book->gIrrep( orb_j ), book->gIrrep( orb_k )) == Irreps::directProd(book->gIrrep( orb_l ), book->gIrrep( orb_i )) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
            #endif
            {
               const double d10 = diagram10( denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin;
               const double d11 = (( orb_j == orb_k ) ? 0.0 :
                                  diagram11( denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin );
               set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_i, d10 + d11 );
               set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_i, d10 - d11 );
               set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_i, orb_l, - 2 * d10 );
            }
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
            #endif
            {
               {
                  const double d12 = diagram12( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  const double d13 = diagram13( denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_i, - 2 * d12 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_i, d12 - d13 );
                  set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_i, orb_k, d12 + d13 );
               }
               if ( orb_j < orb_k ){
                  const double d14 = diagram14( denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  const double d15 = diagram15( denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], Ltensors[orb_i][orb_l-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_j, orb_i, orb_i, - 2 * d14 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_j, orb_i, d14 - d15 );
                  set_dmrg_index( orb_k, orb_l, orb_i, orb_i, orb_i, orb_j, d14 + d15 );
               }
            }
         }
      }
      
      const int upperbound16 = upperbound3 * orb_i;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int combined = 0; combined < upperbound16; combined++){
         const int orb_j  = combined % orb_i;
         const int global = combined / orb_i;
         const int row    = trianglefunction( triangle3, global );
         const int col    = global - ((triangle3-row)*(triangle3-1-row))/2;
         const int orb_m  = orb_i + 1 + row;
         const int orb_n  = orb_m + col;
         if ( Irreps::directProd(book->gIrrep( orb_j ), book->gIrrep( orb_i )) == Irreps::directProd(book->gIrrep( orb_m ), book->gIrrep( orb_n )) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_absigma( orb_m, orb_n ) )
            #endif
            {
               const double d16 = diagram16( denT, Ltensors[orb_i-1][orb_i-1-orb_j], S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin;
               const double d17 = (( orb_m == orb_n ) ? 0.0 :
                                  diagram17( denT, Ltensors[orb_i-1][orb_i-1-orb_j], S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin );
               set_dmrg_index( orb_m, orb_n, orb_i, orb_j, orb_i, orb_i, d16 - d17 );
               set_dmrg_index( orb_m, orb_n, orb_i, orb_i, orb_j, orb_i, d16 + d17 );
               set_dmrg_index( orb_m, orb_n, orb_i, orb_i, orb_i, orb_j, - 2 * d16 );
            }
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_m, orb_n ) )
            #endif
            {
               {
                  const double d18 = diagram18( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  const double d19 = diagram19( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_i, orb_i, d18 + d19 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_n, orb_i, - 2 * d18 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_i, orb_n, d18 - d19 );
               }
               if ( orb_m < orb_n ){
                  const double d20 = diagram20( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin;
                  const double d21 = diagram21( denT, Ltensors[orb_i-1][orb_i-1-orb_j], F1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ) * prefactor_spin;
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
            
            const int irrep_mn = Irreps::directProd( book->gIrrep( orb_m ), book->gIrrep( orb_n ) );
            
            #ifdef CHEMPS2_MPI_COMPILATION
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
         
            #ifdef CHEMPS2_MPI_COMPILATION
               #pragma omp for schedule(dynamic) // Wait after this for loop --> with MPI certain terms will be deleted
            #else
               #pragma omp for schedule(static) nowait
            #endif
            for (int global = 0; global < upperbound1; global++){
               const int orb_k    = orb_i - trianglefunction( orb_i, global ) - 1;
               const int orb_j    = global - (orb_k*(orb_k+1))/2;
               const int irrep_jk = Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) );
               if ( irrep_jk == irrep_mn ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
                  #endif
                  {
                     const double d22 = diagram2_22_23(denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                             S0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin;
                     const double d23 = ((( orb_j == orb_k ) || ( orb_m == orb_n )) ? 0.0 :
                                        diagram2_22_23(denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                             S1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin );
                     const double d24 = diagram24_27_28(denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                              S0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin;
                     if (( orb_j != orb_k ) && ( orb_m != orb_n )){
                        diagram25_26(denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k], S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2);
                     }
                     const double d25 = ((( orb_j == orb_k ) || ( orb_m == orb_n )) ? 0.0 : workmem[0]  * prefactor_spin );
                     const double d26 = ((( orb_j == orb_k ) || ( orb_m == orb_n )) ? 0.0 : workmem2[0] * prefactor_spin );
                     const double d27 = (( orb_m == orb_n ) ? 0.0 :
                                        diagram24_27_28(denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                              S1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin );
                     const double d28 = (( orb_j == orb_k ) ? 0.0 :
                                        diagram24_27_28(denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                              S0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_n, orb_i,  2*d22 + 2*d23 - 2*d24 + d25                   );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_n, orb_m, orb_i,  2*d22 - 2*d23 - 2*d24 - d25                   );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_i, orb_n,   -d22 -   d23 +   d24       + d26 + d27 + d28 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_n, orb_i, orb_m,   -d22 +   d23 +   d24       - d26 - d27 + d28 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_m, orb_n,   -d22 +   d23 +   d24       - d26 + d27 - d28 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_n, orb_m,   -d22 -   d23 +   d24       + d26 - d27 - d28 );
                  }
                  
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
                  #endif
                  {
                     const double d29 = diagram29_30_31_32(denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 F0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d31 = (( orb_m == orb_n ) ? d29 : workmem[0] * prefactor_spin );
                     const double d30 = diagram29_30_31_32(denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 F1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d32 = (( orb_m == orb_n ) ? d30 : workmem[0] * prefactor_spin );
                     
                     const double d34 = diagram34_37_38(denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                              F1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin;
                     const double d37 = workmem[0]  * prefactor_spin;
                     const double d38 = workmem2[0] * prefactor_spin;
                     
                     const double d40 = (( orb_m == orb_n ) ? d34 :
                                        diagram40_43_44(denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                              F1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2) * prefactor_spin );
                     const double d43 = (( orb_m == orb_n ) ? d37 : workmem[0]  * prefactor_spin );
                     const double d44 = (( orb_m == orb_n ) ? d38 : workmem2[0] * prefactor_spin );
                     
                     const double d33 = diagram33_39(denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                           F0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d39 = (( orb_m == orb_n ) ? d33 : workmem[0] * prefactor_spin );
                     
                     const double d35 = diagram35_41(denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                           F1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d41 = (( orb_m == orb_n ) ? d35 : workmem[0] * prefactor_spin );
                     
                     const double d36 = diagram36_42(denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                           F0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d42 = (( orb_m == orb_n ) ? d36 : workmem[0] * prefactor_spin );
                     
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_n, orb_i,  4*(d29 + d33) );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_k, orb_i, -2*(d29 + d30 + d33 + d34) );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_i, orb_n, -2*(d29 + d33 + d35) );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_n, orb_i, orb_k,     d29 + d30 + d33 + d35 + d36 + d37 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_k, orb_n,     d29 + d30 + d33 + d35 + d36 + d38 );
                     set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_n, orb_k, -2*(d29 + d33 + d36) );
                     
                     if ( orb_m < orb_n ){
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_k, orb_m, orb_i,  4*(d31 + d39) );
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_m, orb_k, orb_i, -2*(d31 + d32 + d39 + d40) );
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_k, orb_i, orb_m, -2*(d31 + d39 + d41) );
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_m, orb_i, orb_k,     d31 + d32 + d39 + d41 + d42 + d43 );
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_k, orb_m,     d31 + d32 + d39 + d41 + d42 + d44 );
                        set_dmrg_index( orb_j, orb_n, orb_i, orb_i, orb_m, orb_k, -2*(d31 + d39 + d42) );
                     }
                  
                  }
                  
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_absigma( orb_j, orb_k ) )
                  #endif
                  {
                     const double d45 = diagram45_46_47_48(denT, S0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 F0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin;
                     const double d47 = (( orb_m == orb_n ) ? d45 : workmem[0] * prefactor_spin );
                     
                     const double d46 = (( orb_j == orb_k ) ? 0.0 :
                                        diagram45_46_47_48(denT, S1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 F1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_m < orb_n) * prefactor_spin );
                     const double d48 = (( orb_j == orb_k ) ? 0.0 : (( orb_m == orb_n ) ? d46 : workmem[0] * prefactor_spin ));
                     
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_n, orb_i, orb_i, d45 + d46 );
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_n, orb_i, d45 - d46 );
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_i, orb_n, - 2 * d45 );
                     
                     if ( orb_m < orb_n ){
                        set_dmrg_index( orb_j, orb_k, orb_n, orb_m, orb_i, orb_i, d47 + d48 );
                        set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_m, orb_i, d47 - d48 );
                        set_dmrg_index( orb_j, orb_k, orb_n, orb_i, orb_i, orb_m, - 2 * d47 );
                     }
                  }
                  
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIRANK == MPIchemps2::owner_cdf( L, orb_j, orb_k ) )
                  #endif
                  {
                     const double d49 = diagram49_50_51_52(denT, F0tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 S0tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_j < orb_k) * prefactor_spin;
                     const double d51 = (( orb_j == orb_k ) ? d49 : workmem[0] * prefactor_spin );
                     
                     const double d50 = (( orb_m == orb_n ) ? 0.0 :
                                        diagram49_50_51_52(denT, F1tensors[orb_i-1][orb_k-orb_j][orb_i-1-orb_k],
                                                                 S1tensors[orb_i  ][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, orb_j < orb_k) * prefactor_spin );
                     const double d52 = (( orb_m == orb_n ) ? 0.0 : (( orb_j == orb_k ) ? d50 : workmem[0] * prefactor_spin ));
                     
                     set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_i, orb_i, - 2 * d49 );
                     set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_k, orb_i, d49 + d50 );
                     set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_i, orb_k, d49 - d50 );
                     
                     if ( orb_j < orb_k ){
                        set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_i, orb_i, - 2 * d51 );
                        set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_j, orb_i, d51 + d52 );
                        set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_i, orb_j, d51 - d52 );
                     }
                  }
                  
               }
            }
            
            const int irrep_imn = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_i ), book->gIrrep( orb_m ) ), book->gIrrep( orb_n ) );
            int counter_jkl = 0;
            for (int global = 0; global < upperbound4; global++){
               tripletrianglefunction( global, jkl );
               const int orb_j     = jkl[0];
               const int orb_k     = jkl[1];
               const int orb_l     = jkl[2];
               const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
               if ( irrep_jkl == irrep_imn ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK ){ counter_jkl++; }
                  #else
                  counter_jkl++;
                  #endif
               }
            }
            Tensor3RDM * bcd_S0_doublet = NULL;
            Tensor3RDM * bcd_S1_doublet = NULL;
            Tensor3RDM * bcd_S1_quartet = NULL;
            if ( counter_jkl > 0 ){
                                      bcd_S0_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
               if ( orb_m != orb_n ){ bcd_S1_doublet = new Tensor3RDM( orb_i, -2, 1, 1, irrep_imn, true, book );
                                      bcd_S1_quartet = new Tensor3RDM( orb_i, -2, 3, 1, irrep_imn, true, book ); }
                                      fill_bcd_S0( denT, bcd_S0_doublet,                 S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem );
               if ( orb_m != orb_n ){ fill_bcd_S1( denT, bcd_S1_doublet, bcd_S1_quartet, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 ); }
            }
            
            for (int global = 0; global < upperbound4; global++){
               tripletrianglefunction( global, jkl );
               const int orb_j     = jkl[0];
               const int orb_k     = jkl[1];
               const int orb_l     = jkl[2];
               const int irrep_jkl = Irreps::directProd( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ), book->gIrrep( orb_l ) );
               if ( irrep_jkl == irrep_imn ){
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
                  #endif
                  {
                     const int cnt1 = orb_k - orb_j;
                     const int cnt2 = orb_l - orb_k;
                     const int cnt3 = orb_i - 1 - orb_l;
                     Tensor3RDM * left[9];
                     double results[9];
                     if ( cnt1 + cnt2 > 0 ){
                        left[0] = dm3_a_J0_doublet[cnt1][cnt2][cnt3];
                        left[1] = dm3_a_J1_doublet[cnt1][cnt2][cnt3];
                        left[2] = dm3_a_J1_quartet[cnt1][cnt2][cnt3];
                        diagram90_93( denT, left, S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                        const double d90 =                         workmem[0] * prefactor_spin;
                        const double d93 = (( cnt1 == 0 ) ? 0.0 : workmem2[0] * prefactor_spin );
                        if ( orb_n != orb_m ){ diagram91_92_94( denT, left, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, results ); }
                        const double d91 = ((( orb_n == orb_m ) || ( cnt1 == 0      )) ? 0.0 : results[0] * prefactor_spin ); // dm3_a_J1_doublet
                        const double d92 = ((( orb_n == orb_m ) || ( cnt1*cnt2 == 0 )) ? 0.0 : results[1] * prefactor_spin ); // dm3_a_J1_quartet
                        const double d94 = ( ( orb_n == orb_m )                        ? 0.0 : results[2] * prefactor_spin ); // dm3_a_J0_doublet
                        diagram95_98( denT, left, S0tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2 );
                        const double d95 =                         workmem[0] * prefactor_spin;
                        const double d98 = (( cnt1 == 0 ) ? 0.0 : workmem2[0] * prefactor_spin );
                        if ( orb_n != orb_m ){ diagram96_97_99( denT, left, S1tensors[orb_i][orb_n-orb_m][orb_m-1-orb_i], workmem, workmem2, results ); }
                        const double d96 = ((( orb_n == orb_m ) || ( cnt1 == 0      )) ? 0.0 : results[0] * prefactor_spin ); // dm3_a_J1_doublet
                        const double d97 = ((( orb_n == orb_m ) || ( cnt1*cnt2 == 0 )) ? 0.0 : results[1] * prefactor_spin ); // dm3_a_J1_quartet
                        const double d99 = ( ( orb_n == orb_m )                        ? 0.0 : results[2] * prefactor_spin ); // dm3_a_J0_doublet
                        const double d9095 = d90 + d95;
                        const double d9196 = d91 + d96;
                        const double d9297 = d92 + d97;
                        const double d9398 = d93 + d98;
                        const double d9499 = d94 + d99;
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_n, orb_i, -2*d9095 + 2*d9196 + d9297 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_n, orb_m, orb_i, -2*d9095 - 2*d9196 - d9297 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_i, orb_n,    d9095 +   d9196 - d9297 + d9398 + d9499 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_n, orb_i, orb_m,    d9095 -   d9196 + d9297 + d9398 - d9499 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_m, orb_n,    d9095 -   d9196 + d9297 - d9398 + d9499 );
                        set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_n, orb_m,    d9095 +   d9196 - d9297 - d9398 - d9499 );
                     }
                     {
   
   const double sq3 = sqrt( 3.0 );
   const double d120_125 = ((     cnt2==0)                  ?0.0:      bcd_S0_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d124_129 = ((cnt1*cnt2==0)                  ?0.0:  sq3*bcd_S0_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d110_115 = ((cnt1+cnt2==0)                  ?0.0:      bcd_S0_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d112_117 = ((cnt1+cnt2==0)                  ?0.0: -sq3*bcd_S0_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d100_105 =                                            -bcd_S0_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin;
   const double d102_107 =                                        -sq3*bcd_S0_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin;
   const double d123_128 = (((orb_n==orb_m)||(     cnt2==0))?0.0:  sq3*bcd_S1_doublet->contract(dm3_b_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d121_126 = (((orb_n==orb_m)||(cnt1*cnt2==0))?0.0:      bcd_S1_doublet->contract(dm3_b_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d111_116 = (((orb_n==orb_m)||(cnt1+cnt2==0))?0.0: -sq3*bcd_S1_doublet->contract(dm3_c_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d113_118 = (((orb_n==orb_m)||(cnt1+cnt2==0))?0.0:      bcd_S1_doublet->contract(dm3_c_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d101_106 = ( (orb_n==orb_m)                 ?0.0:  sq3*bcd_S1_doublet->contract(dm3_d_J0_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d103_108 = ( (orb_n==orb_m)                 ?0.0:      bcd_S1_doublet->contract(dm3_d_J1_doublet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d122_127 = (((orb_n==orb_m)||(cnt1*cnt2==0))?0.0:      bcd_S1_quartet->contract(dm3_b_J1_quartet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d114_119 = (((orb_n==orb_m)||(cnt1+cnt2==0))?0.0:   -2*bcd_S1_quartet->contract(dm3_c_J1_quartet[cnt1][cnt2][cnt3])*prefactor_spin);
   const double d104_109 = (((orb_n==orb_m)||(     cnt2==0))?0.0:    2*bcd_S1_quartet->contract(dm3_d_J1_quartet[cnt1][cnt2][cnt3])*prefactor_spin);

                        if ( cnt2 > 0 ){ // dm3_b if NOT k==l
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_j, orb_k, orb_i,  -d120_125 + 3*d121_126 +   d123_128 - d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_k, orb_j, orb_i,  -d120_125 - 3*d121_126 +   d123_128 + d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_j, orb_i, orb_k,  -d120_125 - 3*d121_126 -   d123_128 - d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_k, orb_i, orb_j,  -d120_125 + 3*d121_126 -   d123_128 + d124_129 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_i, orb_j, orb_k, 2*d120_125 + 2*d121_126 + 2*d122_127 );
                           set_dmrg_index( orb_l, orb_m, orb_n, orb_i, orb_k, orb_j, 2*d120_125 - 2*d121_126 - 2*d122_127 );
                        }
                        if ( cnt1 + cnt2 > 0 ){ // dm3_c if NOT j==k==l
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_l, orb_i, 2*d110_115 + 2*d111_116 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_l, orb_j, orb_i,  -d110_115 -   d111_116 - d112_117 - 3*d113_118 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_j, orb_i, orb_l, 2*d110_115 - 2*d111_116 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_l, orb_i, orb_j,  -d110_115 +   d111_116 - d112_117 + 3*d113_118 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_j, orb_l,  -d110_115 +   d111_116 + d112_117 +   d113_118 + d114_119 );
                           set_dmrg_index( orb_k, orb_m, orb_n, orb_i, orb_l, orb_j,  -d110_115 -   d111_116 + d112_117 -   d113_118 - d114_119 );
                        }
                        { // dm3_d for all j <= k <= l
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_l, orb_i, 2*d100_105 + 2*d101_106 );
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_l, orb_k, orb_i,  -d100_105 -   d101_106 - d102_107 - 3*d103_108 );
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_k, orb_i, orb_l, 2*d100_105 - 2*d101_106 );
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_l, orb_i, orb_k,  -d100_105 +   d101_106 - d102_107 + 3*d103_108 );
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_k, orb_l,  -d100_105 +   d101_106 + d102_107 +   d103_108 + d104_109 );
                           set_dmrg_index( orb_j, orb_m, orb_n, orb_i, orb_l, orb_k,  -d100_105 -   d101_106 + d102_107 -   d103_108 - d104_109 );
                        }
                     }
                     
                     /* Here come diagrams 130 to 189 */
                     
                  }
               }
            }
            
            if ( bcd_S0_doublet != NULL ){ delete bcd_S0_doublet; }
            if ( bcd_S1_doublet != NULL ){ delete bcd_S1_doublet; }
            if ( bcd_S1_quartet != NULL ){ delete bcd_S1_quartet; }
            
            #ifdef CHEMPS2_MPI_COMPILATION
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
      
      const int upperbound53 = upperbound4 * ( L - orb_i - 1 );
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int combined = 0; combined < upperbound53; combined++){
         const int global = combined % upperbound4;
         const int orb_m  = orb_i + 1 + ( combined / upperbound4 );
         tripletrianglefunction( global, jkl );
         const int orb_j = jkl[0];
         const int orb_k = jkl[1];
         const int orb_l = jkl[2];
         if ( Irreps::directProd( book->gIrrep( orb_j ), book->gIrrep( orb_k ) ) == Irreps::directProd( book->gIrrep( orb_l ), book->gIrrep( orb_m ) ) ){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
            #endif
            {
               const int cnt1 = orb_k - orb_j;
               const int cnt2 = orb_l - orb_k;
               const int cnt3 = orb_i - 1 - orb_l;
               if ( cnt1 + cnt2 > 0 ){
                  diagram53_54( denT, dm3_a_J0_doublet[cnt1][cnt2][cnt3], dm3_a_J1_doublet[cnt1][cnt2][cnt3], Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2 );
                  const double d53 =  workmem[0] * prefactor_spin;
                  const double d54 = workmem2[0] * prefactor_spin;
                  set_dmrg_index( orb_j, orb_k, orb_l, orb_m, orb_i, orb_i, d53 - d54 );
                  set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_m, orb_i, d53 + d54 );
                  set_dmrg_index( orb_j, orb_k, orb_l, orb_i, orb_i, orb_m, - 2 * d53 );
               }
               Tensor3RDM * left[6];
               double results[18];
               left[0] = dm3_b_J0_doublet[cnt1][cnt2][cnt3];
               left[1] = dm3_b_J1_doublet[cnt1][cnt2][cnt3];
               left[2] = dm3_c_J0_doublet[cnt1][cnt2][cnt3];
               left[3] = dm3_c_J1_doublet[cnt1][cnt2][cnt3];
               left[4] = dm3_d_J0_doublet[cnt1][cnt2][cnt3];
               left[5] = dm3_d_J1_doublet[cnt1][cnt2][cnt3];
               {
                  diagram55to60( denT, left, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2, results );
                  const double d55 = results[0] * prefactor_spin;
                  const double d56 = results[1] * prefactor_spin;
                  const double d57 = results[2] * prefactor_spin;
                  const double d58 = results[3] * prefactor_spin;
                  const double d59 = results[4] * prefactor_spin;
                  const double d60 = results[5] * prefactor_spin;
                  if ( cnt2 > 0 ){
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_l, orb_i, orb_i, d55 + d56 );
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_l, orb_i, d55 - d56 );
                     set_dmrg_index( orb_j, orb_k, orb_m, orb_i, orb_i, orb_l, - 2 * d55 );
                  }
                  if ( cnt1 + cnt2 > 0 ){
                     set_dmrg_index( orb_j, orb_l, orb_m, orb_k, orb_i, orb_i, - 2 * d57 );
                     set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_k, orb_i, d57 - d58 );
                     set_dmrg_index( orb_j, orb_l, orb_m, orb_i, orb_i, orb_k, d57 + d58 );
                  }
                  set_dmrg_index( orb_k, orb_l, orb_m, orb_j, orb_i, orb_i, - 2 * d59 );
                  set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_j, orb_i, d59 + d60 );
                  set_dmrg_index( orb_k, orb_l, orb_m, orb_i, orb_i, orb_j, d59 - d60 );
               }
               {
                  diagram61_62_70_71_80_81( denT, left, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2, results );
                  const double d61 = results[0]  * prefactor_spin;
                  const double d62 = results[1]  * prefactor_spin;
                  const double d70 = results[2]  * prefactor_spin;
                  const double d71 = results[3]  * prefactor_spin;
                  const double d80 = results[4]  * prefactor_spin;
                  const double d81 = results[5]  * prefactor_spin;
                  diagramJ2half_3_2_1( denT, left, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2, results );
                  const double d63 = results[0]  * prefactor_spin;
                  const double d64 = results[1]  * prefactor_spin;
                  const double d65 = results[2]  * prefactor_spin;
                  const double d66 = results[3]  * prefactor_spin;
                  const double d67 = results[4]  * prefactor_spin;
                  const double d68 = results[5]  * prefactor_spin;
                  const double d72 = results[6]  * prefactor_spin;
                  const double d73 = results[7]  * prefactor_spin;
                  const double d74 = results[8]  * prefactor_spin;
                  const double d75 = results[9]  * prefactor_spin;
                  const double d76 = results[10] * prefactor_spin;
                  const double d77 = results[11] * prefactor_spin;
                  const double d82 = results[12] * prefactor_spin;
                  const double d83 = results[13] * prefactor_spin;
                  const double d84 = results[14] * prefactor_spin;
                  const double d85 = results[15] * prefactor_spin;
                  const double d86 = results[16] * prefactor_spin;
                  const double d87 = results[17] * prefactor_spin;
                  left[0] = dm3_b_J1_quartet[cnt1][cnt2][cnt3];
                  left[1] = dm3_c_J1_quartet[cnt1][cnt2][cnt3];
                  left[2] = dm3_d_J1_quartet[cnt1][cnt2][cnt3];
                  diagramJ2oneandhalf_3_2_1( denT, left, Ltensors[orb_i][orb_m-1-orb_i], workmem, workmem2, results );
                  const double d69 = results[0] * prefactor_spin;
                  const double d78 = results[1] * prefactor_spin;
                  const double d79 = results[2] * prefactor_spin;
                  const double d88 = results[3] * prefactor_spin;
                  const double d89 = results[4] * prefactor_spin;
                  if ( cnt2 > 0 ){
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_m, orb_i, -2*d61 - 2*d62 + d63 + d64 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_l, orb_i, -2*d61 + 2*d62 + d63 - d64 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_l, orb_i, orb_m, d61 + d62 + d65 + d66 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_m, orb_i, orb_l, d61 - d62 + d67 + d68 + d69 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_m, orb_l, d61 + d62 + d67 - d68 - d69 );
                     set_dmrg_index( orb_j, orb_k, orb_i, orb_i, orb_l, orb_m, d61 - d62 + d65 - d66 );
                  }
                  if ( cnt1 + cnt2 > 0 ){
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_m, orb_i, 4*d70 + 2*d72 );
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_m, orb_k, orb_i, -2*d70 + 2*d71 - d72 + d73 );
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_k, orb_i, orb_m, -2*d70 + 2*d74 );
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_m, orb_i, orb_k, d70 - d71 - d74 + d76 + d78 );
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_k, orb_m, d70 - d71 - d74 + d75 );
                     set_dmrg_index( orb_j, orb_l, orb_i, orb_i, orb_m, orb_k, -2*d70 - d72 + d77 + d79 );
                  }
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_l, orb_i, 4*d80 + 2*d82 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_l, orb_k, orb_i, -2*d80 - 2*d81 - d82 - d83 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_k, orb_i, orb_l, -2*d80 + 2*d84 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_l, orb_i, orb_k, d80 + d81 - d84 - d85 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_k, orb_l, d80 + d81 - d84 + d86 + d88 );
                  set_dmrg_index( orb_j, orb_m, orb_i, orb_i, orb_l, orb_k, -2*d80 - d82 + d87 + d89 );
               }
               
            }
         }
      }
      
      delete [] workmem;
      delete [] workmem2;
   
   }

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
                                                         : ( phase( TwoSL + 1 - TwoSLprime )
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
                                  * phase( TwoSR + 3 + TwoSLprime )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
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
                            * phase( TwoSL + 3 - TwoSR )
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
                                  * phase( 2 * TwoSR + 2 )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
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
                            * phase( TwoSL + 3 - TwoSR )
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
                                  * phase( TwoSL + TwoSLprime )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSLprime, TwoSR )
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
                                  * phase( TwoSL + TwoSLprime + 3 )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
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
                            * phase( TwoSL + 1 - TwoSLprime )
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
                                  * phase( 2 * TwoSR )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
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
                            * phase( TwoSL + 1 - TwoSLprime )
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
                                  * phase( TwoSR + TwoSLprime )
                                  * gsl_sf_coupling_6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL )
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

double CheMPS2::ThreeDM::diagram2_22_23(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2) const{

   assert( left->get_irrep() == right->get_irrep() );
   assert( left->get_2j()    == right->get_2j()    );
   assert( left->get_nelec() == right->get_nelec() );
   const int orb_i = denT->gIndex();
   const int two_j = left->get_2j();
   const int nelec = left->get_nelec();
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int ILxIjxIk = Irreps::directProd( IL, left->get_irrep() );
               
               for ( int TwoSLprime = TwoSL-two_j; TwoSLprime <= TwoSL+two_j; TwoSLprime+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL-nelec,   TwoSLprime, ILxIjxIk );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-nelec+2, TwoSLprime, ILxIjxIk );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup    =  denT->gStorage( NL,         TwoSL,      IL,       NL+2,       TwoSL,      IL       );
                     double * Tdown  =  denT->gStorage( NL-nelec,   TwoSLprime, ILxIjxIk, NL-nelec+2, TwoSLprime, ILxIjxIk );
                     double * Lblock =  left->gStorage( NL-nelec,   TwoSLprime, ILxIjxIk, NL,         TwoSL,      IL       );
                     double * Rblock = right->gStorage( NL-nelec+2, TwoSLprime, ILxIjxIk, NL+2,       TwoSL,      IL       );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                     
                     int length = dimLdown * dimRdown;
                     int inc = 1;
                     total += ( TwoSL + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );

                  }
               }
            }
         }
      }
   }
   
   return total;

}

void CheMPS2::ThreeDM::diagram25_26(TensorT * denT, TensorS1 * left, TensorS1 * right, double * workmem, double * workmem2) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_d25 = 0.0;
   double total_d26 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               const int IRup   = Irreps::directProd( ILup,   book->gIrrep( orb_i ) ); // IL x Ii
               const int ILdown = Irreps::directProd( ILup,   left->get_irrep()     ); // IL x Ij x Ik
               const int IRdown = Irreps::directProd( ILdown, book->gIrrep( orb_i ) ); // IL x Ii x Ij x Ik = IL x Ii x Im x In
               
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
               
                     for ( int TwoSLdown = TwoSLup-2; TwoSLdown <= TwoSLup+2; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL-2, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSLdown-1; TwoSRdown <= TwoSLdown+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSRup - TwoSRdown ) <= 2 )){
                                 
                                 double * Tup    =  denT->gStorage( NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup   );
                                 double * Tdown  =  denT->gStorage( NL-2, TwoSLdown, ILdown, NL-1, TwoSRdown, IRdown );
                                 double * Lblock =  left->gStorage( NL-2, TwoSLdown, ILdown, NL,   TwoSLup,   ILup   );
                                 double * Rblock = right->gStorage( NL-1, TwoSRdown, IRdown, NL+1, TwoSRup,   IRup   );
                                 
                                 char trans   = 'T';
                                 char notrans = 'N';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                                 
                                 int length = dimLdown * dimRdown;
                                 int inc = 1;
                                 const double contracted = sqrt( 1.0 * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 ) 
                                                         * ddot_( &length, workmem2, &inc, Tdown, &inc );
                                 total_d25 += phase( TwoSLup + TwoSRdown + 3 )
                                            * gsl_sf_coupling_6j( TwoSLup, TwoSRup, 1, TwoSRdown, TwoSLdown, 2 )
                                            * contracted;
                                 total_d26 += 3 * phase( 2 + TwoSRup + TwoSRdown + TwoSLup + TwoSLdown )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSLdown, TwoSRdown )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSLup   )
                                            * contracted;

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
   
   workmem[0]  = total_d25;
   workmem2[0] = total_d26;

}

double CheMPS2::ThreeDM::diagram24_27_28(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2) const{

   assert( left->get_irrep() == right->get_irrep() );
   assert( left->get_nelec() == right->get_nelec() );
   const int orb_i = denT->gIndex();
   const int twoJ1 =  left->get_2j();
   const int twoJ2 = right->get_2j();
   const int nelec = left->get_nelec();
   
   double total = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               const int IRup   = Irreps::directProd( ILup,   book->gIrrep( orb_i ) ); // IL x Ii
               const int ILdown = Irreps::directProd( ILup,   left->get_irrep()     ); // IL x Ij x Ik
               const int IRdown = Irreps::directProd( ILdown, book->gIrrep( orb_i ) ); // IL x Ii x Ij x Ik = IL x Ii x Im x In
               
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
               
                     for ( int TwoSLdown = TwoSLup-twoJ1; TwoSLdown <= TwoSLup+twoJ1; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL-nelec, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSLdown-1; TwoSRdown <= TwoSLdown+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL-nelec+1, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSRup - TwoSRdown ) <= twoJ2 )){
                                 
                                 double * Tup    =  denT->gStorage( NL,         TwoSLup,   ILup,   NL+1,       TwoSRup,   IRup   );
                                 double * Tdown  =  denT->gStorage( NL-nelec,   TwoSLdown, ILdown, NL-nelec+1, TwoSRdown, IRdown );
                                 double * Lblock =  left->gStorage( NL-nelec,   TwoSLdown, ILdown, NL,         TwoSLup,   ILup   );
                                 double * Rblock = right->gStorage( NL-nelec+1, TwoSRdown, IRdown, NL+1,       TwoSRup,   IRup   );
                                 
                                 char trans   = 'T';
                                 char notrans = 'N';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                                 
                                 int length = dimLdown * dimRdown;
                                 int inc = 1;
                                 total += sqrt( 1.0 * ( twoJ1 + 1 ) * ( twoJ2 + 1 ) * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 )
                                        * phase( twoJ1 + TwoSRup + TwoSRdown + TwoSLup + TwoSLdown )
                                        * gsl_sf_coupling_6j( 1, 1, twoJ1, TwoSLup, TwoSLdown, TwoSRdown )
                                        * gsl_sf_coupling_6j( 1, 1, twoJ2, TwoSRup, TwoSRdown, TwoSLup   )
                                        * ddot_( &length, workmem2, &inc, Tdown, &inc );

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
   
   return total;

}

double CheMPS2::ThreeDM::diagram29_30_31_32(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const{

   assert( left->get_irrep() == right->get_irrep() );
   assert( left->get_2j()    == right->get_2j()    );
   assert( left->get_nelec() == right->get_nelec() );
   assert( left->get_nelec() == 0                  );
   const int orb_i = denT->gIndex();
   const int two_j = left->get_2j();
   
   double total_29_30 = 0.0;
   double total_31_32 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int ILxIjxIk = Irreps::directProd( IL, left->get_irrep() );
               
               for ( int TwoSLprime = TwoSL-two_j; TwoSLprime <= TwoSL+two_j; TwoSLprime+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLprime, ILxIjxIk );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSLprime, ILxIjxIk );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup     =  denT->gStorage( NL,   TwoSL,      IL,       NL+2, TwoSL,      IL       );
                     double * Lblock  =  left->gStorage( NL,   TwoSLprime, ILxIjxIk, NL,   TwoSL,      IL       );
                     double * Rblock  = right->gStorage( NL+2, TwoSLprime, ILxIjxIk, NL+2, TwoSL,      IL       );
                     double * Tdown   =  denT->gStorage( NL,   TwoSLprime, ILxIjxIk, NL+2, TwoSLprime, ILxIjxIk );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     int length   = dimLdown * dimRdown;
                     int inc      = 1;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                     total_29_30 += sqrt(1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 )) * phase(TwoSL - TwoSLprime) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     if ( doTR ){
                        double * Rblock2 = right->gStorage( NL+2, TwoSL, IL, NL+2, TwoSLprime, ILxIjxIk );
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock2, &dimRup, &beta, workmem2, &dimLdown );
                        total_31_32 += ( TwoSL + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_31_32;
   return total_29_30;

}

double CheMPS2::ThreeDM::diagram34_37_38(TensorT * denT, TensorF1 * left, TensorF1 * right, double * workmem, double * workmem2) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_d34 = 0.0;
   double total_d37 = 0.0;
   double total_d38 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               const int IRup   = Irreps::directProd( ILup,   book->gIrrep( orb_i ) ); // IL x Ii
               const int ILdown = Irreps::directProd( ILup,   left->get_irrep()     ); // IL x Ij x Ik
               const int IRdown = Irreps::directProd( ILdown, book->gIrrep( orb_i ) ); // IL x Ii x Ij x Ik = IL x Ii x Im x In
               
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
               
                     for ( int TwoSLdown = TwoSLup-2; TwoSLdown <= TwoSLup+2; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSLdown-1; TwoSRdown <= TwoSLdown+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSRup - TwoSRdown ) <= 2 )){
                                 
                                 double * Tup    =  denT->gStorage( NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup   );
                                 double * Tdown  =  denT->gStorage( NL,   TwoSLdown, ILdown, NL+1, TwoSRdown, IRdown );
                                 double * Lblock =  left->gStorage( NL,   TwoSLdown, ILdown, NL,   TwoSLup,   ILup   );
                                 double * Rblock = right->gStorage( NL+1, TwoSRdown, IRdown, NL+1, TwoSRup,   IRup   );
                                 
                                 char trans   = 'T';
                                 char notrans = 'N';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                                 
                                 int length = dimLdown * dimRdown;
                                 int inc = 1;
                                 const double contracted = sqrt( 1.0 * ( TwoSLup + 1 ) * ( TwoSRup + 1 ) ) * ( TwoSRdown + 1 )
                                                         * ddot_( &length, workmem2, &inc, Tdown, &inc );
                                 total_d34 += 0.5 * phase( TwoSLup + TwoSRup + 3 )
                                            * gsl_sf_coupling_6j( TwoSLup, TwoSRup, 1, TwoSRdown, TwoSLdown, 2 )
                                            * contracted;
                                 total_d37 += 3 * phase( TwoSRup + TwoSRdown + 2 * TwoSLup )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSLdown, TwoSRup   )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSLdown )
                                            * contracted;
                                 total_d38 += 3 * phase( TwoSLup + TwoSLdown + 2 * TwoSRup )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSLup   )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSLdown, TwoSRdown )
                                            * contracted;

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
   
   workmem[0]  = total_d37;
   workmem2[0] = total_d38;
   return total_d34;

}

double CheMPS2::ThreeDM::diagram40_43_44(TensorT * denT, TensorF1 * left, TensorF1 * right, double * workmem, double * workmem2) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_d40 = 0.0;
   double total_d43 = 0.0;
   double total_d44 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               const int IRup   = Irreps::directProd( ILup,   book->gIrrep( orb_i ) ); // IL x Ii
               const int ILdown = Irreps::directProd( ILup,   left->get_irrep()     ); // IL x Ij x Ik
               const int IRdown = Irreps::directProd( ILdown, book->gIrrep( orb_i ) ); // IL x Ii x Ij x Ik = IL x Ii x Im x In
               
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
               
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
               
                     for ( int TwoSLdown = TwoSLup-2; TwoSLdown <= TwoSLup+2; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSLdown-1; TwoSRdown <= TwoSLdown+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSRup - TwoSRdown ) <= 2 )){
                                 
                                 double * Tup    =  denT->gStorage( NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup   );
                                 double * Tdown  =  denT->gStorage( NL,   TwoSLdown, ILdown, NL+1, TwoSRdown, IRdown );
                                 double * Lblock =  left->gStorage( NL,   TwoSLdown, ILdown, NL,   TwoSLup,   ILup   );
                                 double * Rblock = right->gStorage( NL+1, TwoSRup,   IRup,   NL+1, TwoSRdown, IRdown );
                                 
                                 char notrans = 'N';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup, &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRup, &beta, workmem2, &dimLdown );
                                 
                                 int length = dimLdown * dimRdown;
                                 int inc = 1;
                                 const double contracted = sqrt( 1.0 * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 )
                                                         * ddot_( &length, workmem2, &inc, Tdown, &inc );
                                 total_d40 += 0.5 * phase( TwoSLup + TwoSRdown + 3 )
                                            * gsl_sf_coupling_6j( TwoSLup, TwoSRup, 1, TwoSRdown, TwoSLdown, 2 )
                                            * contracted;
                                 total_d43 += 3 * phase( 2 * TwoSRup + 2 * TwoSLup )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSLdown, TwoSRup   )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSLdown )
                                            * contracted;
                                 total_d44 += 3 * phase( TwoSLup + TwoSLdown + TwoSRup + TwoSRdown )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSLup   )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSLdown, TwoSRdown )
                                            * contracted;

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
   
   workmem[0]  = total_d43;
   workmem2[0] = total_d44;
   return total_d40;

}

double CheMPS2::ThreeDM::diagram33_39(TensorT * denT, TensorF0 * left, TensorF0 * right, double * workmem, double * workmem2, const bool do39) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_33 = 0.0;
   double total_39 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILdown = Irreps::directProd( IL,   left->get_irrep()     );
            const int IRup   = Irreps::directProd( IL,   book->gIrrep( orb_i ) );
            const int IRdown = Irreps::directProd( IRup, right->get_irrep()    );
            
            int dimLup   = book->gCurrentDim( orb_i, NL, TwoSL, IL     );
            int dimLdown = book->gCurrentDim( orb_i, NL, TwoSL, ILdown );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRup   );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRdown );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                     
                     double * Tup     =  denT->gStorage( NL,   TwoSL, IL,     NL+1, TwoSR, IRup   );
                     double * Tdown   =  denT->gStorage( NL,   TwoSL, ILdown, NL+1, TwoSR, IRdown );
                     double * Lblock  =  left->gStorage( NL,   TwoSL, ILdown, NL,   TwoSL, IL     );
                     double * Rblock  = right->gStorage( NL+1, TwoSR, IRdown, NL+1, TwoSR, IRup   );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     int length   = dimLdown * dimRdown;
                     int inc      = 1;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                     total_33 += 0.5 * ( TwoSR + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     
                     if ( do39 ){
                        double * Rblock2 = right->gStorage( NL+1, TwoSR, IRup,   NL+1, TwoSR, IRdown );
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock2, &dimRup, &beta, workmem2, &dimLdown );
                        total_39 += 0.5 * ( TwoSR + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_39;
   return total_33;

}

double CheMPS2::ThreeDM::diagram35_41(TensorT * denT, TensorF0 * left, TensorF1 * right, double * workmem, double * workmem2, const bool do41) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_35 = 0.0;
   double total_41 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILdown = Irreps::directProd( IL,   left->get_irrep()     );
            const int IRup   = Irreps::directProd( IL,   book->gIrrep( orb_i ) );
            const int IRdown = Irreps::directProd( IRup, right->get_irrep()    );
            
            int dimLup   = book->gCurrentDim( orb_i, NL, TwoSL, IL     );
            int dimLdown = book->gCurrentDim( orb_i, NL, TwoSL, ILdown );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 )){
               
               for ( int TwoSRup = TwoSL-1; TwoSRup <= TwoSL+1; TwoSRup+=2 ){
                  
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
                  
                     for ( int TwoSRdown = TwoSL-1; TwoSRdown <= TwoSL+1; TwoSRdown+=2 ){
                  
                        int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSRdown, IRdown );
                        
                        if ( dimRdown > 0 ){
                           
                           double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,     NL+1, TwoSRup,   IRup   );
                           double * Tdown   =  denT->gStorage( NL,   TwoSL,     ILdown, NL+1, TwoSRdown, IRdown );
                           double * Lblock  =  left->gStorage( NL,   TwoSL,     ILdown, NL,   TwoSL,     IL     );
                           double * Rblock  = right->gStorage( NL+1, TwoSRdown, IRdown, NL+1, TwoSRup,   IRup   );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           int length   = dimLdown * dimRdown;
                           int inc      = 1;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                           dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                           const double wigner6j = gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSRdown, TwoSL );
                           total_35 += 0.5 * sqrt( 6.0 * ( TwoSRup + 1 ) ) * ( TwoSRdown + 1 )
                                     * phase( TwoSL + TwoSRdown + 3 ) * wigner6j
                                     * ddot_( &length, workmem2, &inc, Tdown, &inc );
                           
                           if ( do41 ){
                              double * Rblock2 = right->gStorage( NL+1, TwoSRup, IRup, NL+1, TwoSRdown, IRdown );
                              dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock2, &dimRup, &beta, workmem2, &dimLdown );
                              total_41 += 0.5 * sqrt( 6.0 * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 )
                                        * phase( TwoSL + TwoSRup + 3 ) * wigner6j
                                        * ddot_( &length, workmem2, &inc, Tdown, &inc );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_41;
   return total_35;

}

double CheMPS2::ThreeDM::diagram36_42(TensorT * denT, TensorF1 * left, TensorF0 * right, double * workmem, double * workmem2, const bool do42) const{

   assert( left->get_irrep() == right->get_irrep() );
   const int orb_i = denT->gIndex();
   
   double total_36 = 0.0;
   double total_42 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILdown = Irreps::directProd( IL,   left->get_irrep()     );
            const int IRup   = Irreps::directProd( IL,   book->gIrrep( orb_i ) );
            const int IRdown = Irreps::directProd( IRup, right->get_irrep()    );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
               
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
                  
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRup   );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRdown );
                  
                  if (( dimRup > 0 ) && ( dimRdown > 0 )){
                  
                     for ( int TwoSLdown = TwoSR-1; TwoSLdown <= TwoSR+1; TwoSLdown+=2 ){
                  
                        int dimLdown = book->gCurrentDim( orb_i, NL, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                           
                           double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,     NL+1, TwoSR, IRup   );
                           double * Tdown   =  denT->gStorage( NL,   TwoSLdown, ILdown, NL+1, TwoSR, IRdown );
                           double * Lblock  =  left->gStorage( NL,   TwoSLdown, ILdown, NL,   TwoSL, IL     );
                           double * Rblock  = right->gStorage( NL+1, TwoSR,     IRdown, NL+1, TwoSR, IRup   );
                           
                           char trans   = 'T';
                           char notrans = 'N';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           int length   = dimLdown * dimRdown;
                           int inc      = 1;
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                           dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                           const double voorfactor = 0.5 * sqrt( 6.0 * ( TwoSL + 1 ) ) * ( TwoSR + 1 )
                                                         * phase( TwoSLdown + TwoSR + 1 )
                                                         * gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSLdown, TwoSR );
                           total_36 += voorfactor * ddot_( &length, workmem2, &inc, Tdown, &inc );
                           
                           if ( do42 ){
                              double * Rblock2 = right->gStorage( NL+1, TwoSR, IRup, NL+1, TwoSR, IRdown );
                              dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock2, &dimRup, &beta, workmem2, &dimLdown );
                              total_42 += voorfactor * ddot_( &length, workmem2, &inc, Tdown, &inc );
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_42;
   return total_36;

}

double CheMPS2::ThreeDM::diagram45_46_47_48(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const{

   assert(  left->get_irrep() == right->get_irrep() );
   assert(  left->get_2j()    == right->get_2j()    ); // Spin of left and right equal
   assert(  left->get_nelec() == 2 ); //  left is Sigma
   assert( right->get_nelec() == 0 ); // right is F
   const int two_j = left->get_2j();
   const int orb_i = denT->gIndex();
   
   double total_45_46 = 0.0;
   double total_47_48 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int Idown = Irreps::directProd( IL, left->get_irrep() );
               
               for ( int TwoSLdown = TwoSL-two_j; TwoSLdown <= TwoSL+two_j; TwoSLdown+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL-2, TwoSLdown, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL,   TwoSLdown, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup    =  denT->gStorage( NL,   TwoSL,     IL,    NL, TwoSL,     IL    );
                     double * Tdown  =  denT->gStorage( NL-2, TwoSLdown, Idown, NL, TwoSLdown, Idown );
                     double * Lblock =  left->gStorage( NL-2, TwoSLdown, Idown, NL, TwoSL,     IL    );
                     double * Rblock = right->gStorage( NL,   TwoSLdown, Idown, NL, TwoSL,     IL    );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     int length   = dimLdown * dimRdown;
                     int inc      = 1;
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup,   &dimLup, &alpha, Lblock,  &dimLdown, Tup,    &dimLup,   &beta, workmem,  &dimLdown );
                     dgemm_( &notrans,   &trans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock, &dimRdown, &beta, workmem2, &dimLdown );
                     total_45_46 += phase( TwoSL - TwoSLdown + 2 )
                                  * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLdown + 1 ) )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     
                     if ( doTR ){
                        double * Rblock2 = right->gStorage( NL, TwoSL, IL, NL, TwoSLdown, Idown );
                        dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimRup, &alpha, workmem, &dimLdown, Rblock2, &dimRup, &beta, workmem2, &dimLdown );
                        total_47_48 += ( -1 ) * ( TwoSL + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_47_48;
   return total_45_46;

}

double CheMPS2::ThreeDM::diagram49_50_51_52(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const{

   assert(  left->get_irrep() == right->get_irrep() );
   assert(  left->get_2j()    == right->get_2j()    ); // Spin of left and right equal
   assert(  left->get_nelec() == 0 ); //  left is F
   assert( right->get_nelec() == 2 ); // right is Sigma
   const int two_j = left->get_2j();
   const int orb_i = denT->gIndex();
   
   double total_49_50 = 0.0;
   double total_51_52 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int Idown = Irreps::directProd( IL, left->get_irrep() );
               
               for ( int TwoSLdown = TwoSL-two_j; TwoSLdown <= TwoSL+two_j; TwoSLdown+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL,   TwoSLdown, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+2, TwoSLdown, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup    =  denT->gStorage( NL, TwoSL,     IL,    NL,   TwoSL,     IL    );
                     double * Tdown  =  denT->gStorage( NL, TwoSLdown, Idown, NL+2, TwoSLdown, Idown );
                     double * Lblock =  left->gStorage( NL, TwoSLdown, Idown, NL,   TwoSL,     IL    );
                     double * Rblock = right->gStorage( NL, TwoSL,     IL,    NL+2, TwoSLdown, Idown );
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     int length   = dimLdown * dimRdown;
                     int inc      = 1;
                     dgemm_( &notrans, &notrans, &dimLup,   &dimRdown, &dimRup, &alpha, Tup,    &dimLup,   Rblock,  &dimRup, &beta, workmem,  &dimLup   );
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRdown, &dimLup, &alpha, Lblock, &dimLdown, workmem, &dimLup, &beta, workmem2, &dimLdown );
                     total_49_50 += phase( TwoSL - TwoSLdown + 2 )
                                  * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLdown + 1 ) )
                                  * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     
                     if ( doTR ){
                        double * Lblock2 = left->gStorage( NL, TwoSL, IL, NL, TwoSLdown, Idown );
                        char trans = 'T';
                        dgemm_( &trans, &notrans, &dimLdown, &dimRdown, &dimLup, &alpha, Lblock2, &dimLup, workmem, &dimLup, &beta, workmem2, &dimLdown );
                        total_51_52 += ( -1 ) * ( TwoSLdown + 1 ) * ddot_( &length, workmem2, &inc, Tdown, &inc );
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0] = total_51_52;
   return total_49_50;

}

void CheMPS2::ThreeDM::diagram53_54(TensorT * denT, Tensor3RDM * left1, Tensor3RDM * left2, TensorL * right, double * workmem, double * workmem2) const{

   assert( left1->get_irrep()  == right->get_irrep() );
   assert( left1->get_two_j2() == 1                  );
   assert( left1->get_nelec()  == 3                  );
   if ( left2 != NULL ){
      assert( left2->get_irrep()  == right->get_irrep() );
      assert( left2->get_two_j2() == 1                  );
      assert( left2->get_nelec()  == 3                  );
   }
   const int orb_i = denT->gIndex();
   
   double total1 = 0.0;
   double total2 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int Idown = Irreps::directProd( IL, right->get_irrep() );
               
               for ( int TwoSLdown = TwoSL-1; TwoSLdown <= TwoSL+1; TwoSLdown+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL-3, TwoSLdown, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSLdown, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,    NL,   TwoSL,     IL    );
                     double * Tdown   =  denT->gStorage( NL-3, TwoSLdown, Idown, NL-1, TwoSLdown, Idown );
                     double * Rblock  = right->gStorage( NL-1, TwoSLdown, Idown, NL,   TwoSL,     IL    );
                     
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                     
                     int length = dimLup * dimLdown;
                     int inc    = 1;
                     double * Lblock1 = left1->gStorage( NL-3, TwoSLdown, Idown, NL, TwoSL, IL );
                     total1 -= ( TwoSL + 1 ) * ddot_( &length, workmem2, &inc, Lblock1, &inc );
                     if ( left2 != NULL ){
                        double * Lblock2 = left2->gStorage( NL-3, TwoSLdown, Idown, NL, TwoSL, IL );
                        total2 -= ( TwoSL + 1 ) * ddot_( &length, workmem2, &inc, Lblock2, &inc );
                     }
                  }
               }
            }
         }
      }
   }
   
                         total1 *= sqrt( 0.5 * ( left1->get_two_j1() + 1 ) );
   if ( left2 != NULL ){ total2 *= sqrt( 0.5 * ( left2->get_two_j1() + 1 ) ); }
   workmem[0]  = total1;
   workmem2[0] = total2;

}

void CheMPS2::ThreeDM::diagram55to60(TensorT * denT, Tensor3RDM ** left, TensorL * right, double * workmem, double * workmem2, double * results) const{

   for ( int cnt = 0; cnt < 6; cnt++ ){
      if ( left[ cnt ] != NULL ){
         assert( left[ cnt ]->get_irrep()  == right->get_irrep() );
         assert( left[ cnt ]->get_two_j2() == 1 );
         assert( left[ cnt ]->get_nelec()  == 1 );
      }
      results[ cnt ] = 0.0;
   }
   const int orb_i = denT->gIndex();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int Idown = Irreps::directProd( IL, right->get_irrep() );
               
               for ( int TwoSLdown = TwoSL-1; TwoSLdown <= TwoSL+1; TwoSLdown+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLdown, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLdown, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup    =  denT->gStorage( NL,   TwoSL,     IL,    NL,   TwoSL,     IL    );
                     double * Tdown  =  denT->gStorage( NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown );
                     double * Rblock = right->gStorage( NL,   TwoSL,     IL,    NL+1, TwoSLdown, Idown );
                     
                     char trans   = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &trans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRup, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans, &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup, &beta, workmem2, &dimLdown );
                     
                     int size = dimLup * dimLdown;
                     int inc  = 1;
                     const double factorBC = phase( TwoSL + 1 - TwoSLdown ) * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLdown + 1 ) );
                     const double factorD  = ( TwoSLdown + 1 );
                     for ( int cnt = 0; cnt < 6; cnt++ ){
                        if ( left[ cnt ] != NULL ){
                           results[ cnt ] += (( cnt < 4 ) ? factorBC : factorD )
                                           * ddot_( &size, workmem2, &inc, left[ cnt ]->gStorage( NL-1, TwoSLdown, Idown, NL, TwoSL, IL ), &inc );
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   const double factorJ0 = sqrt( 0.5 );
   const double factorJ1 = sqrt( 1.5 );
   for ( int cnt = 0; cnt < 6; cnt+=2 ){ if ( left[ cnt ] != NULL ){ results[ cnt ] *= factorJ0; } }
   for ( int cnt = 1; cnt < 6; cnt+=2 ){ if ( left[ cnt ] != NULL ){ results[ cnt ] *= factorJ1; } }

}

void CheMPS2::ThreeDM::diagram61_62_70_71_80_81(TensorT * denT, Tensor3RDM ** left, TensorL * right, double * workmem, double * workmem2, double * results) const{

   for ( int cnt = 0; cnt < 6; cnt++ ){
      if ( left[ cnt ] != NULL ){
         assert( left[ cnt ]->get_irrep()  == right->get_irrep() );
         assert( left[ cnt ]->get_two_j2() == 1 );
         assert( left[ cnt ]->get_nelec()  == 1 );
      }
      results[ cnt ] = 0.0;
   }
   const int orb_i = denT->gIndex();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            int dimLup = book->gCurrentDim( orb_i,   NL,   TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL+2, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               const int Idown = Irreps::directProd( IL, right->get_irrep() );
               
               for ( int TwoSLdown = TwoSL-1; TwoSLdown <= TwoSL+1; TwoSLdown+=2 ){
                  
                  int dimLdown = book->gCurrentDim( orb_i,   NL-1, TwoSLdown, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL+1, TwoSLdown, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                  
                     double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,    NL+2, TwoSL,     IL    );
                     double * Tdown   =  denT->gStorage( NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown );
                     double * Rblock  = right->gStorage( NL+1, TwoSLdown, Idown, NL+2, TwoSL,     IL    );
                     
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                     
                     int size = dimLup * dimLdown;
                     int inc  = 1;
                     const double factorBC = ( TwoSL + 1 );
                     const double factorD  = phase( TwoSL + 1 - TwoSLdown ) * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLdown + 1 ) );
                     for ( int cnt = 0; cnt < 6; cnt++ ){
                        if ( left[ cnt ] != NULL ){
                           results[ cnt ] += (( cnt < 4 ) ? factorBC : factorD )
                                           * ddot_( &size, workmem2, &inc, left[ cnt ]->gStorage( NL-1, TwoSLdown, Idown, NL, TwoSL, IL ), &inc );
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   const double factorJ0 = sqrt( 0.5 );
   const double factorJ1 = sqrt( 1.5 );
   for ( int cnt = 0; cnt < 6; cnt+=2 ){ if ( left[ cnt ] != NULL ){ results[ cnt ] *= factorJ0; } }
   for ( int cnt = 1; cnt < 6; cnt+=2 ){ if ( left[ cnt ] != NULL ){ results[ cnt ] *= factorJ1; } }

}

void CheMPS2::ThreeDM::diagramJ2half_3_2_1(TensorT * denT, Tensor3RDM ** left, TensorL * right, double * workmem, double * workmem2, double * results) const{

   for ( int cnt = 0; cnt < 6; cnt++ ){
      if ( left[ cnt ] != NULL ){
         assert( left[ cnt ]->get_irrep()  == right->get_irrep() );
         assert( left[ cnt ]->get_two_j2() == 1 );
         assert( left[ cnt ]->get_nelec()  == 1 );
      }
   }
   for ( int cnt = 0; cnt < 18; cnt++ ){ results[ cnt ] = 0.0; }
   const int orb_i = denT->gIndex();
   const double sqrt_half      = sqrt(0.5);
   const double sqrt_anderhalf = sqrt(1.5);

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            const int IRup   = Irreps::directProd( ILup, book->gIrrep( orb_i ) );
            const int ILdown = Irreps::directProd( ILup, right->get_irrep() );
            const int IRdown = Irreps::directProd( IRup, right->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
            
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
                     
                     for ( int TwoSLdown = TwoSLup-1; TwoSLdown <= TwoSLup+1; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSRup-1; TwoSRdown <= TwoSRup+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSLdown - TwoSRdown ) <= 1 )){
                              
                                 double * Tup     =  denT->gStorage( NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup   );
                                 double * Tdown   =  denT->gStorage( NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown );
                                 double * Rblock  = right->gStorage( NL,   TwoSRdown, IRdown, NL+1, TwoSRup,   IRup   );
                                 
                                 char notrans = 'N';
                                 char trans   = 'T';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );
                                 
                                 int size = dimLup * dimLdown;
                                 int inc  = 1;
                                 const double dm3_b_J0 = ((left[0] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[0]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_b_J1 = ((left[1] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[1]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_c_J0 = ((left[2] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[2]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_c_J1 = ((left[3] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[3]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_d_J0 = ((left[4] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[4]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_d_J1 = ((left[5] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[5]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 
                                 const double acommon = sqrt( 1.0 * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 );
                                 const double wign_6j = gsl_sf_coupling_6j( TwoSLup, TwoSRup, 1, TwoSRdown, TwoSLdown, 1 );
                                 const double factor1 = phase( TwoSLup + TwoSRdown + 2 ) * acommon * wign_6j;
                                 const double factor2 = sqrt( 1.0 * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * phase( TwoSLup - TwoSRdown );
                                 const double factor3 = 2 * acommon
                                                      * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSLdown, TwoSRdown )
                                                      * gsl_sf_coupling_6j( 1, 1, 2, TwoSRup, TwoSLdown, TwoSLup   );
                                 const double wigprod = gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSRdown, TwoSRup   )
                                                      * gsl_sf_coupling_6j( 1, 1, 2, TwoSLup, TwoSRdown, TwoSLdown );
                                 const double factor4 = 2 * acommon * wigprod
                                                      * phase( TwoSRup + TwoSRdown + TwoSLup + TwoSLdown );
                                 const double wign_9j = gsl_sf_coupling_9j( 1, 1, 2, TwoSLdown, TwoSLup, 1, TwoSRdown, TwoSRup, 1 );
                                 const double factor5 = ( -2 ) * acommon * wign_9j;
                                 
                                 const double dcommon = sqrt( 1.0 * ( TwoSLdown + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 );
                                 const double factor6 = dcommon * phase( TwoSLdown + TwoSRdown + 3 ) * wign_6j;
                                 const double factor7 = sqrt( ( TwoSLdown + 1.0 ) / ( TwoSLup + 1.0 ) ) * ( TwoSRup + 1.0 ) * phase( TwoSLup + 1 - TwoSLdown );
                                 const double factor8 = 2 * dcommon * wigprod * phase( TwoSRdown + 1 - TwoSRup );
                                 const double factor9 = 2 * dcommon * wign_9j * phase( TwoSRdown + 1 - TwoSRup )
                                                      * phase( TwoSLup + TwoSLdown + TwoSRup + TwoSRdown + 2 ); //Last phase for column swap wigner 9j
                                 
                                 
                                 results[0]  += sqrt_half      * factor1         * dm3_b_J0;   // d63
                                 results[1]  += sqrt_anderhalf * factor1         * dm3_b_J1;   // d64
    if ( TwoSRdown == TwoSLup ){ results[2]  -= sqrt_half      * ( TwoSRup + 1 ) * dm3_b_J0;   // d65
                                 results[3]  -= sqrt_anderhalf * ( TwoSRup + 1 ) * dm3_b_J1; } // d66
    if ( TwoSRup == TwoSLdown ){ results[4]  += sqrt_half      * factor2         * dm3_b_J0; } // d67
                                 results[5]  += sqrt_anderhalf * factor3         * dm3_b_J1;   // d68
                                 
                                 results[6]  -= sqrt_half      * factor1         * dm3_c_J0;   // d72
                                 results[7]  -= sqrt_anderhalf * factor1         * dm3_c_J1;   // d73
    if ( TwoSRdown == TwoSLup ){ results[8]  += sqrt_half      * ( TwoSRup + 1 ) * dm3_c_J0;   // d74
                                 results[9]  += sqrt_anderhalf * ( TwoSRup + 1 ) * dm3_c_J1; } // d75
                                 results[10] += sqrt_anderhalf * factor4         * dm3_c_J1;   // d76
                                 results[11] += sqrt_anderhalf * factor5         * dm3_c_J1;   // d77
                                 
                                 results[12] += sqrt_half      * factor6         * dm3_d_J0;   // d82
                                 results[13] += sqrt_anderhalf * factor6         * dm3_d_J1;   // d83
    if ( TwoSRdown == TwoSLup ){ results[14] += sqrt_half      * factor7         * dm3_d_J0;   // d84
                                 results[15] += sqrt_anderhalf * factor7         * dm3_d_J1; } // d85
                                 results[16] += sqrt_anderhalf * factor8         * dm3_d_J1;   // d86
                                 results[17] += sqrt_anderhalf * factor9         * dm3_d_J1;   // d87
                                 
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

}

void CheMPS2::ThreeDM::diagramJ2oneandhalf_3_2_1(TensorT * denT, Tensor3RDM ** left, TensorL * right, double * workmem, double * workmem2, double * results) const{

   for ( int cnt = 0; cnt < 3; cnt++ ){
      if ( left[ cnt ] != NULL ){
         assert( left[ cnt ]->get_irrep()  == right->get_irrep() );
         assert( left[ cnt ]->get_two_j2() == 3 );
         assert( left[ cnt ]->get_nelec()  == 1 );
      }
   }
   for ( int cnt = 0; cnt < 5; cnt++ ){ results[ cnt ] = 0.0; }
   const int orb_i = denT->gIndex();

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSLup = book->gTwoSmin( orb_i, NL ); TwoSLup <= book->gTwoSmax( orb_i, NL ); TwoSLup+=2 ){
         for ( int ILup = 0; ILup < book->getNumberOfIrreps(); ILup++ ){
         
            const int IRup   = Irreps::directProd( ILup, book->gIrrep( orb_i ) );
            const int ILdown = Irreps::directProd( ILup, right->get_irrep() );
            const int IRdown = Irreps::directProd( IRup, right->get_irrep() );
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSLup, ILup );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSRup = TwoSLup-1; TwoSRup <= TwoSLup+1; TwoSRup+=2 ){
            
                  int dimRup = book->gCurrentDim( orb_i+1, NL+1, TwoSRup, IRup );
                  
                  if ( dimRup > 0 ){
                     
                     for ( int TwoSLdown = TwoSLup-3; TwoSLdown <= TwoSLup+3; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL-1, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                        
                           for ( int TwoSRdown = TwoSRup-1; TwoSRdown <= TwoSRup+1; TwoSRdown+=2 ){
                        
                              int dimRdown = book->gCurrentDim( orb_i+1, NL, TwoSRdown, IRdown );
                              
                              if (( dimRdown > 0 ) && ( abs( TwoSLdown - TwoSRdown ) <= 1 )){
                              
                                 double * Tup     =  denT->gStorage( NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup   );
                                 double * Tdown   =  denT->gStorage( NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown );
                                 double * Rblock  = right->gStorage( NL,   TwoSRdown, IRdown, NL+1, TwoSRup,   IRup   );
                                 
                                 char notrans = 'N';
                                 char trans   = 'T';
                                 double alpha = 1.0;
                                 double beta  = 0.0; //set
                                 dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                                 dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );

                                 int size = dimLup * dimLdown;
                                 int inc  = 1;
                                 const double dm3_b = ((left[0] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[0]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_c = ((left[1] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[1]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));
                                 const double dm3_d = ((left[2] == NULL) ? 0.0 : ddot_(&size, workmem2, &inc, left[2]->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup), &inc));

                                 const double acommon = 2 * sqrt( 3.0 * ( TwoSLup + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 );
                                 const double wign_9j = gsl_sf_coupling_9j( 1, 1, 2, TwoSLdown, TwoSLup, 3, TwoSRdown, TwoSRup, 1 );
                                 const double factor1 = acommon * gsl_sf_coupling_6j( 1, 2, 1, TwoSLdown, TwoSRdown, TwoSRup )
                                                                * gsl_sf_coupling_6j( 1, 2, 3, TwoSLdown, TwoSLup,   TwoSRup );
                                 const double wigprod = gsl_sf_coupling_6j( 1, 2, 1, TwoSLup, TwoSRup,   TwoSRdown )
                                                      * gsl_sf_coupling_6j( 1, 2, 3, TwoSLup, TwoSLdown, TwoSRdown );
                                 const double factor2 = acommon * wigprod * phase( TwoSRup + TwoSRdown + TwoSLup + TwoSLdown + 2 );
                                 const double factor3 = acommon * ( -1 ) * wign_9j;
                                 const double dcommon = 2 * sqrt( 3.0 * ( TwoSLdown + 1 ) * ( TwoSRdown + 1 ) ) * ( TwoSRup + 1 );
                                 const double factor4 = dcommon * wigprod * phase( TwoSRdown + 1 - TwoSRup );
                                 const double factor5 = dcommon * wign_9j * phase( TwoSRdown + 1 - TwoSRup )
                                                      * phase( TwoSRup + TwoSRdown + TwoSLup + TwoSLdown ); //Last phase for column swap wigner 9j

                                 results[0] += factor1 * dm3_b; // d69
                                 results[1] += factor2 * dm3_c; // d78
                                 results[2] += factor3 * dm3_c; // d79
                                 results[3] += factor4 * dm3_d; // d88
                                 results[4] += factor5 * dm3_d; // d89
                                 
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

}

void CheMPS2::ThreeDM::diagram90_93(TensorT * denT, Tensor3RDM ** left, TensorS0 * denS0, double * workmem, double * workmem2) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS0->get_irrep(), book->gIrrep( orb_i ) );
   if ( left[0] != NULL ){ assert( left[0]->get_irrep() == ImxInxIi ); assert( left[0]->get_nelec() == 3 ); assert( left[0]->get_two_j1()==0 ); assert( left[0]->get_two_j2()==1 ); }
   if ( left[1] != NULL ){ assert( left[1]->get_irrep() == ImxInxIi ); assert( left[1]->get_nelec() == 3 ); assert( left[1]->get_two_j1()==2 ); assert( left[1]->get_two_j2()==1 ); }

   double total_90 = 0.0;
   double total_93 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){

            const int ILdown = Irreps::directProd( IL, ImxInxIi           ); // IL x Im x In x Ii
            const int IRdown = Irreps::directProd( IL, denS0->get_irrep() ); // IL x Im x In

            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );

            if (( dimLup > 0 ) && ( dimRup > 0 )){

               for ( int TwoSLdown = TwoSL-1; TwoSLdown <= TwoSL+1; TwoSLdown+=2 ){

                  int dimLdown = book->gCurrentDim( orb_i,   NL-3, TwoSLdown, ILdown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-2, TwoSL,     IRdown );

                  if (( dimLdown > 0 ) && ( dimRdown > 0 )){

                     double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,     NL,   TwoSL,     IL     );
                     double * Tdown   =  denT->gStorage( NL-3, TwoSLdown, ILdown, NL-2, TwoSL,     IRdown );
                     double * Rblock  = denS0->gStorage( NL-2, TwoSL,     IRdown, NL,   TwoSL,     IL     );

                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );

                     int size = dimLup * dimLdown;
                     int inc  = 1;
                     if ( left[0] != NULL ){
                        double * Lblock = left[0]->gStorage( NL-3, TwoSLdown, ILdown, NL, TwoSL, IL );
                        total_90 -= ( TwoSL + 1 ) * ddot_(&size, workmem2, &inc, Lblock, &inc);
                     }
                     if ( left[1] != NULL ){
                        double * Lblock = left[1]->gStorage( NL-3, TwoSLdown, ILdown, NL, TwoSL, IL );
                        total_93 -= ( TwoSL + 1 ) * ddot_(&size, workmem2, &inc, Lblock, &inc);
                     }
                  }
               }
            }
         }
      }
   }

   workmem[0]  = total_90 * 0.5;
   workmem2[0] = total_93 * 0.5 * sqrt( 3.0 );

}

void CheMPS2::ThreeDM::diagram91_92_94(TensorT * denT, Tensor3RDM ** left, TensorS1 * denS1, double * workmem, double * workmem2, double * result) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS1->get_irrep(), book->gIrrep( orb_i ) );
   if ( left[0] != NULL ){ assert( left[0]->get_irrep() == ImxInxIi ); assert( left[0]->get_nelec() == 3 ); assert( left[0]->get_two_j1()==0 ); assert( left[0]->get_two_j2()==1 ); }
   if ( left[1] != NULL ){ assert( left[1]->get_irrep() == ImxInxIi ); assert( left[1]->get_nelec() == 3 ); assert( left[1]->get_two_j1()==2 ); assert( left[1]->get_two_j2()==1 ); }
   if ( left[2] != NULL ){ assert( left[2]->get_irrep() == ImxInxIi ); assert( left[2]->get_nelec() == 3 ); assert( left[2]->get_two_j1()==2 ); assert( left[2]->get_two_j2()==3 ); }
   
   double total_91 = 0.0;
   double total_92 = 0.0;
   double total_94 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int ILdown = Irreps::directProd( IL, ImxInxIi           ); // IL x Im x In x Ii
            const int IRdown = Irreps::directProd( IL, denS1->get_irrep() ); // IL x Im x In
            
            int dimLup = book->gCurrentDim( orb_i,   NL, TwoSL, IL );
            int dimRup = book->gCurrentDim( orb_i+1, NL, TwoSL, IL );
            
            if (( dimLup > 0 ) && ( dimRup > 0 )){
            
               for ( int TwoSRdown = TwoSL-2; TwoSRdown <= TwoSL+2; TwoSRdown+=2 ){
               
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-2, TwoSRdown, IRdown );
                  
                  if ( dimRdown > 0 ){

                     for ( int TwoSLdown = TwoSRdown-1; TwoSLdown <= TwoSRdown+1; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i, NL-3, TwoSLdown, ILdown );
                        
                        if ( dimLdown > 0 ){
                           
                           double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,     NL,   TwoSL,     IL     );
                           double * Tdown   =  denT->gStorage( NL-3, TwoSLdown, ILdown, NL-2, TwoSRdown, IRdown );
                           double * Rblock  = denS1->gStorage( NL-2, TwoSRdown, IRdown, NL,   TwoSL,     IL     );
                           
                           char notrans = 'N';
                           char trans   = 'T';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );

                           int size = dimLup * dimLdown;
                           int inc  = 1;
                           const double doublet_6j = gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSRdown, TwoSLdown );
                           const double common     = sqrt( TwoSRdown + 1.0 ) * ( TwoSL + 1 ) * phase( TwoSL + TwoSLdown + 1 );
                           if (( abs( TwoSL - TwoSLdown ) == 1 ) && ( left[1] != NULL )){
                              double * Lblock = left[1]->gStorage( NL-3, TwoSLdown, ILdown, NL, TwoSL, IL );
                              total_91 += common * doublet_6j * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                           if ( left[2] != NULL ){
                              double * Lblock = left[2]->gStorage( NL-3, TwoSLdown, ILdown, NL, TwoSL, IL );
                              total_92 += common * gsl_sf_coupling_6j( 1, 3, 2, TwoSL, TwoSRdown, TwoSLdown ) * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                           if (( abs( TwoSL - TwoSLdown ) == 1 ) && ( left[0] != NULL )){
                              double * Lblock = left[0]->gStorage( NL-3, TwoSLdown, ILdown, NL, TwoSL, IL );
                              total_94 += common * doublet_6j * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   result[0] = total_91 * sqrt( 0.5 );
   result[1] = total_92 * ( -2 );
   result[2] = total_94 * ( -sqrt( 1.5 ) );

}

void CheMPS2::ThreeDM::diagram95_98(TensorT * denT, Tensor3RDM ** left, TensorS0 * denS0, double * workmem, double * workmem2) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS0->get_irrep(), book->gIrrep( orb_i ) );
   if ( left[0] != NULL ){ assert( left[0]->get_irrep() == ImxInxIi ); assert( left[0]->get_nelec() == 3 ); assert( left[0]->get_two_j1()==0 ); assert( left[0]->get_two_j2()==1 ); }
   if ( left[1] != NULL ){ assert( left[1]->get_irrep() == ImxInxIi ); assert( left[1]->get_nelec() == 3 ); assert( left[1]->get_two_j1()==2 ); assert( left[1]->get_two_j2()==1 ); }
   
   double total_95 = 0.0;
   double total_98 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int IRup  = Irreps::directProd( IL, book->gIrrep( orb_i ) ); // IL x Ii
            const int Idown = Irreps::directProd( IL, ImxInxIi              ); // IL x Im x In x Ii
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
            
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRup  );
                  int dimLdown = book->gCurrentDim( orb_i,   NL-3, TwoSR, Idown );
                  int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSR, Idown );
                  
                  if (( dimLdown > 0 ) && ( dimRdown > 0 ) && ( dimRup > 0 )){
                        
                     double * Tup     =  denT->gStorage( NL,   TwoSL, IL,    NL+1, TwoSR, IRup  );
                     double * Tdown   =  denT->gStorage( NL-3, TwoSR, Idown, NL-1, TwoSR, Idown );
                     double * Rblock  = denS0->gStorage( NL-1, TwoSR, Idown, NL+1, TwoSR, IRup  );
                     
                     char notrans = 'N';
                     char trans   = 'T';
                     double alpha = 1.0;
                     double beta  = 0.0; //set
                     dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                     dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );

                     int size = dimLup * dimLdown;
                     int inc  = 1;
                     const double common = sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSR + 1 ) ) * phase( TwoSL + 1 - TwoSR );
                     if ( left[0] != NULL ){
                        double * Lblock = left[0]->gStorage( NL-3, TwoSR, Idown, NL, TwoSL, IL );
                        total_95 += common * ddot_(&size, workmem2, &inc, Lblock, &inc);
                     }
                     if ( left[1] != NULL ){
                        double * Lblock = left[1]->gStorage( NL-3, TwoSR, Idown, NL, TwoSL, IL );
                        total_98 += common * ddot_(&size, workmem2, &inc, Lblock, &inc);
                     }
                  }
               }
            }
         }
      }
   }
   
   workmem[0]  = total_95 * 0.5;
   workmem2[0] = total_98 * 0.5 * sqrt( 3.0 );

}

void CheMPS2::ThreeDM::diagram96_97_99(TensorT * denT, Tensor3RDM ** left, TensorS1 * denS1, double * workmem, double * workmem2, double * result) const{

   const int orb_i     = denT->gIndex();
   const int ImxInxIi  = Irreps::directProd( denS1->get_irrep(), book->gIrrep( orb_i ) );
   if ( left[0] != NULL ){ assert( left[0]->get_irrep() == ImxInxIi ); assert( left[0]->get_nelec() == 3 ); assert( left[0]->get_two_j1()==0 ); assert( left[0]->get_two_j2()==1 ); }
   if ( left[1] != NULL ){ assert( left[1]->get_irrep() == ImxInxIi ); assert( left[1]->get_nelec() == 3 ); assert( left[1]->get_two_j1()==2 ); assert( left[1]->get_two_j2()==1 ); }
   if ( left[2] != NULL ){ assert( left[2]->get_irrep() == ImxInxIi ); assert( left[2]->get_nelec() == 3 ); assert( left[2]->get_two_j1()==2 ); assert( left[2]->get_two_j2()==3 ); }
   
   double total_96 = 0.0;
   double total_97 = 0.0;
   double total_99 = 0.0;

   for ( int NL = book->gNmin( orb_i ); NL <= book->gNmax( orb_i ); NL++ ){
      for ( int TwoSL = book->gTwoSmin( orb_i, NL ); TwoSL <= book->gTwoSmax( orb_i, NL ); TwoSL+=2 ){
         for ( int IL = 0; IL < book->getNumberOfIrreps(); IL++ ){
         
            const int IRup  = Irreps::directProd( IL, book->gIrrep( orb_i ) ); // IL x Ii
            const int Idown = Irreps::directProd( IL, ImxInxIi              ); // IL x Im x In x Ii
            
            int dimLup = book->gCurrentDim( orb_i, NL, TwoSL, IL );
            
            if ( dimLup > 0 ){
            
               for ( int TwoSR = TwoSL-1; TwoSR <= TwoSL+1; TwoSR+=2 ){
            
                  int dimRup   = book->gCurrentDim( orb_i+1, NL+1, TwoSR, IRup  );
            
                  if ( dimRup > 0 ){
                        
                     for ( int TwoSLdown = TwoSR-2; TwoSLdown <= TwoSR+2; TwoSLdown+=2 ){
                        
                        int dimLdown = book->gCurrentDim( orb_i,   NL-3, TwoSLdown, Idown );
                        int dimRdown = book->gCurrentDim( orb_i+1, NL-1, TwoSLdown, Idown );
                        
                        if (( dimLdown > 0 ) && ( dimRdown > 0 )){
                              
                           double * Tup     =  denT->gStorage( NL,   TwoSL,     IL,    NL+1, TwoSR,     IRup  );
                           double * Tdown   =  denT->gStorage( NL-3, TwoSLdown, Idown, NL-1, TwoSLdown, Idown );
                           double * Rblock  = denS1->gStorage( NL-1, TwoSLdown, Idown, NL+1, TwoSR,     IRup  );
                           
                           char notrans = 'N';
                           char trans   = 'T';
                           double alpha = 1.0;
                           double beta  = 0.0; //set
                           dgemm_( &notrans, &notrans, &dimLdown, &dimRup, &dimRdown, &alpha, Tdown,   &dimLdown, Rblock, &dimRdown, &beta, workmem,  &dimLdown );
                           dgemm_( &notrans, &trans,   &dimLdown, &dimLup, &dimRup,   &alpha, workmem, &dimLdown, Tup,    &dimLup,   &beta, workmem2, &dimLdown );

                           int size = dimLup * dimLdown;
                           int inc  = 1;
                           const double doublet_6j = gsl_sf_coupling_6j( 1, 1, 2, TwoSLdown, TwoSR, TwoSL );
                           const double common     = sqrt( TwoSL + 1.0 ) * ( TwoSR + 1 ) * phase( TwoSR + TwoSLdown );
                           if (( abs( TwoSL - TwoSLdown ) == 1 ) && ( left[1] != NULL )){
                              double * Lblock = left[1]->gStorage( NL-3, TwoSLdown, Idown, NL, TwoSL, IL );
                              total_96 += common * doublet_6j * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                           if ( left[2] != NULL ){
                              double * Lblock = left[2]->gStorage( NL-3, TwoSLdown, Idown, NL, TwoSL, IL );
                              total_97 += common * gsl_sf_coupling_6j( 1, 3, 2, TwoSLdown, TwoSR, TwoSL ) * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                           if (( abs( TwoSL - TwoSLdown ) == 1 ) && ( left[0] != NULL )){
                              double * Lblock = left[0]->gStorage( NL-3, TwoSLdown, Idown, NL, TwoSL, IL );
                              total_99 += common * doublet_6j * ddot_(&size, workmem2, &inc, Lblock, &inc);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   result[0] = total_96 * sqrt( 0.5 );
   result[1] = total_97 * 2;
   result[2] = total_99 * ( -sqrt( 1.5 ) );

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
                     double alpha = 0.5 * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSLprime + 1 ) ) * phase( TwoSL + 1 - TwoSLprime );
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
                                            * phase( TwoSL + TwoSRprime )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSL, TwoSRprime, TwoSLprime );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( 1.0 * ( TwoSLprime + 1 ) ) * ( TwoSRprime + 1 )
                                            * phase( TwoSL + TwoSRprime )
                                            * gsl_sf_coupling_6j( 1, 2, 3, TwoSL, TwoSLprime, TwoSRprime );
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
                                            * phase( TwoSL + TwoSLprime + 3 )
                                            * gsl_sf_coupling_6j( 1, 1, 2, TwoSLprime, TwoSR, TwoSL );
                           int length = dimLup * dimLdown;
                           int inc = 1;
                           daxpy_( &length, &prefactor, workmem2, &inc, Wblock, &inc );
                        
                        }
                        {
                        
                           double * Wblock  = quartet->gStorage( NL, TwoSL, IL, NL+1, TwoSLprime, ILxImxInxIi );
                           double prefactor = sqrt( 1.0 * ( TwoSR + 1 ) ) * ( TwoSLprime + 1 )
                                            * phase( TwoSL + TwoSLprime + 1 )
                                            * gsl_sf_coupling_6j( 1, 3, 2, TwoSLprime, TwoSR, TwoSL );
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


