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

void CheMPS2::ThreeDM::fill_site(TensorT * denT, TensorL *** Ltensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors){

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   const int orb_i = denT->gIndex();
   const int DIM = max(book->gMaxDimAtBound( orb_i ), book->gMaxDimAtBound( orb_i + 1 ));
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




