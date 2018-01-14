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
#include <sys/time.h>
#include <assert.h>

#include "DMRG.h"
#include "MPIchemps2.h"
#include "Special.h"

void CheMPS2::DMRG::update_safe_3rdm_operators(const int boundary){

   /*
      indices 0 <= j <= k <= l < boundary
      
      tensor_3rdm[ boundary - 1 == index ][ k - j ][ l - k ][ boundary - 1 - l ]
      
      **************************
      *   anni / anni / anni   *
      **************************

         1/ j == k == l is forbidden ( no three annihilators on the same site )
         2/ if j == k, J1 must be zero ( Sigma_J1 does not exist then )
         3/ if k == l, J2 must be 1/2
         
         tensor_3rdm_a_J0_doublet  -->  j <= k <  l  &  j < k == l  ( or NOT j == k == l )
         tensor_3rdm_a_J1_doublet  -->  j <  k <= l               
         tensor_3rdm_a_J1_quartet  -->  j <  k <  l

      **************************
      *   anni / anni / crea   *
      **************************

         1/ j <= k < l is allowed ( k == l is part of tensor_3rdm_c )
         2/ if j == k, J1 must be zero ( Sigma_J1 does not exist then )
         
         tensor_3rdm_b_J0_doublet  -->  j <= k < l
         tensor_3rdm_b_J1_doublet  -->  j <  k < l
         tensor_3rdm_b_J1_quartet  -->  j <  k < l

      **************************
      *   anni / crea / anni   *
      **************************

         1/ j == k == l is forbidden ( j == k == l is part of tensor_3rdm_d )

         tensor_3rdm_c_J0_doublet  -->  j <= k < l  &  j < k == l  ( or NOT j == k == l )
         tensor_3rdm_c_J1_doublet  -->  j <= k < l  &  j < k == l  ( or NOT j == k == l )
         tensor_3rdm_c_J1_quartet  -->  j <= k < l  &  j < k == l  ( or NOT j == k == l )

      **************************
      *   crea / anni / anni   *
      **************************

         1/ j <= k <= l is allowed
         2/ if k == l, J2 must be 1/2
         
         tensor_3rdm_d_J0_doublet  -->  j <= k <= l
         tensor_3rdm_d_J1_doublet  -->  j <= k <= l
         tensor_3rdm_d_J1_quartet  -->  j <= k <  l
   */

   allocate_3rdm_operators( boundary );
   update_3rdm_operators( boundary );
   if ( boundary >= 2 ){ delete_3rdm_operators( boundary - 1 ); }
   
}

void CheMPS2::DMRG::update_3rdm_operators(const int boundary){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const int index = boundary - 1;
   const int dimL  = denBK->gMaxDimAtBound(boundary-1);
   const int dimR  = denBK->gMaxDimAtBound(boundary);
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   #pragma omp parallel
   {
      double * workmem = new double[dimL*dimR];
   
#ifdef CHEMPS2_MPI_COMPILATION //######( loop j<=k<=l MPI )######//

/* Strategy for MPI:
     - outer loop is (j,k)
     - everyone has a temporary duplicate of S_jk and F_jk
*/
      for ( int orb_j = 0; orb_j < boundary; orb_j++ ){
         for ( int orb_k = orb_j; orb_k < boundary; orb_k++ ){
            
            const int irrjk = Irreps::directProd( denBK->gIrrep( orb_j ), denBK->gIrrep( orb_k ) );
            const int cnt1  = orb_k - orb_j;
            
            #pragma omp single
            if ( orb_k < index-1 ){ // All processes own Fx/Sx[ index - 1 ][ k - j ][ index - 1 - k == 0 ]
               const int own_S_jk = MPIchemps2::owner_absigma( orb_j, orb_k );
               const int own_F_jk = MPIchemps2::owner_cdf(  L, orb_j, orb_k );
               if ( MPIRANK != own_F_jk ){ F0tensors[index-1][cnt1][index-orb_k-1] = new TensorF0( index, irrjk, true, denBK );
                                           F1tensors[index-1][cnt1][index-orb_k-1] = new TensorF1( index, irrjk, true, denBK ); }
               if ( MPIRANK != own_S_jk ){ S0tensors[index-1][cnt1][index-orb_k-1] = new TensorS0( index, irrjk, true, denBK );
                          if ( cnt1 > 0 ){ S1tensors[index-1][cnt1][index-orb_k-1] = new TensorS1( index, irrjk, true, denBK ); }}
                                MPIchemps2::broadcast_tensor( F0tensors[index-1][cnt1][index-orb_k-1], own_F_jk );
                                MPIchemps2::broadcast_tensor( F1tensors[index-1][cnt1][index-orb_k-1], own_F_jk );
                                MPIchemps2::broadcast_tensor( S0tensors[index-1][cnt1][index-orb_k-1], own_S_jk );
               if ( cnt1 > 0 ){ MPIchemps2::broadcast_tensor( S1tensors[index-1][cnt1][index-orb_k-1], own_S_jk ); }
            }
         
            #pragma omp for schedule(dynamic)
            for ( int orb_l = orb_k; orb_l < boundary; orb_l++ ){
               if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK ){
                  const int cnt2 = orb_l - orb_k;
                  const int cnt3 = index - orb_l;

#else //######( loop j<=k<=l MPI )######//

      const int upperbound = (boundary*(boundary+1)*(boundary+2))/6;
      int jkl[] = { 0, 0, 0 };
      #pragma omp for schedule(static)
      for ( int global = 0; global < upperbound; global++ ){
         Special::invert_triangle_three( global, jkl );
         const int orb_j = jkl[ 0 ];
         const int orb_k = jkl[ 1 ];
         const int orb_l = jkl[ 2 ];
         const int recalculate_global = orb_j + (orb_k*(orb_k+1))/2 + (orb_l*(orb_l+1)*(orb_l+2))/6;
         assert( global == recalculate_global );
         const int cnt1 = orb_k - orb_j;
         const int cnt2 = orb_l - orb_k;
         const int cnt3 = index - orb_l;
         
#endif //######( loop j<=k<=l MPI )######//

         /* PERFORM THE UPDATES */
         if ( cnt3 == 0 ){ // Create tensors
            if ( cnt2 > 0 ){
                            tensor_3rdm_a_J0_doublet[index][cnt1][cnt2][0]->a1(S0tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
               if (cnt1>0){ tensor_3rdm_a_J1_doublet[index][cnt1][cnt2][0]->a1(S1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_a_J1_quartet[index][cnt1][cnt2][0]->a1(S1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem); }
                            tensor_3rdm_b_J0_doublet[index][cnt1][cnt2][0]->b1(S0tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
               if (cnt1>0){ tensor_3rdm_b_J1_doublet[index][cnt1][cnt2][0]->b1(S1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_b_J1_quartet[index][cnt1][cnt2][0]->b1(S1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem); }
                            tensor_3rdm_c_J0_doublet[index][cnt1][cnt2][0]->c1(F0tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_c_J1_doublet[index][cnt1][cnt2][0]->c1(F1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_c_J1_quartet[index][cnt1][cnt2][0]->c1(F1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_d_J0_doublet[index][cnt1][cnt2][0]->d1(F0tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_d_J1_doublet[index][cnt1][cnt2][0]->d1(F1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem);
                            tensor_3rdm_d_J1_quartet[index][cnt1][cnt2][0]->d1(F1tensors[index-1][cnt1][cnt2-1], MPS[index], workmem); 
            } else {
               if ( cnt1 > 0 ){
                  tensor_3rdm_a_J0_doublet[index][cnt1][0][0]->extra2(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_a_J1_doublet[index][cnt1][0][0]->extra2(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_c_J0_doublet[index][cnt1][0][0]->extra4(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_c_J1_doublet[index][cnt1][0][0]->extra4(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_c_J1_quartet[index][cnt1][0][0]->extra4(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_d_J0_doublet[index][cnt1][0][0]->extra3(Ltensors[index-1][cnt1-1], MPS[index], workmem);
                  tensor_3rdm_d_J1_doublet[index][cnt1][0][0]->extra3(Ltensors[index-1][cnt1-1], MPS[index], workmem);
               } else {
                  tensor_3rdm_d_J0_doublet[index][0][0][0]->extra1(MPS[index]);
                  tensor_3rdm_d_J1_doublet[index][0][0][0]->extra1(MPS[index]);
               }
            }
         } else { // Update tensors
            if (cnt1+cnt2>0){ tensor_3rdm_a_J0_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_a_J0_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1>0)     { tensor_3rdm_a_J1_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_a_J1_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1*cnt2>0){ tensor_3rdm_a_J1_quartet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_a_J1_quartet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt2>0)     { tensor_3rdm_b_J0_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_b_J0_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1*cnt2>0){ tensor_3rdm_b_J1_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_b_J1_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1*cnt2>0){ tensor_3rdm_b_J1_quartet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_b_J1_quartet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1+cnt2>0){ tensor_3rdm_c_J0_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_c_J0_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1+cnt2>0){ tensor_3rdm_c_J1_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_c_J1_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
            if (cnt1+cnt2>0){ tensor_3rdm_c_J1_quartet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_c_J1_quartet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
                              tensor_3rdm_d_J0_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_d_J0_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem);
                              tensor_3rdm_d_J1_doublet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_d_J1_doublet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem);
            if (cnt2>0)     { tensor_3rdm_d_J1_quartet[index][cnt1][cnt2][cnt3]->update(tensor_3rdm_d_J1_quartet[index-1][cnt1][cnt2][cnt3-1], MPS[index], MPS[index], workmem); }
         }

#ifdef CHEMPS2_MPI_COMPILATION //######( close loop j<=k<=l MPI )######//

               }
            }

            #pragma omp single
            if ( orb_k < index - 1 ){ // All processes own Fx/Sx[ index - 1 ][ k - j ][ index - 1 - k == 0 ]
               const int own_S_jk = MPIchemps2::owner_absigma( orb_j, orb_k );
               const int own_F_jk = MPIchemps2::owner_cdf(  L, orb_j, orb_k );
               if ( MPIRANK != own_F_jk ){ delete F0tensors[index-1][cnt1][index-orb_k-1]; F0tensors[index-1][cnt1][index-orb_k-1] = NULL;
                                           delete F1tensors[index-1][cnt1][index-orb_k-1]; F1tensors[index-1][cnt1][index-orb_k-1] = NULL; }
               if ( MPIRANK != own_S_jk ){ delete S0tensors[index-1][cnt1][index-orb_k-1]; S0tensors[index-1][cnt1][index-orb_k-1] = NULL;
                          if ( cnt1 > 0 ){ delete S1tensors[index-1][cnt1][index-orb_k-1]; S1tensors[index-1][cnt1][index-orb_k-1] = NULL; }}
            }
         }
      }
      
#else //######( close loop j<=k<=l MPI )######//

      }
      
#endif //######( close loop j<=k<=l MPI )######//

      delete [] workmem;
   }
   
   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_CALC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::allocate_3rdm_operators(const int boundary){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   const int index = boundary - 1;
   
   tensor_3rdm_a_J0_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_a_J1_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_a_J1_quartet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_b_J0_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_b_J1_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_b_J1_quartet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_c_J0_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_c_J1_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_c_J1_quartet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_d_J0_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_d_J1_doublet[ index ] = new Tensor3RDM***[ boundary ];
   tensor_3rdm_d_J1_quartet[ index ] = new Tensor3RDM***[ boundary ];
   
   for ( int cnt1 = 0; cnt1 < boundary; cnt1++ ){ // cnt1 = k - j < boundary
   
      tensor_3rdm_a_J0_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_a_J1_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_a_J1_quartet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_b_J0_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_b_J1_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_b_J1_quartet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_c_J0_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_c_J1_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_c_J1_quartet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_d_J0_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_d_J1_doublet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      tensor_3rdm_d_J1_quartet[ index ][ cnt1 ] = new Tensor3RDM**[ boundary - cnt1 ];
      
      for ( int cnt2 = 0; cnt2 < boundary - cnt1; cnt2++ ){ // cnt2 = l - k < boundary - k = boundary - cnt1 - j <= boundary - cnt1
      
         tensor_3rdm_a_J0_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_a_J1_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_a_J1_quartet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_b_J0_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_b_J1_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_b_J1_quartet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_c_J0_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_c_J1_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_c_J1_quartet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_d_J0_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_d_J1_doublet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         tensor_3rdm_d_J1_quartet[ index ][ cnt1 ][ cnt2 ] = new Tensor3RDM*[ boundary - cnt1 - cnt2 ];
         
         for ( int cnt3 = 0; cnt3 < boundary - cnt1 - cnt2; cnt3++ ){ // cnt3 = boundary - 1 - l < boundary - ( l - k ) - ( k - j ) = boundary - cnt1 - cnt2
            
            const int orb_l = boundary - 1 - cnt3;
            const int orb_k = orb_l - cnt2;
            const int orb_j = orb_k - cnt1;
            const int irr   = Irreps::directProd( Irreps::directProd( denBK->gIrrep( orb_j ), denBK->gIrrep( orb_k ) ), denBK->gIrrep( orb_l ) );
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK ){
            #endif
               tensor_3rdm_a_J0_doublet[index][cnt1][cnt2][cnt3] = (cnt1+cnt2>0) ? new Tensor3RDM(boundary, 0, 1, 3, irr, true,  denBK) : NULL; // NOT j == k == l
               tensor_3rdm_a_J1_doublet[index][cnt1][cnt2][cnt3] = (cnt1>0)      ? new Tensor3RDM(boundary, 2, 1, 3, irr, true,  denBK) : NULL; //     j <  k <= l
               tensor_3rdm_a_J1_quartet[index][cnt1][cnt2][cnt3] = (cnt1*cnt2>0) ? new Tensor3RDM(boundary, 2, 3, 3, irr, true,  denBK) : NULL; //     j <  k <  l
               tensor_3rdm_b_J0_doublet[index][cnt1][cnt2][cnt3] = (cnt2>0)      ? new Tensor3RDM(boundary, 0, 1, 1, irr, true,  denBK) : NULL; //     j <= k <  l
               tensor_3rdm_b_J1_doublet[index][cnt1][cnt2][cnt3] = (cnt1*cnt2>0) ? new Tensor3RDM(boundary, 2, 1, 1, irr, true,  denBK) : NULL; //     j <  k <  l
               tensor_3rdm_b_J1_quartet[index][cnt1][cnt2][cnt3] = (cnt1*cnt2>0) ? new Tensor3RDM(boundary, 2, 3, 1, irr, true,  denBK) : NULL; //     j <  k <  l
               tensor_3rdm_c_J0_doublet[index][cnt1][cnt2][cnt3] = (cnt1+cnt2>0) ? new Tensor3RDM(boundary, 0, 1, 1, irr, true,  denBK) : NULL; // NOT j == k == l
               tensor_3rdm_c_J1_doublet[index][cnt1][cnt2][cnt3] = (cnt1+cnt2>0) ? new Tensor3RDM(boundary, 2, 1, 1, irr, true,  denBK) : NULL; // NOT j == k == l
               tensor_3rdm_c_J1_quartet[index][cnt1][cnt2][cnt3] = (cnt1+cnt2>0) ? new Tensor3RDM(boundary, 2, 3, 1, irr, true,  denBK) : NULL; // NOT j == k == l
               tensor_3rdm_d_J0_doublet[index][cnt1][cnt2][cnt3] =                 new Tensor3RDM(boundary, 0, 1, 1, irr, false, denBK);        //     j <= k <= l
               tensor_3rdm_d_J1_doublet[index][cnt1][cnt2][cnt3] =                 new Tensor3RDM(boundary, 2, 1, 1, irr, false, denBK);        //     j <= k <= l
               tensor_3rdm_d_J1_quartet[index][cnt1][cnt2][cnt3] = (cnt2>0)      ? new Tensor3RDM(boundary, 2, 3, 1, irr, false, denBK) : NULL; //     j <= k <  l
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               tensor_3rdm_a_J0_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_a_J1_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_a_J1_quartet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_b_J0_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_b_J1_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_b_J1_quartet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_c_J0_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_c_J1_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_c_J1_quartet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_d_J0_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_d_J1_doublet[index][cnt1][cnt2][cnt3] = NULL;
               tensor_3rdm_d_J1_quartet[index][cnt1][cnt2][cnt3] = NULL;
            }
            #endif
         
         }
      }
   }

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_ALLOC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::delete_3rdm_operators(const int boundary){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   const int index = boundary - 1;
   
   for ( int cnt1 = 0; cnt1 < boundary; cnt1++ ){ // cnt1 = k - j < boundary
      
      for ( int cnt2 = 0; cnt2 < boundary - cnt1; cnt2++ ){ // cnt2 = l - k < boundary - k = boundary - cnt1 - j <= boundary - cnt1
         
         for ( int cnt3 = 0; cnt3 < boundary - cnt1 - cnt2; cnt3++ ){ // cnt3 = boundary - 1 - l < boundary - ( l - k ) - ( k - j ) = boundary - cnt1 - cnt2
            
            const int orb_l = boundary - 1 - cnt3;
            const int orb_k = orb_l - cnt2;
            const int orb_j = orb_k - cnt1;
            
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_3rdm_diagram( L, orb_j, orb_k, orb_l ) == MPIRANK )
            #endif
            {
               if (cnt1+cnt2>0){ delete tensor_3rdm_a_J0_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1>0)     { delete tensor_3rdm_a_J1_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1*cnt2>0){ delete tensor_3rdm_a_J1_quartet[index][cnt1][cnt2][cnt3]; }
               if (cnt2>0)     { delete tensor_3rdm_b_J0_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1*cnt2>0){ delete tensor_3rdm_b_J1_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1*cnt2>0){ delete tensor_3rdm_b_J1_quartet[index][cnt1][cnt2][cnt3]; }
               if (cnt1+cnt2>0){ delete tensor_3rdm_c_J0_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1+cnt2>0){ delete tensor_3rdm_c_J1_doublet[index][cnt1][cnt2][cnt3]; }
               if (cnt1+cnt2>0){ delete tensor_3rdm_c_J1_quartet[index][cnt1][cnt2][cnt3]; }
                                 delete tensor_3rdm_d_J0_doublet[index][cnt1][cnt2][cnt3];
                                 delete tensor_3rdm_d_J1_doublet[index][cnt1][cnt2][cnt3];
               if (cnt2>0)     { delete tensor_3rdm_d_J1_quartet[index][cnt1][cnt2][cnt3]; }
            }
         }
         
         delete [] tensor_3rdm_a_J0_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_a_J1_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_a_J1_quartet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_b_J0_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_b_J1_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_b_J1_quartet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_c_J0_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_c_J1_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_c_J1_quartet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_d_J0_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_d_J1_doublet[ index ][ cnt1 ][ cnt2 ];
         delete [] tensor_3rdm_d_J1_quartet[ index ][ cnt1 ][ cnt2 ];
         
      }
      
      delete [] tensor_3rdm_a_J0_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_a_J1_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_a_J1_quartet[ index ][ cnt1 ];
      delete [] tensor_3rdm_b_J0_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_b_J1_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_b_J1_quartet[ index ][ cnt1 ];
      delete [] tensor_3rdm_c_J0_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_c_J1_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_c_J1_quartet[ index ][ cnt1 ];
      delete [] tensor_3rdm_d_J0_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_d_J1_doublet[ index ][ cnt1 ];
      delete [] tensor_3rdm_d_J1_quartet[ index ][ cnt1 ];
      
   }
   
   delete [] tensor_3rdm_a_J0_doublet[ index ];
   delete [] tensor_3rdm_a_J1_doublet[ index ];
   delete [] tensor_3rdm_a_J1_quartet[ index ];
   delete [] tensor_3rdm_b_J0_doublet[ index ];
   delete [] tensor_3rdm_b_J1_doublet[ index ];
   delete [] tensor_3rdm_b_J1_quartet[ index ];
   delete [] tensor_3rdm_c_J0_doublet[ index ];
   delete [] tensor_3rdm_c_J1_doublet[ index ];
   delete [] tensor_3rdm_c_J1_quartet[ index ];
   delete [] tensor_3rdm_d_J0_doublet[ index ];
   delete [] tensor_3rdm_d_J1_doublet[ index ];
   delete [] tensor_3rdm_d_J1_quartet[ index ];

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_FREE ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::update_correlations_tensors(const int siteindex){

   struct timeval start, end;

   const int dimL = denBK->gMaxDimAtBound(siteindex-1);
   const int dimR = denBK->gMaxDimAtBound(siteindex);
   double * workmemLR = new double[dimL*dimR];
   
   for ( int previousindex = 0; previousindex < siteindex-1; previousindex++ ){
   
      gettimeofday(&start, NULL);
      TensorGYZ * newG = new TensorGYZ(siteindex, 'G', denBK);
      TensorGYZ * newY = new TensorGYZ(siteindex, 'Y', denBK);
      TensorGYZ * newZ = new TensorGYZ(siteindex, 'Z', denBK);
      TensorKM  * newK = new TensorKM( siteindex, 'K', denBK->gIrrep(previousindex), denBK );
      TensorKM  * newM = new TensorKM( siteindex, 'M', denBK->gIrrep(previousindex), denBK );
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_TENS_ALLOC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

      gettimeofday(&start, NULL);
      newG->update(Gtensors[previousindex], MPS[siteindex-1], MPS[siteindex-1], workmemLR);
      newY->update(Ytensors[previousindex], MPS[siteindex-1], MPS[siteindex-1], workmemLR);
      newZ->update(Ztensors[previousindex], MPS[siteindex-1], MPS[siteindex-1], workmemLR);
      newK->update(Ktensors[previousindex], MPS[siteindex-1], MPS[siteindex-1], workmemLR);
      newM->update(Mtensors[previousindex], MPS[siteindex-1], MPS[siteindex-1], workmemLR);
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_TENS_CALC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      gettimeofday(&start, NULL);
      delete Gtensors[previousindex];
      delete Ytensors[previousindex];
      delete Ztensors[previousindex];
      delete Ktensors[previousindex];
      delete Mtensors[previousindex];
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_TENS_FREE ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      Gtensors[previousindex] = newG;
      Ytensors[previousindex] = newY;
      Ztensors[previousindex] = newZ;
      Ktensors[previousindex] = newK;
      Mtensors[previousindex] = newM;
      
   }
   delete [] workmemLR;
   
   gettimeofday(&start, NULL);
   Gtensors[siteindex-1] = new TensorGYZ(siteindex, 'G', denBK);
   Ytensors[siteindex-1] = new TensorGYZ(siteindex, 'Y', denBK);
   Ztensors[siteindex-1] = new TensorGYZ(siteindex, 'Z', denBK);
   Ktensors[siteindex-1] = new TensorKM( siteindex, 'K', denBK->gIrrep(siteindex-1), denBK );
   Mtensors[siteindex-1] = new TensorKM( siteindex, 'M', denBK->gIrrep(siteindex-1), denBK );
   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_ALLOC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   
   gettimeofday(&start, NULL);
   Gtensors[siteindex-1]->construct(MPS[siteindex-1]);
   Ytensors[siteindex-1]->construct(MPS[siteindex-1]);
   Ztensors[siteindex-1]->construct(MPS[siteindex-1]);
   Ktensors[siteindex-1]->construct(MPS[siteindex-1]);
   Mtensors[siteindex-1]->construct(MPS[siteindex-1]);
   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_CALC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}


