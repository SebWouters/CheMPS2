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

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_sf_coupling.h>

#include "TensorOperator.h"
#include "Lapack.h"
#include "Heff.h"

CheMPS2::TensorOperator::TensorOperator(const int index_in, const int two_j_in, const int n_elec_in, const int n_irrep_in, const bool moving_right_in, const bool prime_last_in, const bool jw_phase_in, const SyBookkeeper * denBK_in) : Tensor(){

   // Copy the variables
   index        = index_in;
   two_j        = two_j_in;
   n_elec       = n_elec_in;
   n_irrep      = n_irrep_in;
   moving_right = moving_right_in;
   prime_last   = prime_last_in;
   jw_phase     = jw_phase_in;
   denBK        = denBK_in;
   
   assert( two_j >= 0 );
   assert( n_irrep < denBK->getNumberOfIrreps() );
   
   nKappa = 0;
   for (int n_up = denBK->gNmin(index); n_up <= denBK->gNmax(index); n_up++){
      for (int two_s_up = denBK->gTwoSmin(index, n_up); two_s_up <= denBK->gTwoSmax(index, n_up); two_s_up += 2){
         for (int irrep_up = 0; irrep_up < denBK->getNumberOfIrreps(); irrep_up++){
            int dim_up = denBK->gCurrentDim(index, n_up, two_s_up, irrep_up);
            if (dim_up > 0){
               const int irrep_down = Irreps::directProd(n_irrep, irrep_up);
               const int n_down = n_up + n_elec;
               for (int two_s_down = two_s_up - two_j; two_s_down <= two_s_up + two_j; two_s_down += 2){
                  if (two_s_down >= 0){
                     int dim_down = denBK->gCurrentDim(index, n_down, two_s_down, irrep_down);
                     if (dim_down > 0){
                        nKappa++;
                     }
                  }
               }
            }
         }
      }
   }
   
   sectorN1 = new int[nKappa];
   sectorTwoS1 = new int[nKappa];
   sectorI1 = new int[nKappa];
   sector_2S_down = (( two_j == 0 ) ? sectorTwoS1 : new int[nKappa]);
   kappa2index = new int[nKappa+1];
   kappa2index[0] = 0;
   
   nKappa = 0;
   for (int n_up = denBK->gNmin(index); n_up <= denBK->gNmax(index); n_up++){
      for (int two_s_up = denBK->gTwoSmin(index, n_up); two_s_up <= denBK->gTwoSmax(index, n_up); two_s_up += 2){
         for (int irrep_up = 0; irrep_up < denBK->getNumberOfIrreps(); irrep_up++){
            int dim_up = denBK->gCurrentDim(index, n_up, two_s_up, irrep_up);
            if (dim_up > 0){
               const int irrep_down = Irreps::directProd(n_irrep, irrep_up);
               const int n_down = n_up + n_elec;
               for (int two_s_down = two_s_up - two_j; two_s_down <= two_s_up + two_j; two_s_down += 2){
                  if (two_s_down >= 0){
                     int dim_down = denBK->gCurrentDim(index, n_down, two_s_down, irrep_down);
                     if (dim_down > 0){
                        sectorN1[nKappa] = n_up;
                        sectorTwoS1[nKappa] = two_s_up;
                        sectorI1[nKappa] = irrep_up;
                        sector_2S_down[nKappa] = two_s_down;
                        nKappa++;
                        kappa2index[nKappa] = kappa2index[nKappa-1] + dim_up*dim_down;
                     }
                  }
               }
            }
         }
      }
   }
   
   storage = new double[kappa2index[nKappa]];

}

CheMPS2::TensorOperator::~TensorOperator(){

   delete [] sectorN1;
   delete [] sectorTwoS1;
   delete [] sectorI1;
   delete [] kappa2index;
   delete [] storage;
   if ( two_j != 0 ){ delete [] sector_2S_down; }

}

int CheMPS2::TensorOperator::gNKappa() const { return nKappa; }

double * CheMPS2::TensorOperator::gStorage() { return storage; }
      
int CheMPS2::TensorOperator::gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const{

   if ( (Irreps::directProd(I1,n_irrep)) != I2 ){ return -1; }
   if ( N2 != N1 + n_elec ){ return -1; }
   if ( abs( TwoS1 - TwoS2 ) > two_j ){ return -1; }

   for (int cnt = 0; cnt < nKappa; cnt++){
      if ((sectorN1[cnt]==N1)&&(sectorTwoS1[cnt]==TwoS1)&&(sectorI1[cnt]==I1)&&(sector_2S_down[cnt]==TwoS2)){ return cnt; }
   }
   
   return -1;

}
      
int CheMPS2::TensorOperator::gKappa2index(const int kappa) const{ return kappa2index[kappa]; }
      
double * CheMPS2::TensorOperator::gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2){

   int kappa = gKappa(N1,TwoS1,I1,N2,TwoS2,I2);
   if (kappa == -1) return NULL;
   return storage + kappa2index[kappa];

}

int CheMPS2::TensorOperator::gIndex() const { return index; }

int CheMPS2::TensorOperator::get_2j() const{ return two_j; }

int CheMPS2::TensorOperator::get_nelec() const{ return n_elec; }

int CheMPS2::TensorOperator::get_irrep() const { return n_irrep; }

void CheMPS2::TensorOperator::clear(){

   for (int cnt = 0; cnt < kappa2index[nKappa]; cnt++){ storage[cnt] = 0.0; }

}

void CheMPS2::TensorOperator::update(TensorOperator * previous, TensorT * mps_tensor, double * workmem){

   clear();

   if ( moving_right ){
      for (int ikappa=0; ikappa<nKappa; ikappa++){
         update_moving_right(ikappa, previous, mps_tensor, workmem);
      }
   } else {
      for (int ikappa=0; ikappa<nKappa; ikappa++){
         update_moving_left(ikappa, previous, mps_tensor, workmem);
      }
   }

}

void CheMPS2::TensorOperator::update_moving_right(const int ikappa, TensorOperator * previous, TensorT * mps_tensor, double * workmem){

   const int n_right_up       = sectorN1[ ikappa ];
   const int n_right_down     = n_right_up + n_elec;
   const int two_s_right_up   = sectorTwoS1[ ikappa ];
   const int two_s_right_down = sector_2S_down[ ikappa ];
   const int irrep_right_up   = sectorI1[ ikappa ];
   const int irrep_right_down = Irreps::directProd( irrep_right_up, n_irrep );

   int dim_right_up   = denBK->gCurrentDim( index, n_right_up,   two_s_right_up,   irrep_right_up   );
   int dim_right_down = denBK->gCurrentDim( index, n_right_down, two_s_right_down, irrep_right_down );
   
   for (int geval = 0; geval < 6; geval++){
      int n_left_up, n_left_down, two_s_left_up, two_s_left_down, irrep_left_up, irrep_left_down;
      switch ( geval ){
         case 0: // MPS tensor sector (I,J,N) = (0,0,0)
            two_s_left_up   = two_s_right_up;
            two_s_left_down = two_s_right_down;
            n_left_up       = n_right_up;
            n_left_down     = n_right_down;
            irrep_left_up   = irrep_right_up;
            irrep_left_down = irrep_right_down;
            break;
         case 1: // MPS tensor sector (I,J,N) = (0,0,2)
            two_s_left_up   = two_s_right_up;
            two_s_left_down = two_s_right_down;
            n_left_up       = n_right_up - 2;
            n_left_down     = n_right_down - 2;
            irrep_left_up   = irrep_right_up;
            irrep_left_down = irrep_right_down;
            break;
         case 2: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_left_up   = two_s_right_up - 1;
            two_s_left_down = two_s_right_down - 1;
            n_left_up       = n_right_up - 1;
            n_left_down     = n_right_down - 1;
            irrep_left_up   = Irreps::directProd( irrep_right_up,   denBK->gIrrep(index-1) );
            irrep_left_down = Irreps::directProd( irrep_right_down, denBK->gIrrep(index-1) );
            break;
         case 3: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_left_up   = two_s_right_up - 1;
            two_s_left_down = two_s_right_down + 1;
            n_left_up       = n_right_up - 1;
            n_left_down     = n_right_down - 1;
            irrep_left_up   = Irreps::directProd( irrep_right_up,   denBK->gIrrep(index-1) );
            irrep_left_down = Irreps::directProd( irrep_right_down, denBK->gIrrep(index-1) );
            break;
         case 4: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_left_up   = two_s_right_up + 1;
            two_s_left_down = two_s_right_down - 1;
            n_left_up       = n_right_up - 1;
            n_left_down     = n_right_down - 1;
            irrep_left_up   = Irreps::directProd( irrep_right_up,   denBK->gIrrep(index-1) );
            irrep_left_down = Irreps::directProd( irrep_right_down, denBK->gIrrep(index-1) );
            break;
         case 5: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_left_up   = two_s_right_up + 1;
            two_s_left_down = two_s_right_down + 1;
            n_left_up       = n_right_up - 1;
            n_left_down     = n_right_down - 1;
            irrep_left_up   = Irreps::directProd( irrep_right_up,   denBK->gIrrep(index-1) );
            irrep_left_down = Irreps::directProd( irrep_right_down, denBK->gIrrep(index-1) );
            break;
      }
      
      if ( abs( two_s_left_up - two_s_left_down ) <= two_j ){
      
         int dim_left_up   = denBK->gCurrentDim( index-1, n_left_up,   two_s_left_up,   irrep_left_up   );
         int dim_left_down = denBK->gCurrentDim( index-1, n_left_down, two_s_left_down, irrep_left_down );
         
         if (( dim_left_up > 0 ) && ( dim_left_down > 0 )){

            double * mps_block_up   = mps_tensor->gStorage( n_left_up,   two_s_left_up,   irrep_left_up,   n_right_up,   two_s_right_up,   irrep_right_up   );
            double * mps_block_down = mps_tensor->gStorage( n_left_down, two_s_left_down, irrep_left_down, n_right_down, two_s_right_down, irrep_right_down );
            double * left_block     =   previous->gStorage( n_left_up,   two_s_left_up,   irrep_left_up,   n_left_down,  two_s_left_down,  irrep_left_down  );

            // Prefactor
            double alpha = 1.0;
            if ( geval >= 2 ){
               if ( two_j == 0 ){
                  alpha = ( ( jw_phase ) ? -1.0 : 1.0 );
               } else {
                  if ( prime_last ){
                     alpha = CheMPS2::Heff::phase( two_s_right_up + two_s_left_down + two_j + ( ( jw_phase ) ? 3 : 1 ) )
                           * sqrt( ( two_s_left_down + 1.0 ) * ( two_s_right_up + 1.0 ) )
                           * gsl_sf_coupling_6j( two_s_left_up, two_s_left_down, two_j, two_s_right_down, two_s_right_up, 1 );
                  } else {
                     alpha = CheMPS2::Heff::phase( two_s_right_down + two_s_left_up + two_j + ( ( jw_phase ) ? 3 : 1 ) )
                           * sqrt( ( two_s_left_up + 1.0 ) * ( two_s_right_down + 1.0 ) )
                           * gsl_sf_coupling_6j( two_s_left_down, two_s_left_up, two_j, two_s_right_up, two_s_right_down, 1 );
                  }
               }
            }
            
            // prefactor * mps_block_up^T * left_block --> mem
            char trans = 'T';
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&trans, &notr, &dim_right_up, &dim_left_down, &dim_left_up,
                   &alpha, mps_block_up, &dim_left_up, left_block, &dim_left_up,
                   &beta, workmem, &dim_right_up);

            // mem * mps_block_down --> storage
            alpha = 1.0;
            beta = 1.0; //add
            dgemm_(&notr, &notr, &dim_right_up, &dim_right_down, &dim_left_down,
                   &alpha, workmem, &dim_right_up, mps_block_down, &dim_left_down,
                   &beta, storage + kappa2index[ikappa], &dim_right_up);
         }
      }
   }

}

void CheMPS2::TensorOperator::update_moving_left(const int ikappa, TensorOperator * previous, TensorT * mps_tensor, double * workmem){

   const int n_left_up       = sectorN1[ ikappa ];
   const int n_left_down     = n_left_up + n_elec;
   const int two_s_left_up   = sectorTwoS1[ ikappa ];
   const int two_s_left_down = sector_2S_down[ ikappa ];
   const int irrep_left_up   = sectorI1[ ikappa ];
   const int irrep_left_down = Irreps::directProd( irrep_left_up, n_irrep );

   int dim_left_up   = denBK->gCurrentDim( index, n_left_up,   two_s_left_up,   irrep_left_up   );
   int dim_left_down = denBK->gCurrentDim( index, n_left_down, two_s_left_down, irrep_left_down );
   
   for (int geval = 0; geval < 6; geval++){
      int n_right_up, n_right_down, two_s_right_up, two_s_right_down, irrep_right_up, irrep_right_down;
      switch ( geval ){
         case 0: // MPS tensor sector (I,J,N) = (0,0,0)
            two_s_right_up   = two_s_left_up;
            two_s_right_down = two_s_left_down;
            n_right_up       = n_left_up;
            n_right_down     = n_left_down;
            irrep_right_up   = irrep_left_up;
            irrep_right_down = irrep_left_down;
            break;
         case 1: // MPS tensor sector (I,J,N) = (0,0,2)
            two_s_right_up   = two_s_left_up;
            two_s_right_down = two_s_left_down;
            n_right_up       = n_left_up + 2;
            n_right_down     = n_left_down + 2;
            irrep_right_up   = irrep_left_up;
            irrep_right_down = irrep_left_down;
            break;
         case 2: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_right_up   = two_s_left_up - 1;
            two_s_right_down = two_s_left_down - 1;
            n_right_up       = n_left_up + 1;
            n_right_down     = n_left_down + 1;
            irrep_right_up   = Irreps::directProd( irrep_left_up,   denBK->gIrrep(index) );
            irrep_right_down = Irreps::directProd( irrep_left_down, denBK->gIrrep(index) );
            break;
         case 3: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_right_up   = two_s_left_up - 1;
            two_s_right_down = two_s_left_down + 1;
            n_right_up       = n_left_up + 1;
            n_right_down     = n_left_down + 1;
            irrep_right_up   = Irreps::directProd( irrep_left_up,   denBK->gIrrep(index) );
            irrep_right_down = Irreps::directProd( irrep_left_down, denBK->gIrrep(index) );
            break;
         case 4: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_right_up   = two_s_left_up + 1;
            two_s_right_down = two_s_left_down - 1;
            n_right_up       = n_left_up + 1;
            n_right_down     = n_left_down + 1;
            irrep_right_up   = Irreps::directProd( irrep_left_up,   denBK->gIrrep(index) );
            irrep_right_down = Irreps::directProd( irrep_left_down, denBK->gIrrep(index) );
            break;
         case 5: // MPS tensor sector (I,J,N) = (Ilocal,1/2,1)
            two_s_right_up   = two_s_left_up + 1;
            two_s_right_down = two_s_left_down + 1;
            n_right_up       = n_left_up + 1;
            n_right_down     = n_left_down + 1;
            irrep_right_up   = Irreps::directProd( irrep_left_up,   denBK->gIrrep(index) );
            irrep_right_down = Irreps::directProd( irrep_left_down, denBK->gIrrep(index) );
            break;
      }
      
      if ( abs( two_s_right_up - two_s_right_down ) <= two_j ){
      
         int dim_right_up   = denBK->gCurrentDim( index+1, n_right_up,   two_s_right_up,   irrep_right_up   );
         int dim_right_down = denBK->gCurrentDim( index+1, n_right_down, two_s_right_down, irrep_right_down );
         
         if (( dim_right_up > 0 ) && ( dim_right_down > 0 )){

            double * mps_block_up   = mps_tensor->gStorage( n_left_up,   two_s_left_up,   irrep_left_up,   n_right_up,   two_s_right_up,   irrep_right_up   );
            double * mps_block_down = mps_tensor->gStorage( n_left_down, two_s_left_down, irrep_left_down, n_right_down, two_s_right_down, irrep_right_down );
            double * right_block    =   previous->gStorage( n_right_up,  two_s_right_up,  irrep_right_up,  n_right_down, two_s_right_down, irrep_right_down );

            // Prefactor
            double alpha = 1.0;
            if ( geval >= 2 ){
               if ( two_j == 0 ){
                  alpha = ( ( jw_phase ) ? -1.0 : 1.0 ) * ( two_s_right_up + 1.0 ) / ( two_s_left_up + 1.0 );
               } else {
                  if ( prime_last ){
                     alpha = CheMPS2::Heff::phase( two_s_right_up + two_s_left_down + two_j + ( ( jw_phase ) ? 3 : 1 ) )
                           * ( two_s_right_down + 1 ) * sqrt( ( two_s_right_up + 1.0 ) / ( two_s_left_down + 1.0 ) )
                           * gsl_sf_coupling_6j( two_s_right_up, two_s_right_down, two_j, two_s_left_down, two_s_left_up, 1 );
                  } else {
                     alpha = CheMPS2::Heff::phase( two_s_right_down + two_s_left_up + two_j + ( ( jw_phase ) ? 3 : 1 ) )
                           * ( two_s_right_up + 1 ) * sqrt( ( two_s_right_down + 1.0 ) / ( two_s_left_up + 1.0 ) )
                           * gsl_sf_coupling_6j( two_s_right_down, two_s_right_up, two_j, two_s_left_up, two_s_left_down, 1 );
                  }
               }
            }
            
            // prefactor * mps_block_up * right_block --> mem
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&notr, &notr, &dim_left_up, &dim_right_down, &dim_right_up,
                   &alpha, mps_block_up, &dim_left_up, right_block, &dim_right_up,
                   &beta, workmem, &dim_left_up);

            // mem * mps_block_down^T --> storage
            char trans = 'T';
            alpha = 1.0;
            beta = 1.0; //add
            dgemm_(&notr, &trans, &dim_left_up, &dim_left_down, &dim_right_down,
                   &alpha, workmem, &dim_left_up, mps_block_down, &dim_left_down,
                   &beta, storage + kappa2index[ikappa], &dim_left_up);
         }
      }
   }

}

void CheMPS2::TensorOperator::daxpy(double alpha, TensorOperator * to_add){

   assert( kappa2index[nKappa] == to_add->gKappa2index(to_add->gNKappa()) );
   int inc = 1;
   daxpy_(kappa2index+nKappa, &alpha, to_add->gStorage(), &inc, storage, &inc);

}

void CheMPS2::TensorOperator::daxpy_transpose_tensorCD(const double alpha, TensorOperator * to_add){

   assert( kappa2index[nKappa] == to_add->gKappa2index(to_add->gNKappa()) );
   assert( n_elec == 0 );
   assert( ( two_j == 0 ) || ( two_j == 2 ) );

   for (int ikappa=0; ikappa<nKappa; ikappa++){

      const int irrep_up   = sectorI1[ikappa];
      const int irrep_down = Irreps::directProd( irrep_up, n_irrep );
      const int two_s_up   = sectorTwoS1[ikappa];
      const int two_s_down = sector_2S_down[ikappa];
      const int n_updown   = sectorN1[ikappa];

      const int dim_up   = denBK->gCurrentDim( index, n_updown, two_s_up,   irrep_up   );
      const int dim_down = denBK->gCurrentDim( index, n_updown, two_s_down, irrep_down );

      double prefactor = alpha;
      /*
         This phase factor comes historically from the TensorD and is not valid in general,
         as it is tightly coupled to the specific change from (for moving_right == true ):
           < 1/2 m1 1/2 -m2 | 1 (m1-m2) > * (-1)^{1/2-m2} * < j_L' j_L^z' 1 (m1-m2) | j_L  j_L^z  >
         = < 1/2 m2 1/2 -m1 | 1 (m2-m1) > * (-1)^{1/2-m1} * < j_L  j_L^z  1 (m1-m2) | j_L' j_L^z' > * prefactor
      */
      if ( two_s_up != two_s_down ){
         prefactor *= CheMPS2::Heff::phase( two_s_up - two_s_down )
                    * sqrt(( moving_right ) ? (( two_s_up + 1.0 ) / ( two_s_down + 1.0 )) : (( two_s_down + 1.0 ) / ( two_s_up + 1.0 )));
      }

      double * block = to_add->gStorage( n_updown, two_s_down, irrep_down, n_updown, two_s_up, irrep_up );
      for (int irow = 0; irow < dim_up; irow++){
         for (int icol = 0; icol < dim_down; icol++){
            storage[ kappa2index[ikappa] + irow + dim_up * icol ] += prefactor * block[ icol + dim_down * irow ];
         }
      }

   }

}

double CheMPS2::TensorOperator::inproduct( TensorOperator * buddy, const char trans ) const{

   if ( buddy == NULL ){ return 0.0; }

   const int orb_i = gIndex();
   assert( get_2j() == buddy->get_2j()    );
   assert( n_elec   == buddy->get_nelec() );
   assert( n_irrep  == buddy->get_irrep() );

   double value = 0.0;

   if ( trans == 'N' ){

      int length = kappa2index[ nKappa ];
      int inc    = 1;
      value = ddot_( &length, storage, &inc, buddy->gStorage(), &inc );
      return value;

   } else {

      assert( n_elec == 0 );
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){

         const int n_updown   = sectorN1[ ikappa ];
         const int two_j_up   = sectorTwoS1[ ikappa ];
         const int two_j_down = sector_2S_down[ ikappa ];
         const int irrep_up   = sectorI1[ ikappa ];
         const int irrep_down = Irreps::directProd( sectorI1[ ikappa ], n_irrep );

         double * my_block    = storage + kappa2index[ ikappa ];
         double * buddy_block = buddy->gStorage( n_updown, two_j_down, irrep_down, n_updown, two_j_up, irrep_up );
         const int dim_up     = denBK->gCurrentDim( index, n_updown, two_j_up,   irrep_up   );
         const int dim_down   = denBK->gCurrentDim( index, n_updown, two_j_down, irrep_down );
         
         double temp = 0.0;
         for ( int row = 0; row < dim_up; row++ ){
            for ( int col = 0; col < dim_down; col++ ){
               temp += my_block[ row + dim_up * col ] * buddy_block[ col + dim_down * row ];
            }
         }
         
         const double prefactor = (( get_2j() == 0 ) ? 1.0 : ( sqrt( ( two_j_up + 1.0 ) / ( two_j_down + 1.0 ) ) * Heff::phase( two_j_up - two_j_down ) ));
         value += prefactor * temp;
      
      }

      return value;

   }

   return value;

}


