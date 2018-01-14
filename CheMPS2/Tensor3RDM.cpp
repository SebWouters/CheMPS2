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
#include <assert.h>
#include <math.h>

#include "Tensor3RDM.h"
#include "Special.h"
#include "Lapack.h"
#include "Wigner.h"

CheMPS2::Tensor3RDM::Tensor3RDM(const int boundary, const int two_j1_in, const int two_j2, const int nelec, const int irrep, const bool prime_last, const SyBookkeeper * book):
TensorOperator(boundary,
               two_j2,
               nelec,
               irrep,
               true, // moving_right
               prime_last,
               true, // jw_phase
               book,
               book){
               
   two_j1 = two_j1_in;

}

CheMPS2::Tensor3RDM::~Tensor3RDM(){ }

int CheMPS2::Tensor3RDM::get_two_j1() const{ return two_j1; }

int CheMPS2::Tensor3RDM::get_two_j2() const{ return get_2j(); }

bool CheMPS2::Tensor3RDM::get_prime_last() const{ return prime_last; }

void CheMPS2::Tensor3RDM::a1(TensorOperator * Sigma, TensorT * denT, double * workmem){

   clear();
   assert( two_j1 == Sigma->get_2j() );
   assert( n_elec == 3 );
   assert( n_irrep == Irreps::directProd( Sigma->get_irrep(), bk_up->gIrrep( index-1 ) ) );
   const int two_j2 = two_j;

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index, nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index, nr_up+3, two_jr_down, ir_down );
      
      { // Contribution 1
         const int il_down = Irreps::directProd( ir_down, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_down = two_jr_down-1; two_jl_down <= two_jr_down+1; two_jl_down+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up,   two_jr_up,   ir_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up+2, two_jl_down, il_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_up - two_jl_down ) <= two_j1 )){
            
               double * Sblock = Sigma->gStorage( nr_up,   two_jr_up,   ir_up,   nr_up+2, two_jl_down, il_down );
               double * Tup    =  denT->gStorage( nr_up,   two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  =  denT->gStorage( nr_up+2, two_jl_down, il_down, nr_up+3, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jl_down + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_up, two_jr_down, two_jl_down )
                            * Special::phase( two_jr_up + two_jr_down + two_j1 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Sblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }        
         }
      }
      { // Contribution 2
         const int il_up = Irreps::directProd( ir_up, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_up = two_jr_up-1; two_jl_up <= two_jr_up+1; two_jl_up+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up-1, two_jl_up,   il_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up+1, two_jr_down, ir_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_down - two_jl_up ) <= two_j1 )){
            
               double * Sblock = Sigma->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up+1, two_jr_down, ir_down );
               double * Tup    =  denT->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  =  denT->gStorage( nr_up+1, two_jr_down, ir_down, nr_up+3, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jr_up + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_down, two_jr_up, two_jl_up )
                            * Special::phase( two_jl_up + two_jr_down + two_j2 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Sblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }
         }
      }
   }
}

void CheMPS2::Tensor3RDM::b1(TensorOperator * Sigma, TensorT * denT, double * workmem){

   clear();
   assert( two_j1 == Sigma->get_2j() );
   assert( n_elec == 1 );
   assert( n_irrep == Irreps::directProd( Sigma->get_irrep(), bk_up->gIrrep( index-1 ) ) );
   const int two_j2 = two_j;

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index, nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index, nr_up+1, two_jr_down, ir_down );
      
      { // Contribution 1
         const int il_up = Irreps::directProd( ir_up, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_up = two_jr_up-1; two_jl_up <= two_jr_up+1; two_jl_up+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up-1, two_jl_up,   il_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up+1, two_jr_down, ir_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_down - two_jl_up ) <= two_j1 )){
            
               double * Sblock = Sigma->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up+1, two_jr_down, ir_down );
               double * Tup    =  denT->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  =  denT->gStorage( nr_up+1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jr_up + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_down, two_jr_up, two_jl_up )
                            * Special::phase( two_jl_up + two_jr_down + two_j2 + 3 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Sblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }
         }
      }
      { // Contribution 2
         const int il_down = Irreps::directProd( ir_down, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_down = two_jr_down-1; two_jl_down <= two_jr_down+1; two_jl_down+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up-2, two_jr_up,   ir_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up,   two_jl_down, il_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_up - two_jl_down ) <= two_j1 )){
            
               double * Sblock = Sigma->gStorage( nr_up-2, two_jr_up,   ir_up,   nr_up,   two_jl_down, il_down );
               double * Tup    =  denT->gStorage( nr_up-2, two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  =  denT->gStorage( nr_up,   two_jl_down, il_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jl_down + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_up, two_jr_down, two_jl_down )
                            * Special::phase( two_jr_up + two_jr_down + two_j1 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Sblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }        
         }
      }
   }
}

void CheMPS2::Tensor3RDM::c1(TensorOperator * denF, TensorT * denT, double * workmem){

   clear();
   assert( two_j1 == denF->get_2j() );
   assert( n_elec == 1 );
   assert( n_irrep == Irreps::directProd( denF->get_irrep(), bk_up->gIrrep( index-1 ) ) );
   const int two_j2 = two_j;

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index, nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index, nr_up+1, two_jr_down, ir_down );
      
      { // Contribution 1
         const int il_down = Irreps::directProd( ir_down, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_down = two_jr_down-1; two_jl_down <= two_jr_down+1; two_jl_down+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up, two_jr_up,   ir_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up, two_jl_down, il_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_up - two_jl_down ) <= two_j1 )){
            
               double * Fblock = denF->gStorage( nr_up, two_jr_up,   ir_up,   nr_up,   two_jl_down, il_down );
               double * Tup    = denT->gStorage( nr_up, two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  = denT->gStorage( nr_up, two_jl_down, il_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jl_down + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_up, two_jr_down, two_jl_down )
                            * Special::phase( two_jr_up + two_jr_down + two_j1 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Fblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }        
         }
      }
      { // Contribution 2
         const int il_up = Irreps::directProd( ir_up, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_up = two_jr_up-1; two_jl_up <= two_jr_up+1; two_jl_up+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up-1, two_jl_up,   il_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up-1, two_jr_down, ir_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_down - two_jl_up ) <= two_j1 )){
            
               double * Fblock = denF->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up-1, two_jr_down, ir_down );
               double * Tup    = denT->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jr_up + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_down, two_jr_up, two_jl_up )
                            * Special::phase( two_jl_up + two_jr_down + two_j2 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Fblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }
         }
      }
   }
}

void CheMPS2::Tensor3RDM::d1(TensorOperator * denF, TensorT * denT, double * workmem){

   clear();
   assert( two_j1 == denF->get_2j() );
   assert( n_elec == 1 );
   assert( n_irrep == Irreps::directProd( denF->get_irrep(), bk_up->gIrrep( index-1 ) ) );
   const int two_j2 = two_j;

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index, nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index, nr_up+1, two_jr_down, ir_down );
      
      { // Contribution 1
         const int il_down = Irreps::directProd( ir_down, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_down = two_jr_down-1; two_jl_down <= two_jr_down+1; two_jl_down+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up, two_jr_up,   ir_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up, two_jl_down, il_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_up - two_jl_down ) <= two_j1 )){
            
               double * Fblock = denF->gStorage( nr_up, two_jl_down, il_down, nr_up,   two_jr_up,   ir_up   );
               double * Tup    = denT->gStorage( nr_up, two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  = denT->gStorage( nr_up, two_jl_down, il_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jr_down + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_up, two_jr_down, two_jl_down )
                            * Special::phase( two_jr_up + two_jl_down + two_j2 + 3 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &trans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Fblock, &dimLdown, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }        
         }
      }
      { // Contribution 2
         const int il_up = Irreps::directProd( ir_up, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_up = two_jr_up-1; two_jl_up <= two_jr_up+1; two_jl_up+=2 ){
         
            int dimLup   = bk_up->gCurrentDim( index-1, nr_up-1, two_jl_up,   il_up   );
            int dimLdown = bk_up->gCurrentDim( index-1, nr_up-1, two_jr_down, ir_down );
            
            if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jr_down - two_jl_up ) <= two_j1 )){
            
               double * Fblock = denF->gStorage( nr_up-1, two_jr_down, ir_down, nr_up-1, two_jl_up,   il_up   );
               double * Tup    = denT->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jr_up,   ir_up   );
               double * Tdown  = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
               
               double alpha = sqrt( 1.0 * ( two_j2 + 1 ) * ( two_jl_up + 1 ) )
                            * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_down, two_jr_up, two_jl_up )
                            * Special::phase( two_jr_up + two_jr_down + two_j1 + 1 );
               double beta  = 0.0; //set
               char trans   = 'T';
               char notrans = 'N';
               dgemm_( &trans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Fblock, &dimLdown, Tdown, &dimLdown, &beta, workmem, &dimLup );
               alpha = 1.0;
               beta  = 1.0; //add
               dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
            
            }
         }
      }
   }
}

void CheMPS2::Tensor3RDM::extra1(TensorT * denT){

   clear();
   assert( n_elec == 1 );
   assert( n_irrep == bk_up->gIrrep( index-1 ) );
   const int two_j2 = two_j;
   assert( two_j2 == 1 );

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index,   nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index,   nr_up+1, two_jr_down, ir_down );
      int dimL     = bk_up->gCurrentDim( index-1, nr_up-1, two_jr_down, ir_down );
      
      if ( dimL > 0 ){
      
         double * Tup   = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up,   two_jr_up,   ir_up   );
         double * Tdown = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
         
         double alpha = sqrt( 0.5 * ( two_j1 + 1 ) ) * Special::phase( two_j1 );
         double beta  = 0.0; //set --> only contribution to this symmetry sector
         char trans   = 'T';
         char notrans = 'N';
         dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimL, &alpha, Tup, &dimL, Tdown, &dimL, &beta, storage + kappa2index[ikappa], &dimRup );
      
      }
   }
}

void CheMPS2::Tensor3RDM::extra2(TensorL * denL, TensorT * denT, double * workmem){

   clear();
   assert( n_elec == 3 );
   assert( n_irrep == denL->get_irrep() );
   const int two_j2 = two_j;
   assert( two_j2 == 1 );

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index,   nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index,   nr_up+3, two_jr_down, ir_down );
      int dimLup   = bk_up->gCurrentDim( index-1, nr_up,   two_jr_up,   ir_up   );
      int dimLdown = bk_up->gCurrentDim( index-1, nr_up+1, two_jr_down, ir_down );
      
      if (( dimLup > 0 ) && ( dimLdown > 0 )){
      
         double * Tup    = denT->gStorage( nr_up,   two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
         double * Tdown  = denT->gStorage( nr_up+1, two_jr_down, ir_down, nr_up+3, two_jr_down, ir_down );
         double * Lblock = denL->gStorage( nr_up,   two_jr_up,   ir_up,   nr_up+1, two_jr_down, ir_down );
         
         double alpha = sqrt( 0.5 * ( two_j1 + 1 ) ) * Special::phase( two_j1 + 2 );
         double beta  = 0.0; //set
         char trans   = 'T';
         char notrans = 'N';
         dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Lblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
         alpha = 1.0;
         dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
      
      }
   }
}

void CheMPS2::Tensor3RDM::extra3(TensorL * denL, TensorT * denT, double * workmem){

   clear();
   assert( n_elec == 1 );
   assert( n_irrep == denL->get_irrep() );
   const int two_j2 = two_j;
   assert( two_j2 == 1 );

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index,   nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index,   nr_up+1, two_jr_down, ir_down );
      int dimLup   = bk_up->gCurrentDim( index-1, nr_up,   two_jr_up,   ir_up   );
      int dimLdown = bk_up->gCurrentDim( index-1, nr_up-1, two_jr_down, ir_down );
      
      if (( dimLup > 0 ) && ( dimLdown > 0 )){
      
         double * Tup    = denT->gStorage( nr_up,   two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
         double * Tdown  = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
         double * Lblock = denL->gStorage( nr_up-1, two_jr_down, ir_down, nr_up,   two_jr_up,   ir_up   );
         
         double alpha = sqrt( 0.5 * ( two_j1 + 1 ) ) * Special::phase( two_j1 );
         double beta  = 0.0; //set
         char trans   = 'T';
         char notrans = 'N';
         dgemm_( &trans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Lblock, &dimLdown, Tdown, &dimLdown, &beta, workmem, &dimLup );
         alpha = 1.0;
         dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
      
      }
   }
}

void CheMPS2::Tensor3RDM::extra4(TensorL * denL, TensorT * denT, double * workmem){

   clear();
   assert( n_elec == 1 );
   assert( n_irrep == denL->get_irrep() );
   const int two_j2 = two_j;

   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
       
      const int two_jr_up   = sector_spin_up[ ikappa ];
      const int nr_up       = sector_nelec_up[ ikappa ];
      const int ir_up       = sector_irrep_up[ ikappa ];
      const int two_jr_down = sector_spin_down[ ikappa ];
      const int ir_down     = Irreps::directProd( ir_up, n_irrep );
      
      int dimRup   = bk_up->gCurrentDim( index,   nr_up,   two_jr_up,   ir_up   );
      int dimRdown = bk_up->gCurrentDim( index,   nr_up+1, two_jr_down, ir_down );
      
      if ( two_j2 == 1 ){ // Contribution 2
      
         int dimLup   = bk_up->gCurrentDim( index-1, nr_up-2, two_jr_up,   ir_up   );
         int dimLdown = bk_up->gCurrentDim( index-1, nr_up-1, two_jr_down, ir_down );
         
         if (( dimLup > 0 ) && ( dimLdown > 0 )){
         
            double * Tup    = denT->gStorage( nr_up-2, two_jr_up,   ir_up,   nr_up,   two_jr_up,   ir_up   );
            double * Tdown  = denT->gStorage( nr_up-1, two_jr_down, ir_down, nr_up+1, two_jr_down, ir_down );
            double * Lblock = denL->gStorage( nr_up-2, two_jr_up,   ir_up,   nr_up-1, two_jr_down, ir_down );
            
            double alpha = sqrt( 0.5 * ( two_j1 + 1 ) ) * Special::phase( two_j1 + 2 );
            double beta  = 0.0; //set
            char trans   = 'T';
            char notrans = 'N';
            dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Lblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
            alpha = 1.0;
            beta  = 1.0; // add
            dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
         
         }
      }
      { // Contribution 1
         const int il_up   = Irreps::directProd( ir_up,   bk_up->gIrrep( index-1 ) );
         const int il_down = Irreps::directProd( ir_down, bk_up->gIrrep( index-1 ) );
         for ( int two_jl_up = two_jr_up-1; two_jl_up <= two_jr_up+1; two_jl_up+=2 ){
            for ( int two_jl_down = two_jr_down-1; two_jl_down <= two_jr_down+1; two_jl_down+=2 ){
      
               int dimLup   = bk_up->gCurrentDim( index-1, nr_up-1, two_jl_up,   il_up   );
               int dimLdown = bk_up->gCurrentDim( index-1, nr_up,   two_jl_down, il_down );
               
               if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( abs( two_jl_up - two_jl_down ) <= 1 )){
               
                  double * Tup    = denT->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jr_up,   ir_up   );
                  double * Tdown  = denT->gStorage( nr_up,   two_jl_down, il_down, nr_up+1, two_jr_down, ir_down );
                  double * Lblock = denL->gStorage( nr_up-1, two_jl_up,   il_up,   nr_up,   two_jl_down, il_down );
                  
                  double alpha = sqrt( 1.0 * ( two_j1 + 1 ) * ( two_j2 + 1 ) * ( two_jr_up + 1 ) * ( two_jl_down + 1 ) )
                               * Special::phase( two_jr_up + two_jr_down - two_jl_up - two_jl_down )
                               * Wigner::wigner6j( 1, 1, two_j1, two_jl_down, two_jr_up, two_jl_up )
                               * Wigner::wigner6j( 1, two_j1, two_j2, two_jr_up, two_jr_down, two_jl_down );
                  double beta  = 0.0; //set
                  char trans   = 'T';
                  char notrans = 'N';
                  dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &alpha, Lblock, &dimLup, Tdown, &dimLdown, &beta, workmem, &dimLup );
                  alpha = 1.0;
                  beta  = 1.0; // add
                  dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, workmem, &dimLup, &beta, storage + kappa2index[ikappa], &dimRup );
               
               }
            }
         }
      }
   }
}

double CheMPS2::Tensor3RDM::contract( Tensor3RDM * buddy ) const{

   if ( buddy == NULL ){ return 0.0; }

   assert( get_two_j2() == buddy->get_two_j2() );
   assert( n_elec       == buddy->get_nelec()  );
   assert( n_irrep      == buddy->get_irrep()  );
   
   double value = 0.0;

   if ( buddy->get_prime_last() ){ // Tensor abc
   
      int length = kappa2index[ nKappa ];
      int inc    = 1;
      value = ddot_( &length, storage, &inc, buddy->gStorage(), &inc );
      return value;
         
   } else { // Tensor d
   
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
      
         int offset = kappa2index[ ikappa ];
         int length = kappa2index[ ikappa + 1 ] - offset;
         int inc    = 1;
         double prefactor = sqrt( ( sector_spin_up[ ikappa ] + 1.0 ) / ( sector_spin_down[ ikappa ] + 1 ) )
                          * Special::phase( sector_spin_up[ ikappa ] + 1 - sector_spin_down[ ikappa ] );
         value += prefactor * ddot_( &length, storage + offset, &inc, buddy->gStorage() + offset, &inc );
      
      }

      return value;

   }

   return value;

}


