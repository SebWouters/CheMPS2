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
#include <math.h>
#include <assert.h>

#include "TensorL.h"
#include "Lapack.h"
#include "Special.h"

CheMPS2::TensorL::TensorL( const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * book_up, const SyBookkeeper * book_down ) :
TensorOperator( boundary_index,
                1, //two_j
                1, //n_elec
                Idiff,
                moving_right,
                true, //prime_last
                true, //jw_phase (one 2nd quantized operator)
                book_up,
                book_down ){ }

CheMPS2::TensorL::~TensorL(){ }

void CheMPS2::TensorL::create( TensorT * mps_tensor ){

   clear();
   assert( bk_up == bk_down );

   if ( moving_right ){
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_right( ikappa, mps_tensor, mps_tensor, NULL, NULL ); }
   } else {
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_left( ikappa, mps_tensor, mps_tensor, NULL, NULL ); }
   }

}

void CheMPS2::TensorL::create( TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem ){

   clear();

   if ( moving_right ){
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_right( ikappa, mps_tensor_up, mps_tensor_down, previous, workmem ); }
   } else {
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_left( ikappa, mps_tensor_up, mps_tensor_down, previous, workmem ); }
   }

}

void CheMPS2::TensorL::create_right( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem ){

   const int NRup      = sector_nelec_up[ ikappa ];
   const int NRdown    = NRup + 1;
   const int IRup      = sector_irrep_up[ ikappa ];
   const int IRdown    = Irreps::directProd( IRup, n_irrep );
   const int TwoSRup   = sector_spin_up[ ikappa ];
   const int TwoSRdown = sector_spin_down[ ikappa ];

   int dimRup   = bk_up  ->gCurrentDim( index, NRup,   TwoSRup,   IRup   );
   int dimRdown = bk_down->gCurrentDim( index, NRdown, TwoSRdown, IRdown );

   for ( int geval = 0; geval < 2; geval++ ){
      int NL, TwoSL, IL;
      switch ( geval ){
         case 0:
            NL    = NRup;
            TwoSL = TwoSRup;
            IL    = IRup;
            break;
         case 1:
            NL    = NRup - 1;
            TwoSL = TwoSRdown;
            IL    = IRdown;
            break;
      }

      int dimLup   = bk_up  ->gCurrentDim( index - 1, NL, TwoSL, IL );
      int dimLdown = bk_down->gCurrentDim( index - 1, NL, TwoSL, IL );

      if ( previous == NULL ){
         assert( dimLup == dimLdown );
         if ( dimLup > 0 ){

            double * Tup   = mps_tensor_up  ->gStorage( NL, TwoSL, IL, NRup,   TwoSRup,   IRup   );
            double * Tdown = mps_tensor_down->gStorage( NL, TwoSL, IL, NRdown, TwoSRdown, IRdown );

            char trans   = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if ( geval == 1 ){
               alpha = Special::phase( TwoSRdown - TwoSRup + 1 ) * sqrt( ( TwoSRup + 1.0 ) / ( TwoSRdown + 1 ) );
            }
            double add = 1.0;
            dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, Tdown, &dimLup, &add, storage + kappa2index[ ikappa ], &dimRup );

         }
      } else {
         if (( dimLup > 0 ) && ( dimLdown > 0 )){

            double * Tup   = mps_tensor_up  ->gStorage( NL, TwoSL, IL, NRup,   TwoSRup,   IRup   );
            double * Tdown = mps_tensor_down->gStorage( NL, TwoSL, IL, NRdown, TwoSRdown, IRdown );
            double * Opart =        previous->gStorage( NL, TwoSL, IL, NL,     TwoSL,     IL     );

            char trans   = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if ( geval == 1 ){
               alpha = Special::phase( TwoSRdown - TwoSRup + 1 ) * sqrt( ( TwoSRup + 1.0 ) / ( TwoSRdown + 1 ) );
            }
            double set = 0.0;
            dgemm_( &trans, &notrans, &dimRup, &dimLdown, &dimLup, &alpha, Tup, &dimLup, Opart, &dimLup, &set, workmem, &dimRup );
            double one = 1.0;
            dgemm_( &notrans, &notrans, &dimRup, &dimRdown, &dimLdown, &one, workmem, &dimRup, Tdown, &dimLdown, &one, storage + kappa2index[ ikappa ], &dimRup );

         }
      }
   }

}

void CheMPS2::TensorL::create_left( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem ){

   const int NLup      = sector_nelec_up[ ikappa ];
   const int NLdown    = NLup + 1;
   const int ILup      = sector_irrep_up[ ikappa ];
   const int ILdown    = Irreps::directProd( ILup, n_irrep );
   const int TwoSLup   = sector_spin_up[ ikappa ];
   const int TwoSLdown = sector_spin_down[ ikappa ];

   int dimLup   = bk_up  ->gCurrentDim( index, NLup,   TwoSLup,   ILup   );
   int dimLdown = bk_down->gCurrentDim( index, NLdown, TwoSLdown, ILdown );

   for ( int geval = 0; geval < 2; geval++ ){
      int NR, TwoSR, IR;
      switch ( geval ){
         case 0:
            NR    = NLdown;
            TwoSR = TwoSLdown;
            IR    = ILdown;
            break;
         case 1:
            NR    = NLup + 2;
            TwoSR = TwoSLup;
            IR    = ILup;
            break;
      }

      int dimRup   = bk_up  ->gCurrentDim( index + 1, NR, TwoSR, IR );
      int dimRdown = bk_down->gCurrentDim( index + 1, NR, TwoSR, IR );

      if ( previous == NULL ){
         assert( dimRup == dimRdown );
         if ( dimRup > 0 ){

            double * Tup   = mps_tensor_up  ->gStorage( NLup,   TwoSLup,   ILup,   NR, TwoSR, IR );
            double * Tdown = mps_tensor_down->gStorage( NLdown, TwoSLdown, ILdown, NR, TwoSR, IR );

            char trans   = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if ( geval == 1 ){
               alpha = Special::phase( TwoSLup - TwoSLdown + 1 ) * sqrt( ( TwoSLup + 1.0 ) / ( TwoSLdown + 1 ) );
            }
            double add = 1.0;
            dgemm_( &notrans, &trans, &dimLup, &dimLdown, &dimRup, &alpha, Tup, &dimLup, Tdown, &dimLdown, &add, storage + kappa2index[ ikappa ], &dimLup );

         }
      } else {
         if (( dimRup > 0 ) && ( dimRdown > 0 )){

            double * Tup   = mps_tensor_up  ->gStorage( NLup,   TwoSLup,   ILup,   NR, TwoSR, IR );
            double * Tdown = mps_tensor_down->gStorage( NLdown, TwoSLdown, ILdown, NR, TwoSR, IR );
            double * Opart =        previous->gStorage( NR,     TwoSR,     IR,     NR, TwoSR, IR );

            char trans   = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if ( geval == 1 ){
               alpha = Special::phase( TwoSLup - TwoSLdown + 1 ) * sqrt( ( TwoSLup + 1.0 ) / ( TwoSLdown + 1 ) );
            }
            double set = 0.0;
            dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimRup, &alpha, Tup, &dimLup, Opart, &dimRup, &set, workmem, &dimLup );
            double one = 1.0;
            dgemm_( &notrans, &trans, &dimLup, &dimLdown, &dimRdown, &one, workmem, &dimLup, Tdown, &dimLdown, &one, storage + kappa2index[ ikappa ], &dimLup );

         }
      }
   }

}


