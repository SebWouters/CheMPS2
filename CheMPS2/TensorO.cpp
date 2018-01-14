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
#include <algorithm>

#include "TensorO.h"
#include "Lapack.h"

CheMPS2::TensorO::TensorO( const int boundary_index, const bool moving_right, const SyBookkeeper * book_up, const SyBookkeeper * book_down ) :
TensorOperator( boundary_index,
                0, //two_j
                0, //n_elec
                0, //n_irrep
                moving_right,
                true,  //prime_last (doesn't matter for spin-0 tensors)
                false, //jw_phase (no operators)
                book_up,
                book_down ){ }

CheMPS2::TensorO::~TensorO(){ }

void CheMPS2::TensorO::update_ownmem( TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous ){

   clear();

   if ( moving_right ){

      const int dimL = std::max( bk_up->gMaxDimAtBound( index - 1 ), bk_down->gMaxDimAtBound( index - 1 ) );
      const int dimR = std::max( bk_up->gMaxDimAtBound( index ),     bk_down->gMaxDimAtBound( index )     );

      #pragma omp parallel
      {
         double * workmem = new double[ dimL * dimR ];

         #pragma omp for schedule(dynamic)
         for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
            update_moving_right( ikappa, previous, mps_tensor_up, mps_tensor_down, workmem );
         }

         delete [] workmem;
      }
   } else {

      const int dimL = std::max( bk_up->gMaxDimAtBound( index ),     bk_down->gMaxDimAtBound( index )     );
      const int dimR = std::max( bk_up->gMaxDimAtBound( index + 1 ), bk_down->gMaxDimAtBound( index + 1 ) );

      #pragma omp parallel
      {
         double * workmem = new double[ dimL * dimR ];

         #pragma omp for schedule(dynamic)
         for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
            update_moving_left( ikappa, previous, mps_tensor_up, mps_tensor_down, workmem );
         }

         delete [] workmem;
      }
   }

}

void CheMPS2::TensorO::create( TensorT * mps_tensor_up, TensorT * mps_tensor_down ){

   clear();

   if ( moving_right ){
      #pragma omp parallel for schedule(dynamic)
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_right( ikappa, mps_tensor_up, mps_tensor_down ); }
   } else {
      #pragma omp parallel for schedule(dynamic)
      for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){ create_left( ikappa, mps_tensor_up, mps_tensor_down ); }
   }

}

void CheMPS2::TensorO::create_right( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down ){

   const int NR    = sector_nelec_up[ ikappa ];
   const int IR    = sector_irrep_up[ ikappa ];
   const int TwoSR = sector_spin_up [ ikappa ];

   int dimRup   =   bk_up->gCurrentDim( index, NR, TwoSR, IR );
   int dimRdown = bk_down->gCurrentDim( index, NR, TwoSR, IR );

   for ( int geval = 0; geval < 4; geval++ ){
      int IL, TwoSL, NL;
      switch ( geval ){
         case 0:
            NL    = NR;
            TwoSL = TwoSR;
            IL    = IR;
            break;
         case 1:
            NL    = NR - 2;
            TwoSL = TwoSR;
            IL    = IR;
            break;
         case 2:
            NL    = NR - 1;
            TwoSL = TwoSR - 1;
            IL    = Irreps::directProd( IR, bk_up->gIrrep( index - 1 ) ); // bk_up and bk_down treat the same orbitals and ordering
            break;
         case 3:
            NL    = NR - 1;
            TwoSL = TwoSR + 1;
            IL    = Irreps::directProd( IR, bk_up->gIrrep( index - 1 ) ); // bk_up and bk_down treat the same orbitals and ordering
            break;
      }

      int dimLup   =   bk_up->gCurrentDim( index - 1, NL, TwoSL, IL );
      int dimLdown = bk_down->gCurrentDim( index - 1, NL, TwoSL, IL );
      if (( dimLup > 0 ) && ( dimLdown > 0 ) && ( dimLup == dimLdown )){

         double alpha = 1.0;
         double beta  = 1.0; //add
         char trans   = 'T';
         char notrans = 'N';
         double * Tup   =   mps_tensor_up->gStorage( NL, TwoSL, IL, NR, TwoSR, IR );
         double * Tdown = mps_tensor_down->gStorage( NL, TwoSL, IL, NR, TwoSR, IR );
         dgemm_( &trans, &notrans, &dimRup, &dimRdown, &dimLup, &alpha, Tup, &dimLup, Tdown, &dimLdown, &beta, storage + kappa2index[ ikappa ], &dimRup );

      }
   }

}

void CheMPS2::TensorO::create_left( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down ){

   const int NL    = sector_nelec_up[ ikappa ];
   const int IL    = sector_irrep_up[ ikappa ];
   const int TwoSL = sector_spin_up [ ikappa ];

   int dimLup   =   bk_up->gCurrentDim( index, NL, TwoSL, IL );
   int dimLdown = bk_down->gCurrentDim( index, NL, TwoSL, IL );

   for ( int geval = 0; geval < 4; geval++ ){
      int IR, TwoSR, NR;
      switch ( geval ){
         case 0:
            NR    = NL;
            TwoSR = TwoSL;
            IR    = IL;
            break;
         case 1:
            NR    = NL + 2;
            TwoSR = TwoSL;
            IR    = IL;
            break;
         case 2:
            NR    = NL + 1;
            TwoSR = TwoSL - 1;
            IR    = Irreps::directProd( IL, bk_up->gIrrep( index ) ); // bk_up and bk_down treat the same orbitals and ordering
            break;
         case 3:
            NR    = NL + 1;
            TwoSR = TwoSL + 1;
            IR    = Irreps::directProd( IL, bk_up->gIrrep( index ) ); // bk_up and bk_down treat the same orbitals and ordering
            break;
      }

      int dimRup   =   bk_up->gCurrentDim( index + 1, NR, TwoSR, IR );
      int dimRdown = bk_down->gCurrentDim( index + 1, NR, TwoSR, IR );
      if (( dimRup > 0 ) && ( dimRdown > 0 ) && ( dimRup == dimRdown )){

         double alpha = (( geval > 1 ) ? (( TwoSR + 1.0 ) / ( TwoSL + 1 )) : 1.0 );
         double beta  = 1.0; //add
         char trans   = 'T';
         char notrans = 'N';
         double * Tup   =   mps_tensor_up->gStorage( NL, TwoSL, IL, NR, TwoSR, IR );
         double * Tdown = mps_tensor_down->gStorage( NL, TwoSL, IL, NR, TwoSR, IR );
         dgemm_( &notrans, &trans, &dimLup, &dimLdown, &dimRup, &alpha, Tup, &dimLup, Tdown, &dimLdown, &beta, storage + kappa2index[ ikappa ], &dimLup );

      }
   }

}


