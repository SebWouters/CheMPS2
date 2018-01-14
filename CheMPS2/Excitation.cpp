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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "Excitation.h"
#include "Lapack.h"
#include "Special.h"
#include "Wigner.h"

double CheMPS2::Excitation::matvec( const SyBookkeeper * book_up, const SyBookkeeper * book_down, const int orb1, const int orb2, const double alpha, const double beta, const double gamma, Sobject * S_up, Sobject * S_down, TensorO ** overlaps, TensorL ** regular, TensorL ** trans ){

   const int indx = S_up->gIndex();
   assert( orb1 < orb2 );
   assert( indx == S_down->gIndex() );
   assert( indx >= orb1 );
   assert( indx <  orb2 );
   const int DIM = std::max( std::max( book_up->gMaxDimAtBound( indx ),
                                       book_up->gMaxDimAtBound( indx + 2 ) ),
                           std::max( book_down->gMaxDimAtBound( indx ),
                                     book_down->gMaxDimAtBound( indx + 2 ) ) );
   assert( book_up->gIrrep( orb1 ) == book_up->gIrrep( orb2 ) );

   S_down->prog2symm();

   double inproduct = 0.0;

   #pragma omp parallel reduction(+:inproduct)
   {
      if ( orb1 + 1 == orb2 ){ // The second quantized operators are neighbours
         #pragma omp for schedule(dynamic)
         for ( int dummy = 0; dummy < S_up->gNKappa(); dummy++ ){
            const int ikappa = S_up->gReorder( dummy );
            clear( ikappa, S_up );
            inproduct += neighbours( ikappa, book_up, book_down, alpha, beta, gamma, S_up, S_down );
         }
      } else {
         double * workmem1 = new double[ DIM * DIM ];
         if ( orb1 == indx ){ // All the way at the left
            #pragma omp for schedule(dynamic)
            for ( int dummy = 0; dummy < S_up->gNKappa(); dummy++ ){
               const int ikappa = S_up->gReorder( dummy );
               clear( ikappa, S_up );
                first_left( ikappa, book_up, book_down, alpha, S_up, S_down,   trans[ indx + 1 ] );
               second_left( ikappa, book_up, book_down, beta,  S_up, S_down, regular[ indx + 1 ] );
               inproduct += third_left( ikappa, book_up, book_down, gamma, S_up, S_down, overlaps[ indx + 1 ], workmem1 );
            }
         } else {
            if ( orb2 == indx + 1 ){ // All the way at the right
               #pragma omp for schedule(dynamic)
               for ( int dummy = 0; dummy < S_up->gNKappa(); dummy++ ){
                  const int ikappa = S_up->gReorder( dummy );
                  clear( ikappa, S_up );
                   first_right( ikappa, book_up, book_down, alpha, S_up, S_down,   trans[ indx - 1 ] );
                  second_right( ikappa, book_up, book_down, beta,  S_up, S_down, regular[ indx - 1 ] );
                  inproduct += third_right( ikappa, book_up, book_down, gamma, S_up, S_down, overlaps[ indx - 1 ], workmem1 );
               }
            } else { // Somewhere in the middle
               double * workmem2 = new double[ DIM * DIM ];
               #pragma omp for schedule(dynamic)
               for ( int dummy = 0; dummy < S_up->gNKappa(); dummy++ ){
                  const int ikappa = S_up->gReorder( dummy );
                  clear( ikappa, S_up );
                   first_middle( ikappa, book_up, book_down, alpha, S_up, S_down,   trans[ indx - 1 ],   trans[ indx + 1 ], workmem1 );
                  second_middle( ikappa, book_up, book_down, beta,  S_up, S_down, regular[ indx - 1 ], regular[ indx + 1 ], workmem1 );
                  inproduct += third_middle( ikappa, book_up, book_down, gamma, S_up, S_down, overlaps[ indx - 1 ], overlaps[ indx + 1 ], workmem1, workmem2 );
               }
               delete [] workmem2;
            }
         }
         delete [] workmem1;
      }
   }

   S_up  ->symm2prog();
   // S_down->symm2prog(); S_down is not used anymore afterwards

   return inproduct;

}

void CheMPS2::Excitation::clear( const int ikappa, Sobject * S_up ){

   const int start  = S_up->gKappa2index( ikappa );
   const int stop   = S_up->gKappa2index( ikappa + 1 );
   double * storage = S_up->gStorage();
   for ( int cnt = start; cnt < stop; cnt++ ){ storage[ cnt ] = 0.0; }

}

double CheMPS2::Excitation::neighbours( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, const double beta, const double gamma, Sobject * S_up, Sobject * S_down ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int dimLup   = book_up  ->gCurrentDim( index,     NL, TwoSL, IL );
   const int dimRup   = book_up  ->gCurrentDim( index + 2, NR, TwoSR, IR );
   const int dimLdown = book_down->gCurrentDim( index,     NL, TwoSL, IL );
   const int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( dimLup == dimLdown );
   assert( dimRup == dimRdown );
   assert( book_up->gIrrep( index ) == book_up->gIrrep( index + 1 ) );

   double * block_up = S_up->gStorage() + S_up->gKappa2index( ikappa );
   int size = dimLup * dimRup;
   int inc1 = 1;

   // Add a^+ a
   if ( fabs( alpha ) > 0.0 ){
      if (( TwoJ == 0 ) && ( N1 == 1 ) && ( N2 == 1 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 0, 2, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = sqrt( 2.0 ) * alpha;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 1 ) && ( N1 == 2 ) && ( N2 == 1 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 1, 2, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = -alpha;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 1 ) && ( N1 == 1 ) && ( N2 == 0 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 0, 1, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = alpha;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 0 ) && ( N1 == 2 ) && ( N2 == 0 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 1, 1, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = sqrt( 2.0 ) * alpha;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
   }

   // Add a a^+
   if ( fabs( beta ) > 0.0 ){
      if (( TwoJ == 0 ) && ( N1 == 1 ) && ( N2 == 1 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 2, 0, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = sqrt( 2.0 ) * beta;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 1 ) && ( N1 == 0 ) && ( N2 == 1 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 1, 0, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = beta;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 1 ) && ( N1 == 1 ) && ( N2 == 2 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 2, 1, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = -beta;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
      if (( TwoJ == 0 ) && ( N1 == 0 ) && ( N2 == 2 )){
         double * block_down = S_down->gStorage( NL, TwoSL, IL, 1, 1, TwoJ, NR, TwoSR, IR );
         assert( block_down != NULL );
         double factor = sqrt( 2.0 ) * beta;
         daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
      }
   }

   // Add the constant part
   double * block_down = S_down->gStorage( NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR );
   assert( block_down != NULL );
   if ( fabs( gamma ) > 0.0 ){
      double factor = gamma;
      daxpy_( &size, &factor, block_down, &inc1, block_up, &inc1 );
   }
   const double inproduct = ddot_( &size, block_down, &inc1, block_up, &inc1 );
   return inproduct;

}

void CheMPS2::Excitation::first_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Rtrans ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int IRdown = Irreps::directProd( IR, book_up->gIrrep( index ) );
   const int TwoS2  = (( N2 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   const int dimLdown = book_down->gCurrentDim( index, NL, TwoSL, IL );
   assert( dimLup == dimLdown );
   assert( book_up->gIrrep( index ) == Rtrans->get_irrep() );

   if (( N1 == 1 ) && ( fabs( alpha ) > 0.0 )){ // Cfr. 3K1A in HeffDiagrams3.cpp
      for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
         if (( abs( TwoSL - TwoSRdown ) <= TwoS2 ) && ( TwoSRdown >= 0 )){
            const int memSkappa = S_down->gKappa( NL, TwoSL, IL, 0, N2, TwoS2, NR - 1, TwoSRdown, IRdown );
            if ( memSkappa != -1 ){
               int dimRdown = book_down->gCurrentDim( index + 2, NR - 1, TwoSRdown, IRdown );
               double factor = alpha
                             * Special::phase( TwoSL + TwoSR + TwoJ + 2 * TwoS2 )
                             * sqrt( 1.0 * ( TwoJ + 1 ) * ( TwoSR + 1 ) )
                             * Wigner::wigner6j( TwoS2, TwoJ, 1, TwoSR, TwoSRdown, TwoSL );
               double add = 1.0;
               char notrans = 'N';
               double * block_right = Rtrans->gStorage( NR - 1, TwoSRdown, IRdown, NR, TwoSR, IR );
               double * block_down  = S_down->gStorage() + S_down->gKappa2index( memSkappa );
               double * block_up    = S_up  ->gStorage() + S_up  ->gKappa2index( ikappa );
               dgemm_( &notrans, &notrans, &dimLup, &dimRup, &dimRdown, &factor, block_down, &dimLup, block_right, &dimRdown, &add, block_up, &dimLup );
            }
         }
      }
   }

   if (( N1 == 2 ) && ( fabs( alpha ) > 0.0 )){ // Cfr. 3K1B in HeffDiagrams3.cpp
      for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
         int dimRdown = book_down->gCurrentDim( index + 2, NR - 1, TwoSRdown, IRdown );
         if (( dimRdown > 0 ) && ( TwoSRdown >= 0 )){
            const int TwoJstart = ((( TwoSRdown != TwoSL ) || ( TwoS2 == 0 )) ? 1 + TwoS2 : 0 );
            for ( int TwoJdown = TwoJstart; TwoJdown <= 1 + TwoS2; TwoJdown += 2 ){
               if ( abs( TwoSL - TwoSRdown ) <= TwoJdown ){
                  const int memSkappa = S_down->gKappa( NL, TwoSL, IL, 1, N2, TwoJdown, NR - 1, TwoSRdown, IRdown );
                  if ( memSkappa != -1 ){
                     double factor = alpha
                                   * Special::phase( TwoSL + TwoSR + TwoJdown + 1 + 2 * TwoS2 )
                                   * sqrt( 1.0 * ( TwoJdown + 1 ) * ( TwoSR + 1 ) )
                                   * Wigner::wigner6j( TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL );
                     double add = 1.0;
                     char notrans = 'N';
                     double * block_right = Rtrans->gStorage( NR - 1, TwoSRdown, IRdown, NR, TwoSR, IR );
                     double * block_down  = S_down->gStorage() + S_down->gKappa2index( memSkappa );
                     double * block_up    = S_up  ->gStorage() + S_up  ->gKappa2index( ikappa );
                     dgemm_( &notrans, &notrans, &dimLup, &dimRup, &dimRdown, &factor, block_down, &dimLup, block_right, &dimRdown, &add, block_up, &dimLup );
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Excitation::second_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Rregular ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int IRdown = Irreps::directProd( IR, book_up->gIrrep( index ) );
   const int TwoS2  = (( N2 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   const int dimLdown = book_down->gCurrentDim( index, NL, TwoSL, IL );
   assert( dimLup == dimLdown );
   assert( book_up->gIrrep( index ) == Rregular->get_irrep() );

   if (( N1 == 0 )  && ( fabs( beta ) > 0.0 )){ // Cfr. 3K2A in HeffDiagrams3.cpp
      for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
         int dimRdown = book_down->gCurrentDim( index + 2, NR + 1, TwoSRdown, IRdown );
         if (( dimRdown > 0 ) && ( TwoSRdown >= 0 )){
            const int TwoJstart = ((( TwoSRdown != TwoSL ) || ( TwoS2 == 0 )) ? 1 + TwoS2 : 0 );
            for ( int TwoJdown = TwoJstart; TwoJdown <= 1 + TwoS2; TwoJdown += 2 ){
               if ( abs( TwoSL - TwoSRdown ) <= TwoJdown ){
                  const int memSkappa = S_down->gKappa( NL, TwoSL, IL, 1, N2, TwoJdown, NR + 1, TwoSRdown, IRdown );
                  if ( memSkappa != -1 ){
                     double factor = beta
                                   * Special::phase( TwoSL + TwoSRdown + TwoJdown + 2 * TwoS2 )
                                   * sqrt( 1.0 * ( TwoJdown + 1 ) * ( TwoSRdown + 1 ) )
                                   * Wigner::wigner6j( TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL );
                     double add = 1.0;
                     char notrans = 'N';
                     char trans = 'T';
                     double * block_right = Rregular->gStorage( NR, TwoSR, IR, NR + 1, TwoSRdown, IRdown );
                     double * block_down  = S_down  ->gStorage() + S_down->gKappa2index( memSkappa );
                     double * block_up    = S_up    ->gStorage() + S_up  ->gKappa2index( ikappa );
                     dgemm_( &notrans, &trans, &dimLup, &dimRup, &dimRdown, &factor, block_down, &dimLup, block_right, &dimRup, &add, block_up, &dimLup );
                  }
               }
            }
         }
      }
   }

   if (( N1 == 1 ) && ( fabs( beta ) > 0.0 )){ // Cfr. 3K2B in HeffDiagrams3.cpp
      for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
         if (( abs( TwoSL - TwoSRdown ) <= TwoS2 ) && ( TwoSRdown >= 0 )){
            const int memSkappa = S_down->gKappa( NL, TwoSL, IL, 2, N2, TwoS2, NR + 1, TwoSRdown, IRdown );
            if ( memSkappa != -1 ){
               int dimRdown = book_down->gCurrentDim( index + 2, NR + 1, TwoSRdown, IRdown );
               double factor = beta
                             * Special::phase( TwoSL + TwoSRdown + TwoJ + 1 + 2 * TwoS2 )
                             * sqrt( 1.0 * ( TwoJ + 1 ) * ( TwoSRdown + 1 ) )
                             * Wigner::wigner6j( TwoS2, TwoJ, 1, TwoSR, TwoSRdown, TwoSL );
               double add = 1.0;
               char notrans = 'N';
               char trans = 'T';
               double * block_right = Rregular->gStorage( NR, TwoSR, IR, NR + 1, TwoSRdown, IRdown );
               double * block_down  = S_down  ->gStorage() + S_down->gKappa2index( memSkappa );
               double * block_up    = S_up    ->gStorage() + S_up  ->gKappa2index( ikappa );
               dgemm_( &notrans, &trans, &dimLup, &dimRup, &dimRdown, &factor, block_down, &dimLup, block_right, &dimRup, &add, block_up, &dimLup );
            }
         }
      }
   }

}

double CheMPS2::Excitation::third_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Rovlp, double * workmem ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   int dimLup   = book_up  ->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup   = book_up  ->gCurrentDim( index + 2, NR, TwoSR, IR );
   int dimLdown = book_down->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( dimLup == dimLdown );

   double inproduct = 0.0;
   if ( dimRdown > 0 ){
      double * block_down  = S_down->gStorage( NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR );
      double * block_right = Rovlp ->gStorage( NR, TwoSR, IR, NR, TwoSR, IR );
      char trans   = 'T';
      char notrans = 'N';
      double one = 1.0;
      double set = 0.0;
      dgemm_( &notrans, &trans, &dimLup, &dimRup, &dimRdown, &one, block_down, &dimLup, block_right, &dimRup, &set, workmem, &dimLup );

      double * block_up = S_up->gStorage() + S_up->gKappa2index( ikappa );
      int size = dimLup * dimRup;
      int inc1 = 1;
      if ( fabs( gamma ) > 0.0 ){
         double factor = gamma;
         daxpy_( &size, &factor, workmem, &inc1, block_up, &inc1 );
      }
      inproduct = ddot_( &size, workmem, &inc1, block_up, &inc1 );
   }
   return inproduct;

}

void CheMPS2::Excitation::first_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Ltrans ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int ILdown = Irreps::directProd( IL, book_up->gIrrep( index + 1 ) );
   const int TwoS1  = (( N1 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   const int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( dimRup == dimRdown );
   assert( book_up->gIrrep( index + 1 ) == Ltrans->get_irrep() );

   if (( N2 == 0 ) && ( fabs( alpha ) > 0.0 )){ // Cfr. 3B2A in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         int dimLdown = book_down->gCurrentDim( index, NL - 1, TwoSLdown, ILdown );
         if (( dimLdown > 0 ) && ( TwoSLdown >= 0 )){
            const int TwoJstart = ((( TwoSR != TwoSLdown ) || ( TwoS1 == 0 )) ? 1 + TwoS1 : 0 );
            for ( int TwoJdown = TwoJstart; TwoJdown <= 1 + TwoS1; TwoJdown += 2 ){
               if ( abs( TwoSLdown - TwoSR ) <= TwoJdown ){
                  const int memSkappa = S_down->gKappa( NL - 1, TwoSLdown, ILdown, N1, 1, TwoJdown, NR, TwoSR, IR );
                  if ( memSkappa != -1 ){
                     double factor = alpha
                                   * Special::phase( TwoSLdown + TwoSR + 2 - TwoJdown )
                                   * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoJdown + 1 ) )
                                   * Wigner::wigner6j( TwoJdown, TwoS1, 1, TwoSL, TwoSLdown, TwoSR );
                     double add = 1.0;
                     char notrans = 'N';
                     char trans = 'T';
                     double * block_left = Ltrans->gStorage( NL - 1, TwoSLdown, ILdown, NL, TwoSL, IL );
                     double * block_down = S_down->gStorage() + S_down->gKappa2index( memSkappa );
                     double * block_up   = S_up  ->gStorage() + S_up  ->gKappa2index( ikappa );
                     dgemm_( &trans, &notrans, &dimLup, &dimRup, &dimLdown, &factor, block_left, &dimLdown, block_down, &dimLdown, &add, block_up, &dimLup );
                  }
               }
            }
         }
      }
   }

   if (( N2 == 1 ) && ( fabs( alpha ) > 0.0 )){ // Cfr. 3B2B in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         if (( abs( TwoSLdown - TwoSR ) <= TwoS1 ) && ( TwoSLdown >= 0 )){
            const int memSkappa = S_down->gKappa( NL - 1, TwoSLdown, ILdown, N1, 2, TwoS1, NR, TwoSR, IR );
            if ( memSkappa != -1 ){
               int dimLdown = book_down->gCurrentDim( index, NL - 1, TwoSLdown, ILdown );
               double factor = alpha
                             * Special::phase( TwoSLdown + TwoSR + 3 - TwoJ )
                             * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoJ + 1 ) )
                             * Wigner::wigner6j( TwoS1, TwoJ, 1, TwoSL, TwoSLdown, TwoSR );
               double add = 1.0;
               char notrans = 'N';
               char trans = 'T';
               double * block_left = Ltrans->gStorage( NL - 1, TwoSLdown, ILdown, NL, TwoSL, IL );
               double * block_down = S_down->gStorage() + S_down->gKappa2index( memSkappa );
               double * block_up   = S_up  ->gStorage() + S_up  ->gKappa2index( ikappa );
               dgemm_( &trans, &notrans, &dimLup, &dimRup, &dimLdown, &factor, block_left, &dimLdown, block_down, &dimLdown, &add, block_up, &dimLup );
            }
         }
      }
   }

}

void CheMPS2::Excitation::second_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Lregular ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int ILdown = Irreps::directProd( IL, book_up->gIrrep( index + 1 ) );
   const int TwoS1  = (( N1 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   const int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( dimRup == dimRdown );
   assert( book_up->gIrrep( index + 1 ) == Lregular->get_irrep() );

   if (( N2 == 2 ) && ( fabs( beta ) > 0.0 )){ // Cfr. 3B1A in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         int dimLdown = book_down->gCurrentDim( index, NL + 1, TwoSLdown, ILdown );
         if (( dimLdown > 0 ) && ( TwoSLdown >= 0 )){         
            const int TwoJstart = ((( TwoSR != TwoSLdown ) || ( TwoS1 == 0 )) ? 1 + TwoS1 : 0 );
            for ( int TwoJdown = TwoJstart; TwoJdown <= 1 + TwoS1; TwoJdown += 2 ){
               if ( abs( TwoSLdown - TwoSR ) <= TwoJdown ){
                  const int memSkappa = S_down->gKappa( NL + 1, TwoSLdown, ILdown, N1, 1, TwoJdown, NR, TwoSR, IR );
                  if ( memSkappa != -1 ){
                     double factor = beta
                                   * Special::phase( TwoSL + TwoSR + 3 - TwoJdown )
                                   * sqrt( 1.0 * ( TwoJdown + 1 ) * ( TwoSLdown + 1 ) )
                                   * Wigner::wigner6j( TwoJdown, TwoS1, 1, TwoSL, TwoSLdown, TwoSR );
                     double add = 1.0;
                     char notrans = 'N';
                     double * block_left = Lregular->gStorage( NL, TwoSL, IL, NL + 1, TwoSLdown, ILdown );
                     double * block_down = S_down  ->gStorage() + S_down->gKappa2index( memSkappa );
                     double * block_up   = S_up    ->gStorage() + S_up  ->gKappa2index( ikappa );
                     dgemm_( &notrans, &notrans, &dimLup, &dimRup, &dimLdown, &factor, block_left, &dimLup, block_down, &dimLdown, &add, block_up, &dimLup );
                  }
               }
            }
         }
      }
   }

   if (( N2 == 1 ) && ( fabs( beta ) > 0.0 )){ // Cfr. 3B1B in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         if (( abs( TwoSLdown - TwoSR ) <= TwoS1 ) && ( TwoSLdown >= 0 )){
            const int memSkappa = S_down->gKappa( NL + 1, TwoSLdown, ILdown, N1, 0, TwoS1, NR, TwoSR, IR );
            if ( memSkappa != -1 ){
               int dimLdown = book_down->gCurrentDim( index, NL + 1, TwoSLdown, ILdown );
               double factor = beta
                             * Special::phase( TwoSL + TwoSR + 2 - TwoJ )
                             * sqrt( 1.0 * ( TwoSLdown + 1 ) * ( TwoJ + 1 ) )
                             * Wigner::wigner6j( TwoS1, TwoJ, 1, TwoSL, TwoSLdown, TwoSR );
               double add = 1.0;
               char notrans = 'N';
               double * block_left = Lregular->gStorage( NL, TwoSL, IL, NL + 1, TwoSLdown, ILdown );
               double * block_down = S_down  ->gStorage() + S_down->gKappa2index( memSkappa );
               double * block_up   = S_up    ->gStorage() + S_up  ->gKappa2index( ikappa );
               dgemm_( &notrans, &notrans, &dimLup, &dimRup, &dimLdown, &factor, block_left, &dimLup, block_down, &dimLdown, &add, block_up, &dimLup );
            }
         }
      }
   }

}

double CheMPS2::Excitation::third_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Lovlp, double * workmem ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   int dimLup   = book_up  ->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup   = book_up  ->gCurrentDim( index + 2, NR, TwoSR, IR );
   int dimLdown = book_down->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( dimRup == dimRdown );

   double inproduct = 0.0;
   if ( dimLdown > 0 ){
      double * block_down = S_down->gStorage( NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR );
      double * block_left = Lovlp ->gStorage( NL, TwoSL, IL, NL, TwoSL, IL );
      char notrans = 'N';
      double one = 1.0;
      double set = 0.0;
      dgemm_( &notrans, &notrans, &dimLup, &dimRup, &dimLdown, &one, block_left, &dimLup, block_down, &dimLdown, &set, workmem, &dimLup );

      double * block_up = S_up->gStorage() + S_up->gKappa2index( ikappa );
      int size = dimLup * dimRup;
      int inc1 = 1;
      if ( fabs( gamma ) > 0.0 ){
         double factor = gamma;
         daxpy_( &size, &factor, workmem, &inc1, block_up, &inc1 );
      }
      inproduct = ddot_( &size, workmem, &inc1, block_up, &inc1 );
   }
   return inproduct;

}

void CheMPS2::Excitation::first_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Ltrans, TensorL * Rtrans, double * workmem ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int ILdown = Irreps::directProd( IL, Ltrans->get_irrep() );
   const int IRdown = Irreps::directProd( IR, Rtrans->get_irrep() );
   const int TwoS1  = (( N1 == 1 ) ? 1 : 0 );
   const int TwoS2  = (( N2 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( Ltrans->get_irrep() == Rtrans->get_irrep() );

   if ( fabs( alpha ) > 0.0 ){ // Cfr. 3C2 in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
            if (( abs( TwoSLdown - TwoSRdown ) <= TwoJ ) && ( TwoSLdown >= 0 ) && ( TwoSRdown >= 0 )){
               const int memSkappa = S_down->gKappa( NL - 1, TwoSLdown, ILdown, N1, N2, TwoJ, NR - 1, TwoSRdown, IRdown );
               if ( memSkappa != -1 ){
                  int dimLdown = book_down->gCurrentDim( index,     NL - 1, TwoSLdown, ILdown );
                  int dimRdown = book_down->gCurrentDim( index + 2, NR - 1, TwoSRdown, IRdown );
                  double factor = alpha
                                * Special::phase( TwoSL + TwoSRdown + TwoJ + 1 + 2 * TwoS1 + 2 * TwoS2 )
                                * sqrt( 1.0 * ( TwoSL + 1 ) * ( TwoSR + 1 ) )
                                * Wigner::wigner6j( TwoSL, TwoSR, TwoJ, TwoSRdown, TwoSLdown, 1 );
                  char trans = 'T';
                  char notrans = 'N';
                  double set = 0.0;
                  double one = 1.0;
                  double * block_left  = Ltrans->gStorage( NL - 1, TwoSLdown, ILdown, NL, TwoSL, IL );
                  double * block_right = Rtrans->gStorage( NR - 1, TwoSRdown, IRdown, NR, TwoSR, IR );
                  double * block_down  = S_down->gStorage() + S_down->gKappa2index( memSkappa );
                  double * block_up    = S_up  ->gStorage() + S_up  ->gKappa2index( ikappa );
                  dgemm_( &trans,   &notrans, &dimLup, &dimRdown, &dimLdown, &factor, block_left, &dimLdown, block_down,  &dimLdown, &set, workmem,  &dimLup );
                  dgemm_( &notrans, &notrans, &dimLup, &dimRup,   &dimRdown, &one,    workmem,    &dimLup,   block_right, &dimRdown, &one, block_up, &dimLup );
               }
            }
         }
      }
   }

}

void CheMPS2::Excitation::second_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Lregular, TensorL * Rregular, double * workmem ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   const int ILdown = Irreps::directProd( IL, Lregular->get_irrep() );
   const int IRdown = Irreps::directProd( IR, Rregular->get_irrep() );
   const int TwoS1  = (( N1 == 1 ) ? 1 : 0 );
   const int TwoS2  = (( N2 == 1 ) ? 1 : 0 );

   int dimLup = book_up->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup = book_up->gCurrentDim( index + 2, NR, TwoSR, IR );
   assert( Lregular->get_irrep() == Rregular->get_irrep() );

   if ( fabs( beta ) > 0.0 ){ // Cfr. 3C1 in HeffDiagrams3.cpp
      for ( int TwoSLdown = TwoSL - 1; TwoSLdown <= TwoSL + 1; TwoSLdown += 2 ){
         for ( int TwoSRdown = TwoSR - 1; TwoSRdown <= TwoSR + 1; TwoSRdown += 2 ){
            if (( abs( TwoSLdown - TwoSRdown ) <= TwoJ ) && ( TwoSLdown >= 0 ) && ( TwoSRdown >= 0 )){
               const int memSkappa = S_down->gKappa( NL + 1, TwoSLdown, ILdown, N1, N2, TwoJ, NR + 1, TwoSRdown, IRdown );
               if ( memSkappa != -1 ){
                  int dimLdown = book_down->gCurrentDim( index,     NL + 1, TwoSLdown, ILdown );
                  int dimRdown = book_down->gCurrentDim( index + 2, NR + 1, TwoSRdown, IRdown );
                  double factor = beta
                                * Special::phase( TwoSLdown + TwoSR + TwoJ + 1 + 2 * TwoS1 + 2 * TwoS2 )
                                * sqrt( 1.0 * ( TwoSLdown + 1 ) * ( TwoSRdown + 1 ) )
                                * Wigner::wigner6j( TwoSL, TwoSR, TwoJ, TwoSRdown, TwoSLdown, 1 );
                  char trans = 'T';
                  char notrans = 'N';
                  double set = 0.0;
                  double one = 1.0;
                  double * block_left  = Lregular->gStorage( NL, TwoSL, IL, NL + 1, TwoSLdown, ILdown );
                  double * block_right = Rregular->gStorage( NR, TwoSR, IR, NR + 1, TwoSRdown, IRdown );
                  double * block_down  = S_down  ->gStorage() + S_down->gKappa2index( memSkappa );
                  double * block_up    = S_up    ->gStorage() + S_up  ->gKappa2index( ikappa );
                  dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &factor, block_left, &dimLup, block_down,  &dimLdown, &set, workmem,  &dimLup );
                  dgemm_( &notrans, &trans,   &dimLup, &dimRup,   &dimRdown, &one,    workmem,    &dimLup, block_right, &dimRup,   &one, block_up, &dimLup );
               }
            }
         }
      }
   }

}

double CheMPS2::Excitation::third_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Lovlp, TensorO * Rovlp, double * workmem1, double * workmem2 ){

   const int index = S_up->gIndex();
   const int TwoSL = S_up->gTwoSL( ikappa );
   const int TwoSR = S_up->gTwoSR( ikappa );
   const int TwoJ = S_up->gTwoJ( ikappa );
   const int NL = S_up->gNL( ikappa );
   const int NR = S_up->gNR( ikappa );
   const int IL = S_up->gIL( ikappa );
   const int IR = S_up->gIR( ikappa );
   const int N1 = S_up->gN1( ikappa );
   const int N2 = S_up->gN2( ikappa );

   int dimLup   = book_up  ->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRup   = book_up  ->gCurrentDim( index + 2, NR, TwoSR, IR );
   int dimLdown = book_down->gCurrentDim( index,     NL, TwoSL, IL );
   int dimRdown = book_down->gCurrentDim( index + 2, NR, TwoSR, IR );

   double inproduct = 0.0;
   if (( dimLdown > 0 ) && ( dimRdown > 0 )){
      double * block_down  = S_down->gStorage( NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR );
      double * block_left  = Lovlp ->gStorage( NL, TwoSL, IL, NL, TwoSL, IL );
      double * block_right = Rovlp ->gStorage( NR, TwoSR, IR, NR, TwoSR, IR );
      char trans = 'T';
      char notrans = 'N';
      double one = 1.0;
      double set = 0.0;
      dgemm_( &notrans, &notrans, &dimLup, &dimRdown, &dimLdown, &one, block_left, &dimLup, block_down,  &dimLdown, &set, workmem1, &dimLup );
      dgemm_( &notrans, &trans,   &dimLup, &dimRup,   &dimRdown, &one, workmem1,   &dimLup, block_right, &dimRup,   &set, workmem2, &dimLup );

      double * block_up = S_up->gStorage() + S_up->gKappa2index( ikappa );
      int size = dimLup * dimRup;
      int inc1 = 1;
      if ( fabs( gamma ) > 0.0 ){
         double factor = gamma;
         daxpy_( &size, &factor, workmem2, &inc1, block_up, &inc1 );
      }
      inproduct = ddot_( &size, workmem2, &inc1, block_up, &inc1 );
   }
   return inproduct;

}


