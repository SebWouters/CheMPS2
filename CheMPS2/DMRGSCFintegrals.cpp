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

#include <assert.h>

#include "DMRGSCFintegrals.h"

CheMPS2::DMRGSCFintegrals::DMRGSCFintegrals(DMRGSCFindices * iHandler){

   numberOfIrreps = iHandler->getNirreps();
   NCORE    = new int[ numberOfIrreps ];
   NVIRTUAL = new int[ numberOfIrreps ];
   NTOTAL   = new int[ numberOfIrreps ];
   
   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
      NCORE[    irrep ] = iHandler->getNOCC(  irrep ) + iHandler->getNDMRG( irrep );
      NVIRTUAL[ irrep ] = iHandler->getNVIRT( irrep );
      NTOTAL[   irrep ] = iHandler->getNORB(  irrep );
   }

   coulomb_size   = calcNumCoulombElements(  true );
   exchange_size  = calcNumExchangeElements( true );
   coulomb_array  = new double[  coulomb_size ];
   exchange_array = new double[ exchange_size ];
   
}

long long CheMPS2::DMRGSCFintegrals::calcNumCoulombElements(const bool allocate){

   // The object sizes
   long long theSize = 0;
   
   if (allocate){ coulomb_ptr = new long long***[ numberOfIrreps ]; }
   for (int I_cc = 0; I_cc < numberOfIrreps; I_cc++){ // Loop the irrep I_cc = I_c1 x I_c2 = I_a1 x I_a2
      if (allocate){ coulomb_ptr[ I_cc ] = new long long**[ numberOfIrreps ]; }
      for (int I_c1 = 0; I_c1 < numberOfIrreps; I_c1++){
         const int I_c2 = Irreps::directProd( I_cc , I_c1 );
         if ( ( NCORE[ I_c1 ] > 0 ) && ( NCORE[ I_c2 ] > 0 ) && ( I_c1 <= I_c2 ) ){
            if (allocate){ coulomb_ptr[ I_cc ][ I_c1 ] = new long long*[ numberOfIrreps ]; }
            for (int I_a1 = 0; I_a1 < numberOfIrreps; I_a1++){
               const int I_a2 = Irreps::directProd( I_cc, I_a1 );
               if ( ( NTOTAL[ I_a1 ] > 0 ) && ( NTOTAL[ I_a2 ] > 0 ) && ( I_a1 <= I_a2 ) ){
                  if ( I_cc == 0 ){ // I_c1 == I_c2 and I_a1 == I_a2
                     if (allocate){
                        const long long coretriangle = ( NCORE[  I_c1 ] * ( NCORE[  I_c1 ] + 1 ) ) / 2;
                        const long long  alltriangle = ( NTOTAL[ I_a1 ] * ( NTOTAL[ I_a1 ] + 1 ) ) / 2;
                        coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ] = new long long[ coretriangle ];
                        for (int combinedcore = 0; combinedcore < coretriangle; combinedcore++){
                           coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ][ combinedcore ] = theSize;
                           theSize += alltriangle;
                        }
                     } else { delete [] coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ]; }
                  } else { // I_c1 < I_c2 and I_a1 < I_a2
                     if (allocate){
                        const long long coresquare = NCORE[  I_c1 ] * NCORE[  I_c2 ];
                        const long long  allsquare = NTOTAL[ I_a1 ] * NTOTAL[ I_a2 ];
                        coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ] = new long long[ coresquare ];
                        for (int combinedcore = 0; combinedcore < coresquare; combinedcore++){
                           coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ][ combinedcore ] = theSize;
                           theSize += allsquare;
                        }
                     } else { delete [] coulomb_ptr[ I_cc ][ I_c1 ][ I_a1 ]; }
                  }
               }
            }
            if (!allocate){ delete [] coulomb_ptr[ I_cc ][ I_c1 ]; }
         }
      }
      if (!allocate){ delete [] coulomb_ptr[ I_cc ]; }
   }
   if (!allocate){ delete [] coulomb_ptr; }
   
   return theSize;
   
}

long long CheMPS2::DMRGSCFintegrals::calcNumExchangeElements(const bool allocate){

   // The object sizes
   long long theSize = 0;
   
   if (allocate){ exchange_ptr = new long long***[ numberOfIrreps ]; }
   for (int I_cc = 0; I_cc < numberOfIrreps; I_cc++){ // Loop the irrep I_cc = I_c1 x I_c2 = I_v1 x I_v2
      if (allocate){ exchange_ptr[ I_cc ] = new long long**[ numberOfIrreps ]; }
      for (int I_c1 = 0; I_c1 < numberOfIrreps; I_c1++){
         const int I_c2 = Irreps::directProd( I_cc , I_c1 );
         if ( ( NCORE[ I_c1 ] > 0 ) && ( NCORE[ I_c2 ] > 0 ) && ( I_c1 <= I_c2 ) ){
            if (allocate){ exchange_ptr[ I_cc ][ I_c1 ] = new long long*[ numberOfIrreps ]; }
            for (int I_v1 = 0; I_v1 < numberOfIrreps; I_v1++){
               const int I_v2 = Irreps::directProd( I_cc, I_v1 );
               if ( ( NTOTAL[ I_v1 ] > 0 ) && ( NTOTAL[ I_v2 ] > 0 ) ){ // Here no I_v1 <= I_v2 !!
                  const long long virtualsquare = NVIRTUAL[ I_v1 ] * NVIRTUAL[ I_v2 ];
                  if ( I_cc == 0 ){ // I_c1 == I_c2 and I_v1 == I_v2
                     if (allocate){
                        const long long coretriangle = ( NCORE[ I_c1 ] * ( NCORE[ I_c1 ] + 1 ) ) / 2;
                        exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ] = new long long[ coretriangle ];
                        for (int combinedcore = 0; combinedcore < coretriangle; combinedcore++){
                           exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ][ combinedcore ] = theSize;
                           theSize += virtualsquare;
                        }
                     } else { delete [] exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ]; }
                  } else { // I_c1 < I_c2 and I_v1 != I_v2
                     if (allocate){
                        const long long coresquare = NCORE[ I_c1 ] * NCORE[ I_c2 ];
                        exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ] = new long long[ coresquare ];
                        for (int combinedcore = 0; combinedcore < coresquare; combinedcore++){
                           exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ][ combinedcore ] = theSize;
                           theSize += virtualsquare;
                        }
                     } else { delete [] exchange_ptr[ I_cc ][ I_c1 ][ I_v1 ]; }
                  }
               }
            }
            if (!allocate){ delete [] exchange_ptr[ I_cc ][ I_c1 ]; }
         }
      }
      if (!allocate){ delete [] exchange_ptr[ I_cc ]; }
   }
   if (!allocate){ delete [] exchange_ptr; }
   
   return theSize;
   
}

CheMPS2::DMRGSCFintegrals::~DMRGSCFintegrals(){

   delete [] coulomb_array;
   delete [] exchange_array;

   calcNumCoulombElements(  false );
   calcNumExchangeElements( false );

   delete [] NCORE;
   delete [] NVIRTUAL;
   delete [] NTOTAL;
   
}

void CheMPS2::DMRGSCFintegrals::clear(){

   for (long long counter = 0; counter < coulomb_size;  counter++){  coulomb_array[ counter ] = 0.0; }
   for (long long counter = 0; counter < exchange_size; counter++){ exchange_array[ counter ] = 0.0; }

}

long long CheMPS2::DMRGSCFintegrals::get_coulomb_ptr( const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2 ) const{

   const int Icc = Irreps::directProd( Ic1, Ic2 );
   assert( Icc == Irreps::directProd( Ia1, Ia2 ) );
   
   if ( Icc == 0 ){ // Ic1 == Ic2 and Ia1 == Ia2
      const int index_c = ( c1 <= c2 ) ? c1 + (c2 * ( c2 + 1 ))/2 : c2 + (c1 * ( c1 + 1 ))/2 ;
      const int index_a = ( a1 <= a2 ) ? a1 + (a2 * ( a2 + 1 ))/2 : a2 + (a1 * ( a1 + 1 ))/2 ;
      return coulomb_ptr[ Icc ][ Ic1 ][ Ia1 ][ index_c ] + index_a ;
   }
   
   // Ic1 != Ic2 and Ia1 != Ia2
   const int irrep_c = ( Ic1 < Ic2 ) ? Ic1 : Ic2 ;
   const int irrep_a = ( Ia1 < Ia2 ) ? Ia1 : Ia2 ;
   const int index_c = ( Ic1 < Ic2 ) ? c1 + NCORE[  Ic1 ] * c2 : c2 + NCORE[  Ic2 ] * c1 ;
   const int index_a = ( Ia1 < Ia2 ) ? a1 + NTOTAL[ Ia1 ] * a2 : a2 + NTOTAL[ Ia2 ] * a1 ;
   return coulomb_ptr[ Icc ][ irrep_c ][ irrep_a ][ index_c ] + index_a ;

}

void CheMPS2::DMRGSCFintegrals::set_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2, const double val){
   
   coulomb_array[ get_coulomb_ptr( Ic1, Ic2, Ia1, Ia2, c1, c2, a1, a2 ) ] = val;

}

void CheMPS2::DMRGSCFintegrals::add_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2, const double val){
   
   coulomb_array[ get_coulomb_ptr( Ic1, Ic2, Ia1, Ia2, c1, c2, a1, a2 ) ] += val;

}

double CheMPS2::DMRGSCFintegrals::get_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2) const{
   
   return coulomb_array[ get_coulomb_ptr( Ic1, Ic2, Ia1, Ia2, c1, c2, a1, a2 ) ];

}

long long CheMPS2::DMRGSCFintegrals::get_exchange_ptr( const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2 ) const{

   const int Icc = Irreps::directProd( Ic1, Ic2 );
   assert( Icc == Irreps::directProd( Iv1, Iv2 ) );
   
   if ( Icc == 0 ){ // Ic1 == Ic2 and Iv1 == Iv2
   
      if ( c1 <= c2 ){
         return exchange_ptr[ Icc ][ Ic1 ][ Iv1 ][ c1 + (c2 * ( c2 + 1 ))/2 ] + v1 - NCORE[ Iv1 ] + NVIRTUAL[ Iv1 ] * ( v2 - NCORE[ Iv2 ] ) ;
      } else {
         return exchange_ptr[ Icc ][ Ic2 ][ Iv2 ][ c2 + (c1 * ( c1 + 1 ))/2 ] + v2 - NCORE[ Iv2 ] + NVIRTUAL[ Iv2 ] * ( v1 - NCORE[ Iv1 ] ) ;
      }
   
   } else { // Ic1 != Ic2 
   
      if ( Ic1 < Ic2 ){
         return exchange_ptr[ Icc ][ Ic1 ][ Iv1 ][ c1 + NCORE[ Ic1 ] * c2 ] + v1 - NCORE[ Iv1 ] + NVIRTUAL[ Iv1 ] * ( v2 - NCORE[ Iv2 ] ) ;
      } else {
         return exchange_ptr[ Icc ][ Ic2 ][ Iv2 ][ c2 + NCORE[ Ic2 ] * c1 ] + v2 - NCORE[ Iv2 ] + NVIRTUAL[ Iv2 ] * ( v1 - NCORE[ Iv1 ] ) ;
      }
   
   }
   
   return -1;

}

void CheMPS2::DMRGSCFintegrals::set_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2, const double val){
   
   exchange_array[ get_exchange_ptr( Ic1, Ic2, Iv1, Iv2, c1, c2, v1, v2 ) ] = val;

}

void CheMPS2::DMRGSCFintegrals::add_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2, const double val){
   
   exchange_array[ get_exchange_ptr( Ic1, Ic2, Iv1, Iv2, c1, c2, v1, v2 ) ] += val;

}

double CheMPS2::DMRGSCFintegrals::get_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2) const{
   
   return exchange_array[ get_exchange_ptr( Ic1, Ic2, Iv1, Iv2, c1, c2, v1, v2 ) ];

}

double CheMPS2::DMRGSCFintegrals::FourIndexAPI(const int I1, const int I2, const int I3, const int I4, const int index1, const int index2, const int index3, const int index4) const{

   assert( Irreps::directProd( I1, I2 ) == Irreps::directProd( I3, I4 ) );
   
   const bool core1 = ( index1 < NCORE[I1] ) ? true : false;
   const bool core2 = ( index2 < NCORE[I2] ) ? true : false;
   const bool core3 = ( index3 < NCORE[I3] ) ? true : false;
   const bool core4 = ( index4 < NCORE[I4] ) ? true : false;
   
   const int numCore = ( ( core1 ) ? 1 : 0 ) + ( ( core2 ) ? 1 : 0 ) + ( ( core3 ) ? 1 : 0 ) + ( ( core4 ) ? 1 : 0 );
   assert( numCore >= 2 );
   
   if ( numCore == 4 ){
      return get_coulomb(I1, I3, I2, I4, index1, index3, index2, index4);
   }
   
   if ( numCore == 3 ){
      if (( !core1 ) || ( !core3 )){ return get_coulomb(I2, I4, I1, I3, index2, index4, index1, index3); }
      if (( !core2 ) || ( !core4 )){ return get_coulomb(I1, I3, I2, I4, index1, index3, index2, index4); }
   }
   
   if ( numCore == 2 ){
      if ( !core1 ){
         if ( !core2 ){ return get_exchange(I3, I4, I1, I2, index3, index4, index1, index2); }
         if ( !core3 ){ return get_coulomb( I2, I4, I1, I3, index2, index4, index1, index3); }
         if ( !core4 ){ return get_exchange(I3, I2, I1, I4, index3, index2, index1, index4); }
      }
      if ( !core2 ){
         if ( !core3 ){ return get_exchange(I4, I1, I2, I3, index4, index1, index2, index3); }
         if ( !core4 ){ return get_coulomb( I1, I3, I2, I4, index1, index3, index2, index4); }
      }
      return get_exchange(I1, I2, I3, I4, index1, index2, index3, index4);
   }

   assert( 0 == 1 );
   return 0.0;

}

