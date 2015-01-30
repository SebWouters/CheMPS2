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

#include <climits>
#include <assert.h>
#include <iostream>
#include <math.h>

using std::cout;
using std::endl;

#include "FCI.h"
#include "Irreps.h"
#include "Lapack.h"
#include "Davidson.h"

CheMPS2::FCI::FCI(Hamiltonian * Ham, const unsigned int theNel_up, const unsigned int theNel_down, const int TargetIrrep_in, const double maxMemWorkMB_in, const int FCIverbose_in){

   // Copy the basic information
   FCIverbose   = FCIverbose_in;
   maxMemWorkMB = maxMemWorkMB_in;
   L = Ham->getL();
   assert( theNel_up    <= L );
   assert( theNel_down  <= L );
   assert( maxMemWorkMB >  0.0 );
   Nel_up   = theNel_up;
   Nel_down = theNel_down;
   
   // Construct the irrep product table and the list with the orbitals irreps
   CheMPS2::Irreps myIrreps( Ham->getNGroup() );
   NumIrreps         = myIrreps.getNumberOfIrreps();
   TargetIrrep       = TargetIrrep_in;
   orb2irrep         = new int[ L ];
   for (unsigned int orb = 0; orb < L; orb++){ orb2irrep[ orb ] = Ham->getOrbitalIrrep( orb ); }

   /* Copy the Hamiltonian over:
         G_ij = T_ij - 0.5 \sum_k <ik|kj> and ERI_{ijkl} = <ij|kl>
         <ij|kl> is the electron repulsion integral, int dr1 dr2 i(r1) j(r1) k(r2) l(r2) / |r1-r2| */
   Econstant = Ham->getEconst();
   Gmat = new double[ L * L ];
   ERI  = new double[ L * L * L * L ];
   for (unsigned int orb1 = 0; orb1 < L; orb1++){
      for (unsigned int orb2 = 0; orb2 < L; orb2++){
         double tempvar = 0.0;
         for (unsigned int orb3 = 0; orb3 < L; orb3++){
            tempvar += Ham->getVmat( orb1, orb3, orb3, orb2 );
            for (unsigned int orb4 = 0; orb4 < L; orb4++){
               // CheMPS2::Hamiltonian uses physics notation ; ERI chemists notation.
               ERI[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ] = Ham->getVmat( orb1 , orb3 , orb2 , orb4 );
            }
         }
         Gmat[ orb1 + L * orb2 ] = Ham->getTmat( orb1 , orb2 ) - 0.5 * tempvar;
      }
   }
   
   // Set all other internal variables
   StartupCountersVsBitstrings();
   StartupLookupTables();
   StartupIrrepCenter();

}

CheMPS2::FCI::~FCI(){
   
   // FCI::FCI
   delete [] orb2irrep;
   delete [] Gmat;
   delete [] ERI;
   
   // FCI::StartupCountersVsBitstrings
   for ( unsigned int irrep=0; irrep<NumIrreps; irrep++ ){
      delete [] str2cnt_up[irrep];
      delete [] str2cnt_down[irrep];
      delete [] cnt2str_up[irrep];
      delete [] cnt2str_down[irrep];
   }
   delete [] str2cnt_up;
   delete [] str2cnt_down;
   delete [] cnt2str_up;
   delete [] cnt2str_down;
   delete [] numPerIrrep_up;
   delete [] numPerIrrep_down;
   
   // FCI::StartupLookupTables
   for ( unsigned int irrep=0; irrep<NumIrreps; irrep++ ){
      delete [] lookup_cnt_alpha[irrep];
      delete [] lookup_cnt_beta[irrep];
      delete [] lookup_irrep_alpha[irrep];
      delete [] lookup_irrep_beta[irrep];
      delete [] lookup_sign_alpha[irrep];
      delete [] lookup_sign_beta[irrep];
   }
   delete [] lookup_cnt_alpha;
   delete [] lookup_cnt_beta;
   delete [] lookup_irrep_alpha;
   delete [] lookup_irrep_beta;
   delete [] lookup_sign_alpha;
   delete [] lookup_sign_beta;
   
   // FCI::StartupIrrepCenter
   for ( unsigned int irrep=0; irrep<NumIrreps; irrep++ ){
      delete [] irrep_center_crea_orb[irrep];
      delete [] irrep_center_anni_orb[irrep];
      delete [] irrep_center_jumps[irrep];
   }
   delete [] irrep_center_crea_orb;
   delete [] irrep_center_anni_orb;
   delete [] irrep_center_jumps;
   delete [] irrep_center_num;
   delete [] HXVworksmall;
   delete [] HXVworkbig1;
   delete [] HXVworkbig2;

}

void CheMPS2::FCI::StartupCountersVsBitstrings(){

   // Can you represent the alpha and beta Slater determinants as unsigned integers?
   assert( L <= CHAR_BIT * sizeof(unsigned int) );
   
   // Variable which is only needed here: 2^L
   unsigned int TwoPowL = 1;
   for (unsigned int orb = 0; orb < L; orb++){ TwoPowL *= 2; }

   // Create the required arrays to perform the conversions between counters and bitstrings
   numPerIrrep_up     = new unsigned int[ NumIrreps ];
   numPerIrrep_down   = new unsigned int[ NumIrreps ];
   str2cnt_up         = new int*[ NumIrreps ];
   str2cnt_down       = new int*[ NumIrreps ];
   cnt2str_up         = new unsigned int*[ NumIrreps ];
   cnt2str_down       = new unsigned int*[ NumIrreps ];
   
   for (unsigned int irrep = 0; irrep < NumIrreps; irrep++){
      numPerIrrep_up  [ irrep ] = 0;
      numPerIrrep_down[ irrep ] = 0;
      str2cnt_up  [ irrep ] = new int[ TwoPowL ];
      str2cnt_down[ irrep ] = new int[ TwoPowL ];
   }
   
   int * bits = new int[ L ]; // Temporary helper array
   
   // Loop over all allowed bit strings in the spinless fermion Fock space
   for (unsigned int bitstring = 0; bitstring < TwoPowL; bitstring++){
   
      // Find the number of particles and the irrep which correspond to each basis vector
      str2bits( L , bitstring , bits );
      unsigned int Nparticles = 0;
      int Irrep = 0;
      for (unsigned int orb=0; orb<L; orb++){
         if ( bits[orb] ){
            Nparticles++;
            Irrep = getIrrepProduct( Irrep , getOrb2Irrep( orb ) );
         }
      }
      
      // If allowed: set the corresponding str2cnt to the correct counter and keep track of the number of allowed vectors
      for ( unsigned int irr = 0; irr < NumIrreps; irr++ ){
         str2cnt_up  [ irr ][ bitstring ] = -1;
         str2cnt_down[ irr ][ bitstring ] = -1;
      }
      if ( Nparticles == Nel_up ){
         str2cnt_up[ Irrep ][ bitstring ] = numPerIrrep_up[ Irrep ];
         numPerIrrep_up[ Irrep ]++;
      }
      if ( Nparticles == Nel_down ){
         str2cnt_down[ Irrep ][ bitstring ] = numPerIrrep_down[ Irrep ];
         numPerIrrep_down[ Irrep ]++;
      }
   
   }
   
   // Fill the reverse info array: cnt2str
   for ( unsigned int irrep = 0; irrep < NumIrreps; irrep++ ){
   
      if ( FCIverbose>1 ){
         cout << "FCI::Startup : For irrep " << irrep << " there are " << numPerIrrep_up  [ irrep ] << " alpha Slater determinants and "
                                                                       << numPerIrrep_down[ irrep ] <<  " beta Slater determinants." << endl;
      }
      
      cnt2str_up  [ irrep ] = new unsigned int[ numPerIrrep_up  [ irrep ] ];
      cnt2str_down[ irrep ] = new unsigned int[ numPerIrrep_down[ irrep ] ];
      for (unsigned int bitstring = 0; bitstring < TwoPowL; bitstring++){
         if ( str2cnt_up  [ irrep ][ bitstring ] != -1 ){ cnt2str_up  [ irrep ][ str2cnt_up  [ irrep ][ bitstring ] ] = bitstring; }
         if ( str2cnt_down[ irrep ][ bitstring ] != -1 ){ cnt2str_down[ irrep ][ str2cnt_down[ irrep ][ bitstring ] ] = bitstring; }
      }
   
   }
   
   delete [] bits; // Delete temporary helper array

}

void CheMPS2::FCI::StartupLookupTables(){

   // Create a bunch of stuff
   lookup_cnt_alpha   = new int*[ NumIrreps ];
   lookup_cnt_beta    = new int*[ NumIrreps ];
   lookup_irrep_alpha = new int*[ NumIrreps ];
   lookup_irrep_beta  = new int*[ NumIrreps ];
   lookup_sign_alpha  = new int*[ NumIrreps ];
   lookup_sign_beta   = new int*[ NumIrreps ];
   
   int * bits = new int[ L ]; // Temporary helper array
   
   // Quick lookup tables for " sign | new > = E^spinproj_{ij} | old >
   for ( unsigned int irrep = 0; irrep < NumIrreps; irrep++ ){
      
      lookup_cnt_alpha  [ irrep ] = new int[ L * L * numPerIrrep_up[ irrep ] ];
      lookup_irrep_alpha[ irrep ] = new int[ L * L * numPerIrrep_up[ irrep ] ];
      lookup_sign_alpha [ irrep ] = new int[ L * L * numPerIrrep_up[ irrep ] ];
      
      for ( unsigned int cnt_new_alpha = 0; cnt_new_alpha < numPerIrrep_up[ irrep ]; cnt_new_alpha++ ){
      
         for ( unsigned int i = 0; i < L; i++ ){
            for ( unsigned int j = 0; j < L; j++){
               // Check for the sign. If no check for sign, you multiply with sign 0 and everything should be OK...
               lookup_cnt_alpha  [irrep][ i + L * ( j + L * cnt_new_alpha ) ] = 0;
               lookup_irrep_alpha[irrep][ i + L * ( j + L * cnt_new_alpha ) ] = 0;
               lookup_sign_alpha [irrep][ i + L * ( j + L * cnt_new_alpha ) ] = 0;
            }
         }
      
         str2bits( L , cnt2str_up[ irrep ][ cnt_new_alpha ] , bits );
      
         int phase_creator = 1;
         for ( unsigned int creator = 0; creator < L; creator++ ){
            if ( bits[ creator ] ){
               bits[ creator ] = 0;
                  
               int phase_annihilator = 1;
               for ( unsigned int annihilator = 0; annihilator < L; annihilator++ ){
                  if ( !(bits[ annihilator ]) ){
                     bits[ annihilator ] = 1;
                     
                     const int irrep_old = getIrrepProduct( irrep , getIrrepProduct( getOrb2Irrep( creator ) , getOrb2Irrep( annihilator ) ) );
                     const int cnt_old = str2cnt_up[ irrep_old ][ bits2str( L , bits ) ];
                     const int phase = phase_creator * phase_annihilator;
                     
                     lookup_cnt_alpha  [irrep][ creator + L * ( annihilator + L * cnt_new_alpha ) ] = cnt_old;
                     lookup_irrep_alpha[irrep][ creator + L * ( annihilator + L * cnt_new_alpha ) ] = irrep_old;
                     lookup_sign_alpha [irrep][ creator + L * ( annihilator + L * cnt_new_alpha ) ] = phase;
                     
                     bits[ annihilator ] = 0;
                  } else {
                     phase_annihilator *= -1;
                  }
               }
               
               bits[ creator ] = 1;
               phase_creator *= -1;
            }
         }
      }
      
      lookup_cnt_beta  [ irrep ] = new int[ L * L * numPerIrrep_down[ irrep ] ];
      lookup_irrep_beta[ irrep ] = new int[ L * L * numPerIrrep_down[ irrep ] ];
      lookup_sign_beta [ irrep ] = new int[ L * L * numPerIrrep_down[ irrep ] ];
      
      for ( unsigned int cnt_new_beta = 0; cnt_new_beta < numPerIrrep_down[ irrep ]; cnt_new_beta++ ){
      
         for ( unsigned int i = 0; i < L; i++ ){
            for ( unsigned int j = 0; j < L; j++ ){
               // Check for the sign. If no check for sign, you multiply with sign and everything should be OK...
               lookup_cnt_beta  [irrep][ i + L * ( j + L * cnt_new_beta ) ] = 0;
               lookup_irrep_beta[irrep][ i + L * ( j + L * cnt_new_beta ) ] = 0;
               lookup_sign_beta [irrep][ i + L * ( j + L * cnt_new_beta ) ] = 0;
            }
         }
      
         str2bits( L , cnt2str_down[ irrep ][ cnt_new_beta ] , bits );
      
         int phase_creator = 1;
         for ( unsigned int creator = 0; creator < L; creator++ ){
            if ( bits[ creator ] ){
               bits[ creator ] = 0;
                  
               int phase_annihilator = 1;
               for ( unsigned int annihilator = 0; annihilator < L; annihilator++ ){
                  if ( !(bits[ annihilator ]) ){
                     bits[ annihilator ] = 1;
                     
                     const int irrep_old = getIrrepProduct( irrep , getIrrepProduct( getOrb2Irrep( creator ) , getOrb2Irrep( annihilator ) ) );
                     const int cnt_old = str2cnt_down[ irrep_old ][ bits2str( L , bits ) ];
                     const int phase = phase_creator * phase_annihilator;
                     
                     lookup_cnt_beta  [irrep][ creator + L * ( annihilator + L * cnt_new_beta ) ] = cnt_old;
                     lookup_irrep_beta[irrep][ creator + L * ( annihilator + L * cnt_new_beta ) ] = irrep_old;
                     lookup_sign_beta [irrep][ creator + L * ( annihilator + L * cnt_new_beta ) ] = phase;
                     
                     bits[ annihilator ] = 0;
                  } else {
                     phase_annihilator *= -1;
                  }
               }
               
               bits[ creator ] = 1;
               phase_creator *= -1;
            }
         }
      }
      
   }
   
   delete [] bits; // Delete temporary helper array

}

void CheMPS2::FCI::StartupIrrepCenter(){

   // Find the orbital combinations which can form a center irrep
   irrep_center_num      = new unsigned int [ NumIrreps ];
   irrep_center_crea_orb = new unsigned int*[ NumIrreps ];
   irrep_center_anni_orb = new unsigned int*[ NumIrreps ];
   
   for ( unsigned int irrep_center = 0; irrep_center < NumIrreps; irrep_center++ ){
      const int irrep_center_const_signed = irrep_center;
   
      irrep_center_num[ irrep_center ] = 0;
      for ( unsigned int creator = 0; creator < L; creator++ ){
         for ( unsigned int annihilator = creator; annihilator < L; annihilator++ ){
            if ( getIrrepProduct( getOrb2Irrep( creator ) , getOrb2Irrep( annihilator ) ) == irrep_center_const_signed ){
               irrep_center_num[ irrep_center ] += 1;
            }
         }
      }
      irrep_center_crea_orb[ irrep_center ] = new unsigned int[ irrep_center_num[ irrep_center ] ];
      irrep_center_anni_orb[ irrep_center ] = new unsigned int[ irrep_center_num[ irrep_center ] ];
      irrep_center_num[ irrep_center ] = 0;
      for ( unsigned int creator = 0; creator < L; creator++ ){
         for ( unsigned int annihilator = creator; annihilator < L; annihilator++){
            if ( getIrrepProduct( getOrb2Irrep( creator ) , getOrb2Irrep( annihilator ) ) == irrep_center_const_signed ){
               irrep_center_crea_orb[ irrep_center ][ irrep_center_num[ irrep_center ] ] = creator;
               irrep_center_anni_orb[ irrep_center ][ irrep_center_num[ irrep_center ] ] = annihilator;
               irrep_center_num[ irrep_center ] += 1;
            }
         }
      }
   
   }
   
   irrep_center_jumps = new unsigned long long*[ NumIrreps ];
   HXVsizeWorkspace = 0;
   for ( unsigned int irrep_center = 0; irrep_center < NumIrreps; irrep_center++ ){
   
      irrep_center_jumps[ irrep_center ] = new unsigned long long[ NumIrreps+1 ];
      const int localTargetIrrep = getIrrepProduct( irrep_center , getTargetIrrep() );
      irrep_center_jumps[ irrep_center ][ 0 ] = 0;
      for ( unsigned int irrep_up = 0; irrep_up < NumIrreps; irrep_up++ ){
         const int irrep_down = getIrrepProduct( irrep_up , localTargetIrrep );
         unsigned long long temp  = numPerIrrep_up  [ irrep_up   ];
                            temp *= numPerIrrep_down[ irrep_down ];
         irrep_center_jumps[ irrep_center ][ irrep_up+1 ] = irrep_center_jumps[ irrep_center ][ irrep_up ] + temp;
      }
      if ( irrep_center_num[ irrep_center ] * irrep_center_jumps[ irrep_center ][ NumIrreps ] > HXVsizeWorkspace ){
         HXVsizeWorkspace = irrep_center_num[ irrep_center ] * irrep_center_jumps[ irrep_center ][ NumIrreps ];
      }
   }
   if ( FCIverbose>0 ){
      cout << "FCI::Startup : Number of variables in the FCI vector = " << getVecLength(0) << endl;
      unsigned long long numberOfBytes = 2 * sizeof(double) * HXVsizeWorkspace;
      
      cout << "FCI::Startup : Without additional loops the FCI matrix-vector product requires a workspace of " << 1e-6 * numberOfBytes << " MB memory." << endl;
      if ( maxMemWorkMB < 1e-6 * numberOfBytes ){
         HXVsizeWorkspace = (unsigned long long) ceil( ( maxMemWorkMB * 1e6 ) / ( 2 * sizeof(double ) ) );
         numberOfBytes = 2 * sizeof(double) * HXVsizeWorkspace;
         cout << "               For practical purposes, the workspace is constrained to " << 1e-6 * numberOfBytes << " MB memory." << endl;
      }
   }
   HXVworksmall = new double[ L * L * L * L ];
   HXVworkbig1  = new double[ HXVsizeWorkspace ];
   HXVworkbig2  = new double[ HXVsizeWorkspace ];
   
   // Check for the lapack routines { dgemm_ , daxpy_ , dscal_ , dcopy_ , ddot_ }
   unsigned long long maxVecLength = 0;
   for ( unsigned int irrep = 0; irrep < NumIrreps; irrep++ ){
      if ( getVecLength( irrep ) > maxVecLength ){ maxVecLength = getVecLength( irrep ); }
   }
   const unsigned int max_integer = INT_MAX;
   assert( max_integer >= maxVecLength );
   
}

void CheMPS2::FCI::str2bits(const unsigned int Lval, const unsigned int bitstring, int * bits){

   for (unsigned int bit = 0; bit < Lval; bit++){ bits[ bit ] = ( bitstring & ( 1 << bit ) ) >> bit; }

}

unsigned int CheMPS2::FCI::bits2str(const unsigned int Lval, int * bits){

   unsigned int factor = 1;
   unsigned int result = 0;
   for (unsigned int bit = 0; bit < Lval; bit++){
      result += bits[ bit ] * factor;
      factor *= 2;
   }
   return result;

}

int CheMPS2::FCI::getUpIrrepOfCounter(const int irrep_center, const unsigned long long counter) const{

   int irrep_up = NumIrreps;
   while ( counter < irrep_center_jumps[ irrep_center ][ irrep_up-1 ] ){ irrep_up--; }
   return irrep_up-1;
   
}

void CheMPS2::FCI::getBitsOfCounter(const int irrep_center, const unsigned long long counter, int * bits_up, int * bits_down) const{

   const int localTargetIrrep = getIrrepProduct( irrep_center , TargetIrrep );
   
   const int irrep_up   = getUpIrrepOfCounter( irrep_center , counter );
   const int irrep_down = getIrrepProduct( irrep_up , localTargetIrrep );
   
   const unsigned int count_up   = ( counter - irrep_center_jumps[ irrep_center ][ irrep_up ] ) % numPerIrrep_up[ irrep_up ];
   const unsigned int count_down = ( counter - irrep_center_jumps[ irrep_center ][ irrep_up ] ) / numPerIrrep_up[ irrep_up ];
   
   const unsigned int string_up   = cnt2str_up  [ irrep_up   ][ count_up   ];
   const unsigned int string_down = cnt2str_down[ irrep_down ][ count_down ];
   
   str2bits( L , string_up   , bits_up   );
   str2bits( L , string_down , bits_down );

}

double CheMPS2::FCI::getFCIcoeff(int * bits_up, int * bits_down, double * vector) const{

   const unsigned string_up   = bits2str(L, bits_up  );
   const unsigned string_down = bits2str(L, bits_down);
   
   int irrep_up   = 0;
   int irrep_down = 0;
   for ( unsigned int orb = 0; orb < L; orb++ ){
      if ( bits_up  [ orb ] ){ irrep_up   = getIrrepProduct( irrep_up   , getOrb2Irrep( orb ) ); }
      if ( bits_down[ orb ] ){ irrep_down = getIrrepProduct( irrep_down , getOrb2Irrep( orb ) ); }
   }
   
   const int counter_up   = str2cnt_up  [ irrep_up   ][ string_up   ];
   const int counter_down = str2cnt_down[ irrep_down ][ string_down ];
   
   if (( counter_up == -1 ) || ( counter_down == -1 )){ return 0.0; }
   
   return vector[ irrep_center_jumps[ 0 ][ irrep_up ] + counter_up + numPerIrrep_up[ irrep_up ] * counter_down ];

}

/*void CheMPS2::FCI::CheckHamDEBUG() const{

   const unsigned long long vecLength = getVecLength( 0 );
   
   // Building Ham by HamTimesVec
   double * HamHXV = new double[ vecLength * vecLength ];
   double * workspace = new double[ vecLength ];
   for (unsigned long long count = 0; count < vecLength; count++){
   
      ClearVector( vecLength , workspace );
      workspace[ count ] = 1.0;
      HamTimesVec( workspace , HamHXV + count*vecLength );
   
   }
   
   // Building Diag by HamDiag
   DiagHam( workspace );
   double RMSdiagdifference = 0.0;
   for (unsigned long long row = 0; row < vecLength; row++){
      double diff = workspace[ row ] - HamHXV[ row + vecLength * row ];
      RMSdiagdifference += diff * diff;
   }
   RMSdiagdifference = sqrt( RMSdiagdifference );
   cout << "The RMS difference of DiagHam() and diag(HamHXV) = " << RMSdiagdifference << endl;
   
   // Building Ham by getMatrixElement
   int * work     = new int[ 8 ];
   int * ket_up   = new int[ L ];
   int * ket_down = new int[ L ];
   int * bra_up   = new int[ L ];
   int * bra_down = new int[ L ];
   double RMSconstructiondiff = 0.0;
   for (unsigned long long row = 0; row < vecLength; row++){
      for (unsigned long long col = 0; col < vecLength; col++){
         getBitsOfCounter( 0 , row , bra_up , bra_down );
         getBitsOfCounter( 0 , col , ket_up , ket_down );
         double tempvar = HamHXV[ row + vecLength * col ] - GetMatrixElement( bra_up , bra_down , ket_up , ket_down , work );
         RMSconstructiondiff += tempvar * tempvar;
      }
   }
   cout << "The RMS difference of HamHXV - HamMXELEM = " << RMSconstructiondiff << endl;
   delete [] work;
   delete [] ket_up;
   delete [] ket_down;
   delete [] bra_up;
   delete [] bra_down;
   
   // Building Ham^2 by HamTimesVec
   double * workspace2 = new double[ vecLength ];
   for (unsigned long long count = 0; count < vecLength; count++){
   
      ClearVector( vecLength , workspace );
      workspace[ count ] = 1.0;
      HamTimesVec( workspace , workspace2 );
      HamTimesVec( workspace2 , HamHXV + count*vecLength );
   
   }
   
   // Building diag( Ham^2 ) by DiagHamSquared
   DiagHamSquared( workspace );
   double RMSdiagdifference2 = 0.0;
   for (unsigned long long row = 0; row < vecLength; row++){
      double diff = workspace[ row ] - HamHXV[ row + vecLength * row ];
      RMSdiagdifference2 += diff * diff;
   }
   RMSdiagdifference2 = sqrt( RMSdiagdifference2 );
   cout << "The RMS difference of DiagHamSquared() and diag(HamSquared by HXV) = " << RMSdiagdifference2 << endl;
   
   delete [] workspace2;
   delete [] workspace;
   delete [] HamHXV;

}*/

void CheMPS2::FCI::HamTimesVec(double * input, double * output) const{

   ClearVector( getVecLength( 0 ) , output );

   // P.J. Knowles and N.C. Handy, A new determinant-based full configuration interaction method, Chemical Physics Letters 111 (4-5), 315-321 (1984)
   
   // irrep_center is the center irrep of the ERI : (ij|kl) --> irrep_center = I_i x I_j = I_k x I_l
   for ( unsigned int irrep_center = 0; irrep_center < NumIrreps; irrep_center++ ){
   
      const unsigned long long localVecLength = getVecLength( irrep_center );
      const int localTargetIrrep = getIrrepProduct( TargetIrrep , irrep_center );
      const unsigned int numPairs = irrep_center_num[ irrep_center ];

      const unsigned int space_per_vectorpiece = (unsigned int) floor(( 1.0 * HXVsizeWorkspace ) / numPairs);
      unsigned int numIterations = localVecLength / space_per_vectorpiece;
      if ( localVecLength > numIterations * space_per_vectorpiece ){ numIterations++; }

      for ( unsigned int iteration = 0; iteration < numIterations; iteration++ ){
      
         const unsigned long long veccounter_start = iteration * space_per_vectorpiece;
         const unsigned long long guess_stop       = ( iteration + 1 ) * space_per_vectorpiece;
         const unsigned long long veccounter_stop  = ( guess_stop > localVecLength ) ? localVecLength : guess_stop;
   
         /******************************************************************************************************************************
          *   First build workbig1[ i<=j + size(i<=j) * veccounter ] = E_{i<=j} + ( 1 - delta_i==j ) E_{j>i} (irrep_center) | input >  *
          ******************************************************************************************************************************/
         const unsigned long long loopsize = numPairs * ( veccounter_stop - veccounter_start );
         #pragma omp parallel for schedule(static)
         for ( unsigned long long loopvariable = 0; loopvariable < loopsize; loopvariable++ ){
         
            const unsigned int pair               = loopvariable % numPairs;
            const unsigned long long veccounter   = veccounter_start + ( loopvariable / numPairs );
            const unsigned int creator            = irrep_center_crea_orb[ irrep_center ][ pair ];
            const unsigned int annihilator        = irrep_center_anni_orb[ irrep_center ][ pair ];
            const int irrep_new_up                = getUpIrrepOfCounter( irrep_center , veccounter );
            const int irrep_new_down              = getIrrepProduct( irrep_new_up , localTargetIrrep );
            const unsigned int count_new_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_new_up ] ) % numPerIrrep_up[ irrep_new_up ];
            const unsigned int count_new_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_new_up ] ) / numPerIrrep_up[ irrep_new_up ];
            
            double myResult = 0.0;
            
            {
               // E^{alpha}_{creator <= annihilator}
               const int entry_up = creator + L * ( annihilator + L * count_new_up );
               const int sign_up  = lookup_sign_alpha [ irrep_new_up ][ entry_up ];
               if ( sign_up != 0 ){ // Required for one-electron calculations
                  const int irrep_old_up = lookup_irrep_alpha[ irrep_new_up ][ entry_up ];
                  const int cnt_old_up   = lookup_cnt_alpha  [ irrep_new_up ][ entry_up ];
                  myResult = sign_up * input[ irrep_center_jumps[ 0 ][ irrep_old_up ] + cnt_old_up + numPerIrrep_up[ irrep_old_up ] * count_new_down ];
               }
               
               // E^{beta}_{creator <= annihilator}
               const int entry_down = creator + L * ( annihilator + L * count_new_down );
               const int sign_down  = lookup_sign_beta[ irrep_new_down ][ entry_down ];
               if ( sign_down != 0 ){ // Required for one-electron calculations
                  const int cnt_old_down = lookup_cnt_beta[ irrep_new_down ][ entry_down ];
                  myResult += sign_down * input[ irrep_center_jumps[ 0 ][ irrep_new_up ] + count_new_up + numPerIrrep_up[ irrep_new_up ] * cnt_old_down ];
               }
            }
            
            if ( annihilator > creator ){
               // E^{alpha}_{annihilator > creator}
               const int entry_up = annihilator + L * ( creator + L * count_new_up );
               const int sign_up  = lookup_sign_alpha [ irrep_new_up ][ entry_up ];
               if ( sign_up != 0 ){ // Required for one-electron calculations
                  const int irrep_old_up = lookup_irrep_alpha[ irrep_new_up ][ entry_up ];
                  const int cnt_old_up   = lookup_cnt_alpha  [ irrep_new_up ][ entry_up ];
                  myResult += sign_up * input[ irrep_center_jumps[ 0 ][ irrep_old_up ] + cnt_old_up + numPerIrrep_up[ irrep_old_up ] * count_new_down ];
               }
               
               // E^{beta}_{annihilator > creator}
               const int entry_down = annihilator + L * ( creator + L * count_new_down );
               const int sign_down  = lookup_sign_beta[ irrep_new_down ][ entry_down ];
               if ( sign_down != 0 ){ // Required for one-electron calculations
                  const int cnt_old_down = lookup_cnt_beta[ irrep_new_down ][ entry_down ];
                  myResult += sign_down * input[ irrep_center_jumps[ 0 ][ irrep_new_up ] + count_new_up + numPerIrrep_up[ irrep_new_up ] * cnt_old_down ];
               }
            }
            
            HXVworkbig1[ loopvariable ] = myResult;
            
         }
         
         /************************************************
          *   If irrep_center==0, do the one-body terms  *
          ************************************************/
         if ( irrep_center==0 ){
            for ( unsigned int pair = 0; pair < numPairs; pair++ ){
               HXVworksmall[ pair ] = getGmat( irrep_center_crea_orb[ irrep_center ][ pair ] , irrep_center_anni_orb[ irrep_center ][ pair ] );
            }
            char trans  = 'T';
            char notran = 'N';
            double one  = 1.0;
            int mdim = veccounter_stop - veccounter_start; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
            int kdim = numPairs;
            int ndim = 1;
            dgemm_( &trans, &notran, &mdim, &ndim, &kdim, &one, HXVworkbig1, &kdim, HXVworksmall, &kdim, &one, output + veccounter_start, &mdim );
         }
         
         /****************************************************************************************************************************
          *   Now build workbig2[ i<=j + size(i<=j) * veccounter] = 0.5 * ( i<=j | k<=l ) * workbig1[ k<=l + size(k<=l) * counter ]  *
          ****************************************************************************************************************************/
         {
            for ( unsigned int pair1 = 0; pair1 < numPairs; pair1++ ){
               for ( unsigned int pair2 = 0; pair2 < numPairs; pair2++ ){
                  HXVworksmall[ pair1 + numPairs * pair2 ]
                     = 0.5 * getERI( irrep_center_crea_orb[ irrep_center ][ pair1 ] , irrep_center_anni_orb[ irrep_center ][ pair1 ] ,
                                     irrep_center_crea_orb[ irrep_center ][ pair2 ] , irrep_center_anni_orb[ irrep_center ][ pair2 ] );
               }
            }
            char notran = 'N';
            double one  = 1.0;
            double zero = 0.0;
            int mdim = numPairs;
            int kdim = numPairs;
            int ndim = veccounter_stop - veccounter_start; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
            dgemm_( &notran , &notran , &mdim , &ndim , &kdim , &one , HXVworksmall , &mdim , HXVworkbig1 , &kdim , &zero , HXVworkbig2 , &mdim );
         }
         
         /*************************************************************************************************************
          *   Finally do output <-- E_{i<=j} + (1 - delta_{i==j}) E_{j>i} workbig2[ i<=j + size(i<=j) * veccounter ]  *
          *************************************************************************************************************/
         for ( unsigned int pair = 0; pair < numPairs; pair++ ){
         
            const unsigned int orbi = irrep_center_crea_orb[ irrep_center ][ pair ];
            const unsigned int orbj = irrep_center_anni_orb[ irrep_center ][ pair ];

            #pragma omp parallel for schedule(static) // The given E_{i<=j}^{alpha} connects exactly one veccounter with one location_new
            for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
               const int irrep_old_up            = getUpIrrepOfCounter( irrep_center , veccounter );
               const unsigned int count_old_up   = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
               const int entry_up                = orbj + L * ( orbi + L * count_old_up );
               const int sign_up                 = lookup_sign_alpha[ irrep_old_up ][ entry_up ];
               if ( sign_up != 0 ){ // Required for thread safety
                  const unsigned int count_old_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
                  const int irrep_new_up                = lookup_irrep_alpha[ irrep_old_up ][ entry_up ];
                  const int cnt_new_up                  = lookup_cnt_alpha  [ irrep_old_up ][ entry_up ];
                  const unsigned long long location_new = irrep_center_jumps[ 0 ][ irrep_new_up ] + cnt_new_up + numPerIrrep_up[ irrep_new_up ] * count_old_down;
                  output[ location_new ] += sign_up * HXVworkbig2[ pair + numPairs * ( veccounter - veccounter_start ) ];
               }
            }
            
            #pragma omp parallel for schedule(static) // The given E_{i<=j}^{beta} connects exactly one veccounter with one location_new
            for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
               const int irrep_old_up            = getUpIrrepOfCounter( irrep_center , veccounter );
               const unsigned int count_old_down = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
               const int entry_down              = orbj + L * ( orbi + L * count_old_down );
               const int irrep_old_down          = getIrrepProduct( irrep_old_up , localTargetIrrep );
               const int sign_down               = lookup_sign_beta[ irrep_old_down ][ entry_down ];
               if ( sign_down != 0 ){ // Required for thread safety
                  const unsigned int count_old_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
                  const int cnt_new_down                = lookup_cnt_beta[ irrep_old_down ][ entry_down ];
                  const unsigned long long location_new = irrep_center_jumps[ 0 ][ irrep_old_up ] + count_old_up + numPerIrrep_up[ irrep_old_up ] * cnt_new_down;
                  output[ location_new ] += sign_down * HXVworkbig2[ pair + numPairs * ( veccounter - veccounter_start ) ];
               }
            }
            
            if ( orbj > orbi ){
            
               #pragma omp parallel for schedule(static) // The given E_{j>i}^{alpha} connects exactly one veccounter with one location_new
               for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
                  const int irrep_old_up            = getUpIrrepOfCounter( irrep_center , veccounter );
                  const unsigned int count_old_up   = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
                  const int entry_up                = orbi + L * ( orbj + L * count_old_up );
                  const int sign_up                 = lookup_sign_alpha[ irrep_old_up ][ entry_up ];
                  if ( sign_up != 0 ){ // Required for thread safety
                     const unsigned int count_old_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
                     const int irrep_new_up                = lookup_irrep_alpha[ irrep_old_up ][ entry_up ];
                     const int cnt_new_up                  = lookup_cnt_alpha  [ irrep_old_up ][ entry_up ];
                     const unsigned long long location_new = irrep_center_jumps[ 0 ][ irrep_new_up ] + cnt_new_up + numPerIrrep_up[ irrep_new_up ] * count_old_down;
                     output[ location_new ] += sign_up * HXVworkbig2[ pair + numPairs * ( veccounter - veccounter_start ) ];
                  }
               }
               
               #pragma omp parallel for schedule(static) // The given E_{j>i}^{beta} connects exactly one veccounter with one location_new
               for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
                  const int irrep_old_up            = getUpIrrepOfCounter( irrep_center , veccounter );
                  const unsigned int count_old_down = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
                  const int entry_down              = orbi + L * ( orbj + L * count_old_down );
                  const int irrep_old_down          = getIrrepProduct( irrep_old_up , localTargetIrrep );
                  const int sign_down               = lookup_sign_beta[ irrep_old_down ][ entry_down ];
                  if ( sign_down != 0 ){ // Required for thread safety
                     const unsigned int count_old_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
                     const int cnt_new_down                = lookup_cnt_beta[ irrep_old_down ][ entry_down ];
                     const unsigned long long location_new = irrep_center_jumps[ 0 ][ irrep_old_up ] + count_old_up + numPerIrrep_up[ irrep_old_up ] * cnt_new_down;
                     output[ location_new ] += sign_down * HXVworkbig2[ pair + numPairs * ( veccounter - veccounter_start ) ];
                  }
               }
            }
            
         }
      }
   }

}

double CheMPS2::FCI::Fill2RDM(double * vector, double * TwoRDM) const{

   ClearVector( L*L*L*L , TwoRDM );
   assert( Nel_up + Nel_down >= 2 );

   // Based HamTimesVec: Gamma^2(i,j,k,l) = sum_sigma,tau < a^+_i,sigma a^+j,tau a_l,tau a_l,sigma > = TwoRDM[ i + L * ( j + L * ( k + L * l ) ) ]
   
   // irrep_center is the center irrep of the ERI : <ij|kl> --> irrep_center = I_i x I_j = I_k x I_l
   for ( unsigned int irrep_center = 0; irrep_center < NumIrreps; irrep_center++ ){
   
      const unsigned long long localVecLength = getVecLength( irrep_center );
      const int localTargetIrrep = getIrrepProduct( TargetIrrep , irrep_center );
      const unsigned int numPairs = irrep_center_num[ irrep_center ];
      
      const unsigned int space_per_vectorpiece = (unsigned int) floor(( 1.0 * HXVsizeWorkspace ) / numPairs);
      unsigned int numIterations = localVecLength / space_per_vectorpiece;
      if ( localVecLength > numIterations * space_per_vectorpiece ){ numIterations++; }

      for ( unsigned int iteration = 0; iteration < numIterations; iteration++ ){
      
         const unsigned long long veccounter_start = iteration * space_per_vectorpiece;
         const unsigned long long guess_stop       = ( iteration + 1 ) * space_per_vectorpiece;
         const unsigned long long veccounter_stop  = ( guess_stop > localVecLength ) ? localVecLength : guess_stop;
   
         /*************************************************************************************************
          *   First build workbig1[ i<=j + size(i<=j) * veccounter ] = E_{i<=j} (irrep_center) | input >  *
          *************************************************************************************************/
         const unsigned long long loopsize = numPairs * ( veccounter_stop - veccounter_start );
         #pragma omp parallel for schedule(static)
         for ( unsigned long long loopvariable = 0; loopvariable < loopsize; loopvariable++ ){
         
            const unsigned int pair               = loopvariable % numPairs;
            const unsigned long long veccounter   = veccounter_start + ( loopvariable / numPairs );
            const unsigned int creator            = irrep_center_crea_orb[ irrep_center ][ pair ];
            const unsigned int annihilator        = irrep_center_anni_orb[ irrep_center ][ pair ];
            const int irrep_new_up                = getUpIrrepOfCounter( irrep_center , veccounter );
            const int irrep_new_down              = getIrrepProduct( irrep_new_up , localTargetIrrep );
            const unsigned int count_new_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_new_up ] ) % numPerIrrep_up[ irrep_new_up ];
            const unsigned int count_new_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_new_up ] ) / numPerIrrep_up[ irrep_new_up ];
            
            double myResult = 0.0;

            // E^{alpha}_{creator <= annihilator}
            const int entry_up = creator + L * ( annihilator + L * count_new_up );
            const int sign_up  = lookup_sign_alpha[ irrep_new_up ][ entry_up ];
            if ( sign_up != 0 ){ // Required for one-electron calculations
               const int irrep_old_up = lookup_irrep_alpha[ irrep_new_up ][ entry_up ];
               const int cnt_old_up   = lookup_cnt_alpha  [ irrep_new_up ][ entry_up ];
               myResult = sign_up * vector[ irrep_center_jumps[ 0 ][ irrep_old_up ] + cnt_old_up + numPerIrrep_up[ irrep_old_up ] * count_new_down ];
            }
            
            // E^{beta}_{creator <= annihilator}
            const int entry_down = creator + L * ( annihilator + L * count_new_down );
            const int sign_down  = lookup_sign_beta[ irrep_new_down ][ entry_down ];
            if ( sign_down != 0 ){ // Required for one-electron calculations
               const int cnt_old_down = lookup_cnt_beta[ irrep_new_down ][ entry_down ];
               myResult += sign_down * vector[ irrep_center_jumps[ 0 ][ irrep_new_up ] + count_new_up + numPerIrrep_up[ irrep_new_up ] * cnt_old_down ];
            }
            
            HXVworkbig1[ loopvariable ] = myResult;
            
         }
         
         /************************************************
          *   If irrep_center==0, do the one-body terms  *
          ************************************************/
         if ( irrep_center==0 ){
            
            int mdim = numPairs;
            int kdim = veccounter_stop - veccounter_start; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
            int ndim = 1;
            
            char notran = 'N';
            double one  = 1.0;
            double zero = 0.0;
            dgemm_( &notran, &notran, &mdim, &ndim, &kdim, &one, HXVworkbig1, &mdim, vector + veccounter_start, &kdim, &zero, HXVworksmall, &mdim );
            
            for ( unsigned int pair = 0; pair < numPairs; pair++ ){
               const unsigned int orb1 = irrep_center_crea_orb[ irrep_center ][ pair ];
               const unsigned int orb2 = irrep_center_anni_orb[ irrep_center ][ pair ];
               for ( unsigned int orbk = 0; orbk < L; orbk++ ){
                  TwoRDM[ orb1 + L * ( orbk + L * ( orbk + L * orb2 ) ) ] -= HXVworksmall[ pair ];
               }
               if ( orb2 > orb1 ){
                  for ( unsigned int orbk = 0; orbk < L; orbk++ ){
                     TwoRDM[ orb2 + L * ( orbk + L * ( orbk + L * orb1 ) ) ] -= HXVworksmall[ pair ];
                  }
               }
            }
            
         }
         
         /*********************************************************************
          *   Finally do E_{i<=k} workbig1[ j<=l + size(j<=l) * veccounter ]  *
          *              E_{k> i} workbig1[ j<=l + size(j<=l) * veccounter ]  *
          *********************************************************************/
         const unsigned int numPairsSquared = numPairs * numPairs;
         #pragma omp parallel for schedule(static) //Thread safe since each thread can only write to its own location in TwoRDM
         for ( unsigned int pairproduct = 0; pairproduct < numPairsSquared; pairproduct++ ){
         
            const unsigned int pair1 = pairproduct % numPairs;
            const unsigned int pair2 = pairproduct / numPairs;
            const unsigned int orb_i = irrep_center_crea_orb[ irrep_center ][ pair1 ];
            const unsigned int orb_k = irrep_center_anni_orb[ irrep_center ][ pair1 ];
            const unsigned int orb_j = irrep_center_crea_orb[ irrep_center ][ pair2 ];
            const unsigned int orb_l = irrep_center_anni_orb[ irrep_center ][ pair2 ];
            // orb_i <= orb_k and orb_j <= orb_l
            
            double myResult = 0.0;
            
            // E_{i<=k} E_{j<=l} vector
            for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
            
               const unsigned long long location_old = pair2 + numPairs * ( veccounter - veccounter_start );
               const int irrep_old_up                = getUpIrrepOfCounter( irrep_center , veccounter );
               const int irrep_old_down              = getIrrepProduct( irrep_old_up , localTargetIrrep );
               const unsigned int count_old_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
               const unsigned int count_old_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
               
               // E^{alpha}_{i<=k}
               const int entry_up = orb_k + L * ( orb_i + L * count_old_up );
               const int sign_up  = lookup_sign_alpha [ irrep_old_up ][ entry_up ];
               if ( sign_up != 0 ){
                  const int irrep_new_up = lookup_irrep_alpha[ irrep_old_up ][ entry_up ];
                  const int cnt_new_up   = lookup_cnt_alpha  [ irrep_old_up ][ entry_up ];
                  myResult += sign_up * HXVworkbig1[ location_old ] * vector[ irrep_center_jumps[ 0 ][ irrep_new_up ] + cnt_new_up + numPerIrrep_up[ irrep_new_up ] * count_old_down ];
               }
               
               // E^{beta}_{i<=k}
               const int entry_down = orb_k + L * ( orb_i + L * count_old_down );
               const int sign_down  = lookup_sign_beta[ irrep_old_down ][ entry_down ];
               if ( sign_down != 0 ){
                  const int cnt_new_down = lookup_cnt_beta[ irrep_old_down ][ entry_down ];
                  myResult += sign_down * HXVworkbig1[ location_old ] * vector[ irrep_center_jumps[ 0 ][ irrep_old_up ] + count_old_up + numPerIrrep_up[ irrep_old_up ] * cnt_new_down ];
               }
            }
            
            TwoRDM[ orb_i + L * ( orb_j + L * ( orb_k + L * orb_l ) ) ] += myResult; // i<=k and j<=l
            
            if ( orb_k > orb_i ){
            
               myResult = 0.0;
               
               // E_{k>i} E_{j<=l} vector
               for ( unsigned long long veccounter = veccounter_start; veccounter < veccounter_stop; veccounter++ ){
               
                  const unsigned long long location_old = pair2 + numPairs * ( veccounter - veccounter_start );
                  const int irrep_old_up                = getUpIrrepOfCounter( irrep_center , veccounter );
                  const int irrep_old_down              = getIrrepProduct( irrep_old_up , localTargetIrrep );
                  const unsigned int count_old_up       = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) % numPerIrrep_up[ irrep_old_up ];
                  const unsigned int count_old_down     = ( veccounter - irrep_center_jumps[ irrep_center ][ irrep_old_up ] ) / numPerIrrep_up[ irrep_old_up ];
                  
                  // E^{alpha}_{k>i}
                  const int entry_up = orb_i + L * ( orb_k + L * count_old_up );
                  const int sign_up  = lookup_sign_alpha [ irrep_old_up ][ entry_up ];
                  if ( sign_up != 0 ){
                     const int irrep_new_up = lookup_irrep_alpha[ irrep_old_up ][ entry_up ];
                     const int cnt_new_up   = lookup_cnt_alpha  [ irrep_old_up ][ entry_up ];
                     myResult += sign_up * HXVworkbig1[ location_old ] * vector[ irrep_center_jumps[ 0 ][ irrep_new_up ] + cnt_new_up + numPerIrrep_up[ irrep_new_up ] * count_old_down ];
                  }
                  
                  // E^{beta}_{k>i}
                  const int entry_down = orb_i + L * ( orb_k + L * count_old_down );
                  const int sign_down  = lookup_sign_beta[ irrep_old_down ][ entry_down ];
                  if ( sign_down != 0 ){
                     const int cnt_new_down = lookup_cnt_beta[ irrep_old_down ][ entry_down ];
                     myResult += sign_down * HXVworkbig1[ location_old ] * vector[ irrep_center_jumps[ 0 ][ irrep_old_up ] + count_old_up + numPerIrrep_up[ irrep_old_up ] * cnt_new_down ];
                  }
               }
            
               TwoRDM[ orb_k + L * ( orb_j + L * ( orb_i + L * orb_l ) ) ] += myResult; // "i>k" and j<=l
               
            }
         }
      }
   }
   
   // For Gamma^2_{ijkl} the elements with i,k anything and j<=l are set: fill the rest
   for ( unsigned int orb1 = 0; orb1 < L; orb1++ ){
      for ( unsigned int orb3 = 0; orb3 < L; orb3++ ){
         for ( unsigned int orb2 = 0; orb2 < L-1; orb2++ ){
            for ( unsigned int orb4 = orb2+1; orb4 < L; orb4++ ){
               TwoRDM[ orb3 + L * ( orb4 + L * ( orb1 + L * orb2 ) ) ] = TwoRDM[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ];
            }
         }
      }
   }

   // Calculate the FCI energy
   double FCIenergy = getEconst();
   for ( unsigned int orb1 = 0; orb1 < L; orb1++ ){
      for ( unsigned int orb2 = 0; orb2 < L; orb2++ ){
         double tempvar = 0.0;
         double tempvar2 = 0.0;
         for ( unsigned int orb3 = 0; orb3 < L; orb3++ ){
            tempvar  += getERI( orb1 , orb3 , orb3 , orb2 );
            tempvar2 += TwoRDM[ orb1 + L * ( orb3 + L * ( orb2 + L * orb3 ) ) ];
            for ( unsigned int orb4 = 0; orb4 < L; orb4++ ){
               FCIenergy += 0.5 * TwoRDM[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ] * getERI( orb1 , orb3 , orb2 , orb4 );
            }
         }
         FCIenergy += ( getGmat( orb1 , orb2 ) + 0.5 * tempvar ) * tempvar2 / ( Nel_up + Nel_down - 1.0);
      }
   }
   if ( FCIverbose > 0 ){ cout << "FCI::Fill2RDM : Energy based on 2-RDM = " << FCIenergy << endl; }
   return FCIenergy;

}

double CheMPS2::FCI::CalcSpinSquared(double * vector) const{

   const unsigned long long vecLength = getVecLength( 0 );
   double result = 0.0;
      
   #pragma omp parallel for schedule(static) reduction(+:result)
   for ( unsigned long long counter = 0; counter < vecLength; counter++ ){
      for ( unsigned int orbi = 0; orbi < L; orbi++ ){
         
         const int irrep_up     = getUpIrrepOfCounter( 0 , counter );
         const int irrep_down   = getIrrepProduct( irrep_up , TargetIrrep );
         const int count_up     = ( counter - irrep_center_jumps[ 0 ][ irrep_up ] ) % numPerIrrep_up[ irrep_up ];
         const int count_down   = ( counter - irrep_center_jumps[ 0 ][ irrep_up ] ) / numPerIrrep_up[ irrep_up ];
         
         // Diagonal terms
         const int diff_ii = lookup_sign_alpha[ irrep_up   ][ orbi + L * ( orbi + L * count_up   ) ]
                           - lookup_sign_beta [ irrep_down ][ orbi + L * ( orbi + L * count_down ) ]; //Signed integers so subtracting is OK
         const double vector_at_counter_squared = vector[ counter ] * vector[ counter ];
         result += 0.75 * diff_ii * diff_ii * vector_at_counter_squared;
         
         for ( unsigned int orbj = orbi+1; orbj < L; orbj++ ){
         
            // Sz Sz
            const int diff_jj = lookup_sign_alpha[ irrep_up   ][ orbj + L * ( orbj + L * count_up   ) ]
                              - lookup_sign_beta [ irrep_down ][ orbj + L * ( orbj + L * count_down ) ]; //Signed integers so subtracting is OK
            result += 0.5 * diff_ii * diff_jj * vector_at_counter_squared;
            
            const int irrep_up_bis = getIrrepProduct( irrep_up , getIrrepProduct( getOrb2Irrep( orbi ) , getOrb2Irrep( orbj ) ) );
            
            // - ( a_i,up^+ a_j,up )( a_j,down^+ a_i,down )
            const int entry_down_ji = orbj + L * ( orbi + L * count_down );
            const int sign_down_ji  = lookup_sign_beta [ irrep_down ][ entry_down_ji ];
            const int entry_up_ij   = orbi + L * ( orbj + L * count_up );
            const int sign_up_ij    = lookup_sign_alpha[ irrep_up ][ entry_up_ij ];
            const int sign_product1 = sign_up_ij * sign_down_ji;
            if ( sign_product1 != 0 ){
               const int cnt_down_ji = lookup_cnt_beta[ irrep_down ][ entry_down_ji ];
               const int cnt_up_ij   = lookup_cnt_alpha[ irrep_up ][ entry_up_ij ];
               result -= sign_product1 * vector[ irrep_center_jumps[ 0 ][ irrep_up_bis ] + cnt_up_ij + numPerIrrep_up[ irrep_up_bis ] * cnt_down_ji ] * vector[ counter ];
            }

            // - ( a_j,up^+ a_i,up )( a_i,down^+ a_j,down )
            const int entry_down_ij = orbi + L * ( orbj + L * count_down );
            const int sign_down_ij  = lookup_sign_beta[ irrep_down ][ entry_down_ij ];
            const int entry_up_ji   = orbj + L * ( orbi + L * count_up );
            const int sign_up_ji    = lookup_sign_alpha[ irrep_up ][ entry_up_ji ];
            const int sign_product2 = sign_up_ji * sign_down_ij;
            if ( sign_product2 != 0 ){
               const int cnt_down_ij = lookup_cnt_beta[ irrep_down ][ entry_down_ij ];
               const int cnt_up_ji   = lookup_cnt_alpha[ irrep_up ][ entry_up_ji ];
               result -= sign_product2 * vector[ irrep_center_jumps[ 0 ][ irrep_up_bis ] + cnt_up_ji + numPerIrrep_up[ irrep_up_bis ] * cnt_down_ij ] * vector[ counter ];
            }
         
         }
      }
   }
   
   if ( FCIverbose > 0 ){
      const double intendedS = fabs( 0.5 * Nel_up - 0.5 * Nel_down ); // Be careful with subtracting unsigned integers...
      cout << "FCI::CalcSpinSquared : For intended spin " << intendedS
           << " the measured S(S+1) = " << result << " and intended S(S+1) = " << intendedS * (intendedS + 1.0) << endl;
   }
   return result;

}

void CheMPS2::FCI::DiagHam(double * diag) const{

   const unsigned long long vecLength = getVecLength( 0 );

   #pragma omp parallel
   {

      int * bits_up   = new int[ L ];
      int * bits_down = new int[ L ];
      
      #pragma omp for schedule(static)
      for ( unsigned long long counter = 0; counter < vecLength; counter++ ){
      
         double myResult = 0.0;
         getBitsOfCounter( 0 , counter , bits_up , bits_down ); // Fetch the corresponding bits
         
         for ( unsigned int orb1 = 0; orb1 < L; orb1++ ){
            const int n_tot_orb1 = bits_up[ orb1 ] + bits_down[ orb1 ];
            myResult += n_tot_orb1 * getGmat( orb1 , orb1 );
            for ( unsigned int orb2 = 0; orb2 < L; orb2++ ){
               myResult += 0.5 * n_tot_orb1 * ( bits_up[ orb2 ] + bits_down[ orb2 ] ) * getERI( orb1 , orb1 , orb2 , orb2 );
               myResult += 0.5 * ( n_tot_orb1 - bits_up[ orb1 ] * bits_up[ orb2 ] - bits_down[ orb1 ] * bits_down[ orb2 ] ) * getERI( orb1 , orb2 , orb2 , orb1 );
            }
         }
         
         diag[ counter ] = myResult;
         
      }
      
      delete [] bits_up;
      delete [] bits_down;
      
   }

}

void CheMPS2::FCI::DiagHamSquared(double * output) const{

   const unsigned long long vecLength = getVecLength( 0 );
   
   /*
      Completeness relation strategy to calculate the diagonal elements of Ham*Ham:
      Loop over the FCI determinants |det>
         < det | Ham * Ham | det > = 0.0
         Loop over the FCI determinants which are connected with |det> by Ham |det2>
            < det | Ham * Ham | det > += ( < det2 | Ham | det > )^2
   */
   
   #pragma omp parallel
   {

      int * bits_up_ket    = new int[ L ];
      int * bits_down_ket  = new int[ L ];
      int * bits_up_bra    = new int[ L ];
      int * bits_down_bra  = new int[ L ];
      int * work = new int[ 8 ];
      bool * workspace = new bool[ vecLength ];
      
      #pragma omp for schedule(static)
      for (unsigned long long outcounter = 0; outcounter < vecLength; outcounter++){
      
         double myResult = 0.0;
         
         //Workspace keeps track that a det2 is used only once for the completeness relation
         for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ workspace[ cnt ] = true; }
      
         //Find the irrep_up and irrep_down
         const int irrep_out_up   = getUpIrrepOfCounter( 0 , outcounter );
         const int irrep_out_down = getIrrepProduct( irrep_out_up , TargetIrrep );
         
         //Find the counters for the alpha and the beta electrons
         const unsigned int outcount_up   = ( outcounter - irrep_center_jumps[ 0 ][ irrep_out_up ] ) % numPerIrrep_up[ irrep_out_up ];
         const unsigned int outcount_down = ( outcounter - irrep_center_jumps[ 0 ][ irrep_out_up ] ) / numPerIrrep_up[ irrep_out_up ];
         const unsigned int string_out_up   = cnt2str_up  [ irrep_out_up   ][ outcount_up   ];
         const unsigned int string_out_down = cnt2str_down[ irrep_out_down ][ outcount_down ];
         str2bits(L, string_out_up,   bits_up_ket  );
         str2bits(L, string_out_down, bits_down_ket);
         
         for (unsigned int cnt = 0; cnt < L; cnt++){
            bits_up_bra  [ cnt ] = bits_up_ket  [ cnt ];
            bits_down_bra[ cnt ] = bits_down_ket[ cnt ];
         }
         
         // If only 1 electron: Acting with a^+_i,up a_k,up
         if (( Nel_up == 1 ) && ( Nel_down == 0 )){
            const unsigned long long incounter_offset = irrep_center_jumps[ 0 ][ irrep_out_up ] + numPerIrrep_up[ irrep_out_up ] * outcount_down;
            
            //Loop i to remove from the vector
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_up_bra[ orbi ] ){
                  bits_up_bra[ orbi ] = 0;

                  //Loop k to add again to the vector
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if (( orb2irrep[ orbi ] == orb2irrep[ orbk ] ) && ( !(bits_up_bra[ orbk ]) )){
                        bits_up_bra[ orbk ] = 1;

                        //find the corresponding string and counter
                        const unsigned int string_in_up = bits2str(L, bits_up_bra);
                        const unsigned long long incounter = incounter_offset + str2cnt_up[ irrep_out_up ][ string_in_up ];
                        
                        if ( workspace[ incounter ] ){
                           workspace[ incounter ] = false;
                           const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work);
                           myResult += factor * factor;
                        }

                        bits_up_bra[ orbk ] = 0;
                     }
                  }
                  bits_up_bra[ orbi ] = 1;
               }
            }
         
         }
         
         // If only 1 electron: Acting with a^+_i,down a_k,down
         if (( Nel_up == 0 ) && ( Nel_down == 1 )){
            const unsigned long long incounter_offset = irrep_center_jumps[ 0 ][ irrep_out_up ] + outcount_up;
            
            //Loop i to remove from the vector
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_down_bra[ orbi ] ){
                  bits_down_bra[ orbi ] = 0;
                        
                  //Loop k to add again to the vector
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if (( orb2irrep[ orbi ] == orb2irrep[ orbk ] ) && ( !(bits_down_bra[ orbk ]) )){
                        bits_down_bra[ orbk ] = 1;
                                                      
                        //find the corresponding string and counter
                        const unsigned int string_in_down = bits2str(L, bits_down_bra);
                        const unsigned long long incounter = incounter_offset + numPerIrrep_up[ irrep_out_up ] * str2cnt_down[ irrep_out_down ][ string_in_down ];
                        
                        if ( workspace[ incounter ] ){
                           workspace[ incounter ] = false;
                           const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work);
                           myResult += factor * factor;
                        }
                        
                        bits_down_bra[ orbk ] = 0;
                     }
                  }

                  bits_down_bra[ orbi ] = 1;
               }
            }
            
         }
         
         // Case 1: Acting with a^+_i,up a^+_j,up a_l,up a_k,up
         if ( Nel_up >= 2 ){
            const unsigned long long incounter_offset = irrep_center_jumps[ 0 ][ irrep_out_up ] + numPerIrrep_up[ irrep_out_up ] * outcount_down;
            
            //Loop i<j to remove from the vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_up_bra[ orbi ] ){
                  bits_up_bra[ orbi ] = 0;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_up_bra[ orbj ] ){
                        bits_up_bra[ orbj ] = 0;
                        const int IrrepProd_ij = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_up_bra[ orbk ]) ){
                              bits_up_bra[ orbk ] = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_up_bra[ orbl ]) ) && ( IrrepProd_ij == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_up_bra[ orbl ] = 1;
                                    
                                    //find the corresponding string and counter
                                    const unsigned int string_in_up = bits2str(L, bits_up_bra);
                                    const unsigned long long incounter = incounter_offset + str2cnt_up[ irrep_out_up ][ string_in_up ];
                                    
                                    if ( workspace[ incounter ] ){
                                       workspace[ incounter ] = false;
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work);
                                       myResult += factor * factor;
                                    }
                                    
                                    bits_up_bra[ orbl ] = 0;
                                 }
                              }
                              bits_up_bra[ orbk ] = 0;
                           }
                        }
                        bits_up_bra[ orbj ] = 1;
                     }
                  }
                  bits_up_bra[ orbi ] = 1;
               }
            }
         
         }
         
         // Case 2: Acting with a^+_i,down a^+_j,down a_l,down a_k,down
         if ( Nel_down >=2 ){
            const unsigned long long incounter_offset = irrep_center_jumps[ 0 ][ irrep_out_up ] + outcount_up;
            
            //Loop i<j to remove from the vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_down_bra[ orbi ] ){
                  bits_down_bra[ orbi ] = 0;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_down_bra[ orbj ] ){
                        bits_down_bra[ orbj ] = 0;
                        const int IrrepProd12 = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_down_bra[ orbk ]) ){
                              bits_down_bra[ orbk ] = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_down_bra[ orbl ]) ) && ( IrrepProd12 == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_down_bra[ orbl ] = 1;
                                    
                                    //find the corresponding string and counter
                                    const unsigned int string_in_down = bits2str(L, bits_down_bra);
                                    const unsigned long long incounter = incounter_offset + numPerIrrep_up[ irrep_out_up ] * str2cnt_down[ irrep_out_down ][ string_in_down ];
                                    
                                    if ( workspace[ incounter ] ){
                                       workspace[ incounter ] = false;
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work);
                                       myResult += factor * factor;
                                    }
                                    
                                    bits_down_bra[ orbl ] = 0;
                                 }
                              }
                              bits_down_bra[ orbk ] = 0;
                           }
                        }
                        bits_down_bra[ orbj ] = 1;
                     }
                  }
                  bits_down_bra[ orbi ] = 1;
               }
            }
            
         }
         
         // Case 3: Acting with a^+_i,up a_k,up a^+_j,down a_l,down
         if (( Nel_up >= 1 ) && ( Nel_down >= 1 )){
         
            //Loop orbi to REMOVE an UP electron from the vector
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_up_bra[ orbi ] ){
                  bits_up_bra[ orbi ] = 0;
                  
                  //Loop orbk to ADD an UP electron again to the vector
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if ( !(bits_up_bra[ orbk ]) ){
                        bits_up_bra[ orbk ] = 1;
                        const int IrrepProd_ik = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbk ] );
                        const int irrep_in_up = getIrrepProduct( irrep_out_up , IrrepProd_ik );
                        const unsigned int string_in_up = bits2str(L, bits_up_bra);
                        const int incount_up = str2cnt_up[ irrep_in_up ][ string_in_up ];
                        
                        //Loop orbj to REMOVE a DOWN electron from the vector
                        for (unsigned int orbj = 0; orbj < L; orbj++){
                           if ( bits_down_bra[ orbj ] ){
                              bits_down_bra[ orbj ] = 0;

                              //Loop orbl to ADD a DOWN electron again to the vector
                              for (unsigned int orbl = 0; orbl < L; orbl++){
                                 const int IrrepProd_jl = getIrrepProduct( orb2irrep[ orbj ] , orb2irrep[ orbl ] );
                                 if (( bits_down_bra[ orbl ] == 0 ) && ( IrrepProd_ik == IrrepProd_jl )){
                                    bits_down_bra[ orbl ] = 1;
                                    
                                    const int irrep_in_down = getIrrepProduct( irrep_out_down , IrrepProd_jl );
                                    const unsigned int string_in_down = bits2str(L, bits_down_bra);
                                    const int incount_down = str2cnt_down[ irrep_in_down ][ string_in_down ];
                                    const unsigned long long incounter = irrep_center_jumps[ 0 ][ irrep_in_up ] + incount_up + numPerIrrep_up[ irrep_in_up ] * incount_down;
                                    
                                    if ( workspace[ incounter ] ){
                                       workspace[ incounter ] = false;
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work);
                                       myResult += factor * factor;
                                    }
                                    
                                    bits_down_bra[ orbl ] = 0;
                                 }
                              }
                              bits_down_bra[ orbj ] = 1;
                           }
                        }
                        bits_up_bra[ orbk ] = 0;
                     }
                  }
                  bits_up_bra[ orbi ] = 1;
               }
            }
         
         }
         
         output[ outcounter ] = myResult;
      
      }
      
      delete [] bits_up_ket;
      delete [] bits_down_ket;
      delete [] bits_up_bra;
      delete [] bits_down_bra;
      delete [] work;
      delete [] workspace;
   
   }

}


unsigned long long CheMPS2::FCI::LowestEnergyDeterminant() const{

   const unsigned long long vecLength = getVecLength( 0 );
   double * energies = new double[ vecLength ];

   // Fetch the Slater determinant energies
   DiagHam( energies );
   
   // Find the determinant with minimum energy
   unsigned long long minEindex = 0;
   for ( unsigned long long count = 1; count < vecLength; count++ ){
      if ( energies[ count ] < energies[ minEindex ] ){
         minEindex = count;
      }
   }
   
   delete [] energies;
   
   return minEindex;

}

double CheMPS2::FCI::GetMatrixElement(int * bits_bra_up, int * bits_bra_down, int * bits_ket_up, int * bits_ket_down, int * work) const{
   
   int count_annih_up   = 0;
   int count_creat_up   = 0;
   int count_annih_down = 0;
   int count_creat_down = 0;
   
   int * annih_up   = work;
   int * creat_up   = work+2;
   int * annih_down = work+4;
   int * creat_down = work+6;
   
   // Find the differences between bra and ket, and store them in the arrays annih/creat
   for (unsigned int orb = 0; orb < L; orb++){
      if ( bits_bra_up[ orb ] != bits_ket_up[ orb ] ){
         if ( bits_ket_up[ orb ] ){ // Electron is annihilated in the ket
            if ( count_annih_up == 2 ){ return 0.0; }
            annih_up[ count_annih_up ] = orb;
            count_annih_up++;
         } else { // Electron is created in the ket
            if ( count_creat_up == 2 ){ return 0.0; }
            creat_up[ count_creat_up ] = orb;
            count_creat_up++;
         }
      }
      if ( bits_bra_down[ orb ] != bits_ket_down[ orb ] ){
         if ( bits_ket_down[ orb ] ){ // Electron is annihilated in the ket
            if ( count_annih_down == 2 ){ return 0.0; }
            annih_down[ count_annih_down ] = orb;
            count_annih_down++;
         } else { // Electron is created in the ket
            if ( count_creat_down == 2 ){ return 0.0; }
            creat_down[ count_creat_down ] = orb;
            count_creat_down++;
         }
      }
   }
   
   // Sanity check: spin symmetry
   if ( count_annih_up   != count_creat_up   ){ return 0.0; }
   if ( count_annih_down != count_creat_down ){ return 0.0; }
   
   // Sanity check: At most 2 annihilators and 2 creators can connect the ket and bra
   if ( count_annih_up + count_annih_down > 2 ){ return 0.0; }
   if ( count_creat_up + count_creat_down > 2 ){ return 0.0; }
   
   if (( count_annih_up == 0 ) && ( count_annih_down == 0 )){ // |bra> == |ket>  -->  copy from DiagHam( )
   
      double result = 0.0;
      for ( unsigned int orb1 = 0; orb1 < L; orb1++ ){
         const int n_tot_orb1 = bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ];
         result += n_tot_orb1 * getGmat( orb1 , orb1 );
         for ( unsigned int orb2 = 0; orb2 < L; orb2++ ){
            result += 0.5 * n_tot_orb1 * ( bits_ket_up[ orb2 ] + bits_ket_down[ orb2 ] ) * getERI( orb1 , orb1 , orb2 , orb2 )
            + 0.5 * ( n_tot_orb1 - bits_ket_up[ orb1 ] * bits_ket_up[ orb2 ] - bits_ket_down[ orb1 ] * bits_ket_down[ orb2 ] ) * getERI( orb1 , orb2 , orb2 , orb1 );
         }
      }
      return result;
      
   }
   
   if (( count_annih_up == 1 ) && ( count_annih_down == 0 )){ // |bra> = a^+_j,up a_l,up |ket>
   
      const int orbj = creat_up[ 0 ];
      const int orbl = annih_up[ 0 ];
   
      double result = getGmat( orbj , orbl );
      for (unsigned int orb1 = 0; orb1 < L; orb1++){
         result += getERI(orbj, orb1, orb1, orbl) * ( 0.5 - bits_ket_up[ orb1 ] ) + getERI(orb1, orb1, orbj, orbl) * ( bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ] );
      }
      int phase = 1;
      if ( orbj < orbl ){ for (int orbital = orbj+1; orbital < orbl; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } }
      if ( orbl < orbj ){ for (int orbital = orbl+1; orbital < orbj; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } }
      return ( result * phase );
   
   }
   
   if (( count_annih_up == 0 ) && ( count_annih_down == 1 )){ // |bra> = a^+_j,down a_l,down |ket>
   
      const int orbj = creat_down[ 0 ];
      const int orbl = annih_down[ 0 ];
   
      double result = getGmat( orbj , orbl );
      for (unsigned int orb1 = 0; orb1 < L; orb1++){
         result += getERI(orbj, orb1, orb1, orbl) * ( 0.5 - bits_ket_down[ orb1 ] ) + getERI(orb1, orb1, orbj, orbl) * ( bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ] );
      }
      int phase = 1;
      if ( orbj < orbl ){ for (int orbital = orbj+1; orbital < orbl; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      if ( orbl < orbj ){ for (int orbital = orbl+1; orbital < orbj; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      return ( result * phase );
   
   }
   
   if (( count_annih_up == 2 ) && ( count_annih_down == 0 )){
   
      // creat and annih are filled in increasing orbital index
      const int orbi = creat_up[ 0 ];
      const int orbj = creat_up[ 1 ];
      const int orbk = annih_up[ 0 ];
      const int orbl = annih_up[ 1 ];
      
      double result = getERI(orbi, orbk, orbj, orbl) - getERI(orbi, orbl, orbj, orbk);
      int phase = 1;
      for (int orbital = orbk+1; orbital < orbl; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } // Fermion phases orbk and orbl measured in the ket
      for (int orbital = orbi+1; orbital < orbj; orbital++){ if ( bits_bra_up[ orbital ] ){ phase *= -1; } } // Fermion phases orbi and orbj measured in the bra
      return ( result * phase );
   
   }
   
   if (( count_annih_up == 0 ) && ( count_annih_down == 2 )){
   
      // creat and annih are filled in increasing orbital index
      const int orbi = creat_down[ 0 ];
      const int orbj = creat_down[ 1 ];
      const int orbk = annih_down[ 0 ];
      const int orbl = annih_down[ 1 ];
      
      double result = getERI(orbi, orbk, orbj, orbl) - getERI(orbi, orbl, orbj, orbk);
      int phase = 1;
      for (int orbital = orbk+1; orbital < orbl; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } // Fermion phases orbk and orbl measured in the ket
      for (int orbital = orbi+1; orbital < orbj; orbital++){ if ( bits_bra_down[ orbital ] ){ phase *= -1; } } // Fermion phases orbi and orbj measured in the bra
      return ( result * phase );
   
   }
   
   if (( count_annih_up == 1 ) && ( count_annih_down == 1 )){
   
      const int orbi = creat_up  [ 0 ];
      const int orbj = creat_down[ 0 ];
      const int orbk = annih_up  [ 0 ];
      const int orbl = annih_down[ 0 ];
      
      double result = getERI(orbi, orbk, orbj, orbl);
      int phase = 1;
      if ( orbi < orbk ){ for (int orbital = orbi+1; orbital < orbk; orbital++){ if ( bits_ket_up  [ orbital ] ){ phase *= -1; } } }
      if ( orbk < orbi ){ for (int orbital = orbk+1; orbital < orbi; orbital++){ if ( bits_ket_up  [ orbital ] ){ phase *= -1; } } }
      if ( orbj < orbl ){ for (int orbital = orbj+1; orbital < orbl; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      if ( orbl < orbj ){ for (int orbital = orbl+1; orbital < orbj; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      return ( result * phase );
   
   }
   
   return 0.0;

}

void CheMPS2::FCI::FCIdcopy(const unsigned long long vecLength, double * origin, double * target){

   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   dcopy_( &length , origin , &inc , target , &inc );

}

double CheMPS2::FCI::FCIddot(const unsigned long long vecLength, double * vec1, double * vec2){

   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   return ddot_( &length , vec1 , &inc , vec2 , &inc );

}

double CheMPS2::FCI::FCIfrobeniusnorm(const unsigned long long vecLength, double * vec){

   return sqrt( FCIddot( vecLength , vec , vec ) );

}

void CheMPS2::FCI::FCIdaxpy(const unsigned long long vecLength, const double alpha, double * vec_x, double * vec_y){

   double factor = alpha;
   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   daxpy_( &length , &factor , vec_x , &inc , vec_y , &inc );

}

void CheMPS2::FCI::FCIdscal(const unsigned long long vecLength, const double alpha, double * vec){

   double factor = alpha;
   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   dscal_( &length , &factor , vec , &inc );

}

void CheMPS2::FCI::ClearVector(const unsigned long long vecLength, double * vec){

   for ( unsigned long long cnt = 0; cnt < vecLength; cnt++ ){ vec[cnt] = 0.0; }

}

void CheMPS2::FCI::FillRandom(const unsigned long long vecLength, double * vec){

   for ( unsigned long long cnt = 0; cnt < vecLength; cnt++ ){ vec[cnt] = ( ( 2.0 * rand() ) / RAND_MAX ) - 1.0; }

}

double CheMPS2::FCI::GSDavidson(double * inoutput, const int DAVIDSON_NUM_VEC) const{

   const int veclength = getVecLength( 0 ); // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   const double RTOL   = CheMPS2::HEFF_DAVIDSON_RTOL_BASE * sqrt( 1.0 * veclength );
   
   Davidson deBoskabouter( veclength, DAVIDSON_NUM_VEC, CheMPS2::HEFF_DAVIDSON_NUM_VEC_KEEP, RTOL, CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF, false ); // No debug printing for FCI
   double ** whichpointers = new double*[2];
   
   char instruction = deBoskabouter.FetchInstruction( whichpointers );
   assert( instruction == 'A' );
   if ( inoutput != NULL ){ FCIdcopy( veclength, inoutput, whichpointers[0] ); }
   else { FillRandom( veclength, whichpointers[0] ); }
   DiagHam( whichpointers[1] );
   
   instruction = deBoskabouter.FetchInstruction( whichpointers );
   while ( instruction == 'B' ){
      HamTimesVec( whichpointers[0], whichpointers[1] );
      instruction = deBoskabouter.FetchInstruction( whichpointers );
   }
   
   assert( instruction == 'C' );
   if ( inoutput != NULL ){ FCIdcopy( veclength, whichpointers[0], inoutput ); }
   const double FCIenergy = whichpointers[1][0] + getEconst();
   if ( FCIverbose > 1 ){ cout << "FCI::GSDavidson : Required number of matrix-vector multiplications = " << deBoskabouter.GetNumMultiplications() << endl; }
   if ( FCIverbose > 0 ){ cout << "FCI::GSDavidson : Converged ground state energy = " << FCIenergy << endl; }
   delete [] whichpointers;
   return FCIenergy;

}

/*********************************************************************************
 *                                                                               *
 *   Below this block all functions are for the Green's function calculations.   *
 *                                                                               *
 *********************************************************************************/

void CheMPS2::FCI::ActWithNumberOperator(const unsigned int orbIndex, double * resultVector, double * sourceVector) const{

   assert( orbIndex<L );

   int * bits_up    = new int[ L ];
   int * bits_down  = new int[ L ];

   const unsigned long long vecLength = getVecLength( 0 );
   for (unsigned long long counter = 0; counter < vecLength; counter++){
      getBitsOfCounter( 0 , counter , bits_up , bits_down );
      resultVector[ counter ] = ( bits_up[ orbIndex ] + bits_down[ orbIndex ] ) * sourceVector[ counter ];
   }
   
   delete [] bits_up;
   delete [] bits_down;

}

void CheMPS2::FCI::ActWithSecondQuantizedOperator(const char whichOperator, const bool isUp, const unsigned int orbIndex, double * thisVector, const FCI * otherFCI, double * otherVector) const{

   assert( ( whichOperator=='C' ) || ( whichOperator=='A' ) ); //Operator should be a (C) Creator, or (A) Annihilator
   assert( orbIndex<L  );
   assert( L==otherFCI->getL() );

   const unsigned long long vecLength = getVecLength( 0 );

   if ( getTargetIrrep() != getIrrepProduct( otherFCI->getTargetIrrep() , getOrb2Irrep( orbIndex ) )){
      ClearVector( vecLength , thisVector );
      return;
   }

   int * bits_up    = new int[ L ];
   int * bits_down  = new int[ L ];
   
   if (( whichOperator=='C') && ( isUp )){
      for (unsigned long long counter = 0; counter < vecLength; counter++){
         
         getBitsOfCounter( 0 , counter , bits_up , bits_down );
         
         if ( bits_up[ orbIndex ] == 1 ){ // Operator = creator_up
            bits_up[ orbIndex ] = 0;
            int phase = 1;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_up[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         } else {
            thisVector[ counter ] = 0.0;
         }
         
      }
   }
   
   if (( whichOperator=='C') && ( !(isUp) )){
      const int startphase = (( Nel_up % 2 ) == 0) ? 1 : -1;
      for (unsigned long long counter = 0; counter < vecLength; counter++){

         getBitsOfCounter( 0 , counter , bits_up , bits_down );

         if ( bits_down[ orbIndex ] == 1 ){ // Operator = creator_down
            bits_down[ orbIndex ] = 0;
            int phase = startphase;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_down[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         } else {
            thisVector[ counter ] = 0.0;
         }

      }
   }
   
   if (( whichOperator=='A') && ( isUp )){
      for (unsigned long long counter = 0; counter < vecLength; counter++){

         getBitsOfCounter( 0 , counter , bits_up , bits_down );
         
         if ( bits_up[ orbIndex ] == 0 ){ // Operator = annihilator_up
            bits_up[ orbIndex ] = 1;
            int phase = 1;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_up[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         } else {
            thisVector[ counter ] = 0.0;
         }

      }
   }
   
   if (( whichOperator=='A') && ( !(isUp) )){
      const int startphase = (( Nel_up % 2 ) == 0) ? 1 : -1;
      for (unsigned long long counter = 0; counter < vecLength; counter++){

         getBitsOfCounter( 0 , counter , bits_up , bits_down );
         
         if ( bits_down[ orbIndex ] == 0 ){ // Operator = annihilator_down
            bits_down[ orbIndex ] = 1;
            int phase = startphase;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_down[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         } else {
            thisVector[ counter ] = 0.0;
         }
         
      }
   }
   
   delete [] bits_up;
   delete [] bits_down;

}

void CheMPS2::FCI::CGSolveSystem(const double alpha, const double beta, const double eta, double * RHS, double * RealSol, double * ImagSol, const bool checkError) const{

   const unsigned long long vecLength = getVecLength( 0 );
   
   // Create a few helper arrays
   double * RESID  = new double[ vecLength ];
   double * PVEC   = new double[ vecLength ];
   double * OxPVEC = new double[ vecLength ];
   double * temp   = new double[ vecLength ];
   double * temp2  = new double[ vecLength ];
   double * precon = new double[ vecLength ];
   CGDiagPrecond( alpha , beta , eta , precon , temp );
   
   assert( RealSol != NULL );
   assert( ImagSol != NULL );
   assert( fabs( eta ) > 0.0 );

   /* 
         ( alpha + beta H + I eta ) Solution = RHS
      
      is solved with the conjugate gradient (CG) method. Solution = RealSol + I * ImagSol.
      CG requires a symmetric positive definite operator. Therefore:
      
         precon * [ ( alpha + beta H )^2 + eta^2 ] * precon * SolutionTilde = precon * ( alpha + beta H - I eta ) * RHS
         Solution = precon * SolutionTilde
         
      Clue: Solve for ImagSol first. RealSol is then simply
      
         RealSol = - ( alpha + beta H ) / eta * ImagSol
   */
   
   /**** Solve for ImagSol ****/
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RESID[ cnt ] = - eta * precon[ cnt ] * RHS[ cnt ]; } // RESID = - eta * precon * RHS
   if ( FCIverbose > 1 ){ cout << "FCI::CGSolveSystem : Two-norm of the RHS for the imaginary part = " << FCIfrobeniusnorm( vecLength , RESID ) << endl; }
   FCIdcopy( vecLength, RESID, ImagSol ); // Well educated initial guess for the imaginary part (guess is exact if operator is diagonal)
   CGCoreSolver( alpha, beta, eta, precon, ImagSol, RESID, PVEC, OxPVEC, temp, temp2 ); // RESID contains the RHS of ( precon * Op * precon ) * |x> = |b>
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ ImagSol[ cnt ] = precon[ cnt ] * ImagSol[ cnt ]; }
   
   /**** Solve for RealSol ****/
   CGAlphaPlusBetaHAM( -alpha/eta, -beta/eta, ImagSol, RealSol ); // Initial guess RealSol can be obtained from ImagSol
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
      if ( fabs( precon[cnt] ) > CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF ){
         RealSol[cnt] = RealSol[cnt] / precon[cnt];
      } else {
         RealSol[cnt] = RealSol[cnt] / CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF;
      }
   }
   CGAlphaPlusBetaHAM( alpha, beta, RHS, RESID ); // RESID = ( alpha + beta * H ) * RHS
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RESID[ cnt ] = precon[ cnt ] * RESID[ cnt ]; } // RESID = precon * ( alpha + beta * H ) * RHS
   if ( FCIverbose > 1 ){ cout << "FCI::CGSolveSystem : Two-norm of the RHS for the real part = " << FCIfrobeniusnorm( vecLength , RESID ) << endl; }
   CGCoreSolver( alpha, beta, eta, precon, RealSol, RESID, PVEC, OxPVEC, temp, temp2 ); // RESID contains the RHS of ( precon * Op * precon ) * |x> = |b>
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RealSol[ cnt ] = precon[ cnt ] * RealSol[ cnt ]; }
   
   if (( checkError ) && ( FCIverbose > 0 )){
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ precon[ cnt ] = 1.0; }
      CGOperator( alpha , beta , eta , precon , RealSol , temp , temp2 , OxPVEC );
      CGAlphaPlusBetaHAM( alpha , beta , RHS , RESID );
      FCIdaxpy( vecLength , -1.0 , RESID , OxPVEC );
      double RMSerror = FCIddot( vecLength , OxPVEC , OxPVEC );
      CGOperator( alpha , beta , eta , precon , ImagSol , temp , temp2 , OxPVEC );
      FCIdaxpy( vecLength , eta , RHS , OxPVEC );
      RMSerror += FCIddot( vecLength , OxPVEC , OxPVEC );
      RMSerror = sqrt( RMSerror );
      cout << "FCI::CGSolveSystem : RMS error when checking the solution (without preconditioner) = " << RMSerror << endl;
   }
   
   // Clean up
   delete [] temp;
   delete [] temp2;
   delete [] RESID;
   delete [] PVEC;
   delete [] OxPVEC;
   delete [] precon;

}

void CheMPS2::FCI::CGCoreSolver(const double alpha, const double beta, const double eta, double * precon, double * Sol, double * RESID, double * PVEC, double * OxPVEC, double * temp, double * temp2) const{

   const unsigned long long vecLength = getVecLength( 0 );
   const double CGRESIDUALTHRESHOLD = 100.0 * CheMPS2::HEFF_DAVIDSON_RTOL_BASE * sqrt( 1.0 * vecLength );
   if ( FCIverbose>1 ){ cout << "FCI::CGCoreSolver : The residual norm for convergence = " << CGRESIDUALTHRESHOLD << endl; }

   /*
      Operator = precon * [ ( alpha + beta * H )^2 + eta^2 ] * precon : positive definite and symmetric
      x_0 = Sol
      p_0 (PVEC) = r_0 (RESID) = b - Operator * x_0
      k (count_k) = 0
   */

   int count_k = 0;
   CGOperator( alpha , beta , eta , precon , Sol , temp , temp2 , OxPVEC ); // O_p_k = Operator * Sol
   FCIdaxpy( vecLength , -1.0 , OxPVEC , RESID );                           // r_0 = b - Operator * x_0
   FCIdcopy( vecLength , RESID , PVEC );                                    // p_0 = r_0
   double rkT_rk = FCIddot( vecLength , RESID , RESID );
   double residualNorm = sqrt( rkT_rk );
   
   while ( residualNorm >= CGRESIDUALTHRESHOLD ){
   
      CGOperator( alpha , beta , eta , precon, PVEC , temp , temp2 , OxPVEC ); // O_p_k = Operator * p_k
      const double alpha_k = rkT_rk / FCIddot( vecLength , PVEC , OxPVEC );    // alpha_k = r_k^T * r_k / ( p_k^T * O_p_k )
      FCIdaxpy( vecLength ,  alpha_k , PVEC   , Sol   );                       // x_{k+1} = x_k + alpha_k * p_k
      FCIdaxpy( vecLength , -alpha_k , OxPVEC , RESID );                       // r_{k+1} = r_k - alpha_k * A * p_k
      const double rkplus1T_rkplus1 = FCIddot( vecLength , RESID , RESID );
      const double beta_k = rkplus1T_rkplus1 / rkT_rk ;                        // beta_k = r_{k+1}^T * r_{k+1} / ( r_k^T * r_k )
      for ( unsigned long long cnt = 0; cnt < vecLength; cnt++ ){
         PVEC[ cnt ] = RESID[ cnt ] + beta_k * PVEC[ cnt ];                    // p_{k+1} = r_{k+1} + beta_k * p_k
      }
      count_k++;
      rkT_rk = rkplus1T_rkplus1;
      residualNorm = sqrt( rkT_rk );
      if ( FCIverbose > 1 ){ cout << "FCI::CGCoreSolver : At step " << count_k << " the residual norm is " << residualNorm << endl; }
      
   }

}

void CheMPS2::FCI::CGAlphaPlusBetaHAM(const double alpha, const double beta, double * in, double * out) const{

   HamTimesVec( in , out );
   const unsigned long long vecLength = getVecLength( 0 );
   const double prefactor = alpha + beta * getEconst(); // HamTimesVec does only the parts with second quantized operators
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
      out[ cnt ] = prefactor * in[ cnt ] + beta * out[ cnt ]; // out = ( alpha + beta * H ) * in
   }

}

void CheMPS2::FCI::CGOperator(const double alpha, const double beta, const double eta, double * precon, double * in, double * temp, double * temp2, double * out) const{

   const unsigned long long vecLength = getVecLength( 0 );
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
      temp[ cnt ] = precon[ cnt ] * in[ cnt ];                 // temp  = precon * in
   }
   CGAlphaPlusBetaHAM( alpha , beta , temp  , temp2 );         // temp2 = ( alpha + beta * H )   * precon * in
   CGAlphaPlusBetaHAM( alpha , beta , temp2 , out   );         // out   = ( alpha + beta * H )^2 * precon * in
   FCIdaxpy( vecLength , eta*eta , temp , out );               // out   = [ ( alpha + beta * H )^2 + eta*eta ] * precon * in
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
      out[ cnt ] = precon[ cnt ] * out[ cnt ];                 // out   = precon * [ ( alpha + beta * H )^2 + eta*eta ] * precon * in
   }

}

void CheMPS2::FCI::CGDiagPrecond(const double alpha, const double beta, const double eta, double * precon, double * workspace) const{

   // With operator = [ ( alpha + beta * H )^2 + eta*eta ] ; precon becomes 1 / sqrt( diag ( operator ) ).
   
   DiagHam( precon );
   DiagHamSquared( workspace );
   
   const unsigned long long vecLength = getVecLength( 0 );
   const double alpha_bis = alpha + beta * getEconst();
   const double factor1 = alpha_bis * alpha_bis + eta * eta;
   const double factor2 = 2 * alpha_bis * beta;
   const double factor3 = beta * beta;
   for (unsigned long long row = 0; row < vecLength; row++){
      const double diagElement = factor1 + factor2 * precon[ row ] + factor3 * workspace[ row ];
      precon[ row ] = 1.0 / sqrt( diagElement );
   }
   
   if ( FCIverbose>1 ){
      double minval = precon[0];
      double maxval = precon[0];
      for (unsigned long long cnt = 1; cnt < vecLength; cnt++){
         if ( precon[ cnt ] > maxval ){ maxval = precon[ cnt ]; }
         if ( precon[ cnt ] < minval ){ minval = precon[ cnt ]; }
      }
      cout << "FCI::CGDiagPrecond : Minimum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << 1.0/(maxval*maxval) << endl;
      cout << "FCI::CGDiagPrecond : Maximum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << 1.0/(minval*minval) << endl;
   }

}

void CheMPS2::FCI::RetardedGF(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF) const{

   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );

   // G( omega, alpha, beta, eta ) = < 0 | a_{alpha,spin}  [ omega - Ham + E_0 + I*eta ]^{-1} a^+_{beta,spin} | 0 > (addition amplitude)
   //                              + < 0 | a^+_{beta,spin} [ omega + Ham - E_0 + I*eta ]^{-1} a_{alpha,spin}  | 0 > (removal  amplitude)

   double Realpart, Imagpart;
   RetardedGF_addition(omega, eta, orb_alpha, orb_beta, isUp, GSenergy, GSvector, Ham, &Realpart, &Imagpart, NULL, NULL);
   RePartGF[0] = Realpart; // Set
   ImPartGF[0] = Imagpart; // Set
   
   RetardedGF_removal( omega, eta, orb_alpha, orb_beta, isUp, GSenergy, GSvector, Ham, &Realpart, &Imagpart, NULL, NULL);
   RePartGF[0] += Realpart; // Add
   ImPartGF[0] += Imagpart;
   
   if ( FCIverbose>0 ){
      cout << "FCI::RetardedGF : G( omega = " << omega << " ; eta = " << eta << " ; i = " << orb_alpha << " ; j = " << orb_beta << " ) = " << RePartGF[0] << " + I * " << ImPartGF[0] << endl;
      cout << "                  Local density of states (LDOS) = " << - ImPartGF[0] / M_PI << endl;
   }

}

void CheMPS2::FCI::RetardedGF_addition(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag) const{

   // Addition amplitude < 0 | a_{alpha, spin} [ omega - Ham + E_0 + I*eta ]^{-1} a^+_{beta, spin} | 0 >

   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const bool isOK = ( isUp ) ? ( getNel_up() < L ) : ( getNel_down() < L ); // The electron can be added
   if (( getOrb2Irrep( orb_alpha ) != getOrb2Irrep( orb_beta ) ) || ( !isOK )){
      RePartGF[0] = 0.0;
      ImPartGF[0] = 0.0;
      const unsigned int Lpow4 = L*L*L*L;
      if ( TwoRDMreal != NULL ){ for ( unsigned int cnt = 0; cnt < Lpow4; cnt++ ){ TwoRDMreal[ cnt ] = 0.0; } }
      if ( TwoRDMimag != NULL ){ for ( unsigned int cnt = 0; cnt < Lpow4; cnt++ ){ TwoRDMimag[ cnt ] = 0.0; } }
      return;
   }

   const unsigned int addNelUP   = getNel_up()   + ((isUp) ? 1 : 0);
   const unsigned int addNelDOWN = getNel_down() + ((isUp) ? 0 : 1);
   const int addIrrep = getIrrepProduct( getTargetIrrep() , getOrb2Irrep( orb_beta ) );

   CheMPS2::FCI additionFCI( Ham , addNelUP , addNelDOWN , addIrrep , maxMemWorkMB , FCIverbose );
   const unsigned long long addVecLength = additionFCI.getVecLength( 0 );
   double * addBetaVector  = new double[ addVecLength ];
   double * addAlphaVector = ( orb_alpha == orb_beta ) ? addBetaVector : new double[ addVecLength ];
   additionFCI.ActWithSecondQuantizedOperator( 'C' , isUp , orb_beta , addBetaVector , this , GSvector ); // | addBetaVector > = a^+_beta,spin | GSvector >
   if ( orb_alpha != orb_beta ){
      additionFCI.ActWithSecondQuantizedOperator( 'C' , isUp , orb_alpha , addAlphaVector , this , GSvector ); // | addAlphaVector > = a^+_alpha,spin | GSvector >
   }
   
   double * RealPartSolution = new double[ addVecLength ];
   double * ImagPartSolution = new double[ addVecLength ];
   additionFCI.CGSolveSystem( omega + GSenergy , -1.0 , eta , addBetaVector , RealPartSolution , ImagPartSolution );
   if ( TwoRDMreal != NULL ){ additionFCI.Fill2RDM( RealPartSolution , TwoRDMreal ); } // Sets the TwoRDMreal
   RePartGF[0] = FCIddot( addVecLength , addAlphaVector , RealPartSolution );
   delete [] RealPartSolution;
   if ( TwoRDMimag != NULL ){ additionFCI.Fill2RDM( ImagPartSolution , TwoRDMimag ); } // Sets the TwoRDMimag
   ImPartGF[0] = FCIddot( addVecLength , addAlphaVector , ImagPartSolution );
   delete [] ImagPartSolution;

   if ( orb_alpha != orb_beta ){ delete [] addAlphaVector; }
   delete [] addBetaVector;

}

void CheMPS2::FCI::RetardedGF_removal(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag) const{

   // Removal amplitude < 0 | a^+_{beta, spin} [ omega + Ham - E_0 + I*eta ]^{-1} a_{alpha, spin} | 0 >

   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const bool isOK = ( isUp ) ? ( getNel_up() > 0 ) : ( getNel_down() > 0 ); // The electron can be removed
   if (( getOrb2Irrep( orb_alpha ) != getOrb2Irrep( orb_beta ) ) || ( !isOK )){
      RePartGF[0] = 0.0;
      ImPartGF[0] = 0.0;
      const unsigned int Lpow4 = L*L*L*L;
      if ( TwoRDMreal != NULL ){ for ( unsigned int cnt = 0; cnt < Lpow4; cnt++ ){ TwoRDMreal[ cnt ] = 0.0; } }
      if ( TwoRDMimag != NULL ){ for ( unsigned int cnt = 0; cnt < Lpow4; cnt++ ){ TwoRDMimag[ cnt ] = 0.0; } }
      return;
   }
   
   const unsigned int removeNelUP   = getNel_up()   - ((isUp) ? 1 : 0);
   const unsigned int removeNelDOWN = getNel_down() - ((isUp) ? 0 : 1);
   const int removeIrrep  = getIrrepProduct( getTargetIrrep() , getOrb2Irrep( orb_alpha ) );
   
   CheMPS2::FCI removalFCI( Ham , removeNelUP , removeNelDOWN , removeIrrep , maxMemWorkMB , FCIverbose );
   const unsigned long long removeVecLength = removalFCI.getVecLength( 0 );
   double * removeAlphaVector = new double[ removeVecLength ];
   double * removeBetaVector  = ( orb_alpha == orb_beta ) ? removeAlphaVector : new double[ removeVecLength ];
   removalFCI.ActWithSecondQuantizedOperator( 'A' , isUp , orb_alpha , removeAlphaVector , this , GSvector ); // | removeAlphaVector > = a_alpha,spin | GSvector >
   if ( orb_alpha != orb_beta ){
      removalFCI.ActWithSecondQuantizedOperator( 'A' , isUp , orb_beta , removeBetaVector , this , GSvector ); // | removeBetaVector > = a_beta,spin | GSvector >
   }
   
   double * RealPartSolution = new double[ removeVecLength ];
   double * ImagPartSolution = new double[ removeVecLength ];
   removalFCI.CGSolveSystem( omega - GSenergy , 1.0 , eta , removeAlphaVector , RealPartSolution , ImagPartSolution );
   if ( TwoRDMreal != NULL ){ removalFCI.Fill2RDM( RealPartSolution , TwoRDMreal ); } // Sets the TwoRDMreal
   RePartGF[0] = FCIddot( removeVecLength , removeBetaVector , RealPartSolution );
   delete [] RealPartSolution;
   if ( TwoRDMimag != NULL ){ removalFCI.Fill2RDM( ImagPartSolution , TwoRDMimag ); } // Sets the TwoRDMimag
   ImPartGF[0] = FCIddot( removeVecLength , removeBetaVector , ImagPartSolution );
   delete [] ImagPartSolution;

   if ( orb_alpha != orb_beta ){ delete [] removeBetaVector; }
   delete [] removeAlphaVector;

}

void CheMPS2::FCI::DensityResponseGF(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF) const{

   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );

   // X( omega, alpha, beta, eta ) = < 0 | ( n_alpha - <0| n_alpha |0> ) [ omega - Ham + E_0 + I*eta ]^{-1} ( n_beta  - <0| n_beta  |0> ) | 0 > (forward  amplitude)
   //                              - < 0 | ( n_beta  - <0| n_beta  |0> ) [ omega + Ham - E_0 + I*eta ]^{-1} ( n_alpha - <0| n_alpha |0> ) | 0 > (backward amplitude)
   
   double Realpart, Imagpart;
   DensityResponseGF_forward( omega, eta, orb_alpha, orb_beta, GSenergy, GSvector, &Realpart, &Imagpart, NULL, NULL);
   RePartGF[0] = Realpart; // Set
   ImPartGF[0] = Imagpart; // Set
   
   DensityResponseGF_backward(omega, eta, orb_alpha, orb_beta, GSenergy, GSvector, &Realpart, &Imagpart, NULL, NULL);
   RePartGF[0] -= Realpart; // Subtract !!!
   ImPartGF[0] -= Imagpart; // Subtract !!!

   if ( FCIverbose>0 ){
      cout << "FCI::DensityResponseGF : X( omega = " << omega << " ; eta = " << eta << " ; i = " << orb_alpha << " ; j = " << orb_beta << " ) = " << RePartGF[0] << " + I * " << ImPartGF[0] << endl;
      cout << "                         Local density-density response (LDDR) = " << - ImPartGF[0] / M_PI << endl;
   }

}

void CheMPS2::FCI::DensityResponseGF_forward(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag) const{

   // Forward amplitude: < 0 | ( n_alpha - <0| n_alpha |0> ) [ omega - Ham + E_0 + I*eta ]^{-1} ( n_beta  - <0| n_beta  |0> ) | 0 >
   
   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const unsigned long long vecLength = getVecLength( 0 );
   double * densityAlphaVector = new double[ vecLength ];
   double * densityBetaVector  = ( orb_alpha == orb_beta ) ? densityAlphaVector : new double[ vecLength ];
   ActWithNumberOperator( orb_alpha , densityAlphaVector , GSvector );             // densityAlphaVector = n_alpha |0>
   const double n_alpha_0 = FCIddot( vecLength , densityAlphaVector , GSvector );  // <0| n_alpha |0>
   FCIdaxpy( vecLength , -n_alpha_0 , GSvector , densityAlphaVector );             // densityAlphaVector = ( n_alpha - <0| n_alpha |0> ) |0>
   if ( orb_alpha != orb_beta ){
      ActWithNumberOperator( orb_beta , densityBetaVector , GSvector );            // densityBetaVector = n_beta |0>
      const double n_beta_0 = FCIddot( vecLength , densityBetaVector , GSvector ); // <0| n_beta |0>
      FCIdaxpy( vecLength , -n_beta_0 , GSvector , densityBetaVector );            // densityBetaVector = ( n_beta - <0| n_beta |0> ) |0>
   }

   double * RealPartSolution = new double[ vecLength ];
   double * ImagPartSolution = new double[ vecLength ];
   CGSolveSystem( omega + GSenergy , -1.0 , eta , densityBetaVector , RealPartSolution , ImagPartSolution );
   if ( TwoRDMreal != NULL ){ Fill2RDM( RealPartSolution , TwoRDMreal ); } // Sets the TwoRDMreal
   RePartGF[0] = FCIddot( vecLength , densityAlphaVector , RealPartSolution );
   delete [] RealPartSolution;
   if ( TwoRDMimag != NULL ){ Fill2RDM( ImagPartSolution , TwoRDMimag ); } // Sets the TwoRDMimag
   ImPartGF[0] = FCIddot( vecLength , densityAlphaVector , ImagPartSolution );
   delete [] ImagPartSolution;

   if ( orb_alpha != orb_beta ){ delete [] densityBetaVector; }
   delete [] densityAlphaVector;

}

void CheMPS2::FCI::DensityResponseGF_backward(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag) const{

   // Backward amplitude: < 0 | ( n_beta  - <0| n_beta  |0> ) [ omega + Ham - E_0 + I*eta ]^{-1} ( n_alpha - <0| n_alpha |0> ) | 0 >
   
   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const unsigned long long vecLength = getVecLength( 0 );
   double * densityAlphaVector = new double[ vecLength ];
   double * densityBetaVector  = ( orb_alpha == orb_beta ) ? densityAlphaVector : new double[ vecLength ];
   ActWithNumberOperator( orb_alpha , densityAlphaVector , GSvector );             // densityAlphaVector = n_alpha |0>
   const double n_alpha_0 = FCIddot( vecLength , densityAlphaVector , GSvector );  // <0| n_alpha |0>
   FCIdaxpy( vecLength , -n_alpha_0 , GSvector , densityAlphaVector );             // densityAlphaVector = ( n_alpha - <0| n_alpha |0> ) |0>
   if ( orb_alpha != orb_beta ){
      ActWithNumberOperator( orb_beta , densityBetaVector , GSvector );            // densityBetaVector = n_beta |0>
      const double n_beta_0 = FCIddot( vecLength , densityBetaVector , GSvector ); // <0| n_beta |0>
      FCIdaxpy( vecLength , -n_beta_0 , GSvector , densityBetaVector );            // densityBetaVector = ( n_beta - <0| n_beta |0> ) |0>
   }

   double * RealPartSolution = new double[ vecLength ];
   double * ImagPartSolution = new double[ vecLength ];
   CGSolveSystem( omega - GSenergy , 1.0 , eta , densityAlphaVector , RealPartSolution , ImagPartSolution );
   if ( TwoRDMreal != NULL ){ Fill2RDM( RealPartSolution , TwoRDMreal ); } // Sets the TwoRDMreal
   RePartGF[0] = FCIddot( vecLength , densityBetaVector , RealPartSolution );
   delete [] RealPartSolution;
   if ( TwoRDMimag != NULL ){ Fill2RDM( ImagPartSolution , TwoRDMimag ); } // Sets the TwoRDMimag
   ImPartGF[0] = FCIddot( vecLength , densityBetaVector , ImagPartSolution );
   delete [] ImagPartSolution;
   
   if ( orb_alpha != orb_beta ){ delete [] densityBetaVector; }
   delete [] densityAlphaVector;

}


