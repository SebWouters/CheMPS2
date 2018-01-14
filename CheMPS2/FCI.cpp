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

#include <climits>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <algorithm>

using std::cout;
using std::endl;

#include "FCI.h"
#include "Irreps.h"
#include "Lapack.h"
#include "Davidson.h"
#include "ConjugateGradient.h"

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
   num_irreps  = Irreps::getNumberOfIrreps( Ham->getNGroup() );
   TargetIrrep = TargetIrrep_in;
   orb2irrep   = new int[ L ];
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
   for ( unsigned int irrep=0; irrep<num_irreps; irrep++ ){
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
   for ( unsigned int irrep = 0; irrep < num_irreps; irrep++ ){
      for ( unsigned int ij = 0; ij < L * L; ij++ ){
         delete [] lookup_cnt_alpha[irrep][ij];
         delete [] lookup_cnt_beta[irrep][ij];
         delete [] lookup_sign_alpha[irrep][ij];
         delete [] lookup_sign_beta[irrep][ij];
      }
      delete [] lookup_cnt_alpha[irrep];
      delete [] lookup_cnt_beta[irrep];
      delete [] lookup_sign_alpha[irrep];
      delete [] lookup_sign_beta[irrep];
   }
   delete [] lookup_cnt_alpha;
   delete [] lookup_cnt_beta;
   delete [] lookup_sign_alpha;
   delete [] lookup_sign_beta;

   // FCI::StartupIrrepCenter
   for ( unsigned int irrep=0; irrep<num_irreps; irrep++ ){
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
   numPerIrrep_up   = new unsigned int[ num_irreps ];
   numPerIrrep_down = new unsigned int[ num_irreps ];
   str2cnt_up       = new int*[ num_irreps ];
   str2cnt_down     = new int*[ num_irreps ];
   cnt2str_up       = new unsigned int*[ num_irreps ];
   cnt2str_down     = new unsigned int*[ num_irreps ];

   for (unsigned int irrep = 0; irrep < num_irreps; irrep++){
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
      int irrep = 0;
      for (unsigned int orb=0; orb<L; orb++){
         if ( bits[orb] ){
            Nparticles++;
            irrep = Irreps::directProd( irrep, getOrb2Irrep( orb ) );
         }
      }
      
      // If allowed: set the corresponding str2cnt to the correct counter and keep track of the number of allowed vectors
      for ( unsigned int irr = 0; irr < num_irreps; irr++ ){
         str2cnt_up  [ irr ][ bitstring ] = -1;
         str2cnt_down[ irr ][ bitstring ] = -1;
      }
      if ( Nparticles == Nel_up ){
         str2cnt_up[ irrep ][ bitstring ] = numPerIrrep_up[ irrep ];
         numPerIrrep_up[ irrep ]++;
      }
      if ( Nparticles == Nel_down ){
         str2cnt_down[ irrep ][ bitstring ] = numPerIrrep_down[ irrep ];
         numPerIrrep_down[ irrep ]++;
      }
   
   }
   
   // Fill the reverse info array: cnt2str
   for ( unsigned int irrep = 0; irrep < num_irreps; irrep++ ){
   
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
   lookup_cnt_alpha  = new int**[ num_irreps ];
   lookup_cnt_beta   = new int**[ num_irreps ];
   lookup_sign_alpha = new int**[ num_irreps ];
   lookup_sign_beta  = new int**[ num_irreps ];

   int * bits = new int[ L ]; // Temporary helper array

   // Quick lookup tables for " sign | new > = E^spinproj_{ij} | old >
   for ( unsigned int irrep = 0; irrep < num_irreps; irrep++ ){

      const unsigned int num_up   = numPerIrrep_up  [ irrep ];
      const unsigned int num_down = numPerIrrep_down[ irrep ];

      lookup_cnt_alpha [ irrep ] = new int*[ L * L ];
      lookup_cnt_beta  [ irrep ] = new int*[ L * L ];
      lookup_sign_alpha[ irrep ] = new int*[ L * L ];
      lookup_sign_beta [ irrep ] = new int*[ L * L ];

      for ( unsigned int ij = 0; ij < L * L; ij++ ){

         lookup_cnt_alpha [ irrep ][ ij ] = new int[ num_up ];
         lookup_cnt_beta  [ irrep ][ ij ] = new int[ num_down ];
         lookup_sign_alpha[ irrep ][ ij ] = new int[ num_up ];
         lookup_sign_beta [ irrep ][ ij ] = new int[ num_down ];

         for ( unsigned int cnt_new_alpha = 0; cnt_new_alpha < num_up; cnt_new_alpha++ ){
            // Check for the sign. If no check for sign, you multiply with sign 0 and everything should be OK...
            lookup_cnt_alpha [ irrep ][ ij ][ cnt_new_alpha ] = 0;
            lookup_sign_alpha[ irrep ][ ij ][ cnt_new_alpha ] = 0;
         }
         for ( unsigned int cnt_new_beta = 0; cnt_new_beta < num_down; cnt_new_beta++ ){
            // Check for the sign. If no check for sign, you multiply with sign and everything should be OK...
            lookup_cnt_beta [ irrep ][ ij ][ cnt_new_beta ] = 0;
            lookup_sign_beta[ irrep ][ ij ][ cnt_new_beta ] = 0;
         }
      }

      for ( unsigned int cnt_new_alpha = 0; cnt_new_alpha < num_up; cnt_new_alpha++ ){

         str2bits( L , cnt2str_up[ irrep ][ cnt_new_alpha ] , bits );

         int phase_crea = 1;
         for ( unsigned int crea = 0; crea < L; crea++ ){
            if ( bits[ crea ] ){
               bits[ crea ] = 0;

               int phase_anni = 1;
               for ( unsigned int anni = 0; anni < L; anni++ ){
                  if ( !(bits[ anni ]) ){
                     bits[ anni ] = 1;

                     const int irrep_old = Irreps::directProd( irrep , Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) ) );
                     const int cnt_old = str2cnt_up[ irrep_old ][ bits2str( L , bits ) ];
                     const int phase = phase_crea * phase_anni;

                     lookup_cnt_alpha [ irrep ][ crea + L * anni ][ cnt_new_alpha ] = cnt_old;
                     lookup_sign_alpha[ irrep ][ crea + L * anni ][ cnt_new_alpha ] = phase;

                     bits[ anni ] = 0;
                  } else {
                     phase_anni *= -1;
                  }
               }

               bits[ crea ] = 1;
               phase_crea *= -1;
            }
         }
      }

      for ( unsigned int cnt_new_beta = 0; cnt_new_beta < num_down; cnt_new_beta++ ){

         str2bits( L , cnt2str_down[ irrep ][ cnt_new_beta ] , bits );

         int phase_crea = 1;
         for ( unsigned int crea = 0; crea < L; crea++ ){
            if ( bits[ crea ] ){
               bits[ crea ] = 0;

               int phase_anni = 1;
               for ( unsigned int anni = 0; anni < L; anni++ ){
                  if ( !(bits[ anni ]) ){
                     bits[ anni ] = 1;

                     const int irrep_old = Irreps::directProd( irrep , Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) ) );
                     const int cnt_old = str2cnt_down[ irrep_old ][ bits2str( L , bits ) ];
                     const int phase = phase_crea * phase_anni;

                     lookup_cnt_beta  [ irrep ][ crea + L * anni ][ cnt_new_beta ] = cnt_old;
                     lookup_sign_beta [ irrep ][ crea + L * anni ][ cnt_new_beta ] = phase;

                     bits[ anni ] = 0;
                  } else {
                     phase_anni *= -1;
                  }
               }

               bits[ crea ] = 1;
               phase_crea *= -1;
            }
         }
      }

   }

   delete [] bits; // Delete temporary helper array

}

void CheMPS2::FCI::StartupIrrepCenter(){

   // Find the orbital combinations which can form a center irrep
   irrep_center_num      = new unsigned int [ num_irreps ];
   irrep_center_crea_orb = new unsigned int*[ num_irreps ];
   irrep_center_anni_orb = new unsigned int*[ num_irreps ];
   
   for ( unsigned int irrep_center = 0; irrep_center < num_irreps; irrep_center++ ){
      const int irrep_center_const_signed = irrep_center;
   
      irrep_center_num[ irrep_center ] = 0;
      for ( unsigned int crea = 0; crea < L; crea++ ){
         for ( unsigned int anni = crea; anni < L; anni++ ){
            if ( Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) ) == irrep_center_const_signed ){
               irrep_center_num[ irrep_center ] += 1;
            }
         }
      }
      irrep_center_crea_orb[ irrep_center ] = new unsigned int[ irrep_center_num[ irrep_center ] ];
      irrep_center_anni_orb[ irrep_center ] = new unsigned int[ irrep_center_num[ irrep_center ] ];
      irrep_center_num[ irrep_center ] = 0;
      for ( unsigned int creator = 0; creator < L; creator++ ){
         for ( unsigned int annihilator = creator; annihilator < L; annihilator++){
            if ( Irreps::directProd( getOrb2Irrep( creator ) , getOrb2Irrep( annihilator ) ) == irrep_center_const_signed ){
               irrep_center_crea_orb[ irrep_center ][ irrep_center_num[ irrep_center ] ] = creator;
               irrep_center_anni_orb[ irrep_center ][ irrep_center_num[ irrep_center ] ] = annihilator;
               irrep_center_num[ irrep_center ] += 1;
            }
         }
      }
   
   }

   irrep_center_jumps = new unsigned int*[ num_irreps ];
   HXVsizeWorkspace = 0;
   for ( unsigned int irrep_center = 0; irrep_center < num_irreps; irrep_center++ ){
      unsigned long long check = 0;
      irrep_center_jumps[ irrep_center ] = new unsigned int[ num_irreps + 1 ];
      const int localTargetIrrep = Irreps::directProd( irrep_center, getTargetIrrep() );
      irrep_center_jumps[ irrep_center ][ 0 ] = 0;
      for ( unsigned int irrep_up = 0; irrep_up < num_irreps; irrep_up++ ){
         const int irrep_down = Irreps::directProd( irrep_up, localTargetIrrep );
         check += ((unsigned long long) numPerIrrep_up[ irrep_up ] ) * ((unsigned long long) numPerIrrep_down[ irrep_down ] );
         unsigned int temp = numPerIrrep_up[ irrep_up ] * numPerIrrep_down[ irrep_down ];
         irrep_center_jumps[ irrep_center ][ irrep_up + 1 ] = irrep_center_jumps[ irrep_center ][ irrep_up ] + temp;
         HXVsizeWorkspace = std::max( HXVsizeWorkspace, ((unsigned long long) irrep_center_num[ irrep_center ] ) * ((unsigned long long) temp ) );
      }
      assert( check <= ((unsigned int) INT_MAX ) ); // Length of FCI vectors should be less then the SIGNED integer size (to be able to call lapack)
   }
   if ( FCIverbose > 0 ){
      cout << "FCI::Startup : Number of variables in the FCI vector = " << getVecLength( 0 ) << endl;
      double num_megabytes = ( 2.0 * sizeof(double) * HXVsizeWorkspace ) / 1048576;
      cout << "FCI::Startup : Without additional loops the FCI matrix-vector product requires a workspace of " << num_megabytes << " MB memory." << endl;
      if ( maxMemWorkMB < num_megabytes ){
         HXVsizeWorkspace = (unsigned int) ceil( ( maxMemWorkMB * 1048576 ) / ( 2 * sizeof(double) ) );
         num_megabytes = ( 2.0 * sizeof(double) * HXVsizeWorkspace ) / 1048576;
         cout << "               For practical purposes, the workspace is constrained to " << num_megabytes << " MB memory." << endl;
      }
   }
   HXVworksmall = new double[ L * L * L * L ];
   HXVworkbig1  = new double[ HXVsizeWorkspace ];
   HXVworkbig2  = new double[ HXVsizeWorkspace ];

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

int CheMPS2::FCI::getUpIrrepOfCounter(const int irrep_center, const unsigned int counter) const{

   int irrep_up = num_irreps;
   while ( counter < irrep_center_jumps[ irrep_center ][ irrep_up-1 ] ){ irrep_up--; }
   return irrep_up-1;
   
}

void CheMPS2::FCI::getBitsOfCounter(const int irrep_center, const unsigned int counter, int * bits_up, int * bits_down) const{

   const int localTargetIrrep = Irreps::directProd( irrep_center, TargetIrrep );
   
   const int irrep_up   = getUpIrrepOfCounter( irrep_center, counter );
   const int irrep_down = Irreps::directProd( irrep_up, localTargetIrrep );
   
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
      if ( bits_up  [ orb ] ){ irrep_up   = Irreps::directProd( irrep_up   , getOrb2Irrep( orb ) ); }
      if ( bits_down[ orb ] ){ irrep_down = Irreps::directProd( irrep_down , getOrb2Irrep( orb ) ); }
   }
   
   const int counter_up   = str2cnt_up  [ irrep_up   ][ string_up   ];
   const int counter_down = str2cnt_down[ irrep_down ][ string_down ];
   
   if (( counter_up == -1 ) || ( counter_down == -1 )){ return 0.0; }
   
   return vector[ irrep_center_jumps[ 0 ][ irrep_up ] + counter_up + numPerIrrep_up[ irrep_up ] * counter_down ];

}

/*void CheMPS2::FCI::CheckHamDEBUG() const{

   const unsigned int vecLength = getVecLength( 0 );
   
   // Building Ham by matvec
   double * HamHXV = new double[ vecLength * vecLength ];
   double * workspace = new double[ vecLength ];
   for (unsigned int count = 0; count < vecLength; count++){
   
      ClearVector( vecLength , workspace );
      workspace[ count ] = 1.0;
      matvec( workspace , HamHXV + count*vecLength );
   
   }
   
   // Building Diag by HamDiag
   DiagHam( workspace );
   double RMSdiagdifference = 0.0;
   for (unsigned int row = 0; row < vecLength; row++){
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
   for (unsigned int row = 0; row < vecLength; row++){
      for (unsigned int col = 0; col < vecLength; col++){
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
   
   // Building Ham^2 by matvec
   double * workspace2 = new double[ vecLength ];
   for (unsigned int count = 0; count < vecLength; count++){
   
      ClearVector( vecLength , workspace );
      workspace[ count ] = 1.0;
      matvec( workspace , workspace2 );
      matvec( workspace2 , HamHXV + count*vecLength );
   
   }
   
   // Building diag( Ham^2 ) by DiagHamSquared
   DiagHamSquared( workspace );
   double RMSdiagdifference2 = 0.0;
   for (unsigned int row = 0; row < vecLength; row++){
      double diff = workspace[ row ] - HamHXV[ row + vecLength * row ];
      RMSdiagdifference2 += diff * diff;
   }
   RMSdiagdifference2 = sqrt( RMSdiagdifference2 );
   cout << "The RMS difference of DiagHamSquared() and diag(HamSquared by HXV) = " << RMSdiagdifference2 << endl;
   
   delete [] workspace2;
   delete [] workspace;
   delete [] HamHXV;

}*/

void CheMPS2::FCI::excite_alpha_omp( const unsigned int dim_new_up, const unsigned int dim_old_up, const unsigned int dim_down, double * origin, double * result, int * signmap, int * countmap ){

   #pragma omp parallel for schedule(static)
   for ( unsigned int cnt_new_up = 0; cnt_new_up < dim_new_up; cnt_new_up++ ){
      const int sign_up = signmap[ cnt_new_up ];
      if ( sign_up != 0 ){
         const int cnt_old_up = countmap[ cnt_new_up ];
         for ( unsigned int cnt_down = 0; cnt_down < dim_down; cnt_down++ ){
            result[ cnt_new_up + dim_new_up * cnt_down ] += sign_up * origin[ cnt_old_up + dim_old_up * cnt_down ];
         }
      }
   }

}

void CheMPS2::FCI::excite_beta_omp( const unsigned int dim_up, const unsigned int dim_new_down, double * origin, double * result, int * signmap, int * countmap ){

   #pragma omp parallel for schedule(static)
   for ( unsigned int cnt_new_down = 0; cnt_new_down < dim_new_down; cnt_new_down++ ){
      const int sign_down = signmap[ cnt_new_down ];
      if ( sign_down != 0 ){
         const int cnt_old_down = countmap[ cnt_new_down ];
         for ( unsigned int cnt_up = 0; cnt_up < dim_up; cnt_up++ ){
            result[ cnt_up + dim_up * cnt_new_down ] += sign_down * origin[ cnt_up + dim_up * cnt_old_down ];
         }
      }
   }

}

void CheMPS2::FCI::excite_alpha_first( const unsigned int dim_new_up, const unsigned int dim_old_up, const unsigned int start_down, const unsigned int stop_down, double * origin, double * result, int * signmap, int * countmap ){

   for ( unsigned int cnt_new_up = 0; cnt_new_up < dim_new_up; cnt_new_up++ ){
      const int sign_up = signmap[ cnt_new_up ];
      if ( sign_up != 0 ){
         const int cnt_old_up = countmap[ cnt_new_up ];
         for ( unsigned int cnt_down = start_down; cnt_down < stop_down; cnt_down++ ){
            result[ cnt_new_up + dim_new_up * ( cnt_down - start_down ) ] += sign_up * origin[ cnt_old_up + dim_old_up * cnt_down ];
         }
      }
   }

}

void CheMPS2::FCI::excite_beta_first( const unsigned int dim_up, const unsigned int start_down, const unsigned int stop_down, double * origin, double * result, int * signmap, int * countmap ){

   for ( unsigned int cnt_new_down = start_down; cnt_new_down < stop_down; cnt_new_down++ ){
      const int sign_down = signmap[ cnt_new_down ];
      if ( sign_down != 0 ){
         const int cnt_old_down = countmap[ cnt_new_down ];
         for ( unsigned int cnt_up = 0; cnt_up < dim_up; cnt_up++ ){
            result[ cnt_up + dim_up * ( cnt_new_down - start_down ) ] += sign_down * origin[ cnt_up + dim_up * cnt_old_down ];
         }
      }
   }

}

void CheMPS2::FCI::excite_alpha_second_omp( const unsigned int dim_new_up, const unsigned int dim_old_up, const unsigned int start_down, const unsigned int stop_down, double * origin, double * result, int * signmap, int * countmap ){

   #pragma omp parallel for schedule(static)
   for ( unsigned int cnt_old_up = 0; cnt_old_up < dim_old_up; cnt_old_up++ ){
      const int sign_up = signmap[ cnt_old_up ];
      if ( sign_up != 0 ){ // Required for thread safety
         const int cnt_new_up = countmap[ cnt_old_up ];
         for ( unsigned int cnt_down = start_down; cnt_down < stop_down; cnt_down++ ){
            result[ cnt_new_up + dim_new_up * cnt_down ] += sign_up * origin[ cnt_old_up + dim_old_up * ( cnt_down - start_down ) ];
         }
      }
   }

}

void CheMPS2::FCI::excite_beta_second_omp( const unsigned int dim_up, const unsigned int start_down, const unsigned int stop_down, double * origin, double * result, int * signmap, int * countmap ){

   #pragma omp parallel for schedule(static)
   for ( unsigned int cnt_old_down = start_down; cnt_old_down < stop_down; cnt_old_down++ ){
      const int sign_down = signmap[ cnt_old_down ];
      if ( sign_down != 0 ){ // Required for thread safety
         const int cnt_new_down = countmap[ cnt_old_down ];
         for ( unsigned int cnt_up = 0; cnt_up < dim_up; cnt_up++ ){
            result[ cnt_up + dim_up * cnt_new_down ] += sign_down * origin[ cnt_up + dim_up * ( cnt_old_down - start_down ) ];
         }
      }
   }

}

void CheMPS2::FCI::matvec( double * input, double * output ) const{

   struct timeval start, end;
   gettimeofday( &start, NULL );

   ClearVector( getVecLength( 0 ), output );

   // P.J. Knowles and N.C. Handy, A new determinant-based full configuration interaction method, Chemical Physics Letters 111 (4-5), 315-321 (1984)

   // irrep_center is the center irrep of the ERI : (ij|kl) --> irrep_center = I_i x I_j = I_k x I_l
   for ( unsigned int irrep_center = 0; irrep_center < num_irreps; irrep_center++ ){

      const int irrep_target_center = Irreps::directProd( TargetIrrep, irrep_center );
      const unsigned int num_pairs  = irrep_center_num[ irrep_center ];
      const unsigned int * center_crea_orb = irrep_center_crea_orb[ irrep_center ];
      const unsigned int * center_anni_orb = irrep_center_anni_orb[ irrep_center ];
      const unsigned int * zero_jumps = irrep_center_jumps[ 0 ];

      for ( unsigned int irrep_center_up = 0; irrep_center_up < num_irreps; irrep_center_up++ ){
         const int irrep_center_down = Irreps::directProd( irrep_target_center, irrep_center_up );
         const unsigned int dim_center_up   = numPerIrrep_up  [ irrep_center_up   ];
         const unsigned int dim_center_down = numPerIrrep_down[ irrep_center_down ];
         if ( dim_center_up * dim_center_down > 0 ){
            const unsigned int blocksize_beta  = HXVsizeWorkspace / std::max( (unsigned int) 1, dim_center_up * num_pairs );
            assert( blocksize_beta > 0 ); // At least one full column should fit in the workspaces...
            unsigned int num_block_beta = dim_center_down / blocksize_beta;
            while ( blocksize_beta * num_block_beta < dim_center_down ){ num_block_beta++; }
            for ( unsigned int block = 0; block < num_block_beta; block++ ){
               const unsigned int start_center_down = block * blocksize_beta;
               const unsigned int  stop_center_down = std::min( ( block + 1 ) * blocksize_beta, dim_center_down );
               const unsigned int size_center = dim_center_up * ( stop_center_down - start_center_down );
               if ( size_center > 0 ){

                  // First build workbig1[ veccounter + size_center * pair ] = E_{i<=j} + ( 1 - delta_i==j ) E_{j>i} (irrep_center) | input >  */
                  #pragma omp parallel for schedule(static)
                  for ( unsigned int pair = 0; pair < num_pairs; pair++ ){
                     double * target_space   = HXVworkbig1 + size_center * pair;
                     const unsigned int crea = center_crea_orb[ pair ];
                     const unsigned int anni = center_anni_orb[ pair ];
                     const int irrep_excited = Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) );
                     const int irrep_zero_up = Irreps::directProd( irrep_excited, irrep_center_up );
                     const unsigned int dim_zero_up = numPerIrrep_up[ irrep_zero_up ];
                     for ( unsigned int count = 0; count < size_center; count++ ){ target_space[ count ] = 0.0; }

                     excite_alpha_first( dim_center_up, dim_zero_up, start_center_down, stop_center_down,
                                         input + zero_jumps[ irrep_zero_up ],
                                         target_space,
                                         lookup_sign_alpha[ irrep_center_up ][ crea + L * anni ],
                                         lookup_cnt_alpha [ irrep_center_up ][ crea + L * anni ] );

                     excite_beta_first( dim_center_up, start_center_down, stop_center_down,
                                        input + zero_jumps[ irrep_center_up ],
                                        target_space,
                                        lookup_sign_beta[ irrep_center_down ][ crea + L * anni ],
                                        lookup_cnt_beta [ irrep_center_down ][ crea + L * anni ] );

                     if ( anni > crea ){

                        excite_alpha_first( dim_center_up, dim_zero_up, start_center_down, stop_center_down,
                                            input + zero_jumps[ irrep_zero_up ],
                                            target_space,
                                            lookup_sign_alpha[ irrep_center_up ][ anni + L * crea ],
                                            lookup_cnt_alpha [ irrep_center_up ][ anni + L * crea ] );

                        excite_beta_first( dim_center_up, start_center_down, stop_center_down,
                                           input + zero_jumps[ irrep_center_up ],
                                           target_space,
                                           lookup_sign_beta[ irrep_center_down ][ anni + L * crea ],
                                           lookup_cnt_beta [ irrep_center_down ][ anni + L * crea ] );

                     }
                  }

                  // If irrep_center == 0, do the one-body terms
                  if ( irrep_center == 0 ){
                     for ( unsigned int pair = 0; pair < num_pairs; pair++ ){
                        HXVworksmall[ pair ] = getGmat( center_crea_orb[ pair ], center_anni_orb[ pair ] );
                     }
                     char notrans = 'N';
                     double one = 1.0;
                     int mdim = size_center;
                     int kdim = num_pairs;
                     int ndim = 1;
                     double * target = output + zero_jumps[ irrep_center_up ] + dim_center_up * start_center_down;
                     dgemm_( &notrans, &notrans, &mdim, &ndim, &kdim, &one, HXVworkbig1, &mdim, HXVworksmall, &kdim, &one, target, &mdim );
                  }

                  // Now build workbig2[ veccounter + size_center * new_pair ] = 0.5 * ( new_pair | old_pair ) * workbig1[ veccounter + size_center * old_pair ]
                  {
                     for ( unsigned int pair1 = 0; pair1 < num_pairs; pair1++ ){
                        for ( unsigned int pair2 = 0; pair2 < num_pairs; pair2++ ){
                           HXVworksmall[ pair1 + num_pairs * pair2 ]
                              = 0.5 * getERI( center_crea_orb[ pair1 ], center_anni_orb[ pair1 ] ,
                                              center_crea_orb[ pair2 ], center_anni_orb[ pair2 ] );
                        }
                     }
                     char notrans = 'N';
                     double one = 1.0;
                     double set = 0.0;
                     int mdim = size_center;
                     int kdim = num_pairs;
                     int ndim = num_pairs;
                     dgemm_( &notrans, &notrans, &mdim, &ndim, &kdim, &one, HXVworkbig1, &mdim, HXVworksmall, &kdim, &set, HXVworkbig2, &mdim );
                  }

                  // Finally do output <-- E_{i<=j} + (1 - delta_{i==j}) E_{j>i} workbig2[ veccounter + size_center * pair ]
                  for ( unsigned int pair = 0; pair < num_pairs; pair++ ){
                     double * origin_space   = HXVworkbig2 + size_center * pair;
                     const unsigned int crea = center_crea_orb[ pair ];
                     const unsigned int anni = center_anni_orb[ pair ];
                     const int irrep_excited = Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) );
                     const int irrep_zero_up = Irreps::directProd( irrep_excited, irrep_center_up );
                     const unsigned int dim_zero_up = numPerIrrep_up[ irrep_zero_up ];

                     excite_alpha_second_omp( dim_zero_up, dim_center_up, start_center_down, stop_center_down,
                                              origin_space,
                                              output + zero_jumps[ irrep_zero_up ],
                                              lookup_sign_alpha[ irrep_center_up ][ anni + L * crea ],
                                              lookup_cnt_alpha [ irrep_center_up ][ anni + L * crea ] );

                     excite_beta_second_omp( dim_center_up, start_center_down, stop_center_down,
                                             origin_space,
                                             output + zero_jumps[ irrep_center_up ],
                                             lookup_sign_beta[ irrep_center_down ][ anni + L * crea ],
                                             lookup_cnt_beta [ irrep_center_down ][ anni + L * crea ] );

                     if ( anni > crea ){

                        excite_alpha_second_omp( dim_zero_up, dim_center_up, start_center_down, stop_center_down,
                                                 origin_space,
                                                 output + zero_jumps[ irrep_zero_up ],
                                                 lookup_sign_alpha[ irrep_center_up ][ crea + L * anni ],
                                                 lookup_cnt_alpha [ irrep_center_up ][ crea + L * anni ] );

                        excite_beta_second_omp( dim_center_up, start_center_down, stop_center_down,
                                                origin_space,
                                                output + zero_jumps[ irrep_center_up ],
                                                lookup_sign_beta[ irrep_center_down ][ crea + L * anni ],
                                                lookup_cnt_beta [ irrep_center_down ][ crea + L * anni ] );

                     }
                  }
               }
            }
         }
      }
   }

   gettimeofday( &end, NULL );
   const double elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   if ( FCIverbose >= 1 ){ cout << "FCI::matvec : Wall time = " << elapsed << " seconds" << endl; }

}

void CheMPS2::FCI::apply_excitation( double * orig_vector, double * result_vector, const int crea, const int anni, const int orig_target_irrep ) const{

   const int    excitation_irrep = Irreps::directProd( getOrb2Irrep( crea ), getOrb2Irrep( anni ) );
   const int result_target_irrep = Irreps::directProd( excitation_irrep, orig_target_irrep );
   const int   orig_irrep_center = Irreps::directProd( TargetIrrep,   orig_target_irrep );
   const int result_irrep_center = Irreps::directProd( TargetIrrep, result_target_irrep );

   ClearVector( getVecLength( result_irrep_center ) , result_vector );

   for ( unsigned int result_irrep_up = 0; result_irrep_up < num_irreps; result_irrep_up++ ){

      const int result_irrep_down = Irreps::directProd( result_irrep_up, result_target_irrep );
      const int orig_irrep_up     = Irreps::directProd( excitation_irrep, result_irrep_up );

      excite_alpha_omp( numPerIrrep_up  [ result_irrep_up   ], // dim_new_up
                        numPerIrrep_up  [   orig_irrep_up   ], // dim_old_up
                        numPerIrrep_down[ result_irrep_down ], // dim_down
                        orig_vector   + irrep_center_jumps[   orig_irrep_center ][   orig_irrep_up ], // origin
                        result_vector + irrep_center_jumps[ result_irrep_center ][ result_irrep_up ], // result
                        lookup_sign_alpha[ result_irrep_up ][ crea + L * anni ],   // signmap
                        lookup_cnt_alpha [ result_irrep_up ][ crea + L * anni ] ); // countmap

      excite_beta_omp( numPerIrrep_up  [ result_irrep_up   ], // dim_up
                       numPerIrrep_down[ result_irrep_down ], // dim_new_down
                       orig_vector   + irrep_center_jumps[   orig_irrep_center ][ result_irrep_up ], // origin
                       result_vector + irrep_center_jumps[ result_irrep_center ][ result_irrep_up ], // result
                       lookup_sign_beta[ result_irrep_down ][ crea + L * anni ],   // signmap
                       lookup_cnt_beta [ result_irrep_down ][ crea + L * anni ] ); // countmap

   }

}

double CheMPS2::FCI::Fill2RDM(double * vector, double * two_rdm) const{

   assert( Nel_up + Nel_down >= 2 );

   struct timeval start, end;
   gettimeofday(&start, NULL);
   
   ClearVector( L*L*L*L, two_rdm );
   const unsigned int orig_length = getVecLength( 0 );
   unsigned int max_length = 0;
   for ( unsigned int irrep = 0; irrep < num_irreps; irrep++ ){
      if ( getVecLength( irrep ) > max_length ){ max_length = getVecLength( irrep ); }
   }
   double * workspace1 = new double[ max_length  ];
   double * workspace2 = new double[ orig_length ];
   
   // Gamma_{ijkl} = < E_ik E_jl > - delta_jk < E_il >
   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){ // anni1 = l
      for ( unsigned int crea1 = anni1; crea1 < L; crea1++ ){ // crea1 = j >= l
      
         const int irrep_center1 = Irreps::directProd( getOrb2Irrep( crea1 ), getOrb2Irrep( anni1 ) );
         const int target_irrep1 = Irreps::directProd( TargetIrrep, irrep_center1 );
         apply_excitation( vector, workspace1, crea1, anni1, TargetIrrep );
         
         if ( irrep_center1 == 0 ){
            const double value = FCIddot( orig_length, workspace1, vector ); // < E_{crea1,anni1} >
            for ( unsigned int jk = anni1; jk < L; jk++ ){
               two_rdm[ crea1 + L * ( jk + L * ( jk + L * anni1 ) ) ] -= value;
            }
         }
         
         for ( unsigned int crea2 = anni1; crea2 < L; crea2++ ){ // crea2 = i >= l
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){ // anni2 = k >= l
            
               const int irrep_center2 = Irreps::directProd( getOrb2Irrep( crea2 ), getOrb2Irrep( anni2 ) );
               if ( irrep_center2 == irrep_center1 ){
               
                  apply_excitation( workspace1, workspace2, crea2, anni2, target_irrep1 );
                  const double value = FCIddot( orig_length, workspace2, vector ); // < E_{crea2,anni2} E_{crea1,anni1} >
                  two_rdm[ crea2 + L * ( crea1 + L * ( anni2 + L * anni1 ) ) ] += value;
                  
               }
            }
         }
      }
   }
   delete [] workspace1;
   delete [] workspace2;
   
   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){
      for ( unsigned int crea1 = anni1; crea1 < L; crea1++ ){
         const int irrep_center1 = Irreps::directProd( getOrb2Irrep( crea1 ) , getOrb2Irrep( anni1 ) );
         for ( unsigned int crea2 = anni1; crea2 < L; crea2++ ){
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){
               const int irrep_center2 = Irreps::directProd( getOrb2Irrep( crea2 ) , getOrb2Irrep( anni2 ) );
               if ( irrep_center2 == irrep_center1 ){
                  const double value = two_rdm[ crea2 + L * ( crea1 + L * ( anni2 + L * anni1 ) ) ];
                                       two_rdm[ crea1 + L * ( crea2 + L * ( anni1 + L * anni2 ) ) ] = value;
                                       two_rdm[ anni2 + L * ( anni1 + L * ( crea2 + L * crea1 ) ) ] = value;
                                       two_rdm[ anni1 + L * ( anni2 + L * ( crea1 + L * crea2 ) ) ] = value;
               }
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
            tempvar2 += two_rdm[ orb1 + L * ( orb3 + L * ( orb2 + L * orb3 ) ) ];
            for ( unsigned int orb4 = 0; orb4 < L; orb4++ ){
               FCIenergy += 0.5 * two_rdm[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ] * getERI( orb1 , orb3 , orb2 , orb4 );
            }
         }
         FCIenergy += ( getGmat( orb1 , orb2 ) + 0.5 * tempvar ) * tempvar2 / ( Nel_up + Nel_down - 1.0);
      }
   }
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   if ( FCIverbose > 0 ){ cout << "FCI::Fill2RDM : Wall time = " << elapsed << " seconds" << endl; }
   if ( FCIverbose > 0 ){ cout << "FCI::Fill2RDM : Energy (Ham * 2-RDM) = " << FCIenergy << endl; }
   return FCIenergy;

}

void CheMPS2::FCI::Fill4RDM(double * vector, double * four_rdm) const{

   assert( Nel_up + Nel_down >= 4 );

   struct timeval start, end;
   gettimeofday(&start, NULL);
   
   /*
      Gamma_{ijkl,pqrt} = < E_ip E_jq E_kr E_lt >
                        - delta_lr < E_ip E_jq E_kt >
                        - delta_lq < E_ip E_kr E_jt >
                        - delta_lp < E_jq E_kr E_it >
                        - delta_kq < E_ip E_jr E_lt >
                        - delta_kp < E_ir E_jq E_lt >
                        - delta_jp < E_iq E_kr E_lt >
                        + delta_lr delta_kq < E_ip E_jt >
                        + delta_lr delta_kp < E_jq E_it >
                        + delta_lr delta_jp < E_iq E_kt >
                        + delta_lq delta_jr < E_ip E_kt >
                        + delta_lq delta_kp < E_ir E_jt >
                        + delta_lq delta_jp < E_kr E_it >
                        + delta_lp delta_kq < E_jr E_it >
                        + delta_lp delta_iq < E_kr E_jt >
                        + delta_lp delta_ir < E_jq E_kt >
                        + delta_kq delta_jp < E_ir E_lt >
                        + delta_kp delta_jr < E_iq E_lt >
                        - delta_jp delta_kq delta_lr < E_it >
                        - delta_kp delta_lq delta_jr < E_it >
                        - delta_kp delta_iq delta_lr < E_jt >
                        - delta_lp delta_kq delta_ir < E_jt >
                        - delta_jp delta_lq delta_ir < E_kt >
                        - delta_lp delta_iq delta_jr < E_kt >
   */
   
   ClearVector( L*L*L*L*L*L*L*L, four_rdm );
   const unsigned int orig_length = getVecLength( 0 );
   unsigned int max_length = getVecLength( 0 );
   for ( unsigned int irrep = 1; irrep < num_irreps; irrep++ ){
      if ( getVecLength( irrep ) > max_length ){ max_length = getVecLength( irrep ); }
   }
   double * workspace1 = new double[ max_length  ];
   double * workspace2 = new double[ max_length  ];
   double * workspace3 = new double[ max_length  ];
   double * workspace4 = new double[ orig_length ];

   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){ // anni1 = t
      for ( unsigned int crea1 = anni1; crea1 < L; crea1++ ){ // crea1 = l >= t

         const int irrep_center1 = Irreps::directProd( getOrb2Irrep( crea1 ), getOrb2Irrep( anni1 ) );
         const int target_irrep1 = Irreps::directProd( TargetIrrep, irrep_center1 );
         apply_excitation( vector, workspace1, crea1, anni1, TargetIrrep );

         if ( irrep_center1 == 0 ){

            // value = < E_{crea1,anni1} >
            const double value = FCIddot( orig_length, workspace1, vector );
            for ( unsigned int p = anni1; p < L; p++ ){
               for ( unsigned int q = anni1; q < L; q++ ){
                  for ( unsigned int r = anni1; r < L; r++ ){
                     //        i           j           k           l       p       q       r       t
                     four_rdm[ crea1 + L*( p     + L*( q     + L*( r + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_jp delta_kq delta_lr < E_it >
                     four_rdm[ crea1 + L*( r     + L*( p     + L*( q + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_kp delta_lq delta_jr < E_it >
                     four_rdm[ q     + L*( crea1 + L*( p     + L*( r + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_kp delta_iq delta_lr < E_jt >
                     four_rdm[ r     + L*( crea1 + L*( q     + L*( p + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_lp delta_kq delta_ir < E_jt >
                     four_rdm[ r     + L*( p     + L*( crea1 + L*( q + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_jp delta_lq delta_ir < E_kt >
                     four_rdm[ q     + L*( r     + L*( crea1 + L*( p + L*( p + L*( q + L*( r + L * anni1 ))))))] -= value; // - delta_lp delta_iq delta_jr < E_kt >
                  }
               }
            }

         }

         for ( unsigned int crea2 = anni1; crea2 < L; crea2++ ){ // crea2 = k >= t
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){ // anni2 = r >= t

               const int irrep_center2 = Irreps::directProd( getOrb2Irrep( crea2 ), getOrb2Irrep( anni2 ) );
               const int target_irrep2 = Irreps::directProd( target_irrep1, irrep_center2 );
               apply_excitation( workspace1, workspace2, crea2, anni2, target_irrep1 );

               if ( irrep_center1 == irrep_center2 ){

                  // value = < E_{crea2,anni2} E_{crea1,anni1} >
                  const double value = FCIddot( orig_length, workspace2, vector );
                  for ( unsigned int orb = anni1; orb < L; orb++ ){
                     for ( unsigned int xyz = anni1; xyz < L; xyz++ ){
                        //        i           j           k           l           p           q           r           t
                        four_rdm[ crea2 + L*( crea1 + L*( xyz   + L*( orb   + L*( anni2 + L*( xyz   + L*( orb   + L * anni1 ))))))] += value; // + delta_lr delta_kq < E_ip E_jt >
                        four_rdm[ crea1 + L*( crea2 + L*( xyz   + L*( orb   + L*( xyz   + L*( anni2 + L*( orb   + L * anni1 ))))))] += value; // + delta_lr delta_kp < E_jq E_it >
                        four_rdm[ crea2 + L*( xyz   + L*( crea1 + L*( orb   + L*( xyz   + L*( anni2 + L*( orb   + L * anni1 ))))))] += value; // + delta_lr delta_jp < E_iq E_kt >
                        four_rdm[ crea2 + L*( xyz   + L*( crea1 + L*( orb   + L*( anni2 + L*( orb   + L*( xyz   + L * anni1 ))))))] += value; // + delta_lq delta_jr < E_ip E_kt >
                        four_rdm[ crea2 + L*( crea1 + L*( xyz   + L*( orb   + L*( xyz   + L*( orb   + L*( anni2 + L * anni1 ))))))] += value; // + delta_lq delta_kp < E_ir E_jt >
                        four_rdm[ crea1 + L*( xyz   + L*( crea2 + L*( orb   + L*( xyz   + L*( orb   + L*( anni2 + L * anni1 ))))))] += value; // + delta_lq delta_jp < E_kr E_it >
                        four_rdm[ crea1 + L*( crea2 + L*( xyz   + L*( orb   + L*( orb   + L*( xyz   + L*( anni2 + L * anni1 ))))))] += value; // + delta_lp delta_kq < E_jr E_it >
                        four_rdm[ xyz   + L*( crea1 + L*( crea2 + L*( orb   + L*( orb   + L*( xyz   + L*( anni2 + L * anni1 ))))))] += value; // + delta_lp delta_iq < E_kr E_jt >
                        four_rdm[ xyz   + L*( crea2 + L*( crea1 + L*( orb   + L*( orb   + L*( anni2 + L*( xyz   + L * anni1 ))))))] += value; // + delta_lp delta_ir < E_jq E_kt >
                        four_rdm[ crea2 + L*( xyz   + L*( orb   + L*( crea1 + L*( xyz   + L*( orb   + L*( anni2 + L * anni1 ))))))] += value; // + delta_kq delta_jp < E_ir E_lt >
                        four_rdm[ crea2 + L*( xyz   + L*( orb   + L*( crea1 + L*( orb   + L*( anni2 + L*( xyz   + L * anni1 ))))))] += value; // + delta_kp delta_jr < E_iq E_lt >
                     }
                  }

               }

               for ( unsigned int crea3 = anni1; crea3 < L; crea3++ ){ // crea3 = j >= t = anni1
                  for ( unsigned int anni3 = (( anni1 == anni2 ) ? anni1 + 1 : anni1 ); anni3 < L; anni3++ ){ // anni3 = q >= t

                     const int irrep_center3 = Irreps::directProd( getOrb2Irrep( crea3 ), getOrb2Irrep( anni3 ) );
                     const int target_irrep3 = Irreps::directProd( target_irrep2, irrep_center3 );
                     const int irrep_center4 = Irreps::directProd( Irreps::directProd( irrep_center1 , irrep_center2 ), irrep_center3 );
                     apply_excitation( workspace2, workspace3, crea3, anni3, target_irrep2 );

                     if ( irrep_center4 == 0 ){

                        // value = < E_{crea3,anni3} E_{crea2,anni2} E_{crea1,anni1} >
                        const double value = FCIddot( orig_length, workspace3, vector );
                        for ( unsigned int orb = anni1; orb < L; orb++ ){
                           //        i           j           k           l           p           q           r           t
                           four_rdm[ crea3 + L*( crea2 + L*( crea1 + L*( orb   + L*( anni3 + L*( anni2 + L*( orb   + L * anni1 ))))))] -= value; // - delta_lr < E_ip E_jq E_kt >
                           four_rdm[ crea3 + L*( crea1 + L*( crea2 + L*( orb   + L*( anni3 + L*( orb   + L*( anni2 + L * anni1 ))))))] -= value; // - delta_lq < E_ip E_kr E_jt >
                           four_rdm[ crea1 + L*( crea3 + L*( crea2 + L*( orb   + L*( orb   + L*( anni3 + L*( anni2 + L * anni1 ))))))] -= value; // - delta_lp < E_jq E_kr E_it >
                           four_rdm[ crea3 + L*( crea2 + L*( orb   + L*( crea1 + L*( anni3 + L*( orb   + L*( anni2 + L * anni1 ))))))] -= value; // - delta_kq < E_ip E_jr E_lt >
                           four_rdm[ crea3 + L*( crea2 + L*( orb   + L*( crea1 + L*( orb   + L*( anni2 + L*( anni3 + L * anni1 ))))))] -= value; // - delta_kp < E_ir E_jq E_lt >
                           four_rdm[ crea3 + L*( orb   + L*( crea2 + L*( crea1 + L*( orb   + L*( anni3 + L*( anni2 + L * anni1 ))))))] -= value; // - delta_jp < E_iq E_kr E_lt >
                        }

                     }

                     if ( crea3 >= (( crea1 == crea2 ) ? crea2 + 1 : crea2 ) ){ // crea3 = j >= k = crea2 >= t = anni1

                        for ( unsigned int crea4 = (( crea2 == crea3 ) ? crea3 + 1 : crea3 ); crea4 < L; crea4++ ){ // crea4 = i >= j = crea3 >= k = crea2 >= t = anni1
                           for ( unsigned int anni4 = (( anni1 == anni2 ) ? anni1 + 1 : anni1 ); anni4 < L; anni4++ ){ // anni4 = p >= t

                              if ( (( anni2 == anni3 ) && ( anni3 == anni4 )) == false ){

                                 const int irrep_product4 = Irreps::directProd( getOrb2Irrep( crea4 ), getOrb2Irrep( anni4 ) );
                                 if ( irrep_product4 == irrep_center4 ){

                                    apply_excitation( workspace3, workspace4, crea4, anni4, target_irrep3 );
                                    // value = < E_{crea4,anni4} E_{crea3,anni3} E_{crea2,anni2} E_{crea1,anni1} >
                                    const double value = FCIddot( orig_length, workspace4, vector );
                                    four_rdm[ crea4 + L*( crea3 + L*( crea2 + L*( crea1 + L*( anni4 + L*( anni3 + L*( anni2 + L * anni1 ))))))] += value;

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
   delete [] workspace1;
   delete [] workspace2;
   delete [] workspace3;
   delete [] workspace4;
   
   // Make 48-fold permutation symmetric
   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){ // anni1 = t
      for ( unsigned int crea1 = anni1; crea1 < L; crea1++ ){ // crea1 = l >= t = anni1
      const int irrep_center1 = Irreps::directProd( getOrb2Irrep( crea1 ), getOrb2Irrep( anni1 ) );
         for ( unsigned int crea2 = anni1; crea2 < L; crea2++ ){ // crea2 = k >= t = anni1
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){ // anni2 = r >= t = anni1
               const int irrep_center2 = Irreps::directProd( getOrb2Irrep( crea2 ), getOrb2Irrep( anni2 ) );
               const int irrep_12 = Irreps::directProd( irrep_center1, irrep_center2 );
               for ( unsigned int crea3 = crea2; crea3 < L; crea3++ ){ // crea3 = j >= k = crea2 >= t = anni1
                  for ( unsigned int anni3 = anni1; anni3 < L; anni3++ ){ // anni3 = q >= t = anni1
                     const int irrep_center3 = Irreps::directProd( getOrb2Irrep( crea3 ), getOrb2Irrep( anni3 ) );
                     const int irrep_123 = Irreps::directProd( irrep_12, irrep_center3 );
                     for ( unsigned int crea4 = crea3; crea4 < L; crea4++ ){ // crea4 = i >= j = crea3 >= k = crea2 >= t = anni1
                        for ( unsigned int anni4 = anni1; anni4 < L; anni4++ ){ // anni4 = p >= t = anni1
                           const int irrep_center4 = Irreps::directProd( getOrb2Irrep( crea4 ), getOrb2Irrep( anni4 ) );
                           if ( irrep_123 == irrep_center4 ){

                              /* crea4 >= crea3 >= crea2 >= anni1
                                 crea1, anni2, anni3, anni4 >= anni1 */

         const double value = four_rdm[ crea4 + L*( crea3 + L*( crea2 + L*( crea1 + L*( anni4 + L*( anni3 + L*( anni2 + L * anni1 )))))) ];
                              four_rdm[ crea4 + L*( crea2 + L*( crea1 + L*( crea3 + L*( anni4 + L*( anni2 + L*( anni1 + L * anni3 )))))) ] = value;
                              four_rdm[ crea4 + L*( crea1 + L*( crea3 + L*( crea2 + L*( anni4 + L*( anni1 + L*( anni3 + L * anni2 )))))) ] = value;
                              four_rdm[ crea4 + L*( crea1 + L*( crea2 + L*( crea3 + L*( anni4 + L*( anni1 + L*( anni2 + L * anni3 )))))) ] = value;
                              four_rdm[ crea4 + L*( crea2 + L*( crea3 + L*( crea1 + L*( anni4 + L*( anni2 + L*( anni3 + L * anni1 )))))) ] = value;
                              four_rdm[ crea4 + L*( crea3 + L*( crea1 + L*( crea2 + L*( anni4 + L*( anni3 + L*( anni1 + L * anni2 )))))) ] = value;
                              
                              four_rdm[ crea3 + L*( crea4 + L*( crea2 + L*( crea1 + L*( anni3 + L*( anni4 + L*( anni2 + L * anni1 )))))) ] = value;
                              four_rdm[ crea3 + L*( crea2 + L*( crea1 + L*( crea4 + L*( anni3 + L*( anni2 + L*( anni1 + L * anni4 )))))) ] = value;
                              four_rdm[ crea3 + L*( crea1 + L*( crea4 + L*( crea2 + L*( anni3 + L*( anni1 + L*( anni4 + L * anni2 )))))) ] = value;
                              four_rdm[ crea3 + L*( crea1 + L*( crea2 + L*( crea4 + L*( anni3 + L*( anni1 + L*( anni2 + L * anni4 )))))) ] = value;
                              four_rdm[ crea3 + L*( crea2 + L*( crea4 + L*( crea1 + L*( anni3 + L*( anni2 + L*( anni4 + L * anni1 )))))) ] = value;
                              four_rdm[ crea3 + L*( crea4 + L*( crea1 + L*( crea2 + L*( anni3 + L*( anni4 + L*( anni1 + L * anni2 )))))) ] = value;
                              
                              four_rdm[ crea2 + L*( crea3 + L*( crea4 + L*( crea1 + L*( anni2 + L*( anni3 + L*( anni4 + L * anni1 )))))) ] = value;
                              four_rdm[ crea2 + L*( crea4 + L*( crea1 + L*( crea3 + L*( anni2 + L*( anni4 + L*( anni1 + L * anni3 )))))) ] = value;
                              four_rdm[ crea2 + L*( crea1 + L*( crea3 + L*( crea4 + L*( anni2 + L*( anni1 + L*( anni3 + L * anni4 )))))) ] = value;
                              four_rdm[ crea2 + L*( crea1 + L*( crea4 + L*( crea3 + L*( anni2 + L*( anni1 + L*( anni4 + L * anni3 )))))) ] = value;
                              four_rdm[ crea2 + L*( crea4 + L*( crea3 + L*( crea1 + L*( anni2 + L*( anni4 + L*( anni3 + L * anni1 )))))) ] = value;
                              four_rdm[ crea2 + L*( crea3 + L*( crea1 + L*( crea4 + L*( anni2 + L*( anni3 + L*( anni1 + L * anni4 )))))) ] = value;
                              
                              four_rdm[ crea1 + L*( crea3 + L*( crea2 + L*( crea4 + L*( anni1 + L*( anni3 + L*( anni2 + L * anni4 )))))) ] = value;
                              four_rdm[ crea1 + L*( crea2 + L*( crea4 + L*( crea3 + L*( anni1 + L*( anni2 + L*( anni4 + L * anni3 )))))) ] = value;
                              four_rdm[ crea1 + L*( crea4 + L*( crea3 + L*( crea2 + L*( anni1 + L*( anni4 + L*( anni3 + L * anni2 )))))) ] = value;
                              four_rdm[ crea1 + L*( crea4 + L*( crea2 + L*( crea3 + L*( anni1 + L*( anni4 + L*( anni2 + L * anni3 )))))) ] = value;
                              four_rdm[ crea1 + L*( crea2 + L*( crea3 + L*( crea4 + L*( anni1 + L*( anni2 + L*( anni3 + L * anni4 )))))) ] = value;
                              four_rdm[ crea1 + L*( crea3 + L*( crea4 + L*( crea2 + L*( anni1 + L*( anni3 + L*( anni4 + L * anni2 )))))) ] = value;
                              
                              four_rdm[ anni4 + L*( anni3 + L*( anni2 + L*( anni1 + L*( crea4 + L*( crea3 + L*( crea2 + L * crea1 )))))) ] = value;
                              four_rdm[ anni4 + L*( anni2 + L*( anni1 + L*( anni3 + L*( crea4 + L*( crea2 + L*( crea1 + L * crea3 )))))) ] = value;
                              four_rdm[ anni4 + L*( anni1 + L*( anni3 + L*( anni2 + L*( crea4 + L*( crea1 + L*( crea3 + L * crea2 )))))) ] = value;
                              four_rdm[ anni4 + L*( anni1 + L*( anni2 + L*( anni3 + L*( crea4 + L*( crea1 + L*( crea2 + L * crea3 )))))) ] = value;
                              four_rdm[ anni4 + L*( anni2 + L*( anni3 + L*( anni1 + L*( crea4 + L*( crea2 + L*( crea3 + L * crea1 )))))) ] = value;
                              four_rdm[ anni4 + L*( anni3 + L*( anni1 + L*( anni2 + L*( crea4 + L*( crea3 + L*( crea1 + L * crea2 )))))) ] = value;
                              
                              four_rdm[ anni3 + L*( anni4 + L*( anni2 + L*( anni1 + L*( crea3 + L*( crea4 + L*( crea2 + L * crea1 )))))) ] = value;
                              four_rdm[ anni3 + L*( anni2 + L*( anni1 + L*( anni4 + L*( crea3 + L*( crea2 + L*( crea1 + L * crea4 )))))) ] = value;
                              four_rdm[ anni3 + L*( anni1 + L*( anni4 + L*( anni2 + L*( crea3 + L*( crea1 + L*( crea4 + L * crea2 )))))) ] = value;
                              four_rdm[ anni3 + L*( anni1 + L*( anni2 + L*( anni4 + L*( crea3 + L*( crea1 + L*( crea2 + L * crea4 )))))) ] = value;
                              four_rdm[ anni3 + L*( anni2 + L*( anni4 + L*( anni1 + L*( crea3 + L*( crea2 + L*( crea4 + L * crea1 )))))) ] = value;
                              four_rdm[ anni3 + L*( anni4 + L*( anni1 + L*( anni2 + L*( crea3 + L*( crea4 + L*( crea1 + L * crea2 )))))) ] = value;
                              
                              four_rdm[ anni2 + L*( anni3 + L*( anni4 + L*( anni1 + L*( crea2 + L*( crea3 + L*( crea4 + L * crea1 )))))) ] = value;
                              four_rdm[ anni2 + L*( anni4 + L*( anni1 + L*( anni3 + L*( crea2 + L*( crea4 + L*( crea1 + L * crea3 )))))) ] = value;
                              four_rdm[ anni2 + L*( anni1 + L*( anni3 + L*( anni4 + L*( crea2 + L*( crea1 + L*( crea3 + L * crea4 )))))) ] = value;
                              four_rdm[ anni2 + L*( anni1 + L*( anni4 + L*( anni3 + L*( crea2 + L*( crea1 + L*( crea4 + L * crea3 )))))) ] = value;
                              four_rdm[ anni2 + L*( anni4 + L*( anni3 + L*( anni1 + L*( crea2 + L*( crea4 + L*( crea3 + L * crea1 )))))) ] = value;
                              four_rdm[ anni2 + L*( anni3 + L*( anni1 + L*( anni4 + L*( crea2 + L*( crea3 + L*( crea1 + L * crea4 )))))) ] = value;
                              
                              four_rdm[ anni1 + L*( anni3 + L*( anni2 + L*( anni4 + L*( crea1 + L*( crea3 + L*( crea2 + L * crea4 )))))) ] = value;
                              four_rdm[ anni1 + L*( anni2 + L*( anni4 + L*( anni3 + L*( crea1 + L*( crea2 + L*( crea4 + L * crea3 )))))) ] = value;
                              four_rdm[ anni1 + L*( anni4 + L*( anni3 + L*( anni2 + L*( crea1 + L*( crea4 + L*( crea3 + L * crea2 )))))) ] = value;
                              four_rdm[ anni1 + L*( anni4 + L*( anni2 + L*( anni3 + L*( crea1 + L*( crea4 + L*( crea2 + L * crea3 )))))) ] = value;
                              four_rdm[ anni1 + L*( anni2 + L*( anni3 + L*( anni4 + L*( crea1 + L*( crea2 + L*( crea3 + L * crea4 )))))) ] = value;
                              four_rdm[ anni1 + L*( anni3 + L*( anni4 + L*( anni2 + L*( crea1 + L*( crea3 + L*( crea4 + L * crea2 )))))) ] = value;

                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   for ( unsigned int ann = 0; ann < L; ann++ ){
      for ( unsigned int orb = 0; orb < L; orb++ ){
         for ( unsigned int combo = 0; combo < L*L*L*L; combo++ ){
            four_rdm[ combo + L*L*L*L*( ann + L * ( ann + L * ( ann + L * orb ))) ] = 0.0;
            four_rdm[ combo + L*L*L*L*( ann + L * ( ann + L * ( orb + L * ann ))) ] = 0.0;
            four_rdm[ combo + L*L*L*L*( ann + L * ( orb + L * ( ann + L * ann ))) ] = 0.0;
            four_rdm[ combo + L*L*L*L*( orb + L * ( ann + L * ( ann + L * ann ))) ] = 0.0;
            four_rdm[ L*L*L*L * combo + ann + L * ( ann + L * ( ann + L * orb ))  ] = 0.0;
            four_rdm[ L*L*L*L * combo + ann + L * ( ann + L * ( orb + L * ann ))  ] = 0.0;
            four_rdm[ L*L*L*L * combo + ann + L * ( orb + L * ( ann + L * ann ))  ] = 0.0;
            four_rdm[ L*L*L*L * combo + orb + L * ( ann + L * ( ann + L * ann ))  ] = 0.0;
         }
      }
   }
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   if ( FCIverbose > 0 ){ cout << "FCI::Fill4RDM : Wall time = " << elapsed << " seconds" << endl; }

}

void CheMPS2::FCI::Fock4RDM( double * vector, double * three_rdm, double * fock, double * output ) const{

   assert( Nel_up + Nel_down >= 4 );
   const double elapsed = Driver3RDM( vector, output, three_rdm, fock, L + 1 );
   if ( FCIverbose > 0 ){ cout << "FCI::Fock4RDM : Wall time = " << elapsed << " seconds" << endl; }

}

void CheMPS2::FCI::Fill3RDM( double * vector, double * output ) const{

   assert( Nel_up + Nel_down >= 3 );
   const double elapsed = Driver3RDM( vector, output, NULL, NULL, L + 1 );
   if ( FCIverbose > 0 ){ cout << "FCI::Fill3RDM : Wall time = " << elapsed << " seconds" << endl; }

}

void CheMPS2::FCI::Diag4RDM( double * vector, double * three_rdm, const unsigned int orbz, double * output ) const{

   assert( Nel_up + Nel_down >= 4 );
   const double elapsed = Driver3RDM( vector, output, three_rdm, NULL, orbz );
   if ( FCIverbose > 0 ){ cout << "FCI::Diag4RDM : Wall time = " << elapsed << " seconds" << endl; }

}

double CheMPS2::FCI::Driver3RDM( double * vector, double * output, double * three_rdm, double * fock, const unsigned int orbz ) const{

   struct timeval start, end;
   gettimeofday(&start, NULL);

   ClearVector( L*L*L*L*L*L, output );
   const unsigned int orig_length = getVecLength( 0 );
   unsigned int max_length = getVecLength( 0 );
   for ( unsigned int irrep = 1; irrep < num_irreps; irrep++ ){
      if ( getVecLength( irrep ) > max_length ){ max_length = getVecLength( irrep ); }
   }
   double * workspace1 = new double[ max_length  ];
   double * workspace2 = new double[ max_length  ];
   double * workspace3 = new double[ orig_length ];

   double * chi = NULL;
   const bool task_fock =  ( fock != NULL );
   const bool task_E_zz = (( fock == NULL ) && ( three_rdm != NULL ));
   const bool task_3rdm = (( fock == NULL ) && ( three_rdm == NULL ));

   if ( task_3rdm ){

      assert( orbz == L + 1 );
      chi = vector; // | Chi > = | 0 >

      /* Calculate the 3-RDM:
            Gamma_{ijk,pqr} = < 0 | E_ip E_jq E_kr | Chi >
                            - delta_kq < 0 | E_ip E_jr | Chi >
                            - delta_kp < 0 | E_ir E_jq | Chi >
                            - delta_jp < 0 | E_iq E_kr | Chi >
                            + delta_kq delta_jp < 0 | E_ir | Chi >
                            + delta_kp delta_jr < 0 | E_iq | Chi >
      */
   }

   if ( task_E_zz ){

      assert( orbz < L );
      chi = new double[ orig_length ];
      apply_excitation( vector, chi, orbz, orbz, TargetIrrep ); // | Chi > = E_zz | 0 >

      /* Calculate the 4-RDM elements with orbital z fixed:
            Gamma_{ijkz,pqrz} = < 0 | E_ip E_jq E_kr | Chi >
                              - delta_kq < 0 | E_ip E_jr | Chi >
                              - delta_kp < 0 | E_ir E_jq | Chi >
                              - delta_jp < 0 | E_iq E_kr | Chi >
                              + delta_kq delta_jp < 0 | E_ir | Chi >
                              + delta_kp delta_jr < 0 | E_iq | Chi >
                              - ( delta_pz + delta_qz + delta_rz ) Gamma_{ijk,pqr}
      */
   }

   if ( task_fock ){

      assert( orbz == L + 1 );
      chi = new double[ orig_length ];
      ClearVector( orig_length, chi );
      for ( unsigned int anni = 0; anni < L; anni++ ){
         for ( unsigned int crea = 0; crea < L; crea++ ){
            if ( getOrb2Irrep( crea ) == getOrb2Irrep( anni ) ){
               apply_excitation( vector, workspace1, crea, anni, TargetIrrep );
               FCIdaxpy( orig_length, fock[ crea + L * anni ], workspace1, chi ); // | Chi > = sum_{l,t} fock[l,t] E_lt | 0 >
            }
         }
      }

      /* Calculate the 4-RDM contracted with the Fock operator:
            sum_{l,t} fock[l,t] Gamma_{ijkl,pqrt} = < 0 | E_ip E_jq E_kr | Chi >
                                                  - delta_kq < 0 | E_ip E_jr | Chi >
                                                  - delta_kp < 0 | E_ir E_jq | Chi >
                                                  - delta_jp < 0 | E_iq E_kr | Chi >
                                                  + delta_kq delta_jp < 0 | E_ir | Chi >
                                                  + delta_kp delta_jr < 0 | E_iq | Chi >
                                                  - sum_{t} ( fock[r,t] Gamma_{ijk,pqt}
                                                            + fock[q,t] Gamma_{ijk,ptr}
                                                            + fock[p,t] Gamma_{ijk,tqr} )
      */
   }

   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){ // anni1 = i ( works in on the bra ) ( smaller than j, k )
      for ( unsigned int crea1 = 0; crea1 < L; crea1++ ){ // crea1 = p ( can be anything )

         const int irrep_center1 = Irreps::directProd( getOrb2Irrep( crea1 ), getOrb2Irrep( anni1 ) );
         const int target_irrep1 = Irreps::directProd( TargetIrrep, irrep_center1 );
         apply_excitation( vector, workspace1, crea1, anni1, TargetIrrep );

         if ( irrep_center1 == 0 ){

            // value = < Chi | E_{crea1,anni1} | 0 >
            const double value = FCIddot( orig_length, workspace1, chi );
            for ( unsigned int j = anni1; j < L; j++ ){
               for ( unsigned int k = anni1; k < L; k++ ){
                  //      i           j       k       p       q       r
                  output[ anni1 + L*( j + L*( k + L*( j + L*( k     + L * crea1 )))) ] += value; // + delta_kq delta_jp < 0 | E_ir | Chi >
                  output[ anni1 + L*( j + L*( k + L*( k + L*( crea1 + L * j     )))) ] += value; // + delta_kp delta_jr < 0 | E_iq | Chi >
               }
            }

         }

         for ( unsigned int crea2 = 0; crea2 < L; crea2++ ){ // crea2 = q
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){ // anni2 = j >= ( i = anni1 )

               const int irrep_center2 = Irreps::directProd( getOrb2Irrep( crea2 ), getOrb2Irrep( anni2 ) );
               const int target_irrep2 = Irreps::directProd( target_irrep1, irrep_center2 );
               const int irrep_center3 = Irreps::directProd( irrep_center1, irrep_center2 );
               apply_excitation( workspace1, workspace2, crea2, anni2, target_irrep1 );

               if ( irrep_center1 == irrep_center2 ){

                  // value = < Chi | E_{crea2,anni2} E_{crea1,anni1} | 0 >
                  const double value = FCIddot( orig_length, workspace2, chi );
                  for ( unsigned int orb = anni1; orb < L; orb++ ){
                     //      i           j           k           p           q           r
                     output[ anni1 + L*( anni2 + L*( orb   + L*( crea1 + L*( orb   + L * crea2 )))) ] -= value; // - delta_kq < 0 | E_ip E_jr | Chi >
                     output[ anni1 + L*( anni2 + L*( orb   + L*( orb   + L*( crea2 + L * crea1 )))) ] -= value; // - delta_kp < 0 | E_ir E_jq | Chi >
                     output[ anni1 + L*( orb   + L*( anni2 + L*( orb   + L*( crea1 + L * crea2 )))) ] -= value; // - delta_jp < 0 | E_iq E_kr | Chi >
                  }

               }

               if (( crea1 >= anni1 ) && ( crea2 >= anni1 )){

                  for ( unsigned int crea3 = (( crea1 == crea2 ) ? crea2 + 1 : crea2 ); crea3 < L; crea3++ ){ // crea3 = r >= ( q = crea2 ) >= ( i = anni1 )
                     for ( unsigned int anni3 = (( anni1 == anni2 ) ? anni1 + 1 : anni1 ); anni3 < L; anni3++ ){ // anni3 = k >= ( i = anni1 )

                        const int irrep_product3 = Irreps::directProd( getOrb2Irrep( crea3 ), getOrb2Irrep( anni3 ) );

                        if ( irrep_center3 == irrep_product3 ){ // I1 x I2 x I3 = Itrivial

                           apply_excitation( workspace2, workspace3, crea3, anni3, target_irrep2 );

                           // value = < Chi | E_{crea3,anni3} E_{crea2,anni2} E_{crea1,anni1} | 0 >
                           double value = FCIddot( orig_length, workspace3, chi );

                           if ( task_fock ){
                              for ( unsigned int t = 0; t < L; t++ ){
                                 // Irrep diagonality of fock is checked by three_rdm values being zero
                                 value -= ( fock[ crea3 + L * t ] * three_rdm[ anni1 + L*( anni2 + L*( anni3 + L*( crea1 + L*( crea2 + L * t     )))) ]
                                          + fock[ crea2 + L * t ] * three_rdm[ anni1 + L*( anni2 + L*( anni3 + L*( crea1 + L*( t     + L * crea3 )))) ]
                                          + fock[ crea1 + L * t ] * three_rdm[ anni1 + L*( anni2 + L*( anni3 + L*( t     + L*( crea2 + L * crea3 )))) ] );
                              }
                           }

                           if ( task_E_zz ){
                              const int number = (( orbz == crea1 ) ? 1 : 0 ) + (( orbz == crea2 ) ? 1 : 0 ) + (( orbz == crea3 ) ? 1 : 0 );
                              if ( number > 0 ){
                                 value -= number * three_rdm[ anni1 + L*( anni2 + L*( anni3 + L*( crea1 + L*( crea2 + L * crea3 )))) ];
                              }
                           }

                           output[ anni1 + L*( anni2 + L*( anni3 + L*( crea1 + L*( crea2 + L * crea3 )))) ] += value;

                        }
                     }
                  }
               }
            }
         }
      }
   }
   delete [] workspace1;
   delete [] workspace2;
   delete [] workspace3;
   if (( task_fock ) || ( task_E_zz )){ delete [] chi; }

   // Make 12-fold permutation symmetric
   for ( unsigned int anni1 = 0; anni1 < L; anni1++ ){
      for ( unsigned int crea1 = anni1; crea1 < L; crea1++ ){
         const int irrep_prod1 = Irreps::directProd( getOrb2Irrep( crea1 ) , getOrb2Irrep( anni1 ) ); // Ic1 x Ia1
         for ( unsigned int crea2 = anni1; crea2 < L; crea2++ ){
            const int irrep_prod2 = Irreps::directProd( irrep_prod1 , getOrb2Irrep( crea2 ) ); // Ic1 x Ia1 x Ic2
            for ( unsigned int anni2 = anni1; anni2 < L; anni2++ ){
               const int irrep_prod3 = Irreps::directProd( irrep_prod2 , getOrb2Irrep( anni2 ) ); // Ic1 x Ia1 x Ic2 x Ia2
               for ( unsigned int crea3 = crea2; crea3 < L; crea3++ ){
                  const int irrep_prod4 = Irreps::directProd( irrep_prod3 , getOrb2Irrep( crea3 ) ); // Ic1 x Ia1 x Ic2 x Ia2 x Ic3
                  for ( unsigned int anni3 = anni1; anni3 < L; anni3++ ){
                     if ( irrep_prod4 == getOrb2Irrep( anni3 )){ // Ic1 x Ia1 x Ic2 x Ia2 x Ic3 == Ia3
                     
                        /*      crea3 >= crea2 >= anni1
                           crea1, anni3, anni2 >= anni1  */
                  
   const double value = output[ anni1 + L * ( anni2 + L * ( anni3 + L * ( crea1 + L * ( crea2 + L * crea3 ) ) ) ) ];
                        output[ anni1 + L * ( anni3 + L * ( anni2 + L * ( crea1 + L * ( crea3 + L * crea2 ) ) ) ) ] = value;
                        
                        output[ anni3 + L * ( anni2 + L * ( anni1 + L * ( crea3 + L * ( crea2 + L * crea1 ) ) ) ) ] = value;
                        output[ anni2 + L * ( anni3 + L * ( anni1 + L * ( crea2 + L * ( crea3 + L * crea1 ) ) ) ) ] = value;
                        
                        output[ anni2 + L * ( anni1 + L * ( anni3 + L * ( crea2 + L * ( crea1 + L * crea3 ) ) ) ) ] = value;
                        output[ anni3 + L * ( anni1 + L * ( anni2 + L * ( crea3 + L * ( crea1 + L * crea2 ) ) ) ) ] = value;
                  
                        output[ crea3 + L * ( crea2 + L * ( crea1 + L * ( anni3 + L * ( anni2 + L * anni1 ) ) ) ) ] = value;
                        output[ crea2 + L * ( crea3 + L * ( crea1 + L * ( anni2 + L * ( anni3 + L * anni1 ) ) ) ) ] = value;
                        
                        output[ crea2 + L * ( crea1 + L * ( crea3 + L * ( anni2 + L * ( anni1 + L * anni3 ) ) ) ) ] = value;
                        output[ crea3 + L * ( crea1 + L * ( crea2 + L * ( anni3 + L * ( anni1 + L * anni2 ) ) ) ) ] = value;
                        
                        output[ crea1 + L * ( crea3 + L * ( crea2 + L * ( anni1 + L * ( anni3 + L * anni2 ) ) ) ) ] = value;
                        output[ crea1 + L * ( crea2 + L * ( crea3 + L * ( anni1 + L * ( anni2 + L * anni3 ) ) ) ) ] = value;
                        
                     }
                  }
               }
            }
         }
      }
   }
   
   for ( unsigned int anni = 0; anni < L; anni++ ){
      for ( unsigned int combo = 0; combo < L*L*L; combo++ ){
         output[ combo + L * L * L * anni * ( 1 + L + L * L ) ] = 0.0;
         output[ anni * ( 1 + L + L * L ) + L * L * L * combo ] = 0.0;
      }
   }
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   return elapsed;

}

double CheMPS2::FCI::CalcSpinSquared(double * vector) const{

   const unsigned int vecLength = getVecLength( 0 );
   double result = 0.0;
      
   #pragma omp parallel for schedule(static) reduction(+:result)
   for ( unsigned int counter = 0; counter < vecLength; counter++ ){
      for ( unsigned int orbi = 0; orbi < L; orbi++ ){
         
         const int irrep_up     = getUpIrrepOfCounter( 0 , counter );
         const int irrep_down   = Irreps::directProd( irrep_up , TargetIrrep );
         const int count_up     = ( counter - irrep_center_jumps[ 0 ][ irrep_up ] ) % numPerIrrep_up[ irrep_up ];
         const int count_down   = ( counter - irrep_center_jumps[ 0 ][ irrep_up ] ) / numPerIrrep_up[ irrep_up ];
         
         // Diagonal terms
         const int diff_ii = lookup_sign_alpha[ irrep_up   ][ orbi + L * orbi ][ count_up   ]
                           - lookup_sign_beta [ irrep_down ][ orbi + L * orbi ][ count_down ]; //Signed integers so subtracting is OK
         const double vector_at_counter_squared = vector[ counter ] * vector[ counter ];
         result += 0.75 * diff_ii * diff_ii * vector_at_counter_squared;
         
         for ( unsigned int orbj = orbi+1; orbj < L; orbj++ ){
         
            // Sz Sz
            const int diff_jj = lookup_sign_alpha[ irrep_up   ][ orbj + L * orbj ][ count_up   ]
                              - lookup_sign_beta [ irrep_down ][ orbj + L * orbj ][ count_down ]; //Signed integers so subtracting is OK
            result += 0.5 * diff_ii * diff_jj * vector_at_counter_squared;
            
            const int irrep_up_bis = Irreps::directProd( irrep_up , Irreps::directProd( getOrb2Irrep( orbi ) , getOrb2Irrep( orbj ) ) );
            
            // - ( a_i,up^+ a_j,up )( a_j,down^+ a_i,down )
            const int sign_down_ji  = lookup_sign_beta [ irrep_down ][ orbj + L * orbi ][ count_down ];
            const int sign_up_ij    = lookup_sign_alpha[ irrep_up   ][ orbi + L * orbj ][ count_up ];
            const int sign_product1 = sign_up_ij * sign_down_ji;
            if ( sign_product1 != 0 ){
               const int cnt_down_ji = lookup_cnt_beta [ irrep_down ][ orbj + L * orbi ][ count_down ];
               const int cnt_up_ij   = lookup_cnt_alpha[ irrep_up   ][ orbi + L * orbj ][ count_up ];
               result -= sign_product1 * vector[ irrep_center_jumps[ 0 ][ irrep_up_bis ] + cnt_up_ij + numPerIrrep_up[ irrep_up_bis ] * cnt_down_ji ] * vector[ counter ];
            }

            // - ( a_j,up^+ a_i,up )( a_i,down^+ a_j,down )
            const int sign_down_ij  = lookup_sign_beta [ irrep_down ][ orbi + L * orbj ][ count_down ];
            const int sign_up_ji    = lookup_sign_alpha[ irrep_up   ][ orbj + L * orbi ][ count_up ];
            const int sign_product2 = sign_up_ji * sign_down_ij;
            if ( sign_product2 != 0 ){
               const int cnt_down_ij = lookup_cnt_beta [ irrep_down ][ orbi + L * orbj ][ count_down ];
               const int cnt_up_ji   = lookup_cnt_alpha[ irrep_up   ][ orbj + L * orbi ][ count_up ];
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

   const unsigned int vecLength = getVecLength( 0 );

   #pragma omp parallel
   {

      int * bits_up   = new int[ L ];
      int * bits_down = new int[ L ];
      
      #pragma omp for schedule(static)
      for ( unsigned int counter = 0; counter < vecLength; counter++ ){
      
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

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const unsigned int vecLength = getVecLength( 0 );
   
   //
   //   Wick's theorem to evaluate the Hamiltonian squared:
   //
   //      H = g_ij E_ij + 0.5 * (ij|kl) E_ij E_kl
   //
   //      H^2 = g_ij g_kl E_ij E_kl
   //          + 0.5 * [ g_ab (ij|kl) + (ab|ij) g_kl ] E_ab E_ij E_kl
   //          + 0.25 * (ab|cd) * (ij|kl) E_ab E_cd E_ij E_kl
   //
   //      Short illustration of what is being done:
   //
   //          E_ij E_kl = (i,s1)^+ (j,s1) (k,s2)^+ (l,s2)
   //
   //                    = (i,s1)^+ (j,s1) (k,s2)^+ (l,s2)
   //                        |        |      |        |
   //                        ----------      ----------
   //                    + (i,s1)^+ (j,s1) (k,s2)^+ (l,s2)
   //                        |        |      |        |
   //                        |        --------        |
   //                        --------------------------
   //                    = num(i,s1) delta(i,j) num(k,s2) delta(k,l)
   //                    + num(i,s1) delta(i,l) delta(s1,s2) [1-num(k,s2)] delta(k,j)
   //          
   //          g_ij g_kl E_ij E_kl = g_ii g_kk num(i,s1) num(k,s2) + g_ik g_ki num(i,s1) [1-num(k,s2)] delta(s1,s2)
   //
   
   #pragma omp parallel
   {

      int * bits_up   = new int[ L ];
      int * bits_down = new int[ L ];
      
      double * Jmat       = new double[ L * L ]; // (ij|kk)( n_k,up + n_k,down )
      double * K_reg_up   = new double[ L * L ]; // (ik|kj)( n_k,up )
      double * K_reg_down = new double[ L * L ]; // (ik|kj)( n_k,down )
      double * K_bar_up   = new double[ L * L ]; // (ik|kj)( 1 - n_k,up )
      double * K_bar_down = new double[ L * L ]; // (ik|kj)( 1 - n_k,down )
      
      int * specific_orbs_irrep = new int[ num_irreps * ( L + 1 ) ];
      for ( unsigned int irrep = 0; irrep < num_irreps; irrep++ ){
         int count = 0;
         for ( unsigned int orb = 0; orb < L; orb++){
            specific_orbs_irrep[ orb + ( L + 1 ) * irrep ] = 0;
            if ( getOrb2Irrep(orb) == irrep ){
               specific_orbs_irrep[ count + ( L + 1 ) * irrep ] = orb;
               count++;
            }
         }
         specific_orbs_irrep[ L + ( L + 1 ) * irrep ] = count;
      }
      
      #pragma omp for schedule(static)
      for (unsigned int counter = 0; counter < vecLength; counter++){
      
         getBitsOfCounter( 0 , counter , bits_up , bits_down ); // Fetch the corresponding bits
         
         // Construct the J and K matrices properly
         for ( unsigned int i = 0; i < L; i++ ){
            for ( unsigned int j = i; j < L; j++ ){
               
               double val_J        = 0.0;
               double val_KregUP   = 0.0;
               double val_KregDOWN = 0.0;
               double val_KbarUP   = 0.0;
               double val_KbarDOWN = 0.0;
               
               if ( getOrb2Irrep(i) == getOrb2Irrep(j) ){
                  for ( unsigned int k = 0; k < L; k++ ){
                     const double temp = getERI(i,k,k,j);
                     val_J        += getERI(i,j,k,k) * ( bits_up[k] + bits_down[k] );
                     val_KregUP   += temp * bits_up[k];
                     val_KregDOWN += temp * bits_down[k];
                     val_KbarUP   += temp * ( 1 - bits_up[k] );
                     val_KbarDOWN += temp * ( 1 - bits_down[k] );
                  }
               }
               
               Jmat[ i + L * j ] = val_J;
               Jmat[ j + L * i ] = val_J;
               K_reg_up[ i + L * j ] = val_KregUP;
               K_reg_up[ j + L * i ] = val_KregUP;
               K_reg_down[ i + L * j ] = val_KregDOWN;
               K_reg_down[ j + L * i ] = val_KregDOWN;
               K_bar_up[ i + L * j ] = val_KbarUP;
               K_bar_up[ j + L * i ] = val_KbarUP;
               K_bar_down[ i + L * j ] = val_KbarDOWN;
               K_bar_down[ j + L * i ] = val_KbarDOWN;
               
            }
         }
         
         double temp = 0.0;
         // G[i,i] (n_i,up + n_i,down) + 0.5 * ( J[i,i] (n_i,up + n_i,down) + K_bar_up[i,i] * n_i,up + K_bar_down[i,i] * n_i,down )
         for ( unsigned int i = 0; i < L; i++ ){
            const int num_i = bits_up[i] + bits_down[i];
            temp += getGmat(i, i) * num_i + 0.5 * ( Jmat[ i + L * i ] * num_i + K_bar_up[ i + L * i ]   * bits_up[i]
                                                                              + K_bar_down[ i + L * i ] * bits_down[i] );
         }
         double myResult = temp*temp;
         
         for ( unsigned int p = 0; p < L; p++ ){
            for ( unsigned int q = 0; q < L; q++ ){
               if ( getOrb2Irrep(p) == getOrb2Irrep(q) ){
            
                  const int special_pq         = bits_up[p] * ( 1 - bits_up[q] ) + bits_down[p] * ( 1 - bits_down[q] );
                  const double GplusJ_pq       = getGmat(p, q) + Jmat[ p + L * q ];
                  const double K_cross_pq_up   = ( K_bar_up[ p + L * q ] - K_reg_up[ p + L * q ]     ) * bits_up[p]   * ( 1 - bits_up[q]   );
                  const double K_cross_pq_down = ( K_bar_down[ p + L * q ] - K_reg_down[ p + L * q ] ) * bits_down[p] * ( 1 - bits_down[q] );
                  
                  myResult += ( GplusJ_pq * ( special_pq * GplusJ_pq + K_cross_pq_up + K_cross_pq_down )
                              + 0.25 * ( K_cross_pq_up * K_cross_pq_up + K_cross_pq_down * K_cross_pq_down ) );
                  
               }
            }
         }
         
         /*
         
            Part which can be optimized most still --> For H2O/6-31G takes 82.8 % of time with Intel compiler
            Optimization can in principle be done with lookup tables + dgemm_ as matvec product --> bit of work :-(
         
            For future reference: the quantity which is computed here:
            
               0.5 * (ak|ci) * (ak|ci) * [ n_a,up * (1-n_k,up) + n_a,down * (1-n_k,down) ] * [ n_c,up * (1-n_i,up) + n_c,down * (1-n_i,down) ]
             - 0.5 * (ak|ci) * (ai|ck) * [ n_a,up * n_c,up * (1-n_i,up) * (1-n_k,up) + n_a,down * n_c,down * (1-n_i,down) * (1-n_k,down) ]
         
         */
         for ( unsigned int k = 0; k < L; k++ ){
            if ( bits_up[k] + bits_down[k] < 2 ){
               for ( unsigned int a = 0; a < L; a++ ){
               
                  const int special_ak    = ( bits_up[a]   * ( 1 - bits_up[k]   )
                                            + bits_down[a] * ( 1 - bits_down[k] ) );
                  const int local_ak_up   = bits_up[a]   * ( 1 - bits_up[k]   );
                  const int local_ak_down = bits_down[a] * ( 1 - bits_down[k] );
               
                  if ( ( special_ak > 0 ) || ( local_ak_up > 0 ) || ( local_ak_down > 0 ) ){
                  
                     const int irrep_ak = Irreps::directProd( getOrb2Irrep(a), getOrb2Irrep(k) );
                        
                     for ( unsigned int i = 0; i < L; i++ ){
                        if ( bits_up[i] + bits_down[i] < 2 ){
                  
                           const int offset     = Irreps::directProd( irrep_ak, getOrb2Irrep(i) ) * ( L + 1 );
                           const int bar_i_up   = 1 - bits_up[i];
                           const int bar_i_down = 1 - bits_down[i];
                           const int max_c_cnt  = specific_orbs_irrep[ L + offset ];
                              
                           for ( int c_cnt = 0; c_cnt < max_c_cnt; c_cnt++ ){
                              const int c            = specific_orbs_irrep[ c_cnt + offset ];
                              const int fact_ic_up   = bits_up[c]   * bar_i_up;
                              const int fact_ic_down = bits_down[c] * bar_i_down;
                              const int prefactor1   = ( fact_ic_up + fact_ic_down ) * special_ak;
                              const int prefactor2   = local_ak_up * fact_ic_up + local_ak_down * fact_ic_down;
                              const double eri_akci  = getERI(a, k, c, i);
                              const double eri_aick  = getERI(a, i, c, k);
                              myResult += 0.5 * eri_akci * ( prefactor1 * eri_akci - prefactor2 * eri_aick );
                           }
                        }
                     }
                  }
               }
            }
         }
         
         output[ counter ] = myResult;
      
      }
      
      delete [] bits_up;
      delete [] bits_down;
      
      delete [] Jmat;
      delete [] K_reg_up;
      delete [] K_reg_down;
      delete [] K_bar_up;
      delete [] K_bar_down;
      
      delete [] specific_orbs_irrep;
   
   }
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   if ( FCIverbose > 0 ){ cout << "FCI::DiagHamSquared : Wall time = " << elapsed << " seconds" << endl; }

}


unsigned int CheMPS2::FCI::LowestEnergyDeterminant() const{

   const unsigned int vecLength = getVecLength( 0 );
   double * energies = new double[ vecLength ];

   // Fetch the Slater determinant energies
   DiagHam( energies );
   
   // Find the determinant with minimum energy
   unsigned int minEindex = 0;
   for ( unsigned int count = 1; count < vecLength; count++ ){
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

void CheMPS2::FCI::FCIdcopy(const unsigned int vecLength, double * origin, double * target){

   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   dcopy_( &length , origin , &inc , target , &inc );

}

double CheMPS2::FCI::FCIddot(const unsigned int vecLength, double * vec1, double * vec2){

   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   return ddot_( &length , vec1 , &inc , vec2 , &inc );

}

double CheMPS2::FCI::FCIfrobeniusnorm(const unsigned int vecLength, double * vec){

   return sqrt( FCIddot( vecLength , vec , vec ) );

}

void CheMPS2::FCI::FCIdaxpy(const unsigned int vecLength, const double alpha, double * vec_x, double * vec_y){

   double factor = alpha;
   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   daxpy_( &length , &factor , vec_x , &inc , vec_y , &inc );

}

void CheMPS2::FCI::FCIdscal(const unsigned int vecLength, const double alpha, double * vec){

   double factor = alpha;
   int length = vecLength; // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   int inc = 1;
   dscal_( &length , &factor , vec , &inc );

}

void CheMPS2::FCI::ClearVector(const unsigned int vecLength, double * vec){

   for ( unsigned int cnt = 0; cnt < vecLength; cnt++ ){ vec[cnt] = 0.0; }

}

void CheMPS2::FCI::FillRandom(const unsigned int vecLength, double * vec){

   for ( unsigned int cnt = 0; cnt < vecLength; cnt++ ){ vec[cnt] = ( ( 2.0 * rand() ) / RAND_MAX ) - 1.0; }

}

double CheMPS2::FCI::GSDavidson(double * inoutput, const int DVDSN_NUM_VEC) const{

   const int veclength = getVecLength( 0 ); // Checked "assert( max_integer >= maxVecLength );" at FCI::StartupIrrepCenter()
   Davidson deBoskabouter( veclength, DVDSN_NUM_VEC,
                                      CheMPS2::DAVIDSON_NUM_VEC_KEEP,
                                      CheMPS2::DAVIDSON_FCI_RTOL,
                                      CheMPS2::DAVIDSON_PRECOND_CUTOFF, false ); // No debug printing for FCI
   double ** whichpointers = new double*[2];

   char instruction = deBoskabouter.FetchInstruction( whichpointers );
   assert( instruction == 'A' );
   if ( inoutput != NULL ){ FCIdcopy( veclength, inoutput, whichpointers[0] ); }
   else { FillRandom( veclength, whichpointers[0] ); }
   DiagHam( whichpointers[1] );

   instruction = deBoskabouter.FetchInstruction( whichpointers );
   while ( instruction == 'B' ){
      matvec( whichpointers[0], whichpointers[1] );
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

   const unsigned int vecLength = getVecLength( 0 );
   for (unsigned int counter = 0; counter < vecLength; counter++){
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

   const unsigned int vecLength = getVecLength( 0 );

   if ( getTargetIrrep() != Irreps::directProd( otherFCI->getTargetIrrep() , getOrb2Irrep( orbIndex ) )){
      ClearVector( vecLength , thisVector );
      return;
   }

   int * bits_up    = new int[ L ];
   int * bits_down  = new int[ L ];
   
   if (( whichOperator=='C') && ( isUp )){
      for (unsigned int counter = 0; counter < vecLength; counter++){
         
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
      for (unsigned int counter = 0; counter < vecLength; counter++){

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
      for (unsigned int counter = 0; counter < vecLength; counter++){

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
      for (unsigned int counter = 0; counter < vecLength; counter++){

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

   const unsigned int vecLength = getVecLength( 0 );

   // Calculate the diagonal of the CG operator
   double * temp = new double[ vecLength ];
   double * diag = new double[ vecLength ];
   CGdiagonal( alpha, beta, eta, diag, temp );

   assert( RealSol != NULL );
   assert( ImagSol != NULL );
   assert( fabs( eta ) > 0.0 );

   /* 
         ( alpha + beta H + I eta ) Solution = RHS

      is solved with the conjugate gradient (CG) method. Solution = RealSol + I * ImagSol.
      CG requires a symmetric positive definite operator. Therefore:

         [ ( alpha + beta H )^2 + eta^2 ] * Solution = ( alpha + beta H - I eta ) * RHS

      Clue: Solve for ImagSol first. RealSol is then simply

         RealSol = - ( alpha + beta H ) / eta * ImagSol
   */

   /**** Solve for ImagSol ****/
   double ** pointers = new double*[ 3 ];
   double RMSerror = 0.0;
   {
      ConjugateGradient CG( vecLength, CheMPS2::CONJ_GRADIENT_RTOL, CheMPS2::CONJ_GRADIENT_PRECOND_CUTOFF, false );
      char instruction = CG.step( pointers );
      assert( instruction == 'A' );
      for (unsigned int cnt = 0; cnt < vecLength; cnt++){ pointers[ 0 ][ cnt ] = - eta * RHS[ cnt ] / diag[ cnt ]; } // Initial guess
      for (unsigned int cnt = 0; cnt < vecLength; cnt++){ pointers[ 1 ][ cnt ] = diag[ cnt ]; }                      // Diagonal of the operator
      for (unsigned int cnt = 0; cnt < vecLength; cnt++){ pointers[ 2 ][ cnt ] = - eta * RHS[ cnt ]; }               // RHS of the problem
      instruction = CG.step( pointers );
      assert( instruction == 'B' );
      while ( instruction == 'B' ){
         CGoperator( alpha, beta, eta, pointers[ 0 ], temp, pointers[ 1 ] );
         instruction = CG.step( pointers );
      }
      assert( instruction == 'C' );
      FCIdcopy( vecLength, pointers[ 0 ], ImagSol );
      RMSerror += pointers[ 1 ][ 0 ] * pointers[ 1 ][ 0 ];
   }

   /**** Solve for RealSol ****/
   {
      ConjugateGradient CG( vecLength, CheMPS2::CONJ_GRADIENT_RTOL, CheMPS2::CONJ_GRADIENT_PRECOND_CUTOFF, false );
      char instruction = CG.step( pointers );
      assert( instruction == 'A' );
      CGAlphaPlusBetaHAM( - alpha / eta, - beta / eta, ImagSol, pointers[ 0 ] );                      // Initial guess real part can be obtained from the imaginary part
      for (unsigned int cnt = 0; cnt < vecLength; cnt++){ pointers[ 1 ][ cnt ] = diag[ cnt ]; } // Diagonal of the operator
      CGAlphaPlusBetaHAM( alpha, beta, RHS, pointers[ 2 ] );                                          // RHS of the problem
      instruction = CG.step( pointers );
      assert( instruction == 'B' );
      while ( instruction == 'B' ){
         CGoperator( alpha, beta, eta, pointers[ 0 ], temp, pointers[ 1 ] );
         instruction = CG.step( pointers );
      }
      assert( instruction == 'C' );
      FCIdcopy( vecLength, pointers[ 0 ], RealSol );
      RMSerror += pointers[ 1 ][ 0 ] * pointers[ 1 ][ 0 ];
   }
   RMSerror = sqrt( RMSerror );
   delete [] pointers;
   delete [] temp;
   delete [] diag;

   if (( checkError ) && ( FCIverbose > 0 )){
      cout << "FCI::CGSolveSystem : RMS error when checking the solution = " << RMSerror << endl;
   }

}

void CheMPS2::FCI::CGAlphaPlusBetaHAM(const double alpha, const double beta, double * in, double * out) const{

   matvec( in , out );
   const unsigned int vecLength = getVecLength( 0 );
   const double prefactor = alpha + beta * getEconst(); // matvec does only the parts with second quantized operators
   for (unsigned int cnt = 0; cnt < vecLength; cnt++){
      out[ cnt ] = prefactor * in[ cnt ] + beta * out[ cnt ]; // out = ( alpha + beta * H ) * in
   }

}

void CheMPS2::FCI::CGoperator(const double alpha, const double beta, const double eta, double * in, double * temp, double * out) const{

   const unsigned int vecLength = getVecLength( 0 );
   CGAlphaPlusBetaHAM( alpha, beta, in,   temp ); // temp  = ( alpha + beta * H )   * in
   CGAlphaPlusBetaHAM( alpha, beta, temp, out  ); // out   = ( alpha + beta * H )^2 * in
   FCIdaxpy( vecLength, eta*eta, in, out );       // out   = [ ( alpha + beta * H )^2 + eta*eta ] * in

}

void CheMPS2::FCI::CGdiagonal(const double alpha, const double beta, const double eta, double * diagonal, double * workspace) const{

   // diagonal becomes diag[ ( alpha + beta * H )^2 + eta*eta ]
   
   DiagHam( diagonal );
   DiagHamSquared( workspace );
   
   const unsigned int vecLength = getVecLength( 0 );
   const double alpha_bis = alpha + beta * getEconst();
   const double factor1 = alpha_bis * alpha_bis + eta * eta;
   const double factor2 = 2 * alpha_bis * beta;
   const double factor3 = beta * beta;
   for (unsigned int row = 0; row < vecLength; row++){
      diagonal[ row ] = factor1 + factor2 * diagonal[ row ] + factor3 * workspace[ row ];
   }
   
   if ( FCIverbose>1 ){
      double minval = diagonal[0];
      double maxval = diagonal[0];
      for (unsigned int cnt = 1; cnt < vecLength; cnt++){
         if ( diagonal[ cnt ] > maxval ){ maxval = diagonal[ cnt ]; }
         if ( diagonal[ cnt ] < minval ){ minval = diagonal[ cnt ]; }
      }
      cout << "FCI::CGdiagonal : Minimum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << minval << endl;
      cout << "FCI::CGdiagonal : Maximum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << maxval << endl;
   }

}

void CheMPS2::FCI::RetardedGF(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF) const{

   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );

   // G( omega, alpha, beta, eta ) = < 0 | a_{alpha,spin}  [ omega - Ham + E_0 + I*eta ]^{-1} a^+_{beta,spin} | 0 > (addition amplitude)
   //                              + < 0 | a^+_{beta,spin} [ omega + Ham - E_0 + I*eta ]^{-1} a_{alpha,spin}  | 0 > (removal  amplitude)

   double Realpart, Imagpart;
   RetardedGF_addition(omega, eta, orb_alpha, orb_beta, isUp, GSenergy, GSvector, Ham, &Realpart, &Imagpart);
   RePartGF[0] = Realpart; // Set
   ImPartGF[0] = Imagpart; // Set
   
   RetardedGF_removal( omega, eta, orb_alpha, orb_beta, isUp, GSenergy, GSvector, Ham, &Realpart, &Imagpart);
   RePartGF[0] += Realpart; // Add
   ImPartGF[0] += Imagpart;
   
   if ( FCIverbose>0 ){
      cout << "FCI::RetardedGF : G( omega = " << omega << " ; eta = " << eta << " ; i = " << orb_alpha << " ; j = " << orb_beta << " ) = " << RePartGF[0] << " + I * " << ImPartGF[0] << endl;
      cout << "                  Local density of states (LDOS) = " << - ImPartGF[0] / M_PI << endl;
   }

}

void CheMPS2::FCI::GFmatrix_addition(const double alpha, const double beta, const double eta, int * orbsLeft, const unsigned int numLeft, int * orbsRight, const unsigned int numRight, const bool isUp, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartsGF, double * ImPartsGF, double ** TwoRDMreal, double ** TwoRDMimag, double ** TwoRDMadd) const{

   /*
                                                                         1
       GF[i + numLeft * j] = < 0 | a_{orbsLeft[i], spin} -------------------------------- a^+_{orbsRight[j], spin} | 0 >
                                                          [ alpha + beta * Ham + I*eta ]
   */

   // Check whether some stuff is OK
   assert( numLeft  > 0 );
   assert( numRight > 0 );
   for (unsigned int cnt = 0; cnt < numLeft;  cnt++){ int orbl = orbsLeft[  cnt ]; assert((orbl < L) && (orbl >= 0)); }
   for (unsigned int cnt = 0; cnt < numRight; cnt++){ int orbr = orbsRight[ cnt ]; assert((orbr < L) && (orbr >= 0)); }
   assert( RePartsGF != NULL );
   assert( ImPartsGF != NULL );
   for ( unsigned int counter = 0; counter < numLeft * numRight; counter++ ){
       RePartsGF[ counter ] = 0.0;
       ImPartsGF[ counter ] = 0.0;
   }
   const unsigned int Lpow4 = L*L*L*L;
   for ( unsigned int cnt = 0; cnt < numRight; cnt++ ){
      if ( TwoRDMreal != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){ TwoRDMreal[ cnt ][ elem ] = 0.0; } }
      if ( TwoRDMimag != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){ TwoRDMimag[ cnt ][ elem ] = 0.0; } }
      if ( TwoRDMadd  != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){  TwoRDMadd[ cnt ][ elem ] = 0.0; } }
   }
   
   const bool isOK = ( isUp ) ? ( getNel_up() < L ) : ( getNel_down() < L ); // The electron can be added
   for ( unsigned int cnt_right = 0; cnt_right < numRight; cnt_right++ ){
   
      const int orbitalRight = orbsRight[ cnt_right ];
      bool matchingIrrep = false;
      for ( unsigned int cnt_left = 0; cnt_left < numLeft; cnt_left++ ){
         if ( getOrb2Irrep( orbsLeft[ cnt_left] ) == getOrb2Irrep( orbitalRight ) ){ matchingIrrep = true; }
      }
      
      if ( isOK && matchingIrrep ){
      
         const unsigned int addNelUP   = getNel_up()   + ((isUp) ? 1 : 0);
         const unsigned int addNelDOWN = getNel_down() + ((isUp) ? 0 : 1);
         const int addIrrep = Irreps::directProd( getTargetIrrep(), getOrb2Irrep( orbitalRight ) );
         
         CheMPS2::FCI additionFCI( Ham, addNelUP, addNelDOWN, addIrrep, maxMemWorkMB, FCIverbose );
         const unsigned int addVecLength = additionFCI.getVecLength( 0 );
         double * addVector = new double[ addVecLength ];
         additionFCI.ActWithSecondQuantizedOperator( 'C', isUp, orbitalRight, addVector, this, GSvector ); // | addVector > = a^+_right,spin | GSvector >
         
         double * RealPartSolution = new double[ addVecLength ];
         double * ImagPartSolution = new double[ addVecLength ];
         additionFCI.CGSolveSystem( alpha, beta, eta, addVector, RealPartSolution, ImagPartSolution );
         
         if ( TwoRDMreal != NULL ){ additionFCI.Fill2RDM( RealPartSolution, TwoRDMreal[ cnt_right ] ); }
         if ( TwoRDMimag != NULL ){ additionFCI.Fill2RDM( ImagPartSolution, TwoRDMimag[ cnt_right ] ); }
         if ( TwoRDMadd  != NULL ){ additionFCI.Fill2RDM( addVector,         TwoRDMadd[ cnt_right ] ); }
         
         for ( unsigned int cnt_left = 0; cnt_left < numLeft; cnt_left++ ){
            const int orbitalLeft = orbsLeft[ cnt_left ];
            if ( getOrb2Irrep( orbitalLeft ) == getOrb2Irrep( orbitalRight ) ){
               additionFCI.ActWithSecondQuantizedOperator( 'C', isUp, orbitalLeft, addVector, this, GSvector ); // | addVector > = a^+_left,spin | GSvector >
               RePartsGF[ cnt_left + numLeft * cnt_right ] = FCIddot( addVecLength, addVector, RealPartSolution );
               ImPartsGF[ cnt_left + numLeft * cnt_right ] = FCIddot( addVecLength, addVector, ImagPartSolution );
            }
         }
         
         delete [] RealPartSolution;
         delete [] ImagPartSolution;
         delete [] addVector;
       
      }
   }

}

void CheMPS2::FCI::RetardedGF_addition(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag, double * TwoRDMadd) const{

   // Addition amplitude < 0 | a_{alpha, spin} [ omega - Ham + E_0 + I*eta ]^{-1} a^+_{beta, spin} | 0 >
   
   double ** TwoRDMreal_wrap = NULL; if ( TwoRDMreal != NULL ){ TwoRDMreal_wrap = new double*[1]; TwoRDMreal_wrap[0] = TwoRDMreal; }
   double ** TwoRDMimag_wrap = NULL; if ( TwoRDMimag != NULL ){ TwoRDMimag_wrap = new double*[1]; TwoRDMimag_wrap[0] = TwoRDMimag; }
   double **  TwoRDMadd_wrap = NULL; if (  TwoRDMadd != NULL ){  TwoRDMadd_wrap = new double*[1];  TwoRDMadd_wrap[0] = TwoRDMadd;  }
   
   int orb_left  = orb_alpha;
   int orb_right = orb_beta;
   
   GFmatrix_addition( omega + GSenergy, -1.0, eta, &orb_left, 1, &orb_right, 1, isUp, GSvector, Ham, RePartGF, ImPartGF, TwoRDMreal_wrap, TwoRDMimag_wrap, TwoRDMadd_wrap );
   
   if ( TwoRDMreal != NULL ){ delete [] TwoRDMreal_wrap; }
   if ( TwoRDMimag != NULL ){ delete [] TwoRDMimag_wrap; }
   if (  TwoRDMadd != NULL ){ delete [] TwoRDMadd_wrap;  }

}

void CheMPS2::FCI::GFmatrix_removal(const double alpha, const double beta, const double eta, int * orbsLeft, const unsigned int numLeft, int * orbsRight, const unsigned int numRight, const bool isUp, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartsGF, double * ImPartsGF, double ** TwoRDMreal, double ** TwoRDMimag, double ** TwoRDMrem) const{

   /*
                                                                           1
       GF[i + numLeft * j] = < 0 | a^+_{orbsLeft[i], spin} -------------------------------- a_{orbsRight[j], spin} | 0 >
                                                            [ alpha + beta * Ham + I*eta ]
   */
   
   // Check whether some stuff is OK
   assert( numLeft  > 0 );
   assert( numRight > 0 );
   for (unsigned int cnt = 0; cnt < numLeft;  cnt++){ int orbl = orbsLeft [ cnt ]; assert((orbl < L) && (orbl >= 0)); }
   for (unsigned int cnt = 0; cnt < numRight; cnt++){ int orbr = orbsRight[ cnt ]; assert((orbr < L) && (orbr >= 0)); }
   assert( RePartsGF != NULL );
   assert( ImPartsGF != NULL );
   for ( unsigned int counter = 0; counter < numLeft * numRight; counter++ ){
       RePartsGF[ counter ] = 0.0;
       ImPartsGF[ counter ] = 0.0;
   }
   const unsigned int Lpow4 = L*L*L*L;
   for ( unsigned int cnt = 0; cnt < numRight; cnt++ ){
      if ( TwoRDMreal != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){ TwoRDMreal[ cnt ][ elem ] = 0.0; } }
      if ( TwoRDMimag != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){ TwoRDMimag[ cnt ][ elem ] = 0.0; } }
      if ( TwoRDMrem  != NULL ){ for ( unsigned int elem = 0; elem < Lpow4; elem++ ){  TwoRDMrem[ cnt ][ elem ] = 0.0; } }
   }
   
   const bool isOK = ( isUp ) ? ( getNel_up() > 0 ) : ( getNel_down() > 0 ); // The electron can be removed
   for ( unsigned int cnt_right = 0; cnt_right < numRight; cnt_right++ ){
   
      const int orbitalRight = orbsRight[ cnt_right ];
      bool matchingIrrep = false;
      for ( unsigned int cnt_left = 0; cnt_left < numLeft; cnt_left++ ){
         if ( getOrb2Irrep( orbsLeft[ cnt_left] ) == getOrb2Irrep( orbitalRight ) ){ matchingIrrep = true; }
      }
      
      if ( isOK && matchingIrrep ){
      
         const unsigned int removeNelUP   = getNel_up()   - ((isUp) ? 1 : 0);
         const unsigned int removeNelDOWN = getNel_down() - ((isUp) ? 0 : 1);
         const int removeIrrep = Irreps::directProd( getTargetIrrep(), getOrb2Irrep( orbitalRight ) );
         
         CheMPS2::FCI removalFCI( Ham, removeNelUP, removeNelDOWN, removeIrrep, maxMemWorkMB, FCIverbose );
         const unsigned int removeVecLength = removalFCI.getVecLength( 0 );
         double * removeVector = new double[ removeVecLength ];
         removalFCI.ActWithSecondQuantizedOperator( 'A', isUp, orbitalRight, removeVector, this, GSvector ); // | removeVector > = a_right,spin | GSvector >
         
         double * RealPartSolution = new double[ removeVecLength ];
         double * ImagPartSolution = new double[ removeVecLength ];
         removalFCI.CGSolveSystem( alpha, beta, eta, removeVector, RealPartSolution, ImagPartSolution );
         
         if ( TwoRDMreal != NULL ){ removalFCI.Fill2RDM( RealPartSolution, TwoRDMreal[ cnt_right ] ); }
         if ( TwoRDMimag != NULL ){ removalFCI.Fill2RDM( ImagPartSolution, TwoRDMimag[ cnt_right ] ); }
         if ( TwoRDMrem  != NULL ){ removalFCI.Fill2RDM( removeVector,      TwoRDMrem[ cnt_right ] ); }
         
         for ( unsigned int cnt_left = 0; cnt_left < numLeft; cnt_left++ ){
            const int orbitalLeft = orbsLeft[ cnt_left ];
            if ( getOrb2Irrep( orbitalLeft ) == getOrb2Irrep( orbitalRight ) ){
               removalFCI.ActWithSecondQuantizedOperator( 'A', isUp, orbitalLeft, removeVector, this, GSvector ); // | removeVector > = a_left,spin | GSvector >
               RePartsGF[ cnt_left + numLeft * cnt_right ] = FCIddot( removeVecLength, removeVector, RealPartSolution );
               ImPartsGF[ cnt_left + numLeft * cnt_right ] = FCIddot( removeVecLength, removeVector, ImagPartSolution );
            }
         }
         
         delete [] RealPartSolution;
         delete [] ImagPartSolution;
         delete [] removeVector;
       
      }
   }

}

void CheMPS2::FCI::RetardedGF_removal(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const bool isUp, const double GSenergy, double * GSvector, CheMPS2::Hamiltonian * Ham, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag, double * TwoRDMrem) const{

   // Removal amplitude < 0 | a^+_{beta, spin} [ omega + Ham - E_0 + I*eta ]^{-1} a_{alpha, spin} | 0 >
   
   double ** TwoRDMreal_wrap = NULL; if ( TwoRDMreal != NULL ){ TwoRDMreal_wrap = new double*[1]; TwoRDMreal_wrap[0] = TwoRDMreal; }
   double ** TwoRDMimag_wrap = NULL; if ( TwoRDMimag != NULL ){ TwoRDMimag_wrap = new double*[1]; TwoRDMimag_wrap[0] = TwoRDMimag; }
   double **  TwoRDMrem_wrap = NULL; if (  TwoRDMrem != NULL ){  TwoRDMrem_wrap = new double*[1];  TwoRDMrem_wrap[0] = TwoRDMrem;  }
   
   int orb_left  = orb_beta;
   int orb_right = orb_alpha;
   
   // orb_alpha = orb_right in this case !!
   GFmatrix_removal( omega - GSenergy, 1.0, eta, &orb_left, 1, &orb_right, 1, isUp, GSvector, Ham, RePartGF, ImPartGF, TwoRDMreal_wrap, TwoRDMimag_wrap, TwoRDMrem_wrap );
   
   if ( TwoRDMreal != NULL ){ delete [] TwoRDMreal_wrap; }
   if ( TwoRDMimag != NULL ){ delete [] TwoRDMimag_wrap; }
   if (  TwoRDMrem != NULL ){ delete [] TwoRDMrem_wrap;  }

}

void CheMPS2::FCI::DensityResponseGF(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF) const{

   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );

   // X( omega, alpha, beta, eta ) = < 0 | ( n_alpha - <0| n_alpha |0> ) [ omega - Ham + E_0 + I*eta ]^{-1} ( n_beta  - <0| n_beta  |0> ) | 0 > (forward  amplitude)
   //                              - < 0 | ( n_beta  - <0| n_beta  |0> ) [ omega + Ham - E_0 + I*eta ]^{-1} ( n_alpha - <0| n_alpha |0> ) | 0 > (backward amplitude)
   
   double Realpart, Imagpart;
   DensityResponseGF_forward( omega, eta, orb_alpha, orb_beta, GSenergy, GSvector, &Realpart, &Imagpart);
   RePartGF[0] = Realpart; // Set
   ImPartGF[0] = Imagpart; // Set
   
   DensityResponseGF_backward(omega, eta, orb_alpha, orb_beta, GSenergy, GSvector, &Realpart, &Imagpart);
   RePartGF[0] -= Realpart; // Subtract !!!
   ImPartGF[0] -= Imagpart; // Subtract !!!

   if ( FCIverbose>0 ){
      cout << "FCI::DensityResponseGF : X( omega = " << omega << " ; eta = " << eta << " ; i = " << orb_alpha << " ; j = " << orb_beta << " ) = " << RePartGF[0] << " + I * " << ImPartGF[0] << endl;
      cout << "                         Local density-density response (LDDR) = " << - ImPartGF[0] / M_PI << endl;
   }

}

void CheMPS2::FCI::DensityResponseGF_forward(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag, double * TwoRDMdens) const{

   // Forward amplitude: < 0 | ( n_alpha - <0| n_alpha |0> ) [ omega - Ham + E_0 + I*eta ]^{-1} ( n_beta  - <0| n_beta  |0> ) | 0 >
   
   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const unsigned int vecLength = getVecLength( 0 );
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

   if ( TwoRDMdens != NULL ){ Fill2RDM( densityBetaVector , TwoRDMdens ); } // Sets the TwoRDMdens
   if ( orb_alpha != orb_beta ){ delete [] densityBetaVector; }
   delete [] densityAlphaVector;

}

void CheMPS2::FCI::DensityResponseGF_backward(const double omega, const double eta, const unsigned int orb_alpha, const unsigned int orb_beta, const double GSenergy, double * GSvector, double * RePartGF, double * ImPartGF, double * TwoRDMreal, double * TwoRDMimag, double * TwoRDMdens) const{

   // Backward amplitude: < 0 | ( n_beta  - <0| n_beta  |0> ) [ omega + Ham - E_0 + I*eta ]^{-1} ( n_alpha - <0| n_alpha |0> ) | 0 >
   
   assert( ( orb_alpha<L ) && ( orb_beta<L ) ); // Orbital indices within bound
   assert( RePartGF != NULL );
   assert( ImPartGF != NULL );
   
   const unsigned int vecLength = getVecLength( 0 );
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

   if ( TwoRDMdens != NULL ){ Fill2RDM( densityAlphaVector , TwoRDMdens ); } // Sets the TwoRDMdens
   if ( orb_alpha != orb_beta ){ delete [] densityBetaVector; }
   delete [] densityAlphaVector;

}


