/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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

CheMPS2::FCI::FCI(Problem * Prob, const int FCIverbose_in){

   // Amount of output
   FCIverbose = FCIverbose_in;

   // Number of orbitals
   assert( Prob->checkConsistency() );
   L = Prob->gL();
   assert( L <= CHAR_BIT * sizeof(unsigned int) );

   // Everything which has to do with irreps
   CheMPS2::Irreps myIrreps( Prob->gSy() );
   NumIrreps         = myIrreps.getNumberOfIrreps();
   TargetIrrep       = Prob->gIrrep();
   IrrepProductTable = new int[ NumIrreps * NumIrreps ];
   orb2irrep         = new int[ L ];
   for (unsigned int row = 0; row < NumIrreps; row++){
      for (unsigned int col = 0; col < NumIrreps; col++){
         IrrepProductTable[ row + NumIrreps * col ] = myIrreps.directProd( row , col );
      }
   }
   for (unsigned int orb = 0; orb < L; orb++){ orb2irrep[ orb ] = Prob->gIrrep( orb ); }
   
   // Number of up and down electrons: only spin projection symmetry is used
   Nel_up   = ( Prob->gN() + Prob->gTwoS() )/2;
   Nel_down = ( Prob->gN() - Prob->gTwoS() )/2;
   
   // Copy the Hamiltonian over for speed
   Econstant = Prob->gEconst();
   Hmat = new double[ L * L * L * L ];
   for (unsigned int orb1 = 0; orb1 < L; orb1++){
      for (unsigned int orb2 = 0; orb2 < L; orb2++){
         for (unsigned int orb3 = 0; orb3 < L; orb3++){
            for (unsigned int orb4 = 0; orb4 < L; orb4++){
               Hmat[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ] = Prob->gMxElement( orb1 , orb2 , orb3 , orb4 );
            }
         }
      }
   }
   
   // Set all other internal variables
   Startup();
   if ( FCIverbose > 0 ){ cout << "FCI::FCI : Number of variables in the FCI vector = " << getVecLength() << endl; }

}

CheMPS2::FCI::FCI(Problem * Prob, const unsigned int theNel_up, const unsigned int theNel_down, const int FCIverbose_in){

   // Amount of output
   FCIverbose = FCIverbose_in;

   // Number of orbitals
   assert( Prob->checkConsistency() );
   L = Prob->gL();
   assert( L <= CHAR_BIT * sizeof(unsigned int) );
   assert( theNel_up>=0 );
   assert( theNel_down>=0 );
   assert( theNel_up<=L );
   assert( theNel_down<=L );
   Nel_up = theNel_up;
   Nel_down = theNel_down;

   // Everything which has to do with irreps
   CheMPS2::Irreps myIrreps( Prob->gSy() );
   NumIrreps         = myIrreps.getNumberOfIrreps();
   TargetIrrep       = Prob->gIrrep();
   IrrepProductTable = new int[ NumIrreps * NumIrreps ];
   orb2irrep         = new int[ L ];
   for (unsigned int row = 0; row < NumIrreps; row++){
      for (unsigned int col = 0; col < NumIrreps; col++){
         IrrepProductTable[ row + NumIrreps * col ] = myIrreps.directProd( row , col );
      }
   }
   for (unsigned int orb = 0; orb < L; orb++){ orb2irrep[ orb ] = Prob->gIrrep( orb ); }
   
   // Copy the Hamiltonian over
   Econstant = Prob->gEconst();
   Hmat = new double[ L * L * L * L ];
   for (unsigned int orb1 = 0; orb1 < L; orb1++){
      for (unsigned int orb2 = 0; orb2 < L; orb2++){
         for (unsigned int orb3 = 0; orb3 < L; orb3++){
            for (unsigned int orb4 = 0; orb4 < L; orb4++){
               Hmat[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ] = Prob->gMxElement( orb1 , orb2 , orb3 , orb4 );
            }
         }
      }
   }
   
   // Set all other internal variables
   Startup();
   if ( FCIverbose > 0 ){ cout << "FCI::FCI : Number of variables in the FCI vector = " << getVecLength() << endl; }

}

CheMPS2::FCI::~FCI(){

   delete [] IrrepProductTable;
   delete [] orb2irrep;
   delete [] Hmat;
   for (unsigned int irrep=0; irrep<NumIrreps; irrep++){
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
   delete [] jumps;

}

unsigned int CheMPS2::FCI::getL() const{ return L; }

unsigned int CheMPS2::FCI::getNel_up() const{ return Nel_up; }

unsigned int CheMPS2::FCI::getNel_down() const{ return Nel_down; }

unsigned long long CheMPS2::FCI::getVecLength() const{ return jumps[ NumIrreps ]; }

unsigned int CheMPS2::FCI::getNumIrreps() const{ return NumIrreps; }

int CheMPS2::FCI::getTargetIrrep() const{ return TargetIrrep; }

int CheMPS2::FCI::getIrrepProduct(const int Irrep_row, const int Irrep_col) const{ return IrrepProductTable[ Irrep_row + NumIrreps * Irrep_col ]; }

int CheMPS2::FCI::getOrb2Irrep(const int orb) const{ return orb2irrep[ orb ]; }

double CheMPS2::FCI::getEconst() const{ return Econstant; }

double CheMPS2::FCI::getHmat(const int orb1, const int orb2, const int orb3, const int orb4) const{ return Hmat[ orb1 + L * ( orb2 + L * ( orb3 + L * orb4 ) ) ]; }

void CheMPS2::FCI::Startup(){

   // Variable which is only needed in Startup: 2^L
   unsigned int TwoPowL = 1;
   for (unsigned int orb = 0; orb < L; orb++){ TwoPowL *= 2; }

   // Create a bunch of stuff
   numPerIrrep_up   = new unsigned int[ NumIrreps];
   numPerIrrep_down = new unsigned int[ NumIrreps];
   str2cnt_up       = new int*[NumIrreps];
   str2cnt_down     = new int*[NumIrreps];
   cnt2str_up       = new unsigned int*[NumIrreps];
   cnt2str_down     = new unsigned int*[NumIrreps];

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
         if (bits[orb]){
            Nparticles++;
            Irrep = getIrrepProduct( Irrep , orb2irrep[ orb ] );
         }
      }
      
      // If allowed: set the corresponding str2cnt to the correct counter and keep track of the number of allowed vectors
      for (unsigned int irr = 0; irr < NumIrreps; irr++){
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
   
   delete [] bits; // Delete temporary helper array
   
   // Fill the reverse info array: cnt2str
   for (unsigned int irrep = 0; irrep < NumIrreps; irrep++){
   
      if ( FCIverbose > 1 ){
         cout << "FCI::Startup : For irrep " << irrep << " there are " << numPerIrrep_up  [ irrep ] << " alpha Slater determinants and "
                                                                       << numPerIrrep_down[ irrep ] <<  " beta Slater determinants." << endl;
      }
                                                                              
      cnt2str_up  [ irrep ] = new unsigned int[ numPerIrrep_up  [ irrep ] ];
      cnt2str_down[ irrep ] = new unsigned int[ numPerIrrep_down[ irrep ] ];
      
      for (unsigned int bitstring = 0; bitstring < TwoPowL; bitstring++){
         if ( str2cnt_up[ irrep ][ bitstring ] != -1 ){
            cnt2str_up[ irrep ][ str2cnt_up[ irrep ][ bitstring ] ] = bitstring;
         }
         if ( str2cnt_down[ irrep ][ bitstring ] != -1 ){
            cnt2str_down[ irrep ][ str2cnt_down[ irrep ][ bitstring ] ] = bitstring;
         }
      }
   
   }
   
   // Jumps contains information on irrep of a FCI array index
   jumps = new unsigned long long[ NumIrreps+1 ];
   jumps[ 0 ] = 0;
   for (unsigned int irrep = 0; irrep < NumIrreps; irrep++){
      const int irrep_down = getIrrepProduct( irrep , TargetIrrep );
      unsigned long long temp  = numPerIrrep_up  [ irrep      ];
                         temp *= numPerIrrep_down[ irrep_down ];
      jumps[ irrep+1 ] = jumps[ irrep ] + temp;
   }

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

int CheMPS2::FCI::getUpIrrepOfCounter(const unsigned long long counter) const{

   int irrep_up = NumIrreps;
   while ( counter < jumps[ irrep_up-1 ] ){ irrep_up--; }
   return irrep_up-1;
   
}

void CheMPS2::FCI::getBitsOfCounter(const unsigned long long counter, int * bits_up, int * bits_down) const{

   const int irrep_up   = getUpIrrepOfCounter( counter );
   const int irrep_down = getIrrepProduct( irrep_up , TargetIrrep );
   
   const unsigned int count_up   = ( counter - jumps[ irrep_up ] ) % numPerIrrep_up[ irrep_up ];
   const unsigned int count_down = ( counter - jumps[ irrep_up ] ) / numPerIrrep_up[ irrep_up ];
   
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
   for (unsigned int cnt=0; cnt<L; cnt++){
      if ( bits_up  [ cnt ] ){ irrep_up   = getIrrepProduct( irrep_up   , getOrb2Irrep( cnt ) ); }
      if ( bits_down[ cnt ] ){ irrep_down = getIrrepProduct( irrep_down , getOrb2Irrep( cnt ) ); }
   }
   
   const int counter_up   = str2cnt_up  [ irrep_up   ][ string_up   ];
   const int counter_down = str2cnt_down[ irrep_down ][ string_down ];
   
   if (( counter_up == -1 ) || ( counter_down == -1 )){ return 0.0; }
   
   return vector[ jumps[ irrep_up ] + counter_up + numPerIrrep_up[ irrep_up ] * counter_down ];

}

void CheMPS2::FCI::ActWithSecondQuantizedOperator(const char whichOperator, const bool isUp, const unsigned int orbIndex, double * thisVector, FCI * otherFCI, double * otherVector){

   assert( ( whichOperator=='C' ) || ( whichOperator=='A' ) || ( whichOperator=='N' ) ); //Operator should be a (C) Creator, (A) Annihilator, or (N) Number operator
   assert( orbIndex>=0 );
   assert( orbIndex<L  );
   assert( L==otherFCI->getL() );

   const unsigned long long vecLength = getVecLength();

   if ((( whichOperator=='C' ) || ( whichOperator=='A' )) && ( getTargetIrrep() != getIrrepProduct( otherFCI->getTargetIrrep() , getOrb2Irrep( orbIndex ) ))){
      FCIdclear( vecLength , thisVector );
      return;
   }
   if ((whichOperator=='N' ) && ( ( getTargetIrrep() != otherFCI->getTargetIrrep() ) || ( vecLength != otherFCI->getVecLength() ) )){
      FCIdclear( vecLength , thisVector );
      return;
   }

   int * bits_up    = new int[ L ];
   int * bits_down  = new int[ L ];
   
   for (unsigned long long counter = 0; counter < vecLength; counter++){
   
      // Clear the corresponding entry
      thisVector[ counter ] = 0.0;
      
      // Get the bits corresponding to this counter
      getBitsOfCounter( counter , bits_up , bits_down );
      
      // Fetch the corresponding coefficient
      if ( isUp ){ // Down (beta) electron
      
         if (( whichOperator=='C' ) && ( bits_up[ orbIndex ] == 1 )){ // Operator = creator
            bits_up[ orbIndex ] = 0;
            int phase = 1;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_up[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         }
         if (( whichOperator=='A' ) && ( bits_up[ orbIndex ] == 0 )){ // Operator = annihilator
            bits_up[ orbIndex ] = 1;
            int phase = 1;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_up[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         }
         if ( whichOperator=='N' ){ // Operator = particle number
            thisVector[ counter ] = bits_up[ orbIndex ] * otherVector[ counter ];
         }
         
      } else { // Down (beta) electron
      
         if (( whichOperator=='C' ) && ( bits_down[ orbIndex ] == 1 )){ // Operator = creator
            bits_down[ orbIndex ] = 0;
            int phase = (( Nel_up % 2 ) == 0) ? 1 : -1;
            for (unsigned int cnt = 0; cnt<orbIndex; cnt++){ if ( bits_down[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         }
         if (( whichOperator=='A' ) && ( bits_down[ orbIndex ] == 0 )){ // Operator = annihilator
            bits_down[ orbIndex ] = 1;
            int phase = (( Nel_up % 2 ) == 0) ? 1 : -1;
            for (unsigned int cnt = 0; cnt < orbIndex; cnt++){ if ( bits_down[ cnt ] ){ phase *= -1; } }
            thisVector[ counter ] = phase * otherFCI->getFCIcoeff( bits_up , bits_down , otherVector );
         }
         if ( whichOperator=='N' ){ // Operator = particle number
            thisVector[ counter ] = bits_down[ orbIndex ] * otherVector[ counter ];
         }
      
      }

   }
   
   delete [] bits_up;
   delete [] bits_down;

}

void CheMPS2::FCI::HamTimesVec(double * input, double * output) const{

   #pragma omp parallel
   {

      int * bits_up    = new int[ L ];
      int * bits_down  = new int[ L ];
      
      const unsigned long long vecLength = getVecLength();
      
      #pragma omp for schedule(static)
      for (unsigned long long outcounter = 0; outcounter < vecLength; outcounter++){
      
         //Clear the currently considered entry of the output vector
         output[ outcounter ] = 0.0;
      
         //Find the irrep_up and irrep_down
         const int irrep_out_up   = getUpIrrepOfCounter( outcounter );
         const int irrep_out_down = getIrrepProduct( irrep_out_up , TargetIrrep );
         
         //Find the counters for the alpha and the beta electrons
         const unsigned int outcount_up   = ( outcounter - jumps[ irrep_out_up ] ) % numPerIrrep_up[ irrep_out_up ];
         const unsigned int outcount_down = ( outcounter - jumps[ irrep_out_up ] ) / numPerIrrep_up[ irrep_out_up ];
         
         //Get the bits corresponding to this counter: since we already have most steps of "getBitsOfCounter" : copy the last steps here instead of calling the function
         const unsigned int string_out_up   = cnt2str_up  [ irrep_out_up   ][ outcount_up   ];
         const unsigned int string_out_down = cnt2str_down[ irrep_out_down ][ outcount_down ];
         str2bits(L, string_out_up,   bits_up  );
         str2bits(L, string_out_down, bits_down);
         
         //The Hamiltonian is given under the form H = 0.5 * sum_{ijkl sigma tau} getHmat(i,j,k,l) a^+_{i sigma} a^+_{j tau} a_{l tau} a_{k sigma}
         //The product of the orbital irreps of the electrons which get removed should be equal to the product of the ones which are added again !!
         
         /*
            Case 1: sigma and tau are up
            The Hamiltonian can be rewritten as sum_{ i<j , k<l } ( h_ijkl - h_ijlk ) a^+_{i up} a^+_{j up} a_{l up} a_{k up}
            The Slater determinant for the down electrons does not change. As a result the irreps for the up and down determinants cannot change either!
         */
         {
            const unsigned long long incounter_offset = jumps[ irrep_out_up ] + numPerIrrep_up[ irrep_out_up ] * outcount_down;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_up[ orbi ] ){
                  bits_up[ orbi ] = 0;
                  int phase_ij = 1;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_up[ orbj ] ){
                        bits_up[ orbj ] = 0;
                        const int IrrepProd_ij = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_up[ orbk ]) ){
                              bits_up[ orbk ] = 1;
                              int phase_kl = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_up[ orbl ]) ) && ( IrrepProd_ij == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_up[ orbl ] = 1;
                                    
                                    //find the corresponding string and counter
                                    const unsigned int string_in_up = bits2str(L, bits_up);
                                    const int incount_up = str2cnt_up[ irrep_out_up ][ string_in_up ];
                                    
                                    //Add the contribution
                                    const double factor = phase_ij * phase_kl * ( getHmat(orbi, orbj, orbk, orbl) - getHmat(orbi, orbj, orbl, orbk) );
                                    output[ outcounter ] += factor * input[ incounter_offset + incount_up ];
                                    
                                    bits_up[ orbl ] = 0;
                                 }
                                 if ( bits_up[ orbl ] ){ phase_kl *= -1; }
                              }
                              bits_up[ orbk ] = 0;
                           }
                        }
                        bits_up[ orbj ] = 1;
                        phase_ij *= -1;
                     }
                  }
                  bits_up[ orbi ] = 1;
               }
            }
         
         }
         
         /* 
            Case 2: sigma and tau are down
            The Hamiltonian can be rewritten as sum_{ i<j , k<l } ( h_ijkl - h_ijlk ) a^+_{i down} a^+_{j down} a_{l down} a_{k down}
            The Slater determinant for the up electrons does not change. As a result the irreps for the up and down determinants cannot change either!
         */
         {
            const unsigned long long incounter_offset = jumps[ irrep_out_up ] + outcount_up;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_down[ orbi ] ){
                  bits_down[ orbi ] = 0;
                  int phase_ij = 1;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_down[ orbj ] ){
                        bits_down[ orbj ] = 0;
                        const int IrrepProd12 = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_down[ orbk ]) ){
                              bits_down[ orbk ] = 1;
                              int phase_kl = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_down[ orbl ]) ) && ( IrrepProd12 == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_down[ orbl ] = 1;
                                    
                                    //find the corresponding string and counter
                                    const unsigned int string_in_down = bits2str(L, bits_down);
                                    const int incount_down = str2cnt_down[ irrep_out_down ][ string_in_down ];
                                    
                                    //Add the contribution
                                    const double factor = phase_ij * phase_kl * ( getHmat(orbi, orbj, orbk, orbl) - getHmat(orbi, orbj, orbl, orbk) );
                                    output[ outcounter ] += factor * input[ incounter_offset + numPerIrrep_up[ irrep_out_up ] * incount_down ];
                                    
                                    bits_down[ orbl ] = 0;
                                 }
                                 if ( bits_down[ orbl ] ){ phase_kl *= -1; }
                              }
                              bits_down[ orbk ] = 0;
                           }
                        }
                        bits_down[ orbj ] = 1;
                        phase_ij *= -1;
                     }
                  }
                  bits_down[ orbi ] = 1;
               }
            }
            
         }
         
         /*
            Case 3: sigma and tau are different
            The Hamiltonian can be rewritten as sum_{ ijkl } h_ijkl * ( a^+_{i up} a_{k up} ) * ( a^+_{j down} a_{l down} )
         */
         {
         
            //Loop orbi to REMOVE an UP electron from the output vector
            int phase_i = 1;
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_up[ orbi ] ){
                  bits_up[ orbi ] = 0;
                  
                  //Loop orbk to ADD an UP electron again to the output vector
                  int phase_k = 1;
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if ( !(bits_up[ orbk ]) ){
                        bits_up[ orbk ] = 1;
                        const int IrrepProd_ik = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbk ] );
                        const int irrep_in_up = getIrrepProduct( irrep_out_up , IrrepProd_ik );
                        const unsigned int string_in_up = bits2str(L, bits_up);
                        const int incount_up = str2cnt_up[ irrep_in_up ][ string_in_up ];
                        const int phase_ik = phase_i * phase_k;
                        
                        //Loop orbj to REMOVE a DOWN electron from the output vector
                        int phase_j = 1;
                        for (unsigned int orbj = 0; orbj < L; orbj++){
                           if ( bits_down[ orbj ] ){
                              bits_down[ orbj ] = 0;

                              //Loop orbl to ADD a DOWN electron again to the output vector
                              int phase_l = 1;
                              for (unsigned int orbl = 0; orbl < L; orbl++){
                                 const int IrrepProd_jl = getIrrepProduct( orb2irrep[ orbj ] , orb2irrep[ orbl ] );
                                 if (( bits_down[ orbl ] == 0 ) && ( IrrepProd_ik == IrrepProd_jl )){
                                    bits_down[ orbl ] = 1;
                                    
                                    const int irrep_in_down = getIrrepProduct( irrep_out_down , IrrepProd_jl );
                                    const unsigned int string_in_down = bits2str(L, bits_down);
                                    const int incount_down = str2cnt_down[ irrep_in_down ][ string_in_down ];
                                    
                                    //Add the contribution
                                    const unsigned long long incounter = jumps[ irrep_in_up ] + incount_up + numPerIrrep_up[ irrep_in_up ] * incount_down;
                                    const double factor = phase_ik * phase_j * phase_l * getHmat(orbi, orbj, orbk, orbl);
                                    output[ outcounter ] += factor * input[ incounter ];
                                    
                                    bits_down[ orbl ] = 0;
                                 }
                                 if ( bits_down[ orbl ] ){ phase_l *= -1; }
                              }
                              bits_down[ orbj ] = 1;
                              phase_j *= -1;
                           }
                        }
                        
                        bits_up[ orbk ] = 0;
                        
                     } else { phase_k *= -1; }
                  }
                  
                  bits_up[ orbi ] = 1;
                  phase_i *= -1;
               }
            }
         
         }
      
      }
      
      delete [] bits_up;
      delete [] bits_down;
   
   }

}

double CheMPS2::FCI::Fill2RDM(double * vector, double * TwoRDM) const{

   // Clear the output 2DM
   const unsigned int Lpow4 = L*L*L*L;
   for (unsigned int loop = 0; loop < Lpow4; loop++){ TwoRDM[ loop ] = 0.0; }

   #pragma omp parallel
   {

      int * bits_up    = new int[ L ];
      int * bits_down  = new int[ L ];
      
      // Each thread gets its own local 2DM
      double * local2DM = new double[ Lpow4 ];
      for (unsigned int loop = 0; loop < Lpow4; loop++){ local2DM[ loop ] = 0.0; }
      
      // Now each thread can fill its own local 2DM, based on a slight modification of HamTimesVec
      const unsigned long long vecLength = getVecLength();
      
      #pragma omp for schedule(static)
      for (unsigned long long counter = 0; counter < vecLength; counter++){
      
         //Find the bra up and down irreps, counters, strings, and bit sequences
         const int bra_irrep_up   = getUpIrrepOfCounter( counter );
         const int bra_irrep_down = getIrrepProduct( bra_irrep_up , TargetIrrep );
         const unsigned int bra_cnt_irrep_up   = ( counter - jumps[ bra_irrep_up ] ) % numPerIrrep_up[ bra_irrep_up ];
         const unsigned int bra_cnt_irrep_down = ( counter - jumps[ bra_irrep_up ] ) / numPerIrrep_up[ bra_irrep_up ];
         const unsigned int bra_string_up   = cnt2str_up  [ bra_irrep_up   ][ bra_cnt_irrep_up   ];
         const unsigned int bra_string_down = cnt2str_down[ bra_irrep_down ][ bra_cnt_irrep_down ];
         str2bits( L , bra_string_up   , bits_up   );
         str2bits( L , bra_string_down , bits_down );
         
         // 1: Both spin projections up ==> Down Slater same ==> irreps up and down same for bra and ket
         {
            const unsigned long long ket_offset = jumps[ bra_irrep_up ] + numPerIrrep_up[ bra_irrep_up ] * bra_cnt_irrep_down;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_up[ orbi ] ){
                  bits_up[ orbi ] = 0;
                  int phase_ij = 1;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_up[ orbj ] ){
                        bits_up[ orbj ] = 0;
                        const int IrrepProd_ij = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_up[ orbk ]) ){
                              bits_up[ orbk ] = 1;
                              int phase_kl = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_up[ orbl ]) ) && ( IrrepProd_ij == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_up[ orbl ] = 1;
                                    
                                    //Find the corresponding string and counter
                                    const unsigned int ket_string_up = bits2str( L , bits_up );
                                    const int ket_cnt_irrep_up = str2cnt_up[ bra_irrep_up ][ ket_string_up ];
                                    
                                    //Add the contribution
                                    const double element = phase_ij * phase_kl * vector[ counter ] * vector[ ket_offset + ket_cnt_irrep_up ];
                                    local2DM[ orbi + L * ( orbj + L * ( orbk + L * orbl ) ) ] += element;
                                    local2DM[ orbi + L * ( orbj + L * ( orbl + L * orbk ) ) ] -= element;
                                    local2DM[ orbj + L * ( orbi + L * ( orbl + L * orbk ) ) ] += element;
                                    local2DM[ orbj + L * ( orbi + L * ( orbk + L * orbl ) ) ] -= element;
                                    
                                    bits_up[ orbl ] = 0;
                                 }
                                 if ( bits_up[ orbl ] ){ phase_kl *= -1; }
                              }
                              bits_up[ orbk ] = 0;
                           }
                        }
                        bits_up[ orbj ] = 1;
                        phase_ij *= -1;
                     }
                  }
                  bits_up[ orbi ] = 1;
               }
            }
         
         }
         
         // 2: Both spin projections down ==> Up Slater same ==> irreps up and down same for bra and ket
         {
            const unsigned long long ket_offset = jumps[ bra_irrep_up ] + bra_cnt_irrep_up;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_down[ orbi ] ){
                  bits_down[ orbi ] = 0;
                  int phase_ij = 1;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_down[ orbj ] ){
                        bits_down[ orbj ] = 0;
                        const int IrrepProd12 = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
                        for (unsigned int orbk = 0; orbk < L-1; orbk++){
                           if ( !(bits_down[ orbk ]) ){
                              bits_down[ orbk ] = 1;
                              int phase_kl = 1;
                              for (unsigned int orbl = orbk+1; orbl < L; orbl++){
                                 if (( !(bits_down[ orbl ]) ) && ( IrrepProd12 == getIrrepProduct( orb2irrep[ orbk ] , orb2irrep[ orbl ] ) )){
                                    bits_down[ orbl ] = 1;
                                    
                                    //Find the corresponding string and counter
                                    const unsigned int ket_string_down = bits2str( L , bits_down );
                                    const int ket_cnt_irrep_down = str2cnt_down[ bra_irrep_down ][ ket_string_down ];
                                    
                                    //Add the contribution
                                    const unsigned long long ketcounter = ket_offset + numPerIrrep_up[ bra_irrep_up ] * ket_cnt_irrep_down;
                                    const double element = phase_ij * phase_kl * vector[ counter ] * vector[ ketcounter ];
                                    local2DM[ orbi + L * ( orbj + L * ( orbk + L * orbl ) ) ] += element;
                                    local2DM[ orbi + L * ( orbj + L * ( orbl + L * orbk ) ) ] -= element;
                                    local2DM[ orbj + L * ( orbi + L * ( orbl + L * orbk ) ) ] += element;
                                    local2DM[ orbj + L * ( orbi + L * ( orbk + L * orbl ) ) ] -= element;
                                    
                                    bits_down[ orbl ] = 0;
                                 }
                                 if ( bits_down[ orbl ] ){ phase_kl *= -1; }
                              }
                              bits_down[ orbk ] = 0;
                           }
                        }
                        bits_down[ orbj ] = 1;
                        phase_ij *= -1;
                     }
                  }
                  bits_down[ orbi ] = 1;
               }
            }
            
         }
         
         // 3: sigma and tau are different in a^+_{i,sigma} a^+_{j,tau} a_{l,tau} a_{k,sigma}
         {
         
            //Loop orbi to REMOVE an UP electron from the output vector
            int phase_i = 1;
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_up[ orbi ] ){
                  bits_up[ orbi ] = 0;
                  
                  //Loop orbk to ADD an UP electron again to the output vector
                  int phase_k = 1;
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if ( !(bits_up[ orbk ]) ){
                        bits_up[ orbk ] = 1;
                        const int IrrepProd_ik = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbk ] );
                        const int ket_irrep_up = getIrrepProduct( bra_irrep_up , IrrepProd_ik );
                        const unsigned int ket_string_up = bits2str( L , bits_up );
                        const int ket_cnt_irrep_up = str2cnt_up[ ket_irrep_up ][ ket_string_up ];
                        const int phase_ik = phase_i * phase_k;
                        
                        //Loop orbj to REMOVE a DOWN electron from the output vector
                        int phase_j = 1;
                        for (unsigned int orbj = 0; orbj < L; orbj++){
                           if ( bits_down[ orbj ] ){
                              bits_down[ orbj ] = 0;

                              //Loop orbl to ADD a DOWN electron again to the output vector
                              int phase_l = 1;
                              for (unsigned int orbl = 0; orbl < L; orbl++){
                                 const int IrrepProd_jl = getIrrepProduct( orb2irrep[ orbj ] , orb2irrep[ orbl ] );
                                 if (( bits_down[ orbl ] == 0 ) && ( IrrepProd_ik == IrrepProd_jl )){
                                    bits_down[ orbl ] = 1;
                                    
                                    const int ket_irrep_down = getIrrepProduct( bra_irrep_down , IrrepProd_jl );
                                    const unsigned int ket_string_down = bits2str( L , bits_down );
                                    const int ket_cnt_irrep_down = str2cnt_down[ ket_irrep_down ][ ket_string_down ];
                                    
                                    //Add the contribution
                                    const unsigned long long ketcounter = jumps[ ket_irrep_up ] + ket_cnt_irrep_up
                                                                        + numPerIrrep_up[ ket_irrep_up ] * ket_cnt_irrep_down;
                                    const double element = phase_ik * phase_j * phase_l * vector[ counter ] * vector[ ketcounter ];
                                    local2DM[ orbi + L * ( orbj + L * ( orbk + L * orbl ) ) ] += element;
                                    local2DM[ orbj + L * ( orbi + L * ( orbl + L * orbk ) ) ] += element;
                                    
                                    bits_down[ orbl ] = 0;
                                 }
                                 if ( bits_down[ orbl ] ){ phase_l *= -1; }
                              }
                              bits_down[ orbj ] = 1;
                              phase_j *= -1;
                           }
                        }
                        
                        bits_up[ orbk ] = 0;
                        
                     } else { phase_k *= -1; }
                  }
                  
                  bits_up[ orbi ] = 1;
                  phase_i *= -1;
               }
            }
         
         }
      
      }
      
      delete [] bits_up;
      delete [] bits_down;
      
      #pragma omp critical
      {
         for (unsigned int loop = 0; loop < Lpow4; loop++){ TwoRDM[ loop ] += local2DM[ loop ]; }
      }
      
      delete [] local2DM;
   
   }
   
   // Calculate the FCI energy
   double FCIenergy = getEconst();
   for (unsigned int loop = 0; loop < Lpow4; loop++){ FCIenergy += 0.5 * TwoRDM[ loop ] * Hmat[ loop ]; }
   if ( FCIverbose > 0 ){ cout << "FCI::Fill2DM : Econstant + 0.5 * trace( 2-RDM * Hmat ) = " << FCIenergy << endl; }
   return FCIenergy;

}

double CheMPS2::FCI::CalcSpinSquared(double * vector) const{

   // The total values
   double totalSzSz   = 0.0;
   double totalS2rest = 0.0;

   #pragma omp parallel
   {

      int * bits_up    = new int[ L ];
      int * bits_down  = new int[ L ];
      
      // Each thread gets its own variables
      double mySzSz   = 0.0;
      double myS2rest = 0.0;
      
      const unsigned long long vecLength = getVecLength();
      
      #pragma omp for schedule(static)
      for (unsigned long long counter = 0; counter < vecLength; counter++){
      
         getBitsOfCounter( counter , bits_up , bits_down );
         const double vector_counter_squared = vector[ counter ] * vector[ counter ];
         
         //First the diagonal parts
         for (unsigned int orbi = 0; orbi < L; orbi++){
            const double tempvar = ( bits_up[ orbi ] + bits_down[ orbi ] - 2 * bits_up[ orbi ] * bits_down[ orbi ] ) * vector_counter_squared;
            mySzSz   += 0.25 * tempvar;
            myS2rest += 0.5  * tempvar;
            for (unsigned int orbj = orbi+1; orbj < L; orbj++){
               mySzSz += 0.5 * ( bits_up[ orbi ] - bits_down[ orbi ] ) * ( bits_up[ orbj ] - bits_down[ orbj ] ) * vector_counter_squared;
            }
         }
         
         //And then the off-diagonal parts
         for (unsigned int orbi = 0; orbi < L-1; orbi++){
         
            // a_i,up^+ a_i,down a_j,down^+ a_j,up = a_i,up^+ a_j,up a_i,down a_j,down^+
            if (( bits_up[ orbi ] ) && ( !(bits_down[ orbi ]) )){
               bits_up  [ orbi ] = 0;
               bits_down[ orbi ] = 1;
               int phase = -1; // For a_j,down^+ passing a_i,down^+
               for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                  if (( bits_down[ orbj ] ) && ( !(bits_up[ orbj ]) )){
                     bits_up  [ orbj ] = 1;
                     bits_down[ orbj ] = 0;
                     
                     myS2rest += phase * vector[ counter ] * getFCIcoeff( bits_up , bits_down , vector );
                     
                     bits_up  [ orbj ] = 0;
                     bits_down[ orbj ] = 1;
                  }
                  if ( bits_up  [ orbj ] ){ phase *= -1; }
                  if ( bits_down[ orbj ] ){ phase *= -1; }
               }
               bits_up  [ orbi ] = 1;
               bits_down[ orbi ] = 0;
            }

            // a_i,down^+ a_i,up a_j,up^+ a_j,down = a_i,up a_j,up^+ a_i,down^+ a_j,down
            if (( !(bits_up[ orbi ]) ) && ( bits_down[ orbi ] )){
               bits_up  [ orbi ] = 1;
               bits_down[ orbi ] = 0;
               int phase = -1; // For a_j,up^+ passing a_i,up^+
               for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                  if (( bits_up[ orbj ] ) && ( !(bits_down[ orbj ]) )){
                     bits_up  [ orbj ] = 0;
                     bits_down[ orbj ] = 1;
                     
                     myS2rest += phase * vector[ counter ] * getFCIcoeff( bits_up , bits_down , vector );
                     
                     bits_up  [ orbj ] = 1;
                     bits_down[ orbj ] = 0;
                  }
                  if ( bits_up  [ orbj ] ){ phase *= -1; }
                  if ( bits_down[ orbj ] ){ phase *= -1; }
               }
               bits_up  [ orbi ] = 0;
               bits_down[ orbi ] = 1;
            }
            
         }
      }
      
      #pragma omp critical
      {
         totalSzSz   += mySzSz;
         totalS2rest += myS2rest;
      }
      
      delete [] bits_up;
      delete [] bits_down;
   
   }
   
   const double calcS_squared = totalS2rest + totalSzSz;
   if ( FCIverbose > 0 ){
      const double intendedS         = 0.5 * ( getNel_up() - getNel_down() );
      const double intendedS_squared = intendedS * ( intendedS + 1.0 );
      cout << "FCI::CalcSpinSquared : For intended spin " << intendedS << " the measured S(S+1) = " << calcS_squared
                                                                       << " and intented S(S+1) = " << intendedS_squared << endl;
   }
   return calcS_squared;

}

void CheMPS2::FCI::DiagHam(double * diag) const{

   #pragma omp parallel
   {

      int * bits_up   = new int[ L ];
      int * bits_down = new int[ L ];
      
      const unsigned long long vecLength = getVecLength();
      
      #pragma omp for schedule(static)
      for (unsigned long long counter = 0; counter < vecLength; counter++){
      
         diag[ counter ] = 0.0; // Clear the current entry
         getBitsOfCounter( counter , bits_up , bits_down ); // Fetch the corresponding bits
         
         // Diag( Ham ) = sum_{i} h_{iiii} n_i,up n_i,down + sum_{i<j} h_{ijij} n_i,TOT n_j,TOT - sum_{i<j} h_{ijji} ( n_i,up n_j,up + n_i,down n_j,down )
         for (unsigned int orb1 = 0; orb1 < L; orb1++){
            diag[ counter ] += bits_up[ orb1 ] * bits_down[ orb1 ] * getHmat( orb1 , orb1 , orb1 , orb1 );
            for (unsigned int orb2 = orb1+1; orb2 < L; orb2++){
               diag[ counter ] += ( bits_up[ orb1 ] + bits_down[ orb1 ] ) * ( bits_up[ orb2 ] + bits_down[ orb2 ] ) * getHmat( orb1 , orb2 , orb1 , orb2 );
               diag[ counter ] -= ( bits_up[ orb1 ] * bits_up[ orb2 ] + bits_down[ orb1 ] * bits_down[ orb2 ] ) * getHmat( orb1 , orb2 , orb2 , orb1 );
            
            }
         }
         
      }
      
      delete [] bits_up;
      delete [] bits_down;
      
   }

}

void CheMPS2::FCI::DiagHamSquared(double * output) const{
   
   /*
      Completeness relation strategy to calculate the diagonal elements of Ham*Ham:
      Loop over the FCI determinants |det>
         < det | Ham * Ham | det > = 0.0
         Loop over the FCI determinants which are connected with |det> by Ham |det2>
            < det | Ham * Ham | det > += ( < det2 | Ham | det > )^2
   */
   
   #pragma omp parallel
   {
   
      const unsigned long long vecLength = getVecLength();

      int * bits_up_ket    = new int[ L ];
      int * bits_down_ket  = new int[ L ];
      int * bits_up_bra    = new int[ L ];
      int * bits_down_bra  = new int[ L ];
      int * work1 = new int[ 2 ];
      int * work2 = new int[ 2 ];
      int * work3 = new int[ 2 ];
      int * work4 = new int[ 2 ];
      bool * workspace = new bool[ vecLength ];
      
      #pragma omp for schedule(static)
      for (unsigned long long outcounter = 0; outcounter < vecLength; outcounter++){
      
         //Clear the currently considered entry of the output vector
         output[ outcounter ] = 0.0;
         
         //Workspace keeps track that a det2 is used only once for the completeness relation
         for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ workspace[ cnt ] = true; }
      
         //Find the irrep_up and irrep_down
         const int irrep_out_up   = getUpIrrepOfCounter( outcounter );
         const int irrep_out_down = getIrrepProduct( irrep_out_up , TargetIrrep );
         
         //Find the counters for the alpha and the beta electrons
         const unsigned int outcount_up   = ( outcounter - jumps[ irrep_out_up ] ) % numPerIrrep_up[ irrep_out_up ];
         const unsigned int outcount_down = ( outcounter - jumps[ irrep_out_up ] ) / numPerIrrep_up[ irrep_out_up ];
         
         //Get the bits corresponding to this counter: since we already have most steps of "getBitsOfCounter" : copy the last steps here instead of calling the function
         const unsigned int string_out_up   = cnt2str_up  [ irrep_out_up   ][ outcount_up   ];
         const unsigned int string_out_down = cnt2str_down[ irrep_out_down ][ outcount_down ];
         str2bits(L, string_out_up,   bits_up_ket  );
         str2bits(L, string_out_down, bits_down_ket);
         
         for (unsigned int cnt = 0; cnt < L; cnt++){
            bits_up_bra  [ cnt ] = bits_up_ket  [ cnt ];
            bits_down_bra[ cnt ] = bits_down_ket[ cnt ];
         }
         
         //The Hamiltonian is given under the form H = 0.5 * sum_{ijkl sigma tau} getHmat(i,j,k,l) a^+_{i sigma} a^+_{j tau} a_{l tau} a_{k sigma}
         //The product of the orbital irreps of the electrons which get removed should be equal to the product of the ones which are added again !!
         
         // Case 1: sigma and tau are up
         {
            const unsigned long long incounter_offset = jumps[ irrep_out_up ] + numPerIrrep_up[ irrep_out_up ] * outcount_down;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_up_bra[ orbi ] ){
                  bits_up_bra[ orbi ] = 0;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_up_bra[ orbj ] ){
                        bits_up_bra[ orbj ] = 0;
                        const int IrrepProd_ij = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
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
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work1, work2, work3, work4);
                                       output[ outcounter ] += factor * factor;
                                       workspace[ incounter ] = false;
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
         
         // Case 2: sigma and tau are down
         {
            const unsigned long long incounter_offset = jumps[ irrep_out_up ] + outcount_up;
            
            //Loop i<j to remove from the output vector
            for (unsigned int orbi = 0; orbi < L-1; orbi++){
               if ( bits_down_bra[ orbi ] ){
                  bits_down_bra[ orbi ] = 0;
                  for (unsigned int orbj = orbi+1; orbj < L; orbj++){
                     if ( bits_down_bra[ orbj ] ){
                        bits_down_bra[ orbj ] = 0;
                        const int IrrepProd12 = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbj ] );
                        
                        //Loop k<l to add again to the output vector
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
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work1, work2, work3, work4);
                                       output[ outcounter ] += factor * factor;
                                       workspace[ incounter ] = false;
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
         
         // Case 3: sigma and tau are different
         {
         
            //Loop orbi to REMOVE an UP electron from the output vector
            for (unsigned int orbi = 0; orbi < L; orbi++){
               if ( bits_up_bra[ orbi ] ){
                  bits_up_bra[ orbi ] = 0;
                  
                  //Loop orbk to ADD an UP electron again to the output vector
                  for (unsigned int orbk = 0; orbk < L; orbk++){
                     if ( !(bits_up_bra[ orbk ]) ){
                        bits_up_bra[ orbk ] = 1;
                        const int IrrepProd_ik = getIrrepProduct( orb2irrep[ orbi ] , orb2irrep[ orbk ] );
                        const int irrep_in_up = getIrrepProduct( irrep_out_up , IrrepProd_ik );
                        const unsigned int string_in_up = bits2str(L, bits_up_bra);
                        const int incount_up = str2cnt_up[ irrep_in_up ][ string_in_up ];
                        
                        //Loop orbj to REMOVE a DOWN electron from the output vector
                        for (unsigned int orbj = 0; orbj < L; orbj++){
                           if ( bits_down_bra[ orbj ] ){
                              bits_down_bra[ orbj ] = 0;

                              //Loop orbl to ADD a DOWN electron again to the output vector
                              for (unsigned int orbl = 0; orbl < L; orbl++){
                                 const int IrrepProd_jl = getIrrepProduct( orb2irrep[ orbj ] , orb2irrep[ orbl ] );
                                 if (( bits_down_bra[ orbl ] == 0 ) && ( IrrepProd_ik == IrrepProd_jl )){
                                    bits_down_bra[ orbl ] = 1;
                                    
                                    const int irrep_in_down = getIrrepProduct( irrep_out_down , IrrepProd_jl );
                                    const unsigned int string_in_down = bits2str(L, bits_down_bra);
                                    const int incount_down = str2cnt_down[ irrep_in_down ][ string_in_down ];
                                    const unsigned long long incounter = jumps[ irrep_in_up ] + incount_up + numPerIrrep_up[ irrep_in_up ] * incount_down;
                                    
                                    if ( workspace[ incounter ] ){
                                       const double factor = GetMatrixElement(bits_up_bra, bits_down_bra, bits_up_ket, bits_down_ket, work1, work2, work3, work4);
                                       output[ outcounter ] += factor * factor;
                                       workspace[ incounter ] = false;
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
      
      }
      
      delete [] bits_up_ket;
      delete [] bits_down_ket;
      delete [] bits_up_bra;
      delete [] bits_down_bra;
      delete [] work1;
      delete [] work2;
      delete [] work3;
      delete [] work4;
      delete [] workspace;
   
   }

}


double CheMPS2::FCI::GSSmallCISpace(const unsigned int Nslaters, double * vec) const{

   // Just make sure that we don't do anything stupid for small FCI vector spaces
   const unsigned long long vecLength = getVecLength();
   int numSlaters = ( Nslaters > vecLength ) ? vecLength : Nslaters;

   // Fetch the Slater determinant energies (diagonal elements of the FCI Hamiltonian)
   DiagHam( vec );
   
   // Find the maximum Slater determinant energy
   double max_E = vec[0];
   for (unsigned long long count = 1; count < vecLength; count++){ if ( vec[ count ] > max_E ){ max_E = vec[ count ]; } }
   
   // Find the lowest energy determinants; this destroys vec
   unsigned long long * LowEnergyDeterminants = new unsigned long long[ numSlaters ];
   for (int det = 0; det < numSlaters; det++){
      unsigned long long ID = 0;
      for (unsigned long long count = 1; count < vecLength; count++){
         if ( vec[ count ] < vec[ ID ] ){ ID = count; }
      }
      LowEnergyDeterminants[ det ] = ID;
      vec[ ID ] = max_E + 1.0;
   }
   
   // Construct the CI Hamiltonian
   double * CIham = new double[ numSlaters * numSlaters ];
   int * bits_ket_up   = new int[ L ];
   int * bits_ket_down = new int[ L ];
   int * bits_bra_up   = new int[ L ];
   int * bits_bra_down = new int[ L ];
   int * work1 = new int[ 2 ];
   int * work2 = new int[ 2 ];
   int * work3 = new int[ 2 ];
   int * work4 = new int[ 2 ];
   for (int row = 0; row < numSlaters; row++){
      getBitsOfCounter( LowEnergyDeterminants[ row ] , bits_bra_up , bits_bra_down );
      for (int col = row; col < numSlaters; col++){
         getBitsOfCounter( LowEnergyDeterminants[ col ] , bits_ket_up , bits_ket_down );
         CIham[ row + numSlaters * col ] = GetMatrixElement( bits_bra_up , bits_bra_down , bits_ket_up , bits_ket_down , work1 , work2 , work3 , work4 );
         CIham[ col + numSlaters * row ] = CIham[ row + numSlaters * col ];
      }
   }
   delete [] work1;
   delete [] work2;
   delete [] work3;
   delete [] work4;
   delete [] bits_ket_up;
   delete [] bits_ket_down;
   delete [] bits_bra_up;
   delete [] bits_bra_down;
   
   // Diagonalize CIham: use double * vec to store the eigenvalues
   char jobz = 'V';
   char uplo = 'U';
   int lwork = 3 * numSlaters - 1;
   double * work = new double[ lwork ];
   int info;
   dsyev_( &jobz , &uplo , &numSlaters , CIham , &numSlaters , vec , work , &lwork , &info );
   delete [] work;
   const double CIenergy = vec[0] + getEconst();
   if ( FCIverbose > 0 ){ cout << "FCI::GSSmallCISpace : Small CI space (" << numSlaters << " determinants) diagonalization yields energy = " << CIenergy << endl; }
   
   // Copy lowest energy vector
   FCIdclear( vecLength , vec );
   for (int det = 0; det < numSlaters; det++){
      vec[ LowEnergyDeterminants[ det ] ] = CIham[ det ];
   }
   
   // Clear memory
   delete [] CIham;
   delete [] LowEnergyDeterminants;
   
   return CIenergy;

}

double CheMPS2::FCI::GetMatrixElement(int * bits_bra_up, int * bits_bra_down, int * bits_ket_up, int * bits_ket_down, int * annih_up, int * creat_up, int * annih_down, int * creat_down) const{
   
   int count_annih_up   = 0;
   int count_creat_up   = 0;
   int count_annih_down = 0;
   int count_creat_down = 0;
   
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
   
   // Diag( Ham ) = sum_i h_{iiii} n_i,up n_i,down + sum_{i<j} h_{ijij} n_i,TOT n_j,TOT - sum_{i<j} h_{ijji} ( n_i,up n_j,up + n_i,down n_j,down )
   if (( count_annih_up == 0 ) && ( count_annih_down == 0 )){ // |bra> == |ket>
   
      double result = 0.0;
      for (unsigned int orb1 = 0; orb1 < L; orb1++){
         result += getHmat(orb1, orb1, orb1, orb1) * bits_ket_up[ orb1 ] * bits_ket_down[ orb1 ];
         for (unsigned int orb2 = orb1+1; orb2 < L; orb2++){
            result += ( bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ] ) * ( bits_ket_up[ orb2 ] + bits_ket_down[ orb2 ] ) * getHmat(orb1, orb2, orb1, orb2);
            result -= ( bits_ket_up[ orb1 ] * bits_ket_up[ orb2 ] + bits_ket_down[ orb1 ] * bits_ket_down[ orb2 ] ) * getHmat(orb1, orb2, orb2, orb1);
         }
      }
      return result;
      
   }
   
   // sum_{i ; j!=l} sum_{sigma} ( h_{ijil} - h_{ijli} delta_{sigma,up} ) a^+_j,up a_l,up n_i,sigma
   if (( count_annih_up == 1 ) && ( count_annih_down == 0 )){ // |bra> = a^+_j,up a_l,up |ket>
   
      const int orbj = creat_up[ 0 ];
      const int orbl = annih_up[ 0 ];
   
      double result = 0.0;
      for (unsigned int orb1 = 0; orb1 < L; orb1++){
         result += getHmat(orb1, orbj, orb1, orbl) * ( bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ] )
                 - getHmat(orb1, orbj, orbl, orb1) * bits_ket_up[ orb1 ];
      }
      int phase = 1;
      if ( orbj < orbl ){ for (int orbital = orbj+1; orbital < orbl; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } }
      if ( orbl < orbj ){ for (int orbital = orbl+1; orbital < orbj; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } }
      return ( result * phase );
   
   }
   
   // sum_{i ; j!=l} sum_{sigma} ( h_{ijil} - h_{ijli} delta_{sigma,down} ) a^+_j,down a_l,down n_i,sigma
   if (( count_annih_up == 0 ) && ( count_annih_down == 1 )){ // |bra> = a^+_j,down a_l,down |ket>
   
      const int orbj = creat_down[ 0 ];
      const int orbl = annih_down[ 0 ];
   
      double result = 0.0;
      for (unsigned int orb1 = 0; orb1 < L; orb1++){
         result += getHmat(orb1, orbj, orb1, orbl) * ( bits_ket_up[ orb1 ] + bits_ket_down[ orb1 ] )
                 - getHmat(orb1, orbj, orbl, orb1) * bits_ket_down[ orb1 ];
      }
      int phase = 1;
      if ( orbj < orbl ){ for (int orbital = orbj+1; orbital < orbl; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      if ( orbl < orbj ){ for (int orbital = orbl+1; orbital < orbj; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } }
      return ( result * phase );
   
   }
   
   // sum_{i < j ; k < l} ( h_{ijkl} - h_{ijlk} ) a_^+_i,up a^+_j,up a_l,up a_k,up
   if (( count_annih_up == 2 ) && ( count_annih_down == 0 )){
   
      // creat and annih are filled in increasing orbital index
      const int orbi = creat_up[ 0 ];
      const int orbj = creat_up[ 1 ];
      const int orbk = annih_up[ 0 ];
      const int orbl = annih_up[ 1 ];
      
      double result = getHmat(orbi, orbj, orbk, orbl) - getHmat(orbi, orbj, orbl, orbk);
      int phase = 1;
      for (int orbital = orbk+1; orbital < orbl; orbital++){ if ( bits_ket_up[ orbital ] ){ phase *= -1; } } // Fermion phases orbk and orbl measured in the ket
      for (int orbital = orbi+1; orbital < orbj; orbital++){ if ( bits_bra_up[ orbital ] ){ phase *= -1; } } // Fermion phases orbi and orbj measured in the bra
      return ( result * phase );
   
   }
   
   // sum_{i < j ; k < l} ( h_{ijkl} - h_{ijlk} ) a_^+_i,down a^+_j,down a_l,down a_k,down
   if (( count_annih_up == 0 ) && ( count_annih_down == 2 )){
   
      // creat and annih are filled in increasing orbital index
      const int orbi = creat_down[ 0 ];
      const int orbj = creat_down[ 1 ];
      const int orbk = annih_down[ 0 ];
      const int orbl = annih_down[ 1 ];
      
      double result = getHmat(orbi, orbj, orbk, orbl) - getHmat(orbi, orbj, orbl, orbk);
      int phase = 1;
      for (int orbital = orbk+1; orbital < orbl; orbital++){ if ( bits_ket_down[ orbital ] ){ phase *= -1; } } // Fermion phases orbk and orbl measured in the ket
      for (int orbital = orbi+1; orbital < orbj; orbital++){ if ( bits_bra_down[ orbital ] ){ phase *= -1; } } // Fermion phases orbi and orbj measured in the bra
      return ( result * phase );
   
   }
   
   // sum_{ijkl} h_{ijkl} ( a_^+_i,up a_k,up ) ( a^+_j,down a_l,down )
   if (( count_annih_up == 1 ) && ( count_annih_down == 1 )){
   
      const int orbi = creat_up  [ 0 ];
      const int orbj = creat_down[ 0 ];
      const int orbk = annih_up  [ 0 ];
      const int orbl = annih_down[ 0 ];
      
      double result = getHmat(orbi, orbj, orbk, orbl);
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

   // TODO : Use lapack under the hood

   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ target[cnt] = origin[cnt]; }

}

double CheMPS2::FCI::FCIddot(const unsigned long long vecLength, double * vec1, double * vec2){

   // TODO : Use lapack under the hood

   double result = 0.0;
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ result += vec1[cnt] * vec2[cnt]; }
   return result;

}

double CheMPS2::FCI::FCIfrobeniusnorm(const unsigned long long vecLength, double * vec){

   return sqrt( FCIddot( vecLength , vec , vec ) );

}

void CheMPS2::FCI::FCIdaxpy(const unsigned long long vecLength, const double alpha, double * vec_x, double * vec_y){

   // TODO : Use lapack under the hood

   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ vec_y[cnt] += alpha * vec_x[cnt]; }

}

void CheMPS2::FCI::FCIdscal(const unsigned long long vecLength, const double alpha, double * vec){

   // TODO : Use lapack under the hood

   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ vec[cnt] = alpha * vec[cnt]; }

}

void CheMPS2::FCI::FCIdclear(const unsigned long long vecLength, double * vec){

   // TODO : Use lapack under the hood

   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ vec[cnt] = 0.0; }

}

void CheMPS2::FCI::FillRandom(const unsigned long long vecLength, double * vec){

   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ vec[cnt] = ( 2.0 * rand() - 1.0 ) / RAND_MAX; }

}

double CheMPS2::FCI::GSDavidson(double * inoutput) const{

   const int DAVIDSON_NUM_VEC           = CheMPS2::HEFF_DAVIDSON_NUM_VEC;
   const int DAVIDSON_NUM_VEC_KEEP      = CheMPS2::HEFF_DAVIDSON_NUM_VEC_KEEP;
   const double DAVIDSON_RTOL_BASE      = CheMPS2::HEFF_DAVIDSON_RTOL_BASE;
   const double DAVIDSON_PRECOND_CUTOFF = CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF;

   const unsigned long long length_vec = getVecLength();
   int num_vec = 0;
   double ** vecs  = new double*[DAVIDSON_NUM_VEC];
   double ** Hvecs = new double*[DAVIDSON_NUM_VEC];
   int num_allocated = 0;
   
   double * mxM = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   double * mxM_eigs = new double[DAVIDSON_NUM_VEC];
   double * mxM_vecs = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   int mxM_lwork = 3*DAVIDSON_NUM_VEC-1;
   double * mxM_work = new double[mxM_lwork];
   
   double rtol = DAVIDSON_RTOL_BASE * sqrt( 1.0 * length_vec );
   double rnorm = 10*rtol;
   
   double * t_vec = new double[length_vec];
   double * u_vec = new double[length_vec];
   double * work_vec = new double[length_vec];
   int inc1 = 1;
   
   double Norm = 0.0;
   if ( inoutput != NULL ){
      Norm = FCIfrobeniusnorm(length_vec, inoutput);
      if ( Norm != 0.0 ){ FCIdcopy(length_vec, inoutput, t_vec); }
   }
   if (( Norm == 0.0 ) || ( inoutput == NULL )){
      FillRandom(length_vec, t_vec);
   }
   
   double * Diag = new double[length_vec];
   DiagHam(Diag);
   
   double * Reortho_Lowdin = NULL;
   double * Reortho_Overlap_eigs = NULL;
   double * Reortho_Overlap = NULL;
   double * Reortho_Eigenvecs = NULL;
   bool Reortho_Allocated = false;
   
   int nIterations = 0;
   
   while (rnorm > rtol){

      //1. Orthogonalize the new t_vec w.r.t. the old basis
      for (int cnt=0; cnt<num_vec; cnt++){
         const double minus_overlap = - FCIddot(length_vec, t_vec, vecs[cnt]);
         FCIdaxpy(length_vec, minus_overlap, vecs[cnt], t_vec);
      }
   
      //2. Normalize the t_vec
      double alpha = 1.0/FCIfrobeniusnorm(length_vec, t_vec);
      FCIdscal(length_vec, alpha, t_vec);
      
      //3. t_vec becomes part of vecs
      if (num_vec<num_allocated){
         double * temp = vecs[num_vec];
         vecs[num_vec] = t_vec;
         t_vec = temp;
      } else {
         vecs[num_allocated] = t_vec;
         Hvecs[num_allocated] = new double[length_vec];
         t_vec = new double[length_vec];
         num_allocated++;
      }
      HamTimesVec(vecs[num_vec], Hvecs[num_vec]);
      nIterations++;
      
      //4. mxM contains the Hamiltonian in the basis "vecs"
      for (int cnt=0; cnt<num_vec; cnt++){
         mxM[cnt + DAVIDSON_NUM_VEC * num_vec] = FCIddot(length_vec, vecs[num_vec], Hvecs[cnt]);
         mxM[num_vec + DAVIDSON_NUM_VEC * cnt] = mxM[cnt + DAVIDSON_NUM_VEC * num_vec];
      }
      mxM[num_vec + DAVIDSON_NUM_VEC * num_vec] = FCIddot(length_vec, vecs[num_vec], Hvecs[num_vec]);
      
      //5. When t_vec was added to vecs, the number of vecs was actually increased by one. For convenience (doing 4.), only now the number is incremented.
      num_vec++;
      
      //6. Calculate the eigenvalues and vectors of mxM
      char jobz = 'V';
      char uplo = 'U';
      int info;
      for (int cnt1=0; cnt1<num_vec; cnt1++){
         for (int cnt2=0; cnt2<num_vec; cnt2++){
            mxM_vecs[cnt1 + DAVIDSON_NUM_VEC * cnt2] = mxM[cnt1 + DAVIDSON_NUM_VEC * cnt2];
         }
      }
      int lda = DAVIDSON_NUM_VEC;
      dsyev_(&jobz,&uplo,&num_vec,mxM_vecs,&lda,mxM_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
      if ( FCIverbose > 1 ){ cout << "FCI::GSDavidson : Current lowest eigenvalue = " << mxM_eigs[0] + getEconst() << endl; }
      
      //7. Calculate u and r. r is stored in t_vec, u in u_vec.
      FCIdclear(length_vec, t_vec);
      FCIdclear(length_vec, u_vec);
      for (int cnt=0; cnt<num_vec; cnt++){
         double alpha = mxM_vecs[cnt]; //eigenvector with lowest eigenvalue, hence mxM_vecs[cnt + DAVIDSON_NUM_VEC * 0]
         FCIdaxpy(length_vec, alpha, Hvecs[cnt], t_vec);
         FCIdaxpy(length_vec, alpha,  vecs[cnt], u_vec);
      }
      alpha = -mxM_eigs[0];
      FCIdaxpy(length_vec, alpha, u_vec, t_vec);
      
      //8. Calculate the norm of r
      rnorm = FCIfrobeniusnorm(length_vec, t_vec);
      
      //9. In case convergence is not yet reached: prepare for the following iteration
      if (rnorm > rtol){
      
         //9a. Calculate the new t_vec based on the residual of the lowest eigenvalue, to add to the vecs.
         for (unsigned long long cnt=0; cnt<length_vec; cnt++){
            if (fabs(Diag[cnt] - mxM_eigs[0])> DAVIDSON_PRECOND_CUTOFF ){
               work_vec[cnt] = u_vec[cnt]/(Diag[cnt] - mxM_eigs[0]); // work_vec = K^(-1) u_vec
            } else {
               work_vec[cnt] = u_vec[cnt]/ DAVIDSON_PRECOND_CUTOFF;
            }
         }
         alpha = - FCIddot(length_vec, work_vec, t_vec) / FCIddot(length_vec, work_vec, u_vec); // alpha = - (u^T K^(-1) r) / (u^T K^(-1) u)
         FCIdaxpy(length_vec, alpha, u_vec, t_vec); // t_vec = r - (u^T K^(-1) r) / (u^T K^(-1) u) u
         for (unsigned long long cnt=0; cnt<length_vec; cnt++){
            if (fabs(Diag[cnt] - mxM_eigs[0])> DAVIDSON_PRECOND_CUTOFF ){
               t_vec[cnt] = - t_vec[cnt]/(Diag[cnt] - mxM_eigs[0]); //t_vec = - K^(-1) (r - (u^T K^(-1) r) / (u^T K^(-1) u) u)
            } else {
               t_vec[cnt] = - t_vec[cnt]/ DAVIDSON_PRECOND_CUTOFF ;
            }
         }
         
         // 9b. When the maximum number of vectors is reached: construct the one with lowest eigenvalue & restart
         if (num_vec == DAVIDSON_NUM_VEC){
         
            if (DAVIDSON_NUM_VEC_KEEP<=1){
            
               alpha = 1.0/FCIfrobeniusnorm(length_vec, u_vec);
               FCIdscal(length_vec, alpha, u_vec);
               FCIdcopy(length_vec, u_vec, vecs[0]);
               HamTimesVec(vecs[0], Hvecs[0]);
               nIterations++;
               mxM[0] = FCIddot(length_vec, vecs[0], Hvecs[0]);
            
               num_vec = 1;
            
            } else {
            
               //Construct the lowest DAVIDSON_NUM_VEC_KEEP eigenvectors
               if (!Reortho_Allocated) Reortho_Eigenvecs = new double[length_vec * DAVIDSON_NUM_VEC_KEEP];
               FCIdcopy(length_vec, u_vec, Reortho_Eigenvecs);
               for (int cnt=1; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  for (unsigned long long irow=0; irow<length_vec; irow++){
                     Reortho_Eigenvecs[irow + length_vec * cnt] = 0.0;
                     for (int ivec=0; ivec<DAVIDSON_NUM_VEC; ivec++){
                        Reortho_Eigenvecs[irow + length_vec * cnt] += vecs[ivec][irow] * mxM_vecs[ivec + DAVIDSON_NUM_VEC * cnt];
                     }
                  }
               }
               
               //Reorthonormalize them
               //Reortho: Calculate the overlap matrix
               if (!Reortho_Allocated) Reortho_Overlap = new double[DAVIDSON_NUM_VEC_KEEP * DAVIDSON_NUM_VEC_KEEP];
               for (int row=0; row<DAVIDSON_NUM_VEC_KEEP; row++){
                  for (int col=row; col<DAVIDSON_NUM_VEC_KEEP; col++){
                     Reortho_Overlap[row + DAVIDSON_NUM_VEC_KEEP * col] = FCIddot(length_vec, Reortho_Eigenvecs + row * length_vec, Reortho_Eigenvecs + col * length_vec);
                     Reortho_Overlap[col + DAVIDSON_NUM_VEC_KEEP * row] = Reortho_Overlap[row + DAVIDSON_NUM_VEC_KEEP * col];
                  }
               }
               
               //Reortho: Calculate the Lowdin tfo
               int DVDS_KEEP = DAVIDSON_NUM_VEC_KEEP;
               if (!Reortho_Allocated) Reortho_Overlap_eigs = new double[DVDS_KEEP];
               dsyev_(&jobz,&uplo,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
               for (int icnt=0; icnt<DVDS_KEEP; icnt++){
                  Reortho_Overlap_eigs[icnt] = pow(Reortho_Overlap_eigs[icnt],-0.25);
                  dscal_(&DVDS_KEEP, Reortho_Overlap_eigs+icnt, Reortho_Overlap+DVDS_KEEP*icnt, &inc1);
               }
               if (!Reortho_Allocated) Reortho_Lowdin = new double[DVDS_KEEP*DVDS_KEEP];
               char trans = 'T';
               char notr = 'N';
               double one = 1.0;
               double zero = 0.0; //set
               dgemm_(&notr,&trans,&DVDS_KEEP,&DVDS_KEEP,&DVDS_KEEP,&one,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,&zero,Reortho_Lowdin,&DVDS_KEEP);
               
               //Reortho: Put the Lowdin tfo eigenvecs in vecs
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  FCIdclear(length_vec, vecs[ivec]);
                  for (int ivec2=0; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     FCIdaxpy(length_vec, Reortho_Lowdin[ ivec2 + DAVIDSON_NUM_VEC_KEEP * ivec ], Reortho_Eigenvecs + length_vec * ivec2, vecs[ivec]);
                  }
               }
               
               if (!Reortho_Allocated) Reortho_Allocated = true;
               
               //Construct the H*vecs
               for (int cnt=0; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  HamTimesVec(vecs[cnt], Hvecs[cnt]);
                  nIterations++;
               }
               
               //Build MxM
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  for (int ivec2=ivec; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     mxM[ivec + DAVIDSON_NUM_VEC * ivec2] = FCIddot(length_vec, vecs[ivec], Hvecs[ivec2]);
                     mxM[ivec2 + DAVIDSON_NUM_VEC * ivec] = mxM[ivec + DAVIDSON_NUM_VEC * ivec2];
                  }
               }
               
               //Set num_vec
               num_vec = DAVIDSON_NUM_VEC_KEEP;
            
            }

         }
      }
      
   }
   
   const double FCIenergy = mxM_eigs[0] + getEconst();
   if ( inoutput != NULL ){ FCIdcopy(length_vec, u_vec, inoutput); }
   
   if ( FCIverbose > 1 ){ cout << "FCI::GSDavidson : Number of matrix-vector products which were required = " << nIterations << endl; }
   if ( FCIverbose > 0 ){ cout << "FCI::GSDavidson : Converged ground state energy = " << FCIenergy << endl; }
   
   for (int cnt=0; cnt<num_allocated; cnt++){
      delete [] vecs[cnt];
      delete [] Hvecs[cnt];
   }
   delete [] vecs;
   delete [] Hvecs;
   delete [] t_vec;
   delete [] u_vec;
   delete [] work_vec;
   delete [] mxM;
   delete [] mxM_eigs;
   delete [] mxM_vecs;
   delete [] mxM_work;
   delete [] Diag;
   
   if (Reortho_Allocated){
      delete [] Reortho_Eigenvecs;
      delete [] Reortho_Overlap;
      delete [] Reortho_Overlap_eigs;
      delete [] Reortho_Lowdin;
   }
   
   return FCIenergy;

}

void CheMPS2::FCI::CGSolveSystem(const double alpha, const double beta, const double eta, double * RHS, double * RealSol, double * ImagSol, const bool checkError) const{

   const unsigned long long vecLength = getVecLength();
   
   // Create a few helper arrays
   double * RESID  = new double[ vecLength ];
   double * PVEC   = new double[ vecLength ];
   double * OxPVEC = new double[ vecLength ];
   double * temp   = new double[ vecLength ];
   double * temp2  = new double[ vecLength ];
   double * precon = new double[ vecLength ];
   CGDiagPrecond( alpha , beta , eta , precon , temp );

   /* 
         ( alpha + beta H + I eta ) Solution = RHS
      
      is solved with the conjugate gradient (CG) method. Solution = RealSol + I * ImagSol.
      CG requires a symmetric positive definite operator. Therefore:
      
         precon * [ ( alpha + beta H )^2 + eta^2 ] * precon * SolutionTilde = precon * ( alpha + beta H - I eta ) * RHS
         Solution = precon * SolutionTilde
   */
   
   if ( RealSol != NULL ){
   
      CGAlphaPlusBetaHAM( alpha , beta , RHS , RESID ); // RESID = ( alpha + beta * H ) * RHS
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RESID[ cnt ] = precon[ cnt ] * RESID[ cnt ]; } // RESID = precon * ( alpha + beta * H ) * RHS
      if ( FCIverbose > 1 ){ cout << "FCI::CGSolveSystem : Two-norm of the RHS for the real part = " << FCIfrobeniusnorm( vecLength , RESID ) << endl; }
      //FCIdcopy( vecLength , RESID , RealSol ); // Initial guess for the real part
      FillRandom( vecLength, RealSol );
      CGCoreSolver( alpha , beta , eta , precon , RealSol , RESID , PVEC , OxPVEC , temp , temp2 ); // RESID contains the RHS of ( precon * Op * precon ) * |x> = |b>
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RealSol[ cnt ] = precon[ cnt ] * RealSol[ cnt ]; }
   
   }
   
   if ( ImagSol != NULL ){
   
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ RESID[ cnt ] = - eta * precon[ cnt ] * RHS[ cnt ]; } // RESID = - eta * precon * RHS
      if ( FCIverbose > 1 ){ cout << "FCI::CGSolveSystem : Two-norm of the RHS for the imaginary part = " << FCIfrobeniusnorm( vecLength , RESID ) << endl; }
      //FCIdcopy( vecLength , RESID , ImagSol ); // Initial guess for the imaginary part
      FillRandom( vecLength, ImagSol );
      CGCoreSolver( alpha , beta , eta , precon , ImagSol , RESID , PVEC , OxPVEC , temp , temp2 ); // RESID contains the RHS of ( precon * Op * precon ) * |x> = |b>
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ ImagSol[ cnt ] = precon[ cnt ] * ImagSol[ cnt ]; }
      
   }
   
   if (checkError){
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){ precon[ cnt ] = 1.0; }
      double RMSerror = 0.0;
      if ( RealSol != NULL ){
         CGOperator( alpha , beta , eta , precon , RealSol , temp , temp2 , OxPVEC );
         CGAlphaPlusBetaHAM( alpha , beta , RHS , RESID );
         FCIdaxpy( vecLength , -1.0 , RESID , OxPVEC );
         RMSerror += FCIddot( vecLength , OxPVEC , OxPVEC );
      }
      if ( ImagSol != NULL ){
         CGOperator( alpha , beta , eta , precon , ImagSol , temp , temp2 , OxPVEC );
         FCIdaxpy( vecLength , eta , RHS , OxPVEC );
         RMSerror += FCIddot( vecLength , OxPVEC , OxPVEC );
      }
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

   const unsigned long long vecLength = getVecLength();
   const double CGRESIDUALTHRESHOLD = 10.0 * CheMPS2::HEFF_DAVIDSON_RTOL_BASE * sqrt( 1.0 * vecLength );
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
      for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
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
   const unsigned long long vecLength = getVecLength();
   const double prefactor = alpha + beta * getEconst(); // HamTimesVec does only the parts with second quantized operators
   for (unsigned long long cnt = 0; cnt < vecLength; cnt++){
      out[ cnt ] = prefactor * in[ cnt ] + beta * out[ cnt ]; // out = ( alpha + beta * H ) * in
   }

}

void CheMPS2::FCI::CGOperator(const double alpha, const double beta, const double eta, double * precon, double * in, double * temp, double * temp2, double * out) const{

   const unsigned long long vecLength = getVecLength();
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

   // With operator = [ ( alpha + beta * H )^2 + eta*eta ] the preconditioner becomes 1 / sqrt( diag ( operator ) ).
   
   DiagHam( precon );
   DiagHamSquared( workspace );
   
   const unsigned long long vecLength = getVecLength();
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
      cout << "CGDiagPrecond : Minimum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << 1.0/(maxval*maxval) << endl;
      cout << "CGDiagPrecond : Maximum value of diag[ ( alpha + beta * Ham )^2 + eta^2 ] = " << 1.0/(minval*minval) << endl;
   }

}


