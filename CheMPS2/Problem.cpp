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
#include <iostream>
#include <math.h> // fabs

#include "Problem.h"
#include "Irreps.h"
#include "MPIchemps2.h"

using std::cout;
using std::endl;

CheMPS2::Problem::Problem(const Hamiltonian * Hamin, const int TwoSin, const int Nin, const int Irrepin){

   Ham = Hamin;
   L = Ham->getL();
   TwoS = TwoSin;
   N = Nin;
   Irrep = Irrepin;
   bReorder = false;
   
   checkConsistency();
   mx_elem = NULL;

}

CheMPS2::Problem::~Problem(){

   if (bReorder){
      delete [] f1;
      delete [] f2;
   }
   
   if ( mx_elem != NULL ){ delete [] mx_elem; }

}

void CheMPS2::Problem::SetupReorderD2h(){

   if (bReorder){
      delete [] f1;
      delete [] f2;
      bReorder = false;
   }

   if (gSy()==7){ //Only if D2h of course
   
      bReorder = true;
      f1 = new int[Ham->getL()];
      f2 = new int[Ham->getL()];
      
      int DMRGirrepOrder[8];
      DMRGirrepOrder[0] = 0; //Ag  sigma
      DMRGirrepOrder[1] = 5; //B1u sigma^*
      DMRGirrepOrder[2] = 7; //B3u pi_x
      DMRGirrepOrder[3] = 2; //B2g pi_x^*
      DMRGirrepOrder[4] = 6; //B2u pi_y
      DMRGirrepOrder[5] = 3; //B3g pi_y^*
      DMRGirrepOrder[6] = 1; //B1g
      DMRGirrepOrder[7] = 4; //Au
      
      int DMRGOrb = 0;
      for (int irrep=0; irrep<8; irrep++){
         for (int HamOrb=0; HamOrb<Ham->getL(); HamOrb++){
            if (Ham->getOrbitalIrrep(HamOrb)==DMRGirrepOrder[irrep]){
               f1[HamOrb] = DMRGOrb;
               f2[DMRGOrb] = HamOrb;
               DMRGOrb++;
            }
         }
      }
      assert( DMRGOrb==Ham->getL() );
      
   }

}

void CheMPS2::Problem::SetupReorderC2v(){

   if (bReorder){
      delete [] f1;
      delete [] f2;
      bReorder = false;
   }

   if (gSy()==5){ //Only if C2v of course
   
      bReorder = true;
      f1 = new int[Ham->getL()];
      f2 = new int[Ham->getL()];
      
      int DMRGirrepOrder[4];
      DMRGirrepOrder[0] = 0; //A1
      DMRGirrepOrder[1] = 2; //B1
      DMRGirrepOrder[2] = 3; //B2
      DMRGirrepOrder[3] = 1; //A2
      
      int DMRGOrb = 0;
      {
         const int irrep = 0; // Irrep = 0 = A1 : reverse the order of the orbitals
         for (int HamOrb=Ham->getL()-1; HamOrb>=0; HamOrb--){
            if (Ham->getOrbitalIrrep(HamOrb)==DMRGirrepOrder[irrep]){
               f1[HamOrb] = DMRGOrb;
               f2[DMRGOrb] = HamOrb;
               DMRGOrb++;
            }
         }
      }
      for (int irrep=1; irrep<4; irrep++){
         for (int HamOrb=0; HamOrb<Ham->getL(); HamOrb++){
            if (Ham->getOrbitalIrrep(HamOrb)==DMRGirrepOrder[irrep]){
               f1[HamOrb] = DMRGOrb;
               f2[DMRGOrb] = HamOrb;
               DMRGOrb++;
            }
         }
      }
      assert( DMRGOrb==Ham->getL() );
      
      /*
      cout << "The new orbital order in c2v:" << endl;
      for ( int DMRGOrb=0; DMRGOrb<Ham->getL()-1; DMRGOrb++ ){
          cout << f2[DMRGOrb] << ", ";
      }
      cout << f2[Ham->getL()-1] << endl;
      */
      
   }

}

void CheMPS2::Problem::setup_reorder_custom(int * dmrg2ham){

   if (bReorder){
      delete [] f1;
      delete [] f2;
      bReorder = false;
   }
   
   bReorder = true;
   f1 = new int[Ham->getL()]; // Is going to be the inverse of dmrg2ham
   f2 = new int[Ham->getL()]; // Is going to be dmrg2ham copied

   // Set f1 entries negative to check all elements set
   for ( int  ham_orb = 0;  ham_orb < Ham->getL();  ham_orb++ ){ f1[ ham_orb ] = -2; }
   for ( int dmrg_orb = 0; dmrg_orb < Ham->getL(); dmrg_orb++ ){
   
      assert( dmrg2ham[ dmrg_orb ] >= 0           );
      assert( dmrg2ham[ dmrg_orb ] <  Ham->getL() );
      f2[ dmrg_orb ] = dmrg2ham[ dmrg_orb ];
      f1[ dmrg2ham[ dmrg_orb ] ] = dmrg_orb;
   
   }
   // Check all elements f1 set
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){ assert( f1[ ham_orb ] >= 0 ); }

}

void CheMPS2::Problem::setup_reorder_dinfh(int * docc, const double sp_threshold){

   assert( gSy() == 7 ); // Only for d2h of course
   const int num_irreps = 8;
   const int irrep_ag   = 0;
   const int irrep_b1g  = 1;
   const int irrep_b2g  = 2;
   const int irrep_b3g  = 3;
   const int irrep_au   = 4;
   const int irrep_b1u  = 5;
   const int irrep_b2u  = 6;
   const int irrep_b3u  = 7;
   
   double * sp_energies = new double[ Ham->getL() ];
   int * dmrg2ham = new int[ Ham->getL() ];
   int * partners = new int[ Ham->getL() ];
   
   // Get the single particle energies
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){
      double value = Ham->getTmat( ham_orb, ham_orb );
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         int counter = 0;
         for ( int frozen_orb = 0; frozen_orb < Ham->getL(); frozen_orb++ ){
            if ( Ham->getOrbitalIrrep( frozen_orb ) == irrep ){
               if ( counter < docc[ irrep ] ){
                  value += 2 * Ham->getVmat( ham_orb, frozen_orb, ham_orb, frozen_orb )
                             - Ham->getVmat( ham_orb, ham_orb, frozen_orb, frozen_orb );
               }
               counter += 1;
            }
         }
      }
      sp_energies[ ham_orb ] = value;
   }

   // To check that they have been set, put the partners to a negative value.
   for ( int cnt = 0; cnt < Ham->getL(); cnt++ ){ partners[ cnt ] = -2; }
   int dmrg_orb = 0;
   
   // Copy the b1g ( Delta_g ) orbitals
   for ( int ham_orb = Ham->getL() - 1; ham_orb >= 0; ham_orb-- ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b1g ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         for ( int partner_orb = 0; partner_orb < Ham->getL(); partner_orb++ ){
            if ( Ham->getOrbitalIrrep( partner_orb ) == irrep_ag ){
               if ( fabs( sp_energies[ ham_orb ] - sp_energies[ partner_orb ] ) < sp_threshold ){
                  partners[ dmrg_orb ] = partner_orb;
               }
            }
         }
         dmrg_orb++;
      }
   }
   
   // Copy the au ( Delta_u ) orbitals
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_au ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         for ( int partner_orb = 0; partner_orb < Ham->getL(); partner_orb++ ){
            if ( Ham->getOrbitalIrrep( partner_orb ) == irrep_b1u ){
               if ( fabs( sp_energies[ ham_orb ] - sp_energies[ partner_orb ] ) < sp_threshold ){
                  partners[ dmrg_orb ] = partner_orb;
               }
            }
         }
         dmrg_orb++;
      }
   }
   const int num_delta_ug = dmrg_orb;
   assert( (num_delta_ug % 2) == 0 );
   
   // Copy the partner orbitals
   for ( int cnt = 0; cnt < num_delta_ug; cnt++ ){
      assert( partners[ cnt ] >= 0 ); // Check that each one found a partner
      dmrg2ham[ dmrg_orb ] = partners[ cnt ];
      dmrg_orb++;
   }
   
   // Copy the remaining ag orbitals
   for ( int ham_orb = Ham->getL() - 1; ham_orb >= 0; ham_orb-- ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_ag ){
         bool is_a_partner = false;
         for ( int cnt = 0; cnt < num_delta_ug; cnt++ ){
            if ( ham_orb == partners[ cnt ] ){ is_a_partner = true; }
         }
         if ( is_a_partner == false ){
            dmrg2ham[ dmrg_orb ] = ham_orb;
            dmrg_orb++;
         }
      }
   }
   
   // Copy the remaining b1u orbitals
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b1u ){
         bool is_a_partner = false;
         for ( int cnt = 0; cnt < num_delta_ug; cnt++ ){
            if ( ham_orb == partners[ cnt ] ){ is_a_partner = true; }
         }
         if ( is_a_partner == false ){
            dmrg2ham[ dmrg_orb ] = ham_orb;
            dmrg_orb++;
         }
      }
   }
   
   // Copy the b3u ( Pi_u, pi_x ) orbitals
   for ( int ham_orb = Ham->getL() - 1; ham_orb >= 0; ham_orb-- ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b3u ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         dmrg_orb++;
      }
   }
   
   // Copy the b2g ( Pi_g, pi_x^* ) orbitals
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b2g ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         dmrg_orb++;
      }
   }
   
   // Copy the b2u ( Pi_u, pi_y ) orbitals
   for ( int ham_orb = Ham->getL() - 1; ham_orb >= 0; ham_orb-- ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b2u ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         dmrg_orb++;
      }
   }
   
   // Copy the b3g ( Pi_g, pi_y^* ) orbitals
   for ( int ham_orb = 0; ham_orb < Ham->getL(); ham_orb++ ){
      if ( Ham->getOrbitalIrrep( ham_orb ) == irrep_b3g ){
         dmrg2ham[ dmrg_orb ] = ham_orb;
         dmrg_orb++;
      }
   }
   
   assert( dmrg_orb == Ham->getL() );
   setup_reorder_custom( dmrg2ham );
   
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER )
   #endif
   {
      cout << "Reordered the orbitals according to d(infinity)h:" << endl;
      Irreps myIrreps( gSy() );
      for ( dmrg_orb = 0; dmrg_orb < Ham->getL(); dmrg_orb++ ){
         const int ham_orb = dmrg2ham[ dmrg_orb ];
         cout << "   DMRG orb " << dmrg_orb << " [" << myIrreps.getIrrepName(Ham->getOrbitalIrrep( ham_orb )) << "] has SP energy = " << sp_energies[ ham_orb ] << endl;
      }
   }
   
   delete [] dmrg2ham;
   delete [] partners;
   delete [] sp_energies;

}

int CheMPS2::Problem::gIrrep(const int nOrb) const{
   
   if (!bReorder){
      return Ham->getOrbitalIrrep(nOrb);
   }
   
   return Ham->getOrbitalIrrep(f2[nOrb]);

}

bool CheMPS2::Problem::gReorder() const{ return bReorder; }
int CheMPS2::Problem::gf1(const int HamOrb) const{ return (bReorder)?f1[HamOrb]:-1; }
int CheMPS2::Problem::gf2(const int DMRGOrb) const{ return (bReorder)?f2[DMRGOrb]:-1; }

double CheMPS2::Problem::gMxElement(const int alpha, const int beta, const int gamma, const int delta) const{

   return mx_elem[ alpha + L * ( beta + L * ( gamma + L * delta ) ) ];

}

void CheMPS2::Problem::setMxElement(const int alpha, const int beta, const int gamma, const int delta, const double value){

   mx_elem[ alpha + L * ( beta + L * ( gamma + L * delta ) ) ] = value;

}

void CheMPS2::Problem::construct_mxelem(){

   if ( mx_elem == NULL ){ mx_elem = new double[ L*L*L*L ]; }
   const double prefact = 1.0/(N-1);
   
   for (int orb1 = 0; orb1 < L; orb1++){
      const int map1 = (( !bReorder ) ? orb1 : f2[ orb1 ]);
      for (int orb2 = 0; orb2 < L; orb2++){
         const int map2 = (( !bReorder ) ? orb2 : f2[ orb2 ]);
         for (int orb3 = 0; orb3 < L; orb3++){
            const int map3 = (( !bReorder ) ? orb3 : f2[ orb3 ]);
            for (int orb4 = 0; orb4 < L; orb4++){
               const int map4 = (( !bReorder ) ? orb4 : f2[ orb4 ]);
               setMxElement( orb1, orb2, orb3, orb4, Ham->getVmat(map1,map2,map3,map4)
                                                     + prefact*((orb1==orb3)?Ham->getTmat(map2,map4):0)
                                                     + prefact*((orb2==orb4)?Ham->getTmat(map1,map3):0) );
            }
         }
      }
   }

}

bool CheMPS2::Problem::checkConsistency() const{

   Irreps SymmInfo(gSy());
   if ((gIrrep()<0) || (gIrrep()>=SymmInfo.getNumberOfIrreps())){
      cout << "Problem::Problem() : Irrep out of bound : Irrep = " << gIrrep() << endl;
      return false;
   }
   if (gTwoS()<0){
      cout << "Problem::checkConsistency() : TwoS = " << gTwoS() << endl;
      return false;
   }
   if (gN()<0){
      cout << "Problem::checkConsistency() : N = " << gN() << endl;
      return false;
   }
   if (gL()<0){
      cout << "Problem::checkConsistency() : L = " << gL() << endl;
      return false;
   }
   if (gN()>2*gL()){
      cout << "Problem::checkConsistency() : N > 2*L ; N = " << gN() << " and L = " << gL() << endl;
      return false;
   }
   if ( (gN()%2) != (gTwoS()%2) ){
      cout << "Problem::checkConsistency() : N%2 != TwoS%2 ; N = " << gN() << " and TwoS = " << gTwoS() << endl;
      return false;
   }
   if ( gTwoS() > gL() - abs(gN() - gL()) ){
      cout << "Problem::checkConsistency() : TwoS > L - |N-L| ; N = " << gN() << " and TwoS = " << gTwoS() << " and L = " << gL() << endl;
      return false;
   }
   
   return true;

}

bool CheMPS2::Problem::check_rohf_occ( int * occupancies ){

   int sum_n__tot = 0;
   int sum_2s_tot = 0;
   int sum_i__tot = 0;
   //Irreps SymmInfo(gSy());

   for ( int site = 0; site < gL(); site++ ){
      //cout << "Site " << site << " has irrep psi4 = " << gIrrep( site ) << " = " << SymmInfo.getIrrepName( gIrrep( site ) ) << endl;
      //cout << "Site " << site << " has occupancy  = " << occupancies[ site ] << endl;
      if (( occupancies[ site ] < 0 ) || ( occupancies[ site ] > 2 )){
         cout << "Problem::check_rohf_occ() : occupancies[ " << site << " ] = " << occupancies[ site ] << " and should be 0, 1 or 2." << endl;
         return false;
      }
      sum_n__tot += occupancies[ site ];
      if ( occupancies[ site ] == 1 ){
         sum_2s_tot += 1;
         sum_i__tot = Irreps::directProd( sum_i__tot, gIrrep( site ) );
      }
   }

   if (( sum_n__tot != gN() ) || ( sum_2s_tot != gTwoS() ) || ( sum_i__tot != gIrrep() )){
      cout << "Problem::check_rohf_occ() : occupancies corresponds to ( N, 2S, I ) = ( " << sum_n__tot << ", " << sum_2s_tot << ", " << sum_i__tot << " ), while the DMRG targeted sector is ( N, 2S, I ) = ( " << gN() << ", " << gTwoS() << ", " << gIrrep() << " )." << endl;
      return false;
   }

   return true;

}


