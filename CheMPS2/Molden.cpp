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
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>

#include "Molden.h"
#include "Lapack.h"
#include "DMRGSCFmatrix.h"

using std::ifstream;
using std::ofstream;

CheMPS2::Molden::Molden( const int L, const int group, int * irrep_sizes ){

   this->L = L;
   SymmInfo.setGroup( group );
   this->num_irreps = SymmInfo.getNumberOfIrreps();

   Isizes  = new int[ num_irreps ];
   molden  = new double*[ num_irreps ];
   unitary = new double*[ num_irreps ];
   product = new double*[ num_irreps ];

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      Isizes [ irrep ] = irrep_sizes[ irrep ];
      molden [ irrep ] = new double[ Isizes[ irrep ] * L ];
      unitary[ irrep ] = new double[ Isizes[ irrep ] * Isizes[ irrep ] ];
      product[ irrep ] = new double[ Isizes[ irrep ] * L ];
   }

}

CheMPS2::Molden::~Molden(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      delete [] molden[ irrep ];
      delete [] unitary[ irrep ];
      delete [] product[ irrep ];
   }

   delete [] molden;
   delete [] unitary;
   delete [] product;
   delete [] Isizes;

}

void CheMPS2::Molden::read_molden( const string filename ){

   int * Icounter = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ Icounter[ irrep ] = 0; }

   int * psi2molpro = new int[ num_irreps ];
   CheMPS2::Irreps::symm_psi2molpro( psi2molpro, SymmInfo.getGroupName() );

   string line, part;

   std::ifstream inputfile( filename.c_str() );

   bool stop = false;
   string start = "[MO]";
   do {
      getline( inputfile, line );
      const int pos = line.find( start );
      if ( pos != string::npos ){ stop = true; }
   } while ( stop == false );

   int total_num_orbs = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ total_num_orbs += Isizes[ irrep ]; }

   for ( int orbs = 0; orbs < total_num_orbs; orbs++ ){ // Only read in the alpha orbitals

      // Read in the Symm
      getline( inputfile, line );
      const int pos_equal = line.find( '=' );
      const int pos_dot   = line.find( '.', pos_equal + 1 );
      int psi4_irrep      = -1;
      if ( pos_dot != string::npos ){ // Molpro
         part = line.substr( pos_dot + 1, line.size() - 1 - pos_dot );
         const int molpro_irrep = atoi( part.c_str() );
         for ( int irrep = 0; irrep < num_irreps; irrep++ ){
            if ( psi2molpro[ irrep ] == molpro_irrep ){ psi4_irrep = irrep; }
         }
         //std::cout << "MO " << orbs + 1 << " has molpro_irrep " << molpro_irrep << " and psi4_irrep " << SymmInfo.getIrrepName( psi4_irrep ) << "." << std::endl;
      } else { // Psi4
         part = line.substr( pos_equal + 1, line.size() - 1 - pos_equal );
         for ( int irrep = 0; irrep < num_irreps; irrep++ ){
            if ( part.find( SymmInfo.getIrrepName( irrep ) ) != string::npos ){ psi4_irrep = irrep; }
         }
         if ( SymmInfo.getGroupNumber() == 3 ){ // Cs
            if ( part.find( "A'"  ) != string::npos ){ psi4_irrep = 0; }
            if ( part.find( "A\"" ) != string::npos ){ psi4_irrep = 1; }
         }
         //std::cout << "MO " << orbs + 1 << " has psi4_irrep " << SymmInfo.getIrrepName( psi4_irrep ) << "." << std::endl;
      }
      assert( psi4_irrep != -1 );

      // Read in the Ene, Spin, Occup lines
      for ( int cnt = 0; cnt < 3; cnt++ ){ getline( inputfile, line ); }

      for ( int cnt = 0; cnt < L; cnt++ ){
         getline( inputfile, line );
         const int pos_number = line.find_first_not_of( ' ' ); // Start of the primitive counter number
         const int pos_gap    = line.find( ' ', pos_number );  // Gap between the primitive number and the actual coefficient
         part = line.substr( pos_gap + 1, line.size() - 1 - pos_gap );
         const double value = atof( part.c_str() );
         molden[ psi4_irrep ][ cnt + L * Icounter[ psi4_irrep ] ] = value;
         //std::cout << cnt+1 << " " << value << std::endl;
      }

      Icounter[ psi4_irrep ] += 1;

   }

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      assert( Icounter[ irrep ] == Isizes[ irrep ] );
   }

   delete [] Icounter;
   delete [] psi2molpro;

   inputfile.close();

}

void CheMPS2::Molden::read_unitary( const string filename ){

   CheMPS2::DMRGSCFmatrix::read( filename, num_irreps, unitary );

}

void CheMPS2::Molden::print( const string original, const string output ){

   char trans   = 'T';
   char notrans = 'N';
   double one   = 1.0;
   double set   = 0.0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      if ( Isizes[ irrep ] > 0 ){
         dgemm_( &notrans, &trans, &L, Isizes + irrep, Isizes + irrep, &one, molden[ irrep ], &L, unitary[ irrep ], Isizes + irrep, &set, product[ irrep ], &L );
      }
   }

   string line, part;

   std::ifstream inputfile( original.c_str() );
   std::ofstream outputfile( output.c_str(), std::ofstream::trunc );

   // First go to the start of the orbital coefficient dump
   bool stop = false;
   string start = "[MO]";
   do {
      getline( inputfile,line );
      outputfile << line << std::endl;
      const int pos = line.find( start );
      if ( pos != string::npos ){ stop = true; }
   } while ( stop == false );

   inputfile.close();

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      for ( int orb = 0; orb < Isizes[ irrep ]; orb++ ){
         outputfile << " Sym= " << SymmInfo.getIrrepName( irrep ) << std::endl;
         outputfile << " Ene= N/A" << std::endl;
         outputfile << " Spin= Restricted" << std::endl;
         outputfile << " Occup= N/A" << std::endl;
         for ( int cnt = 0; cnt < L; cnt++ ){
            outputfile << cnt + 1 << " " << product[ irrep ][ cnt + L * orb ] << std::endl;
         }
      }
   }

   outputfile.close();

}


