/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

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

#include "DMRGSCFmatrix.h"

CheMPS2::DMRGSCFmatrix::DMRGSCFmatrix( const DMRGSCFindices * iHandler ){

   this->iHandler = iHandler;
   this->num_irreps = iHandler->getNirreps();

   entries = new double*[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = iHandler->getNORB( irrep );
      entries[ irrep ] = new double[ NORB * NORB ];
   }

}

CheMPS2::DMRGSCFmatrix::~DMRGSCFmatrix(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      delete [] entries[ irrep ];
   }
   delete [] entries;

}

void CheMPS2::DMRGSCFmatrix::set( const int irrep, const int p, const int q, const double val ){

   entries[ irrep ][ p + iHandler->getNORB( irrep ) * q ] = val;

}

void CheMPS2::DMRGSCFmatrix::clear(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int size = iHandler->getNORB( irrep ) * iHandler->getNORB( irrep );
      for ( int counter = 0; counter < size; counter++ ){
         entries[ irrep ][ counter ] = 0.0;
      }
   }

}

void CheMPS2::DMRGSCFmatrix::identity(){

   clear();
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = iHandler->getNORB( irrep );
      for ( int diag = 0; diag < NORB; diag++ ){
         entries[ irrep ][ ( NORB + 1 ) * diag ] = 1.0;
      }
   }

}

double CheMPS2::DMRGSCFmatrix::get( const int irrep, const int p, const int q ) const{

   return entries[ irrep ][ p + iHandler->getNORB( irrep ) * q ];

}

double * CheMPS2::DMRGSCFmatrix::getBlock( const int irrep ){

   return entries[ irrep ];

}

double CheMPS2::DMRGSCFmatrix::rms_deviation( const DMRGSCFmatrix * other ) const{

   // Should in principle check whether iHandler matches
   double rms_diff = 0.0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = iHandler->getNORB( irrep );
      for ( int row = 0; row < NORB; row++ ){
         for ( int col = 0; col < NORB; col++ ){
            const double diff = this->get( irrep, row, col ) - other->get( irrep, row, col );
            rms_diff += diff * diff;
         }
      }
   }
   rms_diff = sqrt( rms_diff );
   return rms_diff;

}

#ifdef CHEMPS2_MPI_COMPILATION
void CheMPS2::DMRGSCFmatrix::broadcast( const int ROOT ){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = iHandler->getNORB( irrep );
      const int size = NORB * NORB;
      if ( size > 0 ){
         MPIchemps2::broadcast_array_double( entries[ irrep ], size, ROOT );
      }
   }

}
#endif

