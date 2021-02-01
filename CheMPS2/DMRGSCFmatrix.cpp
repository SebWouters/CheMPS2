/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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
#include <string>
#include <sstream>

#include <hdf5.h>

#include "DMRGSCFmatrix.h"
#include "MPIchemps2.h"

using std::string;

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

void CheMPS2::DMRGSCFmatrix::write( const string filename, const DMRGSCFindices * idx, double ** storage ){

   hid_t file_id  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   hid_t group_id = H5Gcreate( file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

   for ( int irrep = 0; irrep < idx->getNirreps(); irrep++ ){

      std::stringstream irrepname;
      irrepname << "irrep_" << irrep;

      hsize_t dimarray   = idx->getNORB( irrep ) * idx->getNORB( irrep );
      hid_t dataspace_id = H5Screate_simple( 1, &dimarray, NULL );
      hid_t dataset_id   = H5Dcreate( group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
      H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[ irrep ] );

      H5Dclose( dataset_id );
      H5Sclose( dataspace_id );

   }

   H5Gclose( group_id );
   H5Fclose( file_id );

}

void CheMPS2::DMRGSCFmatrix::read( const string filename, const int n_irreps, double ** storage ){

   hid_t file_id  = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
   hid_t group_id = H5Gopen( file_id, "/Data", H5P_DEFAULT );

   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      std::stringstream irrepname;
      irrepname << "irrep_" << irrep;

      hid_t dataset_id = H5Dopen( group_id, irrepname.str().c_str(), H5P_DEFAULT );
      H5Dread( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[ irrep ] );

      H5Dclose( dataset_id );

   }

   H5Gclose( group_id );
   H5Fclose( file_id );

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

