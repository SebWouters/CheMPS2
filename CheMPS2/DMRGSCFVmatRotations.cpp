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

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <string>

#include "DMRGSCFVmatRotations.h"
#include "Lapack.h"

using std::min;
using std::max;
using std::string;

CheMPS2::DMRGSCFVmatRotations::DMRGSCFVmatRotations(){ }
CheMPS2::DMRGSCFVmatRotations::~DMRGSCFVmatRotations(){ }

void CheMPS2::DMRGSCFVmatRotations::fetch( double * eri, const FourIndex * ORIG_VMAT, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack ){

   if ( pack ){

      assert( irrep3 == irrep4 );
      assert( irrep1 == irrep2 );

      const int NORB12 = idx->getNORB( irrep1 );
      const int NORB34 = idx->getNORB( irrep3 );

      int counter = 0; // counter = cnt3 + ( cnt4 * ( cnt4 + 1 )) / 2
      for ( int cnt4 = 0; cnt4 < NORB34; cnt4++ ){
         for ( int cnt3 = 0; cnt3 <= cnt4; cnt3++ ){
            if (( start <= counter ) && ( counter < stop )){
               for ( int cnt2 = 0; cnt2 < NORB12; cnt2++ ){
                  for ( int cnt1 = 0; cnt1 < NORB12; cnt1++ ){
                     eri[ cnt1 + NORB12 * ( cnt2 + NORB12 * ( counter - start ) ) ]
                        = ORIG_VMAT->get( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, cnt2, cnt4 );
                        // Indices (12) and indices (34) are Coulomb pairs
                  }
               }
            }
            counter++;
         }
      }

   } else {

      assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

      const int NORB1 = idx->getNORB( irrep1 );
      const int NORB2 = idx->getNORB( irrep2 );
      const int NORB3 = idx->getNORB( irrep3 );
      const int NORB4 = idx->getNORB( irrep4 );

      int counter = 0; // counter = cnt3 + NORB3 * cnt4
      for ( int cnt4 = 0; cnt4 < NORB4; cnt4++ ){
         for ( int cnt3 = 0; cnt3 < NORB3; cnt3++ ){
            if (( start <= counter ) && ( counter < stop )){
               for ( int cnt2 = 0; cnt2 < NORB2; cnt2++ ){
                  for ( int cnt1 = 0; cnt1 < NORB1; cnt1++ ){
                     eri[ cnt1 + NORB1 * ( cnt2 + NORB2 * ( counter - start ) ) ]
                        = ORIG_VMAT->get( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, cnt2, cnt4 );
                        // Indices (12) and indices (34) are Coulomb pairs
                  }
               }
            }
            counter++;
         }
      }

   }

}

void CheMPS2::DMRGSCFVmatRotations::write( double * eri, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack ){

   assert(( space == 'F' ) || ( space == 'A' ) || ( space == 'C' ) || ( space == 'E' ));

   if (( space == 'A' ) || ( space =='F' )){

      if ( pack ){

         assert( irrep1 == irrep2 );
         assert( irrep3 == irrep4 );

         const int NEW12 = (( space == 'A' ) ? idx->getNDMRG( irrep1 ) : idx->getNORB( irrep1 ));
         const int NEW34 = (( space == 'A' ) ? idx->getNDMRG( irrep3 ) : idx->getNORB( irrep3 ));
         const int SIZE  = stop - start;

         int counter = 0; // counter = cnt1 + ( cnt2 * ( cnt2 + 1 )) / 2
         for ( int cnt2 = 0; cnt2 < NEW12; cnt2++ ){
            for ( int cnt1 = 0; cnt1 <= cnt2; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW34; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 <= cnt4; cnt3++ ){
                        NEW_VMAT->set( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, cnt2, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW34 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }

      } else {

         assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

         const int NEW1 = (( space == 'A' ) ? idx->getNDMRG( irrep1 ) : idx->getNORB( irrep1 ));
         const int NEW2 = (( space == 'A' ) ? idx->getNDMRG( irrep2 ) : idx->getNORB( irrep2 ));
         const int NEW3 = (( space == 'A' ) ? idx->getNDMRG( irrep3 ) : idx->getNORB( irrep3 ));
         const int NEW4 = (( space == 'A' ) ? idx->getNDMRG( irrep4 ) : idx->getNORB( irrep4 ));
         const int SIZE = stop - start;

         int counter = 0; // counter = cnt1 + NEW1 * cnt2
         for ( int cnt2 = 0; cnt2 < NEW2; cnt2++ ){
            for ( int cnt1 = 0; cnt1 < NEW1; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW4; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 < NEW3; cnt3++ ){
                        NEW_VMAT->set( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, cnt2, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW3 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }
      }
   }

   if ( space == 'C' ){

      if ( pack ){

         assert( irrep1 == irrep2 );
         assert( irrep3 == irrep4 );

         const int NEW_C12 = idx->getNOCC( irrep1 ) + idx->getNDMRG( irrep1 );
         const int NEW_A34 = idx->getNORB( irrep3 );
         const int SIZE    = stop - start;

         int counter = 0; // counter = cnt1 + ( cnt2 * ( cnt2 + 1 )) / 2
         for ( int cnt2 = 0; cnt2 < NEW_C12; cnt2++ ){
            for ( int cnt1 = 0; cnt1 <= cnt2; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW_A34; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 <= cnt4; cnt3++ ){
                        ROT_TEI->set_coulomb( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW_A34 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }

      } else {

         assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

         const int NEW_C1 = idx->getNOCC( irrep1 ) + idx->getNDMRG( irrep1 );
         const int NEW_C2 = idx->getNOCC( irrep2 ) + idx->getNDMRG( irrep2 );
         const int NEW_A3 = idx->getNORB( irrep3 );
         const int NEW_A4 = idx->getNORB( irrep4 );
         const int SIZE   = stop - start;

         int counter = 0; // counter = cnt1 + NEW_C1 * cnt2
         for ( int cnt2 = 0; cnt2 < NEW_C2; cnt2++ ){
            for ( int cnt1 = 0; cnt1 < NEW_C1; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW_A4; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 < NEW_A3; cnt3++ ){
                        ROT_TEI->set_coulomb( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW_A3 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }
      }
   }

   if ( space == 'E' ){

      assert( pack == false );
      assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

      const int  NEW_C1 = idx->getNOCC( irrep1 ) + idx->getNDMRG( irrep1 );
      const int  NEW_C3 = idx->getNOCC( irrep3 ) + idx->getNDMRG( irrep3 );
      const int  NEW_V2 = idx->getNVIRT( irrep2 );
      const int  NEW_V4 = idx->getNVIRT( irrep4 );
      const int JUMP_V2 = idx->getNOCC( irrep2 ) + idx->getNDMRG( irrep2 );
      const int JUMP_V4 = idx->getNOCC( irrep4 ) + idx->getNDMRG( irrep4 );
      const int SIZE    = stop - start;

      int counter = 0; // counter = cnt1 + NEW_C1 * cnt2
      for ( int cnt2 = 0; cnt2 < NEW_V2; cnt2++ ){
         for ( int cnt1 = 0; cnt1 < NEW_C1; cnt1++ ){
            if (( start <= counter ) && ( counter < stop )){
               for ( int cnt4 = 0; cnt4 < NEW_V4; cnt4++ ){
                  for ( int cnt3 = 0; cnt3 < NEW_C3; cnt3++ ){
                     ROT_TEI->set_exchange( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, JUMP_V2 + cnt2, JUMP_V4 + cnt4,
                        eri[ ( counter - start ) + SIZE * ( cnt3 + NEW_C3 * cnt4 ) ] );
                        // Indices (12) and indices (34) are Coulomb pairs
                  }
               }
            }
            counter++;
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_first( double * origin, double * target, int orig1, int dim2, const int dim34, double * umat1, int new1, int lda1 ){

   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int right_dim = dim2 * dim34;
   dgemm_( &notrans, &notrans, &new1, &right_dim, &orig1, &one, umat1, &lda1, origin, &orig1, &set, target, &new1 );

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_second( double * origin, double * target, int dim1, int orig2, const int dim34, double * umat2, int new2, int lda2 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   const int right_dim = dim34;
   const int jump_old  = dim1 * orig2;
   const int jump_new  = dim1 * new2;
   for ( int index = 0; index < right_dim; index++ ){
      dgemm_( &notrans, &trans, &dim1, &new2, &orig2, &one, origin + jump_old * index, &dim1, umat2, &lda2, &set, target + jump_new * index, &dim1 );
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_third( double * origin, double * target, const int dim12, int orig3, int dim4, double * umat3, int new3, int lda3 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim12;
   const int jump_old = dim12 * orig3;
   const int jump_new = dim12 * new3;
   for ( int index = 0; index < dim4; index++ ){
      dgemm_( &notrans, &trans, &left_dim, &new3, &orig3, &one, origin + jump_old * index, &left_dim, umat3, &lda3, &set, target + jump_new * index, &left_dim );
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_fourth( double * origin, double * target, const int dim12, int dim3, int orig4, double * umat4, int new4, int lda4 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim12 * dim3;
   dgemm_( &notrans, &trans, &left_dim, &new4, &orig4, &one, origin, &left_dim, umat4, &lda4, &set, target, &left_dim );

}

void CheMPS2::DMRGSCFVmatRotations::open_file( hid_t * file_id, hid_t * dspc_id, hid_t * dset_id, const int first, const int second, const string filename ){

   file_id[0] = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   hsize_t fdim_h5[] = { second, first }; // C is row major: [ col + ncol * row ] is assumed
   dspc_id[0] = H5Screate_simple( 2, fdim_h5, NULL );
   dset_id[0] = H5Dcreate( file_id[0], "storage", H5T_NATIVE_DOUBLE, dspc_id[0], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

}

void CheMPS2::DMRGSCFVmatRotations::close_file( hid_t file_id, hid_t dspc_id, hid_t dset_id ){

   H5Dclose( dset_id );
   H5Sclose( dspc_id );
   H5Fclose( file_id );

}

void CheMPS2::DMRGSCFVmatRotations::write_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int first_write ){

   hsize_t stride_h5[] = { 1, 1 };
   hsize_t  count_h5[] = { 1, 1 };
   hsize_t  start_h5[] = { start, 0 };
   hsize_t  block_h5[] = { size, first_write };
   H5Sselect_hyperslab( dspc_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
   hsize_t mem_h5 = size * first_write; // Should be OK to multiply as integers as it is smaller than mem_size
   hid_t mem_id = H5Screate_simple( 1, &mem_h5, NULL );
   H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, mem_id, dspc_id, H5P_DEFAULT, eri );
   H5Sclose( mem_id );

}

void CheMPS2::DMRGSCFVmatRotations::read_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int second_read ){

   hsize_t stride_h5[] = { 1, 1 };
   hsize_t  count_h5[] = { 1, 1 };
   hsize_t  start_h5[] = { 0, start };
   hsize_t  block_h5[] = { second_read, size };
   H5Sselect_hyperslab( dspc_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
   hsize_t mem_h5 = second_read * size; // Should be OK to multiply as integers as it is smaller than mem_size
   hid_t mem_id = H5Screate_simple( 1, &mem_h5, NULL );
   H5Dread( dset_id, H5T_NATIVE_DOUBLE, mem_id, dspc_id, H5P_DEFAULT, eri );
   H5Sclose( mem_id );

}

void CheMPS2::DMRGSCFVmatRotations::rotate( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, const char space, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename ){

   assert(( space == 'A' ) || ( space == 'F' ));
   const int num_irreps = idx->getNirreps();

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      for ( int irrep2 = irrep1; irrep2 < num_irreps; irrep2++ ){ // irrep2 >= irrep1
         const int product_symm = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = irrep1; irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( product_symm, irrep3 );
            if ( irrep4 >= irrep3 ){ // irrep4 >= irrep3

               const int NEW1 = (( space == 'A' ) ? idx->getNDMRG( irrep1 ) : idx->getNORB( irrep1 ));
               const int NEW2 = (( space == 'A' ) ? idx->getNDMRG( irrep2 ) : idx->getNORB( irrep2 ));
               const int NEW3 = (( space == 'A' ) ? idx->getNDMRG( irrep3 ) : idx->getNORB( irrep3 ));
               const int NEW4 = (( space == 'A' ) ? idx->getNDMRG( irrep4 ) : idx->getNORB( irrep4 ));

               if (( NEW1 > 0 ) && ( NEW2 > 0 ) && ( NEW3 > 0 ) && ( NEW4 > 0 )){

                  const int ORIG1 = idx->getNORB( irrep1 );
                  const int ORIG2 = idx->getNORB( irrep2 );
                  const int ORIG3 = idx->getNORB( irrep3 );
                  const int ORIG4 = idx->getNORB( irrep4 );

                  double * umat1 = umat->getBlock( irrep1 ) + (( space == 'A' ) ? idx->getNOCC( irrep1 ) : 0 );
                  double * umat2 = umat->getBlock( irrep2 ) + (( space == 'A' ) ? idx->getNOCC( irrep2 ) : 0 );
                  double * umat3 = umat->getBlock( irrep3 ) + (( space == 'A' ) ? idx->getNOCC( irrep3 ) : 0 );
                  double * umat4 = umat->getBlock( irrep4 ) + (( space == 'A' ) ? idx->getNOCC( irrep4 ) : 0 );

                  const int  first_new =  NEW1 * NEW2;
                  const int second_old = ORIG3 * ORIG4;

                  const int block_size1 = mem_size / ( ORIG1 * ORIG2 ); // Floor of amount of times first_old  fits in mem_size
                  const int block_size2 = mem_size / second_old;        // Floor of amount of times second_old fits in mem_size
                  assert( block_size1 > 0 );
                  assert( block_size2 > 0 );

                  const bool io_free = (( block_size1 >= second_old ) && ( block_size2 >= first_new ));
                  hid_t file_id, dspc_id, dset_id;

                  if ( io_free == false ){
                     assert( filename.compare( "edmistonruedenberg" ) != 0 );
                     open_file( &file_id, &dspc_id, &dset_id, first_new, second_old, filename );
                  }

                  // First half transformation
                  int start = 0;
                  while ( start < second_old ){
                     const int stop = min( start + block_size1, second_old );
                     const int size = stop - start;
                     fetch( mem1, ORIG_VMAT, irrep1, irrep2, irrep3, irrep4, idx, start, stop, false );
                     blockwise_first(  mem1, mem2, ORIG1, ORIG2, size, umat1, NEW1, ORIG1 );
                     blockwise_second( mem2, mem1, NEW1,  ORIG2, size, umat2, NEW2, ORIG2 );
                     // pack first potentially
                     if ( io_free == false ){ write_file( dspc_id, dset_id, mem1, start, size, first_new ); }
                     start += size;
                  }
                  assert( start == second_old );

                  // Do the second half transformation
                  start = 0;
                  while ( start < first_new ){
                     const int stop = min( start + block_size2, first_new );
                     const int size = stop - start;
                     if ( io_free == false ){ read_file( dspc_id, dset_id, mem1, start, size, second_old ); }
                     // unpack second potentially
                     blockwise_fourth( mem1, mem2, size, ORIG3, ORIG4, umat4, NEW4, ORIG4 );
                     blockwise_third(  mem2, mem1, size, ORIG3, NEW4,  umat3, NEW3, ORIG3 );
                     write( mem1, NEW_VMAT, NULL, space, irrep1, irrep2, irrep3, irrep4, idx, start, stop, false );
                     start += size;
                  }
                  assert( start == first_new );
                  if ( io_free == false ){ close_file( file_id, dspc_id, dset_id ); }
               }
            }
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::rotate( const FourIndex * ORIG_VMAT, DMRGSCFintegrals * ROT_TEI, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename ){

   const int num_irreps = idx->getNirreps();

   // First do Coulomb object : ( c1 <= c2 | a1 <= a2 )
   for ( int Ic1 = 0; Ic1 < num_irreps; Ic1++){
      for ( int Ic2 = Ic1; Ic2 < num_irreps; Ic2++){ // Ic2 >= Ic1
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for ( int Ia1 = 0; Ia1 < num_irreps; Ia1++ ){
            const int Ia2 = Irreps::directProd( Ia1, Icc );
            if ( Ia1 <= Ia2 ){ // Ia2 >= Ia1

               const int NEW_C1 = idx->getNOCC( Ic1 ) + idx->getNDMRG( Ic1 );
               const int NEW_C2 = idx->getNOCC( Ic2 ) + idx->getNDMRG( Ic2 );
               const int NEW_A1 = idx->getNORB( Ia1 );
               const int NEW_A2 = idx->getNORB( Ia2 );

               if (( NEW_C1 > 0 ) && ( NEW_C2 > 0 ) && ( NEW_A1 > 0 ) && ( NEW_A2 > 0 )){

                  const int NORB_C1 = idx->getNORB( Ic1 );
                  const int NORB_C2 = idx->getNORB( Ic2 );
                  const int NORB_A1 = idx->getNORB( Ia1 );
                  const int NORB_A2 = idx->getNORB( Ia2 );

                  double * umat_C1 = umat->getBlock( Ic1 );
                  double * umat_C2 = umat->getBlock( Ic2 );
                  double * umat_A1 = umat->getBlock( Ia1 );
                  double * umat_A2 = umat->getBlock( Ia2 );

                  const int  first_new =  NEW_C1 * NEW_C2;
                  const int second_old = NORB_A1 * NORB_A2;

                  const int block_size1 = mem_size / ( NORB_C1 * NORB_C2 ); // Floor of amount of times first_old  fits in mem_size
                  const int block_size2 = mem_size / second_old;            // Floor of amount of times second_old fits in mem_size
                  assert( block_size1 > 0 );
                  assert( block_size2 > 0 );

                  const bool io_free = (( block_size1 >= second_old ) && ( block_size2 >= first_new ));
                  hid_t file_id, dspc_id, dset_id;

                  if ( io_free == false ){
                     open_file( &file_id, &dspc_id, &dset_id, first_new, second_old, filename );
                  }

                  // First half transformation
                  int start = 0;
                  while ( start < second_old ){
                     const int stop = min( start + block_size1, second_old );
                     const int size = stop - start;
                     fetch( mem1, ORIG_VMAT, Ic1, Ic2, Ia1, Ia2, idx, start, stop, false );
                     blockwise_first(  mem1, mem2, NORB_C1, NORB_C2, size, umat_C1, NEW_C1, NORB_C1 );
                     blockwise_second( mem2, mem1,  NEW_C1, NORB_C2, size, umat_C2, NEW_C2, NORB_C2 );
                     // pack first potentially
                     if ( io_free == false ){ write_file( dspc_id, dset_id, mem1, start, size, first_new ); }
                     start += size;
                  }
                  assert( start == second_old );

                  // Do the second half transformation
                  start = 0;
                  while ( start < first_new ){
                     const int stop = min( start + block_size2, first_new );
                     const int size = stop - start;
                     if ( io_free == false ){ read_file( dspc_id, dset_id, mem1, start, size, second_old ); }
                     // unpack second potentially
                     blockwise_fourth( mem1, mem2, size, NORB_A1, NORB_A2, umat_A2, NEW_A2, NORB_A2 );
                     blockwise_third(  mem2, mem1, size, NORB_A1,  NEW_A2, umat_A1, NEW_A1, NORB_A1 );
                     write( mem1, NULL, ROT_TEI, 'C', Ic1, Ic2, Ia1, Ia2, idx, start, stop, false );
                     start += size;
                  }
                  assert( start == first_new );
                  if ( io_free == false ){ close_file( file_id, dspc_id, dset_id ); }
               }
            }
         }
      }
   }

   // Now do Exchange object ( c1 v1 | c2 v2 ) with c1 <= c2
   for ( int Ic1 = 0; Ic1 < num_irreps; Ic1++ ){
      for ( int Ic2 = Ic1; Ic2 < num_irreps; Ic2++ ){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for ( int Iv1 = 0; Iv1 < num_irreps; Iv1++ ){
            const int Iv2 = Irreps::directProd( Iv1, Icc );

            const int NEW_C1 = idx->getNOCC( Ic1 ) + idx->getNDMRG( Ic1 );
            const int NEW_C2 = idx->getNOCC( Ic2 ) + idx->getNDMRG( Ic2 );
            const int NEW_V1 = idx->getNVIRT( Iv1 );
            const int NEW_V2 = idx->getNVIRT( Iv2 );

            if (( NEW_C1 > 0 ) && ( NEW_C2 > 0 ) && ( NEW_V1 > 0 ) && ( NEW_V2 > 0 )){

               const int NORB_C1 = idx->getNORB( Ic1 );
               const int NORB_C2 = idx->getNORB( Ic2 );
               const int NORB_V1 = idx->getNORB( Iv1 );
               const int NORB_V2 = idx->getNORB( Iv2 );

               const int JUMP_V1 = idx->getNOCC( Iv1 ) + idx->getNDMRG( Iv1 );
               const int JUMP_V2 = idx->getNOCC( Iv2 ) + idx->getNDMRG( Iv2 );

               double * umat_C1 = umat->getBlock( Ic1 );
               double * umat_C2 = umat->getBlock( Ic2 );
               double * umat_V1 = umat->getBlock( Iv1 ) + JUMP_V1;
               double * umat_V2 = umat->getBlock( Iv2 ) + JUMP_V2;

               const int  first_new =  NEW_C1 * NEW_V1;
               const int second_old = NORB_C2 * NORB_V2;

               const int block_size1 = mem_size / ( NORB_C1 * NORB_V1 ); // Floor of amount of times first_old  fits in mem_size
               const int block_size2 = mem_size / second_old;            // Floor of amount of times second_old fits in mem_size
               assert( block_size1 > 0 );
               assert( block_size2 > 0 );

               const bool io_free = (( block_size1 >= second_old ) && ( block_size2 >= first_new ));
               hid_t file_id, dspc_id, dset_id;

               if ( io_free == false ){
                  open_file( &file_id, &dspc_id, &dset_id, first_new, second_old, filename );
               }

               // First half transformation
               int start = 0;
               while ( start < second_old ){
                  const int stop = min( start + block_size1, second_old );
                  const int size = stop - start;
                  fetch( mem1, ORIG_VMAT, Ic1, Iv1, Ic2, Iv2, idx, start, stop, false );
                  blockwise_first(  mem1, mem2, NORB_C1, NORB_V1, size, umat_C1, NEW_C1, NORB_C1 );
                  blockwise_second( mem2, mem1,  NEW_C1, NORB_V1, size, umat_V1, NEW_V1, NORB_V1 );
                  // do not pack first because EXCHANGE
                  if ( io_free == false ){ write_file( dspc_id, dset_id, mem1, start, size, first_new ); }
                  start += size;
               }
               assert( start == second_old );

               // Do the second half transformation
               start = 0;
               while ( start < first_new ){
                  const int stop = min( start + block_size2, first_new );
                  const int size = stop - start;
                  if ( io_free == false ){ read_file( dspc_id, dset_id, mem1, start, size, second_old ); }
                  // unpack second potentially --> is allowed for EXCHANGE
                  blockwise_fourth( mem1, mem2, size, NORB_C2, NORB_V2, umat_V2, NEW_V2, NORB_V2 );
                  blockwise_third(  mem2, mem1, size, NORB_C2,  NEW_V2, umat_C2, NEW_C2, NORB_C2 );
                  write( mem1, NULL, ROT_TEI, 'E', Ic1, Iv1, Ic2, Iv2, idx, start, stop, false );
                  start += size;
               }
               assert( start == first_new );
               if ( io_free == false ){ close_file( file_id, dspc_id, dset_id ); }
            }
         }
      }
   }

}



