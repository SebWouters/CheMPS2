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
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <string>

#include "DMRGSCFrotations.h"
#include "Lapack.h"

using std::min;
using std::max;
using std::string;

void CheMPS2::DMRGSCFrotations::fetch( double * eri, const FourIndex * ORIG_VMAT, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack ){

   if ( pack ){

      assert( irrep1 == irrep2 );
      assert( irrep3 == irrep4 );

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

void CheMPS2::DMRGSCFrotations::write( double * eri, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space1, const char space2, const char space3, const char space4, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack ){

   bool written = false;
   if (( space1 == space2 ) && ( space1 == space3 ) && ( space1 == space4 )){ // All four spaces equal

      if ( pack ){

         assert( irrep1 == irrep2 );
         assert( irrep3 == irrep4 );

         const int NEW12 = dimension( idx, irrep1, space1 );
         const int NEW34 = dimension( idx, irrep3, space3 );
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
         written = true;

      } else {

         assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

         const int NEW1 = dimension( idx, irrep1, space1 );
         const int NEW2 = dimension( idx, irrep2, space2 );
         const int NEW3 = dimension( idx, irrep3, space3 );
         const int NEW4 = dimension( idx, irrep4, space4 );
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
         written = true;

      }
   }

   if (( space1 == 'C' ) && ( space2 == 'C' ) && ( space3 == 'F' ) && ( space4 == 'F' )){

      if ( pack ){

         assert( irrep1 == irrep2 );
         assert( irrep3 == irrep4 );

         const int NEW12 = dimension( idx, irrep1, space1 );
         const int NEW34 = dimension( idx, irrep3, space3 );
         const int SIZE  = stop - start;

         int counter = 0; // counter = cnt1 + ( cnt2 * ( cnt2 + 1 )) / 2
         for ( int cnt2 = 0; cnt2 < NEW12; cnt2++ ){
            for ( int cnt1 = 0; cnt1 <= cnt2; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW34; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 <= cnt4; cnt3++ ){
                        ROT_TEI->set_coulomb( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW34 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }
         written = true;

      } else {

         assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

         const int NEW1 = dimension( idx, irrep1, space1 );
         const int NEW2 = dimension( idx, irrep2, space2 );
         const int NEW3 = dimension( idx, irrep3, space3 );
         const int NEW4 = dimension( idx, irrep4, space4 );
         const int SIZE = stop - start;

         int counter = 0; // counter = cnt1 + NEW1 * cnt2
         for ( int cnt2 = 0; cnt2 < NEW2; cnt2++ ){
            for ( int cnt1 = 0; cnt1 < NEW1; cnt1++ ){
               if (( start <= counter ) && ( counter < stop )){
                  for ( int cnt4 = 0; cnt4 < NEW4; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 < NEW3; cnt3++ ){
                        ROT_TEI->set_coulomb( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                           eri[ ( counter - start ) + SIZE * ( cnt3 + NEW3 * cnt4 ) ] );
                           // Indices (12) and indices (34) are Coulomb pairs
                     }
                  }
               }
               counter++;
            }
         }
         written = true;

      }
   }

   if (( space1 == 'C' ) && ( space2 == 'V' ) && ( space3 == 'C' ) && ( space4 == 'V' )){ // ( C V | C V )

      assert( pack == false );
      assert( Irreps::directProd( irrep1, irrep2 ) == Irreps::directProd( irrep3, irrep4 ) );

      const int NEW1  = dimension( idx, irrep1, space1 );
      const int NEW2  = dimension( idx, irrep2, space2 );
      const int NEW3  = dimension( idx, irrep3, space3 );
      const int NEW4  = dimension( idx, irrep4, space4 );
      const int JUMP2 = jump( idx, irrep2, space2 );
      const int JUMP4 = jump( idx, irrep4, space4 );
      const int SIZE  = stop - start;

      int counter = 0; // counter = cnt1 + NEW1 * cnt2
      for ( int cnt2 = 0; cnt2 < NEW2; cnt2++ ){
         for ( int cnt1 = 0; cnt1 < NEW1; cnt1++ ){
            if (( start <= counter ) && ( counter < stop )){
               for ( int cnt4 = 0; cnt4 < NEW4; cnt4++ ){
                  for ( int cnt3 = 0; cnt3 < NEW3; cnt3++ ){
                     ROT_TEI->set_exchange( irrep1, irrep3, irrep2, irrep4, cnt1, cnt3, JUMP2 + cnt2, JUMP4 + cnt4,
                        eri[ ( counter - start ) + SIZE * ( cnt3 + NEW3 * cnt4 ) ] );
                        // Indices (12) and indices (34) are Coulomb pairs
                  }
               }
            }
            counter++;
         }
      }
      written = true;

   }

   assert( written == true );

}

void CheMPS2::DMRGSCFrotations::blockwise_first( double * origin, double * target, int orig1, int dim2, const int dim34, double * umat1, int new1, int lda1 ){

   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int right_dim = dim2 * dim34;
   dgemm_( &notrans, &notrans, &new1, &right_dim, &orig1, &one, umat1, &lda1, origin, &orig1, &set, target, &new1 );

}

void CheMPS2::DMRGSCFrotations::blockwise_second( double * origin, double * target, int dim1, int orig2, const int dim34, double * umat2, int new2, int lda2 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   const int jump_old  = dim1 * orig2;
   const int jump_new  = dim1 * new2;
   #pragma omp parallel for schedule(static)
   for ( int index = 0; index < dim34; index++ ){
      dgemm_( &notrans, &trans, &dim1, &new2, &orig2, &one, origin + jump_old * index, &dim1, umat2, &lda2, &set, target + jump_new * index, &dim1 );
   }

}

void CheMPS2::DMRGSCFrotations::blockwise_third( double * origin, double * target, const int dim12, int orig3, const int dim4, double * umat3, int new3, int lda3 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim12;
   const int jump_old = dim12 * orig3;
   const int jump_new = dim12 * new3;
   #pragma omp parallel for schedule(static)
   for ( int index = 0; index < dim4; index++ ){
      dgemm_( &notrans, &trans, &left_dim, &new3, &orig3, &one, origin + jump_old * index, &left_dim, umat3, &lda3, &set, target + jump_new * index, &left_dim );
   }

}

void CheMPS2::DMRGSCFrotations::blockwise_fourth( double * origin, double * target, const int dim12, int dim3, int orig4, double * umat4, int new4, int lda4 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim12 * dim3;
   dgemm_( &notrans, &trans, &left_dim, &new4, &orig4, &one, origin, &left_dim, umat4, &lda4, &set, target, &left_dim );

}

void CheMPS2::DMRGSCFrotations::open_file( hid_t * file_id, hid_t * dspc_id, hid_t * dset_id, const int first, const int second, const string filename ){

   file_id[0] = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
   hsize_t fdim_h5[] = { (hsize_t) second, (hsize_t) first }; // C is row major: [ col + ncol * row ] is assumed
   dspc_id[0] = H5Screate_simple( 2, fdim_h5, NULL );
   dset_id[0] = H5Dcreate( file_id[0], "storage", H5T_NATIVE_DOUBLE, dspc_id[0], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

}

void CheMPS2::DMRGSCFrotations::close_file( hid_t file_id, hid_t dspc_id, hid_t dset_id ){

   H5Dclose( dset_id );
   H5Sclose( dspc_id );
   H5Fclose( file_id );

}

void CheMPS2::DMRGSCFrotations::write_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int first_write ){

   hsize_t stride_h5[] = { 1, 1 };
   hsize_t  count_h5[] = { 1, 1 };
   hsize_t  start_h5[] = { (hsize_t) start, 0 };
   hsize_t  block_h5[] = { (hsize_t) size, (hsize_t) first_write };
   H5Sselect_hyperslab( dspc_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
   hsize_t mem_h5 = size * first_write; // Should be OK to multiply as integers as it is smaller than mem_size
   hid_t mem_id = H5Screate_simple( 1, &mem_h5, NULL );
   H5Dwrite( dset_id, H5T_NATIVE_DOUBLE, mem_id, dspc_id, H5P_DEFAULT, eri );
   H5Sclose( mem_id );

}

void CheMPS2::DMRGSCFrotations::read_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int second_read ){

   hsize_t stride_h5[] = { 1, 1 };
   hsize_t  count_h5[] = { 1, 1 };
   hsize_t  start_h5[] = { 0, (hsize_t) start };
   hsize_t  block_h5[] = { (hsize_t) second_read, (hsize_t) size };
   H5Sselect_hyperslab( dspc_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
   hsize_t mem_h5 = second_read * size; // Should be OK to multiply as integers as it is smaller than mem_size
   hid_t mem_id = H5Screate_simple( 1, &mem_h5, NULL );
   H5Dread( dset_id, H5T_NATIVE_DOUBLE, mem_id, dspc_id, H5P_DEFAULT, eri );
   H5Sclose( mem_id );

}

int CheMPS2::DMRGSCFrotations::dimension( DMRGSCFindices * idx, const int irrep, const char space ){

   if ( space == 'O' ){ return idx->getNOCC( irrep ); }
   if ( space == 'A' ){ return idx->getNDMRG( irrep ); }
   if ( space == 'V' ){ return idx->getNVIRT( irrep ); }
   if ( space == 'C' ){ return ( idx->getNOCC( irrep ) + idx->getNDMRG( irrep ) ); }
   if ( space == 'F' ){ return idx->getNORB( irrep ); }
   return -1;

}

int CheMPS2::DMRGSCFrotations::jump( DMRGSCFindices * idx, const int irrep, const char space ){

   if ( space == 'A' ){ return idx->getNOCC( irrep ); }
   if ( space == 'V' ){ return idx->getNOCC( irrep ) + idx->getNDMRG( irrep ); }
   return 0; // O, F, C

}

void CheMPS2::DMRGSCFrotations::unpackage_second( double * mem1, double * mem2, const int SIZE, const int ORIG ){

   #pragma omp parallel for schedule(static)
   for ( int cnt4 = 0; cnt4 < ORIG; cnt4++ ){
      for ( int cnt3 = 0; cnt3 < ORIG; cnt3++ ){
         const int combined = (( cnt3 < cnt4 ) ? ( cnt3 + ( cnt4 * ( cnt4 + 1 )) / 2 )
                                               : ( cnt4 + ( cnt3 * ( cnt3 + 1 )) / 2 ));
         #pragma omp simd
         for ( int cnt12 = 0; cnt12 < SIZE; cnt12++ ){
            mem2[ cnt12 + SIZE * ( cnt3 + ORIG * cnt4 ) ] = mem1[ cnt12 + SIZE * combined ];
         }
      }
   }

}

void CheMPS2::DMRGSCFrotations::package_first( double * mem1, double * mem2, const int NEW, const int PACKED, const int SIZE ){

   #pragma omp parallel for schedule(static)
   for ( int cnt34 = 0; cnt34 < SIZE; cnt34++ ){
      for ( int cnt2 = 0; cnt2 < NEW; cnt2++ ){
         #pragma omp simd
         for ( int cnt1 = 0; cnt1 <= cnt2; cnt1++ ){
            mem2[ cnt1 + ( cnt2 * ( cnt2 + 1 ))/2 + PACKED * cnt34 ] = mem1[ cnt1 + NEW * ( cnt2 + NEW * cnt34 ) ];
         }
      }
   }

}

void CheMPS2::DMRGSCFrotations::rotate( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space1, const char space2, const char space3, const char space4, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename ){

   /* Matrix elements ( 1 2 | 3 4 ) */

   assert(( space1 == 'O' ) || ( space1 == 'A' ) || ( space1 == 'V' ) || ( space1 == 'C' ) || ( space1 == 'F' ));
   assert(( space2 == 'O' ) || ( space2 == 'A' ) || ( space2 == 'V' ) || ( space2 == 'C' ) || ( space2 == 'F' ));
   assert(( space3 == 'O' ) || ( space3 == 'A' ) || ( space3 == 'V' ) || ( space3 == 'C' ) || ( space3 == 'F' ));
   assert(( space4 == 'O' ) || ( space4 == 'A' ) || ( space4 == 'V' ) || ( space4 == 'C' ) || ( space4 == 'F' ));

   const int num_irreps = idx->getNirreps();
   const bool equal12 = ( space1 == space2 );
   const bool equal34 = ( space3 == space4 );
   const bool eightfold = (( space1 == space3 ) && ( space2 == space4 ));

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      for ( int irrep2 = (( equal12 ) ? irrep1 : 0 ); irrep2 < num_irreps; irrep2++ ){ // irrep2 >= irrep1 if space1 == space2
         const int product_symm = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = (( eightfold ) ? irrep1 : 0 ); irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( product_symm, irrep3 );
            if ( irrep4 >= (( equal34 ) ? irrep3 : 0 ) ){ // irrep4 >= irrep3 if space3 == space4

               const int NEW1 = dimension( idx, irrep1, space1 );
               const int NEW2 = dimension( idx, irrep2, space2 );
               const int NEW3 = dimension( idx, irrep3, space3 );
               const int NEW4 = dimension( idx, irrep4, space4 );

               if (( NEW1 > 0 ) && ( NEW2 > 0 ) && ( NEW3 > 0 ) && ( NEW4 > 0 )){

                  const int ORIG1 = idx->getNORB( irrep1 );
                  const int ORIG2 = idx->getNORB( irrep2 );
                  const int ORIG3 = idx->getNORB( irrep3 );
                  const int ORIG4 = idx->getNORB( irrep4 );

                  double * umat1 = umat->getBlock( irrep1 ) + jump( idx, irrep1, space1 );
                  double * umat2 = umat->getBlock( irrep2 ) + jump( idx, irrep2, space2 );
                  double * umat3 = umat->getBlock( irrep3 ) + jump( idx, irrep3, space3 );
                  double * umat4 = umat->getBlock( irrep4 ) + jump( idx, irrep4, space4 );

                  const int block_size1 = mem_size / ( ORIG1 * ORIG2 ); // Floor of amount of times orig( first  ) fits in mem_size
                  const int block_size2 = mem_size / ( ORIG3 * ORIG4 ); // Floor of amount of times orig( second ) fits in mem_size
                  assert( block_size1 > 0 );
                  assert( block_size2 > 0 );

                  const bool pack_first  = (( equal12 ) && ( irrep1 == irrep2 ));
                  const bool pack_second = (( equal34 ) && ( irrep3 == irrep4 ));
                  const int   first_size = (( pack_first  ) ? (  NEW1 * (  NEW1 + 1 )) / 2 :  NEW1 * NEW2  );
                  const int  second_size = (( pack_second ) ? ( ORIG3 * ( ORIG3 + 1 )) / 2 : ORIG3 * ORIG4 );

                  const bool io_free = (( block_size1 >= second_size ) && ( block_size2 >= first_size ));
                  hid_t file_id, dspc_id, dset_id;

                  if ( io_free == false ){
                     assert( filename.compare( "edmistonruedenberg" ) != 0 );
                     open_file( &file_id, &dspc_id, &dset_id, first_size, second_size, filename );
                  }

                  // First half transformation
                  int start = 0;
                  while ( start < second_size ){
                     const int stop = min( start + block_size1, second_size );
                     const int size = stop - start;
                     fetch( mem1, ORIG_VMAT, irrep1, irrep2, irrep3, irrep4, idx, start, stop, pack_second );
                     blockwise_first(  mem1, mem2, ORIG1, ORIG2, size, umat1, NEW1, ORIG1 );
                     blockwise_second( mem2, mem1, NEW1,  ORIG2, size, umat2, NEW2, ORIG2 );
                     if ( pack_first ){
                        package_first( mem1, mem2, NEW1, first_size, size );
                        double * temp = mem1;
                        mem1 = mem2;
                        mem2 = temp;
                     }
                     if ( io_free == false ){ write_file( dspc_id, dset_id, mem1, start, size, first_size ); }
                     start += size;
                  }
                  assert( start == second_size );

                  // Do the second half transformation
                  start = 0;
                  while ( start < first_size ){
                     const int stop = min( start + block_size2, first_size );
                     const int size = stop - start;
                     if ( io_free == false ){ read_file( dspc_id, dset_id, mem1, start, size, second_size ); }
                     if ( pack_second ){
                        unpackage_second( mem1, mem2, size, ORIG3 );
                        double * temp = mem1;
                        mem1 = mem2;
                        mem2 = temp;
                     }
                     blockwise_fourth( mem1, mem2, size, ORIG3, ORIG4, umat4, NEW4, ORIG4 );
                     blockwise_third(  mem2, mem1, size, ORIG3, NEW4,  umat3, NEW3, ORIG3 );
                     write( mem1, NEW_VMAT, ROT_TEI, space1, space2, space3, space4, irrep1, irrep2, irrep3, irrep4, idx, start, stop, pack_first );
                     start += size;
                  }
                  assert( start == first_size );
                  if ( io_free == false ){ close_file( file_id, dspc_id, dset_id ); }
               }
            }
         }
      }
   }

}


