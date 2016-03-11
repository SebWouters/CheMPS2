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
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <string>

#include "DMRGSCFVmatRotations.h"
#include "Lapack.h"
#include "MPIchemps2.h"
#include "MyHDF5.h"

using std::min;
using std::max;
using std::string;
using std::cout;
using std::endl;

CheMPS2::DMRGSCFVmatRotations::DMRGSCFVmatRotations(){ }
CheMPS2::DMRGSCFVmatRotations::~DMRGSCFVmatRotations(){ }

void CheMPS2::DMRGSCFVmatRotations::blockwise_first( double * origin, double * target, int orig1, int dim2, int dim3, int dim4, double * umat1, int new1, int lda1 ){

   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int right_dim = dim2 * dim3 * dim4;
   dgemm_( &notrans, &notrans, &new1, &right_dim, &orig1, &one, umat1, &lda1, origin, &orig1, &set, target, &new1 );

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_second( double * origin, double * target, int dim1, int orig2, int dim3, int dim4, double * umat2, int new2, int lda2 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   const int right_dim = dim3 * dim4;
   const int jump_old  = dim1 * orig2;
   const int jump_new  = dim1 * new2;
   for ( int index = 0; index < right_dim; index++ ){
      dgemm_( &notrans, &trans, &dim1, &new2, &orig2, &one, origin + jump_old * index, &dim1, umat2, &lda2, &set, target + jump_new * index, &dim1 );
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_third( double * origin, double * target, int dim1, int dim2, int orig3, int dim4, double * umat3, int new3, int lda3 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim1 * dim2;
   const int jump_old = dim1 * dim2 * orig3;
   const int jump_new = dim1 * dim2 * new3;
   for ( int index = 0; index < dim4; index++ ){
      dgemm_( &notrans, &trans, &left_dim, &new3, &orig3, &one, origin + jump_old * index, &left_dim, umat3, &lda3, &set, target + jump_new * index, &left_dim );
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_fourth( double * origin, double * target, int dim1, int dim2, int dim3, int orig4, double * umat4, int new4, int lda4 ){

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   int left_dim = dim1 * dim2 * dim3;
   dgemm_( &notrans, &trans, &left_dim, &new4, &orig4, &one, origin, &left_dim, umat4, &lda4, &set, target, &left_dim );

}

void CheMPS2::DMRGSCFVmatRotations::full_base( double * eri, double * work, double * umat1, int new1, int orig1,
                                                                            double * umat2, int new2, int orig2,
                                                                            double * umat3, int new3, int orig3,
                                                                            double * umat4, int new4, int orig4 ){

   blockwise_first(  eri, work, orig1, orig2, orig3, orig4, umat1, new1, orig1 ); // (ijkl) -> (ajkl)
   blockwise_fourth( work, eri, new1,  orig2, orig3, orig4, umat4, new4, orig4 ); // (ajkl) -> (ajkd)
   blockwise_third(  eri, work, new1,  orig2, orig3, new4,  umat3, new3, orig3 ); // (ajkd) -> (ajcd)
   blockwise_second( work, eri, new1,  orig2, new3,  new4,  umat2, new2, orig2 ); // (ajcd) -> (abcd)

}

void CheMPS2::DMRGSCFVmatRotations::full( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, const char space, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * temp1, double * temp2 ){

   assert(( space == 'A' ) || ( space == 'F' )); // Active of full
   const int num_irreps = idx->getNirreps();

   //Two-body terms --> use eightfold permutation symmetry in the irreps
   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      for ( int irrep2 = irrep1; irrep2 < num_irreps; irrep2++ ){
         const int product_symm = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = irrep1; irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( product_symm, irrep3 );
            if ( irrep4 >= irrep2 ){

               const int NEW1 = (( space == 'A' ) ? idx->getNDMRG( irrep1 ) : idx->getNORB( irrep1 ) );
               const int NEW2 = (( space == 'A' ) ? idx->getNDMRG( irrep2 ) : idx->getNORB( irrep2 ) );
               const int NEW3 = (( space == 'A' ) ? idx->getNDMRG( irrep3 ) : idx->getNORB( irrep3 ) );
               const int NEW4 = (( space == 'A' ) ? idx->getNDMRG( irrep4 ) : idx->getNORB( irrep4 ) );

               if (( NEW1 > 0 ) && ( NEW2 > 0 ) && ( NEW3 > 0 ) && ( NEW4 > 0 )){

                  const int NORB1 = idx->getNORB( irrep1 );
                  const int NORB2 = idx->getNORB( irrep2 );
                  const int NORB3 = idx->getNORB( irrep3 );
                  const int NORB4 = idx->getNORB( irrep4 );

                  for ( int cnt4 = 0; cnt4 < NORB4; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 < NORB3; cnt3++ ){
                        for ( int cnt2 = 0; cnt2 < NORB2; cnt2++ ){
                           for ( int cnt1 = 0; cnt1 < NORB1; cnt1++ ){
                              temp1[ cnt1 + NORB1 * ( cnt2 + NORB2 * ( cnt3 + NORB3 * cnt4 ) ) ]
                                = ORIG_VMAT->get( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4 );
                           }
                        }
                     }
                  }

                  double * umat1 = umat->getBlock( irrep1 ) + (( space == 'A' ) ? idx->getNOCC( irrep1 ) : 0 );
                  double * umat2 = umat->getBlock( irrep2 ) + (( space == 'A' ) ? idx->getNOCC( irrep2 ) : 0 );
                  double * umat3 = umat->getBlock( irrep3 ) + (( space == 'A' ) ? idx->getNOCC( irrep3 ) : 0 );
                  double * umat4 = umat->getBlock( irrep4 ) + (( space == 'A' ) ? idx->getNOCC( irrep4 ) : 0 );

                  full_base( temp1, temp2, umat1, NEW1, NORB1,
                                           umat2, NEW2, NORB2,
                                           umat3, NEW3, NORB3,
                                           umat4, NEW4, NORB4 );

                  for ( int cnt4 = 0; cnt4 < NEW4; cnt4++ ){
                     for ( int cnt3 = 0; cnt3 < NEW3; cnt3++ ){
                        for ( int cnt2 = 0; cnt2 < NEW2; cnt2++ ){
                           for ( int cnt1 = 0; cnt1 < NEW1; cnt1++ ){
                              NEW_VMAT->set( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                                             temp1[ cnt1 + NEW1 * ( cnt2 + NEW2 * ( cnt3 + NEW3 * cnt4 ) ) ] );
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

void CheMPS2::DMRGSCFVmatRotations::full( const FourIndex * ORIG_VMAT, DMRGSCFintegrals * ROT_TEI, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * temp1, double * temp2 ){

   const int num_irreps = idx->getNirreps();

   // First do Coulomb object : ( c1 <= c2 | a1 <= a2 )
   for ( int Ic1 = 0; Ic1 < num_irreps; Ic1++){
      for ( int Ic2 = Ic1; Ic2 < num_irreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for ( int Ia1 = 0; Ia1 < num_irreps; Ia1++ ){
            const int Ia2 = Irreps::directProd( Ia1, Icc );
            if ( Ia1 <= Ia2 ){

               const int NEW_C1 = idx->getNOCC( Ic1 ) + idx->getNDMRG( Ic1 );
               const int NEW_C2 = idx->getNOCC( Ic2 ) + idx->getNDMRG( Ic2 );
               const int NEW_A1 = idx->getNORB( Ia1 );
               const int NEW_A2 = idx->getNORB( Ia2 );

               if (( NEW_C1 > 0 ) && ( NEW_C2 > 0 ) && ( NEW_A1 > 0 ) && ( NEW_A2 > 0 )){

                  const int NORB_C1 = idx->getNORB( Ic1 );
                  const int NORB_C2 = idx->getNORB( Ic2 );
                  const int NORB_A1 = idx->getNORB( Ia1 );
                  const int NORB_A2 = idx->getNORB( Ia2 );

                  for ( int c2 = 0; c2 < NORB_C2; c2++ ){
                     for ( int a2 = 0; a2 < NORB_A2; a2++ ){
                        for ( int a1 = 0; a1 < NORB_A1; a1++ ){
                           for ( int c1 = 0; c1 < NORB_C1; c1++ ){
                              // We try to make the Coulomb elements !
                              temp1[ c1 + NORB_C1 * ( a1 + NORB_A1 * ( a2 + NORB_A2 * c2 ) ) ]
                                = ORIG_VMAT->get( Ic1, Ia1, Ic2, Ia2, c1, a1, c2, a2 );
                           }
                        }
                     }
                  }

                  full_base( temp1, temp2, umat->getBlock( Ic1 ), NEW_C1, NORB_C1,
                                           umat->getBlock( Ia1 ), NEW_A1, NORB_A1,
                                           umat->getBlock( Ia2 ), NEW_A2, NORB_A2,
                                           umat->getBlock( Ic2 ), NEW_C2, NORB_C2 );

                  for ( int c2 = 0; c2 < NEW_C2; c2++ ){
                     for ( int a2 = 0; a2 < NEW_A2; a2++ ){
                        for ( int a1 = 0; a1 < NEW_A1; a1++ ){
                           for ( int c1 = 0; c1 < NEW_C1; c1++ ){
                              ROT_TEI->set_coulomb( Ic1, Ic2, Ia1, Ia2, c1, c2, a1, a2,
                                  temp1[ c1 + NEW_C1 * ( a1 + NEW_A1 * ( a2 + NEW_A2 * c2 ) ) ] );
                           }
                        }
                     }
                  }
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

               for ( int c2 = 0; c2 < NORB_C2; c2++ ){
                  for ( int v2 = 0; v2 < NORB_V2; v2++ ){
                     for ( int v1 = 0; v1 < NORB_V1; v1++ ){
                        for ( int c1 = 0; c1 < NORB_C1; c1++ ){
                           // We try to make the Exchange elements !
                           temp1[ c1 + NORB_C1 * ( v1 + NORB_V1 * ( v2 + NORB_V2 * c2 ) ) ]
                             = ORIG_VMAT->get( Ic1, Ic2, Iv1, Iv2, c1, c2, v1, v2 );
                        }
                     }
                  }
               }

               const int JUMP_V1 = idx->getNOCC( Iv1 ) + idx->getNDMRG( Iv1 );
               const int JUMP_V2 = idx->getNOCC( Iv2 ) + idx->getNDMRG( Iv2 );

               full_base( temp1, temp2, umat->getBlock( Ic1 ),           NEW_C1, NORB_C1,
                                        umat->getBlock( Iv1 ) + JUMP_V1, NEW_V1, NORB_V1,
                                        umat->getBlock( Iv2 ) + JUMP_V2, NEW_V2, NORB_V2,
                                        umat->getBlock( Ic2 ),           NEW_C2, NORB_C2 );

               for ( int c2 = 0; c2 < NEW_C2; c2++ ){
                  for ( int v2 = 0; v2 < NEW_V2; v2++ ){
                     for ( int v1 = 0; v1 < NEW_V1; v1++ ){
                        for ( int c1 = 0; c1 < NEW_C1; c1++ ){
                           ROT_TEI->set_exchange( Ic1, Ic2, Iv1, Iv2, c1, c2, JUMP_V1 + v1, JUMP_V2 + v2,
                               temp1[ c1 + NEW_C1 * ( v1 + NEW_V1 * ( v2 + NEW_V2 * c2 ) ) ] );
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_disk( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, const char space, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename ){

   assert(( space == 'A' ) || ( space == 'F' ));
   const int num_irreps = idx->getNirreps();

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      for ( int irrep2 = irrep1; irrep2 < num_irreps; irrep2++ ){
         const int product_symm = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = irrep1; irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( product_symm, irrep3 );
            if ( irrep4 >= irrep2 ){

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

                  {  // Do the first half transformation
                     hid_t   file_id      = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                     hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                     hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                     hid_t   dataset_id   = H5Dcreate( file_id, "storage", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

                     const int block_size = mem_size / ( ORIG1 * ORIG2 ); // Floor of amount of times first_old fits in mem_size
                     assert( block_size > 0 );

                     int start = 0;
                     while ( start < second_old ){

                        const int stop = min( start + block_size, second_old );
                        const int size = stop - start;
                        assert( size > 0 );

                        // count = ( cnt3 + ORIG3 * cnt4 )
                        for ( int count = start; count < stop; count++ ){
                           const int cnt3 = count % ORIG3;
                           const int cnt4 = count / ORIG3;
                           const int rel  = count - start;
                           for ( int cnt2 = 0; cnt2 < ORIG2; cnt2++ ){
                              for ( int cnt1 = 0; cnt1 < ORIG1; cnt1++ ){
                                 mem1[ cnt1 + ORIG1 * ( cnt2 + ORIG2 * rel ) ]
                                    = ORIG_VMAT->get( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4 );
                              }
                           }
                        }

                        blockwise_first(  mem1, mem2, ORIG1, ORIG2, size, 1, umat1, NEW1, ORIG1 );
                        blockwise_second( mem2, mem1, NEW1,  ORIG2, size, 1, umat2, NEW2, ORIG2 );

                        {
                           hsize_t stride_h5[] = { 1, 1 };
                           hsize_t  count_h5[] = { 1, 1 };
                           hsize_t  start_h5[] = { start, 0         }; // C is row major: [ col + ncol * row ] is assumed
                           hsize_t  block_h5[] = { size,  first_new }; // C is row major: [ col + ncol * row ] is assumed
                           H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                        }

                        hsize_t memsize_h5 = size * fdim_h5[ 1 ];
                        hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                        H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                        H5Sclose( memspace_id );

                        start += size;
                     }
                     assert( start == second_old );

                     H5Dclose( dataset_id );
                     H5Sclose( dataspace_id );
                     H5Fclose( file_id );
                  }
                  {  // Do the second half transformation
                     hid_t   file_id      = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
                     hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                     hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                     hid_t   dataset_id   = H5Dopen( file_id, "storage", H5P_DEFAULT );

                     const int block_size = mem_size / second_old; // Floor of amount of times second_old fits in mem_size
                     assert( block_size > 0 );

                     int start = 0;
                     while ( start < first_new ){

                        const int stop = min( start + block_size, first_new );
                        const int size = stop - start;
                        assert( size > 0 );

                        {
                           hsize_t stride_h5[] = { 1, 1 };
                           hsize_t  count_h5[] = { 1, 1 };
                           hsize_t  start_h5[] = { 0,          start }; // C is row major: [ col + ncol * row ] is assumed
                           hsize_t  block_h5[] = { second_old, size  }; // C is row major: [ col + ncol * row ] is assumed
                           H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                        }

                        hsize_t memsize_h5 = fdim_h5[ 0 ] * size;
                        hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                        H5Dread( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                        H5Sclose( memspace_id );

                        blockwise_fourth( mem1, mem2, 1, size, ORIG3, ORIG4, umat4, NEW4, ORIG4 );
                        blockwise_third(  mem2, mem1, 1, size, ORIG3, NEW4,  umat3, NEW3, ORIG3 );

                        // count = ( cnt1 + NEW1 * cnt2 )
                        for ( int count = start; count < stop; count++ ){
                           const int cnt1 = count % NEW1;
                           const int cnt2 = count / NEW1;
                           const int rel  = count - start;
                           for ( int cnt4 = 0; cnt4 < NEW4; cnt4++ ){
                              for ( int cnt3 = 0; cnt3 < NEW3; cnt3++ ){
                                 NEW_VMAT->set( irrep1, irrep2, irrep3, irrep4, cnt1, cnt2, cnt3, cnt4,
                                                mem1[ rel + size * ( cnt3 + NEW3 * cnt4 ) ] );
                              }
                           }
                        }

                        start += size;
                     }
                     assert( start == first_new );

                     H5Dclose( dataset_id );
                     H5Sclose( dataspace_id );
                     H5Fclose( file_id );
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::DMRGSCFVmatRotations::blockwise_disk( const FourIndex * ORIG_VMAT, DMRGSCFintegrals * ROT_TEI, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename ){

   const int num_irreps = idx->getNirreps();

   // First do Coulomb object : ( c1 <= c2 | a1 <= a2 )
   for ( int Ic1 = 0; Ic1 < num_irreps; Ic1++){
      for ( int Ic2 = Ic1; Ic2 < num_irreps; Ic2++){
         const int Icc = Irreps::directProd( Ic1, Ic2 );
         for ( int Ia1 = 0; Ia1 < num_irreps; Ia1++ ){
            const int Ia2 = Irreps::directProd( Ia1, Icc );
            if ( Ia1 <= Ia2 ){

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

                  {  // Do the first half transformation
                     hid_t   file_id      = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                     hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                     hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                     hid_t   dataset_id   = H5Dcreate( file_id, "storage", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

                     const int block_size = mem_size / ( NORB_C1 * NORB_C2 ); // Floor of amount of times first_old fits in mem_size
                     assert( block_size > 0 );

                     int start = 0;
                     while ( start < second_old ){

                        const int stop = min( start + block_size, second_old );
                        const int size = stop - start;
                        assert( size > 0 );

                        // count = ( cnt_A1 + NORB_A1 * cnt_A2 )
                        for ( int count = start; count < stop; count++ ){
                           const int cnt_A1 = count % NORB_A1;
                           const int cnt_A2 = count / NORB_A1;
                           const int rel  = count - start;
                           for ( int cnt_C2 = 0; cnt_C2 < NORB_C2; cnt_C2++ ){
                              for ( int cnt_C1 = 0; cnt_C1 < NORB_C1; cnt_C1++ ){
                                 mem1[ cnt_C1 + NORB_C1 * ( cnt_C2 + NORB_C2 * rel ) ]
                                    = ORIG_VMAT->get( Ic1, Ia1, Ic2, Ia2, cnt_C1, cnt_A1, cnt_C2, cnt_A2 );
                              }
                           }
                        }

                        blockwise_first(  mem1, mem2, NORB_C1, NORB_C2, size, 1, umat_C1, NEW_C1, NORB_C1 );
                        blockwise_second( mem2, mem1,  NEW_C1, NORB_C2, size, 1, umat_C2, NEW_C2, NORB_C2 );

                        {
                           hsize_t stride_h5[] = { 1, 1 };
                           hsize_t  count_h5[] = { 1, 1 };
                           hsize_t  start_h5[] = { start, 0         }; // C is row major: [ col + ncol * row ] is assumed
                           hsize_t  block_h5[] = { size,  first_new }; // C is row major: [ col + ncol * row ] is assumed
                           H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                        }

                        hsize_t memsize_h5 = size * fdim_h5[ 1 ];
                        hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                        H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                        H5Sclose( memspace_id );

                        start += size;
                     }
                     assert( start == second_old );

                     H5Dclose( dataset_id );
                     H5Sclose( dataspace_id );
                     H5Fclose( file_id );
                  }
                  {  // Do the second half transformation
                     hid_t   file_id      = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
                     hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                     hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                     hid_t   dataset_id   = H5Dopen( file_id, "storage", H5P_DEFAULT );

                     const int block_size = mem_size / second_old; // Floor of amount of times second_old fits in mem_size
                     assert( block_size > 0 );

                     int start = 0;
                     while ( start < first_new ){

                        const int stop = min( start + block_size, first_new );
                        const int size = stop - start;
                        assert( size > 0 );

                        {
                           hsize_t stride_h5[] = { 1, 1 };
                           hsize_t  count_h5[] = { 1, 1 };
                           hsize_t  start_h5[] = { 0,          start }; // C is row major: [ col + ncol * row ] is assumed
                           hsize_t  block_h5[] = { second_old, size  }; // C is row major: [ col + ncol * row ] is assumed
                           H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                        }

                        hsize_t memsize_h5 = fdim_h5[ 0 ] * size;
                        hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                        H5Dread( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                        H5Sclose( memspace_id );

                        blockwise_fourth( mem1, mem2, 1, size, NORB_A1, NORB_A2, umat_A2, NEW_A2, NORB_A2 );
                        blockwise_third(  mem2, mem1, 1, size, NORB_A1,  NEW_A2, umat_A1, NEW_A1, NORB_A1 );

                        // count = ( cnt_C1 + NEW_C1 * cnt_C2 )
                        for ( int count = start; count < stop; count++ ){
                           const int cnt_C1 = count % NEW_C1;
                           const int cnt_C2 = count / NEW_C1;
                           const int rel  = count - start;
                           for ( int cnt_A2 = 0; cnt_A2 < NEW_A2; cnt_A2++ ){
                              for ( int cnt_A1 = 0; cnt_A1 < NEW_A1; cnt_A1++ ){
                                 ROT_TEI->set_coulomb( Ic1, Ic2, Ia1, Ia2, cnt_C1, cnt_C2, cnt_A1, cnt_A2,
                                    mem1[ rel + size * ( cnt_A1 + NEW_A1 * cnt_A2 ) ] );
                              }
                           }
                        }

                        start += size;
                     }
                     assert( start == first_new );

                     H5Dclose( dataset_id );
                     H5Sclose( dataspace_id );
                     H5Fclose( file_id );
                  }
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

               const int  first_new =  NEW_C1 * NEW_C2;
               const int second_old = NORB_V1 * NORB_V2;

               {  // Do the first half transformation
                  hid_t   file_id      = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                  hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                  hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                  hid_t   dataset_id   = H5Dcreate( file_id, "storage", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

                  const int block_size = mem_size / ( NORB_C1 * NORB_C2 ); // Floor of amount of times first_old fits in mem_size
                  assert( block_size > 0 );

                  int start = 0;
                  while ( start < second_old ){

                     const int stop = min( start + block_size, second_old );
                     const int size = stop - start;
                     assert( size > 0 );

                     // count = ( cnt_V1 + NORB_V1 * cnt_V2 )
                     for ( int count = start; count < stop; count++ ){
                        const int cnt_V1 = count % NORB_V1;
                        const int cnt_V2 = count / NORB_V1;
                        const int rel  = count - start;
                        for ( int cnt_C2 = 0; cnt_C2 < NORB_C2; cnt_C2++ ){
                           for ( int cnt_C1 = 0; cnt_C1 < NORB_C1; cnt_C1++ ){
                              mem1[ cnt_C1 + NORB_C1 * ( cnt_C2 + NORB_C2 * rel ) ]
                                 = ORIG_VMAT->get( Ic1, Ic2, Iv1, Iv2, cnt_C1, cnt_C2, cnt_V1, cnt_V2 );
                           }
                        }
                     }

                     blockwise_first(  mem1, mem2, NORB_C1, NORB_C2, size, 1, umat_C1, NEW_C1, NORB_C1 );
                     blockwise_second( mem2, mem1,  NEW_C1, NORB_C2, size, 1, umat_C2, NEW_C2, NORB_C2 );

                     {
                        hsize_t stride_h5[] = { 1, 1 };
                        hsize_t  count_h5[] = { 1, 1 };
                        hsize_t  start_h5[] = { start, 0         }; // C is row major: [ col + ncol * row ] is assumed
                        hsize_t  block_h5[] = { size,  first_new }; // C is row major: [ col + ncol * row ] is assumed
                        H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                     }

                     hsize_t memsize_h5 = size * fdim_h5[ 1 ];
                     hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                     H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                     H5Sclose( memspace_id );

                     start += size;
                  }
                  assert( start == second_old );

                  H5Dclose( dataset_id );
                  H5Sclose( dataspace_id );
                  H5Fclose( file_id );
               }
               {  // Do the second half transformation
                  hid_t   file_id      = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
                  hsize_t fdim_h5[]    = { second_old, first_new }; // C is row major: [ col + ncol * row ] is assumed
                  hid_t   dataspace_id = H5Screate_simple( 2, fdim_h5, NULL );
                  hid_t   dataset_id   = H5Dopen( file_id, "storage", H5P_DEFAULT );

                  const int block_size = mem_size / second_old; // Floor of amount of times second_old fits in mem_size
                  assert( block_size > 0 );

                  int start = 0;
                  while ( start < first_new ){

                     const int stop = min( start + block_size, first_new );
                     const int size = stop - start;
                     assert( size > 0 );

                     {
                        hsize_t stride_h5[] = { 1, 1 };
                        hsize_t  count_h5[] = { 1, 1 };
                        hsize_t  start_h5[] = { 0,          start }; // C is row major: [ col + ncol * row ] is assumed
                        hsize_t  block_h5[] = { second_old, size  }; // C is row major: [ col + ncol * row ] is assumed
                        H5Sselect_hyperslab( dataspace_id, H5S_SELECT_SET, start_h5, stride_h5, count_h5, block_h5 );
                     }

                     hsize_t memsize_h5 = fdim_h5[ 0 ] * size;
                     hid_t memspace_id  = H5Screate_simple( 1, &memsize_h5, NULL );

                     H5Dread( dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, mem1 );
                     H5Sclose( memspace_id );

                     blockwise_fourth( mem1, mem2, 1, size, NORB_V1, NORB_V2, umat_V2, NEW_V2, NORB_V2 );
                     blockwise_third(  mem2, mem1, 1, size, NORB_V1,  NEW_V2, umat_V1, NEW_V1, NORB_V1 );

                     // count = ( cnt_C1 + NEW_C1 * cnt_C2 )
                     for ( int count = start; count < stop; count++ ){
                        const int cnt_C1 = count % NEW_C1;
                        const int cnt_C2 = count / NEW_C1;
                        const int rel  = count - start;
                        for ( int cnt_V2 = 0; cnt_V2 < NEW_V2; cnt_V2++ ){
                           for ( int cnt_V1 = 0; cnt_V1 < NEW_V1; cnt_V1++ ){
                              ROT_TEI->set_exchange( Ic1, Ic2, Iv1, Iv2, cnt_C1, cnt_C2, JUMP_V1 + cnt_V1, JUMP_V2 + cnt_V2,
                                 mem1[ rel + size * ( cnt_V1 + NEW_V1 * cnt_V2 ) ] );
                           }
                        }
                     }

                     start += size;
                  }
                  assert( start == first_new );

                  H5Dclose( dataset_id );
                  H5Sclose( dataspace_id );
                  H5Fclose( file_id );
               }
            }
         }
      }
   }

}



