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
#include <math.h>

#include "Lapack.h"
#include "DMRGSCFunitary.h"

using std::cout;
using std::endl;

CheMPS2::DMRGSCFunitary::DMRGSCFunitary( const DMRGSCFindices * iHandler ) : DMRGSCFmatrix( iHandler ){

   this->identity();

   //Find the unique indices for OCC-ACT, OCC-VIRT, and ACT-VIRT rotations
   x_linearlength = 0;
   jumper = new int*[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      jumper[ irrep ] = new int[ 3 ];
      const int NOCC = iHandler->getNOCC( irrep );
      const int NACT = iHandler->getNDMRG( irrep );
      const int NVIR = iHandler->getNVIRT( irrep );
      jumper[ irrep ][ 0 ] = x_linearlength;
      x_linearlength += NOCC * NACT;
      jumper[ irrep ][ 1 ] = x_linearlength;
      x_linearlength += NACT * NVIR;
      jumper[ irrep ][ 2 ] = x_linearlength;
      x_linearlength += NOCC * NVIR;
   }

}

CheMPS2::DMRGSCFunitary::~DMRGSCFunitary(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ delete [] jumper[ irrep ]; }
   delete [] jumper;

}

int CheMPS2::DMRGSCFunitary::getNumVariablesX() const{ return x_linearlength; }

void CheMPS2::DMRGSCFunitary::buildSkewSymmX( const int irrep, double * result, double * Xelem, const bool compact ) const{

   const int linsize = iHandler->getNORB( irrep );
   const int NOCC = iHandler->getNOCC( irrep );
   const int NACT = iHandler->getNDMRG( irrep );
   const int NVIR = iHandler->getNVIRT( irrep );
   for ( int cnt = 0; cnt < linsize * linsize; cnt++ ){ result[ cnt ] = 0.0; }

   if ( compact ){

      for ( int occ = 0; occ < NOCC; occ++ ){
         for ( int act = 0; act < NACT; act++ ){
            const int xsolindex = jumper[ irrep ][ 0 ] + act + NACT * occ;
            const int index1 = NOCC + act;
            result[ index1 + linsize * occ ] =   Xelem[ xsolindex ];
            result[ occ + linsize * index1 ] = - Xelem[ xsolindex ];
         }
      }
      for ( int act = 0; act < NACT; act++ ){
         for ( int vir = 0; vir < NVIR; vir++ ){
            const int xsolindex = jumper[ irrep ][ 1 ] + vir + NVIR * act;
            const int index1 = NOCC + NACT + vir;
            const int index2 = NOCC + act;
            result[ index1 + linsize * index2 ] =   Xelem[ xsolindex ];
            result[ index2 + linsize * index1 ] = - Xelem[ xsolindex ];
         }
      }
      for ( int occ = 0; occ < NOCC; occ++ ){
         for ( int vir = 0; vir < NVIR; vir++ ){
            const int xsolindex = jumper[ irrep ][ 2 ] + vir + NVIR * occ;
            const int index1 = NOCC + NACT + vir;
            result[ index1 + linsize * occ ] =   Xelem[ xsolindex ];
            result[ occ + linsize * index1 ] = - Xelem[ xsolindex ];
         }
      }

   } else { //NOT compact

      int jump = 0;
      for ( int cnt = 0; cnt < irrep; cnt++ ){
         const int NORBx = iHandler->getNORB( cnt );
         jump += ( NORBx * ( NORBx - 1 ) ) / 2;
      }

      for ( int row = 0; row < linsize; row++ ){
         for ( int col = row+1; col < linsize; col++ ){
            result[ row + linsize * col ] =   Xelem[ jump + row + ( col * ( col - 1 ) ) / 2 ];
            result[ col + linsize * row ] = - Xelem[ jump + row + ( col * ( col - 1 ) ) / 2 ];
         }
      }

   }

}

void CheMPS2::DMRGSCFunitary::updateUnitary( double * temp1, double * temp2, double * vector, const bool multiply, const bool compact ){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;

      if ( linsize > 1 ){ //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok

         double * xblock    = temp1;              //linsize*linsize
         double * Bmat      = temp1 +     size;   //linsize*linsize
         double * work1     = temp1 + 2 * size;   //linsize*linsize
         double * work2     = temp1 + 3 * size;   //linsize*linsize

         double * workLARGE = temp2;    //4*size
         int     lworkLARGE = 4 * size; //4*size = 4*linsize*linsize > 3*linsize-1

         // Construct the antisymmetric x-matrix
         buildSkewSymmX( irrep, xblock, vector, compact );

         //Bmat <= xblock * xblock
         char notrans = 'N';
         double one = 1.0;
         double set = 0.0;
         dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &one, xblock, &linsize, xblock, &linsize, &set, Bmat, &linsize );

         //Bmat * work1 * Bmat^T <= xblock * xblock
         char uplo = 'U';
         char jobz = 'V';
         int info;
         dsyev_( &jobz, &uplo, &linsize, Bmat, &linsize, work1, workLARGE, &lworkLARGE, &info );

         //work2 <= Bmat^T * xblock * Bmat
         dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &one, xblock, &linsize, Bmat, &linsize, &set, work1, &linsize );
         char trans = 'T';
         dgemm_( &trans, &notrans, &linsize, &linsize, &linsize, &one, Bmat, &linsize, work1, &linsize, &set, work2, &linsize );

         if (CheMPS2::DMRGSCF_debugPrint){
            cout << "   DMRGSCFunitary::updateUnitary : Lambdas of irrep block " << irrep << " : " << endl;
            for (int cnt=0; cnt<linsize/2; cnt++){
               cout << "      block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
               cout << "              [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
            }
         }

         //work1 <= values of the antisymmetric 2x2 blocks
         for ( int cnt = 0; cnt < linsize/2; cnt++ ){
            work1[cnt] = 0.5 * ( work2[ 2 * cnt     + linsize * ( 2 * cnt + 1 ) ]
                               - work2[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] );
         }

         if ( CheMPS2::DMRGSCF_debugPrint ){
            for ( int cnt = 0; cnt < linsize/2; cnt++ ){
               work2[ 2 * cnt     + linsize * ( 2 * cnt + 1 ) ] -= work1[ cnt ];
               work2[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] += work1[ cnt ];
            }
            double rms_diff = 0.0;
            for ( int cnt = 0; cnt < size; cnt++ ){ rms_diff += work2[ cnt ] * work2[ cnt ]; }
            rms_diff = sqrt( rms_diff );
            cout << "   DMRGSCFunitary::updateUnitary : RMSdeviation of irrep block " << irrep << " (should be 0.0) = " << rms_diff << endl;
         }

         //work2 <= exp(Bmat^T * xblock * Bmat)
         for ( int cnt = 0; cnt < size; cnt++ ){ work2[ cnt ] = 0.0; }
         for ( int cnt = 0; cnt < linsize/2; cnt++ ){
            double cosine = cos( work1[ cnt ] );
            double sine   = sin( work1[ cnt ] );
            work2[ 2 * cnt     + linsize * ( 2 * cnt    ) ] = cosine;
            work2[ 2 * cnt + 1 + linsize * ( 2 * cnt + 1) ] = cosine;
            work2[ 2 * cnt     + linsize * ( 2 * cnt + 1) ] = sine;
            work2[ 2 * cnt + 1 + linsize * ( 2 * cnt    ) ] = -sine;
         }
         for ( int cnt = 2*(linsize/2); cnt < linsize; cnt++ ){ work2[ cnt * (linsize + 1) ] = 1.0; }

         //work2 <= Bmat * exp(Bmat^T * xblock * Bmat) * Bmat^T = exp(xblock)
         dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &one, Bmat, &linsize, work2, &linsize, &set, work1, &linsize );
         dgemm_( &notrans, &trans, &linsize, &linsize, &linsize, &one, work1, &linsize, Bmat, &linsize, &set, work2, &linsize );

         int inc1 = 1;
         if ( multiply ){ //U <-- exp(x) * U
            dgemm_( &notrans, &notrans, &linsize, &linsize, &linsize, &one, work2, &linsize, entries[ irrep ], &linsize, &set, work1, &linsize );
            dcopy_( &size, work1, &inc1, entries[ irrep ], &inc1 );
         } else { //U <-- exp(x)
            dcopy_( &size, work2, &inc1, entries[ irrep ], &inc1 );
         }

      }
   }

   if (CheMPS2::DMRGSCF_debugPrint){ CheckDeviationFromUnitary( temp2 ); }

}

void CheMPS2::DMRGSCFunitary::rotateActiveSpaceVectors( double * eigenvecs, double * work ){

   int passed = 0;
   int tot_dmrg = iHandler->getDMRGcumulative( num_irreps );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int NDMRG = iHandler->getNDMRG( irrep );
      int NOCC  = iHandler->getNOCC( irrep );
      int NORB  = iHandler->getNORB( irrep );

      if ( NDMRG > 1 ){

         double * rotation = eigenvecs + passed * ( tot_dmrg + 1 );

         char trans = 'T';
         char notrans = 'N';
         double one = 1.0;
         double set = 0.0;
         dgemm_( &trans, &notrans, &NDMRG, &NORB, &NDMRG, &one, rotation, &tot_dmrg, entries[ irrep ] + NOCC, &NORB, &set, work, &NDMRG );

         for ( int row = 0; row < NDMRG; row++ ){
            for ( int col = 0; col < NORB; col++ ){
               entries[ irrep ][ NOCC + row + NORB * col ] = work[ row + NDMRG * col ];
            }
         }

      }
      passed += NDMRG;
   }

   if (CheMPS2::DMRGSCF_debugPrint){ CheckDeviationFromUnitary( work ); }

}

void CheMPS2::DMRGSCFunitary::CheckDeviationFromUnitary( double * work ) const{

   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int linsize = iHandler->getNORB( irrep );
      if ( linsize > 0 ){
         dgemm_( &trans, &notrans, &linsize, &linsize, &linsize, &one, entries[ irrep ], &linsize, entries[ irrep ], &linsize, &set, work, &linsize );
         for ( int diag = 0; diag < linsize; diag++ ){
            work[ diag * ( 1 + linsize ) ] -= 1.0;
         }
         double rms_diff = 0.0;
         for ( int cnt = 0; cnt < linsize * linsize; cnt++ ){
            rms_diff += work[ cnt ] * work[ cnt ];
         }
         rms_diff = sqrt( rms_diff );
         cout << "   DMRGSCFunitary::CheckDeviationFromUnitary : 2-norm of U[" << irrep << "]^T * U[" << irrep << "] - I = " << rms_diff << endl;
      }
   }

}

double CheMPS2::DMRGSCFunitary::get_determinant( const int irrep, double * work1, double * work2, double * work_eig, int lwork_eig ) const{

   int NORB = iHandler->getNORB( irrep );

   // S = U + U^T --> work1
   for ( int row = 0; row < NORB; row++ ){
      for ( int col = 0; col < NORB; col++ ){
         work1[ row + NORB * col ] = ( entries[ irrep ][ row + NORB * col ]
                                     + entries[ irrep ][ col + NORB * row ] );
      }
   }

   // Get the eigenvectors of S = U + U^T: eigvals in work2, eigvecs in work1
   char jobz = 'V';
   char uplo = 'U';
   int info;
   dsyev_( &jobz, &uplo, &NORB, work1, &NORB, work2, work_eig, &lwork_eig, &info );

   // Transform the original orthogonal matrix with the real orthonormal eigenbasis of S --> the result V^T U V is stored in work2
   char trans = 'T';
   char notrans = 'N';
   double one = 1.0;
   double set = 0.0;
   dgemm_( &trans,   &notrans, &NORB, &NORB, &NORB, &one, work1,    &NORB, entries[ irrep ], &NORB, &set, work_eig, &NORB );
   dgemm_( &notrans, &notrans, &NORB, &NORB, &NORB, &one, work_eig, &NORB, work1,            &NORB, &set, work2,    &NORB );

   /*  Work2 should contain
          * 1x1 blocks containing [-1]
          * 2x2 blocks containing [[cos(u) sin(u)][-sin(u) cos(u)]]
          * 1x1 blocks containing [+1]
       and is hence TRIDIAGONAL. Its determinant is easy to compute:
   */
   double f_low  = 1.0;
   double f_high = work2[ 0 ];
   for ( int diag = 1; diag < NORB; diag++ ){
      double temp = work2[ diag * ( 1 + NORB ) ] * f_high - work2[ diag + NORB * ( diag - 1 ) ] * work2[ diag - 1 + NORB * diag ] * f_low;
      f_low = f_high;
      f_high = temp;
   }
   const double determinant = f_high;

   if ( CheMPS2::DMRGSCF_debugPrint ){
      cout << "   DMRGSCFunitary::get_determinant : det( U[" << irrep << "] ) = " << determinant << endl;
   }

   return determinant;

}

void CheMPS2::DMRGSCFunitary::getLog( double * vector, double * temp1, double * temp2 ) const{

   int jump = 0;

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int linsize = iHandler->getNORB( irrep );
      int size = linsize * linsize;

      if ( linsize > 1 ){
         /* linsize >= 2; hence temp1 is at least of size 4*linsize*linsize
            if linsize <= 1; there corresponds no block in xmatrix to it */

         double * work1 = temp1;            //linsize * linsize
         double * work2 = temp1 +     size; //linsize * linsize
         double * work3 = temp1 + 2 * size; //linsize * linsize
         double * work4 = temp1 + 3 * size; //linsize * linsize

         double * workLARGE = temp2;
         int lworkLARGE     = 4 * size;

         // work1 contains eigvec( U + U^T )
         // work3 contains eigvec( U + U^T )^T * U * eigvec( U + U^T ) ( TRIDIAGONAL matrix )
         const double determinant = get_determinant( irrep, work1, work3, workLARGE, lworkLARGE );
         assert( determinant > 0.0 );

         //Fill work2 with ln(V^T U V) = ln(work3).
         for ( int cnt = 0; cnt < size; cnt++ ){ work2[ cnt ] = 0.0; }
         for ( int cnt = 0; cnt < linsize/2; cnt++ ){ //2x2 blocks
            const double cosinus = 0.5 * ( work3[ 2 * cnt + linsize * ( 2 * cnt     ) ] + work3[ 2 * cnt + 1 + linsize * ( 2 * cnt + 1 ) ] );
            const double sinus   = 0.5 * ( work3[ 2 * cnt + linsize * ( 2 * cnt + 1 ) ] - work3[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] );
            const double theta   = atan2( sinus, cosinus );
            if (CheMPS2::DMRGSCF_debugPrint){
               cout << "   DMRGSCFunitary::getLog : block = [ " << work3[ 2 * cnt     + linsize * ( 2 * cnt     ) ] << " , "
                                                                << work3[ 2 * cnt     + linsize * ( 2 * cnt + 1 ) ] << " ]" << endl;
               cout << "                                    [ " << work3[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] << " , "
                                                                << work3[ 2 * cnt + 1 + linsize * ( 2 * cnt + 1 ) ] << " ] and corresponds to theta = " << theta << endl;
            }
            work3[ 2 * cnt     + linsize * ( 2 * cnt     ) ] -= cosinus;
            work3[ 2 * cnt + 1 + linsize * ( 2 * cnt + 1 ) ] -= cosinus;
            work3[ 2 * cnt     + linsize * ( 2 * cnt + 1 ) ] -= sinus;
            work3[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] += sinus;
            work2[ 2 * cnt     + linsize * ( 2 * cnt + 1 ) ] = theta;
            work2[ 2 * cnt + 1 + linsize * ( 2 * cnt     ) ] = -theta;
         } //The rest are 1x1 blocks, corresponding to eigenvalue 1 --> ln(1) = 0 --> work2 does not need to be updated.
         for ( int cnt = 2*(linsize/2); cnt < linsize; cnt++ ){ work3[ cnt * ( linsize + 1 ) ] -= 1; }

         //Calculate the 2-norm of updated work3 (should be 0.0)
         if (CheMPS2::DMRGSCF_debugPrint){
            double RMSdeviation = 0.0;
            for ( int cnt = 0; cnt < size; cnt++ ){ RMSdeviation += work3[ cnt ] * work3[ cnt ]; }
            RMSdeviation = sqrt( RMSdeviation );
            cout << "   DMRGSCFunitary::getLog : 2-norm of [ V^T*U*V - assumed blocks ] for irrep " << irrep << " (should be 0.0) = " << RMSdeviation << endl;
         }

         //Calculate V * ln(blocks) * V^T --> work4
         char trans = 'T';
         char notrans = 'N';
         double one = 1.0;
         double set = 0.0;
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &one, work1,     &linsize, work2, &linsize, &set, workLARGE, &linsize ); //V*ln(blocks)
         dgemm_(&notrans, &trans,   &linsize, &linsize, &linsize, &one, workLARGE, &linsize, work1, &linsize, &set, work4,     &linsize ); //V*ln(blocks)*V^T

         //Fill the vector with the just calculated ln(U) = work4
         for ( int row = 0; row < linsize; row++ ){
            for ( int col = row + 1; col < linsize; col++ ){
               vector[ jump + row + ( col * ( col - 1 ) ) / 2 ] = 0.5 * ( work4[ row + linsize * col ] - work4[ col + linsize * row ] );
            }
         }

         jump += ( linsize * ( linsize - 1 ) ) / 2;

      }
   }

   if ( true ){
      DMRGSCFunitary tempU = DMRGSCFunitary( iHandler );
      tempU.updateUnitary( temp1, temp2, vector, false, false ); //multiply = compact = false
      const double rms_diff = rms_deviation( &tempU );
      cout << "   DMRGSCFunitary::getLog : 2-norm of [ U - exp(ln(U)) ] (should be 0.0) = " << rms_diff << endl;
   }

}

void CheMPS2::DMRGSCFunitary::makeSureAllBlocksDetOne( double * temp1, double * temp2 ){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int linsize = iHandler->getNORB( irrep );
      int size = linsize * linsize;

      if ( linsize > 1 ){
         /* linsize >= 2; hence temp1 is at least of size 4*linsize*linsize
            if linsize <= 1; there corresponds no block in xmatrix to it */

         double * work1 = temp1;        //linsize * linsize
         double * work2 = temp1 + size; //linsize * linsize

         double * workLARGE = temp2;
         int lworkLARGE     = 4 * size;

         const double determinant = get_determinant( irrep, work1, work2, workLARGE, lworkLARGE );
         if ( determinant < 0.0 ){ // determinant = -1
            // U <-- diag(-1, 1, 1, 1, ...) * U : First row of U changes sign!
            for ( int cnt = 0; cnt < linsize; cnt++ ){ entries[ irrep ][ 0 + linsize * cnt ] *= -1; }
         }

      }
   }

}

void CheMPS2::DMRGSCFunitary::saveU( const string filename ) const{

   CheMPS2::DMRGSCFmatrix::write( filename, iHandler, entries );

}

void CheMPS2::DMRGSCFunitary::loadU( const string filename ){

   CheMPS2::DMRGSCFmatrix::read( filename, num_irreps, entries );

}


