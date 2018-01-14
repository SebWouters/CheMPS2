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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include "Davidson.h"
#include "Lapack.h"

using std::cout;
using std::endl;

CheMPS2::Davidson::Davidson( const int veclength, const int MAX_NUM_VEC, const int NUM_VEC_KEEP, const double RTOL, const double DIAG_CUTOFF, const bool debug_print, const char problem_type ){

   assert( ( problem_type == 'E' ) || ( problem_type == 'L' ) );

   this->debug_print  = debug_print;
   this->veclength    = veclength;
   this->problem_type = problem_type;
   this->MAX_NUM_VEC  = MAX_NUM_VEC;
   this->NUM_VEC_KEEP = NUM_VEC_KEEP;
   this->DIAG_CUTOFF  = DIAG_CUTOFF;
   this->RTOL         = RTOL;

   state = 'I'; // <I>nitialized Davidson
   nMultiplications = 0;

   // To store the vectors and the matrix x vectors
   num_vec = 0;
   vecs  = new double*[ MAX_NUM_VEC ];
   Hvecs = new double*[ MAX_NUM_VEC ];
   num_allocated = 0;

   // The projected problem
   mxM       = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];
   mxM_eigs  = new double[ MAX_NUM_VEC ];
   mxM_vecs  = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];
   mxM_lwork = 3 * MAX_NUM_VEC - 1;
   mxM_work  = new double[ mxM_lwork ];
   mxM_rhs   = (( problem_type == 'L' ) ? new double[ MAX_NUM_VEC ] : NULL );

   // Vector spaces
   diag     = new double[ veclength ];
   t_vec    = new double[ veclength ];
   u_vec    = new double[ veclength ];
   work_vec = new double[ veclength ];
   RHS      = (( problem_type == 'L' ) ? new double[ veclength ] : NULL );

   // For the deflation
   Reortho_Lowdin       = NULL;
   Reortho_Overlap_eigs = NULL;
   Reortho_Overlap      = NULL;
   Reortho_Eigenvecs    = NULL;

}

CheMPS2::Davidson::~Davidson(){

   for (int cnt = 0; cnt < num_allocated; cnt++){
      delete [] vecs[cnt];
      delete [] Hvecs[cnt];
   }
   delete [] vecs;
   delete [] Hvecs;

   delete [] mxM;
   delete [] mxM_eigs;
   delete [] mxM_vecs;
   delete [] mxM_work;
   if ( mxM_rhs != NULL ){ delete [] mxM_rhs; }

   delete [] diag;
   delete [] t_vec;
   delete [] u_vec;
   delete [] work_vec;
   if ( RHS != NULL ){ delete [] RHS; }

   if ( Reortho_Lowdin       != NULL ){ delete [] Reortho_Lowdin; }
   if ( Reortho_Overlap_eigs != NULL ){ delete [] Reortho_Overlap_eigs; }
   if ( Reortho_Overlap      != NULL ){ delete [] Reortho_Overlap; }
   if ( Reortho_Eigenvecs    != NULL ){ delete [] Reortho_Eigenvecs; }

}

int CheMPS2::Davidson::GetNumMultiplications() const{ return nMultiplications; }

char CheMPS2::Davidson::FetchInstruction( double ** pointers ){

   /* 
      Possible states:
       - I : just initialized
       - U : just before the big loop, the initial guess and the diagonal are set
       - N : a new vector has just been added to the list and a matrix-vector multiplication has been performed
       - F : the space has been deflated and a few matrix-vector multiplications are required
       - C : convergence was reached

      Possible instructions:
       - A : copy the initial guess to pointers[0], the diagonal to pointers[1], and if problem_type=='L' the right-hand side of the linear problem to pointers[2]
       - B : perform pointers[1] = symmetric matrix times pointers[0]
       - C : copy the converged solution from pointers[0] back; pointers[1][0] contains the converged energy if problem_type=='E' and the residual norm if problem_type=='L'
       - D : there was an error
   */

   if ( state == 'I' ){
      pointers[ 0 ] = t_vec;
      pointers[ 1 ] = diag;
      if ( problem_type == 'L' ){ pointers[ 2 ] = RHS; }
      state = 'U';
      return 'A';
   }

   if ( state == 'U' ){
      SafetyCheckGuess();
      AddNewVec();
      pointers[ 0 ] =  vecs[ num_vec ];
      pointers[ 1 ] = Hvecs[ num_vec ];
      nMultiplications++;
      state = 'N';
      return 'B';
   }

   if ( state == 'N' ){
      const double rnorm = DiagonalizeSmallMatrixAndCalcResidual();
      // if ( debug_print ){ cout << "WARNING AT DAVIDSON : Current residual norm = " << rnorm << endl; }
      if ( rnorm > RTOL ){ // Not yet converged
         CalculateNewVec();
         if ( num_vec == MAX_NUM_VEC ){
            Deflation();
            pointers[ 0 ] =  vecs[ num_vec ];
            pointers[ 1 ] = Hvecs[ num_vec ];
            nMultiplications++;
            num_vec++;
            state = 'F';
            return 'B';
         }
         AddNewVec();
         pointers[ 0 ] =  vecs[ num_vec ];
         pointers[ 1 ] = Hvecs[ num_vec ];
         nMultiplications++;
         state = 'N';
         return 'B';
      } else { // Converged
         state = 'C';
         pointers[ 0 ] = u_vec;
         pointers[ 1 ] = work_vec;
         if ( problem_type == 'E' ){ work_vec[ 0 ] = mxM_eigs[ 0 ]; }
         if ( problem_type == 'L' ){ work_vec[ 0 ] = rnorm; }
         return 'C';
      }
   }

   if ( state == 'F' ){
      if ( num_vec == NUM_VEC_KEEP ){
         MxMafterDeflation();
         AddNewVec();
         pointers[ 0 ] =  vecs[ num_vec ];
         pointers[ 1 ] = Hvecs[ num_vec ];
         nMultiplications++;
         state = 'N';
         return 'B';
      } else {
         pointers[ 0 ] =  vecs[ num_vec ];
         pointers[ 1 ] = Hvecs[ num_vec ];
         nMultiplications++;
         num_vec++;
         state = 'F';
         return 'B';
      }
   }

   return 'D';

}

double CheMPS2::Davidson::FrobeniusNorm( double * current_vector ){

   char frobenius = 'F';
   int inc1 = 1;
   const double twonorm = dlange_( &frobenius, &veclength, &inc1, current_vector, &veclength, NULL ); // Work is not referenced for Frobenius norm
   return twonorm;

}

void CheMPS2::Davidson::SafetyCheckGuess(){

   const double twonorm = FrobeniusNorm( t_vec );
   if ( twonorm == 0.0 ){
      for ( int cnt = 0; cnt < veclength; cnt++ ){ t_vec[ cnt ] = ( (double) rand() ) / RAND_MAX; }
      if ( debug_print ){
         cout << "WARNING AT DAVIDSON : Initial guess was a zero-vector. Now it is overwritten with random numbers." << endl;
      }
   }

}

void CheMPS2::Davidson::AddNewVec(){

   int inc1 = 1;

   // Orthogonalize the new vector w.r.t. the old basis
   for ( int cnt = 0; cnt < num_vec; cnt++ ){
      double minus_overlap = - ddot_( &veclength, t_vec, &inc1, vecs[ cnt ], &inc1 );
      daxpy_( &veclength, &minus_overlap, vecs[ cnt ], &inc1, t_vec, &inc1 );
   }

   // Normalize the new vector
   double alpha = 1.0 / FrobeniusNorm( t_vec );
   dscal_( &veclength, &alpha, t_vec, &inc1 );

   // The new vector becomes part of vecs
   if ( num_vec < num_allocated ){
      double * temp = vecs[ num_vec ];
      vecs[ num_vec ] = t_vec;
      t_vec = temp;
   } else {
      vecs[ num_allocated ] = t_vec;
      Hvecs[ num_allocated ] = new double[ veclength ];
      t_vec = new double[ veclength ];
      num_allocated++;
   }

}

double CheMPS2::Davidson::DiagonalizeSmallMatrixAndCalcResidual(){

   int inc1 = 1;

   if ( problem_type == 'E' ){ // EIGENVALUE PROBLEM
      // mxM contains V^T . A . V
      for ( int cnt = 0; cnt < num_vec; cnt++ ){
         mxM[ cnt + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, vecs[ num_vec ], &inc1, Hvecs[ cnt ], &inc1 );
         mxM[ num_vec + MAX_NUM_VEC * cnt ] = mxM[ cnt + MAX_NUM_VEC * num_vec ];
      }
      mxM[ num_vec + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, vecs[ num_vec ], &inc1, Hvecs[ num_vec ], &inc1 );
   } else { // LINEAR PROBLEM
      // mxM contains V^T . A^T . A . V
      for ( int cnt = 0; cnt < num_vec; cnt++ ){
         mxM[ cnt + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, Hvecs[ num_vec ], &inc1, Hvecs[ cnt ], &inc1 );
         mxM[ num_vec + MAX_NUM_VEC * cnt ] = mxM[ cnt + MAX_NUM_VEC * num_vec ];
      }
      mxM[ num_vec + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, Hvecs[ num_vec ], &inc1, Hvecs[ num_vec ], &inc1 );
      // mxM_rhs contains V^T . A^T . RHS
      mxM_rhs[ num_vec ] = ddot_( &veclength, Hvecs[ num_vec ], &inc1, RHS, &inc1 );
   }

   // When t-vec was added to vecs, the number of vecs was actually increased by one. Now the number is incremented.
   num_vec++;

   // Diagonalize mxM ( always )
   char jobz = 'V';
   char uplo = 'U';
   int info;
   for ( int cnt1 = 0; cnt1 < num_vec; cnt1++ ){
      for ( int cnt2 = 0; cnt2 < num_vec; cnt2++ ){
         mxM_vecs[ cnt1 + MAX_NUM_VEC * cnt2 ] = mxM[ cnt1 + MAX_NUM_VEC * cnt2 ];
      }
   }
   dsyev_( &jobz, &uplo, &num_vec, mxM_vecs, &MAX_NUM_VEC, mxM_eigs, mxM_work, &mxM_lwork, &info ); // Ascending order of eigenvalues

   // If problem_type == linear, solve the linear problem in least-squares sense!
   if ( problem_type == 'L' ){
      double one = 1.0;
      double set = 0.0;
      char trans = 'T';
      char notra = 'N';
      dgemm_( &trans, &notra, &num_vec, &inc1, &num_vec, &one, mxM_vecs, &MAX_NUM_VEC, mxM_rhs, &MAX_NUM_VEC, &set, mxM_work, &MAX_NUM_VEC ); // mxM_work = mxM_vecs^T * mxM_rhs
      for ( int cnt = 0; cnt < num_vec; cnt++ ){
         double current_eigenvalue = mxM_eigs[ cnt ];
         if ( fabs( current_eigenvalue ) < DIAG_CUTOFF ){
            current_eigenvalue = DIAG_CUTOFF * (( current_eigenvalue < 0.0 ) ? -1 : 1 );
            if ( debug_print ){
               cout << "WARNING AT DAVIDSON : The eigenvalue " << mxM_eigs[ cnt ] << " to solve Ax = b has been overwritten with " << current_eigenvalue << "." << endl;
            }
         }
         mxM_work[ cnt ] = mxM_work[ cnt ] / current_eigenvalue;
      }
      dgemm_( &notra, &notra, &num_vec, &inc1, &num_vec, &one, mxM_vecs, &MAX_NUM_VEC, mxM_work, &MAX_NUM_VEC, &set, mxM_work + MAX_NUM_VEC, &MAX_NUM_VEC );
      for ( int cnt = 0; cnt < num_vec; cnt++ ){ mxM_vecs[ cnt ] = mxM_work[ MAX_NUM_VEC + cnt ]; } // mxM_vecs = U * eigs^{-1} * U^T * RHS
   }

   // Calculate u and r. r is stored in t_vec, u in u_vec.
   for ( int cnt = 0; cnt < veclength; cnt++ ){ t_vec[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < veclength; cnt++ ){ u_vec[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < num_vec; cnt++ ){
      double alpha = mxM_vecs[ cnt ]; // Eigenvector with lowest eigenvalue, hence mxM_vecs[ cnt + MAX_NUM_VEC * 0 ]
      daxpy_( &veclength, &alpha, Hvecs[ cnt ], &inc1, t_vec, &inc1 );
      daxpy_( &veclength, &alpha,  vecs[ cnt ], &inc1, u_vec, &inc1 );
   }
   if ( problem_type == 'E' ){ // t_vec = H * x - lambda * x
      double alpha = - mxM_eigs[ 0 ];
      daxpy_( &veclength, &alpha, u_vec, &inc1, t_vec, &inc1 );
   } else { // t_vec = H * x - RHS
      double alpha = -1.0;
      daxpy_( &veclength, &alpha, RHS, &inc1, t_vec, &inc1 );
   }

   // Calculate the norm of r
   const double rnorm = FrobeniusNorm( t_vec );
   //cout << "    Davidson :: rnorm = " << rnorm << endl;
   return rnorm;

}

void CheMPS2::Davidson::CalculateNewVec(){

   int inc1 = 1;
   const double shift = (( problem_type == 'E' ) ? mxM_eigs[ 0 ] : 0.0 );

   // Calculate the new t_vec based on the residual of the lowest eigenvalue, to add to the vecs.
   for ( int cnt = 0; cnt < veclength; cnt++ ){
      const double difference = diag[ cnt ] - shift;
      const double fabsdiff   = fabs( difference );
      if ( fabsdiff > DIAG_CUTOFF ){
         work_vec[ cnt ] = u_vec[ cnt ] / difference; // work_vec = K^(-1) u_vec
      } else {
         work_vec[ cnt ] = u_vec[ cnt ] / DIAG_CUTOFF;
         if ( debug_print ){ cout << "WARNING AT DAVIDSON : fabs( precon[" << cnt << "] ) = " << fabsdiff << endl; }
      }
   }
   double alpha = - ddot_( &veclength, work_vec, &inc1, t_vec, &inc1 ) / ddot_( &veclength, work_vec, &inc1, u_vec, &inc1 ); // alpha = - (u^T K^(-1) r) / (u^T K^(-1) u)
   daxpy_( &veclength, &alpha, u_vec, &inc1, t_vec, &inc1 ); // t_vec = r - (u^T K^(-1) r) / (u^T K^(-1) u) u
   for ( int cnt = 0; cnt < veclength; cnt++ ){
      const double difference = diag[ cnt ] - shift;
      const double fabsdiff   = fabs( difference );
       if ( fabsdiff > DIAG_CUTOFF ){
         t_vec[ cnt ] = - t_vec[ cnt ] / difference; // t_vec = - K^(-1) (r - (u^T K^(-1) r) / (u^T K^(-1) u) u)
      } else {
         t_vec[ cnt ] = - t_vec[ cnt ] / DIAG_CUTOFF;
      }
   }

}

void CheMPS2::Davidson::Deflation(){

   int inc1 = 1;

   // When the maximum number of vectors is reached: construct the one with lowest eigenvalue & restart
   if ( NUM_VEC_KEEP <= 1 ){

      double alpha = 1.0 / FrobeniusNorm( u_vec );
      dscal_( &veclength, &alpha, u_vec, &inc1 );
      dcopy_( &veclength, u_vec, &inc1, vecs[ 0 ], &inc1 );

   } else {

      if ( problem_type == 'L' ){ SolveLinearSystemDeflation( NUM_VEC_KEEP ); }

      if ( Reortho_Eigenvecs    == NULL ){ Reortho_Eigenvecs    = new double[    veclength * NUM_VEC_KEEP ]; }
      if ( Reortho_Overlap      == NULL ){ Reortho_Overlap      = new double[ NUM_VEC_KEEP * NUM_VEC_KEEP ]; }
      if ( Reortho_Overlap_eigs == NULL ){ Reortho_Overlap_eigs = new double[ NUM_VEC_KEEP                ]; }
      if ( Reortho_Lowdin       == NULL ){ Reortho_Lowdin       = new double[ NUM_VEC_KEEP * NUM_VEC_KEEP ]; }
   
      // Construct the lowest NUM_VEC_KEEP eigenvectors
      dcopy_( &veclength, u_vec, &inc1, Reortho_Eigenvecs, &inc1 );
      for ( int cnt = 1; cnt < NUM_VEC_KEEP; cnt++ ){
         for ( int irow = 0; irow < veclength; irow++ ){
            Reortho_Eigenvecs[ irow + veclength * cnt ] = 0.0;
            for ( int ivec = 0; ivec < MAX_NUM_VEC; ivec++ ){
               Reortho_Eigenvecs[ irow + veclength * cnt ] += vecs[ ivec ][ irow ] * mxM_vecs[ ivec + MAX_NUM_VEC * cnt ];
            }
         }
      }

      // Calculate the overlap matrix
      char trans  = 'T';
      char notr   = 'N';
      double one  = 1.0;
      double zero = 0.0; //set
      dgemm_( &trans, &notr, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &veclength, &one, Reortho_Eigenvecs, &veclength, Reortho_Eigenvecs, &veclength, &zero, Reortho_Overlap, &NUM_VEC_KEEP );

      // Calculate the Lowdin tfo
      char jobz = 'V';
      char uplo = 'U';
      int info;
      dsyev_( &jobz, &uplo, &NUM_VEC_KEEP, Reortho_Overlap, &NUM_VEC_KEEP, Reortho_Overlap_eigs, mxM_work, &mxM_lwork, &info ); // Ascending order of eigs
      for ( int icnt = 0; icnt < NUM_VEC_KEEP; icnt++ ){
         Reortho_Overlap_eigs[ icnt ] = pow( Reortho_Overlap_eigs[ icnt ], -0.25 );
         dscal_( &NUM_VEC_KEEP, Reortho_Overlap_eigs + icnt, Reortho_Overlap + NUM_VEC_KEEP * icnt, &inc1 );
      }
      dgemm_( &notr, &trans, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &one, Reortho_Overlap, &NUM_VEC_KEEP, Reortho_Overlap, &NUM_VEC_KEEP, &zero, Reortho_Lowdin, &NUM_VEC_KEEP );

      // Reortho: Put the Lowdin tfo eigenvecs in vecs
      for ( int ivec = 0; ivec < NUM_VEC_KEEP; ivec++ ){
         for ( int loop = 0; loop < veclength; loop++ ){ vecs[ ivec ][ loop ] = 0.0; }
         for ( int ivec2 = 0; ivec2 < NUM_VEC_KEEP; ivec2++ ){
            daxpy_( &veclength, Reortho_Lowdin + ivec2 + NUM_VEC_KEEP * ivec, Reortho_Eigenvecs + veclength * ivec2, &inc1, vecs[ ivec ], &inc1 );
         }
      }
   }

   num_vec = 0;

}

void CheMPS2::Davidson::SolveLinearSystemDeflation( const int NUM_SOLUTIONS ){

   assert( problem_type == 'L' );
   assert( num_vec == MAX_NUM_VEC );
   assert( NUM_SOLUTIONS <= MAX_NUM_VEC );
   assert( 2 <= NUM_SOLUTIONS );

   double * work1 = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];   // projector
   double * work3 = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];   // projector times mxM_rhs
   double * work2 = new double[ MAX_NUM_VEC * NUM_SOLUTIONS ]; // solutions

   for ( int solution = 0; solution < NUM_SOLUTIONS; solution++ ){

      // work1 = ( 1 - sum_j v_j v_j^T ) = projector
      for ( int cntr = 0; cntr < MAX_NUM_VEC * MAX_NUM_VEC; cntr++ ){ work1[ cntr ] = 0.0; }
      for ( int diag = 0; diag < MAX_NUM_VEC; diag++ ){ work1[ diag * ( 1 + MAX_NUM_VEC ) ] = 1.0; }
      for ( int prev = 0; prev < solution; prev++ ){
         for ( int col = 0; col < MAX_NUM_VEC; col++ ){
            for ( int row = 0; row < MAX_NUM_VEC; row++ ){
               work1[ row + MAX_NUM_VEC * col ] -= work2[ row + MAX_NUM_VEC * prev ] * work2[ col + MAX_NUM_VEC * prev ];
            }
         }
      }

      // work3    = ( 1 - sum_j v_j v_j^T ) * [ U^T * A^T * b ]                               = work1 * mxM_rhs
      // mxM_vecs = ( 1 - sum_j v_j v_j^T ) * [ U^T * A^T * A * U ] * ( 1 - sum_j v_j v_j^T ) = work1 * mxM * work1
      {
         double one = 1.0;
         double set = 0.0;
         char notrans = 'N';
         int inc1 = 1;
         dgemm_( &notrans, &notrans, &MAX_NUM_VEC, &MAX_NUM_VEC, &MAX_NUM_VEC, &one, work1, &MAX_NUM_VEC, mxM,     &MAX_NUM_VEC, &set, work3,    &MAX_NUM_VEC );
         dgemm_( &notrans, &notrans, &MAX_NUM_VEC, &MAX_NUM_VEC, &MAX_NUM_VEC, &one, work3, &MAX_NUM_VEC, work1,   &MAX_NUM_VEC, &set, mxM_vecs, &MAX_NUM_VEC );
         dgemm_( &notrans, &notrans, &MAX_NUM_VEC, &inc1,        &MAX_NUM_VEC, &one, work1, &MAX_NUM_VEC, mxM_rhs, &MAX_NUM_VEC, &set, work3,    &MAX_NUM_VEC );
      }

      // Diagonalize mxM_vecs = V * lambda * V^T   ===>   Ascending order of eigenvalues   &    ( V, lambda ) = ( mxM_vecs, mxM_eigs )
      {
         char jobz = 'V';
         char uplo = 'U';
         int info;
         dsyev_( &jobz, &uplo, &MAX_NUM_VEC, mxM_vecs, &MAX_NUM_VEC, mxM_eigs, mxM_work, &mxM_lwork, &info );
      }

      // work2[ :, solution ] = [ ( 1 - sum_j v_j v_j^T ) * [ U^T * A^T * A * U ] * ( 1 - sum_j v_j v_j^T ) ]^{-1} * ( 1 - sum_j v_j v_j^T ) * [ U^T * A^T * b ] = V * lambda^{-1} * V^T * work3
      {
         double one = 1.0;
         double set = 0.0;
         char trans = 'T';
         char notrans = 'N';
         int inc1 = 1;
         dgemm_( &trans, &notrans, &MAX_NUM_VEC, &inc1, &MAX_NUM_VEC, &one, mxM_vecs, &MAX_NUM_VEC, work3, &MAX_NUM_VEC, &set, mxM_work, &MAX_NUM_VEC );
         for ( int diag = 0; diag < MAX_NUM_VEC; diag++ ){
            if ( diag < solution ){
               mxM_work[ diag ] = 0.0; // PSEUDOINVERSE
            } else {
               double current_eigenvalue = mxM_eigs[ diag ];
               if ( fabs( current_eigenvalue ) < DIAG_CUTOFF ){
                  current_eigenvalue = DIAG_CUTOFF * (( current_eigenvalue < 0.0 ) ? -1 : 1 );
                  if ( debug_print ){
                     cout << "WARNING AT DAVIDSON : The eigenvalue " << mxM_eigs[ diag ] << " to solve Ax = b has been overwritten with " << current_eigenvalue << "." << endl;
                  }
               }
               mxM_work[ diag ] = mxM_work[ diag ] / current_eigenvalue;
            }
         }
         dgemm_( &notrans, &notrans, &MAX_NUM_VEC, &inc1, &MAX_NUM_VEC, &one, mxM_vecs, &MAX_NUM_VEC, mxM_work, &MAX_NUM_VEC, &set, work2 + MAX_NUM_VEC * solution, &MAX_NUM_VEC );
      }

      // Normalize work2[ :, solution ]
      {
         int inc1 = 1;
         double * ptr = work2 + MAX_NUM_VEC * solution;
         const double twonorm = sqrt( ddot_( &MAX_NUM_VEC, ptr, &inc1, ptr, &inc1 ) );
         /*if ( debug_print ){
            cout << "Davidson :: Deflation :: Norm of solution " << solution << " = " << twonorm << endl;
         }*/
         double factor = 1.0 / twonorm;
         dscal_( &MAX_NUM_VEC, &factor, ptr, &inc1 );
      }

   }

   // Copy over work2 to mxM_vecs
   {
      int inc1 = 1;
      int size = MAX_NUM_VEC * NUM_SOLUTIONS;
      dcopy_( &size, work2, &inc1, mxM_vecs, &inc1 );
   }

   delete [] work1;
   delete [] work2;
   delete [] work3;

}

void CheMPS2::Davidson::MxMafterDeflation(){

   int inc1 = 1;

   if ( problem_type == 'E' ){ // EIGENVALUE PROBLEM
      // mxM contains V^T . A . V
      for ( int ivec = 0; ivec < NUM_VEC_KEEP; ivec++ ){
         for ( int ivec2 = ivec; ivec2 < NUM_VEC_KEEP; ivec2++ ){
            mxM[ ivec + MAX_NUM_VEC * ivec2 ] = ddot_( &veclength, vecs[ ivec ], &inc1, Hvecs[ ivec2 ], &inc1 );
            mxM[ ivec2 + MAX_NUM_VEC * ivec ] = mxM[ ivec + MAX_NUM_VEC * ivec2 ];
         }
      }
   } else { // LINEAR PROBLEM
      // mxM contains V^T . A^T . A . V
      for ( int ivec = 0; ivec < NUM_VEC_KEEP; ivec++ ){
         for ( int ivec2 = ivec; ivec2 < NUM_VEC_KEEP; ivec2++ ){
            mxM[ ivec + MAX_NUM_VEC * ivec2 ] = ddot_( &veclength, Hvecs[ ivec ], &inc1, Hvecs[ ivec2 ], &inc1 );
            mxM[ ivec2 + MAX_NUM_VEC * ivec ] = mxM[ ivec + MAX_NUM_VEC * ivec2 ];
         }
      }
      // mxM_rhs contains V^T . A^T . RHS
      for ( int ivec = 0; ivec < NUM_VEC_KEEP; ivec++ ){
         mxM_rhs[ ivec ] = ddot_( &veclength, Hvecs[ ivec ], &inc1, RHS, &inc1 );
      }
   }

}


