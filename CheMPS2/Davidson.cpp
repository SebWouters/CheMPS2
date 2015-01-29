/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2015 Sebastian Wouters

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

#include "Davidson.h"
#include "Lapack.h"

using std::cout;
using std::endl;

CheMPS2::Davidson::Davidson(const int veclength_in, const int MAX_NUM_VEC_in, const int NUM_VEC_KEEP_in, const double RTOL_in, const double DIAG_CUTOFF_in, const bool debugPrint_in){

   debugPrint = debugPrint_in;
   veclength = veclength_in;
   state = 'I'; // <I>nitialized Davidson
   nMultiplications = 0;
   
   MAX_NUM_VEC = MAX_NUM_VEC_in;
   NUM_VEC_KEEP = NUM_VEC_KEEP_in;
   DIAG_CUTOFF = DIAG_CUTOFF_in;
   RTOL = RTOL_in;
   
   // To store the vectors and the matrix x vectors
   num_vec = 0;
   vecs  = new double*[ MAX_NUM_VEC ];
   Hvecs = new double*[ MAX_NUM_VEC ];
   num_allocated = 0;
   
   // The effective diagonalization problem
   mxM       = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];
   mxM_eigs  = new double[ MAX_NUM_VEC ];
   mxM_vecs  = new double[ MAX_NUM_VEC * MAX_NUM_VEC ];
   mxM_lwork = 3 * MAX_NUM_VEC - 1;
   mxM_work  = new double[ mxM_lwork ];
   
   // Vector spaces
   diag     = new double[ veclength ];
   t_vec    = new double[ veclength ];
   u_vec    = new double[ veclength ];
   work_vec = new double[ veclength ];
   
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
   
   delete [] diag;
   delete [] t_vec;
   delete [] u_vec;
   delete [] work_vec;
   
   if ( Reortho_Lowdin       != NULL ){ delete [] Reortho_Lowdin; }
   if ( Reortho_Overlap_eigs != NULL ){ delete [] Reortho_Overlap_eigs; }
   if ( Reortho_Overlap      != NULL ){ delete [] Reortho_Overlap; }
   if ( Reortho_Eigenvecs    != NULL ){ delete [] Reortho_Eigenvecs; }

}

int CheMPS2::Davidson::GetNumMultiplications() const{ return nMultiplications; }

char CheMPS2::Davidson::FetchInstruction(double ** whichpointers){

   /* 
      Possible states:
       - I : just initialized
       - U : just before the big loop, the initial guess and the diagonal are set
       - N : a new vector has just been added to the list and a matrix-vector multiplication has been performed
       - F : the space has been deflated and a few matrix-vector multiplications are required
       - C : convergence was reached
   
      Possible instructions:
       - A : copy the initial guess to whichpointers[0] and the diagonal to whichpointers[1]
       - B : perform whichpointers[1] = symmetric matrix times whichpointers[0]
       - C : copy the converged solution from whichpointers[0] back; whichpointers[1][0] contains the converged energy
       - D : there was an error
   */

   if ( state == 'I' ){
      whichpointers[0] = t_vec;
      whichpointers[1] = diag;
      state = 'U';
      return 'A';
   }
   
   if ( state == 'U' ){
      SafetyCheckGuess();
      AddNewVec();
      whichpointers[0] =  vecs[ num_vec ];
      whichpointers[1] = Hvecs[ num_vec ];
      nMultiplications++;
      state = 'N';
      return 'B';
   }
   
   if ( state == 'N' ){
      const double rnorm = DiagonalizeSmallMatrixAndCalcResidual();
      if ( rnorm > RTOL ){ // Not yet converged
         CalculateNewVec();
         if ( num_vec == MAX_NUM_VEC ){
            Deflation();
            whichpointers[0] =  vecs[ num_vec ];
            whichpointers[1] = Hvecs[ num_vec ];
            nMultiplications++;
            num_vec++;
            state = 'F';
            return 'B';
         }
         AddNewVec();
         whichpointers[0] =  vecs[ num_vec ];
         whichpointers[1] = Hvecs[ num_vec ];
         nMultiplications++;
         state = 'N';
         return 'B';
      } else { // Converged
         state = 'C';
         whichpointers[0] = u_vec;
         whichpointers[1] = mxM_eigs;
         return 'C';
      }
   }
   
   if ( state == 'F' ){
      if ( num_vec == NUM_VEC_KEEP ){
         MxMafterDeflation();
         AddNewVec();
         whichpointers[0] =  vecs[ num_vec ];
         whichpointers[1] = Hvecs[ num_vec ];
         nMultiplications++;
         state = 'N';
         return 'B';
      } else {
         whichpointers[0] =  vecs[ num_vec ];
         whichpointers[1] = Hvecs[ num_vec ];
         nMultiplications++;
         num_vec++;
         state = 'F';
         return 'B';
      }
   }

   return 'D';

}

void CheMPS2::Davidson::SafetyCheckGuess(){

   char frobenius = 'F';
   int inc1 = 1;
   const double twonorm = dlange_( &frobenius, &veclength, &inc1, t_vec, &veclength, NULL ); // Work is not referenced for Frobenius norm
   if ( twonorm == 0.0 ){
      for (int cnt = 0; cnt < veclength; cnt++){ t_vec[ cnt ] = ((double) rand())/RAND_MAX; }
      if ( debugPrint ){
         cout << "WARNING AT DAVIDSON : Initial guess was a zero-vector. Now it is overwritten with random numbers." << endl;
      }
   }

}

void CheMPS2::Davidson::AddNewVec(){

   int inc1 = 1;
   
   //1. Orthogonalize the new vector w.r.t. the old basis
   for (int cnt = 0; cnt < num_vec; cnt++){
      double min_overlap = - ddot_( &veclength, t_vec, &inc1, vecs[ cnt ], &inc1 );
      daxpy_( &veclength, &min_overlap, vecs[ cnt ], &inc1, t_vec, &inc1 );
   }

   //2. Normalize the new vector
   char frobenius = 'F';
   double alpha = 1.0 / dlange_( &frobenius, &veclength, &inc1, t_vec, &veclength, NULL ); // Work is not referenced for Frobenius norm
   dscal_( &veclength, &alpha, t_vec, &inc1 );
   
   //3. The new vector becomes part of vecs
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

   //4. mxM contains the Hamiltonian in the basis "vecs"
   for (int cnt = 0; cnt < num_vec; cnt++){
      mxM[ cnt + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, vecs[ num_vec ], &inc1, Hvecs[ cnt ], &inc1 );
      mxM[ num_vec + MAX_NUM_VEC * cnt ] = mxM[ cnt + MAX_NUM_VEC * num_vec ];
   }
   mxM [ num_vec + MAX_NUM_VEC * num_vec ] = ddot_( &veclength, vecs[ num_vec ], &inc1, Hvecs[ num_vec ], &inc1 );
   
   //5. When t-vec was added to vecs, the number of vecs was actually increased by one. For convenience (doing 4.), only now the number is incremented.
   num_vec++;
   
   //6. Calculate the eigenvalues and vectors of mxM
   char jobz = 'V';
   char uplo = 'U';
   int info;
   for (int cnt1 = 0; cnt1 < num_vec; cnt1++){
      for (int cnt2 = 0; cnt2 < num_vec; cnt2++){
         mxM_vecs[ cnt1 + MAX_NUM_VEC * cnt2 ] = mxM[ cnt1 + MAX_NUM_VEC * cnt2 ];
      }
   }
   dsyev_( &jobz, &uplo, &num_vec, mxM_vecs, &MAX_NUM_VEC, mxM_eigs, mxM_work, &mxM_lwork, &info ); // Ascending order of eigenvalues
   
   //7. Calculate u and r. r is stored in t_vec, u in u_vec.
   for (int cnt = 0; cnt < veclength; cnt++){ t_vec[ cnt ] = 0.0; }
   for (int cnt = 0; cnt < veclength; cnt++){ u_vec[ cnt ] = 0.0; }
   for (int cnt = 0; cnt < num_vec; cnt++){
      double alpha = mxM_vecs[ cnt ]; // Eigenvector with lowest eigenvalue, hence mxM_vecs[ cnt + MAX_NUM_VEC * 0 ]
      daxpy_( &veclength, &alpha, Hvecs[ cnt ], &inc1, t_vec, &inc1 );
      daxpy_( &veclength, &alpha,  vecs[ cnt ], &inc1, u_vec, &inc1 );
   }
   double theEigenvalue = -mxM_eigs[0];
   daxpy_( &veclength, &theEigenvalue, u_vec, &inc1, t_vec, &inc1 );
   
   //8. Calculate the norm of r
   char frobenius = 'F';
   const double rnorm = dlange_( &frobenius, &veclength, &inc1, t_vec, &veclength, NULL ); // Work is not referenced for Frobenius norm
   return rnorm;

}

void CheMPS2::Davidson::CalculateNewVec(){

   int inc1 = 1;

   //9a. Calculate the new t_vec based on the residual of the lowest eigenvalue, to add to the vecs.
   for (int cnt = 0; cnt < veclength; cnt++){
      const double difference = diag[ cnt ] - mxM_eigs[0];
      const double fabsdiff   = fabs( difference );
      if ( fabsdiff > DIAG_CUTOFF ){
         work_vec[ cnt ] = u_vec[ cnt ] / difference; // work_vec = K^(-1) u_vec
      } else {
         work_vec[ cnt ] = u_vec[ cnt ] / DIAG_CUTOFF;
         if ( debugPrint ){ cout << "WARNING AT DAVIDSON : | (diag[" << cnt << "] - mxM_eigs[0]) | = " << fabsdiff << endl; }
      }
   }
   double alpha = - ddot_( &veclength, work_vec, &inc1, t_vec, &inc1 ) / ddot_( &veclength, work_vec, &inc1, u_vec, &inc1 ); // alpha = - (u^T K^(-1) r) / (u^T K^(-1) u)
   daxpy_( &veclength, &alpha, u_vec, &inc1, t_vec, &inc1 ); // t_vec = r - (u^T K^(-1) r) / (u^T K^(-1) u) u
   for (int cnt = 0; cnt < veclength; cnt++){
      const double difference = diag[ cnt ] - mxM_eigs[0];
      const double fabsdiff   = fabs( difference );
       if ( fabsdiff > DIAG_CUTOFF ){
         t_vec[ cnt ] = - t_vec[ cnt ] / difference; //t_vec = - K^(-1) (r - (u^T K^(-1) r) / (u^T K^(-1) u) u)
      } else {
         t_vec[ cnt ] = - t_vec[ cnt ] / DIAG_CUTOFF;
      }
   }

}

void CheMPS2::Davidson::Deflation(){

   int inc1 = 1;

   // 9b. When the maximum number of vectors is reached: construct the one with lowest eigenvalue & restart
   if ( NUM_VEC_KEEP <= 1 ){
   
      char frobenius = 'F';
      double alpha = 1.0 / dlange_( &frobenius, &veclength, &inc1, u_vec, &veclength, NULL ); // Work is not referenced for Frobenius norm
      dscal_( &veclength, &alpha, u_vec, &inc1 );
      dcopy_( &veclength, u_vec, &inc1, vecs[0], &inc1 );
   
   } else {
   
      if ( Reortho_Eigenvecs    == NULL ){ Reortho_Eigenvecs    = new double[    veclength * NUM_VEC_KEEP ]; }
      if ( Reortho_Overlap      == NULL ){ Reortho_Overlap      = new double[ NUM_VEC_KEEP * NUM_VEC_KEEP ]; }
      if ( Reortho_Overlap_eigs == NULL ){ Reortho_Overlap_eigs = new double[ NUM_VEC_KEEP                ]; }
      if ( Reortho_Lowdin       == NULL ){ Reortho_Lowdin       = new double[ NUM_VEC_KEEP * NUM_VEC_KEEP ]; }
   
      //Construct the lowest NUM_VEC_KEEP eigenvectors
      dcopy_( &veclength, u_vec, &inc1, Reortho_Eigenvecs, &inc1 );
      for (int cnt = 1; cnt < NUM_VEC_KEEP; cnt++){
         for (int irow = 0; irow < veclength; irow++){
            Reortho_Eigenvecs[ irow + veclength * cnt ] = 0.0;
            for (int ivec = 0; ivec < MAX_NUM_VEC; ivec++){
               Reortho_Eigenvecs[ irow + veclength * cnt ] += vecs[ ivec ][ irow ] * mxM_vecs[ ivec + MAX_NUM_VEC * cnt ];
            }
         }
      }
      
      //Calculate the overlap matrix
      char trans  = 'T';
      char notr   = 'N';
      double one  = 1.0;
      double zero = 0.0; //set
      dgemm_( &trans, &notr, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &veclength, &one, Reortho_Eigenvecs, &veclength, Reortho_Eigenvecs, &veclength, &zero, Reortho_Overlap, &NUM_VEC_KEEP );
      
      //Calculate the Lowdin tfo
      char jobz = 'V';
      char uplo = 'U';
      int info;
      dsyev_( &jobz, &uplo, &NUM_VEC_KEEP, Reortho_Overlap, &NUM_VEC_KEEP, Reortho_Overlap_eigs, mxM_work, &mxM_lwork, &info ); //Ascending order of eigs
      for (int icnt = 0; icnt < NUM_VEC_KEEP; icnt++){
         Reortho_Overlap_eigs[ icnt ] = pow( Reortho_Overlap_eigs[ icnt ], -0.25 );
         dscal_( &NUM_VEC_KEEP, Reortho_Overlap_eigs + icnt, Reortho_Overlap + NUM_VEC_KEEP * icnt, &inc1 );
      }
      dgemm_( &notr, &trans, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &NUM_VEC_KEEP, &one, Reortho_Overlap, &NUM_VEC_KEEP, Reortho_Overlap, &NUM_VEC_KEEP, &zero, Reortho_Lowdin, &NUM_VEC_KEEP );
      
      //Reortho: Put the Lowdin tfo eigenvecs in vecs
      for (int ivec = 0; ivec < NUM_VEC_KEEP; ivec++){
         for (int loop = 0; loop < veclength; loop++){ vecs[ ivec ][ loop ] = 0.0; }
         for (int ivec2 = 0; ivec2 < NUM_VEC_KEEP; ivec2++){
            daxpy_( &veclength, Reortho_Lowdin + ivec2 + NUM_VEC_KEEP * ivec, Reortho_Eigenvecs + veclength * ivec2, &inc1, vecs[ ivec ], &inc1 );
         }
      }
   }

   num_vec = 0;

}

void CheMPS2::Davidson::MxMafterDeflation(){

   int inc1 = 1;
   for (int ivec = 0; ivec < NUM_VEC_KEEP; ivec++){
      for (int ivec2 = ivec; ivec2 < NUM_VEC_KEEP; ivec2++){
         mxM[ ivec + MAX_NUM_VEC * ivec2 ] = ddot_( &veclength, vecs[ ivec ], &inc1, Hvecs[ ivec2 ], &inc1 );
         mxM[ ivec2 + MAX_NUM_VEC * ivec ] = mxM[ ivec + MAX_NUM_VEC * ivec2 ];
      }
   }

}


