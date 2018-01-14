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

#ifndef DAVIDSON_CHEMPS2_H
#define DAVIDSON_CHEMPS2_H

namespace CheMPS2{
/** Davidson class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date January 29, 2015
    
    The Davidson class implements Davidson's algorithm to find the lowest eigenvalue and corresponding eigenvector of a symmetric operator.
    Information can be found in \n
     
     [1] E.R. Davidson, J. Comput. Phys. 17 (1), 87-94 (1975). http://dx.doi.org/10.1016/0021-9991(75)90065-0 \n
     [2] http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter11.pdf (In this class algorithm 11.1 is implemented, with equation (11.3) instead of line (16).)
*/
   class Davidson{

      public:

         //! Constructor
         /** \param veclength    Linear dimension of the symmetric matrix, or the length of the vectors
             \param MAX_NUM_VEC  The maximum number of vectors in which the symmetric matrix is approximately diagonalized
             \param NUM_VEC_KEEP The number of vectors to keep during deflation
             \param RTOL         The tolerance for the two-norm of the residual ( for convergence )
             \param DIAG_CUTOFF  Cutoff value for the diagonal preconditioner
             \param debug_print  Whether or not to debug print
             \param problem_type 'E' for eigenvalue or 'L' for linear problem. */
         Davidson( const int veclength, const int MAX_NUM_VEC, const int NUM_VEC_KEEP, const double RTOL, const double DIAG_CUTOFF, const bool debug_print, const char problem_type = 'E' );

         //! Destructor
         virtual ~Davidson();

         //! The iterator to converge the ground state vector
         /** \param pointers Array of double* of length 2 when problem_type=='E' or length 3 when problem_type=='L'.
             \return Instruction character. 'A' means copy the initial guess to pointers[0] and the diagonal of the symmetric matrix to pointers[1]. If 'A' and problem_type=='E', the right-hand side of the problem should be copied to pointers[2]. 'B' means calculate pointers[1] as the result of multiplying the symmetric matrix with pointers[0]. 'C' means that the converged solution can be copied back from pointers[0], and pointers[1][0] contains the ground-state energy if problem_type=='E' or the residual norm if problem_type=='L'. 'D' means that an error has occurred. */
         char FetchInstruction( double ** pointers );

         //! Get the number of matrix vector multiplications which have been performed
         /** \return The number of matrix vector multiplications which have been performed */
         int GetNumMultiplications() const;

      private:

         int veclength; // The vector length
         int nMultiplications; // Current number of requested matrix-vector multiplications
         char state; // Current state of the algorithm --> based on this parameter the next instruction is given
         bool debug_print;
         char problem_type;

         // Davidson parameters
         int MAX_NUM_VEC;
         int NUM_VEC_KEEP;
         double DIAG_CUTOFF;
         double RTOL;

         // To store the vectors and the matrix x vectors
         int num_vec;
         double ** vecs;
         double ** Hvecs;
         int num_allocated;

         // The effective diagonalization problem
         double * mxM;
         double * mxM_eigs;
         double * mxM_vecs;
         int mxM_lwork;
         double * mxM_work;
         double * mxM_rhs;

         // Vector spaces
         double * t_vec;
         double * u_vec;
         double * work_vec;
         double * diag;
         double * RHS;

         // For the deflation
         double * Reortho_Lowdin;
         double * Reortho_Overlap_eigs;
         double * Reortho_Overlap;
         double * Reortho_Eigenvecs;

         // Control script functions
         double FrobeniusNorm( double * current_vector );
         void SafetyCheckGuess();
         void AddNewVec();
         double DiagonalizeSmallMatrixAndCalcResidual(); // Returns the residual norm
         void CalculateNewVec();
         void Deflation();
         void MxMafterDeflation();
         void SolveLinearSystemDeflation( const int NUM_SOLUTIONS );

   };
}

#endif

