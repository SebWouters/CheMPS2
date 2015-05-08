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
         /** \param veclength_in Linear dimension of the symmetric matrix, or the length of the vectors
             \param MAX_NUM_VEC_in The maximum number of vectors in which the symmetric matrix is approximately diagonalized
             \param NUM_VEC_KEEP_in The number of vectors to keep on deflation
             \param RTOL_in The tolerance for the two-norm of the residual (for convergence)
             \param DIAG_CUTOFF_in Cutoff value for the diagonal preconditioner
             \param debugPrint_in Whether or not to debug print */
         Davidson(const int veclength_in, const int MAX_NUM_VEC_in, const int NUM_VEC_KEEP_in, const double RTOL_in, const double DIAG_CUTOFF_in, const bool debugPrint_in);
         
         //! Destructor
         virtual ~Davidson();
         
         //! The iterator to converge the ground state vector
         /** \param whichpointers Array of double* of length 2 to return pointers to vectors to the caller
             \return Instruction character. 'A' means copy the initial guess to whichpointers[0] and the diagonal of the symmetric matrix to whichpointers[1]. 'B' means calculate whichpointers[1] as the result of multiplying the symmetric matrix with whichpointers[0]. 'C' means that the converged solution can be copied back from whichpointers[0], and the ground-state energy from whichpointers[1][0]. 'D' means that an error has occurred. */
         char FetchInstruction(double ** whichpointers);
         
         //! Get the number of matrix vector multiplications which have been performed
         /** \return The number of matrix vector multiplications which have been performed */
         int GetNumMultiplications() const;
         
      private:
      
         int veclength; // The vector length
         int nMultiplications; // Current number of requested matrix-vector multiplications
         char state; // Current state of the algorithm --> based on this parameter the next instruction is given
         bool debugPrint;
         
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
         
         // Vector spaces
         double * t_vec;
         double * u_vec;
         double * work_vec;
         double * diag;
         
         // For the deflation
         double * Reortho_Lowdin;
         double * Reortho_Overlap_eigs;
         double * Reortho_Overlap;
         double * Reortho_Eigenvecs;
         
         // Control script functions
         void SafetyCheckGuess();
         void AddNewVec();
         double DiagonalizeSmallMatrixAndCalcResidual(); // Returns the residual norm
         void CalculateNewVec();
         void Deflation();
         void MxMafterDeflation();
         
   };
}

#endif

