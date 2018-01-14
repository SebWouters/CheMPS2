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

#ifndef CONJUGATEGRADIENT_CHEMPS2_H
#define CONJUGATEGRADIENT_CHEMPS2_H

namespace CheMPS2{
/** Conjugate gradient class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date January 27, 2016

    The ConjugateGradient class implements the conjugate gradient algorithm to solve the symmetric linear problem \n

        \f$ operator * x = b \f$. \n

    With precon = 1 / sqrt( diag( operator ) ), the problem is turned into \n

        \f$ precon * operator * precon * xtilde = precon * b \f$ \n
        \f$ x = precon * xtilde \f$.
*/
   class ConjugateGradient{

      public:

         //! Constructor
         /** \param veclength_in Linear dimension of the symmetric matrix
             \param RTOL_in The tolerance for the two-norm of the residual
             \param DIAG_CUTOFF_in The cutoff to truncate the diagonal elements of operator
             \param print_in Whether or not to print */
         ConjugateGradient(const int veclength_in, const double RTOL_in, const double DIAG_CUTOFF_in, const bool print_in);

         //! Destructor
         virtual ~ConjugateGradient();

         //! The iterator to converge the ground state vector
         /** \param pointers Array of double* of length 3 to return pointers to vectors to the caller
             \return Instruction character. 'A' means copy the initial guess to pointers[0], the diagonal of the symmetric matrix to pointers[1], and the right-hand side of the problem to pointers[2]. 'B' means calculate pointers[1] = symmetric matrix times pointers[0]. 'C' means that the converged solution can be copied back from pointers[0], and the residual norm from pointers[1][0]. 'D' means that an error has occurred. */
         char step( double ** pointers );

         //! Get the number of matrix vector multiplications which have been performed
         /** \return The number of matrix vector multiplications which have been performed */
         int get_num_matvec() const;

      private:

         int veclength;
         double RTOL;
         double DIAG_CUTOFF;
         bool print;

         char state;     // Current state of the algorithm
         int num_matvec; // Current number of matvec multiplications

         // Helper arrays
         double * XVEC;
         double * PRECON;
         double * RHS;
         double * WORK;
         double * RESID;
         double * PVEC;
         double * OPVEC;

         // Helper variables
         double rnorm;
         double rdotr;

         // Internal functions to hop between states
         void stepL2K();
         void stepY2Z();
         void stepJ2K();
         void stepG2H();
         double inprod( double * vector );
         double inprod( double * vector, double * othervector );
         void apply_precon( double * vector );
         void apply_precon( double * vector, double * result );

   };
}

#endif

