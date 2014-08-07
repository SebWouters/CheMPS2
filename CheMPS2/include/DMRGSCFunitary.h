/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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

#ifndef DMRGSCFUNITARY_H
#define DMRGSCFUNITARY_H

#include "Options.h"
#include "DMRGSCFindices.h"
#include "DIIS.h"

namespace CheMPS2{
/** DMRGSCFunitary class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 11, 2014
    
    The DMRGSCFunitary class is a storage and manipulation class for the DMRGSCF unitary matrix. This matrix is blockdiagonal in the irreducible representations, and is formed by stepwise multiplying in new unitary rotations due to the Newton-Raphson algorithm.
*/
   class DMRGSCFunitary{

      public:
      
         //! Constructor
         /** \param iHandlerIn The DMRGSCF indices */
         DMRGSCFunitary(DMRGSCFindices * iHandlerIn);
         
         //! Destructor
         virtual ~DMRGSCFunitary();
         
         //! Get the number of variables in the x-parametrization of the unitary update
         /** \return The number of unique variables in the x-matrix */
         int getNumVariablesX() const;
         
         //! Get the first Hamiltonian index corresponding to linearindex
         /** \param linearindex The linear index of the x-parametrization
             \return The first Hamiltonian index corresponding to linearindex */
         int getFirstIndex(const int linearindex) const;
         
         //! Get the second Hamiltonian index corresponding to linearindex
         /** \param linearindex The linear index of the x-parametrization
             \return The second Hamiltonian index corresponding to linearindex */
         int getSecondIndex(const int linearindex) const;
         
         //! Get the unitary rotation for block irrep
         /** \param irrep The irreducible representation
             \return Pointer to the desired unitary block */
         double * getBlock(const int irrep);
         
         //! Update the unitary transformation
         /** \param workmem1 Work memory
             \param workmem2 Work memory
             \param vector The elements in X
             \param multiply Boolean whether exp(X)*U or exp(X) should become the new U. If multiply==true, U <-- exp(X)*U. If multiply==false, U <-- exp(X).
             \param compact Boolean which indicates how the elements X are stored */
         void updateUnitary(double * workmem1, double * workmem2, double * vector, const bool multiply, const bool compact);
         
         //! Rotate the unitary matrix to the NO eigenbasis
         /** \param eigenvecs The NO eigenbasis
             \param work Work memory */
         void rotateUnitaryNOeigenvecs(double * eigenvecs, double * work);
         
         //! Calculate the two-norm of U^T*U - I
         /** \param work Work memory */
         void CheckDeviationFromUnitary(double * work) const;
         
         //! Obtain the logarithm of the unitary matrix
         /** \param vector Where the logarithm should be stored
             \param temp1 Work memory
             \param temp2 Work memory */
         void getLog(double * vector, double * temp1, double * temp2) const;
         
         //! Obtain the logarithm of the current unitary matrix based on the BCH formula
         /** \param Xprev The logarithm of the previous unitary matrix
             \param step The update based on the gradient and the Hessian
             \param Xnew The approximated logarithm of the current unitary matrix (write)
             \param temp1 Work memory
             \param temp2 Work memory */
         void BCH(double * Xprev, double * step, double * Xnew, double * temp1, double * temp2) const;
         
         //! Orbitals are defined up to a phase factor. Make sure that the logarithm of each block of the unitary has determinant 1.
         /** \param temp1 Work memory
             \param temp2 Work memory */
         void makeSureAllBlocksDetOne(double * temp1, double * temp2);
         
         //! Save the unitary to disk
         void saveU() const;
         
         //! Load the unitary from disk
         void loadU();
         
         //! Delete the stored unitary (on disk)
         void deleteStoredUnitary() const;

         
      private:
      
         //Externally created and destroyed index handler
         DMRGSCFindices * iHandler;
         
         //Number of variables in the x-matrix
         int x_linearlength;
         
         //Helper arrays to jump from linear x-matrix index to orbital indices and back
         int * x_firstindex;
         int * x_secondindex;
         
         //The unitary matrix (e^x * previous unitary): unitary[irrep][row + size_irrep * col]
         double ** unitary;
         
         // Find the linear index corresponding to p and q
         /** \param p_index The first Hamiltonian index
             \param q_index The second Hamiltonian index
             \return The linear index corresponding to (p,q). If no index is found -1 is returned. */
         int getLinearIndex(const int p_index, const int q_index) const;
         
         // Build in result the skew symmetric matrix X for irrep block irrep based on the elements in Xelem. If compact==true, they are stored in gradient form.
         void buildSkewSymmX(const int irrep, double * result, double * Xelem, const bool compact) const;
         
   };
}

#endif
