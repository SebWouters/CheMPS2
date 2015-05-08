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

#ifndef TENSORT_CHEMPS2_H
#define TENSORT_CHEMPS2_H

#include "Tensor.h"
#include "TensorDiag.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorT class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 18, 2013
    
    The TensorT class is a storage class for MPS tensors. Extra functionalities are
     - left- and right-orthonormalization
     - left- and right-multiplication with the R-/L- matrices from the orthonormalization */
   class TensorT : public Tensor{

      public:
      
         //! Constructor
         /** \param indexIn The site index
             \param IlocalIn The local irrep index
             \param denBKIn The problem to be solved */
         TensorT(const int indexIn, const int IlocalIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorT();
         
         //! Get the number of symmetry blocks
         /** \return The number of symmetry blocks */
         int gNKappa() const;
         
         //! Get the pointer to the storage
         /** return pointer to the storage */
         double * gStorage();
         
         //! Get the index corresponding to a certain tensor block
         /** \param N1 The left particle number sector
             \param TwoS1 The left spin symmetry sector
             \param I1 The left irrep sector
             \param N2 The right particle number sector
             \param TwoS2 The right spin symmetry sector
             \param I2 The right irrep sector
             \return The kappa corresponding to the input parameters; -1 means no such block */
         int gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const;
         
         //! Get the storage jump corresponding to a certain tensor block
         /** \param kappa The symmetry block
             \return kappa2index[kappa], the memory jumper to a certain block */
         int gKappa2index(const int kappa) const;
         
         //! Get the pointer to the storage of a certain tensor block
         /** \param N1 The left particle number sector
             \param TwoS1 The left spin symmetry sector
             \param I1 The left irrep sector
             \param N2 The right particle number sector
             \param TwoS2 The right spin symmetry sector
             \param I2 The right irrep sector
             \return Pointer to the storage of the specified tensor block; NULL means no such block */
         double * gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2);
         
         //! Get the location index
         /** \return the index */
         int gIndex() const;
         
         //! Fill storage with random numbers 0 < val < 1.
         void random();
         
         //! Left-normalization
         /** \param Rstorage Where the R-part of the QR-decomposition can be stored. */
         void QR(TensorDiag * Rstorage);
         
         //! Right-normalization
         /** \param Lstorage Where the L-part of the LQ-decomposition can be stored. */
         void LQ(TensorDiag * Lstorage);
         
         //! Multiply at the left with a TensorDiag
         /** \param Mx The TensorDiag with which the current TensorT should be multiplied at the left */
         void LeftMultiply(TensorDiag * Mx);
         
         //! Multiply at the right with a TensorDiag
         /** \param Mx The TensorDiag with which the current TensorT should be multiplied at the right */
         void RightMultiply(TensorDiag * Mx);
         
         //! Reset the TensorT (if virtual dimensions are changed)
         void Reset();
         
         //! Check whether the TensorT is left-normal
         /** \return Whether TensorT is left-normal */
         bool CheckLeftNormal() const;
         
         //! Check whether the TensorT is right-normal
         /** \return Whether TensorT is right-normal */
         bool CheckRightNormal() const;
         
      private:
      
         //The local irrep
         int Ilocal;
      
         //right particle number sector
         int * sectorNR;
         
         //right spin sector
         int * sectorTwoSR;
         
         //right irrep sector
         int * sectorIR;
         
         //Delete all arrays
         void DeleteAllArrays();
         
         //Allocate all arrays
         void AllocateAllArrays();
         
   };
}

#endif
