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

#ifndef TENSORDIAG_CHEMPS2_H
#define TENSORDIAG_CHEMPS2_H

#include "Tensor.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorDiag class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 15, 2013
    
    The TensorDiag class is a storage class for boundary tensors that are entirely block-diagonal:
     - the L and R parts of the LQ and QR decompositions of TensorT
     - the storage part of TensorX and TensorGYZ */
   class TensorDiag : public Tensor{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param denBKIn The problem to be solved */
         TensorDiag(const int indexIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorDiag();
         
         //! Get the number of symmetry blocks
         /** \return The number of symmetry blocks */
         int gNKappa() const;
         
         //! Get the pointer to the storage
         /** return pointer to the storage */
         double * gStorage();
         
         //! Get the index corresponding to a certain tensor block
         /** \param N1 The left or up particle number sector
             \param TwoS1 The left or up spin symmetry sector
             \param I1 The left or up irrep sector
             \param N2 The right or down particle number sector
             \param TwoS2 The right or down spin symmetry sector
             \param I2 The right or down irrep sector
             \return The kappa corresponding to the input parameters; -1 means no such block */
         int gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const;
         
         //! Get the storage jump corresponding to a certain tensor block
         /** \param kappa The symmetry block
             \return kappa2index[kappa], the memory jumper to a certain block */
         int gKappa2index(const int kappa) const;
         
         //! Get the pointer to the storage of a certain tensor block
         /** \param N1 The left or up particle number sector
             \param TwoS1 The left or up spin symmetry sector
             \param I1 The left or up irrep sector
             \param N2 The right or down particle number sector
             \param TwoS2 The right or down spin symmetry sector
             \param I2 The right or down irrep sector
             \return Pointer to the storage of the specified tensor block; NULL means no such block */
         double * gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2);
         
         //! Get the location index
         /** \return the index */
         int gIndex() const;
         
      protected:
      
         //! Update the previous TensorDiag with a TensorT and add the result to this TensorDiag, when moving right.
         /** \param ikappa Symmetry block of this TensorDiag which is updated
             \param denT TensorT to update the TensorDiag
             \param diagPrevious The previous TensorDiag
             \param workmemLR Work memory */
         void updateRight(const int ikappa, Tensor * denT, TensorDiag * diagPrevious, double * workmemLR);
         
         //! Update the previous TensorDiag with a TensorT and add the result to this TensorDiag, when moving left.
         /** \param ikappa Symmetry block of this TensorDiag which is updated
             \param denT TensorT to update the TensorDiag
             \param diagPrevious The previous TensorDiag
             \param workmemLR Work memory */
         void updateLeft(const int ikappa, Tensor * denT, TensorDiag * diagPrevious, double * workmemLR);
         
   };
}

#endif
