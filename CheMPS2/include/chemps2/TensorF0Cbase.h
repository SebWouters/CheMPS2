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

#ifndef TENSORF0CBASE_CHEMPS2_H
#define TENSORF0CBASE_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorF0Cbase class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 1, 2013
    
    The TensorF0Cbase class is a storage class for TensorF0 and TensorC. It also performs the update of a previous TensorF0Cbase based on movingRight.*/
   class TensorF0Cbase : public Tensor{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn Direct product of irreps of the two 2nd quantized operators; both sandwiched & to sandwich
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn The problem to be solved */
         TensorF0Cbase(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorF0Cbase();
         
         //! Get the number of symmetry blocks
         /** \return The number of symmetry blocks */
         int gNKappa() const;
         
         //! Get the pointer to the storage
         /** return pointer to the storage */
         double * gStorage();
         
         //! Get the index corresponding to a certain tensor block
         /** \param N1 The up particle number sector
             \param TwoS1 The up spin symmetry sector
             \param I1 The up irrep sector
             \param N2 The down particle number sector
             \param TwoS2 The down spin symmetry sector
             \param I2 The down irrep sector
             \return The kappa corresponding to the input parameters; -1 means no such block */
         int gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const;
         
         //! Get the storage jump corresponding to a certain tensor block
         /** \param kappa The symmetry block
             \return kappa2index[kappa], the memory jumper to a certain block */
         int gKappa2index(const int kappa) const;
         
         //! Get the pointer to the storage of a certain tensor block
         /** \param N1 The up particle number sector
             \param TwoS1 The up spin symmetry sector
             \param I1 The up irrep sector
             \param N2 The down particle number sector
             \param TwoS2 The down spin symmetry sector
             \param I2 The down irrep sector
             \return Pointer to the storage of the specified tensor block; NULL means no such block */
         double * gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2);
         
         //! Get the boundary index
         /** \return the index */
         int gIndex() const;
         
         //! Get Idiff
         /** \return Idiff */
         int gIdiff() const;
         
         //! Clear and update
         /** \param F0CbasePrevious TensorF0Cbase needed for the update
             \param denT TensorT needed for the update.
             \param workmem Work memory */
         void update(TensorF0Cbase * F0CbasePrevious, TensorT * denT, double * workmem);
         
      protected:
      
         //! The irrep difference (direct product of the two 2nd quantized operators; both sandwiched & to sandwich)
         int Idiff;
        
         //! Whether or not moving right
         bool movingRight;
         
         //! Set all storage variables to 0.0
         void Clear();
         
      private:
         
         //update when movingright
         void updateMovingRight(TensorF0Cbase * F0CbasePrevious, TensorT * denT, double * workmem);
         
         //update when movingleft
         void updateMovingLeft(TensorF0Cbase * F0CbasePrevious, TensorT * denT, double * workmem);

   };
}
#endif
