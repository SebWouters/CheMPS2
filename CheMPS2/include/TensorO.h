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

#ifndef TENSORO_CHEMPS2_H
#define TENSORO_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "Problem.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorO class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 31, 2013
    
    The TensorO class is a storage class for overlaps between different MPSs. */
   class TensorO : public Tensor{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKupIn The symmetry bookkeeper with the upper symmetry sector virtual dimensions (old MPS)
             \param denBKdownIn The symmetry bookkeeper with the lower symmetry sector virtual dimensions (current MPS)
             \param ProbIn The Problem containing the Hamiltonian matrix elements */
         TensorO(const int indexIn, const bool movingRightIn, const SyBookkeeper * denBKupIn, const SyBookkeeper * denBKdownIn, const Problem * ProbIn);
         
         //! Destructor
         virtual ~TensorO();
         
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
         
         //! Clear and add the relevant terms to the TensorO
         /** \param denTup Upper TensorT from which the new TensorO should be made (old MPS)
             \param denTdown Lower TensorT from which the new TensorO should be made (current MPS)
             \param denO The previous TensorO */
         void update(TensorT * denTup, TensorT * denTdown, TensorO * denO);
         
         //! Clear and add the relevant terms to the TensorO
         /** \param denTup Upper TensorT from which the new TensorO should be made (old MPS)
             \param denTdown Lower TensorT from which the new TensorO should be made (current MPS) */
         void update(TensorT * denTup, TensorT * denTdown);
         
      private:
      
         //Externally created and destroyed BK for the old MPS
         const SyBookkeeper * denBKup;
      
         //whether moving right or not
         bool movingRight;
         
         //Problem containing the matrix elements
         const Problem * Prob;
         
         //Clear the memory
         void Clear();
      
         //helper functions
         void updateRight(const int ikappa, TensorT * denTup, TensorT * denTdown);
         void updateLeft(const int ikappa, TensorT * denTup, TensorT * denTdown);
         void updateRight(const int ikappa, TensorT * denTup, TensorT * denTdown, TensorO * denO, double * workmem);
         void updateLeft(const int ikappa, TensorT * denTup, TensorT * denTdown, TensorO * denO, double * workmem);
         
   };
}

#endif
