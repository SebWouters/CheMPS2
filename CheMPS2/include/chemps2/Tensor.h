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

#ifndef TENSOR_CHEMPS2_H
#define TENSOR_CHEMPS2_H

#include "SyBookkeeper.h"

namespace CheMPS2{
/** Pure virtual Tensor class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 15, 2013

    The Tensor class defines parameters and functions which all Tensors must have. */
   class Tensor{

      public:

         //! Get the number of tensor blocks
         /** return The number of tensor blocks */
         virtual int gNKappa() const = 0;

         //! Get the pointer to the storage
         /** return pointer to the storage */
         virtual double * gStorage() = 0;

         //! Get the index corresponding to a certain tensor block
         /** \param N1 The left or up particle number sector
             \param TwoS1 The left or up spin symmetry sector
             \param I1 The left or up irrep sector
             \param N2 The right or down particle number sector
             \param TwoS2 The right or down spin symmetry sector
             \param I2 The right or down irrep sector
             \return The kappa corresponding to the input parameters; -1 means no such block */
         virtual int gKappa( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 ) const = 0;

         //! Get the storage jump corresponding to a certain tensor block
         /** \param kappa The symmetry block
             \return kappa2index[ kappa ], the memory jumper to a certain block */
         virtual int gKappa2index( const int kappa ) const = 0;

         //! Get the pointer to the storage of a certain tensor block
         /** \param N1 The left or up particle number sector
             \param TwoS1 The left or up spin symmetry sector
             \param I1 The left or up irrep sector
             \param N2 The right or down particle number sector
             \param TwoS2 The right or down spin symmetry sector
             \param I2 The right or down irrep sector
             \return Pointer to the storage of the specified tensor block; NULL means no such block */
         virtual double * gStorage( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 ) = 0;

         //! Get the location index
         /** \return the index */
         virtual int gIndex() const = 0;

      protected:

         //! Index of the Tensor object. For TensorT: a site index; for other tensors: a boundary index
         int index;

         //! The actual variables. Tensor block kappa begins at storage+kappa2index[kappa] and ends at storage+kappa2index[kappa+1].
         double * storage;

         //! Number of Tensor blocks.
         int nKappa;

         //! kappa2index[kappa] indicates the start of tensor block kappa in storage. kappa2index[nKappa] gives the size of storage.
         int * kappa2index;

   };
}

#endif
