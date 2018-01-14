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

#ifndef TENSORO_CHEMPS2_H
#define TENSORO_CHEMPS2_H

#include "TensorOperator.h"
#include "TensorT.h"
#include "Problem.h"

namespace CheMPS2{
/** TensorO class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 31, 2013

    The TensorO class is a storage class for overlaps between different MPSs. */
   class TensorO : public TensorOperator{

      public:

         //! Constructor
         /** \param boundary_index The boundary index
             \param moving_right If true: sweep from left to right. If false: sweep from right to left
             \param book_up   The symmetry bookkeeper with the upper symmetry sector virtual dimensions
             \param book_down The symmetry bookkeeper with the lower symmetry sector virtual dimensions */
         TensorO( const int boundary_index, const bool moving_right, const SyBookkeeper * book_up, const SyBookkeeper * book_down );

         //! Destructor
         virtual ~TensorO();

         //! Clear and add the relevant terms to the TensorO
         /** \param mps_tensor_up   Upper MPS tensor from which the new TensorO should be made
             \param mps_tensor_down Lower MPS tensor from which the new TensorO should be made */
         void create( TensorT * mps_tensor_up, TensorT * mps_tensor_down );

         //! Update the previous TensorO
         /** \param mps_tensor_up   Upper MPS tensor from which the update should be made
             \param mps_tensor_down Lower MPS tensor from which the update should be made
             \param previous        Previous TensorO from which the update should be made */
         void update_ownmem( TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous );

      private:

         //helper functions
         void create_right( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down );
         void create_left(  const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down );

   };
}

#endif
