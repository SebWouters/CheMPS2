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

#ifndef TENSORL_CHEMPS2_H
#define TENSORL_CHEMPS2_H

#include "TensorT.h"
#include "TensorO.h"

namespace CheMPS2{
/** TensorL class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 20, 2013

    The TensorL class is a storage and manipulation class for a single contracted creator/annihilitor. */
   class TensorL : public TensorOperator{

      public:

         //! Constructor
         /** \param boundary_index The boundary index
             \param Idiff          The irrep of the one creator ( sandwiched if TensorL ; to sandwich if TensorQ )
             \param moving_right   If true: sweep from left to right. If false: sweep from right to left
             \param book_up        Symmetry bookkeeper of the upper MPS
             \param book_down      Symmetry bookkeeper of the lower MPS */
         TensorL( const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * book_up, const SyBookkeeper * book_down );

         //! Destructor
         virtual ~TensorL();

         //! Create a new TensorL
         /** \param mps_tensor TensorT from which the new TensorL should be made. Please not that this function assumes book_up == book_down. */
         void create( TensorT * mps_tensor );

         //! Create a new TensorL
         /** \param mps_tensor_up   Upper TensorT from which the new TensorL should be made
             \param mps_tensor_down Lower TensorT from which the new TensorL should be made
             \param previous        Overlap matrix on the previous edge
             \param workmem         Work memory of size max(dimLup,down) * max(dimRup,down) */
         void create( TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem );

      private:

         //! Make new when moving_right == true
         void create_right( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem );

         //! Make new when moving_right == false
         void create_left( const int ikappa, TensorT * mps_tensor_up, TensorT * mps_tensor_down, TensorO * previous, double * workmem );

   };
}

#endif
