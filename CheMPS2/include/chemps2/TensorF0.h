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

#ifndef TENSORF0_CHEMPS2_H
#define TENSORF0_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "TensorL.h"
#include "TensorOperator.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorF0 class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 1, 2013
    
    The TensorF0 class is a storage and manipulation class for the spin-0 component of a contracted creator & annihilator. */
   class TensorF0 : public TensorOperator{

      public:
      
         //! Constructor
         /** \param boundary_index The boundary index
             \param Idiff Direct product of irreps of the two 2nd quantized operators; both sandwiched & to sandwich
             \param moving_right If true: sweep from left to right. If false: sweep from right to left
             \param denBK The symmetry sector bookkeeper */
         TensorF0( const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * denBK );
         
         //! Destructor
         virtual ~TensorF0();
         
         //Make new TensorF0 (vs update)
         /** \param denT TensorT from which the new TensorF0 should be made. */
         void makenew(TensorT * denT);
         
         //Make new TensorF0 (vs update)
         /** \param denL TensorL from which the new TensorF0 should be made.
             \param denT TensorT from which the new TensorF0 should be made.
             \param workmem Work memory */
         void makenew(TensorL * denL, TensorT * denT, double * workmem);
         
      private:
         
         //makenew when movingright
         void makenewRight(TensorL * denL, TensorT * denT, double * workmem);
      
         //makenew when movingleft
         void makenewLeft(TensorL * denL, TensorT * denT, double * workmem);
         
         //makenew when movingright
         void makenewRight(TensorT * denT);
      
         //makenew when movingleft
         void makenewLeft(TensorT * denT);
         
   };
}

#endif
