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

#ifndef TENSORS1_CHEMPS2_H
#define TENSORS1_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "TensorL.h"
#include "TensorS1Bbase.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorS1 class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 1, 2013
    
    The TensorS1 class is a storage and manipulation class for the spin-1 component of two contracted creators or two contracted annihilators. */
   class TensorS1 : public TensorS1Bbase{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn Direct product of irreps of the two 2nd quantized operators; both sandwiched & to sandwich
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn The problem to be solved */
         TensorS1(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorS1();
         
         //Make new TensorS1 (vs update)
         /** \param denL TensorL from which the new TensorS1 should be made.
             \param denT TensorT from which the new TensorS1 should be made.
             \param workmem Work memory */
         void makenew(TensorL * denL, TensorT * denT, double * workmem);
         
      private:
         
         //makenew when movingright
         void makenewRight(TensorL * denL, TensorT * denT, double * workmem);
      
         //makenew when movingleft
         void makenewLeft(TensorL * denL, TensorT * denT, double * workmem);
         
   };
}

#endif
