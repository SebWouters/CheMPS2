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

#ifndef TENSORL_CHEMPS2_H
#define TENSORL_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "TensorSwap.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorL class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 20, 2013
    
    The TensorL class is a storage and manipulation class for a single contracted creator/annihilitor. */
   class TensorL : public TensorSwap{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn The irrep of the one creator ( sandwiched if TensorL ; to sandwich if TensorQ )
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn The problem to be solved */
         TensorL(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorL();
         
         //Make new TensorL (vs update)
         /** \param denT TensorT from which the new TensorL should be made. */
         void makenew(TensorT * denT);
         
      private:
      
         //makenew when movingright
         void makenewRight(TensorT * denT);
      
         //makenew when movingleft
         void makenewLeft(TensorT * denT);
         
         
   };
}

#endif
