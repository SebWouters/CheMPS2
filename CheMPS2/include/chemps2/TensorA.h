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

#ifndef TENSORA_CHEMPS2_H
#define TENSORA_CHEMPS2_H

#include "TensorS0Abase.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorA class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 4, 2013
    
    The TensorA class is a storage class for the complementary operator of the spin-0 component of two contracted creators or two contracted annihilators. */
   class TensorA : public TensorS0Abase{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn Direct product of irreps of the two 2nd quantized operators; both sandwiched & to sandwich
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn The problem to be solved */
         TensorA(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorA();
         
         //! Clear the tensor when update is not yet possible
         void ClearStorage();
         
         //! Add a term
         /** \param alpha prefactor
             \param TermToAdd The TensorS0Abase to add */
         void AddATerm(double alpha, TensorS0Abase * TermToAdd);

   };
}

#endif
