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

#ifndef TENSORM_CHEMPS2_H
#define TENSORM_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "TensorSwap.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorM class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 9, 2014
    
    The TensorM class is a storage and manipulation class for a single contracted creator/annihilitor for the two-orbital mutual information. It only exists moving left to right. */
   class TensorM : public TensorSwap{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn The irrep of the one creator ( sandwiched since TensorM )
             \param denBKIn The symmetry bookkeeper of the problem to be solved */
         TensorM(const int indexIn, const int IdiffIn, const SyBookkeeper * denBKIn);
         
         //! Destructor
         virtual ~TensorM();
         
         //Construct new TensorM (vs update)
         /** \param denT TensorT from which the new TensorM should be made. */
         void construct(TensorT * denT);
         
   };
}

#endif
