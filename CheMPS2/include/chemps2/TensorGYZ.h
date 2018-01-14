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

#ifndef TENSORGYZ_CHEMPS2_H
#define TENSORGYZ_CHEMPS2_H

#include "TensorT.h"
#include "Problem.h"
#include "TensorOperator.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorGYZ class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 9, 2014
    
    The TensorGYZ class is a storage and manipulation class for the two-orbital mutual information. It only exists moving left to right. */
   class TensorGYZ : public TensorOperator{

      public:
      
         //! Constructor
         /** \param boundary_index The boundary index
             \param identity The identity: G, Y, Z (capitals!)
             \param denBK The symmetry bookkeeper with symmetry sector virtual dimensions */
         TensorGYZ(const int boundary_index, const char identity, const SyBookkeeper * denBK);
         
         //! Destructor
         virtual ~TensorGYZ();
         
         //! Construct a new TensorGYZ
         /** \param denT TensorT from which the new TensorGYZ should be made */
         void construct(TensorT * denT);
         
      private:
         
         //! The identity
         char identity;
         
   };
}

#endif
