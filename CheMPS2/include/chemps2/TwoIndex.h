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

#ifndef TWOINDEX_CHEMPS2_H
#define TWOINDEX_CHEMPS2_H

#include <string>
#include "Irreps.h"

namespace CheMPS2{
/** TwoIndex class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 8, 2013
    
    Container class for symmetric two-index tensors with Abelian point group symmetry (real character table; see Irreps.h): 1DMs and 1-particle matrix elements. */
   class TwoIndex{

      public:
      
         //! Constructor
         /** \param nGroup The symmetry group number (see Irreps.h)
             \param IrrepSizes Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep */
         TwoIndex(const int nGroup, const int * IrrepSizes);
         
         //! Destructor
         virtual ~TwoIndex();
         
         //! Set all one-body matrix elements to zero
         void Clear();
         
         //! Set an element
         /** \param irrep The irrep number (see Irreps.h)
             \param i The first index (within the symmetry block)
             \param j The second index (within the symmetry block)
             \param val The value to which the element of the matrix should be set */
         void set(const int irrep, const int i, const int j, const double val);

         //! Get an element
         /** \param irrep The irrep number (see Irreps.h)
             \param i The first index (within the symmetry block)
             \param j The second index (within the symmetry block)
             \return The value of the matrix element */
         double get(const int irrep, const int i, const int j) const;
         
         //! Save the TwoIndex object
         /** \param name filename */
         void save(const std::string name) const;
         
         //! Load the TwoIndex object
         /** \param name filename */
         void read(const std::string name);
      
      private:
      
         //Contains the group number, the number of irreps, and the multiplication table
         Irreps SymmInfo;
         
         //Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep
         int * Isizes;
         
         //storage[I_i][i+j*(j+1)/2] = mx_ij = mx_ji
         double ** storage;

   };
}

#endif

