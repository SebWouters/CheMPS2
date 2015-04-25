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

#ifndef TWODMSTORAGE_CHEMPS2_H
#define TWODMSTORAGE_CHEMPS2_H

#include "Irreps.h"

namespace CheMPS2{
/** TwoDMstorage class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date June 20, 2014
    
    Container class for 2DM storage. The 2DM-A and 2DM-B have 4-fold permutation symmetry. Currently only the sparsity due to the Abelian point group symmetry (real character table; see Irreps.h) is taken into account.
*/
   class TwoDMstorage{

      public:
      
         //! Constructor
         /** \param nGroup The symmetry group number (see Irreps.h)
             \param IrrepSizes Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep */
         TwoDMstorage(const int nGroup, const int * IrrepSizes);
         
         //! Destructor
         virtual ~TwoDMstorage();
         
         //! Set an element
         /** \param irrep_i The irrep number of the first orbital (see Irreps.h)
             \param irrep_j The irrep number of the second orbital
             \param irrep_k The irrep number of the third orbital
             \param irrep_l The irrep number of the fourth orbital
             \param i The first index (within the symmetry block)
             \param j The second index (within the symmetry block)
             \param k The third index (within the symmetry block)
             \param l The fourth index (within the symmetry block)
             \param val The value to which the element of the matrix should be set */
         void set(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val);

         //! Get an element
         /** \param irrep_i The irrep number of the first orbital (see Irreps.h)
             \param irrep_j The irrep number of the second orbital
             \param irrep_k The irrep number of the third orbital
             \param irrep_l The irrep number of the fourth orbital
             \param i The first index (within the symmetry block)
             \param j The second index (within the symmetry block)
             \param k The third index (within the symmetry block)
             \param l The fourth index (within the symmetry block) */
         double get(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
         
         //! Save the TwoDMstorage object
         /** \param name filename */
         void save(const std::string name) const;
         
         //! Load the TwoDMstorage object
         /** \param name filename */
         void read(const std::string name);
      
      private:
      
         //Contains the group number, the number of irreps, and the multiplication table
         Irreps SymmInfo;
         
         //Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep
         int * Isizes;
         
         /*The following conventions are used for storage:
            - 4-fold permutation symmetry: 2DM_ijkl = 2DM_jilk = 2DM_klij = 2DM_lkji
            - Abelian point group symmetry: I_i x I_j = I_k x I_l
            - Reorder indices untill irreps I_i <= {I_j , I_k , I_l}
            - Icenter = I_i x I_j = I_k x I_l
            - storage[Icenter][I_i][I_k - I_i][ i + dim_i * ( j + dim_j * ( k + dim_k * l ) ) ]
                   with   * I_i given
                          * I_j = I_i x I_center >= I_i
                          * I_k >= I_i given
                          * I_l = Icenter x I_k with I_l >= I_i   */
         long long *** storage;
         
         //Calculate the number of elements --> true means allocate the storage and false means deallocate the storage!
         long long calcNumberOfElements(const bool allocateStorage);
         
         //The number of TwoDMstorage elements
         long long arrayLength;
         
         //The actual TwoDM elements
         double * theElements;
         
         //Functions to get the correct pointer to memory
         long long getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;

   };
}

#endif

