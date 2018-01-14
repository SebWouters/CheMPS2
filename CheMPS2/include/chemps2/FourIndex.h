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

#ifndef FOURINDEX_CHEMPS2_H
#define FOURINDEX_CHEMPS2_H

#include "Irreps.h"

namespace CheMPS2{
/** FourIndex class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 8, 2013
    
    Container class for four-index tensors with Abelian point group symmetry (real character table; see Irreps.h): 2-particle matrix elements. The four-index tensor element V_ijkl has 8-fold permutation symmetry and is only nonzero when I_i x I_j = I_k x I_l. To clarify the convention, the potential energy is given: \n
    \f$\frac{1}{2} \sum\limits_{ijkl\sigma\tau} V_{ijkl} \delta_{I_i \otimes I_j \otimes I_k \otimes I_l, I_{trivial}} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{l \tau} \hat{a}_{k \sigma} \f$.\n
    Hence \f$ V_{ijkl} = ( ij \mid V \mid kl ) \f$ in physics notation, i.e. with i and k the same integration coordinate for the electron repulsion integrals.
*/
   class FourIndex{

      public:
      
         //! Constructor
         /** \param nGroup The symmetry group number (see Irreps.h)
             \param IrrepSizes Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep */
         FourIndex(const int nGroup, const int * IrrepSizes);
         
         //! Destructor
         virtual ~FourIndex();
         
         //! Set all two-body matrix elements to zero
         void Clear();
         
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
         
         //! Add a double to an element
         /** \param irrep_i The irrep number of the first orbital (see Irreps.h)
             \param irrep_j The irrep number of the second orbital
             \param irrep_k The irrep number of the third orbital
             \param irrep_l The irrep number of the fourth orbital
             \param i The first index (within the symmetry block)
             \param j The second index (within the symmetry block)
             \param k The third index (within the symmetry block)
             \param l The fourth index (within the symmetry block)
             \param val The value which should be added to the matrixelement */
         void add(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val);

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
         
         //! Get a given irrep size
         /** \param irrep The irrep for which you want to know the irrep size
             \return The corresponding irrep size */
         int get_irrep_size( const int irrep ) const;
         
         //! Save the FourIndex object
         /** \param name filename */
         void save(const std::string name) const;

         //! Load the FourIndex object
         /** \param name filename */
         void read(const std::string name);

         #ifdef CHEMPS2_MPI_COMPILATION
         //! Broadcast the ERI to all processes
         /** \param ROOT The process which should broadcast */
         void broadcast( const int ROOT );
         #endif

      private:
      
         //Contains the group number, the number of irreps, and the multiplication table
         Irreps SymmInfo;
         
         //Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep
         int * Isizes;
         
         /*The following conventions are used for storage:
            - 8-fold permutation symmetry: V_ijkl = V_jilk = V_ilkj = V_lijk = V_kjil = V_jkli = V_klij = V_lkji
            - Reorder indices until V_ijkl (I_i x I_j = I_k x I_l --> stored per center symm block):
                  - I_i <= I_j <= I_l and I_i <= I_k
                  - Icenter == Itriv : I_i == I_j and I_k == I_l
                  - Icenter >  Itriv : I_i <  I_j and I_k != I_l
                  - Vmat[Icenter][I_i][I_k] --> only created if ordering of all sectors is ok (I_i <= Icent x I_i ; I_k >= I_i ; Icent x I_k >= I_cent x I_i)
            - Once the order is established based on symmetry sectors, the order within symmetry sectors has to be set too:
                  - Case Icenter == Itrivial :
                        - If I_i == I_j == I_k == I_l : index (within symm block) i smallest; l>=j>=i and k>=i ; if i==j then l>=k
                        - Vmat[Itriv][I_i][I_i][i + k(k+1)/2]
                              - if i==j --> [j-i][l-k]
                              - if i< j --> [j-i][l-j]
                        - If I_i == I_j <  I_k == I_l : index i<=j ; if i<j then k and l fixed ; if i==j then l>=k
                        - Vmat[Itriv][I_i][I_k][i + nOrbWithIrrepI_i * k]
                              - if i==j --> [j-i][l-k]
                              - if i< j --> [j-i][l]
                  - Case Icenter > Itrivial (I_i < I_j and I_k != I_l) :
                        - If I_i == I_k and hence I_j == I_l : index k>=i and index l>=j
                        - Vmat[Icent][I_i][I_i][i + k*(k+1)/2][j][l-j]
                        - If I_i <  I_k and hence I_j <  I_l : fixed by block order
                        - Vmat[Icent][I_i][I_k][i + nOrbWithIrrepI_i * k][j][l] */
         long long ***** storage;
         
         //Calculate the number of unique FourIndex elements --> true means allocate the storage and false means deallocate the storage!
         long long calcNumberOfUniqueElements(const bool allocateStorage);
         
         //The number of unique FourIndex elements
         long long arrayLength;
         
         //The actual two-body matrix elements
         double * theElements;
         
         //Functions to get the correct pointer to memory
         long long getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
         long long getPtrIrrepOrderOK(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const;
         long long getPtrAllOK1(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;
         long long getPtrAllOK2(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;
         long long getPtrAllOK5(const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const;

   };
}

#endif

