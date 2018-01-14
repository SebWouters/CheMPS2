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

#ifndef DMRGSCFINTEGRALS_CHEMPS2_H
#define DMRGSCFINTEGRALS_CHEMPS2_H

#include "DMRGSCFindices.h"

namespace CheMPS2{
/** DMRGSCFintegrals class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date January 27, 2015
    
    Container class for rotated DMRGSCF integrals. At most two virtual indices are needed to compute the gradient and Hessian, which allows for computational as well as memory savings. The integrals are stored in two objects:
    \f{eqnarray*}{
       Coulomb & : & ( c_1 c_2 | a_1 a_2 ) \\
      Exchange & : & ( c_1 v_1 | c_2 v_2 )
    \f}
    where \f$c\f$ denote core (occupied + active) orbitals, \f$v\f$ virtual orbitals, and \f$a\f$ both core and virtual orbitals. For the Coulomb integrals the fourfold permutation symmetry \f$c_1 \leq c_2\f$ and \f$a_1 \leq a_2\f$ is taken into account and for the exchange integrals the twofold permutation symmetry \f$c_1 \leq c_2\f$.
    
    Consider for example a molecule in C1 symmetry with 100 core (occupied + active) and 300 virtual orbitals. The FourIndex class needs to store \f$\frac{1}{8}(N_{core} + N_{virt})^4\f$ doubles or 25.6 GB. This class needs to store \f$\frac{1}{4} N_{core}^2 (N_{core}+N_{virt})^2 + \frac{1}{2} N_{core}^2 N_{virt}^2\f$ doubles or 6.8 GB. The eightfold permutation symmetry of the \f$( c_1 c_2 | c_3 c_4 )\f$ part of the Coulomb object is not used in this class. For the example, this would however only result in a reduction of \f$\frac{1}{8} N_{core}^4\f$ doubles or 0.1 GB.
*/
   class DMRGSCFintegrals{

      public:
      
         //! Constructor
         /** \param iHandler The DMRGSCFindices which contain information on the core (occupied + active) and virtual spaces */
         DMRGSCFintegrals(DMRGSCFindices * iHandler);
         
         //! Destructor
         virtual ~DMRGSCFintegrals();
         
         //! Set the storage objects to zero
         void clear();
         
         //! Set an element of the Coulomb object ( c1 c2 | a1 a2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Ia1 The irrep of the first core+virtual index
             \param Ia2 The irrep of the second core+virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param a1 The first core+virtual index (within the irrep block)
             \param a2 The second core+virtual index (within the irrep block)
             \param val The value to which the element of the object should be set */
         void set_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2, const double val);

         //! Add a double to an element of the Coulomb object ( c1 c2 | a1 a2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Ia1 The irrep of the first core+virtual index
             \param Ia2 The irrep of the second core+virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param a1 The first core+virtual index (within the irrep block)
             \param a2 The second core+virtual index (within the irrep block)
             \param val The value with which the element should be augmented */
         void add_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2, const double val);

         //! Get an element of the Coulomb object ( c1 c2 | a1 a2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Ia1 The irrep of the first core+virtual index
             \param Ia2 The irrep of the second core+virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param a1 The first core+virtual index (within the irrep block)
             \param a2 The second core+virtual index (within the irrep block)
             \return The desired value of the Coulomb object */
         double get_coulomb(const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2) const;

         //! Set an element of the Exchange object ( c1 v1 | c2 v2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Iv1 The irrep of the first virtual index
             \param Iv2 The irrep of the second virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param v1 The first virtual index (within the irrep block; counting starts at NCORE[ Iv1 ])
             \param v2 The second virtual index (within the irrep block; counting starts at NCORE[ Iv2 ])
             \param val The value to which the element of the object should be set */
         void set_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2, const double val);
         
         //! Add a double to an element of the Exchange object ( c1 v1 | c2 v2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Iv1 The irrep of the first virtual index
             \param Iv2 The irrep of the second virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param v1 The first virtual index (within the irrep block; counting starts at NCORE[ Iv1 ])
             \param v2 The second virtual index (within the irrep block; counting starts at NCORE[ Iv2 ])
             \param val The value with which the element should be augmented */
         void add_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2, const double val);
         
         //! Get an element of the Exchange object ( c1 v1 | c2 v2 )
         /** \param Ic1 The irrep of the first core index
             \param Ic2 The irrep of the second core index
             \param Iv1 The irrep of the first virtual index
             \param Iv2 The irrep of the second virtual index
             \param c1 The first core index (within the irrep block)
             \param c2 The second core index (within the irrep block)
             \param v1 The first virtual index (within the irrep block; counting starts at NCORE[ Iv1 ])
             \param v2 The second virtual index (within the irrep block; counting starts at NCORE[ Iv2 ])
             \return The desired value of the Exchange object */
         double get_exchange(const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2) const;
         
         //! Get a two-body matrix element with at most 2 virtual indices, using the FourIndex API
         /** \param I1 The irrep of the first index
             \param I2 The irrep of the second index
             \param I3 The irrep of the third index
             \param I4 The irrep of the fourth index
             \param index1 The first index (within the irrep block)
             \param index2 The second index (within the irrep block)
             \param index3 The third index (within the irrep block)
             \param index4 The fourth index (within the irrep block)
             \return The desired two-body matrix element */
         double FourIndexAPI(const int I1, const int I2, const int I3, const int I4, const int index1, const int index2, const int index3, const int index4) const;
      
      private:
      
         // The number of irreps
         int numberOfIrreps;
         
         // The number of core orbitals per irrep
         int * NCORE;
         
         // The number of virtual orbitals per irrep
         int * NVIRTUAL;
         
         // The total number of orbitals per irrep
         int * NTOTAL;
         
         // Pointers to the Coulomb array:
         //    - coulomb_ptr[ Ic1 x Ic2 != 0 ][ I_c1 ][ I_a1 ][ c1 + size(c1) * c2 ] points to the relevant part of coulomb_array
         //    - coulomb_ptr[ Ic1 x Ic2 == 0 ][ I_c1 ][ I_a1 ][ c1 + (c2*(c2+1))/2 ] points to the relevant part of coulomb_array
         long long **** coulomb_ptr;
         long long coulomb_size;
         double * coulomb_array;
         long long calcNumCoulombElements(const bool allocate);
         long long get_coulomb_ptr( const int Ic1, const int Ic2, const int Ia1, const int Ia2, const int c1, const int c2, const int a1, const int a2 ) const;
         
         // Pointers to the Exchange array:
         //    - exchange_ptr[ Ic1 x Ic2 != 0 ][ I_c1 ][ I_v1 ][ c1 + size(c1) * c2 ] points to the relevant part of exchange_array
         //    - exchange_ptr[ Ic1 x Ic2 == 0 ][ I_c1 ][ I_v1 ][ c1 + (c2*(c2+1))/2 ] points to the relevant part of exchange_array
         long long **** exchange_ptr;
         long long exchange_size;
         double * exchange_array;
         long long calcNumExchangeElements(const bool allocate);
         long long get_exchange_ptr( const int Ic1, const int Ic2, const int Iv1, const int Iv2, const int c1, const int c2, const int v1, const int v2 ) const;

   };
}

#endif

