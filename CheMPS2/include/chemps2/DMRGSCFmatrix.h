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

#ifndef DMRGSCFMATRIX_CHEMPS2_H
#define DMRGSCFMATRIX_CHEMPS2_H

#include "DMRGSCFindices.h"

namespace CheMPS2{
/** DMRGSCFmatrix class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date January 26, 2015
    
    Container class for DMRGSCF matrices which are blockdiagonal in the irreps.
*/
   class DMRGSCFmatrix{

      public:
      
         //! Constructor
         /** \param iHandler_in The DMRGSCFindices which contain information on the occupied, active, and virtual spaces */
         DMRGSCFmatrix(DMRGSCFindices * iHandler_in);
         
         //! Destructor
         virtual ~DMRGSCFmatrix();
         
         //! Clear the matrix
         void clear();
         
         //! Set an element
         /** \param irrep The irrep number of the indices
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \param val The value to which the element of the matrix should be set */
         void set(const int irrep, const int p, const int q, const double val);

         //! Get an element
         /** \param irrep The irrep number of the indices
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \return The requested matrix element */
         double get(const int irrep, const int p, const int q) const;
         
         //! Get a matrix block
         /** \param irrep The irrep number of the matrix block
             \return Pointer to the requested block */
         double * getBlock(const int irrep);
      
      private:
      
         // The information on the occupied, active, and virtual spaces
         DMRGSCFindices * iHandler;
         
         // The matrix entries
         double ** entries;

   };
}

#endif

