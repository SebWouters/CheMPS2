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
         /** \param iHandler The DMRGSCFindices which contain information on the occupied, active, and virtual spaces */
         DMRGSCFmatrix( const DMRGSCFindices * iHandler );

         //! Destructor
         virtual ~DMRGSCFmatrix();

         //! Clear the matrix
         void clear();

         //! Make this matrix the identity matrix
         void identity();

         //! Set an element
         /** \param irrep The irrep number of the indices
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \param val The value to which the element of the matrix should be set */
         void set( const int irrep, const int p, const int q, const double val );

         //! Get an element
         /** \param irrep The irrep number of the indices
             \param p The first index (within the symmetry block)
             \param q The second index (within the symmetry block)
             \return The requested matrix element */
         double get( const int irrep, const int p, const int q ) const;

         //! Get a matrix block
         /** \param irrep The irrep number of the matrix block
             \return Pointer to the requested block */
         double * getBlock( const int irrep );

         //! Get the RMS deviation with another DMRGSCFmatrix
         /** \param other DMRGSCFmatrix to compare to
             \return The RMS deviation */
         double rms_deviation( const DMRGSCFmatrix * other ) const;

         //! Write a DMRGSCFmatrix to disk
         /** \param filename Filename to store the DMRGSCFmatrix to
             \param idx      DMRGSCFindices which contain the number of orbitals per irrep
             \param storage  Where the entries of the DMRGSCFmatrix are currently stored */
         static void write( const string filename, const DMRGSCFindices * idx, double ** storage );

         //! Read the DMRGSCFmatrix from disk
         /** \param filename Filename where the DMRGSCFmatrix is stored
             \param n_irreps The number of irreps for the symmetry group under consideration
             \param storage  Where the entries of the DMRGSCFmatrix should be loaded to */
         static void read( const string filename, const int n_irreps, double ** storage );

         #ifdef CHEMPS2_MPI_COMPILATION
         //! Broadcast this matrix to all processes
         /** \param ROOT The process which should broadcast */
         void broadcast( const int ROOT );
         #endif

      protected:

         //! The information on the occupied, active, and virtual spaces
         const DMRGSCFindices * iHandler;

         //! The matrix entries
         double ** entries;

         //! The number of irreps
         int num_irreps;

   };
}

#endif

