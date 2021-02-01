/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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

#ifndef DMRGSCFROTATIONS_CHEMPS2_H
#define DMRGSCFROTATIONS_CHEMPS2_H

#include "Options.h"
#include "Hamiltonian.h"
#include "DMRGSCFunitary.h"
#include "DMRGSCFintegrals.h"

#include <hdf5.h>

namespace CheMPS2{
/** DMRGSCFrotations class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 14, 2014
    
    The DMRGSCFrotations class performs the two-body matrix element rotations for the DMRGSCF and Edmiston-Ruedenberg classes.
*/
   class DMRGSCFrotations{

      public:

         //! Fill the rotated two-body matrix elements for the space. If the blocks become too large, disk is used.
         /** \param ORIG_VMAT The FourIndex object with the original ERI.
             \param NEW_VMAT The FourIndex object where the new ERI should be stored.
             \param ROT_TEI The rotated two-body matrix elements are stored here.
             \param space1 Orbital space 1 (O, A, V, C, or F).
             \param space2 Orbital space 2 (O, A, V, C, or F).
             \param space3 Orbital space 3 (O, A, V, C, or F).
             \param space4 Orbital space 4 (O, A, V, C, or F).
             \param idx The DMRGSCF indices.
             \param umat The unitary matrix to rotate ORIG_VMAT to NEW_VMAT.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4.
             \param mem_size Sizes of the work memories.
             \param filename Where to store the temporary intermediate objects. */
         static void rotate( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space1, const char space2, const char space3, const char space4, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename );

      private:

         // Blockwise rotations
         static void blockwise_first(  double * origin, double * target, int orig1, int dim2, const int dim34, double * umat1, int new1, int lda1 );
         static void blockwise_second( double * origin, double * target, int dim1, int orig2, const int dim34, double * umat2, int new2, int lda2 );
         static void blockwise_third(  double * origin, double * target, const int dim12, int orig3, const int dim4, double * umat3, int new3, int lda3 );
         static void blockwise_fourth( double * origin, double * target, const int dim12, int dim3, int orig4, double * umat4, int new4, int lda4 );

         // Space sizes
         static int dimension( DMRGSCFindices * idx, const int irrep, const char space );
         static int      jump( DMRGSCFindices * idx, const int irrep, const char space );

         // Copy the required integrals from ORIG_VMAT to eri
         static void fetch( double * eri, const FourIndex * ORIG_VMAT, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack );

         // Unpack the second coulomb pair
         static void unpackage_second( double * mem1, double * mem2, const int SIZE, const int ORIG );

         // Pack the first coulomb pair
         static void package_first( double * mem1, double * mem2, const int NEW, const int PACKED, const int SIZE );

         // Copy the rotated integrals from eri to NEW_VMAT or ROT_TEI, depending on 'space'
         static void write( double * eri, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space1, const char space2, const char space3, const char space4, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack );

         // HDF5 file handling
         static void open_file( hid_t * file_id, hid_t * dspc_id, hid_t * dset_id, const int first, const int second, const string filename );
         static void write_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int first_write );
         static void  read_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int second_read );
         static void close_file( hid_t file_id, hid_t dspc_id, hid_t dset_id );

   };
}

#endif
