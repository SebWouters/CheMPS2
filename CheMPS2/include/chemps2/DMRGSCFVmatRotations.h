/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

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

#ifndef DMRGSCFVMATROTATIONS_CHEMPS2_H
#define DMRGSCFVMATROTATIONS_CHEMPS2_H

#include "Options.h"
#include "Hamiltonian.h"
#include "DMRGSCFunitary.h"
#include "DMRGSCFintegrals.h"
#include "MyHDF5.h"

namespace CheMPS2{
/** DMRGSCFVmatRotations class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 14, 2014
    
    The DMRGSCFVmatRotations class performs the two-body matrix element rotations for the DMRGSCF and Edmiston-Ruedenberg classes.
*/
   class DMRGSCFVmatRotations{

      public:

         //! Constructor
         DMRGSCFVmatRotations();

         //! Destructor
         virtual ~DMRGSCFVmatRotations();

         //! Fill the rotated two-body matrix elements for the space. If the blocks become too large, disk is used.
         /** \param ORIG_VMAT The FourIndex object with the original ERI.
             \param NEW_VMAT The FourIndex object where the new ERI should be stored.
             \param space Which orbital space NEW_VMAT corresponds to. Should be 'A' (active) or 'F' (full).
             \param idx The DMRGSCF indices.
             \param umat The unitary matrix to rotate ORIG_VMAT to NEW_VMAT.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4.
             \param mem_size Sizes of the work memories.
             \param fileame Where to store the temporary intermediate objects. */
         static void rotate( const FourIndex * ORIG_VMAT, FourIndex * NEW_VMAT, const char space, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename );

         //! Fill the rotated two-body matrix elements needed for CASSCF and CASPT2. If the blocks become too large, disk is used.
         /** \param ORIG_VMAT The FourIndex object with the original ERI.
             \param ROT_TEI The rotated two-body matrix elements are stored here.
             \param idx The DMRGSCF indices.
             \param umat The unitary matrix to rotate ORIG_VMAT to NEW_VMAT.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4.
             \param mem_size Sizes of the work memories.
             \param fileame Where to store the temporary intermediate objects. */
         static void rotate( const FourIndex * ORIG_VMAT, DMRGSCFintegrals * ROT_TEI, DMRGSCFindices * idx, DMRGSCFunitary * umat, double * mem1, double * mem2, const int mem_size, const string filename );

      private:

         // Blockwise rotations
         static void blockwise_first(  double * origin, double * target, int orig1, int dim2,  int dim3,  int dim4,  double * umat1, int new1, int lda1 );
         static void blockwise_second( double * origin, double * target, int dim1,  int orig2, int dim3,  int dim4,  double * umat2, int new2, int lda2 );
         static void blockwise_third(  double * origin, double * target, int dim1,  int dim2,  int orig3, int dim4,  double * umat3, int new3, int lda3 );
         static void blockwise_fourth( double * origin, double * target, int dim1,  int dim2,  int dim3,  int orig4, double * umat4, int new4, int lda4 );

         // Copy the required integrals from ORIG_VMAT to eri
         static void fetch( double * eri, const FourIndex * ORIG_VMAT, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack );
         
         // Copy the rotated integrals from eri to NEW_VMAT or ROT_TEI, depending on 'space'
         static void write( double * eri, FourIndex * NEW_VMAT, DMRGSCFintegrals * ROT_TEI, const char space, const int irrep1, const int irrep2, const int irrep3, const int irrep4, DMRGSCFindices * idx, const int start, const int stop, const bool pack );
         
         // HDF5 file handling
         static void open_file( hid_t * file_id, hid_t * dspc_id, hid_t * dset_id, const int first, const int second, const string filename );
         static void write_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int first_write );
         static void  read_file( hid_t dspc_id, hid_t dset_id, double * eri, const int start, const int size, const int second_read );
         static void close_file( hid_t file_id, hid_t dspc_id, hid_t dset_id );

   };
}

#endif
