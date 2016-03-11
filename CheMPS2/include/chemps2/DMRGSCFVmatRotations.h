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

         // Perform the rotations 1,4,3,2 when it does not happen blockwise
         static void full_base( double * eri, double * work, double * umat1, int new1, int orig1,
                                                             double * umat2, int new2, int orig2,
                                                             double * umat3, int new3, int orig3,
                                                             double * umat4, int new4, int orig4 );

         //Blockwise rotations
         static void blockwise_first(  double * origin, double * target, int orig1, int dim2,  int dim3,  int dim4,  double * umat1, int new1, int lda1 );
         static void blockwise_second( double * origin, double * target, int dim1,  int orig2, int dim3,  int dim4,  double * umat2, int new2, int lda2 );
         static void blockwise_third(  double * origin, double * target, int dim1,  int dim2,  int orig3, int dim4,  double * umat3, int new3, int lda3 );
         static void blockwise_fourth( double * origin, double * target, int dim1,  int dim2,  int dim3,  int orig4, double * umat4, int new4, int lda4 );

         //Determine the block sizes
         static int blocksize( const int total_size, const int max_block_size, int * num_blocks );

   };
}

#endif
