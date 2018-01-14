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

#ifndef SYBOOKKEEPER_CHEMPS2_H
#define SYBOOKKEEPER_CHEMPS2_H

#include "Problem.h"

namespace CheMPS2{
/** SyBookkeeper class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 14, 2013

    The SyBookkeeper class keeps track of all the symmetry at the boundaries. This includes:
     - the FCI virtual dimensions per symmetry sector
     - an extra consistency check (next to the one in Problem.cpp) to check whether the desired symmetry in the Problem class is possible (non-zero FCI dimensions)
     - the current virtual dimensions per symmetry sector, so everyone can check the current dimensions here. */
   class SyBookkeeper{

      public:

         //! Constructor
         /** \param Prob The problem to be solved
             \param D    The initial number of reduced renormalized DMRG basis states */
         SyBookkeeper( const Problem * Prob, const int D );

         //! Copy constructor
         /** \param tocopy The SyBookkeeper to be copied */
         SyBookkeeper( const SyBookkeeper & tocopy );

         //! Destructor
         virtual ~SyBookkeeper();

         //! Get the problem
         /** \return The Problem of the SyBookkeeper */
         const Problem * gProb() const;

         //! Get the number of orbitals
         /** \return The number of orbitals */
         int gL() const;

         //! Get an orbital irrep
         /** \param orbital The DMRG orbital index (because the Problem class is asked)
             \return The irrep of the orbital */
         int gIrrep( const int orbital ) const;

         //! Get twice the targeted spin
         /** \return Twice the targeted spin */
         int gTwoS() const;

         //! Get the targeted particle number
         /** \return The targeted particle number */
         int gN() const;

         //! Get the targeted irrep
         /** \return The targeted irrep */
         int gIrrep() const;

         //! Get the total number of irreps
         /** \return The number of irreps */
         int getNumberOfIrreps() const;

         //! Get the min. possible particle number for a certain boundary
         /** \param boundary The boundary index ( from 0 to L ( included ) )
             \return Nmin[ bound ] */
         int gNmin( const int boundary ) const;

         //! Get the max. possible particle number for a certain boundary
         /** \param boundary The boundary index
             \return Nmax[ bound ] */
         int gNmax( const int boundary ) const;

         //! Get the minimum possible spin value for a certain boundary and particle number
         /** \param boundary The boundary index
             \param N The particle number
             \return The corresponding minimum spin value */
         int gTwoSmin( const int boundary, const int N ) const;

         //! Get the maximum possible spin value for a certain boundary and particle number
         /** \param boundary The boundary index
             \param N The particle number
             \return The corresponding maximum spin value */
         int gTwoSmax( const int boundary, const int N ) const;

         //! Get the FCI virtual dimensions ( bound by SYBK_dimensionCutoff )
         /** \param boundary The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param irrep The irrep
             \return The corresponding FCI virtual dimension */
         int gFCIdim( const int boundary, const int N, const int TwoS, const int irrep ) const;

         //! Get the current virtual dimensions
         /** \param boundary The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param irrep The irrep
             \return The corresponding current virtual dimension */
         int gCurrentDim( const int boundary, const int N, const int TwoS, const int irrep ) const;

         //! Get whether the desired symmetry sector is possible
         /** \return Whether the desired symmetry sector is possible */
         bool IsPossible() const;

         //! Get the current virtual dimensions
         /** \param boundary The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param irrep The irrep
             \param value The new dimension size */
         void SetDim( const int boundary, const int N, const int TwoS, const int irrep, const int value );

         //! Get the maximum virtual dimension at a certain boundary
         /** \param boundary The boundary index
             \return The maximum virtual dimension at the boundary */
         int gMaxDimAtBound( const int boundary ) const;

         //! Get the total reduced virtual dimension at a certain boundary
         /** \param boundary The boundary index
             \return The total reduced virtual dimension at the boundary */
         int gTotDimAtBound( const int boundary ) const;

         //! Restart by setting the virtual dimensions from boundary start to boundary stop to FCI virtual dimensions based on the environment
         /** \param start Start boundary index to create FCI virtual dimensions based on the environment
             \param stop  Stop  boundary index to create FCI virtual dimensions based on the environment
             \param virtual_dim The total virtual dimension to rescale the newly created symmetry sectors */
         void restart( const int start, const int stop, const int virtual_dim );

      private:

         // Pointer to the Problem --> constructed and destructed outside of this class
         const Problem * Prob;

         // The number of irreps
         int num_irreps;

         // Contains the minimum particle number possible at a boundary ( array length = L + 1 )
         int * Nmin;

         // Contains the maximum particle number possible at a boundary ( array length = L + 1 )
         int * Nmax;

         // Contains twice the minimum spin projection possible at a boundary & particle number ( array dimensions = ( L + 1 ) x ( Nmax[ boundary ] - Nmin[ boundary ] + 1 )
         int ** TwoSmin;

         // Contains twice the maximum spin projection possible at a boundary & particle number ( array dimensions = ( L + 1 ) x ( Nmax[ boundary ] - Nmin[ boundary ] + 1 )
         int ** TwoSmax;

         // FCI dimensions ( array access: FCIdim[ boundary ][ N - gNmin( boundary ) ][ ( TwoS - gTwoSmin( boundary, N ) ) / 2 ][ irrep ] )
         int **** FCIdim;

         // CUR dimensions ( array access: CURdim[ boundary ][ N - gNmin( boundary ) ][ ( TwoS - gTwoSmin( boundary, N ) ) / 2 ][ irrep ] )
         int **** CURdim;

         // Allocate the arrays
         void allocate_arrays();

         // Construct the FCI dimensions
         void fillFCIdim();
         void fill_fci_dim_right( int **** storage, const int start, const int stop );
         void fill_fci_dim_left(  int **** storage, const int start, const int stop );

         // Get a dimension of FCIdim / CURdim
         int gDimPrivate( int **** storage, const int boundary, const int N, const int TwoS, const int irrep ) const;

         // Scale CURdim with virtual_dim from boundary start to boundary stop ( both included )
         void ScaleCURdim( const int virtual_dim, const int start, const int stop );

         // Copy dimension arrays
         void CopyDim( int **** origin, int **** target );

   };
}

#endif
