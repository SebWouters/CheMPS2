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

#ifndef MOLDEN_CHEMPS2_H
#define MOLDEN_CHEMPS2_H

#include <string>
#include "Irreps.h"

using std::string;

namespace CheMPS2{
/** Molden class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 24, 2016

    This class allows to rotate an molpro and psi4 molden files to a new molden file based on a CheMPS2 unitary. */
   class Molden{

      public:

         //! Constructor
         /** \param L The number of primitives
             \param group The psi4 group number
             \param irrep_sizes The number of orbitals per irrep */
         Molden( const int L, const int group, int * irrep_sizes );

         //! Destructor
         ~Molden();

         //! Read a molden file
         /** \param filename Filename of the original molden file */
         void read_molden( const string filename );

         //! Read a unitary matrix
         /** \param filename Filename of the unitary rotation */
         void read_unitary( const string filename );

         //! Multiply and print the new molden file
         /** \param original Filename of the original molden file
             \param output Filename of the new (rotated) molden file */
         void print( const string original, const string output );

      private:

         // Symmetry information
         Irreps SymmInfo;

         // Number of primitives
         int L;

         // Number of orbitals per irrep
         int num_irreps;

         // Irrep sizes
         int * Isizes;

         // molden[ irrep ][ gausnr + L * mo_orb ]
         double ** molden;

         // unitary[ irrep ][ cas_orb + Isizes[ irrep ] * mo_orb ]
         double ** unitary;

         // product[ irrep ][ gausnr + L * cas_orb ]
         double ** product;

   };
}

#endif
