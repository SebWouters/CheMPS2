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

#ifndef SYBOOKKEEPER_CHEMPS2_H
#define SYBOOKKEEPER_CHEMPS2_H

#include "Problem.h"
#include "Irreps.h"
#include "Options.h"

namespace CheMPS2{
/** SyBookkeeper class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 14, 2013
    
    The SyBookkeeper class keeps track of all the symmetry at the boundaries. This includes:
     - the FCI virtual dimensions per symmetry sector
     - an extra consistency check (next to the one in Problem.cpp) to check whether the desired symmetry in the Problem class is possible (non-zero FCI dimensions)
     - the current virtual dimensions per symmetry sector, so everyone can check the current dimensions here. */
   class SyBookkeeper : public Irreps{

      public:
      
         //! Constructor
         /** \param Probin The problem to be solved
             \param Din The initial number of reduced renormalized DMRG basis states */
         SyBookkeeper(const Problem * Probin, const int Din);
         
         //! Destructor
         virtual ~SyBookkeeper();
         
         //! Get the number of orbitals
         /** \return The number of orbitals */
         int gL() const;
         
         //! Get an orbital irrep
         /** \param nOrb The orbital index
             \return The irrep of the orbital with index nOrb */
         int gIrrep(const int nOrb) const;
         
         //! Get twice the targeted spin
         /** \return Twice the targeted spin */
         int gTwoS() const;
         
         //! Get the targeted particle number
         /** \return The targeted particle number */
         int gN() const;
         
         //! Get the targeted irrep
         /** \return The targeted irrep */
         int gIrrep() const;
         
         //! Get the min. possible particle number for a certain boundary
         /** \param bound The boundary index (from 0 to L (included))
             \return Nmin[bound] */
         int gNmin(const int bound) const;

         //! Get the max. possible particle number for a certain boundary
         /** \param bound The boundary index
             \return Nmax[bound] */
         int gNmax(const int bound) const;
         
         //! Get the min. possible spin value for a certain boundary and particle number
         /** \param bound The boundary index
             \param N The particle number
             \return TwoSmin[bound][N-Nmin[bound]] */
         int gTwoSmin(const int bound, const int N) const;

         //! Get the max. possible spin value for a certain boundary and particle number
         /** \param bound The boundary index
             \param N The particle number
             \return TwoSmax[bound][N-Nmin[bound]] */
         int gTwoSmax(const int bound, const int N) const;
         
         //! Get the FCI virtual dimensions (bound by cutoff)
         /** \param bound The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param Icnt The irrep
             \return FCIdim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] */
         int gFCIdim(const int bound, const int N, const int TwoS, const int Icnt) const;
         
         //! Get the total (reduced) virtual dimension
         /** \return The virtual dimension (reduced) */
         int gD() const;
         
         //! Get the current virtual dimensions
         /** \param bound The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param Icnt The irrep
             \return CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] */
         int gCurrentDim(const int bound, const int N, const int TwoS, const int Icnt) const;
         
         //! Get whether the desired symmetry sector is possible
         /** \return The virtual dimension (reduced) */
         bool IsPossible() const;
         
         //! Get the current virtual dimensions
         /** \param bound The boundary index
             \param N The particle number
             \param TwoS Twice the spin sector
             \param Icnt The irrep
             \param val The new dimension size */
         void SetDim(const int bound, const int N, const int TwoS, const int Icnt, const int val);
         
         //! Get the max. virtual dimension at a certain boundary. Useful function to preallocate memory when constructing Heff.
         /** \param iBound The boundary index
             \return The max. virtual dimension at iBound */
         int gMaxDimAtBound(const int iBound) const;
         
      private:
      
         //Pointer to the Problem --> constructed and destructed outside of this class
         const Problem * Prob;
         
         //Contains the min. particle number possible at a boundary: length of array = L+1
         int * Nmin;
         
         //Contains the max. particle number possible at a boundary: length of array = L+1
         int * Nmax;
         
         //Contains twice the min. spin projection possible at a boundary & particle number: length of array = (L+1) x (Nmax[bound] - Nmin[bound] + 1)
         int ** TwoSmin;
         
         //Contains twice the max. spin projection possible at a boundary & particle number: length of array = (L+1) x (Nmax[bound] - Nmin[bound] + 1)
         int ** TwoSmax;
         
         //FCI dimensions FCIdim[bound][N][TwoS][Irrep]
         int **** FCIdim;
         
         //Current dimensions CurrentDim[bound][N][TwoS][Irrep]
         int **** CurrentDim;
         
         //Internal helpers
         void fillFCIdim(); //Fill the FCIdim
         int gDimPrivate(int **** storage, const int bound, const int N, const int TwoS, const int Icnt) const; //--> same as gFCIdim
         void ScaleCurrentDim(const int virtualD); //Do a scaling reduction from FCIdim to CurrentDim, with gD() as bounds.
         void print() const;
         
         
   };
}

#endif
