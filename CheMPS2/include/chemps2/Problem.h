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

#ifndef PROBLEM_CHEMPS2_H
#define PROBLEM_CHEMPS2_H

#include "Hamiltonian.h"

namespace CheMPS2{
/** Problem class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 14, 2013
    
    Setup of the problem that is fed to the DMRG class. It contains
     - the Hamiltonian
     - the targeted spin
     - the targeted particle number
     - the targeted wavefunction irrep */
   class Problem{

      public:
      
         //! Constructor
         /** \param Hamin Pointer to the Hamiltonian
             \param TwoSin Twice the targeted spin
             \param Nin The targeted particle number
             \param Irrepin The targeted irrep  */
         Problem(const Hamiltonian * Hamin, const int TwoSin, const int Nin, const int Irrepin);
         
         //! Destructor
         virtual ~Problem();
         
         //! Get the number of orbitals
         /** \return The number of orbitals */
         int gL() const;
         
         //! Get the point group symmetry
         /** \return The point group symmetry */
         int gSy() const;
         
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
         
         //! Get the constant part of the Hamiltonian
         /** \return The constant part of the Hamiltonian */
         double gEconst() const;
         
         //! Get a specific interaction matrix element
         /** \param alpha The first index (0 <= alpha < L)
             \param beta The second index
             \param gamma The third index
             \param delta The fourth index
             \return \f$ h_{\alpha \beta ; \gamma \delta} = \left(\alpha \beta \mid V \mid \gamma \delta \right) + \frac{1}{N-1} \left( \left( \alpha \mid T \mid \gamma \right) \delta_{\beta \delta} + \delta_{\alpha \gamma} \left( \beta \mid T \mid \delta \right) \right) \f$ */
         double gMxElement(const int alpha, const int beta, const int gamma, const int delta) const;
         
         //! Check whether the given parameters L, N, and TwoS are not inconsistent and whether 0<=Irrep<nIrreps. A more thorough test will be done when the FCI virtual dimensions are constructed.
         /** \return True if consistent, else false */
         bool checkConsistency() const;
         
         //! Get whether the Hamiltonian orbitals are reordered for the DMRG calculation
         /** \return Whether the Hamiltonian orbitals are reordered */
         bool gReorderD2h() const;
         
         //! Get the DMRG index corresponding to a Ham index
         /** \param HamOrb The Ham index
             \return The DMRG index */
         int gf1(const int HamOrb) const;
         
         //! Get the Ham index corresponding to a DMRG index
         /** \param DMRGOrb The DMRG index
             \return The Ham index */
         int gf2(const int DMRGOrb) const;
         
         //! Reorder the orbitals, so that they form irrep blocks, with order of irreps Ag B1u B3u B2g B2u B3g B1g Au
         void SetupReorderD2h();
         
      private:
      
         //Pointer to the Hamiltonian --> constructed and destructed outside of this class
         const Hamiltonian * Ham;
         
         //Twice the targeted spin
         int TwoS;
         
         //The targeted particle number
         int N;
         
         //1/(N-1)
         double OneOverNMinusOne;
         
         //The targeted irrep
         int Irrep;
         
         //Whether or not to reorder D2h
         bool bReorderD2h;
         
         //f1[HamiltonianIndex] = DMRGindex
         int * f1;
         
         //f2[DMRGIndex] = HamiltonianIndex
         int * f2;
         
   };
}

#endif
