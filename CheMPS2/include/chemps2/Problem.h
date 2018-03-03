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
         int gL() const{ return L; }
         
         //! Get the point group symmetry
         /** \return The point group symmetry */
         int gSy() const{ return Ham->getNGroup(); }
         
         //! Get an orbital irrep
         /** \param nOrb The orbital index
             \return The irrep of the orbital with index nOrb */
         int gIrrep(const int nOrb) const;
         
         //! Get twice the targeted spin
         /** \return Twice the targeted spin */
         int gTwoS() const{ return TwoS; }
         
         //! Get the targeted particle number
         /** \return The targeted particle number */
         int gN() const{ return N; }
         
         //! Get the targeted irrep
         /** \return The targeted irrep */
         int gIrrep() const{ return Irrep; }
         
         //! Get the constant part of the Hamiltonian
         /** \return The constant part of the Hamiltonian */
         double gEconst() const{ return Ham->getEconst(); }
         
         //! Get a specific interaction matrix element
         /** \param alpha The first index (0 <= alpha < L)
             \param beta The second index
             \param gamma The third index
             \param delta The fourth index
             \return \f$ h_{\alpha \beta ; \gamma \delta} = \left(\alpha \beta \mid V \mid \gamma \delta \right) + \frac{1}{N-1} \left( \left( \alpha \mid T \mid \gamma \right) \delta_{\beta \delta} + \delta_{\alpha \gamma} \left( \beta \mid T \mid \delta \right) \right) \f$ */
         double gMxElement(const int alpha, const int beta, const int gamma, const int delta) const;
         
         //! Set the matrix elements: Note that each time you create a DMRG object, they will be overwritten with the eightfold permutation symmetric Hamiltonian again!!!
         /** \param alpha The first index (0 <= alpha < L)
             \param beta The second index
             \param gamma The third index
             \param delta The fourth index
             \param value The value to set the matrix element to */
         void setMxElement(const int alpha, const int beta, const int gamma, const int delta, const double value);
         
         //! Construct a table with the h-matrix elements (two-body augmented with one-body). Remember to recall this function each time you change the Hamiltonian!
         void construct_mxelem();
         
         //! Check whether the given parameters L, N, and TwoS are not inconsistent and whether 0<=Irrep<nIrreps. A more thorough test will be done when the FCI virtual dimensions are constructed.
         /** \return True if consistent, else false */
         bool checkConsistency() const;
         
         //! Get whether the Hamiltonian orbitals are reordered for the DMRG calculation
         /** \return Whether the Hamiltonian orbitals are reordered */
         bool gReorder() const;
         
         //! Get the DMRG index corresponding to a Ham index
         /** \param HamOrb The Ham index
             \return The DMRG index */
         int gf1(const int HamOrb) const;
         
         //! Get the Ham index corresponding to a DMRG index
         /** \param DMRGOrb The DMRG index
             \return The Ham index */
         int gf2(const int DMRGOrb) const;
         
         //! Reorder the orbitals, so that they form irrep blocks, with order of irreps Ag B1u B3u B2g B2u B3g B1g Au. Previous reorderings are cleared.
         void SetupReorderD2h();
         
         //! Reorder the orbitals, so that they form irrep blocks, with order of irreps A1 B1 B2 A2. Previous reorderings are cleared.
         void SetupReorderC2v();
         
         //! Reorder the orbitals to a custom ordering. Previous reorderings are cleared.
         /** \param dmrg2ham Array which contains the reordering: dmrg2ham[ dmrg_lattice_site ] = hamiltonian_index. */
         void setup_reorder_custom(int * dmrg2ham);
         
         //! Reorder the orbitals for d(infinity)h. Previous reorderings are cleared.
         /** \param docc Array which contains for each irrep the number of doubly occupied orbitals
             \param sp_threshold Threshold to detect Delta_g and Delta_u partners based on single-particle energies */
         void setup_reorder_dinfh(int * docc, const double sp_threshold=1e-5);

         //! Check that ROHF-style occupancies are compatible with the currently targeted symmetry sector.
         /** \param occupancies Array which contains per DMRG orbital (not HAM orbital ordering!) the ROHF-style occupancy (0, 1 or 2) */
         bool check_rohf_occ( int * occupancies );

      private:
      
         //Pointer to the Hamiltonian --> constructed and destructed outside of this class
         const Hamiltonian * Ham;
         
         //The number of orbitals
         int L;
         
         //Twice the targeted spin
         int TwoS;
         
         //The targeted particle number
         int N;
         
         //The targeted irrep
         int Irrep;
         
         //Whether or not to reorder
         bool bReorder;
         
         //f1[HamiltonianIndex] = DMRGindex
         int * f1;
         
         //f2[DMRGIndex] = HamiltonianIndex
         int * f2;
         
         //Matrix element table
         double * mx_elem;
         
   };
}

#endif
