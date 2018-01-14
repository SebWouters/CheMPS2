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

#ifndef DMRGSCFINDICES_CHEMPS2_H
#define DMRGSCFINDICES_CHEMPS2_H

#include "Irreps.h"

namespace CheMPS2{
/** DMRGSCFindices class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 14, 2014
    
    The DMRGSCFindices class contains the index conversion conventions for the CASSCF class. This class assumes that for the given molecular point group, the L orbitals are ordered in increasing irrep number: I_1 <= I_2 <= I_3 <= ... <= I_L. The integer arrays NOCC, NDMRG and NVIRT contain the number of occupied, active and virtual orbitals per point group irrep.
*/
   class DMRGSCFindices{

      public:
      
         //! Constructor
         /** \param L The total number of orbitals
             \param Group The group number, which defines the number of irreps
             \param NOCCin The number of occupied orbitals per irrep
             \param NDMRGin The number of active orbitals per irrep
             \param NVIRTin The number of virtual orbitals per irrep */
         DMRGSCFindices(const int L, const int Group, int * NOCCin, int * NDMRGin, int * NVIRTin);
         
         //! Destructor
         virtual ~DMRGSCFindices();
         
         //! Get the number of orbitals
         /** \return The number of orbitals */
         int getL() const;
         
         //! Get the group number
         /** \return The group number */
         int getGroupNumber() const;
         
         //! Get the number of irreps
         /** \return The number of irreps */
         int getNirreps() const;
         
         //! Get the number of orbitals for an irrep
         /** \param irrep The irreducible representation
             \return The number of orbitals for irrep */
         int getNORB(const int irrep) const;
         
         //! Get the number of occupied orbitals for an irrep
         /** \param irrep The irreducible representation
             \return The number of occupied orbitals for irrep */
         int getNOCC(const int irrep) const;
         
         //! Get the number of active orbitals for an irrep
         /** \param irrep The irreducible representation
             \return The number of active orbitals for irrep */
         int getNDMRG(const int irrep) const;
         
         //! Get the number of virtual orbitals for an irrep
         /** \param irrep The irreducible representation
             \return The number of virtual orbitals for irrep */
         int getNVIRT(const int irrep) const;
         
         //! Get the cumulative number of active orbitals for an irrep
         /** \param irrep The irreducible representation
             \return The cumulative number of active orbitals in the irreps < irrep */
         int getDMRGcumulative(const int irrep) const;
         
         //! Get in the original Hamiltonian index the start orbital for the occupied orbitals with a certain irrep
         /** \param irrep The irreducible representation
             \return The start orbital in the original Hamiltonian order for the occupied orbitals of irrep */
         int getOrigNOCCstart(const int irrep) const;
         
         //! Get in the original Hamiltonian index the start orbital for the active orbitals with a certain irrep
         /** \param irrep The irreducible representation
             \return The start orbital in the original Hamiltonian order for the active orbitals of irrep */
         int getOrigNDMRGstart(const int irrep) const;
         
         //! Get in the original Hamiltonian index the start orbital for the virtual orbitals with a certain irrep
         /** \param irrep The irreducible representation
             \return The start orbital in the original Hamiltonian order for the virtual orbitals of irrep */
         int getOrigNVIRTstart(const int irrep) const;
         
         //! Get an array with the irreps of each DMRG orbital
         /** \return pointer to irrepOfEachDMRGorbital */
         int * getIrrepOfEachDMRGorbital();
         
         //! Get the irrep corresponding to a global orbital index
         /** \param index The global orbital index
             \return The irrep of the corresponding orbital */
         int getOrbitalIrrep(const int index) const;
         
         //! Get the total number of occupied orbitals
         /** \return The total number of occupied orbitals */
         int getNOCCsum() const;
         
         //! Get the maximum NORB
         /** \return The maximum NORB */
         int getNORBmax() const;
         
         //! Get the orbital rotation parameter space size
         /** \return The orbital rotation parameter space size */
         int getROTparamsize() const;
         
         //! Print my contents
         void Print() const;

         
      private:
         
         //Number of orbitals
         int L;
         
         //Irreps controller which contains the group number
         Irreps SymmInfo;
         
         //Number of irreps (follows from Group number)
         int Nirreps;
         
         //Number of orbitals per irrep
         int * NORB;
         
         //Number of occupied orbitals per irrep
         int * NOCC;
         
         //Number of active orbitals per irrep
         int * NDMRG;
         
         //Number of virtual orbitals per irrep
         int * NVIRT;
         
         //Cumulative number of orbitals per irrep
         int * NORBcumulative;
         
         //Cumulative number of active orbitals per irrep
         int * NDMRGcumulative;
         
         //Irrep of each DMRG orbital (to construct HamDMRG)
         int * irrepOfEachDMRGorbital;
         
         //Irrep of each global orbital
         int * irrepOfEachOrbital;
         
   };
}

#endif
