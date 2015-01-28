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
         /** \param HamOrigIn The original Hamiltonian
             \param iHandlerIn The DMRGSCF indices */
         DMRGSCFVmatRotations(Hamiltonian * HamOrigIn, DMRGSCFindices * iHandlerIn);
         
         //! Destructor
         virtual ~DMRGSCFVmatRotations();
         
         //! Fill the rotated two-body matrix elements of HamDMRG, based on HamOrig and unitary. Do entire blocks at once.
         /** \param HamDMRG The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to VmatRotated.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4. */
         void fillVmatDMRG(Hamiltonian * HamDMRG, DMRGSCFunitary * unitary, double * mem1, double * mem2) const;
         
         //! Fill the rotated two-body matrix elements, based on HamOrig and unitary. Do entire blocks at once.
         /** \param VmatRotated The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to VmatRotated.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4. */
         void fillVmatRotated(FourIndex * VmatRotated, DMRGSCFunitary * unitary, double * mem1, double * mem2) const;
         
         //! Fill the rotated two-body matrix elements with max. two virtual indices, based on HamOrig and unitary. Do entire blocks at once.
         /** \param theRotatedTEI The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to theRotatedTEI.
             \param mem1 Work memory with at least the size max(linsize of irreps)^4.
             \param mem2 Work memory with at least the size max(linsize of irreps)^4. */
         void fillRotatedTEI(DMRGSCFintegrals * theRotatedTEI, DMRGSCFunitary * unitary, double * mem1, double * mem2) const;
         
         //! Fill the rotated two-body matrix elements of HamDMRG, based on HamOrig and unitary. Cut the blocks into chunks with linear size maxBlockSize.
         /** \param HamDMRG The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to Vmat(HamDMRG).
             \param mem1 Work memory with at least the size maxBlockSize^4.
             \param mem2 Work memory with at least the size maxBlockSize^4.
             \param mem3 Work memory with at least the size maxBlockSize^4.
             \param maxBlockSize Parameter which indicates the size of the work memories. */
         void fillVmatDMRGBlockWise(Hamiltonian * HamDMRG, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize) const;
         
         //! Fill the rotated two-body matrix elements, based on HamOrig and unitary. Cut the blocks into chunks with linear size maxBlockSize.
         /** \param VmatRotated The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to VmatRotated.
             \param mem1 Work memory with at least the size maxBlockSize^4.
             \param mem2 Work memory with at least the size maxBlockSize^4.
             \param mem3 Work memory with at least the size maxBlockSize^4.
             \param maxBlockSize Parameter which indicates the size of the work memories.
             \param cutCorners If false, all rotated two-body matrix elements are calculated. If true, at most two virtual indices are considered. */
         void fillVmatRotatedBlockWise(FourIndex * VmatRotated, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize, const bool cutCorners) const;
         
         //! Fill the rotated two-body matrix elements with max. two virtual indices, based on HamOrig and unitary. Cut the blocks into chunks with linear size maxBlockSize.
         /** \param theRotatedTEI The rotated two-body matrix elements are stored here.
             \param unitary The unitary matrix to rotate Vmat(HamOrig) to theRotatedTEI.
             \param mem1 Work memory with at least the size maxBlockSize^4.
             \param mem2 Work memory with at least the size maxBlockSize^4.
             \param mem3 Work memory with at least the size maxBlockSize^4.
             \param maxBlockSize Parameter which indicates the size of the work memories. */
         void fillRotatedTEIBlockWise(DMRGSCFintegrals * theRotatedTEI, DMRGSCFunitary * unitary, double * mem1, double * mem2, double * mem3, const int maxBlockSize) const;
         
      private:
      
         //The original Hamiltonian
         Hamiltonian * HamOrig;
         
         //The indices bookkeeper
         DMRGSCFindices * iHandler;
      
         //Symmetry information
         Irreps SymmInfo;
         
         //The number of irreps
         int numberOfIrreps;
         
   };
}

#endif
