/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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

#ifndef EDMISTONRUEDENBERG_H
#define EDMISTONRUEDENBERG_H

#include "Options.h"
#include "DMRGSCFunitary.h"
#include "Hamiltonian.h"

namespace CheMPS2{
/** EdmistonRuedenberg class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 14, 2014
    
    The EdmistonRuedenberg class localizes the orbitals of the active space based on the maximalization of \f$ O = \sum_i V_{ii;ii} \f$.
*/
   class EdmistonRuedenberg{

      public:
      
         //! Constructor
         /** \param HamIn The active space Hamiltonian
             \param printLevelIn If 0: nothing is printed. If 1: info at beginning and end. If >1: intermediary info as well. */
         EdmistonRuedenberg(Hamiltonian * HamIn, const int printLevelIn=1);
         
         //! Destructor
         virtual ~EdmistonRuedenberg();
         
         //! Maximize the Edmiston-Ruedenberg cost function
         /** \param temp1 Work memory of at least max(dim(irrep(Ham)))^4
             \param temp2 Work memory of at least max(dim(irrep(Ham)))^4
             \param gradThreshold Stop if the norm of the gradient is smaller than this value
             \param maxIter Stop if maxIter iterations have been performed (when gradThreshold is not reached)
             \return The value of the optimized cost function */
         double Optimize(double * temp1, double * temp2, const double gradThreshold=EDMISTONRUED_gradThreshold, const int maxIter=EDMISTONRUED_maxIter);
         
         //! Get the pointer to the unitary to use in DMRGSCF
         /** \return Pointer to the unitary which defines the localized orbitals */
         DMRGSCFunitary * getUnitary();
         
      private:
      
         //Pointer to the active space Hamiltonian
         Hamiltonian * Ham;
         
         //The print level
         int printLevel;
         
         //The symmetry information object
         Irreps SymmInfo;
         
         //The DMRGSCF index handler (in order to be able to recycle the DMRGSCFunitary object)
         DMRGSCFindices * iHandler;
         
         //DMRGSCF unitary
         DMRGSCFunitary * unitary;
         
         //The rotated two-body matrix elements
         FourIndex * VmatRotated;
         
         //Calculate the gradient, hessian and update
         double augmentedHessianNewtonRaphson(double * gradient, double * temp1, double * temp2) const;
         double calcGradientValue(const int irrep, const int p, const int q) const;
         double calcHessianValue( const int irrep, const int p, const int q, const int r, const int s) const;
         
         //Calculate the cost function
         double costFunction() const;
         
   };
}

#endif
