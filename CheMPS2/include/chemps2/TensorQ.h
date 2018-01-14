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

#ifndef TENSORQ_CHEMPS2_H
#define TENSORQ_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorOperator.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "SyBookkeeper.h"
#include "Problem.h"

namespace CheMPS2{
/** TensorQ class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 5, 2013
    
    The TensorQ class is a storage and manipulation class for the complementary operator of three contracted creators/annihilitors. */
   class TensorQ : public TensorOperator{

      public:

         //! Constructor
         /** \param boundary_index The boundary index
             \param Idiff The irrep of the one creator ( sandwiched if TensorL ; to sandwich if TensorQ )
             \param moving_right If true: sweep from left to right. If false: sweep from right to left
             \param denBK Symmetry bookkeeper of the problem at hand
             \param Prob Problem containing the matrix elements
             \param site The site on which the last crea/annih should work */
         TensorQ(const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * denBK, const Problem * Prob, const int site);

         //! Destructor
         virtual ~TensorQ();
         
         //! Add terms after update/clear without previous tensors
         /** \param denT TensorT to construct the Q-term without previous tensors */
         void AddTermSimple(TensorT * denT);
         
         //! Add terms after update/clear with previous TensorL's
         /** \param Ltensors The TensorL's to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsL(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2);
         
         //! Add terms after update/clear with previous A-tensors and B-tensors
         /** \param denA The A-tensor to construct the Q-term
             \param denB The B-tensor to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsAB(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2);
         
         //! Add terms after update/clear with previous C-tensors and D-tensors
         /** \param denC The C-tensor to construct the Q-term
             \param denD The D-tensor to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsCD(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2);
         
      private:
      
         //! Pointer to the problem (contains the matrix elements)
         const Problem * Prob;
         
         //! Site on which the last crea/annih works
         int site;
         
         //Internal stuff
         void AddTermSimpleRight(TensorT * denT);
         void AddTermSimpleLeft(TensorT * denT);
         void AddTermsLRight(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsLLeft(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsABRight(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsABLeft(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsCDRight(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsCDLeft(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2);
         
   };
}

#endif
