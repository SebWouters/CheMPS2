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

#ifndef TENSORQ_CHEMPS2_H
#define TENSORQ_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorA.h"
#include "TensorB.h"
#include "TensorC.h"
#include "TensorD.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "TensorSwap.h"
#include "SyBookkeeper.h"
#include "Problem.h"

namespace CheMPS2{
/** TensorQ class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 5, 2013
    
    The TensorQ class is a storage and manipulation class for the complementary operator of three contracted creators/annihilitors. */
   class TensorQ : public TensorSwap{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param IdiffIn The irrep of the one creator ( sandwiched if TensorL ; to sandwich if TensorQ )
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn Symmetry bookkeeper of the problem at hand
             \param ProbIn Problem containing the matrix elements
             \param siteIn The site on which the last crea/annih should work */
         TensorQ(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn, const Problem * ProbIn, const int siteIn);
         
         //! Destructor
         virtual ~TensorQ();
         
         //! Clear storage in case an update is not possible
         void ClearStorage();
         
         //! Add terms after update/clear without previous tensors
         /** \param denT TensorT to construct the Q-term without previous tensors */
         void AddTermSimple(TensorT * denT);
         
         //! Add terms after update/clear with previous TensorL's
         /** \param Ltensors The TensorL's to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsL(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2);
         
         //! Add terms after update/clear with previous TensorA's and TensorB's
         /** \param denA The TensorA to construct the Q-term
             \param denB The TensorB to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsAB(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2);
         
         //! Add terms after update/clear with previous TensorC's, TensorD's, TensorF0's, TensorF1's
         /** \param denC The TensorC to construct the Q-term
             \param denD The TensorD to construct the Q-term
             \param denT TensorT to construct the Q-term with previous TensorL's
             \param workmem Work memory
             \param workmem2 Work memory */
         void AddTermsCD(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2);
         
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
         void AddTermsABRight(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsABLeft(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsCDRight(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2);
         void AddTermsCDLeft(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2);
         
   };
}

#endif
