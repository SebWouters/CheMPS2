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

#ifndef TENSORX_CHEMPS2_H
#define TENSORX_CHEMPS2_H

#include "Tensor.h"
#include "TensorL.h"
#include "TensorQ.h"
#include "TensorA.h"
#include "TensorC.h"
#include "TensorD.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "Problem.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorX class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date March 6, 2013
    
    The TensorX class is a storage and manipulation class for completely contracted Hamiltonian terms. */
   class TensorX : public TensorDiag{

      public:
      
         //! Constructor
         /** \param indexIn The boundary index
             \param movingRightIn If true: sweep from left to right. If false: sweep from right to left
             \param denBKIn The symmetry bookkeeper with symmetry sector virtual dimensions
             \param ProbIn The Problem containing the Hamiltonian matrix elements */
         TensorX(const int indexIn, const bool movingRightIn, const SyBookkeeper * denBKIn, const Problem * ProbIn);
         
         //! Destructor
         virtual ~TensorX();
         
         //! Clear and add the relevant terms to the TensorX
         /** \param denT TensorT from which the new TensorX should be made
             \param Ltensors Array with the TensorL's
             \param Xtensor The previous TensorX
             \param Qtensor The previous TensorQ
             \param Atensor The previous TensorA
             \param Ctensor The previous TensorC
             \param Dtensor The previous TensorD */
         void update(TensorT * denT, TensorL ** Ltensors, TensorX * Xtensor, TensorQ * Qtensor, TensorA * Atensor, TensorC * Ctensor, TensorD * Dtensor);
         
         //! Clear and add the relevant terms to the TensorX
         /** \param denT TensorT from which the new TensorX should be made */
         void update(TensorT * denT);
         
      private:
      
         //whether moving right or not
         bool movingRight;
         
         //Problem containing the matrix elements
         const Problem * Prob;
      
         //helper functions
         void makenewRight(const int ikappa, TensorT * denT);
         void makenewLeft(const int ikappa, TensorT * denT);
         void addTermQLRight(const int ikappa, TensorT * denT, TensorL ** Lprev, TensorQ * Qprev, double * workmemRR, double * workmemLR, double * workmemLL);
         void addTermQLLeft(const int ikappa, TensorT * denT, TensorL ** Lprev, TensorQ * Qprev, double * workmemLL, double * workmemLR, double * workmemRR);
         void addTermALeft(const int ikappa, TensorT * denT, TensorA * Aprev, double * workmemLR, double * workmemLL);
         void addTermARight(const int ikappa, TensorT * denT, TensorA * Aprev, double * workmemRR, double * workmemLR);
         void addTermCRight(const int ikappa, TensorT * denT, TensorC * denC, double * workmemLR);
         void addTermCLeft(const int ikappa, TensorT * denT, TensorC * denC, double * workmemLR);
         void addTermDRight(const int ikappa, TensorT * denT, TensorD * denD, double * workmemLR);
         void addTermDLeft(const int ikappa, TensorT * denT, TensorD * denD, double * workmemLR);
         
   };
}

#endif
