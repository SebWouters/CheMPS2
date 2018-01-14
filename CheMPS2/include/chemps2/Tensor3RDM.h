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

#ifndef TENSOR3RDM_CHEMPS2_H
#define TENSOR3RDM_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorOperator.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** Tensor3RDM class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 18, 2015
    
    The Tensor3RDM class contains the functions to calculate new renormalized operators for three second quantized operators. */
   class Tensor3RDM : public TensorOperator{

      public:
      
         //! Constructor
         /** \param boundary   The boundary index
             \param two_j1_in  Twice the spin-value of the intermediary coupling of the first two second quantized operators
             \param two_j2     Twice the spin-value of coupling of the three second quantized operators
             \param nelec      How many electrons there are more in the symmetry sector of the lower leg compared to the upper leg
             \param irrep      Direct product of irreps of the three second quantized operators
             \param prime_last Convention in which the tensor operator is stored (see class information)
             \param book       The symmetry sector partitioning */
         Tensor3RDM(const int boundary, const int two_j1_in, const int two_j2, const int nelec, const int irrep, const bool prime_last, const SyBookkeeper * book);
         
         //! Destructor
         virtual ~Tensor3RDM();
         
         //! Get the intermediary spin coupling value of the Tensor3RDM
         /** \return Twice J1, the intermediary spin coupling value of the Tensor3RDM */
         int get_two_j1() const;
         
         //! Get the spin value of the Tensor3RDM
         /** \return Twice J2, the spin value of the Tensor3RDM */
         int get_two_j2() const;
         
         //! Make diagram a1
         /** \param Sigma   The TensorS0 or TensorS1 to make diagram a1
             \param denT    The TensorT to make diagram a1
             \param workmem Work memory */
         void a1(TensorOperator * Sigma, TensorT * denT, double * workmem);
         
         //! Make diagram b1
         /** \param Sigma   The TensorS0 or TensorS1 to make diagram b1
             \param denT    The TensorT to make diagram b1
             \param workmem Work memory */
         void b1(TensorOperator * Sigma, TensorT * denT, double * workmem);
         
         //! Make diagram c1
         /** \param denF    The TensorF0 or TensorF1 to make diagram c1
             \param denT    The TensorT to make diagram c1
             \param workmem Work memory */
         void c1(TensorOperator * denF, TensorT * denT, double * workmem);
         
         //! Make diagram d1
         /** \param denF    The TensorF0 or TensorF1 to make diagram d1
             \param denT    The TensorT to make diagram d1
             \param workmem Work memory */
         void d1(TensorOperator * denF, TensorT * denT, double * workmem);
         
         //! Make diagram extra1
         /** \param denT The TensorT to make diagram extra1 */
         void extra1(TensorT * denT);
         
         //! Make diagram extra2
         /** \param denL    The TensorL to make diagram extra2
             \param denT    The TensorT to make diagram extra2
             \param workmem Work memory */
         void extra2(TensorL * denL, TensorT * denT, double * workmem);
         
         //! Make diagram extra3
         /** \param denL    The TensorL to make diagram extra3
             \param denT    The TensorT to make diagram extra3
             \param workmem Work memory */
         void extra3(TensorL * denL, TensorT * denT, double * workmem);
         
         //! Make diagram extra4
         /** \param denL    The TensorL to make diagram extra4
             \param denT    The TensorT to make diagram extra4
             \param workmem Work memory */
         void extra4(TensorL * denL, TensorT * denT, double * workmem);
         
         //! Make the in-product of two Tensor3RDMs
         /** \param buddy The second tensor
             \return The in-product */
         double contract( Tensor3RDM * buddy ) const;
         
         //! Get whether the tensor convention is prime last
         /** \return Whether the tensor convention is prime last */
         bool get_prime_last() const;
         
      private:
      
         int two_j1;
         
   };
}

#endif
