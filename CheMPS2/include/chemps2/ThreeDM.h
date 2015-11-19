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

#ifndef THREEDM_CHEMPS2_H
#define THREEDM_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "TensorS0.h"
#include "TensorS1.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** ThreeDM class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 16, 2015
    
    The ThreeDM class stores the spin-summed three-particle reduced density matrix (3-RDM) of a converged DMRG calculation: \n
    \f$ \Gamma_{ijk;lmn} = \sum_{\sigma \tau s} \braket{ a^{\dagger}_{i \sigma} a^{\dagger}_{j \tau} a^{\dagger}_{k s} a_{n s} a_{m \tau} a_{l \sigma}} \f$\n
    Because the wave-function belongs to a certain Abelian irrep, \f$ I_{i} \otimes I_{j} \otimes I_{k} = I_{l} \otimes I_{m} \otimes I_{n} \f$ must be valid before the corresponding element \f$ \Gamma_{ijk;lmn} \f$ is non-zero.
*/
   class ThreeDM{

      public:
      
         //! Constructor
         /** \param denBKIn Symmetry sector bookkeeper
             \param ProbIn The problem to be solved */
         ThreeDM(const SyBookkeeper * book_in, const Problem * prob_in);
         
         //! Destructor
         virtual ~ThreeDM();
         
         //! Get a 3-RDM term, using the DMRG indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \param cnt5 the fifth index
             \param cnt6 the sixth index
             \return the desired value */
         double get_dmrg_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6) const;
         
         //! Get a 3-RDM, using the HAM indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \param cnt5 the fifth index
             \param cnt6 the sixth index
             \return the desired value */
         double get_ham_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6) const;
         
         //! Fill the 3-RDM terms corresponding to site denT->gIndex()
         /** \param denT DMRG site-matrices
             \param Ltens Ltensors
             \param F0tens F0tensors
             \param F1tens F1tensors
             \param S0tens S0tensors
             \param S1tens S1tensors*/
         void fill_site(TensorT * denT, TensorL *** Ltens, TensorF0 **** F0tens, TensorF1 **** F1tens, TensorS0 **** S0tens, TensorS1 **** S1tens);
             
         //! Return the trace (should be N(N-1)(N-2))
         /** \return Trace of the 3-RDM */
         double trace() const;
         
         //! Save the 3-RDM to disk
         void save() const;
         
         //! Load the 3-RDM from disk
         void read();
         
         //! Add the 3-RDM elements of all MPI processes
         void mpi_allreduce();
         
      private:
      
         //The BK containing all the irrep information
         const SyBookkeeper * book;
         
         //The problem containing orbital reshuffling and symmetry information
         const Problem * prob;
         
         //The DMRG chain length
         int L;
         
         //The 3-RDM elements (L*L*L*L*L*L array of doubles)
         double * elements;
         
         //Set all twelve-fold permutation symmetric 3-RDM elements, using the DMRG indices.
         void set_dmrg_index(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const int cnt5, const int cnt6, const double value);
         
         //Helper functions
         static int trianglefunction(const int k, const int glob);
         static int phase(const int two_times_power){ return ((((two_times_power/2)%2)!=0)?-1:1); }
         
         //Partitioning 2-4-0
         double diagram1(TensorT * denT, TensorF0 * denF0, double * workmem) const;
         
         //Partitioning 0-4-2
         double diagram3(TensorT * denT, TensorF0 * denF0, double * workmem) const;
         
         //Partitioning 3-3-0 (d4 --> d9)
         
         //Partitioning 2-3-1
         double diagram10(TensorT * denT, TensorS0 * denS0, TensorL * denL, double * workmem, double * workmem2) const;
         double diagram11(TensorT * denT, TensorS1 * denS1, TensorL * denL, double * workmem, double * workmem2) const;
         double diagram12(TensorT * denT, TensorF0 * denF0, TensorL * denL, double * workmem, double * workmem2) const;
         double diagram13(TensorT * denT, TensorF1 * denF1, TensorL * denL, double * workmem, double * workmem2) const;
         double diagram14(TensorT * denT, TensorF0 * denF0, TensorL * denL, double * workmem, double * workmem2) const;
         double diagram15(TensorT * denT, TensorF1 * denF1, TensorL * denL, double * workmem, double * workmem2) const;
         
         //Partitioning 1-3-2
         double diagram16(TensorT * denT, TensorL * denL, TensorS0 * denS0, double * workmem, double * workmem2) const;
         double diagram17(TensorT * denT, TensorL * denL, TensorS1 * denS1, double * workmem, double * workmem2) const;
         double diagram18(TensorT * denT, TensorL * denL, TensorF0 * denF0, double * workmem, double * workmem2) const;
         double diagram19(TensorT * denT, TensorL * denL, TensorF1 * denF1, double * workmem, double * workmem2) const;
         double diagram20(TensorT * denT, TensorL * denL, TensorF0 * denF0, double * workmem, double * workmem2) const;
         double diagram21(TensorT * denT, TensorL * denL, TensorF1 * denF1, double * workmem, double * workmem2) const;
         
         //Partitioning 1-4-1 and 2-2-2
         double     diagram2_22_23(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2) const;
         double    diagram24_27_28(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2) const;
         void         diagram25_26(TensorT * denT, TensorS1       * left, TensorS1       * right, double * workmem, double * workmem2) const;
         double diagram29_30_31_32(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const;
         double       diagram33_39(TensorT * denT, TensorF0       * left, TensorF0       * right, double * workmem, double * workmem2, const bool do39) const;
         double    diagram34_37_38(TensorT * denT, TensorF1       * left, TensorF1       * right, double * workmem, double * workmem2) const;
         double       diagram35_41(TensorT * denT, TensorF0       * left, TensorF1       * right, double * workmem, double * workmem2, const bool do41) const;
         double       diagram36_42(TensorT * denT, TensorF1       * left, TensorF0       * right, double * workmem, double * workmem2, const bool do42) const;
         double    diagram40_43_44(TensorT * denT, TensorF1       * left, TensorF1       * right, double * workmem, double * workmem2) const;
         double diagram45_46_47_48(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const;
         double diagram49_50_51_52(TensorT * denT, TensorOperator * left, TensorOperator * right, double * workmem, double * workmem2, const bool doTR) const;
         
         //Partitioning 3-2-1 (d53 --> d89)
         
         //Partitioning 3-1-2 (d90 --> d189)
         
   };
}

#endif
