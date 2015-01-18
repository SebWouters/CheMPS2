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

#ifndef TWODM_CHEMPS2_H
#define TWODM_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "TensorS0.h"
#include "TensorS1.h"
#include "SyBookkeeper.h"
#include "TwoDMstorage.h"

namespace CheMPS2{
/** TwoDM class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date June 13, 2013
    
    The TwoDM class stores the result of a converged DMRG calculation. With the 2DM \n
    \f$ \Gamma_{(\alpha \sigma) (\beta \tau) ; (\gamma \sigma) (\delta \tau)} = \braket{ a^{\dagger}_{\alpha \sigma} a^{\dagger}_{\beta \tau} a_{\delta \tau} a_{\gamma \sigma} } \f$\n
    we can define two spin-reduced versions of interest:\n
    \f$ \Gamma^A_{\alpha \beta ; \gamma \delta} = \sum_{\sigma \tau} \Gamma_{(\alpha \sigma) (\beta \tau) ; (\gamma \sigma) (\delta \tau)} \f$ \n
    \f$ \Gamma^B_{\alpha \beta ; \gamma \delta} = \sum_{\sigma} \left( \Gamma_{(\alpha \sigma) (\beta \sigma) ; (\gamma \sigma) (\delta \sigma)} - \Gamma_{(\alpha \sigma) (\beta -\sigma) ; (\gamma \sigma) (\delta -\sigma)} \right) \f$. \n
    Because the wave-function belongs to a certain Abelian irrep, \f$ I_{\alpha} \otimes I_{\beta} \otimes I_{\gamma} \otimes I_{\delta} = I_{trivial} \f$ must be valid before the corresponding element \f$ \Gamma^{A,B}_{\alpha \beta ; \gamma \delta} \f$ is non-zero.
*/
   class TwoDM{

      public:
      
         //! Constructor
         /** \param denBKIn Symmetry sector bookkeeper
             \param ProbIn The problem to be solved */
         TwoDM(const SyBookkeeper * denBKIn, const Problem * ProbIn);
         
         //! Destructor
         virtual ~TwoDM();
         
         //! Get a 2DM_A term, using the DMRG indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \return the desired value */
         double getTwoDMA_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const;
         
         //! Get a 2DM_B term, using the DMRG indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \return the desired value */
         double getTwoDMB_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const;
         
         //! Get a 2DM_A term, using the HAM indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \return the desired value */
         double getTwoDMA_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const;
         
         //! Get a 2DM_B term, using the HAM indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \param cnt3 the third index
             \param cnt4 the fourth index
             \return the desired value */
         double getTwoDMB_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const;
         
         //! Fill the 2DM terms with as second site index denT->gIndex()
         /** \param denT DMRG site-matrices
             \param Ltens Ltensors
             \param F0tens F0tensors
             \param F1tens F1tensors
             \param S0tens S0tensors
             \param S1tens S1tensors*/
         void FillSite(TensorT * denT, TensorL *** Ltens, TensorF0 **** F0tens, TensorF1 **** F1tens, TensorS0 **** S0tens, TensorS1 **** S1tens);
             
         //! Return the double trace of 2DM-A (should be N(N-1))
         /** \return Double trace of 2DM-A */
         double doubletrace2DMA();
         
         //! Calculate the energy based on the 2DM-A
         /** \return The energy calculated as 0.5*trace(2DM-A * Ham) */
         double calcEnergy();
         
         //! Save the TwoDMs to disk
         void save();
         
         //! Load the TwoDMs from disk
         void read();
         
      private:
      
         //The BK containing all the irrep information
         const SyBookkeeper * denBK;
         
         //The problem containing orbital reshuffling and symmetry information
         const Problem * Prob;
         
         //The chain length
         int L;
         
         //Two 2DM^{A,B} objects
         TwoDMstorage * TwoDMA;
         TwoDMstorage * TwoDMB;
         
         //number of orbitals per irrep
         int * irrep2num_orb;
         
         //index of an orbital within irrep block
         int * orb2IndexSy;
         
         // Set 2DM terms, using the DMRG indices
         void setTwoDMA_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value);
         void setTwoDMB_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value);
         
         //Helper functions
         static int trianglefunction(const int k, const int glob);
         double doD1(TensorT * denT);
         double doD2(TensorT * denT, TensorL * Lright, double * workmem);
         double doD3(TensorT * denT, TensorS0 * S0right, double * workmem);
         double doD4(TensorT * denT, TensorF0 * F0right, double * workmem);
         double doD5(TensorT * denT, TensorF0 * F0right, double * workmem);
         double doD6(TensorT * denT, TensorF1 * F1right, double * workmem);
         double doD7(TensorT * denT, TensorL * Lleft, double * workmem);
         double doD8(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g);
         void doD9andD10andD11(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, double * d9, double * d10, double * d11, int Irrep_g);
         double doD12(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g);
         double doD13(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g);
         double doD14(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g);
         double doD15(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g);
         double doD16(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g);
         double doD17orD21(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD17);
         double doD18orD22(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD18);
         double doD19orD23(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD19);
         double doD20orD24(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD20);
         
   };
}

#endif
