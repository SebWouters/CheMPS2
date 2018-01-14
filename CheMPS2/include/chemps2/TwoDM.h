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

#ifndef TWODM_CHEMPS2_H
#define TWODM_CHEMPS2_H

#include "TensorT.h"
#include "TensorL.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "TensorS0.h"
#include "TensorS1.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TwoDM class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date June 13, 2013
    
    The TwoDM class stores the result of a converged DMRG calculation. With the 2DM \n
    \f$ \Gamma_{(\alpha \sigma) (\beta \tau) ; (\gamma \sigma) (\delta \tau)} = \braket{ a^{\dagger}_{\alpha \sigma} a^{\dagger}_{\beta \tau} a_{\delta \tau} a_{\gamma \sigma} } \f$\n
    we can define two spin-reduced versions of interest:\n
    \f$ \Gamma^{2A}_{\alpha \beta ; \gamma \delta} = \sum_{\sigma \tau} \Gamma_{(\alpha \sigma) (\beta \tau) ; (\gamma \sigma) (\delta \tau)} \f$ \n
    \f$ \Gamma^{2B}_{\alpha \beta ; \gamma \delta} = \sum_{\sigma} \left( \Gamma_{(\alpha \sigma) (\beta \sigma) ; (\gamma \sigma) (\delta \sigma)} - \Gamma_{(\alpha \sigma) (\beta -\sigma) ; (\gamma \sigma) (\delta -\sigma)} \right) \f$. \n
    Because the wave-function belongs to a certain Abelian irrep, \f$ I_{\alpha} \otimes I_{\beta} \otimes I_{\gamma} \otimes I_{\delta} = I_{trivial} \f$ must be valid before the corresponding element \f$ \Gamma^{A,B}_{\alpha \beta ; \gamma \delta} \f$ is non-zero. \n
\n
    We can also define spin-densities in the spin-ensemble as:\n
    \f{eqnarray*}{
       \Gamma^{spin}_{ij} & = & \frac{3(1 - \delta_{S,0})}{(S+1)(2S+1)} \sum_{S^z} S^z \braket{ S S^z \mid a^{\dagger}_{i \uparrow} a_{j \uparrow} - a^{\dagger}_{i \downarrow} a_{j \downarrow} \mid S S^z } \\
                          & = & \frac{3(1 - \delta_{S,0})}{2(S+1)} \left[ ( 2 - N_{elec} ) \Gamma^1_{ij} - \sum_r \Gamma^{2A}_{ir;rj} - \sum_r \Gamma^{2B}_{ir;rj} \right]
    \f} \n
    The normalization factor is chosen so that \f$ trace( \Gamma^{spin} ) = 2 S \f$.
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
         
         //! Get a 1-RDM term, using the DMRG indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \return the desired value */
         double get1RDM_DMRG(const int cnt1, const int cnt2) const;
         
         //! Get a spin-density term, using the DMRG indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \return the desired value  */
         double spin_density_dmrg( const int cnt1, const int cnt2 ) const;
         
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
         
         //! Get a 1-RDM term, using the HAM indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \return the desired value */
         double get1RDM_HAM(const int cnt1, const int cnt2) const;
         
         //! Get a spin-density term, using the HAM indices
         /** \param cnt1 the first index
             \param cnt2 the second index
             \return the desired value  */
         double spin_density_ham( const int cnt1, const int cnt2 ) const;
         
         //! Fill the 2DM terms with as second site index denT->gIndex()
         /** \param denT DMRG site-matrices
             \param Ltens Ltensors
             \param F0tens F0tensors
             \param F1tens F1tensors
             \param S0tens S0tensors
             \param S1tens S1tensors*/
         void FillSite(TensorT * denT, TensorL *** Ltens, TensorF0 **** F0tens, TensorF1 **** F1tens, TensorS0 **** S0tens, TensorS1 **** S1tens);
         
         //! After the whole 2-RDM is filled, a prefactor for higher multiplicities should be applied
         void correct_higher_multiplicities();
             
         //! Return the double trace of 2DM-A (should be N(N-1))
         /** \return Double trace of 2DM-A */
         double trace() const;
         
         //! Calculate the energy based on the 2DM-A
         /** \return The energy calculated as 0.5*trace(2DM-A * Ham) */
         double energy() const;
         
         //! Print the natural orbital occupation numbers
         void print_noon() const;
         
         //! Save the TwoDMs to disk
         void save() const;
         
         //! Load the TwoDMs from disk
         void read();

         //! Save the 2-RDM-A to disk in Hamiltonian indices
         /** \param filename The filename to store the 2-RDM at */
         void save_HAM( const string filename ) const;

         //! Write the 2-RDM-A to a file
         /** param filename where to write the 2-RDM-A elements to */
         void write2DMAfile(const string filename) const;
         
         //! Add the 2-RDM elements of all MPI processes
         void mpi_allreduce();
         
      private:
      
         //The BK containing all the irrep information
         const SyBookkeeper * denBK;
         
         //The problem containing orbital reshuffling and symmetry information
         const Problem * Prob;
         
         //The chain length
         int L;
         
         //Two 2DM^{A,B} objects
         double * two_rdm_A;
         double * two_rdm_B;
         
         // Set 2DM terms, using the DMRG indices
         void set_2rdm_A_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value);
         void set_2rdm_B_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value);
         
         //Helper functions
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
