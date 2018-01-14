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

#ifndef CORRELATIONS_CHEMPS2_H
#define CORRELATIONS_CHEMPS2_H

#include "SyBookkeeper.h"
#include "Problem.h"
#include "TwoDM.h"
#include "TensorGYZ.h"
#include "TensorKM.h"

namespace CheMPS2{
/** Correlations class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 9, 2014
    
    \section introCorr Introduction
    
    Orbital correlation can be quantified by means of the two-orbital mutual information, as well as several correlation functions. The two-orbital mutual information was introduced in quantum chemistry by Rissler [COR1] in order to find the optimal orbital ordering for DMRG calculations. Later, the groups of Legeza and Reiher used this quantum information measure to study both the orbital ordering, as well as the molecular electronic structure [COR2, COR3, COR4, COR5]. The two-orbital mutual information can be interpreted in terms of generalized correlation functions, which reveal the detailed nature of the electron correlation [COR6].

    \section twoorbmutualinfo Two-orbital mutual information

    In quantum chemistry, molecular electronic structure is often interpreted in terms of single-particle states, called orbitals. For correlated calculations, the orbital occupations become entangled, which is in contrast to Hartree-Fock theory, where the orbitals have definite occupations. Several measures from quantum information theory allow to quantify orbital correlation. From the total Hilbert space, a subset \f$A\f$ of orbitals can be selected. Denote by \f$B\f$ the complement of \f$A\f$, i.e. the other orbitals which form together with \f$A\f$ the total Hilbert space. From a wavefunction in the total Hilbert space \f$\ket{\Psi}\f$, the (positive semidefinite) reduced density matrix of \f$A\f$ can be calculated as
    \f[
    \mathbf{\rho}_A = trace_{B} ~ \ket{\Psi}\bra{\Psi},
    \f]
    where \f$trace_{B}\f$ denotes the summation over the occupations of the orbitals in \f$B\f$. The (nonnegative) eigenspectrum of \f$\mathbf{\rho}_A\f$ directly reflects the quantum entanglement between the orbitals in \f$A\f$ and the ones in \f$B\f$. If there is only one nonzero eigenvalue, \f$A\f$ and \f$B\f$ are not entangled. The total wavefunction can then be factorized: \f$\ket{\Psi} = \ket{\Psi_A}\ket{\Psi_B}\f$. A measurement of an occupation in \f$B\f$ does not influence the outcome of a measurement in \f$A\f$. If there are several nonzero eigenvalues, \f$A\f$ and \f$B\f$ are entangled, and measurements in \f$A\f$ and \f$B\f$ are correlated. The total wavefunction can not be factorized in that case.

    The von Neumann entropy is a measure which quantifies the amount of entanglement between \f$A\f$ and \f$B\f$:
    \f[
    S_{A \mid B} = -trace_A ~ \mathbf{\rho}_A \log( \mathbf{\rho}_A ) = S_{B \mid A} \geq 0.
    \f]
    If \f$A\f$ and \f$B\f$ are unentangled, \f$S_{A \mid B} = 0\f$. With increasing entanglement between \f$A\f$ and \f$B\f$, \f$S_{A \mid B}\f$ becomes larger. When \f$A\f$ consists of a single orbital \f$i\f$, the von Neumann entropy is called the single-orbital entropy \f$S_1(i)\f$. When \f$A\f$ consists of two orbitals \f$i\f$ and \f$j\f$, the von Neumann entropy is called the two-orbital entropy \f$S_2(ij)\f$. The following inequality holds between them:
    \f[
    S_2(ij) \leq S_1(i) + S_1(j),
    \f]
    the so-called subadditivity property. Any entanglement between orbitals \f$i\f$ and \f$j\f$ reduces \f$S_2(ij)\f$ with respect to \f$S_1(i) + S_1(j)\f$. The amount of entanglement between two orbitals can then be quantified by means of the two-orbital mutual information
    \f[
    \left[ \mathbf{I} \right]_{ij} = \frac{1}{2} \left( S_1(i) + S_1(j) - S_2(ij) \right) \left( 1 - \delta_{ij} \right) = \left[ \mathbf{I} \right]_{ji} \geq 0.
    \f]
    With increasing entanglement between orbitals \f$i\f$ and \f$j\f$, \f$\left[ \mathbf{I} \right]_{ij}\f$ becomes larger. In DMRG, this measure can be calculated efficiently [COR1, COR5, COR6]. In a spin-adapted DMRG code, there are only 16 unique elements instead of the 26 which are described in Ref. [COR5].
    
    \section corrfunc Correlation functions

    In addition to the two-orbital mutual information, CheMPS2 also provides access to four other correlation functions. Three of them can be obtained from the reduced 2-RDMs \f$\Gamma^A\f$ and \f$\Gamma^B\f$, and from the reduced 1-RDM \f$\Gamma^1\f$:
    \f{eqnarray*}{
    \Gamma^A_{ij;kl} & = & \sum\limits_{\sigma \tau} \braket{ \hat{a}_{i\sigma}^{\dagger} \hat{a}_{j\tau}^{\dagger} \hat{a}_{l\tau} \hat{a}_{k\sigma} }, \\
    \Gamma^B_{ij;kl} & = & \sum\limits_{\sigma \tau} (-1)^{\sigma-\tau} \braket{ \hat{a}_{i\sigma}^{\dagger} \hat{a}_{j\tau}^{\dagger} \hat{a}_{l\tau} \hat{a}_{k\sigma} }, \\
    \Gamma^1_{i;j} & = & \sum\limits_{\sigma} \braket{ \hat{a}_{i\sigma}^{\dagger}  \hat{a}_{j\sigma} },
    \f}
    where \f$\sigma\f$ and \f$\tau\f$ can be \f$\pm\frac{1}{2}\f$. These correlation functions are the spin correlation function
    \f[
    C_{spin}(i,j) = 4 \left( \braket{\hat{S}_i^z \hat{S}_j^z} - \braket{\hat{S}_i^z} \braket{\hat{S}_j^z} \right) = \Gamma^B_{ij;ij} + \delta_{ij} \Gamma^1_{i;i},
    \f]
    (where \f$\braket{\hat{S}_i^z}=0\f$ due to the spin-adaptation), the spin-flip correlation function
    \f[
    C_{spinflip}(i,j) = \braket{\hat{S}_i^+ \hat{S}_j^-} + \braket{\hat{S}_i^- \hat{S}_j^+} = \frac{1}{2} \left( \Gamma^B_{ij;ji} - \Gamma^A_{ij;ji} \right) + \delta_{ij} \Gamma^1_{i;i},
    \f]
    and the density correlation function
    \f[
    C_{dens}(i,j) = \braket{\hat{n}_i \hat{n}_j} - \braket{\hat{n}_i} \braket{\hat{n}_j} = \Gamma^A_{ij;ij} + \Gamma^1_{i;i} \left( \delta_{ij} - \Gamma^1_{j;j} \right).
    \f]
    The fourth one is the singlet diradical correlation function [COR7]
    \f[
    C_{dirad}(i,j) = \braket{\hat{d}_{i\uparrow} \hat{d}_{j\downarrow}} + \braket{\hat{d}_{i\downarrow} \hat{d}_{j\uparrow}} - \braket{\hat{d}_{i\uparrow}} \braket{\hat{d}_{j\downarrow}} - \braket{\hat{d}_{i\downarrow}} \braket{\hat{d}_{j\uparrow}},
    \f]
    where \f$\hat{d}_{i\sigma} = \hat{n}_{i\sigma} (1 - \hat{n}_{i~-\sigma})\f$.
    
    The two-orbital mutual information and these correlation functions are calculated in the function CheMPS2::DMRG::calc2DMandCorrelations.

    \section biblioCorr References
    
    [COR1] J. Rissler, R. M. Noack and S.R. White, Chemical Physics 323, 519-531 (2006). http://dx.doi.org/10.1016/j.chemphys.2005.10.018 \n
    [COR2] G. Barcza, O. Legeza, K.H. Marti, M. Reiher, Physical Review A 83, 012508 (2011). http://dx.doi.org/10.1103/PhysRevA.83.012508 \n
    [COR3] K. Boguslawski, K.H. Marti, O. Legeza and M. Reiher, Journal of Chemical Theory and Computation 8, 1970-1982 (2012). http://dx.doi.org/10.1021/ct300211j \n
    [COR4] K. Boguslawski, P. Tecmer, O. Legeza and M. Reiher, Journal of Physical Chemistry Letters 3, 3129-3135 (2012). http://dx.doi.org/10.1021/jz301319v \n
    [COR5] K. Boguslawski, P. Tecmer, G. Barcza, O. Legeza and M. Reiher, Journal of Chemical Theory and Computation 9, 2959-2973 (2013). http://dx.doi.org/10.1021/ct400247p \n
    [COR6] G. Barcza, R.M. Noack, J. Solyom, O. Legeza, http://arxiv.org/abs/1406.6643 (2014). \n
    [COR7] J. Hachmann, J.J. Dorando, M. Aviles and G.K.-L. Chan, Journal of Chemical Physics 127, 134309 (2007). http://dx.doi.org/10.1063/1.2768362\n
    */
   class Correlations{

      public:
      
         //! Constructor
         /** \param denBKIn Symmetry sector bookkeeper
             \param ProbIn The problem to be solved
             \param the2DMin The 2-RDM of the active space */
         Correlations(const SyBookkeeper * denBKIn, const Problem * ProbIn, TwoDM * the2DMin);
         
         //! Destructor
         virtual ~Correlations();
         
         //! Get a Cspin term, using the DMRG indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCspin_DMRG(const int row, const int col) const;

         //! Get a Cspin term, using the HAM indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCspin_HAM(const int row, const int col) const;

         //! Get a Cdens term, using the DMRG indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCdens_DMRG(const int row, const int col) const;
         
         //! Get a Cdens term, using the HAM indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCdens_HAM(const int row, const int col) const;

         //! Get a Cspinflip term, using the DMRG indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCspinflip_DMRG(const int row, const int col) const;

         //! Get a Cspinflip term, using the HAM indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCspinflip_HAM(const int row, const int col) const;
         
         //! Get a Cdirad term, using the DMRG indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCdirad_DMRG(const int row, const int col) const;

         //! Get a Cdirad term, using the HAM indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getCdirad_HAM(const int row, const int col) const;

         //! Get a mutual information term, using the DMRG indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getMutualInformation_DMRG(const int row, const int col) const;

         //! Get a mutual information term, using the HAM indices
         /** \param row the first index
             \param col the second index
             \return the desired value */
         double getMutualInformation_HAM(const int row, const int col) const;
         
         //! Get the single-orbital entropy for a certain site, using the DMRG indices
         /** \param index The DMRG index
             \return The single-orbital entropy for this site */
         double SingleOrbitalEntropy_DMRG(const int index) const;
         
         //! Get the single-orbital entropy for a certain site, using hte HAM indices
         /** \param index The HAM index
             \return The single-orbital entropy for this site */
         double SingleOrbitalEntropy_HAM(const int index) const;
         
         //! Fill at the current step of the iterations the two-orbital mutual information and the remaining part of Cdirad
         /** \param denT DMRG site-matrices
             \param Gtensors Tensors required for the calculation
             \param Ytensors Tensors required for the calculation
             \param Ztensors Tensors required for the calculation
             \param Ktensors Tensors required for the calculation
             \param Mtensors Tensors required for the calculation*/
         void FillSite(TensorT * denT, TensorGYZ ** Gtensors, TensorGYZ ** Ytensors, TensorGYZ ** Ztensors, TensorKM ** Ktensors, TensorKM ** Mtensors);
         
         //! Return Idistance(power) (see return for definition)
         /** \param power The power used in Idistance
             \return \f$ Idistance(power) = sum_{ij} I(i,j) * \mid i-j \mid^{power} \f$ */
         double MutualInformationDistance(const double power) const;
         
         //! Print the correlation functions and two-orbital mutual information
         /** \param precision The number of digits to be printed
             \param columnsPerLine Rarara: The number of columns per line */
         void Print(const int precision=6, const int columnsPerLine=8) const;
         
         //! Broadcast the diradical correlation function and the two-orbital mutual information
         void mpi_broadcast();
         
      private:
      
         //The BK containing all the irrep information
         const SyBookkeeper * denBK;
         
         //The problem containing orbital reshuffling and symmetry information
         const Problem * Prob;
         
         //The 2-RDM of the active space
         TwoDM * the2DM;
         
         //The number of active space orbitals
         int L;
         
         //The spin correlation function
         double * Cspin;
         
         //The density correlation function
         double * Cdens;
         
         //The spin-flip correlation function
         double * Cspinflip;
         
         //The singlet diradical correlation function
         double * Cdirad;
         
         //The two-orbital mutual information
         double * MutInfo;
         
         //Helper function: fills Cspin, Cdens, Cspinflip, and Cdirad (the latter only partially)
         void FillSpinDensSpinflip();
         
         //Helper functions for FillSite
         double diagram1(TensorT * denT, TensorGYZ * denY, double * workmem) const;
         double diagram2(TensorT * denT, TensorGYZ * denZ, double * workmem) const;
         double diagram3(TensorT * denT, TensorGYZ * denG, double * workmem) const;
         double diagram4(TensorT * denT, TensorKM * denK, double * workmem) const;
         double diagram5(TensorT * denT, TensorKM * denM, double * workmem) const;
         
         //Helper function to print tables in a nice format
         void PrintTableNice(const double * table, const int sPrecision, const int columnsPerLine) const;
         
   };
}

#endif
