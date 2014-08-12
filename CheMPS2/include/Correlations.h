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

#ifndef CORRELATIONS_H
#define CORRELATIONS_H

#include <iostream>
#include <sstream>

#include "SyBookkeeper.h"
#include "Problem.h"
#include "TwoDM.h"
#include "TensorGYZ.h"
#include "TensorK.h"
#include "TensorM.h"

namespace CheMPS2{
/** Correlations class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 9, 2014
    
    The correlations class contains the following spin-summed correlation functions:
       - The spin correlation function \f$ C_{spin}(i,j) = 4 \left( \braket{ \hat{S}_i^z \hat{S}_j^z } - \braket{ \hat{S}_i^z } \braket{ \hat{S}_j^z } \right) \f$
       - The density correlation function \f$ C_{dens}(i,j) = \braket{ \hat{n}_i \hat{n}_j } - \braket{ \hat{n}_i } \braket{ \hat{n}_j } \f$
       - The spin-flip correlation function \f$ C_{spinflip}(i,j) = \braket{ \hat{S}_i^+ \hat{S}_j^- } + \braket{ \hat{S}_i^- \hat{S}_j^+ } \f$
    
    as well as the two-orbital mutual information:
    \f$ I(i,j) = \frac{1}{2} \left( S_1(i) + S_1(j) - S_2(ij) \right) \left( 1 - \delta_{ij} \right) \f$. \n
    For the latter see:
       - Rissler et al., Chem. Phys. 323, 519 (2006)
       - Boguslawski, Legeza et al., JCTC 9, 2959 (2013)
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
         
         //! Fill at the current step of the iterations the two-orbital mutual information
         /** \param denT DMRG site-matrices
             \param Gtensors Tensors required for the calculation
             \param Ytensors Tensors required for the calculation
             \param Ztensors Tensors required for the calculation
             \param Ktensors Tensors required for the calculation
             \param Mtensors Tensors required for the calculation*/
         void FillSite(TensorT * denT, TensorGYZ ** Gtensors, TensorGYZ ** Ytensors, TensorGYZ ** Ztensors, TensorK ** Ktensors, TensorM ** Mtensors);
         
         //! Return Idistance(power) (see return for definition)
         /** \param power The power used in Idistance
             \return \f$ Idistance(power) = sum_{ij} I(i,j) * \mid i-j \mid^{power} \f$ */
         double MutualInformationDistance(const double power) const;
         
         //! Print the correlation functions and two-orbital mutual information
         /** \param precision The number of digits to be printed
             \param columnsPerLine Rarara: The number of columns per line */
         void Print(const int precision=6, const int columnsPerLine=8) const;
         
      private:
      
         //The BK containing all the irrep information
         const SyBookkeeper * denBK;
         
         //The problem containing orbital reshuffling and symmetry information
         const Problem * Prob;
         
         //The 2-RDM of the active space
         TwoDM * the2DM;
         
         //The number of active space orbitals
         int L;
         
         //The one-body RDM
         double * OneRDM;
         
         //The spin correlation function
         double * Cspin;
         
         //The density correlation function
         double * Cdens;
         
         //The spin-flip correlation function
         double * Cspinflip;
         
         //The two-orbital mutual information
         double * MutInfo;
         
         //Helper function: fills OneRDM, Cspin, Cdens, Cspinflip, MutInfo
         void FillOneRDMSpinDensSpinflip();
         
         //Helper functions for FillSite
         double diagram1(TensorT * denT, TensorGYZ * denY, double * workmem) const;
         double diagram2(TensorT * denT, TensorGYZ * denZ, double * workmem) const;
         double diagram3(TensorT * denT, TensorGYZ * denG, double * workmem) const;
         double diagram4(TensorT * denT, TensorSwap * denK, double * workmem) const;
         double diagram5(TensorT * denT, TensorSwap * denM, double * workmem) const;
         
         //Helper function to print tables in a nice format
         void PrintTableNice(const double * table, const int sPrecision, const int columnsPerLine) const;
         
   };
}

#endif
