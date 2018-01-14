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

#ifndef DMRGSCFOPTIONS_CHEMPS2_H
#define DMRGSCFOPTIONS_CHEMPS2_H

#include "Options.h"

namespace CheMPS2{
/** DMRGSCFoptions class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 19, 2014
    
    The DMRGSCFoptions class contains the options to control the DMRGSCF iterations. The default values are taken from Options.h. \n
    
    For DIIS: \n
    (1)  DoDIIS (bool) : Whether or not to do DIIS \n
    (2)  DIISGradientBranch (double) : Start DIIS when the 2-norm of the update vector (i.e. NOT the gradient vector!) is smaller than this value \n
    (3)  NumDIISVecs (int) : Number of previous updates to keep during DIIS \n
    (4)  StoreDIIS (bool) : Whether or not to store the DIIS checkpoint file \n
    (5)  DIISStorageName (string) : The filename to store the DIIS checkpoint \n
    
    General DMRGSCF control: \n
    (6)  MaxIterations (int) : The maximum number of DMRGSCF iterations \n
    (7)  GradientThreshold (double) : Stop the DMRGSCF iterations when the gradient for orbital rotation has a 2-norm smaller than this value \n
    (8)  StoreUnitary (bool) : Whether or not to store the Orbital Rotation checkpoint file \n
    (9)  UnitaryStorageName (string) : The filename to store the Orbital Rotation checkpoint \n
    (10) StateAveraging (bool) : Whether to do state-averaged or state-specific DMRGSCF \n
    
    DMRG active space options: \n
    (11) WhichActiveSpace (int) : Determines which active space is used for the DMRG (FCI replacement) calculations. If 1: NO, sorted within each irrep by NOON. If 2: Localized Orbitals (Edmiston-Ruedenberg), sorted within each irrep by the exchange matrix (Fiedler vector). If 3: Not localized, but only sorted within each irrep by the Fiedler vector of the exchange matrix. If other value: No additional active space rotations (the ones from DMRGSCF are of course performed). \n
    (12) DumpCorrelations (bool) : Whether or not to print the correlation functions and two-orbital mutual information of the active space \n
    (13) StartLocRandom (bool) : When localized orbitals are used, it is sometimes beneficial to start the localization procedure from a random unitary. A specific example is the reduction of the d2h point group of graphene nanoribbons to the cs point group, in order to make use of locality in the DMRG calculations. Since molecular orbitals will still belong to the full point group d2h, a random unitary helps in constructing localized orbitals which belong to the cs point group.
*/
   class DMRGSCFoptions{

      public:
      
         //! Constructor
         DMRGSCFoptions();
         
         //! Destructor
         virtual ~DMRGSCFoptions();
         
         //! Get whether DIIS should be performed
         /** \return Whether DIIS should be performed */
         bool getDoDIIS() const;
         
         //! Get the threshold for when DIIS should start
         /** \return The threshold for the 2-norm of the update vector (NOT the gradient vector) for starting DIIS */
         double getDIISGradientBranch() const;
         
         //! Get the number of DIIS update vectors which should be kept
         /** \return The number of DIIS update vectors which should be kept */
         int getNumDIISVecs() const;
         
         //! Get whether the DIIS checkpoint should be stored to disk
         /** \return Whether the DIIS checkpoint should be stored to disk */
         bool getStoreDIIS() const;
         
         //! Get the DIIS checkpoint filename
         /** \return The filename for the DIIS checkpoint */
         string getDIISStorageName() const;
         
         //! Get the maximum number of DMRGSCF iterations
         /** \return The maximum number of DMRGSCF iterations */
         int getMaxIterations() const;
         
         //! Get the threshold for DMRGSCF convergence
         /** \return The threshold for the 2-norm of the gradient vector for DMRGSCF convergence */
         double getGradientThreshold() const;
         
         //! Get whether the Orbital Rotation checkpoint should be stored to disk
         /** \return Whether the Orbital Rotation checkpoint should be stored to disk */
         bool getStoreUnitary() const;
         
         //! Get the Orbital Rotation checkpoint filename
         /** \return The filename for the Orbital Rotation checkpoint */
         string getUnitaryStorageName() const;
         
         //! Get whether state-averaging or state-specific DMRGSCF should be performed
         /** \return Whether state-averaging DMRGSCF should be performed */
         bool getStateAveraging() const;
         
         //! Get which active space should be considered in the DMRG routine
         /** \return Which active space should be considered in the DMRG routine. If 1: NO, sorted within each irrep by NOON. If 2: Localized Orbitals (Edmiston-Ruedenberg), sorted within each irrep by the exchange matrix (Fiedler vector). If other value: No additional active space rotations (the ones from DMRGSCF are of course performed). */
         int getWhichActiveSpace() const;
         
         //! Get whether the correlations and two-orbital mutual information should be printed
         /** \return Whether the correlations and two-orbital mutual information should be printed */
         bool getDumpCorrelations() const;
         
         //! Get whether the localization procedure should start from a random unitary
         /** \return Whether the localization procedure should start from a random unitary */
         bool getStartLocRandom() const;

         //! Set whether DIIS should be performed
         /** \param DoDIIS_in Whether DIIS should be performed */
         void setDoDIIS(const bool DoDIIS_in);
         
         //! Set the threshold for when DIIS should start
         /** \param DIISGradientBranch_in The threshold for the 2-norm of the update vector (NOT the gradient vector) for starting DIIS */
         void setDIISGradientBranch(const double DIISGradientBranch_in);
         
         //! Set the number of DIIS update vectors which should be kept
         /** \param NumDIISVecs_in The number of DIIS update vectors which should be kept */
         void setNumDIISVecs(const int NumDIISVecs_in);
         
         //! Set whether the DIIS checkpoint should be stored to disk
         /** \param StoreDIIS_in Whether the DIIS checkpoint should be stored to disk */
         void setStoreDIIS(const bool StoreDIIS_in);
         
         //! Set the DIIS checkpoint filename
         /** \param DIISStorageName_in The filename for the DIIS checkpoint */
         void setDIISStorageName(const string DIISStorageName_in);
         
         //! Set the maximum number of DMRGSCF iterations
         /** \param MaxIterations_in The maximum number of DMRGSCF iterations */
         void setMaxIterations(const int MaxIterations_in);
         
         //! Set the threshold for DMRGSCF convergence
         /** \param GradientThreshold_in The threshold for the 2-norm of the gradient vector for DMRGSCF convergence */
         void setGradientThreshold(const double GradientThreshold_in);
         
         //! Set whether the Orbital Rotation checkpoint should be stored to disk
         /** \param StoreUnitary_in Whether the Orbital Rotation checkpoint should be stored to disk */
         void setStoreUnitary(const bool StoreUnitary_in);
         
         //! Set the Orbital Rotation checkpoint filename
         /** \param UnitaryStorageName_in The filename for the Orbital Rotation checkpoint */
         void setUnitaryStorageName(const string UnitaryStorageName_in);
         
         //! Set whether state-averaging or state-specific DMRGSCF should be performed
         /** \param StateAveraging_in Whether state-averaging DMRGSCF should be performed */
         void setStateAveraging(const bool StateAveraging_in);
         
         //! Set which active space should be considered in the DMRG routine
         /** \param WhichActiveSpace_in Which active space should be considered in the DMRG routine. If 1: NO, sorted within each irrep by NOON. If 2: Localized Orbitals (Edmiston-Ruedenberg), sorted within each irrep by the exchange matrix (Fiedler vector). If other value: No additional active space rotations (the ones from DMRGSCF are of course performed). */
         void setWhichActiveSpace(const int WhichActiveSpace_in);
         
         //! Set whether the correlations and two-orbital mutual information should be printed
         /** \param DumpCorrelations_in Whether the correlations and two-orbital mutual information should be printed */
         void setDumpCorrelations(const bool DumpCorrelations_in);
         
         //! Set whether the localization procedure should start from a random unitary
         /** \param StartLocRandom_in Whether the localization procedure should start from a random unitary */
         void setStartLocRandom(const bool StartLocRandom_in);
         
      private:
      
         //See class information
         bool   DoDIIS;
         double DIISGradientBranch;
         int    NumDIISVecs;
         bool   StoreDIIS;
         string DIISStorageName;
         
         int    MaxIterations;
         double GradientThreshold;
         bool   StoreUnitary;
         string UnitaryStorageName;
         bool   StateAveraging;
         
         int    WhichActiveSpace;
         bool   DumpCorrelations;
         bool   StartLocRandom;
         
   };
}

#endif
