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

#ifndef CASSCF_CHEMPS2_H
#define CASSCF_CHEMPS2_H

#include "Hamiltonian.h"
#include "Irreps.h"
#include "TwoDM.h"
#include "DMRG.h"
#include "Problem.h"
#include "Options.h"
#include "ConvergenceScheme.h"
#include "DMRGSCFindices.h"
#include "DMRGSCFunitary.h"
#include "DIIS.h"
#include "DMRGSCFoptions.h"

namespace CheMPS2{
/** CASSCF class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date June 18, 2013
    
    The CASSCF class provides an implementation of DMRGSCF [Zgid et al., J. Chem. Phys. 128, 144116 (2008)] with as DMRG solver the spin-adapted DMRG code CheMPS2. For the CASSCF update, the following implementations are included:\n
    - the quadratically convergent Newton-Raphson approach [Siegbahn et al., J. Chem. Phys. 74, 2384 (1981)]
    - the augmented Hessian Newton-Raphson approach [Lengsfield, J. Chem. Phys. 73, 382 (1980)]; which is ORCA's prefered choice
*/
   class CASSCF{

      public:
      
         //! Constructor
         /** \param filename The file containing the Psi output to start the DMRGSCF calculations. */
         CASSCF(const string filename);
         
         //! Constructor
         /** \param HamIn Hamiltonian containing the matrix elements of the Hamiltonian for which a CASSCF calculation is desired
             \param DOCCin Array containing the number of doubly occupied HF orbitals per irrep
             \param SOCCin Array containing the number of singly occupied HF orbitals per irrep */
         CASSCF(Hamiltonian * HamIn, int * DOCCin, int * SOCCin);
         
         //! Destructor
         virtual ~CASSCF();
         
         //! Get the number of irreps
         /** \return The number of irreps */
         int getNumberOfIrreps();
         
         //! Set the start of the CASSCF calculation
         /** \param NoccIn Array of length numberOfIrreps containing the number of double occupied HF orbitals per irrep for the CASSCF loop.
             \param NDMRGIn Array of length numberOfIrreps containing the number of active orbitals per irrep for the CASSCF loop.
             \param NvirtIn Array of length numberOfIrreps containing the number of empty orbitals per irrep for the CASSCF loop. */
         void setupStart(int * NoccIn, int * NDMRGIn, int * NvirtIn);
         
         //! Does the state-specific CASSCF cycle with the augmented Hessian Newton-Raphson method
         /** \param Nelectrons Total number of electrons in the system: occupied HF orbitals + active space
             \param TwoS Twice the targeted spin
             \param Irrep Desired wave-function irrep
             \param OptScheme The optimization scheme to run the inner DMRG loop
             \param rootNum Denotes the targeted state in state-specific CASSCF; 1 means ground state, 2 first excited state etc.
             \param theDMRGSCFoptions Contains the DMRGSCF options
             \return The converged DMRGSCF energy */
         double doCASSCFnewtonraphson(const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * theDMRGSCFoptions);
         
         //! CASSCF unitary rotation remove call
         void deleteStoredUnitary(const string filename=CheMPS2::DMRGSCF_unitaryStorageName){ unitary->deleteStoredUnitary(filename); }
         
         //! CASSCF DIIS vectors remove call
         void deleteStoredDIIS(const string filename=CheMPS2::DMRGSCF_DIISstorageName){ if (theDIIS!=NULL){ theDIIS->deleteStoredDIIS(filename); }}
         
      private:
      
         //Index convention handler
         DMRGSCFindices * iHandler;
         
         //Unitary matrix storage and manipulator
         DMRGSCFunitary * unitary;
         
         //DIIS object
         DIIS * theDIIS;
      
         //The original Hamiltonian
         Hamiltonian * HamOrig;
         
         //The rotated 2-body matrix elements
         FourIndex * VmatRotated;
         
         //Irreps controller
         Irreps SymmInfo;
         
         //The number of orbitals
         int L;
         
         //The numberOfIrreps;
         int numberOfIrreps;
         
         //Double occupations
         int * DOCC;
         
         //Single occupations
         int * SOCC;
         
         //Boolean whether or not setupStart has been called
         bool setupStartCalled;
         
         //Number of DMRG orbitals
         int nOrbDMRG;

         //Copy theDMRG2DM over to CASSCF::DMRG2DM
         void copy2DMover(TwoDM * theDMRG2DM);

         //Update the unitary, 2DM and 1DM with the given NO basis
         void rotate2DMand1DM(const int N, double * eigenvecs, double * work); 
         
         //Space for the DMRG 1DM
         double * DMRG1DM;

         //Space for the DMRG 2DM
         double * DMRG2DM;
         
         //Set the DMRG 1DM
         void setDMRG1DM(const int N);

         //The NO in terms of the active space orbitals are stored in the nOrbDMRG*nOrbDMRG array eigenvecs
         void calcNOON(double * eigenvecs, double * workmem);
         
         //Copy the localized orbitals over from unitary to the nOrbDMRG*nOrbDMRG array eigenvecs
         void fillLocalizedOrbitalRotations(DMRGSCFunitary * unitary, double * eigenvecs);
         
         //Helper function to fetch DOCC and SOCC from filename
         void allocateAndFillOCC(const string filename);
         
         //Helper function to copy the DOCC and SOCC arrays
         void allocateAndFillOCC(int * DOCCin, int * SOCCin);
         
         //Helper function to check HF
         void checkHF();
         
         //Fill Econst and Tmat of HamDMRG
         void fillConstAndTmatDMRG(Hamiltonian * HamDMRG) const;
         
         //Calculate the gradient, return function is the gradient 2-norm
         double calcGradient(double * gradient);
         
         //Calculate the hessian
         void calcHessian(double * hessian, const int rowjump);
         
         //Do Augmented Hessian form of NR. On return, gradient contains the RESCALED gradient.
         double augmentedHessianNR(double * gradient, double * updateNorm);
         
         //Fmat function as defined by Eq. (11) in the Siegbahn paper.
         double FmatHelper(const int index1, const int index2) const;
         double Fmat(const int index1, const int index2) const;
         double ** Fmatrix;
         void buildFmat();
         
         //The Coulomb and exchange interaction with the occupied and active electrons respectively
         double ** QmatrixOCC;
         double ** QmatrixACT;
         double ** QmatrixWORK;
         double ** OneBodyMatrixElements;
         void buildQmatrixOCC();
         void buildQmatrixACT();
         void buildOneBodyMatrixElements();
         void rotateOldToNew(double ** matrix);
         void constructCoulombAndExchangeMatrixInOrigIndices(double ** densityMatrix, double ** result);
         double QmatOCC(const int index1, const int index2) const;
         double QmatACT(const int index1, const int index2) const;
         double TmatRotated(const int index1, const int index2) const;
         
         //Wmat function as defined by Eq.(21b) in the Siegbahn paper.
         double Wmat(const int index1, const int index2, const int index3, const int index4) const;
         
         //Function to get the occupancies to obtain coefficients of certain Slater determinants for neutral C2. Important to figure out diatomic D(inf)h symmetries when calculating them in D2h symmetry. The function is not basis set and active space dependent (at least if no B2g, B3g, B2u and B3u orbitals are condensed).
         void PrintCoeff_C2(DMRG * theDMRG);
         
   };
}

#endif
