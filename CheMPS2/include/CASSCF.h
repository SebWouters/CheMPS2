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
    
    \section intro Introduction
    
    In methods which use a FCI solver, this solver can be replaced by DMRG. DMRG allows for an efficient extraction of the 2-RDM [CAS1, CAS2]. The 2-RDM of the active space is required in the complete active space self-consistent field (CASSCF) method to compute the gradient and the Hessian with respect to orbital rotations [CAS3]. It is therefore natural to introduce a CASSCF variant with DMRG as active space solver, called DMRG-SCF [CAS2, CAS4, CAS5], which allows to treat static correlation in large active spaces. In CheMPS2, we have implemented the augmented Hessian Newton-Raphson DMRG-SCF method, with exact Hessian [CAS5, CAS6]. It can be called with the function CheMPS2::CASSCF::doCASSCFnewtonraphson.

    \section origalgo Augmented Hessian Newton-Raphson DMRG-SCF

    Orbital rotations are represented by orthogonal matrices \f$\mathbf{U}\f$. Since orbitals are defined up to a phase factor, one can always choose the orbital gauges so that \f$\mathbf{U}\f$ is special orthogonal: \f$\det(\mathbf{U})=+1\f$. \f$\mathbf{U}\f$ can be parametrized by a real-valued skew-symmetric matrix \f$\mathbf{X} = -\mathbf{X}^T\f$, see CheMPS2::DMRGSCFunitary. The CASSCF energy is a function of \f$\mathbf{U}\f$ and hence of \f$\mathbf{X}\f$. Up to second order, the energy is given by 
    \f[
    E(\vec{x}) = E(0) + \vec{x}^T \vec{g} + \frac{1}{2} \vec{x}^T \mathbf{H} \vec{x},
    \f]
    where \f$\vec{x}\f$ contains the independent variables of \f$\mathbf{X}\f$. The variables \f$\vec{x}\f$ parametrize an additional orbital rotation \f$\mathbf{U}_{add} = \exp(\mathbf{X}(\vec{x}))\f$ with respect to the current orbitals \f$\mathbf{U}(n)\f$. The vector \f$\vec{g}\f$ is the gradient, and the matrix \f$\mathbf{H}\f$ the Hessian, for orbital rotations [CAS3]. The minimum of \f$E(\vec{x})\f$ is found at \f$\vec{x} = - \mathbf{H}^{-1} \vec{g}\f$. The additional orbital rotation \f$\mathbf{U}_{add}\f$ then allows to update the current orbitals \f$\mathbf{U}(n)\f$ to the new orbitals
    \f[
    \mathbf{U}(n+1) = \mathbf{U}_{add} \mathbf{U}(n) = \exp(\mathbf{X}(\vec{x}(n))) \mathbf{U}(n).
    \f]
    This updating scheme is called the Newton-Raphson method [CAS3]. If the Hessian is positive definite, these updates are stable. For saddle points in the energy landscape, the Hessian has negative eigenvalues, and these updates can be unstable. It is therefore better to use the augmented Hessian Newton-Raphson method [CAS7]:
    \f[
    \left[ \begin{array}{cc} \mathbf{H} & \vec{g} \\ \vec{g}^T & 0 \end{array} \right] \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right] = \alpha \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right].
    \f]
    The eigenvector with smallest algebraic eigenvalue determines a stable update \f$\vec{x}\f$, as is well explained in Ref. [CAS7].

    \section diis Direct inversion of the iterative subspace (DIIS)

    When the update norm \f$\|\vec{x}\|_2\f$ is small enough, the convergence can be accelerated by the direct inversion of the iterative subspace (DIIS) [CAS5, CAS8, CAS9, CAS10]. For a given set of orbitals \f$\mathbf{U}(n)\f$, the update \f$\vec{x}(n)\f$ is calculated with the augmented Hessian Newton-Rapshon method. This update defines the next set of orbitals:
    \f[
    \mathbf{U}(n+1) = \mathbf{U}_{add} \mathbf{U}(n) = \exp(\mathbf{X}(\vec{x}(n))) \mathbf{U}(n).
    \f]
    In DIIS, the error vector \f$\vec{x}(n)\f$ and the state vector \f$\mathbf{Y}(n+1) = \log(\mathbf{U}(n+1))\f$ are added to a list. The error
    \f[
    e = \left\| \sum\limits_{i = 1}^{\kappa} c_i \vec{x}(n - \kappa + i) \right\|_2
    \f]
    is minimized under the constraint \f$\sum_{i} c_i = 1\f$. \f$\kappa\f$ is the size of the list memory, i.e. the number of retained vectors. The minimization of the error \f$e\f$ can be performed with Lagrangian calculus:
    \f[
    \left[ \begin{array}{cc} \mathbf{B} & \vec{1} \\ \vec{1}^T & 0 \end{array} \right] \left[ \begin{array}{c} \vec{c} \\ \lambda \end{array} \right] = \left[ \begin{array}{c} \vec{0} \\ 1 \end{array} \right],
    \f]
    where \f$2\lambda\f$ is the Lagrangian multiplier and
    \f[
    \left[\mathbf{B}\right]_{ij} = \vec{x}^T(n - \kappa + i) \vec{x}(n - \kappa + j) = \left[\mathbf{B}\right]_{ji}.
    \f]
    The new state vector is then defined as
    \f[
    \mathbf{Y}_{new} = \sum\limits_{i = 1}^{\kappa} c_i \mathbf{Y}(n+1 - \kappa + i).
    \f]
    The new state vector \f$\mathbf{Y}_{new}\f$ is calculated by the function CheMPS2::DIIS::calculateParam. The current orbitals are then set to \f$\mathbf{U}(n+1) = \exp(\mathbf{Y}_{new})\f$.
    
    \section biblio References

    [CAS1]  D. Zgid and M. Nooijen, Journal of Chemical Physics 128, 144115 (2008). http://dx.doi.org/10.1063/1.2883980 \n
    [CAS2]  D. Ghosh, J. Hachmann, T. Yanai and G.K.-L. Chan, Journal of Chemical Physics 128, 144117 (2008). http://dx.doi.org/10.1063/1.2883976 \n
    [CAS3]  P.E.M. Siegbahn, J. Almlof, A. Heiberg and B.O. Roos, Journal of Chemical Physics 74, 2384-2396 (1981). http://dx.doi.org/10.1063/1.441359 \n
    [CAS4]  D. Zgid and M. Nooijen, Journal of Chemical Physics 128, 144116 (2008). http://dx.doi.org/10.1063/1.2883981 \n
    [CAS5]  T. Yanai, Y. Kurashige, D. Ghosh and G.K.-L. Chan, International Journal of Quantum Chemistry 109, 2178-2190 (2009). http://dx.doi.org/10.1002/qua.22099 \n
    [CAS6]  S. Wouters, W. Poelmans, P.W. Ayers and D. Van Neck, Computer Physics Communications 185, 1501-1514 (2014). http://dx.doi.org/10.1016/j.cpc.2014.01.019 \n
    [CAS7]  A. Banerjee, N. Adams, J. Simons and R. Shepard, Journal of Physical Chemistry 89, 52-57 (1985). http://dx.doi.org/10.1021/j100247a015 \n
    [CAS8]  P. Pulay, Chemical Physics Letters 73, 393-398 (1980). http://dx.doi.org/10.1016/0009-2614(80)80396-4 \n
    [CAS9]  C.D. Sherrill, Programming DIIS, http://vergil.chemistry.gatech.edu/notes/diis/node3.html (2000). \n
    [CAS10] T. Rohwedder and R. Schneider, Journal of Mathematical Chemistry 49, 1889-1914 (2011). http://dx.doi.org/10.1007/s10910-011-9863-y \n
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
