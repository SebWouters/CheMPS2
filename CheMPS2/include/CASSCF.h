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
#include "DMRGSCFwtilde.h"
#include "DMRGSCFmatrix.h"
#include "DMRGSCFintegrals.h"

namespace CheMPS2{
/** CASSCF class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date June 18, 2013
    
    \section intro Introduction
    
    In methods which use a FCI solver, this solver can be replaced by DMRG. DMRG allows for an efficient extraction of the 2-RDM [CAS1, CAS2]. The 2-RDM of the active space is required in the complete active space self-consistent field (CASSCF) method to compute the gradient and the Hessian with respect to orbital rotations [CAS3]. It is therefore natural to introduce a CASSCF variant with DMRG as active space solver, called DMRG-SCF [CAS2, CAS4, CAS5], which allows to treat static correlation in large active spaces. In CheMPS2, I have implemented the augmented Hessian Newton-Raphson DMRG-SCF method, with exact Hessian [CAS5, CAS6]. It can be called with the function CheMPS2::CASSCF::doCASSCFnewtonraphson.
    
    \section equations_2step_dmrgscf Orbital gradient and Hessian

    The calculation of the orbital gradient and Hessian for DMRG-SCF is based on [CAS3]. The basic idea is to express the energy with the unitary group generators:
    \f{eqnarray*}{
      \hat{E}_{pq} & = & \sum\limits_{\sigma} \hat{a}^{\dagger}_{p \sigma} \hat{a}_{q \sigma} \\
      \left[ \hat{E}_{pq} , \hat{E}_{rs} \right] & = & \delta_{qr} \hat{E}_{ps} - \delta_{ps} \hat{E}_{rq} \\
      \hat{E}^{-}_{pq} & = & \hat{E}_{pq} - \hat{E}_{qp} \\
      \hat{T} & = & \sum\limits_{p<q} x_{pq} \hat{E}^{-}_{pq} \\
      E(\vec{x}) & = & \braket{0 | e^{\hat{T}} \hat{H} e^{-\hat{T}} | 0 } \\
      \left. \frac{\partial E(\vec{x})}{\partial x_{ij}} \right|_{0} & = & \braket{ 0 | \left[ \hat{E}_{ij}^{-}, \hat{H} \right] | 0 } \\
      \left. \frac{\partial^2 E(\vec{x})}{\partial x_{ij} \partial x_{kl}} \right|_{0} & = & \frac{1}{2} \braket{ 0 |  \left[ \hat{E}_{ij}^{-}, \left[ \hat{E}_{kl}^{-}, \hat{H} \right] \right] | 0 } + \frac{1}{2} \braket{ 0 |  \left[ \hat{E}_{kl}^{-}, \left[ \hat{E}_{ij}^{-}, \hat{H} \right] \right] | 0 }
    \f}
    The variables \f$x_{pq}\f$ only connect orbitals with the same irrep (\f$I_p=I_q\f$). Assuming that DMRG is exact, \f$x_{pq}\f$ in addition only connects orbitals when they belong to different occupation blocks: occupied, active, virtual. With some algebra, the derivatives can be rewritten. Real-valued symmetric one-electron integrals \f$h_{ij}\f$ and real-valued eightfold permutation symmetric two-electron integrals \f$(ij | kl)\f$ are assumed (chemical notation for the latter).
    \f{eqnarray*}{
      \Gamma^{2A}_{ijkl} & = & \sum\limits_{\sigma \tau} \braket{ 0 | \hat{a}^{\dagger}_{i \sigma} \hat{a}_{j \tau}^{\dagger} \hat{a}_{l \tau} \hat{a}_{k \sigma} | 0} \\
      \Gamma^1_{ij} & = & \sum\limits_{\sigma} \braket{ 0 | \hat{a}^{\dagger}_{i \sigma} \hat{a}_{j \sigma} | 0} \\
      \left. \frac{\partial E(\vec{x})}{\partial x_{ij}} \right|_{0} & = & 2 \left( F_{ij} - F_{ji} \right) \\
      F_{pq} & = & \sum\limits_{r} \Gamma_{pr}^{1} h_{qr} + \sum\limits_{rst} \Gamma^{2A}_{psrt} (qr | st) \\
      \left. \frac{\partial^2 E(\vec{x})}{\partial x_{ij} \partial x_{kl}} \right|_{0} & = & w_{ijkl} - w_{jikl} - w_{ijlk} + w_{jilk} \\
      w_{pqrs} & = & \delta_{qr} \left( F_{ps} + F_{sp} \right) + \tilde{w}_{pqrs} \\
      \tilde{w}_{pqrs} & = & 2 \Gamma^1_{pr} h_{qs} + 2 \sum\limits_{\alpha \beta} \left( \Gamma^{2A}_{r \alpha p \beta} (qs | \alpha \beta ) + \left( \Gamma^{2A}_{r \alpha \beta p} + \Gamma^{2A}_{r p \beta \alpha} \right) (q \alpha | s \beta ) \right)
    \f}
    In the calculation of \f$F_{pq}\f$, the indices \f$prst\f$ can only be occupied or active due to their appearance in the density matrices, and the only index which can be virtual is hence \f$q\f$. Moreover, due to the irrep symmetry of the integrals and density matrices, \f$F_{pq}\f$ is diagonal in the irreps: \f$I_p = I_q\f$. Alternatively, this can be understood by the fact that \f$x_{pq}\f$ only connects orbitals with the same irrep.

    In the calculation of \f$\tilde{w}_{pqrs}\f$, the indices \f$pr\alpha\beta\f$ can only be occupied or active due to their appearance in the density matrices, and the only indices which can be virtual are hence \f$qs\f$. Together with the remark for \f$F_{pq}\f$, this can save time for the two-electron integral rotation. Moreover, as \f$x_{pq}\f$ only connects orbitals with the same irrep, \f$I_p = I_q\f$ and \f$I_r = I_s\f$ in \f$\tilde{w}_{pqrs}\f$.

    By rewriting the density matrices, the calculation of \f$F_{pq}\f$ and \f$\tilde{w}_{pqrs}\f$ can be simplified. In the following, occ and act denote the doubly occupied and active orbital spaces, respectively.
    \f{eqnarray*}{
      \Gamma^1_{ij} & = & 2 \delta_{ij}^{occ} + \Gamma^{1,act}_{ij} \\
      \Gamma^{2A}_{ijkl} & = & \Gamma^{2A,act}_{ijkl} + 2 \delta_{ik}^{occ} \Gamma^{1,act}_{jl} - \delta_{il}^{occ} \Gamma^{1,act}_{jk} \\
      & + & 2 \delta_{jl}^{occ} \Gamma^{1,act}_{ik} - \delta_{jk}^{occ} \Gamma^{1,act}_{il} + 4 \delta_{ik}^{occ} \delta_{jl}^{occ} - 2 \delta_{il}^{occ} \delta_{jk}^{occ}
    \f}
    Define the following symmetric charge (Coulomb + exchange) matrices:
    \f{eqnarray*}{
      Q^{occ}_{ij} & = & \sum\limits_{s \in occ} \left[ 2 (ij | ss) - (is | js ) \right] \\
      Q^{act}_{ij} & = & \sum\limits_{st \in act} \frac{1}{2} \Gamma^{1,act}_{st} \left[ 2 (ij | st) - (is | jt ) \right] 
    \f}
    They can be calculated efficiently by (1) rotating the occupied and active density matrices from the current basis to the original basis, (2) contracting the rotated density matrices with the two-electron integrals in the original basis, and (3) rotating these contractions to the current basis. The constant part and the one-electron integrals of the active space Hamiltonian are:
    \f{eqnarray*}{
      \tilde{E}_{const} & = & E_{const} + \sum\limits_{s \in occ} \left( 2 h_{ss} + Q_{ss}^{occ} \right) \\
      \tilde{h}_{ij} & = & h_{ij} + Q_{ij}^{occ} 
    \f}
    The calculation of \f$F_{pq}\f$ boils down to:
    \f{eqnarray*}{
      p \in occ & : & F_{pq} = 2 \left( h_{qp} + Q^{occ}_{qp} + Q^{act}_{qp} \right) \\
      p \in act & : & F_{pq} = \sum\limits_{r \in act} \Gamma^{1,act}_{pr} \left[ h_{qr} + Q^{occ}_{qr} \right] +  \sum\limits_{rst \in act} \Gamma^{2A,act}_{psrt} (qr | st)
    \f}
    And the calculation of \f$\tilde{w}_{pqrs}\f$ (remember that \f$I_p = I_q\f$ and \f$I_r = I_s\f$):
    \f{eqnarray*}{
      (p,r) & \in & (occ,occ) : \tilde{w}_{pqrs} \\
                          & = & 4 \delta_{pr}^{occ} \left[ h_{qs} + Q^{occ}_{qs} + Q^{act}_{qs} \right] \\
                          & + & 4 \left[ 4 (qp | sr) - ( qs | pr ) - ( qr | sp ) \right] \\
      (p,r) & \in & (act,act) : \tilde{w}_{pqrs} = 2 \Gamma^{1,act}_{rp} \left[ h_{qs} + Q^{occ}_{qs} \right] \\
                          & + & 2 \sum\limits_{\alpha\beta \in act} \left[ \Gamma^{2A,act}_{r \alpha p \beta} (qs | \alpha \beta ) + \left( \Gamma^{2A,act}_{r \alpha \beta p} + \Gamma^{2A,act}_{r p \beta \alpha} \right) (q \alpha | s \beta ) \right] \\
      (p,r) & \in & (act,occ) : \tilde{w}_{pqrs} \\
                          & = & 2 \sum\limits_{\alpha \in act} \Gamma^{1,act}_{\alpha p} \left[ 4 (q \alpha | s r) - (qs | \alpha r) - (qr | s \alpha) \right] \\
      (p,r) & \in & (occ,act) : \tilde{w}_{pqrs} \\
                          & = & 2 \sum\limits_{\beta \in act} \Gamma^{1,act}_{r \beta} \left[ 4 (q p | s \beta) -  (qs | p \beta) - (q \beta | sp) \right]
    \f}

    \section origalgo Augmented Hessian Newton-Raphson DMRG-SCF
    
    The CASSCF energy is a function of \f$\vec{x}\f$. Up to second order, the energy is given by 
    \f[
    E(\vec{x}) = \braket{0 | e^{\hat{T}(\vec{x})} \hat{H} e^{-\hat{T}(\vec{x})} | 0 } \approx E(0) + \vec{x}^T \vec{g} + \frac{1}{2} \vec{x}^T \mathbf{H} \vec{x}.
    \f]
    The vector \f$\vec{g}\f$ is the gradient and the matrix \f$\mathbf{H}\f$ the Hessian for orbital rotations [CAS3]. They have been described in the previous section. The minimum of \f$E(\vec{x})\f$ is found at \f$\vec{x} = - \mathbf{H}^{-1} \vec{g}\f$. The variables \f$\vec{x}\f$ parametrize an additional orbital rotation \f$\mathbf{U}_{add} = \exp(\mathbf{X}(\vec{x}))\f$, with \f$\mathbf{X}(\vec{x}) = -\mathbf{X}^T(\vec{x})\f$ a real-valued skew-symmetric matrix. The additional orbital rotation \f$\mathbf{U}_{add}\f$ transforms the current orbitals \f$\mathbf{U}(n)\f$ to the new orbitals
    \f[
    \mathbf{U}(n+1) = \mathbf{U}_{add} \mathbf{U}(n) = \exp(\mathbf{X}(\vec{x}(n))) \mathbf{U}(n).
    \f]
    This updating scheme is called the Newton-Raphson method [CAS3]. If the Hessian is positive definite, these updates are stable. For saddle points in the energy landscape, the Hessian has negative eigenvalues, and these updates can be unstable. It is therefore better to use the augmented Hessian Newton-Raphson method [CAS7]:
    \f[
    \left[ \begin{array}{cc} \mathbf{H} & \vec{g} \\ \vec{g}^T & 0 \end{array} \right] \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right] = \alpha \left[ \begin{array}{c} \vec{x} \\ 1 \end{array} \right].
    \f]
    The eigenvector with smallest algebraic eigenvalue determines a stable update \f$\vec{x}\f$, as is well explained in Ref. [CAS7].
    
    As a final remark in this section, I would like to say that orbitals have gauge freedom. One can always multiply them with a phase factor. It is therefore possible to choose the orbital gauges so that all \f$\mathbf{U}\f$ are always special orthogonal: \f$\det(\mathbf{U})=+1\f$.

    \section diis Direct inversion of the iterative subspace (DIIS)

    When the update norm \f$\|\vec{x}\|_2\f$ is small enough, the convergence can be accelerated by the direct inversion of the iterative subspace (DIIS) [CAS5, CAS8, CAS9, CAS10]. For a given set of orbitals \f$\mathbf{U}(n)\f$, the update \f$\vec{x}(n)\f$ is calculated with the augmented Hessian Newton-Raphson method. This update defines the next set of orbitals:
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
         /** \param filename The file containing the Psi4 Hamiltonian text dump (mointegrals.so_PRINT) output to start the DMRGSCF calculations. */
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
         
         //! Build the F-matrix (Eq. (11) in the Siegbahn paper [CAS3])
         /** \param localFmat Matrix where the result should be stored
             \param localTmat Matrix which contains the one-electron integrals
             \param localJKocc Matrix which contains the Coulomb and exchange interaction due to the frozen core orbitals
             \param localJKact Matrix which contains the Coulomb and exchange interaction due to the active space
             \param localIdx Orbital index bookkeeper for the CASSCF calculations
             \param theInts The rotated two-electron integrals (at most 2 virtual indices)
             \param local2DM The DMRG 2-RDM
             \param local1DM The DMRG 1-RDM */
         static void buildFmat(DMRGSCFmatrix * localFmat, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM);
         
         //! Build the Wtilde-matrix (Eq. (20b) in the Siegbahn paper [CAS3])
         /** \param localwtilde Where the result should be stored
             \param localTmat Matrix which contains the one-electron integrals
             \param localJKocc Matrix which contains the Coulomb and exchange interaction due to the frozen core orbitals
             \param localJKact Matrix which contains the Coulomb and exchange interaction due to the active space
             \param localIdx Orbital index bookkeeper for the CASSCF calculations
             \param theInts The rotated two-electron integrals (at most 2 virtual indices)
             \param local2DM The DMRG 2-RDM
             \param local1DM The DMRG 1-RDM */
         static void buildWtilde(DMRGSCFwtilde * localwtilde, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM);
         
         //! Calculate the augmented Hessian Newton-Raphson update for the orthogonal orbital rotation matrix
         /** \param localFmat Matrix which contains the Fock operator (Eq. (11) in the Siegbahn paper [CAS3])
             \param localwtilde Object which contains the second order derivative of the energy with respect to the unitary (Eq. (20b) in the Siegbahn paper [CAS3])
             \param localIdx Orbital index bookkeeper for the CASSCF calculations
             \param localUmat The unitary matrix for CASSCF calculations (in this function it is used to fetch the orbital ordering convention of the skew-symmetric parametrization)
             \param theupdate Where the augmented Hessian Newton-Raphson update will be stored
             \param updateNorm Pointer to one double to store the update norm
             \param gradNorm Pointer to one double to store the gradient norm */
         static void augmentedHessianNR(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * theupdate, double * updateNorm, double * gradNorm);
         
         //! Copy over the DMRG 2-RDM
         /** \param theDMRG2DM The 2-RDM from the DMRG object
             \param totOrbDMRG The total number of DMRG orbitals
             \param localDMRG2DM The CASSCF 2-RDM */
         static void copy2DMover(TwoDM * theDMRG2DM, const int totOrbDMRG, double * localDMRG2DM);
         
         //! Construct the 1-RDM from the 2-RDM
         /** \param nDMRGelectrons The number of DMRG active space electrons
             \param totOrbDMRG The total number of DMRG orbitals
             \param localDMRG1DM The CASSCF 1-RDM
             \param localDMRG2DM The CASSCF 2-RDM */
         static void setDMRG1DM(const int nDMRGelectrons, const int totOrbDMRG, double * localDMRG1DM, double * localDMRG2DM);
         
         //! Calculate the natural orbitals and their occupation numbers
         /** \param localIdx Object which handles the index conventions for CASSCF
             \param eigenvecs Where to store the natural orbitals
             \param workmem Work memory
             \param localDMRG1DM The CASSCF 1-RDM */
         static void calcNOON(CheMPS2::DMRGSCFindices * localIdx, double * eigenvecs, double * workmem, double * localDMRG1DM);
         
         //! Rotate the CASSCF 1-RDM and 2-RDM to a new basis (to calculate the gradient and Hessian)
         /** \param nDMRGelectrons The number of DMRG active space electrons
             \param totOrbDMRG The total number of DMRG orbitals
             \param eigenvecs Where the eigenvectors are stored
             \param work Work memory
             \param localDMRG1DM The CASSCF 1-RDM
             \param localDMRG2DM The CASSCF 2-RDM */
         static void rotate2DMand1DM(const int nDMRGelectrons, int totOrbDMRG, double * eigenvecs, double * work, double * localDMRG1DM, double * localDMRG2DM);
         
         //! From an Edmiston-Ruedenberg active space rotation, fetch the eigenvectors and store them in eigenvecs
         /** \param unitary The Edmiston-Ruedenberg active space rotation
             \param localIdx Object which handles the index conventions for CASSCF
             \param eigenvecs Where the eigenvectors are stored */
         static void fillLocalizedOrbitalRotations(CheMPS2::DMRGSCFunitary * unitary, CheMPS2::DMRGSCFindices * localIdx, double * eigenvecs);
         
      private:
      
         //Index convention handler
         DMRGSCFindices * iHandler;
         
         //Unitary matrix storage and manipulator
         DMRGSCFunitary * unitary;
         
         //DIIS object
         DIIS * theDIIS;
      
         //The original Hamiltonian
         Hamiltonian * HamOrig;
         bool shouldHamOrigBeDeleted;
         
         //The rotated 2-body matrix elements with at most two virtual indices
         DMRGSCFintegrals * theRotatedTEI;
         
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
         
         //Space for the DMRG 1DM
         double * DMRG1DM;

         //Space for the DMRG 2DM
         double * DMRG2DM;
         
         //Helper function to fetch DOCC and SOCC from filename
         void allocateAndFillOCC(const string filename);
         
         //Helper function to copy the DOCC and SOCC arrays
         void allocateAndFillOCC(int * DOCCin, int * SOCCin);
         
         //Helper function to check HF
         void checkHF();
         
         //Fill Econst and Tmat of HamDMRG
         void fillConstAndTmatDMRG(Hamiltonian * HamDMRG) const;
         
         //Calculate the gradient, return function is the gradient 2-norm
         static double calcGradient(const DMRGSCFmatrix * localFmat, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * gradient);
         
         //Calculate the hessian
         static void calcHessian(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * hessian, const int rowjump);
         
         //Fmat function as defined by Eq. (11) in the Siegbahn paper.
         DMRGSCFmatrix * theFmatrix;
         
         //The Coulomb and exchange interaction with the occupied and active electrons respectively
         DMRGSCFmatrix * theQmatOCC;
         DMRGSCFmatrix * theQmatACT;
         DMRGSCFmatrix * theQmatWORK;
         DMRGSCFmatrix * theTmatrix;
         void rotateOldToNew(DMRGSCFmatrix * myMatrix);
         void buildTmatrix();
         void constructCoulombAndExchangeMatrixInOrigIndices(DMRGSCFmatrix * densityMatrix, DMRGSCFmatrix * resultMatrix);
         void buildQmatOCC();
         void buildQmatACT();
         
         //The Wmat_tilde function as defined by Eq. (20b) in the Siegbahn paper (see class header for specific definition)
         DMRGSCFwtilde * wmattilde;
         
         //The Wmat function as defined by Eq. (21b) in the Siegbahn paper
         static double Wmat(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const int irrep_pq, const int irrep_rs, const int relindexP, const int relindexQ, const int relindexR, const int relindexS);
         
         //Function to get the occupancies to obtain coefficients of certain Slater determinants for neutral C2. Important to figure out diatomic D(inf)h symmetries when calculating them in D2h symmetry. The function is not basis set and active space dependent (at least if no B2g, B3g, B2u and B3u orbitals are condensed).
         void PrintCoeff_C2(DMRG * theDMRG);
         
   };
}

#endif
