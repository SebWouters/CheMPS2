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

#ifndef EDMISTONRUEDENBERG_CHEMPS2_H
#define EDMISTONRUEDENBERG_CHEMPS2_H

#include "Options.h"
#include "DMRGSCFunitary.h"
#include "Hamiltonian.h"

namespace CheMPS2{
/** EdmistonRuedenberg class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date August 14, 2014
    
    \section introEdRued Introduction
    
    DMRG is not invariant with respect to orbital choice and ordering [LOC1]. An MPS introduces an artificial correlation length. Highly correlated orbitals should be placed close to each other for optimal DMRG performance. The correlation between localized orbitals in extended molecules decreases with increasing spatial separation. The exchange matrix directly reflects orbital overlap and separation. Minimizing its bandwidth yields a good ordering for extended molecules [LOC2, LOC3].
    
    \section locEdRued Localization
    
    Physics notation is assumed for the two-body matrix elements (the electron repulsion integrals with eightfold permutation symmetry):
    \f[
    v_{ij;kl} = \left( i(\vec{r}_1) j(\vec{r}_2) \left| \frac{1}{\mid \vec{r}_1 - \vec{r}_2 \mid} \right| k(\vec{r}_1) l(\vec{r}_2) \right).
    \f]
    Edmiston and Ruedenberg proposed to maximize the cost function
    \f[
    I = \sum\limits_i v_{ii;ii}
    \f]
    in order to obtain localized orthonormal orbitals [LOC4]. When the spatial extent of an orbital is small, its self-repulsion will indeed be high. CheMPS2 obtains localized orbitals (per irrep) with an augmented Hessian Newton-Raphson maximization of \f$I\f$, which can be called with the function CheMPS2::EdmistonRuedenberg::Optimize.
    
    The orbitals which maximize the cost function \f$ I = \sum\limits_i v_{ii;ii} \f$ can be parametrized with an orthogonal matrix \f$\mathbf{U} = \exp(\mathbf{X}(\vec{x}))\f$ (see CheMPS2::CASSCF):
    \f[
    I(\mathbf{U}) = \sum\limits_{ijklm} \left[ \mathbf{U} \right]_{ij} \left[ \mathbf{U} \right]_{ik} \left[ \mathbf{U} \right]_{il} \left[ \mathbf{U} \right]_{im} v_{jk;lm}.
    \f]
    The following derivatives from Ref. [LOC5] are useful:
    \f[
    \left[ \frac{\partial \left[ \mathbf{U} \right]_{\alpha \beta}}{\partial x_{pq}} \right]_0 = \delta_{p \alpha} \delta_{q \beta} - \delta_{q \alpha} \delta_{p \beta},
    \f]
    \f{eqnarray*}{
    2 \left[ \frac{\partial^2 \left[ \mathbf{U} \right]_{\alpha \beta}}{\partial x_{pq} \partial x_{rs}} \right]_0 & = & \delta_{qr} (\delta_{p \alpha} \delta_{s \beta} + \delta_{s \alpha} \delta_{p \beta} ) \\ & - & \delta_{pr} (\delta_{q \alpha} \delta_{s \beta} + \delta_{s \alpha} \delta_{q \beta} ) \\ & - & \delta_{qs} (\delta_{p \alpha} \delta_{r \beta} + \delta_{r \alpha} \delta_{p \beta} ) \\ & + & \delta_{ps} (\delta_{q \alpha} \delta_{r \beta} + \delta_{r \alpha} \delta_{q \beta} ),
    \f}
    in which \f$x_{pq}\f$ is the element in \f$\vec{x}\f$ which corresponds to \f$\left[ \mathbf{X} \right]_{pq}\f$. Keeping the eightfold permutation symmetry of \f$v_{ij;kl}\f$ in mind, the gradient of \f$I\f$ with respect to orbital rotations is
    \f[
    \left[ \frac{\partial I}{\partial x_{pq}} \right]_0 = 4 (v_{pp;pq} - v_{qq;qp}),
    \f]
    and the Hessian
    \f{eqnarray*}{
    \left[ \frac{\partial^2 I}{\partial x_{pq} \partial x_{rs}} \right]_0 & = & \delta_{pr} (8v_{pp;qs} + 4v_{pq;ps} - 2v_{qq;qs} - 2v_{ss;sq}) \\ & + & \delta_{qs} (8v_{qq;pr} + 4v_{qp;qr} - 2v_{pp;pr} - 2v_{rr;rp}) \\ & - & \delta_{ps} (8v_{pp;qr} + 4v_{pq;pr} - 2v_{qq;qr} - 2v_{rr;rq})\\ & - & \delta_{qr} (8v_{qq;ps} + 4v_{qp;qs} - 2v_{pp;ps} - 2v_{ss;sp}).
    \f}
    
    \section fiedlerEdRued Minimizing the bandwidth of the exchange matrix

    Once localized orbitals are obtained, they can be ordered so that the bandwidth of the exchange matrix \f$\mathbf{K}\f$ is minimal. This matrix has only nonnegative entries [LOC6]:
    \f[
    \left[ \mathbf{K} \right]_{ij} = v_{ij;ji} \geq 0.
    \f]
    A fast approximate bandwidth minimization is obtained by the so-called Fiedler vector [LOC7, LOC8, LOC9, LOC10]. Consider the minimization of the cost function
    \f[
    F(\vec{z}) ~ = ~ \frac{1}{2} \sum\limits_{k \neq l} \left[ \mathbf{K} \right]_{kl} (z_k - z_l)^2 ~ \geq 0,
    \f]
    where the (continuous) variable \f$z_k\f$ denotes the optimal position of orbital \f$k\f$. In order to fix translational invariance and normalization of the solution, the constraints \f$\sum_i z_i = 0\f$ and \f$\vec{z}^T \vec{z} = 1\f$ are imposed. The cost function can be rewritten as
    \f[
    F(\vec{z}) ~ = ~ \vec{z}^T \mathbf{L} \vec{z} ~ \geq 0,
    \f]
    with the symmetric and positive semidefinite Laplacian \f$\mathbf{L}\f$:
    \f[
    \left[ \mathbf{L} \right]_{k \neq l} = - \left[ \mathbf{K} \right]_{kl},
    \f]
    \f[
    \left[ \mathbf{L} \right]_{k k} = \sum\limits_{m \neq k} \left[ \mathbf{K} \right]_{km}.
    \f]
    The constant vector \f$\vec{1}\f$ is an eigenvector of \f$\mathbf{L}\f$ with eigenvalue 0, which is discarded due to the translational invariance constraint. When the exchange matrix \f$\mathbf{K}\f$ is not block-diagonal, all other eigenvalues are positive. The eigenvector with smallest nonzero eigenvalue is the Fiedler vector, the solution to the constrained minimization of \f$F(\vec{z})\f$ [LOC7, LOC8, LOC9, LOC10]. When the indices of the orbitals are permuted so that they are ordered according to the continuous variables in the Fiedler vector, a fast approximate bandwidth minimization of the exchange matrix \f$\mathbf{K}\f$ is performed. CheMPS2 orders the localized orbitals according to the Fiedler vector of the exchange matrix, which is performed in the function CheMPS2::EdmistonRuedenberg::FiedlerExchange.

    \section biblioEdRued References
    
    [LOC1]  S. Wouters and D. Van Neck, European Physical Journal D 68, 272 (2014). http://dx.doi.org/10.1140/epjd/e2014-50500-1 \n
    [LOC2]  W. Mizukami, Y. Kurashige and T. Yanai, Journal of Chemical Theory and Computation 9, 401-407 (2013). http://dx.doi.org/10.1021/ct3008974 \n
    [LOC3]  N. Nakatani and G.K.-L. Chan, Journal of Chemical Physics 138, 134113 (2013). http://dx.doi.org/10.1063/1.4798639 \n
    [LOC4]  C. Edmiston and K. Ruedenberg, Reviews of Modern Physics 35, 457-464 (1963). http://dx.doi.org/10.1103/RevModPhys.35.457 \n
    [LOC5]  P.E.M. Siegbahn, J. Almlof, A. Heiberg and B.O. Roos, Journal of Chemical Physics 74, 2384-2396 (1981). http://dx.doi.org/10.1063/1.441359 \n
    [LOC6]  A. Auerbach, Interacting Electrons and Quantum Magnetism, Springer, Heidelberg, 1994. \n
    [LOC7]  M. Fiedler, Czechoslovak Mathematical Journal 23, 298-305 (1973). http://dml.cz/dmlcz/101168 \n
    [LOC8]  M. Fiedler, Czechoslovak Mathematical Journal 25, 619-633 (1975). http://dml.cz/dmlcz/101357 \n
    [LOC9]  G. Barcza, O. Legeza, K. H. Marti, M. Reiher, Physical Review A 83, 012508 (2011). http://dx.doi.org/10.1103/PhysRevA.83.012508 \n
    [LOC10] M.W. Berry, Fiedler ordering, http://web.eecs.utk.edu/~mberry/order/node9.html (1996).
    
    \section boysLocal Boys localization
    
    In the Boys localization method, the cost function
    \f[
      J = \frac{1}{2} \sum\limits_{ij} \left( \braket{ \phi_i \mid \vec{r} \mid \phi_i } - \braket{ \phi_j \mid \vec{r} \mid \phi_j } \right)^2 = \frac{1}{2} \sum\limits_{ij} \left( \vec{r}_{ii} - \vec{r}_{jj} \right)^2
    \f]
    is maximized. With an orthogonal orbital rotation \f$\mathbf{U} = \exp(\mathbf{X}(\vec{x}))\f$, the cost function becomes
    \f[
      J(\mathbf{U}) = \frac{1}{2} \sum\limits_{i,j} \left( \left[ \mathbf{U} \right]_{i\alpha} \left[ \mathbf{U} \right]_{i\beta} \vec{r}_{\alpha\beta} - \left[ \mathbf{U} \right]_{j\alpha} \left[ \mathbf{U} \right]_{j\beta} \vec{r}_{\alpha\beta} \right)^2.
    \f]
    The gradient and hessian of J with respect to orbital rotations are given by:
    \f{eqnarray*}{
      \left[ \frac{\partial J}{\partial x_{pq}} \right]_0 & = & 2 N_{orb} \left[ \left( \vec{r}_{pq} + \vec{r}_{qp} \right).\left( \vec{r}_{pp} - \vec{r}_{qq} \right) \right] \\
      \left[ \frac{\partial^2 J}{\partial x_{pq} \partial x_{rs}} \right]_0 & = & 2 N_{orb} \left( \delta_{pr} + \delta_{qs} - \delta_{qr} - \delta_{ps} \right) \left[ \left( \vec{r}_{pq} + \vec{r}_{qp} \right).\left( \vec{r}_{rs} + \vec{r}_{sr} \right) \right] \\
      & + & N_{orb} \delta_{pr} \left[ \left( \vec{r}_{qs} + \vec{r}_{sq} \right).\left( 2 \vec{r}_{pp} - \vec{r}_{qq} - \vec{r}_{ss} \right) \right] \\
      & + & N_{orb} \delta_{qs} \left[ \left( \vec{r}_{pr} + \vec{r}_{rp} \right).\left( 2 \vec{r}_{qq} - \vec{r}_{pp} - \vec{r}_{rr} \right) \right] \\
      & - & N_{orb} \delta_{qr} \left[ \left( \vec{r}_{ps} + \vec{r}_{sp} \right).\left( 2 \vec{r}_{qq} - \vec{r}_{pp} - \vec{r}_{ss} \right) \right] \\
      & - & N_{orb} \delta_{ps} \left[ \left( \vec{r}_{qr} + \vec{r}_{rq} \right).\left( 2 \vec{r}_{pp} - \vec{r}_{qq} - \vec{r}_{rr} \right) \right].
    \f}
    (Written down here for future reference.)
    
*/
   class EdmistonRuedenberg{

      public:
      
         //! Constructor
         /** \param Vmat The active space two-electron integrals
             \param group The group number
             \param printLevelIn If 0: nothing is printed. If 1: info at beginning and end. If >1: intermediary info as well. */
         EdmistonRuedenberg( const FourIndex * Vmat, const int group, const int printLevelIn=1 );
         
         //! Destructor
         virtual ~EdmistonRuedenberg();
         
         //! Maximize the Edmiston-Ruedenberg cost function
         /** \param temp1 Work memory of at least max(dim(irrep(Ham)))^4
             \param temp2 Work memory of at least max(dim(irrep(Ham)))^4
             \param startFromRandomUnitary The name of the variable says it all. Useful if the point group symmetry is reduced to make use of locality in the DMRG calculations, but when the molecular orbitals still belong to the full point group.
             \param gradThreshold Stop if the norm of the gradient is smaller than this value
             \param maxIter Stop if maxIter iterations have been performed (when gradThreshold is not reached)
             \return The value of the optimized cost function */
         double Optimize(double * temp1, double * temp2, const bool startFromRandomUnitary, const double gradThreshold=EDMISTONRUED_gradThreshold, const int maxIter=EDMISTONRUED_maxIter);
         
         //! Permute the orbitals so that the bandwidth of the exchange matrix is (approximately) minimized (per irrep)
         /** \param maxlinsize max(dim(irrep(Ham)))
             \param temp1 Work memory of at least 4*max(dim(irrep(Ham)))^2
             \param temp2 Work memory of at least 4*max(dim(irrep(Ham)))^2 */
         void FiedlerExchange(const int maxlinsize, double * temp1, double * temp2);

         //! Permute the orbitals so that the bandwidth of the exchange matrix is (approximately) minimized
         /** \param dmrg2ham At exit, dmrg2ham contains the new orbital ordering */
         void FiedlerGlobal( int * dmrg2ham ) const;

         //! Get the pointer to the unitary to use in DMRGSCF
         /** \return Pointer to the unitary which defines the localized orbitals */
         DMRGSCFunitary * getUnitary();
         
      private:
      
         //Pointer to the active space two-electron integrals
         const FourIndex * VMAT_ORIG;
         
         //The print level
         int printLevel;
         
         //The symmetry information object
         Irreps SymmInfo;
         
         //The DMRGSCF index handler (in order to be able to recycle the DMRGSCFunitary object)
         DMRGSCFindices * iHandler;
         
         //DMRGSCF unitary
         DMRGSCFunitary * unitary;
         
         //The rotated two-body matrix elements
         FourIndex * VmatRotated;
         
         //Calculate the gradient, hessian and update
         double augmentedHessianNewtonRaphson(double * gradient, double * temp1, double * temp2) const;
         double calcGradientValue(const int irrep, const int p, const int q) const;
         double calcHessianValue( const int irrep, const int p, const int q, const int r, const int s) const;
         
         //Calculate the cost function
         double costFunction() const;
         
         //Fiedler vector helper function
         double FiedlerExchangeCost() const;
         static double FiedlerGlobalCost( const DMRGSCFindices * idx, const FourIndex * VMAT_LOCAL, int * dmrg2ham );
         void Fiedler(const int irrep, int * reorder, double * laplacian, double * temp2);
         
   };
}

#endif
