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

#ifndef DMRGSCFUNITARY_CHEMPS2_H
#define DMRGSCFUNITARY_CHEMPS2_H

#include "Options.h"
#include "DMRGSCFindices.h"
#include "DMRGSCFmatrix.h"
#include "DIIS.h"

namespace CheMPS2{
/** DMRGSCFunitary class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date July 11, 2014
    
    The DMRGSCFunitary class is a storage and manipulation class for the DMRGSCF orthogonal orbital rotation matrix. This matrix is blockdiagonal in the irreducible representations, and is formed by stepwise multiplying in new unitary rotations due to the augmented Hessian Newton-Raphson algorithm, see CheMPS2::CASSCF.
    
    \section buildexp Exponential of a skew-symmetric matrix
    
    The exponential of a real-valued skew-symmetric matrix \f$\mathbf{X} = -\mathbf{X}^T\f$ is an orthogonal matrix \f$\mathbf{U}\f$:
    \f[
    \mathbf{U}^T \mathbf{U} = \exp(\mathbf{X}^T) \exp(\mathbf{X}) =  \exp(- \mathbf{X}) \exp(\mathbf{X}) = \mathbf{I}.
    \f]
    A real-valued skew-symmetric matrix \f$\mathbf{X}\f$ has purely imaginary eigenvalues, which come in complex conjugate pairs \f$(i \lambda, -i\lambda)\f$. For matrices of odd dimension, there should hence always be one eigenvalue \f$0\f$. \f$\mathbf{B} = \mathbf{X} \mathbf{X}\f$ is then a symmetric matrix with nonpositive real-valued eigenvalues \f$-\lambda^2\f$. Nonzero eigenvalues of \f$\mathbf{B}\f$ occur twice. The eigenvectors \f$\mathbf{V}\f$ of \f$\mathbf{B} = \mathbf{V} diag(-\lambda^2) \mathbf{V}^T\f$ allow to make \f$\mathbf{X}\f$ block-diagonal. \f$\mathbf{C} = \mathbf{V}^T \mathbf{X} \mathbf{V}\f$ is skew-symmetric. Moreover, \f$\mathbf{C}\mathbf{C} = \mathbf{V}^T \mathbf{B} \mathbf{V}\f$ is diagonal: \f$\mathbf{C}\mathbf{C} = diag(-\lambda^2)\f$. \f$\mathbf{C}\f$ is hence block-diagonal with \f$1 \times 1\f$ blocks \f$\left[0\right]\f$ and \f$2 \times 2\f$ blocks
    \f[
    \left[ \begin{array}{cc} 0 & \lambda \\ -\lambda & 0 \end{array} \right].
    \f]
    The exponential of the \f$1 \times 1\f$ block \f$\left[0\right]\f$ is \f$\left[1\right]\f$, and the exponential of the \f$2 \times 2\f$ block is
    \f[
    \exp \left[ \begin{array}{cc} 0 & \lambda \\ -\lambda & 0 \end{array} \right] = \left[ \begin{array}{cc} \cos(\lambda) & \sin(\lambda) \\ -\sin(\lambda) & \cos(\lambda) \end{array} \right].
    \f]
    The matrix \f$\exp(\mathbf{C})\f$ can hence be easily calculated blockwise. The exponential of \f$\mathbf{X}\f$ is then obtained as \f$\exp(\mathbf{X}) = \mathbf{V} \exp(\mathbf{C}) \mathbf{V}^T\f$. It is calculated by the function CheMPS2::DMRGSCFunitary::updateUnitary.
    
    \section buildlog Logarithm of a special orthogonal matrix

    The reverse problem of finding a (nonunique) real-valued logarithm of a special orthogonal matrix \f$\mathbf{U}\f$ can be performed similarly. Since \f$\mathbf{U}\f$ is orthogonal (and hence norm-preserving), its eigenvalues all have norm 1:
    \f[
    \mathbf{U} = \mathbf{V}_{U} diag(e^{i \theta}) \mathbf{V}_{U}^{\dagger} = \mathbf{V}_{U}^* diag(e^{-i \theta}) \mathbf{V}_{U}^{T},
    \f]
    \f[
    \mathbf{U}^T = \mathbf{V}_{U}^* diag(e^{i \theta}) \mathbf{V}_{U}^{T} = \mathbf{V}_{U} diag(e^{-i \theta}) \mathbf{V}_{U}^{\dagger}.
    \f]
    For the second equalities, complex conjugation of the real-valued matrices \f$\mathbf{U}\f$ and \f$\mathbf{U}^T\f$ is used. The eigenvalues of \f$\mathbf{U}\f$ hence come in complex conjugate pairs \f$(e^{i \theta}, e^{-i \theta})\f$. Consider the symmetric matrix \f$ \mathbf{S} = \mathbf{U} + \mathbf{U}^T = \mathbf{V}_{U} diag(2\cos(\theta) ) \mathbf{V}_{U}^{\dagger} \f$. If \f$\cos(\theta) \neq \pm 1\f$, the eigenvalue \f$2\cos(\theta)\f$ occurs twice. For matrices of odd dimension \f$\cos(\theta)=+1\f$ always occurs an odd number of times. Construct the symmetric matrix \f$\mathbf{S}\f$ and diagonalize it (real-valued) as \f$\mathbf{S} = \mathbf{V}_{S} diag(2\cos(\theta)) \mathbf{V}_{S}^T\f$. The matrix \f$\mathbf{D} = \mathbf{V}_{S}^T \mathbf{U} \mathbf{V}_{S}\f$ is also a special orthogonal matrix with \f$1 \times 1\f$ blocks \f$\left[ \pm 1 \right]\f$ and \f$2 \times 2\f$ blocks
    \f[
    \left[ \begin{array}{cc} \cos(\theta) & \sin(\theta) \\ -\sin(\theta) & \cos(\theta) \end{array} \right].
    \f]
    Because we consider special orthogonal matrices, the \f$1 \times 1\f$ blocks \f$\left[ -1 \right]\f$ always occur an even number of times, and they can hence be considered as a special case of the \f$2 \times 2\f$ blocks. If we choose the branchcut for the logarithm on the negative real axis, the logarithm of the \f$1 \times 1\f$ block \f$\left[ 1 \right]\f$ is \f$\left[ 0 \right]\f$ and the logarithm of the \f$2\times 2\f$ block is 
    \f[
    \log \left[ \begin{array}{cc} \cos(\theta) & \sin(\theta) \\ -\sin(\theta) & \cos(\theta) \end{array} \right] = \left[ \begin{array}{cc} 0 & \theta \\ -\theta & 0 \end{array} \right],
    \f]
    with \f$\theta \in \left[ -\pi, \pi \right]\f$. The matrix \f$\log(\mathbf{D})\f$ can hence be easily calculated blockwise. The logarithm of \f$\mathbf{U}\f$ is then obtained as \f$\log(\mathbf{U}) = \mathbf{V}_{S} \log(\mathbf{D}) \mathbf{V}_{S}^T\f$. It is calculated by the function CheMPS2::DMRGSCFunitary::getLog.
*/
   class DMRGSCFunitary : public DMRGSCFmatrix{

      public:
      
         //! Constructor
         /** \param iHandler The DMRGSCF indices */
         DMRGSCFunitary( const DMRGSCFindices * iHandler );
         
         //! Destructor
         virtual ~DMRGSCFunitary();

         //! Get the number of variables in the x-parametrization of the unitary update
         /** \return The number of unique variables in the x-matrix */
         int getNumVariablesX() const;

         //! Update the unitary transformation
         /** \param workmem1 Work memory of at least 4*max(dim(irrep(Ham)))^2
             \param workmem2 Work memory of at least 4*max(dim(irrep(Ham)))^2
             \param vector The elements in X
             \param multiply Boolean whether exp(X)*U or exp(X) should become the new U. If multiply==true, U <-- exp(X)*U. If multiply==false, U <-- exp(X).
             \param compact Boolean which indicates how the elements X are stored */
         void updateUnitary( double * workmem1, double * workmem2, double * vector, const bool multiply, const bool compact );
         
         //! Rotate the unitary matrix
         /** \param eigenvecs The rotation vectors, in a memory block of size nOrbDMRG^2
             \param work Work memory, with size 2*max(dim(irrep(Ham)))^2 */
         void rotateActiveSpaceVectors( double * eigenvecs, double * work );
         
         //! Calculate the two-norm of U^T*U - I
         /** \param work Work memory */
         void CheckDeviationFromUnitary( double * work ) const;
         
         //! Obtain the logarithm of the unitary matrix
         /** \param vector Where the logarithm should be stored
             \param temp1 Work memory of at least 4*max(dim(irrep(Ham)))^2
             \param temp2 Work memory of at least 4*max(dim(irrep(Ham)))^2 */
         void getLog( double * vector, double * temp1, double * temp2 ) const;

         //! Orbitals are defined up to a phase factor. Make sure that the logarithm of each block of the unitary has determinant 1.
         /** \param temp1 Work memory of at least 4*max(dim(irrep(Ham)))^2
             \param temp2 Work memory of at least 4*max(dim(irrep(Ham)))^2 */
         void makeSureAllBlocksDetOne( double * temp1, double * temp2 );

         //! Save the unitary to disk
         /** \param filename Filename to store the unitary to */
         void saveU( const string filename=DMRGSCF_unitary_storage_name ) const;

         //! Load the unitary from disk
         /** \param filename Filename to load the unitary from */
         void loadU( const string filename=DMRGSCF_unitary_storage_name );

      private:

         //Number of variables in the x-matrix
         int x_linearlength;

         //Helper arrays to jump from linear x-matrix index to orbital indices and back
         int ** jumper;

         // Build in result the skew symmetric matrix X for irrep block irrep based on the elements in Xelem. If compact==true, they are stored in gradient form.
         void buildSkewSymmX( const int irrep, double * result, double * Xelem, const bool compact ) const;

         // Get the determinant, and in the process fill work1 with eigvec( U[ irrep ] + U^T[ irrep ] ) and work2 with work1^T U work1 (TRIDIAGONAL matrix).
         double get_determinant( const int irrep, double * work1, double * work2, double * work_eig, int lwork_eig ) const;

   };
}

#endif
