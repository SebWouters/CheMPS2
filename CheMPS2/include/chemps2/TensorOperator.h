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

#ifndef TENSOROPERATOR_CHEMPS2_H
#define TENSOROPERATOR_CHEMPS2_H

#include "Tensor.h"
#include "TensorT.h"
#include "SyBookkeeper.h"

namespace CheMPS2{
/** TensorOperator class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date October 15, 2015
    
    The TensorOperator class is a storage and update class for tensor operators with a given:\n
    - spin (two_j)
    - particle number (n_elec)
    - point group irrep (n_irrep).
    
    It replaces the previous classes TensorDiag, TensorSwap, TensorS0Abase, TensorS1Bbase, TensorF0Cbase, TensorF1Dbase, TensorA, TensorB, TensorC, and TensorD. Their storage and update functions have a common origin. The boolean prime_last denotes whether in which convention the tensor operator is stored:\n
    - \f$ \braket{ j m J M | j' m' } \braket{ j ( N, I ) || J ( n\_elec, n\_irrep ) || j' ( N + n\_elec, I \times n\_irrep ) } \f$ (prime_last == true)
    - \f$ \braket{ j' m' J M | j m } \braket{ j ( N, I ) || J ( n\_elec, n\_irrep ) || j' ( N + n\_elec, I \times n\_irrep ) } \f$ (prime_last == false).
    
    This determines the specific reduced update formulae when contracting with the Clebsch-Gordan coefficients of the reduced MPS tensors. */
   class TensorOperator : public Tensor{

      public:

         //! Constructor
         /** \param boundary_index The boundary index
             \param two_j Twice the spin of the tensor operator
             \param n_elec How many electrons there are more in the symmetry sector of the lower leg compared to the upper leg
             \param n_irrep The (real-valued abelian) point group irrep difference between the symmetry sectors of the lower and upper legs (see Irreps.h)
             \param moving_right If true: sweep from left to right. If false: sweep from right to left.
             \param prime_last Convention in which the tensor operator is stored (see class information)
             \param jw_phase Whether or not to include a Jordan-Wigner phase due to the fermion anti-commutation relations
             \param bk_up   Symmetry bookkeeper of the upper MPS
             \param bk_down Symmetry bookkeeper of the lower MPS */
         TensorOperator( const int boundary_index, const int two_j, const int n_elec, const int n_irrep, const bool moving_right, const bool prime_last, const bool jw_phase, const SyBookkeeper * bk_up, const SyBookkeeper * bk_down );

         //! Destructor
         virtual ~TensorOperator();

         //! Get the number of symmetry blocks
         /** \return The number of symmetry blocks */
         int gNKappa() const;

         //! Get the pointer to the storage
         /** return pointer to the storage */
         double * gStorage();

         //! Get the index corresponding to a certain tensor block
         /** \param N1 The up particle number sector
             \param TwoS1 The up spin symmetry sector
             \param I1 The up irrep sector
             \param N2 The down particle number sector
             \param TwoS2 The down spin symmetry sector
             \param I2 The down irrep sector
             \return The kappa corresponding to the input parameters; -1 means no such block */
         int gKappa( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 ) const;

         //! Get the storage jump corresponding to a certain tensor block
         /** \param kappa The symmetry block
             \return kappa2index[kappa], the memory jumper to a certain block */
         int gKappa2index( const int kappa ) const;

         //! Get the pointer to the storage of a certain tensor block
         /** \param N1 The up particle number sector
             \param TwoS1 The up spin symmetry sector
             \param I1 The up irrep sector
             \param N2 The down particle number sector
             \param TwoS2 The down spin symmetry sector
             \param I2 The down irrep sector
             \return Pointer to the storage of the specified tensor block; NULL means no such block */
         double * gStorage( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 );

         //! Get the boundary index
         /** \return the index */
         int gIndex() const;

         //! Get twice the spin of the tensor operator
         /** \return Twice the spin of the tensor operatorg */
         int get_2j() const;

         //! Get how many electrons there are more in the symmetry sector of the lower leg compared to the upper leg
         /** \return How many electrons there are more in the symmetry sector of the lower leg compared to the upper leg */
         int get_nelec() const;

         //! Get the (real-valued abelian) point group irrep difference between the symmetry sectors of the lower and upper legs (see Irreps.h)
         /** \return The (real-valued abelian) point group irrep difference between the symmetry sectors of the lower and upper legs (see Irreps.h) */
         int get_irrep() const;

         //! Clear and update
         /** \param previous The previous TensorOperator needed for the update
             \param mps_tensor_up   The upper MPS tensor needed for the update 
             \param mps_tensor_down The lower MPS tensor needed for the update
             \param workmem Work memory */
         void update( TensorOperator * previous, TensorT * mps_tensor_up, TensorT * mps_tensor_down, double * workmem );

         //! daxpy for TensorOperator
         /** \param alpha The prefactor
             \param to_add The TensorOperator x which should be added: this <-- this + alpha * to_add */
         void daxpy( double alpha, TensorOperator * to_add );

         //! daxpy_transpose for C- and D-tensors (with special spin-dependent factors)
         /** \param alpha The prefactor
             \param to_add The TensorOperator x which should be added: this <-- this + alpha * special_spin_dependent_factor * to_add^T */
         void daxpy_transpose_tensorCD( const double alpha, TensorOperator * to_add );

         //! Set all storage variables to 0.0
         void clear();

         //! Make the in-product of two TensorOperator
         /** \param buddy The second tensor
             \param trans If trans == 'N' a regular ddot is taken. If trans == 'T' and n_elec==0, the in-product with buddy's transpose is made.
             \return The in-product */
         double inproduct( TensorOperator * buddy, const char trans ) const;

      protected:

         //! The bookkeeper of the upper MPS
         const SyBookkeeper * bk_up;

         //! The bookkeeper of the lower MPS
         const SyBookkeeper * bk_down;

         //! Twice the spin of the tensor operator
         int two_j;

         //! How many electrons there are more in the symmetry sector of the lower leg compared to the upper leg
         int n_elec;

         //! The (real-valued abelian) point group irrep difference between the symmetry sectors of the lower and upper legs (see Irreps.h)
         int n_irrep;

         //! Whether or not moving right
         bool moving_right;

         //! The up particle number sector
         int * sector_nelec_up;

         //! The up spin symmetry sector
         int * sector_irrep_up;

         //! The up spin symmetry sector
         int * sector_spin_up;

         //! The down spin symmetry sector (pointer points to sectorTwoS1 if two_j == 0)
         int * sector_spin_down;

         //! Update moving right
         /** \param ikappa The tensor block which should be updated
             \param previous The previous TensorOperator needed for the update
             \param mps_tensor_up   The upper MPS tensor needed for the update 
             \param mps_tensor_down The lower MPS tensor needed for the update
             \param workmem Work memory */
         void update_moving_right( const int ikappa, TensorOperator * previous, TensorT * mps_tensor_up, TensorT * mps_tensor_down, double * workmem );

         //! Update moving left
         /** \param ikappa The tensor block which should be updated
             \param previous The previous TensorOperator needed for the update
             \param mps_tensor_up   The upper MPS tensor needed for the update 
             \param mps_tensor_down The lower MPS tensor needed for the update
             \param workmem Work memory */
         void update_moving_left( const int ikappa, TensorOperator * previous, TensorT * mps_tensor_up, TensorT * mps_tensor_down, double * workmem );

         //! Convention in which the tensor operator is stored (see class information)
         bool prime_last;

         //! Whether or not to include a Jordan-Wigner phase due to the fermion anti-commutation relations
         bool jw_phase;

      private:

   };
}

#endif
