/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

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

#ifndef CASPT2_CHEMPS2_H
#define CASPT2_CHEMPS2_H

#include "DMRGSCFindices.h"
#include "DMRGSCFintegrals.h"
#include "DMRGSCFmatrix.h"

#define CHEMPS2_CASPT2_A         0
#define CHEMPS2_CASPT2_B_SINGLET 1
#define CHEMPS2_CASPT2_B_TRIPLET 2
#define CHEMPS2_CASPT2_C         3
#define CHEMPS2_CASPT2_D         4
#define CHEMPS2_CASPT2_E_SINGLET 5
#define CHEMPS2_CASPT2_E_TRIPLET 6
#define CHEMPS2_CASPT2_F_SINGLET 7
#define CHEMPS2_CASPT2_F_TRIPLET 8
#define CHEMPS2_CASPT2_G_SINGLET 9
#define CHEMPS2_CASPT2_G_TRIPLET 10
#define CHEMPS2_CASPT2_H_SINGLET 11
#define CHEMPS2_CASPT2_H_TRIPLET 12
#define CHEMPS2_CASPT2_NUM_CASES 13

namespace CheMPS2{
/** CASPT2 class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date December 11, 2015
    
    \section theo_caspt2 Information
    
    The CASPT2 class contains the functions to perform second order multireference perturbation theory on top of a CASSCF wavefuntion [CASPT1, CASPT2]. CASPT2 has recently also been used with DMRG as active space solver: the active space 4-RDM contracted with the Fock operator, together with the 1-, 2- and 3-RDM are required thereto [CASPT3]. Alternatively, cumulant approximations can be used as well [CASPT4]. To mitigate problems with CASPT2 several modifications of the zeroth order Hamiltonian have been introduced: IPEA corrections [CASPT5], the g1 term [CASPT6], real level shifts [CASPT7] and imaginary level shifts [CASPT8].
    
    \section biblio_caspt2 References

    [CASPT1]  K. Andersson, P.-A. Malmqvist, B.O. Roos, A.J. Sadlej and K. Wolinski, Journal of Physical Chemistry 94, 5483-5488 (1990). http://dx.doi.org/10.1021/j100377a012 \n
    [CASPT2]  K. Andersson, P.‐A. Malmqvist and B.O. Roos, Journal of Chemical Physics 96, 1218-1226 (1992). http://dx.doi.org/10.1063/1.462209 \n
    [CASPT3]  Y. Kurashige and T. Yanai, Journal of Chemical Physics 135, 094104 (2011). http://dx.doi.org/10.1063/1.3629454 \n
    [CASPT4]  Y. Kurashige, J. Chalupsky, T.N. Lan and T. Yanai, Journal of Chemical Physics 141, 174111 (2014). http://dx.doi.org/10.1063/1.4900878 \n
    [CASPT5]  G. Ghigo, B.O. Roos and P.-A. Malmqvist, Chemical Physics Letters 396, 142-149 (2004). http://dx.doi.org/10.1016/j.cplett.2004.08.032 \n
    [CASPT6]  K. Andersson, Theoretica Chimica Acta 91, 31-46 (1995). http://dx.doi.org/10.1007/BF01113860 \n
    [CASPT7]  B.O. Roos and K. Andersson, Chemical Physics Letters 245, 215-223 (1995). http://dx.doi.org/10.1016/0009-2614(95)01010-7 \n
    [CASPT8]  N. Forsberg and P.-A. Malmqvist, Chemical Physics Letters 274, 196-204 (1997). http://dx.doi.org/10.1016/S0009-2614(97)00669-6 \n
*/
   class CASPT2{

      public:

         //! Constructor
         /** \param idx      DMRGSCFindices which contain the partitioning into occupied, active, and virtual orbitals per irrep
             \param ints     The two-electron integrals needed for CASSCF and CASPT2
             \param oei      The one-electron integrals
             \param fock     The fock matrix of CASPT2 with the occupied-occupied block and the virtual-virtual block diagonal!
             \param one_dm   The spin-summed one-particle density matrix one_dm[i+L*j] = sum_sigma < a^+_i,sigma a_j,sigma > (with L the number DMRG orbitals)
             \param two_dm   The spin-summed two-particle density matrix two_dm[i+L*(j+L*(k+L*l))] = sum_sigma,tau < a^+_i,sigma a^+_j,tau a_l,tau a_k,sigma > (with L the number DMRG orbitals)
             \param three_dm The spin-summed three-particle density matrix three_dm[i+L*(j+L*(k+L*(l+L*(m+L*n))))] = sum_z,tau,s < a^+_{i,z} a^+_{j,tau} a^+_{k,s} a_{n,s} a_{m,tau} a_{l,z} > (with L the number DMRG orbitals)
             \param contract The spin-summed four-particle density matrix contracted with the fock operator (as performed by Cumulant::gamma4_fock_contract_ham or FCI::Fock4RDM) contract[i+L*(j+L*(k+L*(p+L*(q+L*r))))] = sum_{l,t,sigma,tau,s,z} fock(l,t) < a^+_{i,sigma} a^+_{j,tau} a^+_{k,s} a^+_{l,z} a_{t,z} a_{r,s} a_{q,tau} a_{p,sigma} > (with L the number DMRG orbitals) */
         CASPT2(DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei, DMRGSCFmatrix * fock, double * one_dm, double * two_dm, double * three_dm, double * contract);

         //! Destructor
         virtual ~CASPT2();
         
         //! Solve the CASPT2 energy
         double solve( const bool diag_only = false ) const;

      private:
      
         // The number of occupied, active and virtual orbitals per irrep (externally allocated and deleted)
         const DMRGSCFindices * indices;
         
         // The Fock matrix (externally allocated and deleted)
         const DMRGSCFmatrix * fock;
         
         // The active space 1-RDM (externally allocated and deleted)
         double * one_rdm;
         
         // The active space 2-RDM (externally allocated and deleted)
         double * two_rdm;
         
         // The active space 3-RDM (externally allocated and deleted)
         double * three_rdm;
         
         // The active space 4-RDM contracted with the Fock operator (externally allocated and deleted)
         double * f_dot_4dm;
         
         // The active space 3-RDM contracted with the Fock operator (allocated and deleted in this class)
         double * f_dot_3dm;
         
         // The active space 2-RDM contracted with the Fock operator (allocated and deleted in this class)
         double * f_dot_2dm;
         
         // The active space 1-RDM contracted with the Fock operator
         double f_dot_1dm;
         
         // The sum of the Fock operator corresponding to occupied orbitals
         double sum_f_kk;
         
         // The number of irreps
         int num_irreps;

         // The Fock operator expectation value
         double E_FOCK;

         // Calculate the expectation value of the Fock operator
         void create_f_dots();

         // Calculate the total vector length and the partitioning of the vector in blocks
         int vector_helper();
         long long debug_total_length() const;
         
         // Once make_S**() has been calles, these overlap matrices can be used to contruct the RHS of the linear problem
         void construct_rhs( const DMRGSCFmatrix * oei, const DMRGSCFintegrals * integrals );
         
         // Fill result with the diagonal elements of the Fock operator
         void diagonal( double * result, const double ovlp_prefactor ) const;
         
         // Fill result with Fock operator times vector
         void matvec( double * vector, double * result, double * diag_fock ) const;
         
         // Variables for the partitioning of the vector in blocks
         int * jump;
         int * size_A;
         int * size_C;
         int * size_D;
         int * size_E;
         int * size_G;
         int * size_B_singlet;
         int * size_B_triplet;
         int * size_F_singlet;
         int * size_F_triplet;

         // Functions for the partitioning of the vector in blocks
         int get_maxsize() const;
         static int jump_AC_active( const DMRGSCFindices * idx, const int irrep_t, const int irrep_u, const int irrep_v );
         static int jump_BF_active( const DMRGSCFindices * idx, const int irrep_t, const int irrep_u, const int ST );
         static int jump_D_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_a );
         static int jump_B_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_j, const int ST );
         
         // The RHS of the linear problem
         double * vector_rhs;
         
         // Variables for the overlap
         double ** SAA;
         double ** SCC;
         double ** SDD;
         double ** SEE;
         double ** SGG;
         double ** SBB_singlet;
         double ** SBB_triplet;
         double ** SFF_singlet;
         double ** SFF_triplet;
         
         // Fill some helper variables for the overlap
         void make_SAA_SCC();
         void make_SDD();
         void make_SBB_SFF_singlet();
         void make_SBB_SFF_triplet();
         void make_SEE_SGG();
         
         // Variables for the diagonal part of the Fock operator
         double ** FAA;
         double ** FCC;
         double ** FDD;
         double ** FEE;
         double ** FGG;
         double ** FBB_singlet;
         double ** FBB_triplet;
         double ** FFF_singlet;
         double ** FFF_triplet;
         
         // Variables for the off-diagonal part of the Fock operator: Operator[ IL ][ IR ][ w ][ left + SIZE * right ]
         double **** FAD;
         double **** FCD;
         double *** FEH;
         double *** FGH;
         double **** FAB_singlet;
         double **** FAB_triplet;
         double **** FCF_singlet;
         double **** FCF_triplet;
         double **** FBE_singlet;
         double **** FBE_triplet;
         double **** FFG_singlet;
         double **** FFG_triplet;
         double **** FDE_singlet;
         double **** FDE_triplet;
         double **** FDG_singlet;
         double **** FDG_triplet;

         // Fill helper variables for the Fock operator
         void make_FAA_FCC();
         void make_FDD();
         void make_FBB_FFF_singlet();
         void make_FBB_FFF_triplet();
         void make_FEE_FGG();
         
         void make_FAD_FCD();
         void make_FEH_FGH();
         void make_FAB_FCF_singlet();
         void make_FAB_FCF_triplet();
         void make_FBE_FFG_singlet();
         void make_FBE_FFG_triplet();
         void make_FDE_FDG();

         // Diagonalize the overlap matrices and adjust jump, vector_rhs, and FXX accordingly
         void recreate();
         static int  recreatehelper1( double * FOCK, double * OVLP, int SIZE, double * work, double * eigs, int lwork );
         static void recreatehelper2( double * LEFT, double * RIGHT, double ** matrix, double * work, int OLD_LEFT, int NEW_LEFT, int OLD_RIGHT, int NEW_RIGHT, const int number );
         static void recreatehelper3( double * OVLP, int OLDSIZE, int NEWSIZE, double * rhs_old, double * rhs_new, const int num_rhs );
         
         
   };
}

#endif
