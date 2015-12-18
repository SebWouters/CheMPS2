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
    
    \section biblio References

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
         /** \param ham_in Hamiltonian containing the matrix elements of the Hamiltonian for which a CASSCF calculation is desired
             \param docc_in Array containing the number of doubly occupied HF orbitals per irrep
             \param socc_in Array containing the number of singly occupied HF orbitals per irrep
             \param nocc_in Array containing the number of doubly occupied (inactive) orbitals per irrep
             \param ndmrg_in Array containing the number of active orbitals per irrep
             \param nvirt_in Array containing the number of virtual (secondary) orbitals per irrep */
         CASPT2(DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei, DMRGSCFmatrix * fock, double * one_dm, double * two_dm, double * three_dm, double * contract);
         
         //! Destructor
         virtual ~CASPT2();
         
         //! Get the number of irreps
         /** \return The number of irreps */
         int get_num_irreps();
         
         
         
      private:
      
         // The number of occupied, active and virtual orbitals per irrep (externally allocated and deleted)
         const DMRGSCFindices * indices;
         
         // The one-electron integrals (externally allocated and deleted)
         const DMRGSCFmatrix * oei;
         
         // The Fock matrix (externally allocated and deleted)
         const DMRGSCFmatrix * fock;
         
         // The active space 1-RDM (externally allocated and deleted)
         double * one_rdm;
         
         // The active space 2-RDM (externally allocated and deleted)
         double * two_rdm;
         
         // The active space 3-RDM (externally allocated and deleted)
         double * three_rdm;
         
         // The active space 4-RDM contracted with the Fock operator (externally allocated and deleted)
         double * contracted;
         
         // The number of irreps
         int num_irreps;
         
         // The Fock operator expectation value
         double E_FOCK;
         
         // Calculate the expectation value of the Fock operator
         double calc_fock_expectation() const;
         
         // Calculate the total vector length and the partitioning of the vector in blocks
         int vector_helper();
         long long total_vector_length() const; // For debugging purposes
         
         // Once make_S**() has been calles, these overlap matrices can be used to contruct the RHS of the linear problem
         void construct_rhs( const DMRGSCFintegrals * integrals );
         
         // Variables for the partitioning of the vector in blocks
         int * jump;
         int * size_AC;
         int * size_D;
         int * size_BF_singlet;
         int * size_BF_triplet;
         
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
         
   };
}

#endif
