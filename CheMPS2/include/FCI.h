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

#ifndef FCI_CHEMPS2_H
#define FCI_CHEMPS2_H

#include "Problem.h"

namespace CheMPS2{
/** FCI class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 6, 2014
    
    The FCI class performs the full configuration interaction ground state calculation in a given particle number, irrep, and spin projection sector of a given Hamiltonian.
*/
   class FCI{

      public:
      
         //! Constructor
         /** \param Prob The problem contains the Hamiltonian, as well as the symmetry sector which should be targeted
             \param FCIverbose The FCI verbose level: 0 print nothing, 1 print start and solution, 2 print everything */
         FCI(CheMPS2::Problem * Prob, const int FCIverbose=2);
         
         //! Constructor
         /** \param Prob The problem contains the Hamiltonian, as well as the irrep which should be targeted. Nel_up and Nel_down are specified specifically now.
             \param Nel_up The number of up or alpha electrons of the new FCI object
             \param Nel_down The number of down or beta electrons of the new FCI object
             \param FCIverbose The FCI verbose level: 0 print nothing, 1 print start and solution, 2 print everything */
         FCI(CheMPS2::Problem * Prob, const unsigned int Nel_up, const unsigned int Nel_down, const int FCIverbose=2);
         
         //! Destructor
         virtual ~FCI();
         
         //! Getter for the number of orbitals
         /** \return The number of orbitals */
         unsigned int getL() const;

         //! Getter for the number of up or alpha electrons
         /** \return The number of up or alpha electrons */
         unsigned int getNel_up() const;

         //! Getter for the number of down or beta electrons
         /** \return The number of down or beta electrons */
         unsigned int getNel_down() const;
         
         //! Getter for the number of variables in the FCI vector
         /** \return The number of variables in the FCI vector */
         unsigned long long getVecLength() const;
         
         //! Get the number of irreps
         /** \return The number of irreps */
         unsigned int getNumIrreps() const;
         
         //! Get the target irrep
         /** \return The target irrep */
         int getTargetIrrep() const;
         
         //! Get the direct product of two irreps
         /** \param Irrep_row The first irrep
             \param Irrep_col The second irrep
             \return The direct product Irrep_row x Irrep_col */
         int getIrrepProduct(const int Irrep_row, const int Irrep_col) const;
         
         //! Get an orbital irrep
         /** \param orb The orbital index
             \return The irrep of orbital orb */
         int getOrb2Irrep(const int orb) const;
         
         //! Calculates the FCI ground state with Davidson's algorithm
         /** \param inoutput If inoutput!=NULL, vector with getVecLength() variables which contains on exit the solution of the FCI calculation
             \return The ground state energy */
         double GSDavidson(double * inoutput=NULL) const;
         
         //! Calculates the CI ground state in a small space of low-energy Slater determinants with lapack's dsyev
         /** \param Nslaters The number of low-energy Slater determinants to consider
             \param vec Vector with getVecLength() variables which contains on exit the solution of the small CI calculation
             \return The small CI ground state energy */
         double GSSmallCISpace(const unsigned int Nslaters, double * vec) const;
      
         //! Function which performs the Hamiltonian times Vector product
         /** \param input The vector on which the Hamiltonian should act
             \param output The Hamiltonian times input on exit */
         void HamTimesVec(double * input, double * output) const;
         
         //! Function which returns the diagonal elements of the FCI Hamiltonian
         /** \param diag Vector with getVecLength() variables which contains on exit the diagonal elements of the FCI Hamiltonian */
         void HamDiag(double * diag) const;
         
         //! Function which returns a FCI coefficient
         /** \param bits_up The bit string representation of the up or alpha electron Slater determinant
             \param bits_down The bit string representation of the down or beta electron Slater determinant
             \param vector The FCI vector with getVecLength() variables from which a coefficient is desired
             \return The corresponding FCI coefficient; 0.0 if something is not OK */
         double getFCIcoeff(int * bits_up, int * bits_down, double * vector) const;
         
         //! Set this FCI vector to an operator acting on a previous FCI vector
         /** \param whichOperator With which operator should be acted on the other FCI state: 0 means creator, 1 means annihilator, 2 means particle number
             \param isUp Boolean which denotes if the operator corresponds to an up (alpha) or down (beta) electron
             \param thisVector Vector with length getVecLength() where the result of the operation should be stored
             \param otherFCI FCI instance which corresponds to the FCI vector on which is acted
             \param otherVector Vector with length otherFCI->getVecLength() which contains the FCI vector on which is acted */
         void ActWithOperator(const int whichOperator, const bool isUp, const unsigned int orbIndex, double * thisVector, FCI * otherFCI, double * otherVector);
         
         //! Function which returns the 2-body matrix elements (in physics notation: orb1(r1) * orb2(r2) * orb3(r1) * orb4(r2) / |r1-r2|); the 1-body matrix elements are absorbed in the 2-body matrix elements as explained in the Problem class
         /** \param orb1 First orbital index (electron at position r1)
             \param orb2 Second orbital index (electron at position r2)
             \param orb3 Third orbital index (electron at position r1)
             \param orb4 Fourth orbital index (electron at position r2)
             \return The desired 2-body matrix element, augmented with the 1-body matrix elements */
         double getHmat(const int orb1, const int orb2, const int orb3, const int orb4) const;
         
         //! Function which returns the nuclear repulsion energy
         /** \return The nuclear repulsion energy */
         double getEconst() const;
         
         //! Construct the (spin-summed) 2-RDM of a FCI vector
         /** \param vector The FCI vector
             \param TwoRDM To store the 2-RDM; needs to be of size L^4 (L = number of orbitals); irrep symmetry shows in 2-RDM elements being zero; physics notation (see gHmat)
             \return The energy of the given FCI vector, calculated by contraction of the 2-RDM with Hmat */
         double Fill2RDM(double * vector, double * TwoRDM) const;
         
         //! Fill a vector with random numbers in the interval [-1,1[; used when output for GSDavidson is desired but no specific input can be given (as for the Hubbard model in the site basis)
         /** \param vecLength The length of the vector; normally getVecLength()
             \param vec The vector to fill with random numbers */
         static void FillRandom(const unsigned long long vecLength, double * vec);
         
         //! Copy a very long vector
         /** \param vecLength The vector length
             \param origin Vector to be copied
             \param target Where to copy the vector to */
         static void FCIdcopy(const unsigned long long vecLength, double * origin, double * target);
         
         //! Take the inproduct of two very long vectors
         /** \param vecLength The vector length
             \param vec1 The first vector
             \param vec2 The second vector
             \return The inproduct < vec1 | vec2 > */
         static double FCIddot(const unsigned long long vecLength, double * vec1, double * vec2);
         
         //! Calculate the 2-norm of a very long vector
         /** \param vecLength The vector length
             \param vec The vector
             \return The 2-norm of vec */
         static double FCIfrobeniusnorm(const unsigned long long vecLength, double * vec);
         
         //! Do lapack's daxpy for a very long vector
         /** \param vecLength The vector length
             \param alpha The scalar factor
             \param vec_x The vector which has to be added to vec_y in rescaled form
             \param vec_y The target vector */
         static void FCIdaxpy(const unsigned long long vecLength, const double alpha, double * vec_x, double * vec_y);
         
         //! Do lapack's dscal for a very long vector
         /** \param vecLength The vector length
             \param alpha The scalar factor
             \param vec The vector which has to be rescaled */
         static void FCIdscal(const unsigned long long vecLength, const double alpha, double * vec);
         
         //! Clear a very long vector
         /** \param vecLength The vector length
             \param vec The vector which has to be set to zero */
         static void FCIdclear(const unsigned long long vecLength, double * vec);
         
      private:
      
         //The FCI verbose level
         int FCIverbose;
         
         //Hamiltonian matrix elements
         double Econstant;
         double * Hmat;
         
         //Irrep (point group) information
         unsigned int NumIrreps;
         int TargetIrrep;
         int * IrrepProductTable;
         int * orb2irrep;
         
         //The number of orbitals
         unsigned int L;
         
         //The number of up and down electrons
         unsigned int Nel_up;
         unsigned int Nel_down;
         
         //The number of Slater determinants per irrep, which can be constructed with "Nel" electrons in the "L" orbitals with irreps "orb2irrep"
         unsigned int * numPerIrrep_up;
         unsigned int * numPerIrrep_down;
         
         //Access as str2cnt[ irrep ][ string ]; where string represents the occupation numbers of a Slater determinant; it returns a counter within the irrep block
         int ** str2cnt_up;
         int ** str2cnt_down;
         
         //Reverse functionality compared to the above arrays; cnt2str[ irrep ][ counter ]
         unsigned int ** cnt2str_up;
         unsigned int ** cnt2str_down;
         
         //For jumps[irrep_up] <= globalcounter < jumps[irrep_up+1], globalcounter corresponds to the direct product space Up( irrep_up ) x Down ( irrep_up x TargetIrrep )
         unsigned long long * jumps;
         
         //Initialize the variables above
         void Startup();
         
         //Convert between the integer and bit representations of the occupation numbers
         static void str2bits(const unsigned int L, const unsigned int string, int * bits);
         static unsigned int bits2str(const unsigned int L, int * bits);
         
         //For a certain global counter, get the corresponding irrep_up
         int getUpIrrepOfCounter(const unsigned long long counter) const;
         
         //For a given bra and ket Slater determinant return the FCI Hamiltonian matrix element; only used for GSSmallCISpace and NOT for GSDavidson
         double GetMatrixElement(int * bits_bra_up, int * bits_bra_down, int * bits_ket_up, int * bits_ket_down, int * annih_up, int * creat_up, int * annih_down, int * creat_down) const;
         
         //Convert a certain global counter into the bit representation of the occupation numbers
         void getBitsOfCounter(const unsigned long long counter, int * bits_up, int * bits_down) const;

   };

}

#endif
