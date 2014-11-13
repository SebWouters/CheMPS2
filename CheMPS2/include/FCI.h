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
    
    The FCI class performs the full configuration interaction ground state calculation in a given particle number, irrep, and spin projection sector of a given Hamiltonian. It also contains the functionality to calculate Green's functions.
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
         
//==========> Get basic information on this FCI object
         
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
         
//==========> The core routines which you're most likely going actually going to use
         
         //! Calculates the FCI ground state with Davidson's algorithm
         /** \param inoutput If inoutput!=NULL, vector with getVecLength() variables which contains on exit the solution of the FCI calculation
             \return The ground state energy */
         double GSDavidson(double * inoutput=NULL) const;
         
         //! Calculates the CI ground state in a small space of low-energy Slater determinants with lapack's dsyev
         /** \param Nslaters The number of low-energy Slater determinants to consider
             \param vec Vector with getVecLength() variables which contains on exit the solution of the small CI calculation
             \return The small CI ground state energy */
         double GSSmallCISpace(const unsigned int Nslaters, double * vec) const;
         
         //! Calculates the solution of the equation ( alpha + beta * Hamiltonian + I * eta ) Solution = RHS with conjugate gradient
         /** \param alpha The real part of the scalar in the operator
             \param beta The prefactor of the Hamiltonian in the operator
             \param eta The imaginary part of the scalar in the operator
             \param RHS The real-valued right-hand side of the equation with length getVecLength()
             \param RealSol If not NULL, on exit this array of length getVecLength() contains the real part of the solution
             \param ImagSol If not NULL, on exit this array of length getVecLength() contains the imaginary part of the solution
             \param checkError If true, at the end the RMS error without preconditioner will be printed */
         void CGSolveSystem(const double alpha, const double beta, const double eta, double * RHS, double * RealSol, double * ImagSol, const bool checkError=true) const;
         
         //! Set this FCI vector to an operator acting on a previous FCI vector
         /** \param whichOperator With which operator should be acted on the other FCI state: C means creator, A means annihilator, N means particle number
             \param isUp Boolean which denotes if the operator corresponds to an up (alpha) or down (beta) electron
             \param orbIndex Orbital index on which the operator acts
             \param thisVector Vector with length getVecLength() where the result of the operation should be stored
             \param otherFCI FCI instance which corresponds to the FCI vector on which is acted
             \param otherVector Vector with length otherFCI->getVecLength() which contains the FCI vector on which is acted */
         void ActWithSecondQuantizedOperator(const char whichOperator, const bool isUp, const unsigned int orbIndex, double * thisVector, FCI * otherFCI, double * otherVector);
         
         //! Construct the (spin-summed) 2-RDM of a FCI vector
         /** \param vector The FCI vector
             \param TwoRDM To store the 2-RDM; needs to be of size L^4 (L = number of orbitals); irrep symmetry shows in 2-RDM elements being zero; physics notation (see gHmat)
             \return The energy of the given FCI vector, calculated by contraction of the 2-RDM with Hmat */
         double Fill2RDM(double * vector, double * TwoRDM) const;
         
         //! Measure S(S+1) (spin squared)
         /** \param vector The FCI vector
             \return Measured value of S(S+1) */
         double CalcSpinSquared(double * vector) const;
         
//==========> Functions which involve bitstring representations
         
         //! Function which returns a FCI coefficient
         /** \param bits_up The bit string representation of the up or alpha electron Slater determinant
             \param bits_down The bit string representation of the down or beta electron Slater determinant
             \param vector The FCI vector with getVecLength() variables from which a coefficient is desired
             \return The corresponding FCI coefficient; 0.0 if something is not OK */
         double getFCIcoeff(int * bits_up, int * bits_down, double * vector) const;
         
         //! Find the bit representation of a FCI index
         /** \param counter The given FCI index
             \param bits_up Array of length L to store the bit representation of the up (alpha) electrons in
             \param bits_down Array of length L to store the bit representation of the down (beta) electrons in */
         void getBitsOfCounter(const unsigned long long counter, int * bits_up, int * bits_down) const;
         
         //! Convertor between two representations of a same spin-projection Slater determinant
         /** \param L The number of orbitals
             \param bitstring The input integer, whos bits are the occupation numbers of the orbitals
             \param bits Contains on exit the L bits of bitstring */
         static void str2bits(const unsigned int L, const unsigned int bitstring, int * bits);
         
         //! Convertor between two representations of a same spin-projection Slater determinant
         /** \param L The number of orbitals
             \param bits Array with the L bits which should be combined to a single integer
             \return The single integer which represents the bits in bits */
         static unsigned int bits2str(const unsigned int L, int * bits);
         
         //! Find the irrep of the up Slater determinant to which a certain FCI coefficient belongs
         /** \param counter The FCI coefficient index
             \return The corresponding irrep of the up Slater determinant */
         int getUpIrrepOfCounter(const unsigned long long counter) const;

//==========> Functions involving the FCI Hamiltonian, and its elements

         //! Sandwich the Hamiltonian between two Slater determinants (return a specific element) (without Econstant!!)
         /** \param bits_bra_up Bit representation of the <bra| Slater determinant of the up (alpha) electrons
             \param bits_bra_down Bit representation of the <bra| Slater determinant of the down (beta) electrons
             \param bits_ket_up Bit representation of the |ket> Slater determinant of the up (alpha) electrons
             \param bits_ket_down Bit representation of the |ket> Slater determinant of the down (beta) electrons
             \param work1 Work array of length 2
             \param work2 Work array of length 2
             \param work3 Work array of length 2
             \param work4 Work array of length 2
             \return The FCI Hamiltonian element which connects the given two Slater determinants */
         double GetMatrixElement(int * bits_bra_up, int * bits_bra_down, int * bits_ket_up, int * bits_ket_down, int * work1, int * work2, int * work3, int * work4) const;
      
         //! Function which performs the Hamiltonian times Vector product (without Econstant!!)
         /** \param input The vector on which the Hamiltonian should act
             \param output The Hamiltonian times input on exit */
         void HamTimesVec(double * input, double * output) const;
         
         //! Function which returns the diagonal elements of the FCI Hamiltonian (without Econstant!!) = Slater determinant energies
         /** \param diag Vector with getVecLength() variables which contains on exit the diagonal elements of the FCI Hamiltonian */
         void DiagHam(double * diag) const;
         
         //! Function which returns the diagonal elements of the FCI Hamiltonian squared (without Econstant!!)
         /** \param output Vector with getVecLength() variables which contains on exit the diagonal elements of the FCI Hamiltonian squared */
         void DiagHamSquared(double * output) const;

//==========> Some lapack like routine's which work with unsigned long long's

         //! Fill a vector with random numbers in the interval [-1,1[; used when output for GSDavidson is desired but no specific input can be given (as for the Hubbard model in the site basis)
         /** \param vecLength The length of the vector; normally getVecLength()
             \param vec The vector to fill with random numbers */
         static void FillRandom(const unsigned long long vecLength, double * vec);

         //! Take the inproduct of two very long vectors
         /** \param vecLength The vector length
             \param vec1 The first vector
             \param vec2 The second vector
             \return The inproduct < vec1 | vec2 > */
         static double FCIddot(const unsigned long long vecLength, double * vec1, double * vec2);
         
         //! Clear a very long vector
         /** \param vecLength The vector length
             \param vec The vector which has to be set to zero */
         static void FCIdclear(const unsigned long long vecLength, double * vec);
         
         //! Copy a very long vector
         /** \param vecLength The vector length
             \param origin Vector to be copied
             \param target Where to copy the vector to */
         static void FCIdcopy(const unsigned long long vecLength, double * origin, double * target);
         
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
         
      protected:
         
         //! Calculate the solution of the system Operator |Sol> = |RESID> with Operator = precon * [ ( alpha + beta * H )^2 + eta^2 ] * precon
         /** \param alpha The parameter alpha of the operator
             \param beta The parameter beta of the operator
             \param eta The parameter eta of the operator
             \param precon The diagonal preconditioner
             \param Sol On exit the solution is stored in this array of size getVecLength()
             \param RESID On entry the RHS of the equation is stored in this array of size getVecLength(), on exit is is overwritten
             \param PVEC Workspace of size getVecLength()
             \param OxPVEC Workspace of size getVecLength()
             \param temp Workspace of size getVecLength()
             \param temp2 Workspace of size getVecLength() */
         void CGCoreSolver(const double alpha, const double beta, const double eta, double * precon, double * Sol, double * RESID, double * PVEC, double * OxPVEC, double * temp, double * temp2) const;

         //! Calculate out = (alpha + beta * Hamiltonian) * in (Econstant is taken into account!!)
         /** \param alpha The parameter alpha of the operator
             \param beta The parameter beta of the operator
             \param in On entry the vector of size getVecLength() on which the operator should be applied; unchanged on exit
             \param out Array of size getVecLength(), which contains on exit (alpha + beta * Hamiltonian) * in */
         void CGAlphaPlusBetaHAM(const double alpha, const double beta, double * in, double * out) const;
         
         //! Calculate out = precon * [(alpha + beta * Hamiltonian)^2 + eta^2] * precon * in (Econstant is taken into account!!)
         /** \param alpha The parameter alpha of the operator
             \param beta The parameter beta of the operator
             \param eta The parameter eta of the operator
             \param precon The diagonal preconditioner
             \param in On entry the vector of size getVecLength() on which the operator should be applied; unchanged on exit
             \param out Array of size getVecLength(), which contains on exit precon * [(alpha + beta * Hamiltonian)^2 + eta^2] * precon * in
             \param temp Workspace of size getVecLength()
             \param temp2 Workspace of size getVecLength() */
         void CGOperator(const double alpha, const double beta, const double eta, double * precon, double * in, double * temp, double * temp2, double * out) const;
         
         //! Calculates (without approximation) precon = (diag[(alpha + beta * Hamiltonian)^2 + eta^2])^{ -1/2 } (Econstant is taken into account!!)
         /** \param alpha The parameter alpha of the operator
             \param beta The parameter beta of the operator
             \param eta The parameter eta of the operator
             \param precon Array of size getVecLength(), which contains on exit (diag[(alpha + beta * Hamiltonian)^2 + eta^2])^{ -1/2 }
             \param workspace Workspace of size getVecLength() */
         void CGDiagPrecond(const double alpha, const double beta, const double eta, double * precon, double * workspace) const;
      
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

   };

}

#endif
