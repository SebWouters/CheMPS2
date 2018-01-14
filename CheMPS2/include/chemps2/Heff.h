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

#ifndef HEFF_CHEMPS2_H
#define HEFF_CHEMPS2_H

#include "TensorL.h"
#include "TensorOperator.h"
#include "TensorS0.h"
#include "TensorS1.h"
#include "TensorF0.h"
#include "TensorF1.h"
#include "TensorQ.h"
#include "TensorX.h"
#include "Problem.h"
#include "SyBookkeeper.h"
#include "Sobject.h"
#include "Options.h"

namespace CheMPS2{
/** Heff class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date May 2, 2013
    
    The Heff class contains the sparse eigensolver routines and effective Hamiltonian construction needed for the DMRG algorithm. */
   class Heff{

      public:
      
         //! Constructor
         /** \param denBKIn The SyBookkeeper to get the dimensions
             \param ProbIn The Problem that contains the Hamiltonian
             \param dvdson_rtol_in The residual tolerance for the DMRG Davidson iterations */
         Heff(const SyBookkeeper * denBKIn, const Problem * ProbIn, const double dvdson_rtol_in);
         
         //! Destructor
         virtual ~Heff();
         
         //! Davidson Solver
         /** \param denS Initial guess S-object
             \param Ltensors Pointer to the single contracted 2nd quantized operators
             \param Atensors Spin-0 complementary operators of two creators
             \param Btensors Spin-1 complementary operators of two creators
             \param Ctensors Spin-0 complementary operators of a creator and an annihilator
             \param Dtensors Spin-1 complementary operators of a creator and an annihilator
             \param S0tensors Spin-0 reduction of two creators
             \param S1tensors Spin-1 reduction of two creators
             \param F0tensors Spin-0 reduction of a creator and an annihilator
             \param F1tensors Spin-1 reduction of a creator and an annihilator
             \param Qtensors Complementary operators of three sandwiched 2nd quantized operators
             \param Xtensors Pointer to the completely contracted terms
             \param nLower Number of lower-lying states to project out
             \param VeffTilde The projection operators to project the nLower lower-lying states out */
         double SolveDAVIDSON(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower = 0, double ** VeffTilde = NULL) const;
         
         //! Phase function
         /** \param TwoTimesPower Twice the power of the phase (-1)^{power}
             \return The phase (-1)^{TwoTimesPower/2} */
         static int phase(const int TwoTimesPower){ return (((TwoTimesPower/2)%2)!=0)?-1:1; }
         
      private:
      
         //The SyBookkeeper
         const SyBookkeeper * denBK;
         
         //The Problem (and hence Hamiltonian)
         const Problem * Prob;
         
         //The Davidson residual tolerance
         double dvdson_rtol;
      
         //Do Heff * memS -> memHeff
         void makeHeff(double * memS, double * memHeff, const Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const;
         
         //Fill the diagonal elements
         void fillHeffDiag(double * memHeffDiag, const Sobject * denS, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const;
         
         //Solve Davidson for the MPI_CHEMPS2_MASTER process
         double SolveDAVIDSON_main(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const;
         
         //Solve Davidson for the helper processes
         double SolveDAVIDSON_help(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const;
         
         //The diagrams: Type 1/5
         void addDiagram1A(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorX * Xleft) const;
         void addDiagram1B(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorX * Xright) const;
         void addDiagram1C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, double Helem_links) const;
         void addDiagram1D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, double Helem_rechts) const;
         void addDiagramExcitations(const int ikappa, double * memS, double * memHeff, const Sobject * denS, int nLower, double ** VeffTilde) const;
         
         //The diagrams: Type 2/5
         void addDiagram2a1spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Atensors, TensorS0 **** S0tensors, double * workspace) const;
         void addDiagram2a2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Atensors, TensorS0 **** S0tensors, double * workspace) const;
         void addDiagram2a1spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Btensors, TensorS1 **** S1tensors, double * workspace) const;
         void addDiagram2a2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Btensors, TensorS1 **** S1tensors, double * workspace) const;
         void addDiagram2a3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Ctensors, TensorF0 **** F0tensors, double * workspace) const;
         void addDiagram2a3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Dtensors, TensorF1 **** F1tensors, double * workspace) const;
         void addDiagram2b1and2b2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const;
         void addDiagram2c1and2c2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const;
         void addDiagram2dall(const int ikappa, double * memS, double * memHeff, const Sobject * denS) const;
         void addDiagram2e1and2e2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const;
         void addDiagram2f1and2f2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const;
         void addDiagram2b3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagram2c3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagram2e3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagram2f3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagram2b3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagram2c3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagram2e3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagram2f3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const;
         
         //The diagrams: Type 3/5
         void addDiagram3Aand3D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qleft, TensorL ** Lleft, double * temp) const;
         void addDiagram3Band3I(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qleft, TensorL ** Lleft, double * temp) const;
         void addDiagram3C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ ** Qleft, TensorL ** Lright, double * temp) const;
         void addDiagram3Eand3H(const int ikappa, double * memS, double * memHeff, const Sobject * denS) const;
         void addDiagram3Kand3F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qright, TensorL ** Lright, double * temp) const;
         void addDiagram3Land3G(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qright, TensorL ** Lright, double * temp) const;
         void addDiagram3J(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ ** Qright, TensorL ** Lleft, double * temp) const;
         
         //The diagrams: Type 4/5
         void addDiagram4A1and4A2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atens) const;
         void addDiagram4A1and4A2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Btens) const;
         void addDiagram4A3and4A4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctens) const;
         void addDiagram4A3and4A4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtens) const;
         void addDiagram4B1and4B2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Aleft, TensorL ** Lright, double * temp) const;
         void addDiagram4B1and4B2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Bleft, TensorL ** Lright, double * temp) const;
         void addDiagram4B3and4B4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Cleft, TensorL ** Lright, double * temp) const;
         void addDiagram4B3and4B4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Dleft, TensorL ** Lright, double * temp) const;
         void addDiagram4C1and4C2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Aleft, TensorL ** Lright, double * temp) const;
         void addDiagram4C1and4C2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Bleft, TensorL ** Lright, double * temp) const;
         void addDiagram4C3and4C4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Cleft, TensorL ** Lright, double * temp) const;
         void addDiagram4C3and4C4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Dleft, TensorL ** Lright, double * temp) const;
         void addDiagram4D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, double * temp) const;
         void addDiagram4E(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram4F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lright, double * temp) const;
         void addDiagram4G(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lright, double * temp) const;
         void addDiagram4H(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram4I(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, double * temp) const;
         void addDiagram4J1and4J2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Aright) const;
         void addDiagram4J1and4J2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Bright) const;
         void addDiagram4J3and4J4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Cright) const;
         void addDiagram4J3and4J4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dright) const;
         void addDiagram4K1and4K2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Aright, double * temp) const;
         void addDiagram4L1and4L2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Aright, double * temp) const;
         void addDiagram4K1and4K2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Bright, double * temp) const;
         void addDiagram4L1and4L2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Bright, double * temp) const;
         void addDiagram4K3and4K4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Cright, double * temp) const;
         void addDiagram4L3and4L4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Cright, double * temp) const;
         void addDiagram4K3and4K4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Dright, double * temp) const;
         void addDiagram4L3and4L4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Dright, double * temp) const;
         
         //The diagrams: type 5/5
         void addDiagram5A(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram5B(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram5C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram5D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram5E(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         void addDiagram5F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const;
         
         //All diagonal contributions
         void addDiagonal1A(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorX * Xleft) const;
         void addDiagonal1B(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorX * Xright) const;
         void addDiagonal1C(const int ikappa, double * memHeffDiag, const Sobject * denS, const double Helem_links) const;
         void addDiagonal1D(const int ikappa, double * memHeffDiag, const Sobject * denS, const double Helem_rechts) const;
         void addDiagonal2d3all(const int ikappa, double * memHeffDiag, const Sobject * denS) const;
         void addDiagonal2b3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagonal2c3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagonal2e3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagonal2f3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const;
         void addDiagonal2b3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagonal2c3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagonal2e3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagonal2f3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const;
         void addDiagonal2a3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator **** Ctensors, TensorF0 **** F0tensors) const;
         void addDiagonal2a3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator **** Dtensors, TensorF1 **** F1tensors) const;
         void addDiagonalExcitations(const int ikappa, double * memHeffDiag, const Sobject * denS, int nLower, double ** VeffTilde) const;
         
   };
}

#endif
