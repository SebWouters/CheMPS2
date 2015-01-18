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

#ifndef CONVERGENCESCHEME_CHEMPS2_H
#define CONVERGENCESCHEME_CHEMPS2_H

namespace CheMPS2{
/** ConvergenceScheme class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 7, 2013
    
    The ConvergenceScheme class contains the convergence settings. This is a list of instructions, which are performed in order. Each instruction line contains the information for one particular batch of DMRG sweeps:\n
    (1) the number of renormalized basis states to keep (D)\n
    (2) the energy convergence threshold for energy changes per left- and right-sweep\n
    (3) the maximum number of iterations, in case the energy changes do not drop below the threshold\n
    (4) the noise prefactor f\n
    \n
    The noise level which is added to the Sobject is the product of\n
    (1) f\n
    (2) the maximum discarded weight during the last sweep\n
    (3) a random number in the interval [-0.5,0.5]*/
   class ConvergenceScheme{

      public:
      
         //! Constructor
         /** \param nInstructions the number of instructions */
         ConvergenceScheme(const int nInstructions);
         
         //! Destructor
         virtual ~ConvergenceScheme();
         
         //!Get the number of instructions
         /** return the number of instructions to converge the DMRG calculation */
         int getNInstructions();
         
         //! Set an instruction
         /** \param instruction the number of the instruction
             \param D the number of renormalized states for that instruction
             \param Econv the energy convergence threshold for that instruction
             \param nMax the max. number of sweeps for that instruction
             \param noisePrefactor the noise prefactor for that instruction */
         void setInstruction(const int instruction, const int D, const double Econv, const int nMax, const double noisePrefactor);
         
         //! Get the number of renormalized states for a particular instruction
         /** \param instruction the number of the instruction
             \return the number of renormalized states for this instruction */
         int getD(const int instruction);
         
         //! Get the energy convergence threshold for a particular instruction
         /** \param instruction the number of the instruction
             \return the energy convergence threshold for this instruction */
         double getEconv(const int instruction);
         
         //! Get the maximum number of sweeps for a particular instruction
         /** \param instruction the number of the instruction
             \return the maximum number of sweeps for this instruction */
         int getMaxSweeps(const int instruction);
         
         //! Get the noise prefactor for a particular instruction
         /** \param instruction the number of the instruction
             \return the noise prefactor for this instruction */
         double getNoisePrefactor(const int instruction);
         
      private:
      
         //The number of instructions
         int nInstructions;
      
         //Number of renormalized states
         int * nD;
         
         //Energy convergence thresholds
         double * fEconv;
         
         //Maximum number of sweeps per instruction
         int * nMaxSweeps;
         
         //The noise prefactor for each instruction
         double * fNoisePrefactor;
         
   };
}

#endif
