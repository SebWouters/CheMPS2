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

#ifndef CONVERGENCESCHEME_CHEMPS2_H
#define CONVERGENCESCHEME_CHEMPS2_H

#include "Options.h"

namespace CheMPS2{
/** ConvergenceScheme class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date November 7, 2013

    The ConvergenceScheme class contains the convergence settings. This is a list of instructions, which are performed in order. Each instruction line contains the information for one particular batch of DMRG sweeps:\n
    (1) the number of renormalized basis states to keep (D)\n
    (2) the energy convergence threshold for energy changes per left- and right-sweep\n
    (3) the maximum number of iterations, in case the energy changes do not drop below the threshold\n
    (4) the noise prefactor f\n
    (5) the Davidson residual tolerance\n
    \n
    The noise level which is added to the Sobject is the product of\n
    (1) f\n
    (2) the maximum discarded weight during the last sweep\n
    (3) a random number in the interval [-0.5,0.5]*/
   class ConvergenceScheme{

      public:

         //! Constructor
         /** \param num_instructions the number of instructions */
         ConvergenceScheme(const int num_instructions);

         //! Destructor
         virtual ~ConvergenceScheme();

         //! Get the number of instructions
         /** return the number of instructions */
         int get_number() const;

         //! Set an instruction
         /** \param instruction the number of the instruction
             \param D the number of renormalized states for that instruction
             \param energy_conv the energy convergence threshold for that instruction
             \param max_sweeps the max. number of sweeps for that instruction
             \param noise_prefactor the noise prefactor for that instruction
             \param davidson_rtol the Davidson residual tolerance for that instruction */
         void set_instruction(const int instruction, const int D, const double energy_conv, const int max_sweeps, const double noise_prefactor, const double davidson_rtol);
         
         //! Set an instruction
         /** \param instruction the number of the instruction
             \param D the number of renormalized states for that instruction
             \param energy_conv the energy convergence threshold for that instruction
             \param max_sweeps the max. number of sweeps for that instruction
             \param noise_prefactor the noise prefactor for that instruction */
         void setInstruction(const int instruction, const int D, const double energy_conv, const int max_sweeps, const double noise_prefactor){
            set_instruction( instruction, D, energy_conv, max_sweeps, noise_prefactor, CheMPS2::DAVIDSON_DMRG_RTOL );
         }

         //! Get the number of renormalized states for a particular instruction
         /** \param instruction the number of the instruction
             \return the number of renormalized states for this instruction */
         int get_D(const int instruction) const;

         //! Get the energy convergence threshold for a particular instruction
         /** \param instruction the number of the instruction
             \return the energy convergence threshold for this instruction */
         double get_energy_conv(const int instruction) const;

         //! Get the maximum number of sweeps for a particular instruction
         /** \param instruction the number of the instruction
             \return the maximum number of sweeps for this instruction */
         int get_max_sweeps(const int instruction) const;

         //! Get the noise prefactor for a particular instruction
         /** \param instruction the number of the instruction
             \return the noise prefactor for this instruction */
         double get_noise_prefactor(const int instruction) const;

         //! Get the Davidson residual tolerance for a particular instruction
         /** \param instruction the number of the instruction
             \return the Davidson residual tolerance for this instruction */
         double get_dvdson_rtol(const int instruction) const;

      private:

         //The number of instructions
         int num_instructions;

         //Number of renormalized states
         int * num_D;

         //Energy convergence thresholds
         double * energy_convergence;

         //Maximum number of sweeps per instruction
         int * num_max_sweeps;

         //The noise prefactor for each instruction
         double * noise_prefac;

         //The Davidson residual tolerance for each instruction
         double * dvdson_rtol;

   };
}

#endif
