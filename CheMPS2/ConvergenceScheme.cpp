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

#include <stdlib.h>
#include <assert.h>

#include "ConvergenceScheme.h"

CheMPS2::ConvergenceScheme::ConvergenceScheme( const int num_instructions ){

   this->num_instructions = num_instructions;

   assert( num_instructions > 0 );
   num_D              = new    int[ num_instructions ];
   energy_convergence = new double[ num_instructions ];
   num_max_sweeps     = new    int[ num_instructions ];
   noise_prefac       = new double[ num_instructions ];
   dvdson_rtol        = new double[ num_instructions ];

}

CheMPS2::ConvergenceScheme::~ConvergenceScheme(){

   delete [] num_D;
   delete [] energy_convergence;
   delete [] num_max_sweeps;
   delete [] noise_prefac;
   delete [] dvdson_rtol;

}

int CheMPS2::ConvergenceScheme::get_number() const{ return num_instructions; }
         
void CheMPS2::ConvergenceScheme::set_instruction( const int instruction, const int D, const double energy_conv, const int max_sweeps, const double noise_prefactor, const double davidson_rtol ){

   assert( instruction >= 0 );
   assert( instruction < num_instructions );
   assert( D > 0 );
   assert( energy_conv > 0.0 );
   assert( max_sweeps > 0 );
   assert( davidson_rtol > 0.0 );
   
                num_D[ instruction ] = D;
   energy_convergence[ instruction ] = energy_conv;
       num_max_sweeps[ instruction ] = max_sweeps;
         noise_prefac[ instruction ] = noise_prefactor;
          dvdson_rtol[ instruction ] = davidson_rtol;

}

int CheMPS2::ConvergenceScheme::get_D( const int instruction ) const{ return num_D[ instruction ]; }

double CheMPS2::ConvergenceScheme::get_energy_conv( const int instruction ) const{ return energy_convergence[ instruction ]; }

int CheMPS2::ConvergenceScheme::get_max_sweeps( const int instruction ) const{ return num_max_sweeps[ instruction ]; }

double CheMPS2::ConvergenceScheme::get_noise_prefactor( const int instruction ) const{ return noise_prefac[ instruction ]; }

double CheMPS2::ConvergenceScheme::get_dvdson_rtol( const int instruction ) const{ return dvdson_rtol[ instruction ]; }


