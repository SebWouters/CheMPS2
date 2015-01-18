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

#include <stdlib.h>
#include <assert.h>

#include "ConvergenceScheme.h"

CheMPS2::ConvergenceScheme::ConvergenceScheme(const int nInstructions){

   this->nInstructions = nInstructions;
   
   assert( nInstructions>0 );
   nD              = new int[   nInstructions];
   fEconv          = new double[nInstructions];
   nMaxSweeps      = new int[   nInstructions];
   fNoisePrefactor = new double[nInstructions];

}

CheMPS2::ConvergenceScheme::~ConvergenceScheme(){

   delete [] nD;
   delete [] fEconv;
   delete [] nMaxSweeps;
   delete [] fNoisePrefactor;

}

int CheMPS2::ConvergenceScheme::getNInstructions(){ return nInstructions; }
         
void CheMPS2::ConvergenceScheme::setInstruction(const int instruction, const int D, const double Econv, const int nMax, const double noisePrefactor){

   assert( instruction>=0 );
   assert( instruction<nInstructions );
   assert( D>0 );
   assert( Econv>0.0 );
   assert( nMax>0 );
   assert( noisePrefactor>=0.0 );
   
                nD[instruction] = D;
            fEconv[instruction] = Econv;
        nMaxSweeps[instruction] = nMax;
   fNoisePrefactor[instruction] = noisePrefactor;

}

int CheMPS2::ConvergenceScheme::getD(const int instruction){ return nD[instruction]; }

double CheMPS2::ConvergenceScheme::getEconv(const int instruction){ return fEconv[instruction]; }

int CheMPS2::ConvergenceScheme::getMaxSweeps(const int instruction){ return nMaxSweeps[instruction]; }

double CheMPS2::ConvergenceScheme::getNoisePrefactor(const int instruction){ return fNoisePrefactor[instruction]; }


