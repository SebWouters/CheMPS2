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
#include <math.h>

#include "TensorD.h"
#include "Lapack.h"

CheMPS2::TensorD::TensorD(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : TensorF1Dbase(indexIn, IdiffIn, movingRightIn, denBKIn){

}

CheMPS2::TensorD::~TensorD(){

}

void CheMPS2::TensorD::ClearStorage(){ Clear(); }
      
void CheMPS2::TensorD::AddATerm(double alpha, TensorF1Dbase * TermToAdd){

   int inc = 1;
   daxpy_(kappa2index+nKappa, &alpha, TermToAdd->gStorage(), &inc, storage, &inc);

}

void CheMPS2::TensorD::AddATermTranspose(const double alpha, TensorF1Dbase * TermToAdd){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int fase = ((((sectorTwoS1[ikappa]-sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
      const double prefactor = alpha * fase * ((movingRight)?sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0)):sqrt((sectorTwoSD[ikappa]+1.0)/(sectorTwoS1[ikappa]+1.0)));
      int dimU = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
      int ID = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimD = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoSD[ikappa], ID);
      
      double * BlockToAdd = TermToAdd->gStorage(sectorN1[ikappa], sectorTwoSD[ikappa], ID, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
      for (int irow=0; irow<dimU; irow++){
         for (int icol=0; icol<dimD; icol++){
            storage[kappa2index[ikappa] + irow + dimU * icol] += prefactor * BlockToAdd[icol + dimD * irow];
         }
      }
   }

}

