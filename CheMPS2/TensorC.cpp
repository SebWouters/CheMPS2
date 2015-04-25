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

#include "TensorC.h"
#include "Lapack.h"

CheMPS2::TensorC::TensorC(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : TensorF0Cbase(indexIn, IdiffIn, movingRightIn, denBKIn){

}

CheMPS2::TensorC::~TensorC(){

}

void CheMPS2::TensorC::ClearStorage(){ Clear(); }

void CheMPS2::TensorC::AddATerm(double alpha, TensorF0Cbase * TermToAdd){

   int inc = 1;
   daxpy_(kappa2index+nKappa, &alpha, TermToAdd->gStorage(), &inc, storage, &inc);

}

void CheMPS2::TensorC::AddATermTranspose(const double alpha, TensorF0Cbase * TermToAdd){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimU = denBK->gCurrentDim(index,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
      const int ID = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimD = denBK->gCurrentDim(index,sectorN1[ikappa],sectorTwoS1[ikappa],ID);
      
      double * BlockToAdd = TermToAdd->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],ID,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
      for (int irow=0; irow<dimU; irow++){
         for (int icol=0; icol<dimD; icol++){
            storage[kappa2index[ikappa] + irow + dimU * icol] += alpha * BlockToAdd[icol + dimD * irow];
         }
      }
   }

}

