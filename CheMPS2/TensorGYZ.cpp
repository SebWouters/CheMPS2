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

#include <math.h>

#include "TensorGYZ.h"
#include "Lapack.h"

CheMPS2::TensorGYZ::TensorGYZ(const int indexIn, const char identityIn, const SyBookkeeper * denBKIn) : TensorDiag(indexIn, denBKIn){

   identity = identityIn;

}

CheMPS2::TensorGYZ::~TensorGYZ(){ }

void CheMPS2::TensorGYZ::construct(TensorT * denT){

   for (int ikappa=0; ikappa<nKappa; ikappa++){

      int NL       = -1;
      int IL       = -1;
      int TwoSL    = -1;
      double alpha = 1.0;

      if (identity=='Y'){
         NL    = sectorN1[ikappa];
         TwoSL = sectorTwoS1[ikappa];
         IL    = sectorI1[ikappa];
      }

      if (identity=='Z'){
         NL    = sectorN1[ikappa]-2;
         TwoSL = sectorTwoS1[ikappa];
         IL    = sectorI1[ikappa];
      }

      if (identity=='G'){
         NL    = sectorN1[ikappa]-1;
         TwoSL = sectorTwoS1[ikappa]-1;
         IL    = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
         alpha = sqrt(0.5);
      }

      int dimR = denBK->gCurrentDim(index,   sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimL = denBK->gCurrentDim(index-1, NL,               TwoSL,               IL);

      if (dimL>0){
         double * BlockT = denT->gStorage(NL, TwoSL, IL, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
         char trans = 'T';
         char notr = 'N';
         double beta = 0.0;
         dgemm_(&trans,&notr,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
      } else {
         for (int cnt=kappa2index[ikappa]; cnt<kappa2index[ikappa+1]; cnt++){ storage[cnt] = 0.0; }
      }

      if (identity=='G'){
         TwoSL = sectorTwoS1[ikappa]+1;
         dimL  = denBK->gCurrentDim(index-1, NL, TwoSL, IL);

         if (dimL>0){
            double * BlockT = denT->gStorage(NL, TwoSL, IL, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
            char trans = 'T';
            char notr = 'N';
            double beta = 1.0; //ADD NOW!!!
            dgemm_(&trans,&notr,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
         }
      }

   }

}

void CheMPS2::TensorGYZ::update(TensorT * denT, TensorGYZ * denGYZ, double * workmemLR){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      for (int cnt=kappa2index[ikappa]; cnt<kappa2index[ikappa+1]; cnt++){ storage[cnt] = 0.0; } //Clear
      updateRight(ikappa, denT, denGYZ, workmemLR); //Add the previous TensorGYZ
   }

}


