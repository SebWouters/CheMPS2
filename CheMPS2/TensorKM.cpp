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

#include "TensorKM.h"
#include "Lapack.h"

CheMPS2::TensorKM::TensorKM(const int indexIn, const char identityIn, const int IdiffIn, const SyBookkeeper * denBKIn) :
TensorOperator(indexIn,
               1, // two_j
               1, // n_elec
               IdiffIn,
               true, // TensorKM only exists moving left to right
               true, // prime_last
               false, // No jw_phase when updating (two-orbital mutual information!)
               denBKIn){

   identity = identityIn;

}

CheMPS2::TensorKM::~TensorKM(){ }

void CheMPS2::TensorKM::construct(TensorT * denT){

   clear();
   
   if ( identity == 'K' ){

      for (int ikappa=0; ikappa<nKappa; ikappa++){

         const int IDR = Irreps::directProd( n_irrep, sectorI1[ikappa] );
         int dimUR = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa],    sectorI1[ikappa]);
         int dimDR = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sector_2S_down[ikappa], IDR             );
         int dimL = denBK->gCurrentDim(index-1, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);

         if (dimL>0){

            double * BlockTup   = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa],   sectorTwoS1[ikappa],    sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+1, sector_2S_down[ikappa], IDR);

            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);

         }
      }
   }
   
   if ( identity == 'M' ){
   
      for (int ikappa=0; ikappa<nKappa; ikappa++){

         const int IDR = Irreps::directProd( n_irrep, sectorI1[ikappa] );
         int dimUR = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa],    sectorI1[ikappa]);
         int dimDR = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sector_2S_down[ikappa], IDR             );
         int dimL = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1,  sector_2S_down[ikappa], IDR);

         if (dimL>0){

            double * BlockTup   = denT->gStorage(sectorN1[ikappa]-1, sector_2S_down[ikappa], IDR, sectorN1[ikappa],   sectorTwoS1[ikappa],    sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]-1, sector_2S_down[ikappa], IDR, sectorN1[ikappa]+1, sector_2S_down[ikappa], IDR);

            char trans = 'T';
            char notrans = 'N';
            int fase = ((((sector_2S_down[ikappa] - sectorTwoS1[ikappa] + 1)/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sectorTwoS1[ikappa]+1.0)/(sector_2S_down[ikappa]+1));
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);

         }
      }
   }

}


