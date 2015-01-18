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

#include "TensorL.h"
#include "Lapack.h"
#include "Gsl.h"

CheMPS2::TensorL::TensorL(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : TensorSwap(indexIn, IdiffIn, movingRightIn, denBKIn){

}

CheMPS2::TensorL::~TensorL(){

}

void CheMPS2::TensorL::makenew(TensorT * denT){

   if (movingRight){ makenewRight(denT); }
   else{ makenewLeft( denT); }
   
}

void CheMPS2::TensorL::makenewRight(TensorT * denT){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int IDR = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUR = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDR = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], IDR             );
      
      for (int geval=0; geval<2; geval++){
         int NL,TwoSL,IL;
         switch(geval){
            case 0:
               NL = sectorN1[ikappa];
               TwoSL = sectorTwoS1[ikappa];
               IL = sectorI1[ikappa];
               break;
            case 1:
               NL = sectorN1[ikappa] - 1;
               TwoSL = sectorTwoSD[ikappa];
               IL = IDR;
               break;
         }
         int dimL = denBK->gCurrentDim(index-1, NL, TwoSL, IL);
         if (dimL>0){
         
            double * BlockTup   = denT->gStorage(NL,TwoSL,IL,sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(NL,TwoSL,IL,sectorN1[ikappa]+1,sectorTwoSD[ikappa],IDR);
            
            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if (geval>=1){
               int fase = ((((sectorTwoSD[ikappa] - sectorTwoS1[ikappa] + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
            }
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);
         
         }
      }
   } 

}


void CheMPS2::TensorL::makenewLeft(TensorT * denT){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int IDL = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUL = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDL = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], IDL             );
      
      for (int geval=0; geval<2; geval++){
         int NR,TwoSR,IR;
         switch(geval){
            case 0:
               NR = sectorN1[ikappa]+1;
               TwoSR = sectorTwoSD[ikappa];
               IR = IDL;
               break;
            case 1:
               NR = sectorN1[ikappa]+2;
               TwoSR = sectorTwoS1[ikappa];
               IR = sectorI1[ikappa];
               break;
         }
         int dimR = denBK->gCurrentDim(index+1, NR, TwoSR, IR);
         if (dimR>0){
         
            double * BlockTup   = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], NR, TwoSR, IR);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], IDL,              NR, TwoSR, IR);
            
            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if (geval>=1){
               int fase = ((((sectorTwoS1[ikappa] - sectorTwoSD[ikappa] + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
            }
            double beta = 1.0; //add
            dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimR,&alpha,BlockTup,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
         }
      }
   }

}


