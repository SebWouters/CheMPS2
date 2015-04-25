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

#include "TensorS0.h"
#include "Lapack.h"

CheMPS2::TensorS0::TensorS0(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : TensorS0Abase(indexIn, IdiffIn, movingRightIn, denBKIn){

}

CheMPS2::TensorS0::~TensorS0(){

}

void CheMPS2::TensorS0::makenew(TensorT * denT){

   if (movingRight){ makenewRight(denT); }
   else{ makenewLeft( denT); }

}

void CheMPS2::TensorS0::makenew(TensorL * denL, TensorT * denT, double * workmem){

   if (movingRight){ makenewRight(denL, denT, workmem); }
   else{ makenewLeft( denL, denT, workmem); }
       
}

void CheMPS2::TensorS0::makenewRight(TensorT * denT){ //Idiff = Itrivial

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUR = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDR = denBK->gCurrentDim(index,   sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimL  = denBK->gCurrentDim(index-1, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      if (dimL>0){
         
         double * BlockTup   = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         double * BlockTdown = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);

         char trans = 'T';
         char notrans = 'N';
         double alpha = sqrt(2.0);
         double beta = 1.0; //add
         dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);
         
      }
   }

}

void CheMPS2::TensorS0::makenewLeft(TensorT * denT){ //Idiff = Itrivial

   Clear();

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUL = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDL = denBK->gCurrentDim(index,   sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimR  = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
      if (dimR>0){
         
         double * BlockTup   = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         double * BlockTdown = denT->gStorage(sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);

         char trans = 'T';
         char notrans = 'N';
         double alpha = sqrt(2.0);
         double beta = 1.0; //add
         dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimR,&alpha,BlockTup,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
      }
   }

}

void CheMPS2::TensorS0::makenewRight(TensorL * denL, TensorT * denT, double * workmem){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDR = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUR = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDR = denBK->gCurrentDim(index, sectorN1[ikappa]+2, sectorTwoS1[ikappa], IDR             );
      
      for (int geval=0; geval<4; geval++){
         int NLU,TwoSLU,ILU,TwoSLD,ILD; //NLD = NLU+1
         switch(geval){
            case 0:
               NLU = sectorN1[ikappa];
               TwoSLU = sectorTwoS1[ikappa];
               ILU = sectorI1[ikappa];
               TwoSLD = sectorTwoS1[ikappa]-1;
               ILD = Irreps::directProd( ILU, denL->gIdiff() );
               break;
            case 1:
               NLU = sectorN1[ikappa];
               TwoSLU = sectorTwoS1[ikappa];
               ILU = sectorI1[ikappa];
               TwoSLD = sectorTwoS1[ikappa]+1;
               ILD = Irreps::directProd( ILU, denL->gIdiff() );
               break;
            case 2:
               NLU = sectorN1[ikappa]-1;
               TwoSLU = sectorTwoS1[ikappa]-1;
               ILU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
               TwoSLD = sectorTwoS1[ikappa];
               ILD = IDR;
               break;
            case 3:
               NLU = sectorN1[ikappa]-1;
               TwoSLU = sectorTwoS1[ikappa]+1;
               ILU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
               TwoSLD = sectorTwoS1[ikappa];
               ILD = IDR;
               break;
         }
         int dimLU = denBK->gCurrentDim(index-1, NLU,   TwoSLU, ILU);
         int dimLD = denBK->gCurrentDim(index-1, NLU+1, TwoSLD, ILD);
         if ((dimLU>0) && (dimLD>0)){
         
            double * BlockTup   = denT->gStorage(NLU,   TwoSLU, ILU, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(NLU+1, TwoSLD, ILD, sectorN1[ikappa]+2, sectorTwoS1[ikappa], IDR);
            double * BlockL     = denL->gStorage(NLU,   TwoSLU, ILU, NLU+1,              TwoSLD,              ILD);
            
            //factor * Tup^T * L -> mem
            char trans = 'T';
            char notrans = 'N';
            double alpha;
            if (geval<=1){
               int fase = ((((sectorTwoS1[ikappa] - TwoSLD + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(0.5 * (TwoSLD+1.0) / (sectorTwoS1[ikappa]+1.0) );
            } else {
               alpha = - sqrt(0.5);
            }
            double beta = 0.0; //set
            dgemm_(&trans,&notrans,&dimUR,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,BlockL,&dimLU,&beta,workmem,&dimUR);
            
            //mem * Tdown -> storage
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&notrans,&dimUR,&dimDR,&dimLD,&alpha,workmem,&dimUR,BlockTdown,&dimLD,&beta,storage+kappa2index[ikappa],&dimUR);
         
         }
      }
   } 

}

void CheMPS2::TensorS0::makenewLeft(TensorL * denL, TensorT * denT, double * workmem){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDL = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUL = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDL = denBK->gCurrentDim(index, sectorN1[ikappa]+2, sectorTwoS1[ikappa], IDL             );
      
      for (int geval=0; geval<4; geval++){
         int NRU,TwoSRU,IRU,TwoSRD,IRD; //NRD = NRU+1
         switch(geval){
            case 0:
               NRU = sectorN1[ikappa]+1;
               TwoSRU = sectorTwoS1[ikappa]-1;
               IRU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
               TwoSRD = sectorTwoS1[ikappa];
               IRD = IDL;
               break;
            case 1:
               NRU = sectorN1[ikappa]+1;
               TwoSRU = sectorTwoS1[ikappa]+1;
               IRU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
               TwoSRD = sectorTwoS1[ikappa];
               IRD = IDL;
               break;
            case 2:
               NRU = sectorN1[ikappa]+2;
               TwoSRU = sectorTwoS1[ikappa];
               IRU = sectorI1[ikappa];
               TwoSRD = sectorTwoS1[ikappa]-1;
               IRD = Irreps::directProd( sectorI1[ikappa] , denL->gIdiff() );
               break;
            case 3:
               NRU = sectorN1[ikappa]+2;
               TwoSRU = sectorTwoS1[ikappa];
               IRU = sectorI1[ikappa];
               TwoSRD = sectorTwoS1[ikappa]+1;
               IRD = Irreps::directProd( sectorI1[ikappa] , denL->gIdiff() );
               break;
         }
         int dimRU = denBK->gCurrentDim(index+1, NRU,   TwoSRU, IRU);
         int dimRD = denBK->gCurrentDim(index+1, NRU+1, TwoSRD, IRD);
         if ((dimRU>0) && (dimRD>0)){
         
            double * BlockTup   = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], NRU,   TwoSRU, IRU);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]+2, sectorTwoS1[ikappa], IDL,              NRU+1, TwoSRD, IRD);
            double * BlockL     = denL->gStorage(NRU,                TwoSRU,              IRU,              NRU+1, TwoSRD, IRD);
            
            //factor * Tup * L -> mem
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               int fase = ((((sectorTwoS1[ikappa] - TwoSRU + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(0.5 * (TwoSRU+1.0) / (sectorTwoS1[ikappa]+1.0) );
            } else {
               alpha = - sqrt(0.5) * (TwoSRD+1.0) / (sectorTwoS1[ikappa]+1.0);
            }
            double beta = 0.0; //set
            dgemm_(&notrans,&notrans,&dimUL,&dimRD,&dimRU,&alpha,BlockTup,&dimUL,BlockL,&dimRU,&beta,workmem,&dimUL);
            
            //mem * Tdown^T -> storage
            char trans = 'T';
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimRD,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
         }
      }
   } 

}

