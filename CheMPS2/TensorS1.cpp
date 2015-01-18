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

#include "TensorS1.h"
#include "Lapack.h"
#include "Gsl.h"

CheMPS2::TensorS1::TensorS1(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : TensorS1Bbase(indexIn, IdiffIn, movingRightIn, denBKIn){

}

CheMPS2::TensorS1::~TensorS1(){

}

void CheMPS2::TensorS1::makenew(TensorL * denL, TensorT * denT, double * workmem){

   if (movingRight){ makenewRight(denL, denT, workmem); }
   else{ makenewLeft( denL, denT, workmem); }

}

void CheMPS2::TensorS1::makenewRight(TensorL * denL, TensorT * denT, double * workmem){

   Clear();

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDR = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUR = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDR = denBK->gCurrentDim(index, sectorN1[ikappa]+2, sectorTwoSD[ikappa], IDR             );
      
      for (int geval=0; geval<4; geval++){
         int NLU,TwoSLU,ILU,TwoSLD,ILD; //NLD = NLU+1
         switch(geval){
            case 0:
               NLU = sectorN1[ikappa];
               TwoSLU = sectorTwoS1[ikappa];
               ILU = sectorI1[ikappa];
               TwoSLD = sectorTwoSD[ikappa]-1;
               ILD = Irreps::directProd( ILU, denL->gIdiff() );
               break;
            case 1:
               NLU = sectorN1[ikappa];
               TwoSLU = sectorTwoS1[ikappa];
               ILU = sectorI1[ikappa];
               TwoSLD = sectorTwoSD[ikappa]+1;
               ILD = Irreps::directProd( ILU, denL->gIdiff() );
               break;
            case 2:
               NLU = sectorN1[ikappa]-1;
               TwoSLU = sectorTwoS1[ikappa]-1;
               ILU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
               TwoSLD = sectorTwoSD[ikappa];
               ILD = IDR;
               break;
            case 3:
               NLU = sectorN1[ikappa]-1;
               TwoSLU = sectorTwoS1[ikappa]+1;
               ILU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
               TwoSLD = sectorTwoSD[ikappa];
               ILD = IDR;
               break;
         }
         int dimLU = denBK->gCurrentDim(index-1, NLU,   TwoSLU, ILU);
         int dimLD = denBK->gCurrentDim(index-1, NLU+1, TwoSLD, ILD);
         if ((dimLU>0) && (dimLD>0) && (abs(TwoSLU-TwoSLD)<2)){
         
            double * BlockTup   = denT->gStorage(NLU,   TwoSLU, ILU, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(NLU+1, TwoSLD, ILD, sectorN1[ikappa]+2, sectorTwoSD[ikappa], IDR);
            double * BlockL     = denL->gStorage(NLU,   TwoSLU, ILU, NLU+1,              TwoSLD,              ILD);
            
            //factor * Tup^T * L -> mem
            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               int fase = ((((sectorTwoS1[ikappa] + sectorTwoSD[ikappa] + 2)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0*(TwoSLD+1)) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSLD);
            } else {
               int fase = ((((TwoSLU + sectorTwoSD[ikappa] + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0*(sectorTwoS1[ikappa]+1)) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSLU);
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

void CheMPS2::TensorS1::makenewLeft(TensorL * denL, TensorT * denT, double * workmem){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDL = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimUL = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimDL = denBK->gCurrentDim(index, sectorN1[ikappa]+2, sectorTwoSD[ikappa], IDL             );
      
      for (int geval=0; geval<4; geval++){
         int NRU,TwoSRU,IRU,TwoSRD,IRD; //NRD = NRU+1
         switch(geval){
            case 0:
               NRU = sectorN1[ikappa]+1;
               TwoSRU = sectorTwoS1[ikappa]-1;
               IRU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
               TwoSRD = sectorTwoSD[ikappa];
               IRD = IDL;
               break;
            case 1:
               NRU = sectorN1[ikappa]+1;
               TwoSRU = sectorTwoS1[ikappa]+1;
               IRU = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
               TwoSRD = sectorTwoSD[ikappa];
               IRD = IDL;
               break;
            case 2:
               NRU = sectorN1[ikappa]+2;
               TwoSRU = sectorTwoS1[ikappa];
               IRU = sectorI1[ikappa];
               TwoSRD = sectorTwoSD[ikappa]-1;
               IRD = Irreps::directProd( sectorI1[ikappa] , denL->gIdiff() );
               break;
            case 3:
               NRU = sectorN1[ikappa]+2;
               TwoSRU = sectorTwoS1[ikappa];
               IRU = sectorI1[ikappa];
               TwoSRD = sectorTwoSD[ikappa]+1;
               IRD = Irreps::directProd( sectorI1[ikappa] , denL->gIdiff() );
               break;
         }
         int dimRU = denBK->gCurrentDim(index+1, NRU,   TwoSRU, IRU);
         int dimRD = denBK->gCurrentDim(index+1, NRU+1, TwoSRD, IRD);
         if ((dimRU>0) && (dimRD>0) && (abs(TwoSRD-TwoSRU)<2)){
         
            double * BlockTup   = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], NRU,   TwoSRU, IRU);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]+2, sectorTwoSD[ikappa], IDL,              NRU+1, TwoSRD, IRD);
            double * BlockL     = denL->gStorage(NRU,                TwoSRU,              IRU,              NRU+1, TwoSRD, IRD);
            
            //factor * Tup * L -> mem
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               int fase = ((((sectorTwoS1[ikappa] + sectorTwoSD[ikappa] + 2)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0 * (TwoSRU+1)) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSRU);
            } else {
               int fase = ((((sectorTwoS1[ikappa] + TwoSRD + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0 / (sectorTwoSD[ikappa] + 1.0)) * (TwoSRD + 1) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSRD);
            }
            double beta = 0.0; //set
            dgemm_(&notrans,&notrans,&dimUL,&dimRD,&dimRU,&alpha,BlockTup,&dimUL,BlockL,&dimRU,&beta,workmem,&dimUL);
            
            //mem * Tdown -> storage
            char trans = 'T';
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimRD,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
         }
      }
   }

}

