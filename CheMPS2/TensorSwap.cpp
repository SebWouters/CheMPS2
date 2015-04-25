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

#include "TensorSwap.h"
#include "Lapack.h"
#include "Gsl.h"

CheMPS2::TensorSwap::TensorSwap(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : Tensor(){

   index = indexIn; //boundary = index
   Idiff = IdiffIn;
   movingRight = movingRightIn;
   denBK = denBKIn;
   
   nKappa = 0;
   for (int NU=denBK->gNmin(index); NU<=denBK->gNmax(index); NU++){
      for (int TwoSU=denBK->gTwoSmin(index,NU); TwoSU<=denBK->gTwoSmax(index,NU); TwoSU+=2){
         for (int IU=0; IU<denBK->getNumberOfIrreps(); IU++){
            int dimU = denBK->gCurrentDim(index,NU,TwoSU,IU);
            if (dimU>0){
               int ID = Irreps::directProd(Idiff,IU);
               for (int TwoSD = TwoSU-1; TwoSD<=TwoSU+1; TwoSD+=2){
                  if (TwoSD>=0){
                     int dimD = denBK->gCurrentDim(index,NU+1,TwoSD,ID);
                     if (dimD>0){
                        nKappa++;
                     }
                  }
               }
            }
         }
      }
   }
   
   sectorN1 = new int[nKappa];
   sectorTwoS1 = new int[nKappa];
   sectorI1 = new int[nKappa];
   sectorTwoSD = new int[nKappa];
   kappa2index = new int[nKappa+1];
   kappa2index[0] = 0;
   
   nKappa = 0;
   for (int NU=denBK->gNmin(index); NU<=denBK->gNmax(index); NU++){
      for (int TwoSU=denBK->gTwoSmin(index,NU); TwoSU<=denBK->gTwoSmax(index,NU); TwoSU+=2){
         for (int IU=0; IU<denBK->getNumberOfIrreps(); IU++){
            int dimU = denBK->gCurrentDim(index,NU,TwoSU,IU);
            if (dimU>0){
               int ID = Irreps::directProd(Idiff,IU);
               for (int TwoSD = TwoSU-1; TwoSD<=TwoSU+1; TwoSD+=2){
                  if (TwoSD>=0){
                     int dimD = denBK->gCurrentDim(index,NU+1,TwoSD,ID);
                     if (dimD>0){
                        sectorN1[nKappa] = NU;
                        sectorTwoS1[nKappa] = TwoSU;
                        sectorI1[nKappa] = IU;
                        sectorTwoSD[nKappa] = TwoSD;
                        nKappa++;
                        kappa2index[nKappa] = kappa2index[nKappa-1] + dimU*dimD;
                     }
                  }
               }
            }
         }
      }
   }
   
   storage = new double[kappa2index[nKappa]];

}

CheMPS2::TensorSwap::~TensorSwap(){

   delete [] sectorN1;
   delete [] sectorTwoS1;
   delete [] sectorI1;
   delete [] sectorTwoSD;
   delete [] kappa2index;
   delete [] storage;

}

int CheMPS2::TensorSwap::gNKappa() const { return nKappa; }

double * CheMPS2::TensorSwap::gStorage() { return storage; }
      
int CheMPS2::TensorSwap::gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const{

   if ((Irreps::directProd(I1,Idiff))!=I2) return -1;
   if (N2!=N1+1) return -1;

   for (int cnt=0; cnt<nKappa; cnt++){
      if ((sectorN1[cnt]==N1)&&(sectorTwoS1[cnt]==TwoS1)&&(sectorI1[cnt]==I1)&&(sectorTwoSD[cnt]==TwoS2)) return cnt;
   }
   
   return -1;

}
      
int CheMPS2::TensorSwap::gKappa2index(const int kappa) const{ return kappa2index[kappa]; }
      
double * CheMPS2::TensorSwap::gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2){

   int kappa = gKappa(N1,TwoS1,I1,N2,TwoS2,I2);
   if (kappa == -1) return NULL;
   return storage + kappa2index[kappa];

}

int CheMPS2::TensorSwap::gIndex() const { return index; }

int CheMPS2::TensorSwap::gIdiff() const { return Idiff; }

void CheMPS2::TensorSwap::Clear(){

   for (int cnt=0; cnt<kappa2index[nKappa]; cnt++) storage[cnt] = 0.0;

}

void CheMPS2::TensorSwap::update(TensorSwap * SwapPrevious, TensorT * denT, double * workmem, const bool JordanWigner){

   if (movingRight){ updateMovingRight(SwapPrevious, denT, workmem, JordanWigner); }
   else{ updateMovingLeft( SwapPrevious, denT, workmem, JordanWigner); }

}

void CheMPS2::TensorSwap::updateMovingRight(TensorSwap * SwapPrevious, TensorT * denT, double * workmem, const bool JordanWigner){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUR = denBK->gCurrentDim(index,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
      int IDR = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimDR = denBK->gCurrentDim(index,sectorN1[ikappa]+1,sectorTwoSD[ikappa],IDR);
      
      for (int geval=0; geval<5; geval++){
         int NUL,TwoSUL,IUL,TwoSDL,IDL;
         switch(geval){
            case 0: //N = 0
               NUL = sectorN1[ikappa];
               TwoSUL = sectorTwoS1[ikappa];
               IUL = sectorI1[ikappa];
               TwoSDL = sectorTwoSD[ikappa];
               IDL = IDR;
               break;
            case 1: //N = 2
               NUL = sectorN1[ikappa] - 2;
               TwoSUL = sectorTwoS1[ikappa];
               IUL = sectorI1[ikappa];
               TwoSDL = sectorTwoSD[ikappa];
               IDL = IDR;
               break;
            case 2: //N = 1, both TwoJLs minus one
               NUL = sectorN1[ikappa] - 1;
               TwoSUL = sectorTwoS1[ikappa] - 1;
               IUL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1));
               TwoSDL = sectorTwoSD[ikappa] - 1;
               IDL = Irreps::directProd( IDR , denBK->gIrrep(index-1) );
               break;
            case 3: //N = 1, both TwoJLs plus one
               NUL = sectorN1[ikappa] - 1;
               TwoSUL = sectorTwoS1[ikappa] + 1;
               IUL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1));
               TwoSDL = sectorTwoSD[ikappa] + 1;
               IDL = Irreps::directProd( IDR , denBK->gIrrep(index-1) );
               break;
            case 4: //N = 1, TwoJLs are swap of TwoJRs
               NUL = sectorN1[ikappa] - 1;
               TwoSUL = sectorTwoSD[ikappa];
               IUL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1));
               TwoSDL = sectorTwoS1[ikappa];
               IDL = Irreps::directProd( IDR , denBK->gIrrep(index-1) );
               break;
         }
         int dimUL = denBK->gCurrentDim(index-1,NUL,  TwoSUL,IUL);
         int dimDL = denBK->gCurrentDim(index-1,NUL+1,TwoSDL,IDL);
         
         if ((dimUL>0) && (dimDL>0)){
            
            double * BlockTup   = denT->gStorage(NUL,  TwoSUL,IUL,sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(NUL+1,TwoSDL,IDL,sectorN1[ikappa]+1,sectorTwoSD[ikappa],IDR             );
            double * BlockSwapLeft = SwapPrevious->gStorage(NUL,TwoSUL,IUL,NUL+1,TwoSDL,IDL);
         
            //factor * Tup^T * Swap --> mem
            double alpha = 1.0;
            if (geval>=2){
               int fase = ((((TwoSDL+sectorTwoS1[ikappa])/2)%2)!=0)?-1:1;
               if (!JordanWigner){ fase *= -1; }
               alpha = fase * sqrt((TwoSDL+1.0)*(sectorTwoS1[ikappa]+1.0)) * gsl_sf_coupling_6j(sectorTwoS1[ikappa],sectorTwoSD[ikappa],1,TwoSDL,TwoSUL,1);
            }
            char trans = 'T';
            char notr = 'N';
            double beta = 0.0; //set mem
            dgemm_(&trans,&notr,&dimUR,&dimDL,&dimUL,&alpha,BlockTup,&dimUL,BlockSwapLeft,&dimUL,&beta,workmem,&dimUR);
            
            //mem * Tdown --> storage + kappa2index[ikappa]
            beta = 1.0; //add to storage
            alpha = 1.0; //factor already with previous multiplication
            dgemm_(&notr,&notr,&dimUR,&dimDR,&dimDL,&alpha,workmem,&dimUR,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUR);

         }

      }
   }

}

void CheMPS2::TensorSwap::updateMovingLeft(TensorSwap * SwapPrevious, TensorT * denT, double * workmem, const bool JordanWigner){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUL = denBK->gCurrentDim(index,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
      int IDL = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimDL = denBK->gCurrentDim(index,sectorN1[ikappa]+1,sectorTwoSD[ikappa],IDL);
      
      for (int geval=0; geval<5; geval++){
         int NUR,TwoSUR,IUR,TwoSDR,IDR;
         switch(geval){
            case 0: //N = 0
               NUR = sectorN1[ikappa];
               TwoSUR = sectorTwoS1[ikappa];
               IUR = sectorI1[ikappa];
               TwoSDR = sectorTwoSD[ikappa];
               IDR = IDL;
               break;
            case 1: //N = 2
               NUR = sectorN1[ikappa] + 2;
               TwoSUR = sectorTwoS1[ikappa];
               IUR = sectorI1[ikappa];
               TwoSDR = sectorTwoSD[ikappa];
               IDR = IDL;
               break;
            case 2: //N = 1, both TwoJRs minus one
               NUR = sectorN1[ikappa] + 1;
               TwoSUR = sectorTwoS1[ikappa] - 1;
               IUR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index));
               TwoSDR = sectorTwoSD[ikappa] - 1;
               IDR = Irreps::directProd( IDL , denBK->gIrrep(index) );
               break;
            case 3: //N = 1, both TwoJRs plus one
               NUR = sectorN1[ikappa] + 1;
               TwoSUR = sectorTwoS1[ikappa] + 1;
               IUR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index));
               TwoSDR = sectorTwoSD[ikappa] + 1;
               IDR = Irreps::directProd( IDL , denBK->gIrrep(index) );
               break;
            case 4: //N = 1, TwoJRs are swap of TwoJLs
               NUR = sectorN1[ikappa] + 1;
               TwoSUR = sectorTwoSD[ikappa];
               IUR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index));
               TwoSDR = sectorTwoS1[ikappa];
               IDR = Irreps::directProd( IDL , denBK->gIrrep(index) );
               break;
         }
         int dimUR = denBK->gCurrentDim(index+1,NUR,  TwoSUR,IUR);
         int dimDR = denBK->gCurrentDim(index+1,NUR+1,TwoSDR,IDR);
         
         if ((dimUR>0) && (dimDR>0)){
            
            double * BlockTup   = denT->gStorage(sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa],NUR,  TwoSUR,IUR);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]+1,sectorTwoSD[ikappa],IDL,             NUR+1,TwoSDR,IDR);
            double * BlockSwapRight = SwapPrevious->gStorage(NUR,TwoSUR,IUR,NUR+1,TwoSDR,IDR);
         
            //factor * Tup * Swap --> mem
            double alpha = 1.0;
            if (geval>=2){
               int fase = ((((sectorTwoSD[ikappa]+TwoSUR)/2)%2)!=0)?-1:1;
               if (!JordanWigner){ fase *= -1; }
               alpha = fase*(TwoSDR+1)*sqrt((TwoSUR+1.0)/(sectorTwoSD[ikappa]+1.0))*gsl_sf_coupling_6j(TwoSUR,TwoSDR,1,sectorTwoSD[ikappa],sectorTwoS1[ikappa],1);
            }
            char notr = 'N';
            double beta = 0.0; //set mem
            dgemm_(&notr,&notr,&dimUL,&dimDR,&dimUR,&alpha,BlockTup,&dimUL,BlockSwapRight,&dimUR,&beta,workmem,&dimUL);
            
            //mem * Tdown^T --> storage + kappa2index[ikappa]
            beta = 1.0; //add to storage
            alpha = 1.0; //factor already with previous multiplication
            char trans = 'T';
            dgemm_(&notr,&trans,&dimUL,&dimDL,&dimDR,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);

         }

      }
   }

}



