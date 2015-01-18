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

#include "TensorS0Abase.h"
#include "Lapack.h"

CheMPS2::TensorS0Abase::TensorS0Abase(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn) : Tensor(){

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
               int dimD = denBK->gCurrentDim(index,NU+2,TwoSU,ID);
               if (dimD>0){
                  nKappa++;
               }
            }
         }
      }
   }
   
   sectorN1 = new int[nKappa];
   sectorTwoS1 = new int[nKappa];
   sectorI1 = new int[nKappa];
   kappa2index = new int[nKappa+1];
   kappa2index[0] = 0;
   
   nKappa = 0;
   for (int NU=denBK->gNmin(index); NU<=denBK->gNmax(index); NU++){
      for (int TwoSU=denBK->gTwoSmin(index,NU); TwoSU<=denBK->gTwoSmax(index,NU); TwoSU+=2){
         for (int IU=0; IU<denBK->getNumberOfIrreps(); IU++){
            int dimU = denBK->gCurrentDim(index,NU,TwoSU,IU);
            if (dimU>0){
               int ID = Irreps::directProd(Idiff,IU);
               int dimD = denBK->gCurrentDim(index,NU+2,TwoSU,ID);
               if (dimD>0){
                  sectorN1[nKappa] = NU;
                  sectorTwoS1[nKappa] = TwoSU;
                  sectorI1[nKappa] = IU;
                  nKappa++;
                  kappa2index[nKappa] = kappa2index[nKappa-1] + dimU*dimD;
               }
            }
         }
      }
   }
   
   storage = new double[kappa2index[nKappa]];

}

CheMPS2::TensorS0Abase::~TensorS0Abase(){

   delete [] sectorN1;
   delete [] sectorTwoS1;
   delete [] sectorI1;
   delete [] kappa2index;
   delete [] storage;

}

int CheMPS2::TensorS0Abase::gNKappa() const { return nKappa; }

double * CheMPS2::TensorS0Abase::gStorage() { return storage; }
      
int CheMPS2::TensorS0Abase::gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const{

   if ((Irreps::directProd(I1,Idiff))!=I2) return -1;
   if (N2!=N1+2) return -1;
   if (TwoS1!=TwoS2) return -1;

   for (int cnt=0; cnt<nKappa; cnt++){
      if ((sectorN1[cnt]==N1)&&(sectorTwoS1[cnt]==TwoS1)&&(sectorI1[cnt]==I1)) return cnt;
   }
   
   return -1;

}
      
int CheMPS2::TensorS0Abase::gKappa2index(const int kappa) const{ return kappa2index[kappa]; }
      
double * CheMPS2::TensorS0Abase::gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2){

   int kappa = gKappa(N1,TwoS1,I1,N2,TwoS2,I2);
   if (kappa == -1) return NULL;
   return storage + kappa2index[kappa];

}

int CheMPS2::TensorS0Abase::gIndex() const { return index; }

int CheMPS2::TensorS0Abase::gIdiff() const { return Idiff; }

void CheMPS2::TensorS0Abase::Clear(){

   for (int cnt=0; cnt<kappa2index[nKappa]; cnt++){ storage[cnt] = 0.0; }

}

void CheMPS2::TensorS0Abase::update(TensorS0Abase * S0AbasePrevious, TensorT * denT, double * workmem){

   if (movingRight){ updateMovingRight(S0AbasePrevious, denT, workmem); }
   else{ updateMovingLeft( S0AbasePrevious, denT, workmem); }

}

void CheMPS2::TensorS0Abase::updateMovingRight(TensorS0Abase * S0AbasePrevious, TensorT * denT, double * workmem){

   Clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUR = denBK->gCurrentDim(index,sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa]);
      int IDR = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimDR = denBK->gCurrentDim(index,sectorN1[ikappa]+2,sectorTwoS1[ikappa],IDR);
      
      for (int geval=0; geval<4; geval++){
         int NUL,TwoSL,IUL,IDL; //NDL = NUL+2 ; TwoSDL = TwoSUL
         switch(geval){
            case 0: //N = 0
               NUL = sectorN1[ikappa];
               TwoSL = sectorTwoS1[ikappa];
               IUL = sectorI1[ikappa];
               IDL = IDR;
               break;
            case 1: //N = 2
               NUL = sectorN1[ikappa] - 2;
               TwoSL = sectorTwoS1[ikappa];
               IUL = sectorI1[ikappa];
               IDL = IDR;
               break;
            case 2: //N = 1
               NUL = sectorN1[ikappa] - 1;
               TwoSL = sectorTwoS1[ikappa] - 1;
               IUL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1));
               IDL = Irreps::directProd( IDR , denBK->gIrrep(index-1) );
               break;
            case 3: //N = 1
               NUL = sectorN1[ikappa] - 1;
               TwoSL = sectorTwoS1[ikappa] + 1;
               IUL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1));
               IDL = Irreps::directProd( IDR , denBK->gIrrep(index-1) );
               break;
         }
         int dimUL = denBK->gCurrentDim(index-1,NUL,  TwoSL,IUL);
         int dimDL = denBK->gCurrentDim(index-1,NUL+2,TwoSL,IDL);
         
         if ((dimUL>0) && (dimDL>0)){
            
            double * BlockTup   = denT->gStorage(NUL,  TwoSL,IUL,sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa]);
            double * BlockTdown = denT->gStorage(NUL+2,TwoSL,IDL,sectorN1[ikappa]+2,sectorTwoS1[ikappa],IDR             );
            double * BlockS0AbaseLeft = S0AbasePrevious->gStorage(NUL,TwoSL,IUL,NUL+2,TwoSL,IDL);
         
            //Tup^T * S0Abase --> mem
            double alpha = 1.0; //factor = 1 for Sigma 0 tensor
            char trans = 'T';
            char notr = 'N';
            double beta = 0.0; //set mem
            dgemm_(&trans,&notr,&dimUR,&dimDL,&dimUL,&alpha,BlockTup,&dimUL,BlockS0AbaseLeft,&dimUL,&beta,workmem,&dimUR);
            
            //mem * Tdown --> storage + kappa2index[ikappa]
            beta = 1.0; //add to storage
            dgemm_(&notr,&notr,&dimUR,&dimDR,&dimDL,&alpha,workmem,&dimUR,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUR);

         }

      }
   }

}

void CheMPS2::TensorS0Abase::updateMovingLeft(TensorS0Abase * S0AbasePrevious, TensorT * denT, double * workmem){

   Clear();

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimUL = denBK->gCurrentDim(index,sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa]);
      int IDL = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimDL = denBK->gCurrentDim(index,sectorN1[ikappa]+2,sectorTwoS1[ikappa],IDL);
      
      for (int geval=0; geval<4; geval++){
         int NUR,TwoSR,IUR,IDR; //NDR = NUR + 2; TwoSDR = TwoSUR
         switch(geval){
            case 0: //N = 0
               NUR = sectorN1[ikappa];
               TwoSR = sectorTwoS1[ikappa];
               IUR = sectorI1[ikappa];
               IDR = IDL;
               break;
            case 1: //N = 2
               NUR = sectorN1[ikappa] + 2;
               TwoSR = sectorTwoS1[ikappa];
               IUR = sectorI1[ikappa];
               IDR = IDL;
               break;
            case 2: //N = 1
               NUR = sectorN1[ikappa] + 1;
               TwoSR = sectorTwoS1[ikappa] - 1;
               IUR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index));
               IDR = Irreps::directProd( IDL , denBK->gIrrep(index) );
               break;
            case 3: //N = 1
               NUR = sectorN1[ikappa] + 1;
               TwoSR = sectorTwoS1[ikappa] + 1;
               IUR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index));
               IDR = Irreps::directProd( IDL , denBK->gIrrep(index) );
               break;
         }
         int dimUR = denBK->gCurrentDim(index+1,NUR,  TwoSR,IUR);
         int dimDR = denBK->gCurrentDim(index+1,NUR+2,TwoSR,IDR);
         
         if ((dimUR>0) && (dimDR>0)){
            
            double * BlockTup   = denT->gStorage(sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa],NUR,  TwoSR,IUR);
            double * BlockTdown = denT->gStorage(sectorN1[ikappa]+2,sectorTwoS1[ikappa],IDL,             NUR+2,TwoSR,IDR);
            double * BlockS0AbaseRight = S0AbasePrevious->gStorage(NUR,TwoSR,IUR,NUR+2,TwoSR,IDR);
         
            //factor * Tup * Swap --> mem
            double alpha = 1.0;
            if (geval>=2){
               alpha = (TwoSR+1.0)/(sectorTwoS1[ikappa]+1.0);
            }
            char notr = 'N';
            double beta = 0.0; //set mem
            dgemm_(&notr,&notr,&dimUL,&dimDR,&dimUR,&alpha,BlockTup,&dimUL,BlockS0AbaseRight,&dimUR,&beta,workmem,&dimUL);
            
            //mem * Tdown^T --> storage + kappa2index[ikappa]
            beta = 1.0; //add to storage
            alpha = 1.0; //factor already with previous multiplication
            char trans = 'T';
            dgemm_(&notr,&trans,&dimUL,&dimDL,&dimDR,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);

         }

      }
   }

}



