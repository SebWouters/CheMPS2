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
#include <algorithm>

#include "TensorO.h"
#include "Lapack.h"

CheMPS2::TensorO::TensorO(const int indexIn, const bool movingRightIn, const SyBookkeeper * denBKupIn, const SyBookkeeper * denBKdownIn, const Problem * ProbIn) : Tensor(){

   index = indexIn;
   denBKup = denBKupIn;
   denBK = denBKdownIn;
   Prob = ProbIn;
   movingRight = movingRightIn;
   
   nKappa = 0;
   for (int N=denBK->gNmin(index); N<=denBK->gNmax(index); N++){
      for (int TwoS=denBK->gTwoSmin(index,N); TwoS<=denBK->gTwoSmax(index,N); TwoS+=2){
         for (int Icnt=0; Icnt<denBK->getNumberOfIrreps(); Icnt++){
            int dimDown = denBK->gCurrentDim(index,N,TwoS,Icnt);
            int dimUp = denBKup->gCurrentDim(index,N,TwoS,Icnt);
            if ((dimDown>0) && (dimUp>0)) nKappa++;
         }
      }
   }
   
   sectorN1 = new int[nKappa];
   sectorTwoS1 = new int[nKappa];
   sectorI1 = new int[nKappa];
   kappa2index = new int[nKappa+1];
   kappa2index[0] = 0;
   
   nKappa = 0;
   for (int N=denBK->gNmin(index); N<=denBK->gNmax(index); N++){
      for (int TwoS=denBK->gTwoSmin(index,N); TwoS<=denBK->gTwoSmax(index,N); TwoS+=2){
         for (int Icnt=0; Icnt<denBK->getNumberOfIrreps(); Icnt++){
            int dimDown = denBK->gCurrentDim(index,N,TwoS,Icnt);
            int dimUp = denBKup->gCurrentDim(index,N,TwoS,Icnt);
            if ((dimDown>0) && (dimUp>0)){
               sectorN1[nKappa] = N;
               sectorTwoS1[nKappa] = TwoS;
               sectorI1[nKappa] = Icnt;
               nKappa++;
               kappa2index[nKappa] = kappa2index[nKappa-1] + dimDown*dimUp;
            }
         }
      }
   }
   
   storage = new double[kappa2index[nKappa]];

}

CheMPS2::TensorO::~TensorO(){

   delete [] sectorN1;
   delete [] sectorTwoS1;
   delete [] sectorI1;
   delete [] kappa2index;
   delete [] storage;

}

int CheMPS2::TensorO::gNKappa() const { return nKappa; }

double * CheMPS2::TensorO::gStorage() { return storage; }
      
int CheMPS2::TensorO::gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const{

   if ((N1!=N2) || (TwoS1!=TwoS2) || (I1!=I2)) return -1;

   for (int cnt=0; cnt<nKappa; cnt++){
      if ( (sectorN1[cnt]==N1) && (sectorTwoS1[cnt]==TwoS1) && (sectorI1[cnt]==I1) ) return cnt;
   }
   
   return -1;

}
      
int CheMPS2::TensorO::gKappa2index(const int kappa) const{ return kappa2index[kappa]; }
      
double * CheMPS2::TensorO::gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2){

   int kappa = gKappa(N1,TwoS1,I1,N2,TwoS2,I2);
   if (kappa == -1) return NULL;
   return storage + kappa2index[kappa];

}

int CheMPS2::TensorO::gIndex() const { return index; }

void CheMPS2::TensorO::Clear(){ for (int cnt=0; cnt<kappa2index[nKappa]; cnt++) storage[cnt] = 0.0; }

void CheMPS2::TensorO::update(TensorT * denTup, TensorT * denTdown){

   Clear();

   if (movingRight){
      //PARALLEL
      #pragma omp parallel for schedule(dynamic)
      for (int ikappa=0; ikappa<nKappa; ikappa++){ updateRight(ikappa, denTup, denTdown); }
   } else {
      //PARALLEL
      #pragma omp parallel for schedule(dynamic)
      for (int ikappa=0; ikappa<nKappa; ikappa++){ updateLeft(ikappa, denTup, denTdown); }
   }

}

void CheMPS2::TensorO::update(TensorT * denTup, TensorT * denTdown, TensorO * denO){

   Clear();

   if (movingRight){
   
      if (index>1){
   
         const int dimL = std::max(denBK->gMaxDimAtBound(index-1),denBKup->gMaxDimAtBound(index-1));
         const int dimR = std::max(denBK->gMaxDimAtBound(index),  denBKup->gMaxDimAtBound(index)  );
      
         //PARALLEL
         #pragma omp parallel
         {
         
            double * workmem = new double[dimL*dimR];
         
            #pragma omp for schedule(dynamic)
            for (int ikappa=0; ikappa<nKappa; ikappa++){ updateRight(ikappa, denTup, denTdown, denO, workmem); }
            
            delete [] workmem;
         
         }
      }
   } else {
   
      if (index<Prob->gL()-1){
   
         const int dimL = std::max(denBK->gMaxDimAtBound(index),  denBKup->gMaxDimAtBound(index)  );
         const int dimR = std::max(denBK->gMaxDimAtBound(index+1),denBKup->gMaxDimAtBound(index+1));
         
         //PARALLEL
         #pragma omp parallel
         {
         
            double * workmem = new double[dimL*dimR];
         
            #pragma omp for schedule(dynamic)
            for (int ikappa=0; ikappa<nKappa; ikappa++){ updateLeft(ikappa, denTup, denTdown, denO, workmem); }
            
            delete [] workmem;
      
         }
      }
   }

}

void CheMPS2::TensorO::updateRight(const int ikappa, TensorT * denTup, TensorT * denTdown){
   
   int dimRdown = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   int dimRup = denBKup->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   
   for (int geval=0; geval<4; geval++){
      int IL, TwoSL, NL;
      switch(geval){
         case 0:
            NL = sectorN1[ikappa];
            TwoSL = sectorTwoS1[ikappa];
            IL = sectorI1[ikappa];
            break;
         case 1:
            NL = sectorN1[ikappa]-2;
            TwoSL = sectorTwoS1[ikappa];
            IL = sectorI1[ikappa];
            break;
         case 2:
            NL = sectorN1[ikappa]-1;
            TwoSL = sectorTwoS1[ikappa]-1;
            IL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
            break;
         case 3:
            NL = sectorN1[ikappa]-1;
            TwoSL = sectorTwoS1[ikappa]+1;
            IL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
            break;
      }
      
      int dimLup = denBKup->gCurrentDim(index-1, NL, TwoSL, IL);
      int dimLdown = denBK->gCurrentDim(index-1, NL, TwoSL, IL);
      if ((dimLup>0) && (dimLdown>0) && (dimLup==dimLdown)){
      
         double alpha = 1.0;
         double beta = 1.0; //add
         char trans = 'T';
         char notrans = 'N';
         double * Tup   =   denTup->gStorage(NL,TwoSL,IL,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
         double * Tdown = denTdown->gStorage(NL,TwoSL,IL,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
         dgemm_(&trans,&notrans,&dimRdown,&dimRup,&dimLup,&alpha,Tdown,&dimLdown,Tup,&dimLup,&beta,storage+gKappa2index(ikappa),&dimRdown);
      
      }
   }

}

void CheMPS2::TensorO::updateLeft(const int ikappa, TensorT * denTup, TensorT * denTdown){
   
   int dimLdown = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   int dimLup = denBKup->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   
   for (int geval=0; geval<4; geval++){
      int IR, TwoSR, NR;
      switch(geval){
         case 0:
            NR = sectorN1[ikappa];
            TwoSR = sectorTwoS1[ikappa];
            IR = sectorI1[ikappa];
            break;
         case 1:
            NR = sectorN1[ikappa]+2;
            TwoSR = sectorTwoS1[ikappa];
            IR = sectorI1[ikappa];
            break;
         case 2:
            NR = sectorN1[ikappa]+1;
            TwoSR = sectorTwoS1[ikappa]-1;
            IR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
            break;
         case 3:
            NR = sectorN1[ikappa]+1;
            TwoSR = sectorTwoS1[ikappa]+1;
            IR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
            break;
      }
      
      int dimRup = denBKup->gCurrentDim(index+1, NR, TwoSR, IR);
      int dimRdown = denBK->gCurrentDim(index+1, NR, TwoSR, IR);
      if ((dimRup>0) && (dimRdown>0) && (dimRup==dimRdown)){
      
         double alpha = (geval>1) ? (TwoSR+1.0)/(sectorTwoS1[ikappa]+1.0) : 1.0;
         double beta = 1.0; //add
         char trans = 'T';
         char notrans = 'N';
         double * Tup   =   denTup->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa],NR,TwoSR,IR);
         double * Tdown = denTdown->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa],NR,TwoSR,IR);
         dgemm_(&notrans,&trans,&dimLdown,&dimLup,&dimRup,&alpha,Tdown,&dimLdown,Tup,&dimLup,&beta,storage+gKappa2index(ikappa),&dimLdown);
      
      }
   }

}

void CheMPS2::TensorO::updateRight(const int ikappa, TensorT * denTup, TensorT * denTdown, TensorO * denO, double * workmem){
   
   int dimRdown = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   int dimRup = denBKup->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   
   for (int geval=0; geval<4; geval++){
      int IL, TwoSL, NL;
      switch(geval){
         case 0:
            NL = sectorN1[ikappa];
            TwoSL = sectorTwoS1[ikappa];
            IL = sectorI1[ikappa];
            break;
         case 1:
            NL = sectorN1[ikappa]-2;
            TwoSL = sectorTwoS1[ikappa];
            IL = sectorI1[ikappa];
            break;
         case 2:
            NL = sectorN1[ikappa]-1;
            TwoSL = sectorTwoS1[ikappa]-1;
            IL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
            break;
         case 3:
            NL = sectorN1[ikappa]-1;
            TwoSL = sectorTwoS1[ikappa]+1;
            IL = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index-1) );
            break;
      }
      
      int dimLup = denBKup->gCurrentDim(index-1, NL, TwoSL, IL);
      int dimLdown = denBK->gCurrentDim(index-1, NL, TwoSL, IL);
      if ((dimLup>0) && (dimLdown>0)){
      
         double alpha = 1.0;
         double beta = 0.0; //set
         char trans = 'T';
         char notrans = 'N';
         double * Tup   =   denTup->gStorage(NL,TwoSL,IL,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
         double * Tdown = denTdown->gStorage(NL,TwoSL,IL,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
         double * Opart = denO->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
         
         dgemm_(&trans,&notrans,&dimRdown,&dimLup,&dimLdown,&alpha,Tdown,&dimLdown,Opart,&dimLdown,&beta,workmem,&dimRdown);
         
         beta = 1.0; //add
         
         dgemm_(&notrans,&notrans,&dimRdown,&dimRup,&dimLup,&alpha,workmem,&dimRdown,Tup,&dimLup,&beta,storage+gKappa2index(ikappa),&dimRdown);
      
      }
   }

}

void CheMPS2::TensorO::updateLeft(const int ikappa, TensorT * denTup, TensorT * denTdown, TensorO * denO, double * workmem){
   
   int dimLdown = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   int dimLup = denBKup->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   
   for (int geval=0; geval<4; geval++){
      int IR, TwoSR, NR;
      switch(geval){
         case 0:
            NR = sectorN1[ikappa];
            TwoSR = sectorTwoS1[ikappa];
            IR = sectorI1[ikappa];
            break;
         case 1:
            NR = sectorN1[ikappa]+2;
            TwoSR = sectorTwoS1[ikappa];
            IR = sectorI1[ikappa];
            break;
         case 2:
            NR = sectorN1[ikappa]+1;
            TwoSR = sectorTwoS1[ikappa]-1;
            IR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
            break;
         case 3:
            NR = sectorN1[ikappa]+1;
            TwoSR = sectorTwoS1[ikappa]+1;
            IR = Irreps::directProd( sectorI1[ikappa] , denBK->gIrrep(index) );
            break;
      }
      
      int dimRup = denBKup->gCurrentDim(index+1, NR, TwoSR, IR);
      int dimRdown = denBK->gCurrentDim(index+1, NR, TwoSR, IR);
      if ((dimRup>0) && (dimRdown>0)){
      
         double alpha = (geval>1) ? (TwoSR+1.0)/(sectorTwoS1[ikappa]+1.0) : 1.0;
         double beta = 0.0; //set
         char trans = 'T';
         char notrans = 'N';
         double * Tup   =   denTup->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa],NR,TwoSR,IR);
         double * Tdown = denTdown->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa],NR,TwoSR,IR);
         double * Opart = denO->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
         
         dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,Tdown,&dimLdown,Opart,&dimRdown,&beta,workmem,&dimLdown);
         
         alpha = 1.0;
         beta = 1.0; //add
         
         dgemm_(&notrans,&trans,&dimLdown,&dimLup,&dimRup,&alpha,workmem,&dimLdown,Tup,&dimLup,&beta,storage+gKappa2index(ikappa),&dimLdown);
      
      }
   }

}



