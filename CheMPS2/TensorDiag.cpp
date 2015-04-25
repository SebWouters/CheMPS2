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

#include "TensorDiag.h"
#include "Lapack.h"

CheMPS2::TensorDiag::TensorDiag(const int indexIn, const SyBookkeeper * denBKIn) : Tensor(){

   index = indexIn;
   denBK = denBKIn;
   
   nKappa = 0;
   for (int N=denBK->gNmin(index); N<=denBK->gNmax(index); N++){
      for (int TwoS=denBK->gTwoSmin(index,N); TwoS<=denBK->gTwoSmax(index,N); TwoS+=2){
         for (int Icnt=0; Icnt<denBK->getNumberOfIrreps(); Icnt++){
            int dim = denBK->gCurrentDim(index,N,TwoS,Icnt);
            if (dim>0) nKappa++;
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
            int dim = denBK->gCurrentDim(index,N,TwoS,Icnt);
            if (dim>0){
               sectorN1[nKappa] = N;
               sectorTwoS1[nKappa] = TwoS;
               sectorI1[nKappa] = Icnt;
               nKappa++;
               kappa2index[nKappa] = kappa2index[nKappa-1] + dim*dim;
            }
         }
      }
   }
   
   storage = new double[kappa2index[nKappa]];

}

CheMPS2::TensorDiag::~TensorDiag(){

   delete [] sectorN1;
   delete [] sectorTwoS1;
   delete [] sectorI1;
   delete [] kappa2index;
   delete [] storage;

}

int CheMPS2::TensorDiag::gNKappa() const { return nKappa; }

double * CheMPS2::TensorDiag::gStorage() { return storage; }
      
int CheMPS2::TensorDiag::gKappa(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2) const{

   if ((N1!=N2) || (TwoS1!=TwoS2) || (I1!=I2)) return -1;

   for (int cnt=0; cnt<nKappa; cnt++){
      if ( (sectorN1[cnt]==N1) && (sectorTwoS1[cnt]==TwoS1) && (sectorI1[cnt]==I1) ) return cnt;
   }
   
   return -1;

}
      
int CheMPS2::TensorDiag::gKappa2index(const int kappa) const{ return kappa2index[kappa]; }
      
double * CheMPS2::TensorDiag::gStorage(const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2){

   int kappa = gKappa(N1,TwoS1,I1,N2,TwoS2,I2);
   if (kappa == -1) return NULL;
   return storage + kappa2index[kappa];

}

int CheMPS2::TensorDiag::gIndex() const { return index; }

void CheMPS2::TensorDiag::updateRight(const int ikappa, Tensor * denT, TensorDiag * diagPrevious, double * workmemLR){

   int dimR = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   for (int geval=0; geval<4; geval++){
      int NL,TwoSL,IL;
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
      int dimL = denBK->gCurrentDim(index-1,NL,TwoSL,IL);
      if (dimL>0){

         double * BlockT    = denT->gStorage(NL,TwoSL,IL,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
         double * BlockDiag = diagPrevious->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);

         //T^T * diag --> mem
         char trans = 'T';
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0; //set
         dgemm_(&trans,&notr,&dimR,&dimL,&dimL,&alpha,BlockT,&dimL,BlockDiag,&dimL,&beta,workmemLR,&dimR);

         //mem * T --> storage
         beta = 1.0; //add
         dgemm_(&notr,&notr,&dimR,&dimR,&dimL,&alpha,workmemLR,&dimR,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);

      }
   }

}

void CheMPS2::TensorDiag::updateLeft(const int ikappa, Tensor * denT, TensorDiag * diagPrevious, double * workmemLR){

   int dimL = denBK->gCurrentDim(index, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
   for (int geval=0; geval<4; geval++){
      int NR,TwoSR,IR;
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
      int dimR = denBK->gCurrentDim(index+1,NR,TwoSR,IR);
      if (dimR>0){

         double * BlockT    = denT->gStorage(sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa],NR,TwoSR,IR);
         double * BlockDiag = diagPrevious->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);

         //factor * T * diag --> mem
         char notr = 'N';
         double alpha = (geval>1) ? (TwoSR+1.0)/(sectorTwoS1[ikappa]+1.0) : 1.0;
         double beta = 0.0; //set
         dgemm_(&notr,&notr,&dimL,&dimR,&dimR,&alpha,BlockT,&dimL,BlockDiag,&dimR,&beta,workmemLR,&dimL);

         //mem * T^T --> storage
         char trans = 'T';
         alpha = 1.0;
         beta = 1.0; //add
         dgemm_(&notr,&trans,&dimL,&dimL,&dimR,&alpha,workmemLR,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimL);

      }
   }

}


