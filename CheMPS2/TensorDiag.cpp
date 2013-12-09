/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013 Sebastian Wouters

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


