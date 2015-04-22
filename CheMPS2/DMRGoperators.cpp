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
#include <string.h>
#include <iostream>
#include <sstream>

#include "DMRG.h"

using std::cout;
using std::endl;

void CheMPS2::DMRG::updateMovingRightSafeFirstTime(const int cnt){

   if (isAllocated[cnt]==2){
      deleteTensors(cnt, false);
      isAllocated[cnt]=0;
   }
   if (isAllocated[cnt]==0){
      allocateTensors(cnt, true);
      isAllocated[cnt]=1;
   }
   updateMovingRight(cnt);
   
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){
      if (cnt>0){
         if (isAllocated[cnt-1]==1){
            storeOperators(cnt-1, true);
            deleteTensors(cnt-1, true);
            isAllocated[cnt-1]=0;
         }
      }
   }

}

void CheMPS2::DMRG::updateMovingRightSafe(const int cnt){

   if (isAllocated[cnt]==2){
      deleteTensors(cnt, false);
      isAllocated[cnt]=0;
   }
   if (isAllocated[cnt]==0){
      allocateTensors(cnt, true);
      isAllocated[cnt]=1;
   }
   updateMovingRight(cnt);
   
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){
      if (cnt>0){
         if (isAllocated[cnt-1]==1){
            storeOperators(cnt-1, true);
            deleteTensors(cnt-1, true);
            isAllocated[cnt-1]=0;
         }
      }
      if (cnt+1<L-1){
         if (isAllocated[cnt+1]==2){
            deleteTensors(cnt+1, false);
            isAllocated[cnt+1]=0;
         }
      }
      if (cnt+2<L-1){
         if (isAllocated[cnt+2]==1){
            deleteTensors(cnt+2, true);
            isAllocated[cnt+2]=0;
         }
         if (isAllocated[cnt+2]==0){
            allocateTensors(cnt+2, false);
            isAllocated[cnt+2]=2;
         }
         loadOperators(cnt+2, false);
      }
   }

}

void CheMPS2::DMRG::updateMovingLeftSafe(const int cnt){

   if (isAllocated[cnt]==1){
      deleteTensors(cnt, true);
      isAllocated[cnt]=0;
   }
   if (isAllocated[cnt]==0){
      allocateTensors(cnt, false);
      isAllocated[cnt]=2;
   }
   updateMovingLeft(cnt);
   
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){
      if (cnt+1<L-1){
         if (isAllocated[cnt+1]==2){
            storeOperators(cnt+1, false);
            deleteTensors(cnt+1, false);
            isAllocated[cnt+1]=0;
         }
      }
      if (cnt-1>=0){
         if (isAllocated[cnt-1]==1){
            deleteTensors(cnt-1, true);
            isAllocated[cnt-1]=0;
         }
      }
      if (cnt-2>=0){
         if (isAllocated[cnt-2]==2){
            deleteTensors(cnt-2, false);
            isAllocated[cnt-2]=0;
         }
         if (isAllocated[cnt-2]==0){
            allocateTensors(cnt-2, true);
            isAllocated[cnt-2]=1;
         }
         loadOperators(cnt-2, true);
      }
   }

}

void CheMPS2::DMRG::updateMovingLeftSafe2DM(const int cnt){

   if (isAllocated[cnt]==1){
      deleteTensors(cnt, true);
      isAllocated[cnt]=0;
   }
   if (isAllocated[cnt]==0){
      allocateTensors(cnt, false);
      isAllocated[cnt]=2;
   }
   updateMovingLeft(cnt);
   
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){
      if (cnt+1<L-1){
         if (isAllocated[cnt+1]==2){
            deleteTensors(cnt+1, false);
            isAllocated[cnt+1]=0;
         }
      }
      if (cnt-1>=0){
         if (isAllocated[cnt-1]==2){
            deleteTensors(cnt-1, false);
            isAllocated[cnt-1]=0;
         }
         if (isAllocated[cnt-1]==0){
            allocateTensors(cnt-1, true);
            isAllocated[cnt-1]=1;
         }
         loadOperators(cnt-1, true);
      }
   }

}

void CheMPS2::DMRG::deleteAllBoundaryOperators(){

   for (int cnt=0; cnt<L-1; cnt++){
      if (isAllocated[cnt]==1){ deleteTensors(cnt, true); }
      if (isAllocated[cnt]==2){ deleteTensors(cnt, false); }
      isAllocated[cnt] = 0;
   }

}

int CheMPS2::DMRG::trianglefunction(const int k, const int glob){

   int cnt2tilde = 1;
   while(cnt2tilde*(cnt2tilde+1)/2 <= glob){ cnt2tilde++; }
   return k - cnt2tilde;
   
}

void CheMPS2::DMRG::updateMovingRight(const int index){

   const int dimL = denBK->gMaxDimAtBound(index);
   const int dimR = denBK->gMaxDimAtBound(index+1);
   
   #pragma omp parallel
   {
   
      double * workmem = new double[dimL*dimR];

      //Ltensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<(index+1) ; cnt2++){
         if (cnt2==0){
            Ltensors[index][cnt2]->makenew(MPS[index]);
         } else {
            Ltensors[index][cnt2]->update( Ltensors[index-1][cnt2-1] , MPS[index] , workmem );
         }
      }
      
      //Two-operator tensors
      const int k1 = index+1;
      const int upperbound1 = k1*(k1+1)/2;
      //After this parallel region, WAIT because F0,F1,S0,S1[index][cnt2][cnt3==0] is required for the complementary operators
      #pragma omp for schedule(static)
      for (int glob=0; glob<upperbound1; glob++){
         const int cnt2 = trianglefunction(k1,glob);
         const int cnt3 = glob - (k1-1-cnt2)*(k1-cnt2)/2;
         if (cnt3==0){
            if (cnt2==0){
               F0tensors[index][cnt2][cnt3]->makenew(MPS[index]);
               F1tensors[index][cnt2][cnt3]->makenew(MPS[index]);
               S0tensors[index][cnt2][cnt3]->makenew(MPS[index]);
               //S1[index][0][cnt3] doesn't exist
            } else {
               F0tensors[index][cnt2][cnt3]->makenew(Ltensors[index-1][cnt2-1],MPS[index],workmem);
               F1tensors[index][cnt2][cnt3]->makenew(Ltensors[index-1][cnt2-1],MPS[index],workmem);
               S0tensors[index][cnt2][cnt3]->makenew(Ltensors[index-1][cnt2-1],MPS[index],workmem);
               S1tensors[index][cnt2][cnt3]->makenew(Ltensors[index-1][cnt2-1],MPS[index],workmem);
            }
         } else {
            F0tensors[index][cnt2][cnt3]->update(F0tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
            F1tensors[index][cnt2][cnt3]->update(F1tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
            S0tensors[index][cnt2][cnt3]->update(S0tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
            if (cnt2>0){ S1tensors[index][cnt2][cnt3]->update(S1tensors[index-1][cnt2][cnt3-1],MPS[index],workmem); }
         }
      }
      
      //Complementary two-operator tensors
      const int k2 = L-1-index;
      const int upperbound2 = k2*(k2+1)/2;
      #pragma omp for schedule(static) nowait
      for (int glob=0; glob<upperbound2; glob++){
         const int cnt2 = trianglefunction(k2,glob);
         const int cnt3 = glob - (k2-1-cnt2)*(k2-cnt2)/2;
         if (index==0){
            Atensors[index][cnt2][cnt3]->ClearStorage();
            if (cnt2>0){ Btensors[index][cnt2][cnt3]->ClearStorage(); }
            Ctensors[index][cnt2][cnt3]->ClearStorage();
            Dtensors[index][cnt2][cnt3]->ClearStorage();
         } else {
            Atensors[index][cnt2][cnt3]->update(Atensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
            if (cnt2>0){ Btensors[index][cnt2][cnt3]->update(Btensors[index-1][cnt2][cnt3+1],MPS[index],workmem); }
            Ctensors[index][cnt2][cnt3]->update(Ctensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
            Dtensors[index][cnt2][cnt3]->update(Dtensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
         }
         for (int num=0; num<(index+1); num++){
            if ( Atensors[index][cnt2][cnt3]->gIdiff() == S0tensors[index][num][0]->gIdiff() ){ //Then the matrix elements are not 0 due to symm.
               
               double alpha = Prob->gMxElement(index-num,index,index+1+cnt3,index+1+cnt3+cnt2);
               if ((cnt2==0) && (num==0)) alpha *= 0.5;
               if ((cnt2>0) && (num>0)) alpha += Prob->gMxElement(index-num,index,index+1+cnt2+cnt3,index+1+cnt3);
               Atensors[index][cnt2][cnt3]->AddATerm(alpha,S0tensors[index][num][0]);

               alpha = 2*Prob->gMxElement(index-num,index+1+cnt3,index,index+1+cnt2+cnt3) - Prob->gMxElement(index-num,index,index+1+cnt2+cnt3,index+1+cnt3);
               Ctensors[index][cnt2][cnt3]->AddATerm(alpha,F0tensors[index][num][0]);
               
               alpha = - Prob->gMxElement(index-num,index,index+1+cnt2+cnt3,index+1+cnt3);
               Dtensors[index][cnt2][cnt3]->AddATerm(alpha,F1tensors[index][num][0]);
               
               if (num>0){
                  if (cnt2>0){
                     alpha = Prob->gMxElement(index-num,index,index+1+cnt3,index+1+cnt3+cnt2) - Prob->gMxElement(index-num,index,index+1+cnt2+cnt3,index+1+cnt3);
                     Btensors[index][cnt2][cnt3]->AddATerm(alpha,S1tensors[index][num][0]);
                  }
                  
                  alpha = 2*Prob->gMxElement(index-num,index+1+cnt3,index,index+1+cnt2+cnt3) - Prob->gMxElement(index-num,index,index+1+cnt3,index+1+cnt2+cnt3);
                  Ctensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F0tensors[index][num][0]);
                  
                  alpha = - Prob->gMxElement(index-num,index,index+1+cnt3,index+1+cnt2+cnt3);
                  Dtensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F1tensors[index][num][0]);
               }
            }
         }
      }
      
      //Qtensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<L-1-index ; cnt2++){
         if (index==0){
            Qtensors[index][cnt2]->ClearStorage();
            Qtensors[index][cnt2]->AddTermSimple(MPS[index]);
         } else {
            double * workmemBIS = new double[dimL*dimL];
            Qtensors[index][cnt2]->update(Qtensors[index-1][cnt2+1],MPS[index],workmem);
            Qtensors[index][cnt2]->AddTermSimple(MPS[index]);
            Qtensors[index][cnt2]->AddTermsL(Ltensors[index-1],MPS[index], workmemBIS, workmem);
            Qtensors[index][cnt2]->AddTermsAB(Atensors[index-1][cnt2+1][0], Btensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
            Qtensors[index][cnt2]->AddTermsCD(Ctensors[index-1][cnt2+1][0], Dtensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
            delete [] workmemBIS;
         }
      }
      
      delete [] workmem;
   
   }
   
   //Xtensors
   if (index==0){
      Xtensors[index]->update(MPS[index]);
   } else {
      Xtensors[index]->update(MPS[index], Ltensors[index-1], Xtensors[index-1], Qtensors[index-1][0], Atensors[index-1][0][0], Ctensors[index-1][0][0], Dtensors[index-1][0][0]);
   }
   
   //Otensors
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         if (index==0){
            Exc_Overlaps[state][index]->update(Exc_MPSs[state][index],MPS[index]);
         } else {
            Exc_Overlaps[state][index]->update(Exc_MPSs[state][index],MPS[index],Exc_Overlaps[state][index-1]);
         }
      }
   }
   
}

void CheMPS2::DMRG::updateMovingLeft(const int index){

   const int dimL = denBK->gMaxDimAtBound(index+1);
   const int dimR = denBK->gMaxDimAtBound(index+2);
   
   #pragma omp parallel
   {
   
      double * workmem = new double[dimL*dimR];

      //Ltensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
         if (cnt2==0){
            Ltensors[index][cnt2]->makenew(MPS[index+1]);
         } else {
            Ltensors[index][cnt2]->update( Ltensors[index+1][cnt2-1] , MPS[index+1] , workmem );
         }
      }
      
      //Two-operator tensors
      const int k1 = L-1-index;
      const int upperbound1 = k1*(k1+1)/2;
      //After this parallel region, WAIT because F0,F1,S0,S1[index][cnt2][cnt3==0] is required for the complementary operators
      #pragma omp for schedule(static)
      for (int glob=0; glob<upperbound1; glob++){
         const int cnt2 = trianglefunction(k1,glob);
         const int cnt3 = glob - (k1-1-cnt2)*(k1-cnt2)/2;
         if (cnt3==0){
            if (cnt2==0){
               F0tensors[index][cnt2][cnt3]->makenew(MPS[index+1]);
               F1tensors[index][cnt2][cnt3]->makenew(MPS[index+1]);
               S0tensors[index][cnt2][cnt3]->makenew(MPS[index+1]);
               //S1[index][0] doesn't exist
            } else {
               F0tensors[index][cnt2][cnt3]->makenew(Ltensors[index+1][cnt2-1],MPS[index+1],workmem);
               F1tensors[index][cnt2][cnt3]->makenew(Ltensors[index+1][cnt2-1],MPS[index+1],workmem);
               S0tensors[index][cnt2][cnt3]->makenew(Ltensors[index+1][cnt2-1],MPS[index+1],workmem);
               S1tensors[index][cnt2][cnt3]->makenew(Ltensors[index+1][cnt2-1],MPS[index+1],workmem);
            }
         } else {
            F0tensors[index][cnt2][cnt3]->update(F0tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
            F1tensors[index][cnt2][cnt3]->update(F1tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
            S0tensors[index][cnt2][cnt3]->update(S0tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
            if (cnt2>0){ S1tensors[index][cnt2][cnt3]->update(S1tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem); }
         }
      }
         
      //Complementary two-operator tensors
      const int k2 = index+1;
      const int upperbound2 = k2*(k2+1)/2;
      #pragma omp for schedule(static) nowait
      for (int glob=0; glob<upperbound2; glob++){
         const int cnt2 = trianglefunction(k2,glob);
         const int cnt3 = glob - (k2-1-cnt2)*(k2-cnt2)/2;
         if (index==L-2){
            Atensors[index][cnt2][cnt3]->ClearStorage();
            if (cnt2>0){ Btensors[index][cnt2][cnt3]->ClearStorage(); }
            Ctensors[index][cnt2][cnt3]->ClearStorage();
            Dtensors[index][cnt2][cnt3]->ClearStorage();
         } else {
            Atensors[index][cnt2][cnt3]->update(Atensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
            if (cnt2>0){ Btensors[index][cnt2][cnt3]->update(Btensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem); }
            Ctensors[index][cnt2][cnt3]->update(Ctensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
            Dtensors[index][cnt2][cnt3]->update(Dtensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
         }
         for (int num=0; num<L-index-1; num++){
            if ( Atensors[index][cnt2][cnt3]->gIdiff() == S0tensors[index][num][0]->gIdiff() ){ //Then the matrix elements are not 0 due to symm.
                 
               double alpha = Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1,index+1+num);
               if ((cnt2==0) && (num==0)) alpha *= 0.5;
               if ((cnt2>0) && (num>0)) alpha += Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1+num,index+1);
               Atensors[index][cnt2][cnt3]->AddATerm(alpha,S0tensors[index][num][0]);
               
               alpha = 2*Prob->gMxElement(index-cnt2-cnt3,index+1,index-cnt3,index+1+num) - Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1+num,index+1);
               Ctensors[index][cnt2][cnt3]->AddATerm(alpha,F0tensors[index][num][0]);
               
               alpha = - Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1+num,index+1);
               Dtensors[index][cnt2][cnt3]->AddATerm(alpha,F1tensors[index][num][0]);
               
               if (num>0){
                  if (cnt2>0){
                     alpha = Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1,index+1+num) - Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1+num,index+1);
                     Btensors[index][cnt2][cnt3]->AddATerm(alpha,S1tensors[index][num][0]);
                  }
                  
                  alpha = 2*Prob->gMxElement(index-cnt2-cnt3,index+1,index-cnt3,index+1+num) - Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1,index+1+num);
                  Ctensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F0tensors[index][num][0]);
                  
                  alpha = - Prob->gMxElement(index-cnt2-cnt3,index-cnt3,index+1,index+1+num);
                  Dtensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F1tensors[index][num][0]);
               }
            }
         }
      }
      
      //Qtensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<index+1 ; cnt2++){
         if (index==L-2){
            Qtensors[index][cnt2]->ClearStorage();
            Qtensors[index][cnt2]->AddTermSimple(MPS[index+1]);
         } else {
            double * workmemBIS = new double[dimR*dimR];
            Qtensors[index][cnt2]->update(Qtensors[index+1][cnt2+1],MPS[index+1],workmem);
            Qtensors[index][cnt2]->AddTermSimple(MPS[index+1]);
            Qtensors[index][cnt2]->AddTermsL(Ltensors[index+1],MPS[index+1], workmemBIS, workmem);
            Qtensors[index][cnt2]->AddTermsAB(Atensors[index+1][cnt2+1][0], Btensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
            Qtensors[index][cnt2]->AddTermsCD(Ctensors[index+1][cnt2+1][0], Dtensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
            delete [] workmemBIS;
         }
      }
      
      delete [] workmem;
   
   }
   
   //Xtensors
   if (index==L-2){
      Xtensors[index]->update(MPS[index+1]);
   } else {
      Xtensors[index]->update(MPS[index+1], Ltensors[index+1], Xtensors[index+1], Qtensors[index+1][0], Atensors[index+1][0][0], Ctensors[index+1][0][0], Dtensors[index+1][0][0]);
   }
   
   //Otensors
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         if (index==L-2){
            Exc_Overlaps[state][index]->update(Exc_MPSs[state][index+1],MPS[index+1]);
         } else {
            Exc_Overlaps[state][index]->update(Exc_MPSs[state][index+1],MPS[index+1],Exc_Overlaps[state][index+1]);
         }
      }
   }

}

void CheMPS2::DMRG::allocateTensors(const int index, const bool movingRight){

   if (movingRight){

      //Ltensors
      //To right: Ltens[cnt][cnt2] = operator on site cnt-cnt2; at boundary cnt+1
      Ltensors[index] = new TensorL * [index+1];
      for (int cnt2=0; cnt2<(index+1) ; cnt2++){ Ltensors[index][cnt2] = new TensorL(index+1,denBK->gIrrep(index-cnt2),movingRight,denBK); }
   
      //Two-operator tensors
      //To right: F0tens[cnt][cnt2][cnt3] = operators on sites cnt-cnt3-cnt2 and cnt-cnt3; at boundary cnt+1
      F0tensors[index] = new TensorF0 ** [index+1];
      F1tensors[index] = new TensorF1 ** [index+1];
      S0tensors[index] = new TensorS0 ** [index+1];
      S1tensors[index] = new TensorS1 ** [index+1];
      for (int cnt2=0; cnt2<(index+1) ; cnt2++){
         F0tensors[index][cnt2] = new TensorF0 * [index-cnt2+1];
         F1tensors[index][cnt2] = new TensorF1 * [index-cnt2+1];
         S0tensors[index][cnt2] = new TensorS0 * [index-cnt2+1];
         if (cnt2>0){ S1tensors[index][cnt2] = new TensorS1 * [index-cnt2+1]; }
         for (int cnt3=0; cnt3<(index-cnt2+1); cnt3++){
            const int Iprod = Irreps::directProd(denBK->gIrrep(index-cnt2-cnt3),denBK->gIrrep(index-cnt3));
            F0tensors[index][cnt2][cnt3] = new TensorF0(index+1,Iprod,movingRight,denBK);
            F1tensors[index][cnt2][cnt3] = new TensorF1(index+1,Iprod,movingRight,denBK);
            S0tensors[index][cnt2][cnt3] = new TensorS0(index+1,Iprod,movingRight,denBK);
            if (cnt2>0){ S1tensors[index][cnt2][cnt3] = new TensorS1(index+1,Iprod,movingRight,denBK); }
         }
      }
   
      //Complementary two-operator tensors
      //To right: Atens[cnt][cnt2][cnt3] = operators on sites cnt+1+cnt3 and cnt+1+cnt2+cnt3; at boundary cnt+1
      Atensors[index] = new TensorA ** [L-1-index];
      Btensors[index] = new TensorB ** [L-1-index];
      Ctensors[index] = new TensorC ** [L-1-index];
      Dtensors[index] = new TensorD ** [L-1-index];
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
         Atensors[index][cnt2] = new TensorA * [L-1-index-cnt2];
         if (cnt2>0){ Btensors[index][cnt2] = new TensorB * [L-1-index-cnt2]; }
         Ctensors[index][cnt2] = new TensorC * [L-1-index-cnt2];
         Dtensors[index][cnt2] = new TensorD * [L-1-index-cnt2];
         for (int cnt3=0; cnt3<L-1-index-cnt2; cnt3++){
            const int Idiff = Irreps::directProd(denBK->gIrrep(index+1+cnt2+cnt3),denBK->gIrrep(index+1+cnt3));
            Atensors[index][cnt2][cnt3] = new TensorA(index+1,Idiff,movingRight,denBK);
            if (cnt2>0){ Btensors[index][cnt2][cnt3] = new TensorB(index+1,Idiff,movingRight,denBK); }
            Ctensors[index][cnt2][cnt3] = new TensorC(index+1,Idiff,movingRight,denBK);
            Dtensors[index][cnt2][cnt3] = new TensorD(index+1,Idiff,movingRight,denBK);
         }
      }
   
      //Qtensors
      //To right: Qtens[cnt][cnt2] = operator on site cnt+1+cnt2; at boundary cnt+1
      Qtensors[index] = new TensorQ * [L-1-index];
      for (int cnt2=0; cnt2<L-1-index ; cnt2++){
         Qtensors[index][cnt2] = new TensorQ(index+1,denBK->gIrrep(index+1+cnt2),movingRight,denBK,Prob,index+1+cnt2);
      }
   
      //Xtensors
      Xtensors[index] = new TensorX(index+1,movingRight,denBK,Prob);
      
      //Otensors
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            Exc_Overlaps[state][index] = new TensorO(index+1,movingRight,Exc_BKs[state],denBK,Prob);
         }
      }
   
   } else {
   
      //Ltensors
      //To left: Ltens[cnt][cnt2] = operator on site cnt+1+cnt2; at boundary cnt+1
      Ltensors[index] = new TensorL * [L-1-index];
      for (int cnt2=0; cnt2<L-1-index; cnt2++){ Ltensors[index][cnt2] = new TensorL(index+1,denBK->gIrrep(index+1+cnt2),movingRight,denBK); }
   
      //Two-operator tensors
      //To left: F0tens[cnt][cnt2][cnt3] = operators on sites cnt+1+cnt3 and cnt+1+cnt3+cnt2; at boundary cnt+1
      F0tensors[index] = new TensorF0 ** [L-1-index];
      F1tensors[index] = new TensorF1 ** [L-1-index];
      S0tensors[index] = new TensorS0 ** [L-1-index];
      S1tensors[index] = new TensorS1 ** [L-1-index];
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
         F0tensors[index][cnt2] = new TensorF0 * [L-1-index-cnt2];
         F1tensors[index][cnt2] = new TensorF1 * [L-1-index-cnt2];
         S0tensors[index][cnt2] = new TensorS0 * [L-1-index-cnt2];
         if (cnt2>0){ S1tensors[index][cnt2] = new TensorS1 * [L-1-index-cnt2]; }
         for (int cnt3=0; cnt3<L-1-index-cnt2; cnt3++){
            const int Iprod = Irreps::directProd(denBK->gIrrep(index+1+cnt3),denBK->gIrrep(index+1+cnt2+cnt3));
            F0tensors[index][cnt2][cnt3] = new TensorF0(index+1,Iprod,movingRight,denBK);
            F1tensors[index][cnt2][cnt3] = new TensorF1(index+1,Iprod,movingRight,denBK);
            S0tensors[index][cnt2][cnt3] = new TensorS0(index+1,Iprod,movingRight,denBK);
            if (cnt2>0){ S1tensors[index][cnt2][cnt3] = new TensorS1(index+1,Iprod,movingRight,denBK); }
         }
      }
   
      //Complementary two-operator tensors
      //To left: Atens[cnt][cnt2][cnt3] = operators on sites cnt-cnt2-cnt3 and cnt-cnt3; at boundary cnt+1
      Atensors[index] = new TensorA ** [index+1];
      Btensors[index] = new TensorB ** [index+1];
      Ctensors[index] = new TensorC ** [index+1];
      Dtensors[index] = new TensorD ** [index+1];
      for (int cnt2=0; cnt2<index+1; cnt2++){
         Atensors[index][cnt2] = new TensorA * [index + 1 - cnt2];
         if (cnt2>0){ Btensors[index][cnt2] = new TensorB * [index + 1 - cnt2]; }
         Ctensors[index][cnt2] = new TensorC * [index + 1 - cnt2];
         Dtensors[index][cnt2] = new TensorD * [index + 1 - cnt2];
         for (int cnt3=0; cnt3<index+1-cnt2; cnt3++){
            const int Idiff = Irreps::directProd(denBK->gIrrep(index-cnt2-cnt3),denBK->gIrrep(index-cnt3));
            Atensors[index][cnt2][cnt3] = new TensorA(index+1,Idiff,movingRight,denBK);
            if (cnt2>0){ Btensors[index][cnt2][cnt3] = new TensorB(index+1,Idiff,movingRight,denBK); }
            Ctensors[index][cnt2][cnt3] = new TensorC(index+1,Idiff,movingRight,denBK);
            Dtensors[index][cnt2][cnt3] = new TensorD(index+1,Idiff,movingRight,denBK);
         }
      }
   
      //Qtensors
      //To left: Qtens[cnt][cnt2] = operator on site cnt-cnt2; at boundary cnt+1
      Qtensors[index] = new TensorQ * [index+1];
      for (int cnt2=0; cnt2<index+1 ; cnt2++){ Qtensors[index][cnt2] = new TensorQ(index+1,denBK->gIrrep(index-cnt2),movingRight,denBK,Prob,index-cnt2); }
   
      //Xtensors
      Xtensors[index] = new TensorX(index+1,movingRight,denBK,Prob);
      
      //Otensors
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            Exc_Overlaps[state][index] = new TensorO(index+1,movingRight,Exc_BKs[state],denBK,Prob);
         }
      }
   
   }

}

void CheMPS2::DMRG::MY_HDF5_WRITE(const hid_t file_id, const std::string sPath, Tensor * theTensor){

   const int size = theTensor->gKappa2index(theTensor->gNKappa());
   if (size > 0){

      hid_t group_id          = H5Gcreate(file_id, sPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      hsize_t dimarray        = size;
      hid_t dataspace_id      = H5Screate_simple(1, &dimarray, NULL);
      hid_t dataset_id        = H5Dcreate(group_id, "tensorStorage", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theTensor->gStorage());

      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);

      H5Gclose(group_id);
   
   }

}

void CheMPS2::DMRG::storeOperators(const int index, const bool movingRight){

   const int Nbound = movingRight ? index+1 : L-1-index;
   const int Cbound = movingRight ? L-1-index : index+1;

   std::stringstream thefilename;
   thefilename << tempfolder << "/" << CheMPS2::DMRG_OPERATOR_storage_prefix << thePID << "_index_" << index << ".h5";
   
   //The hdf5 file
   hid_t file_id = H5Fcreate(thefilename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   
      //Ltensors
      for (int cnt2=0; cnt2<Nbound ; cnt2++){
         std::stringstream sstream;
         sstream << "/Ltensor_" << cnt2 ;
         MY_HDF5_WRITE(file_id, sstream.str(), Ltensors[index][cnt2]);
      }
   
      //Two-operator tensors
      for (int cnt2=0; cnt2<Nbound ; cnt2++){
         for (int cnt3=0; cnt3<Nbound-cnt2 ; cnt3++){
            
            std::stringstream sstream1;
            sstream1 << "/F0tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream1.str(), F0tensors[index][cnt2][cnt3]);
            
            std::stringstream sstream2;
            sstream2 << "/F1tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream2.str(), F1tensors[index][cnt2][cnt3]);
            
            std::stringstream sstream3;
            sstream3 << "/S0tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream3.str(), S0tensors[index][cnt2][cnt3]);
            
            if (cnt2>0){
               std::stringstream sstream4;
               sstream4 << "/S1tensor_" << cnt2 << "_" << cnt3 ;
               MY_HDF5_WRITE(file_id, sstream4.str(), S1tensors[index][cnt2][cnt3]);
            }
            
         }
      }
   
      //Complementary two-operator tensors
      for (int cnt2=0; cnt2<Cbound ; cnt2++){
         for (int cnt3=0; cnt3<Cbound-cnt2 ; cnt3++){
         
            std::stringstream sstream1;
            sstream1 << "/Atensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream1.str(), Atensors[index][cnt2][cnt3]);
            
            if (cnt2>0){
               std::stringstream sstream2;
               sstream2 << "/Btensor_" << cnt2 << "_" << cnt3 ;
               MY_HDF5_WRITE(file_id, sstream2.str(), Btensors[index][cnt2][cnt3]);
            }
            
            std::stringstream sstream3;
            sstream3 << "/Ctensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream3.str(), Ctensors[index][cnt2][cnt3]);
            
            std::stringstream sstream4;
            sstream4 << "/Dtensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_WRITE(file_id, sstream4.str(), Dtensors[index][cnt2][cnt3]);
            
         }
      }
   
      //Qtensors
      for (int cnt2=0; cnt2<Cbound ; cnt2++){
         std::stringstream sstream;
         sstream << "/Qtensor_" << cnt2 ;
         MY_HDF5_WRITE(file_id, sstream.str(), Qtensors[index][cnt2]);
      }
   
      //Xtensors
      std::string sPathX = "/Xtensor" ;
      MY_HDF5_WRITE(file_id, sPathX, Xtensors[index]);
      
      //Otensors
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            std::stringstream sstream;
            sstream << "/Otensor_" << state ;
            MY_HDF5_WRITE(file_id, sstream.str(), Exc_Overlaps[state][index]);
         }
      }

   H5Fclose(file_id);
   
}

void CheMPS2::DMRG::MY_HDF5_READ(const hid_t file_id, const std::string sPath, Tensor * theTensor){

   const int size = theTensor->gKappa2index(theTensor->gNKappa());
   if (size > 0){

      hid_t group_id    = H5Gopen(file_id, sPath.c_str(), H5P_DEFAULT);

      hid_t dataset_id  = H5Dopen(group_id, "tensorStorage", H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theTensor->gStorage());
      
      H5Dclose(dataset_id);

      H5Gclose(group_id);
      
   }

}

void CheMPS2::DMRG::loadOperators(const int index, const bool movingRight){

   const int Nbound = movingRight ? index+1 : L-1-index;
   const int Cbound = movingRight ? L-1-index : index+1;

   std::stringstream thefilename;
   thefilename << tempfolder << "/" << CheMPS2::DMRG_OPERATOR_storage_prefix << thePID << "_index_" << index << ".h5";
   
   //The hdf5 file
   hid_t file_id = H5Fopen(thefilename.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   
      //Ltensors
      for (int cnt2=0; cnt2<Nbound ; cnt2++){
         std::stringstream sstream;
         sstream << "/Ltensor_" << cnt2 ;
         MY_HDF5_READ(file_id, sstream.str(), Ltensors[index][cnt2]);
      }
      
      //Two-operator tensors
      for (int cnt2=0; cnt2<Nbound ; cnt2++){
         for (int cnt3=0; cnt3<Nbound-cnt2 ; cnt3++){
         
            std::stringstream sstream1;
            sstream1 << "/F0tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream1.str(), F0tensors[index][cnt2][cnt3]);
            
            std::stringstream sstream2;
            sstream2 << "/F1tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream2.str(), F1tensors[index][cnt2][cnt3]);
            
            std::stringstream sstream3;
            sstream3 << "/S0tensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream3.str(), S0tensors[index][cnt2][cnt3]);
            
            if (cnt2>0){
               std::stringstream sstream4;
               sstream4 << "/S1tensor_" << cnt2 << "_" << cnt3 ;
               MY_HDF5_READ(file_id, sstream4.str(), S1tensors[index][cnt2][cnt3]);
            }
                        
         }
      }
      
      //Complementary two-operator tensors
      for (int cnt2=0; cnt2<Cbound ; cnt2++){
         for (int cnt3=0; cnt3<Cbound-cnt2 ; cnt3++){
         
            std::stringstream sstream1;
            sstream1 << "/Atensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream1.str(), Atensors[index][cnt2][cnt3]);
            
            if (cnt2>0){
               std::stringstream sstream2;
               sstream2 << "/Btensor_" << cnt2 << "_" << cnt3 ;
               MY_HDF5_READ(file_id, sstream2.str(), Btensors[index][cnt2][cnt3]);
            }
            
            std::stringstream sstream3;
            sstream3 << "/Ctensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream3.str(), Ctensors[index][cnt2][cnt3]);
            
            std::stringstream sstream4;
            sstream4 << "/Dtensor_" << cnt2 << "_" << cnt3 ;
            MY_HDF5_READ(file_id, sstream4.str(), Dtensors[index][cnt2][cnt3]);
            
         }
      }
      
      //Qtensors
      for (int cnt2=0; cnt2<Cbound ; cnt2++){
         std::stringstream sstream;
         sstream << "/Qtensor_" << cnt2 ;
         MY_HDF5_READ(file_id, sstream.str(), Qtensors[index][cnt2]);
      }
      
      //Xtensors
      std::string sPathX = "/Xtensor";
      MY_HDF5_READ(file_id, sPathX, Xtensors[index]);
      
      //Otensors
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            std::stringstream sstream;
            sstream << "/Otensor_" << state ;
            MY_HDF5_READ(file_id, sstream.str(), Exc_Overlaps[state][index]);
         }
      }
   
   H5Fclose(file_id);

}

void CheMPS2::DMRG::deleteTensors(const int index, const bool movingRightOfTensors){

   const int upperBoundNormal = (( movingRightOfTensors)?(index+1):(L-1-index));
   const int upperBoundComple = ((!movingRightOfTensors)?(index+1):(L-1-index));
   
   //Ltensors
   for (int cnt2=0; cnt2<upperBoundNormal; cnt2++){ delete Ltensors[index][cnt2]; }
   delete [] Ltensors[index];
   
   //Two-operator tensors
   for (int cnt2=0; cnt2<upperBoundNormal; cnt2++){
      for (int cnt3=0; cnt3<upperBoundNormal-cnt2; cnt3++){
         delete F0tensors[index][cnt2][cnt3];
         delete F1tensors[index][cnt2][cnt3];
         delete S0tensors[index][cnt2][cnt3];
         if (cnt2>0){ delete S1tensors[index][cnt2][cnt3]; }
      }
      delete [] F0tensors[index][cnt2];
      delete [] F1tensors[index][cnt2];
      delete [] S0tensors[index][cnt2];
      if (cnt2>0){ delete [] S1tensors[index][cnt2]; }
   }
   delete [] F0tensors[index];
   delete [] F1tensors[index];
   delete [] S0tensors[index];
   delete [] S1tensors[index];
   
   //Complementary two-operator tensors
   for (int cnt2=0; cnt2<upperBoundComple; cnt2++){
      for (int cnt3=0; cnt3<upperBoundComple-cnt2; cnt3++){
         delete Atensors[index][cnt2][cnt3];
         if (cnt2>0){ delete Btensors[index][cnt2][cnt3]; }
         delete Ctensors[index][cnt2][cnt3];
         delete Dtensors[index][cnt2][cnt3];
      }
      delete [] Atensors[index][cnt2];
      if (cnt2>0){ delete [] Btensors[index][cnt2]; }
      delete [] Ctensors[index][cnt2];
      delete [] Dtensors[index][cnt2];
   }
   delete [] Atensors[index];
   delete [] Btensors[index];
   delete [] Ctensors[index];
   delete [] Dtensors[index];
   
   //Qtensors
   for (int cnt2=0; cnt2<upperBoundComple ; cnt2++){ delete Qtensors[index][cnt2]; }
   delete [] Qtensors[index];
   
   //Xtensors
   delete Xtensors[index];
   
   //Otensors
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){ delete Exc_Overlaps[state][index]; }
   }
   
}

void CheMPS2::DMRG::deleteStoredOperators(){

   std::stringstream temp;
   temp << "rm " << tempfolder << "/" << CheMPS2::DMRG_OPERATOR_storage_prefix << thePID << "*.h5";
   int info = system(temp.str().c_str());
   cout << "Info on DMRG::operators rm call to system: " << info << endl;

}

