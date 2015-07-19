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
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "DMRG.h"
#include "Lapack.h"
#include "MPIchemps2.h"

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
            OperatorsOnDisk(cnt-1, true, true);
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
            OperatorsOnDisk(cnt-1, true, true);
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
         OperatorsOnDisk(cnt+2, false, false);
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
            OperatorsOnDisk(cnt+1, false, true);
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
         OperatorsOnDisk(cnt-2, true, false);
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
         OperatorsOnDisk(cnt-1, true, false);
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

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const int dimL = denBK->gMaxDimAtBound(index);
   const int dimR = denBK->gMaxDimAtBound(index+1);
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   #pragma omp parallel
   {
   
      double * workmem = new double[dimL*dimR];

      //Ltensors : all processes own all Ltensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<(index+1) ; cnt2++){
         if (cnt2==0){
            Ltensors[index][cnt2]->makenew(MPS[index]);
         } else {
            Ltensors[index][cnt2]->update( Ltensors[index-1][cnt2-1] , MPS[index] , workmem );
         }
      }
      
      //Two-operator tensors : certain processes own certain two-operator tensors
      const int k1 = index+1;
      const int upperbound1 = (k1*(k1+1))/2;
      //After this parallel region, WAIT because F0,F1,S0,S1[index][cnt2][cnt3==0] is required for the complementary operators
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static)
      #endif
      for (int glob=0; glob<upperbound1; glob++){
         const int cnt2 = trianglefunction(k1,glob);
         const int cnt3 = glob - (k1-1-cnt2)*(k1-cnt2)/2;
         //Operator[index][cnt2][cnt3] corresponds to site indices (index-cnt3) and (index-cnt3-cnt2)
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = index-cnt3-cnt2;
         const int siteindex2 = index-cnt3;
         #endif
         if (cnt3==0){ //Every MPI process owns the Operator[index][cnt2][cnt3==0]
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
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf(  L, siteindex1, siteindex2 ) == MPIRANK )
            #endif
            {
               F0tensors[index][cnt2][cnt3]->update(F0tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
               F1tensors[index][cnt2][cnt3]->update(F1tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( siteindex1, siteindex2 ) == MPIRANK )
            #endif
            {
               S0tensors[index][cnt2][cnt3]->update(S0tensors[index-1][cnt2][cnt3-1],MPS[index],workmem);
               if (cnt2>0){ S1tensors[index][cnt2][cnt3]->update(S1tensors[index-1][cnt2][cnt3-1],MPS[index],workmem); }
            }
         }
      }
      
      //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
      const int k2 = L-1-index;
      const int upperbound2 = (k2*(k2+1))/2;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int glob=0; glob<upperbound2; glob++){
         const int cnt2 = trianglefunction(k2,glob);
         const int cnt3 = glob - (k2-1-cnt2)*(k2-cnt2)/2;
         //Operator[index][cnt2][cnt3] corresponds to site indices (index+1+cnt3) and (index+1+cnt2+cnt3)
         const int siteindex1 = index+1+cnt3;
         const int siteindex2 = index+1+cnt2+cnt3;
         const int irrep_prod = Irreps::directProd( denBK->gIrrep(siteindex1), denBK->gIrrep(siteindex2) );
         #ifdef CHEMPS2_MPI_COMPILATION
         const bool do_absigma = ( MPIchemps2::owner_absigma( siteindex1, siteindex2 ) == MPIRANK );
         const bool do_cdf     = ( MPIchemps2::owner_cdf(  L, siteindex1, siteindex2 ) == MPIRANK );
         #endif
         if (index==0){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_absigma )
            #endif
            {
               Atensors[index][cnt2][cnt3]->ClearStorage();
               if (cnt2>0){ Btensors[index][cnt2][cnt3]->ClearStorage(); }
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_cdf )
            #endif
            {
               Ctensors[index][cnt2][cnt3]->ClearStorage();
               Dtensors[index][cnt2][cnt3]->ClearStorage();
            }
         } else {
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_absigma )
            #endif
            {
               Atensors[index][cnt2][cnt3]->update(Atensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
               if (cnt2>0){ Btensors[index][cnt2][cnt3]->update(Btensors[index-1][cnt2][cnt3+1],MPS[index],workmem); }
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_cdf )
            #endif
            {
               Ctensors[index][cnt2][cnt3]->update(Ctensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
               Dtensors[index][cnt2][cnt3]->update(Dtensors[index-1][cnt2][cnt3+1],MPS[index],workmem);
            }
         }
         for (int num=0; num<(index+1); num++){
            if ( irrep_prod == S0tensors[index][num][0]->gIdiff() ){ //Then the matrix elements are not 0 due to symm.
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( do_absigma )
               #endif
               {
                  double alpha = Prob->gMxElement(index-num,index,siteindex1,siteindex2);
                  if ((cnt2==0) && (num==0)){ alpha *= 0.5; }
                  if ((cnt2>0) && (num>0)){ alpha += Prob->gMxElement(index-num,index,siteindex2,siteindex1); }
                  Atensors[index][cnt2][cnt3]->AddATerm(alpha, S0tensors[index][num][0]);
                  
                  if ((num>0) && (cnt2>0)){
                     alpha = Prob->gMxElement(index-num,index,siteindex1,siteindex2)
                           - Prob->gMxElement(index-num,index,siteindex2,siteindex1);
                     Btensors[index][cnt2][cnt3]->AddATerm(alpha,S1tensors[index][num][0]);
                  }
               }
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( do_cdf )
               #endif
               {
                  double alpha = 2 * Prob->gMxElement(index-num,siteindex1,index,siteindex2)
                                   - Prob->gMxElement(index-num,siteindex1,siteindex2,index);
                  Ctensors[index][cnt2][cnt3]->AddATerm(alpha,F0tensors[index][num][0]);
                  
                  alpha = - Prob->gMxElement(index-num,siteindex1,siteindex2,index); // Second line for Ctensors
                  Dtensors[index][cnt2][cnt3]->AddATerm(alpha,F1tensors[index][num][0]);
                  
                  if (num>0){
                     alpha = 2 * Prob->gMxElement(index-num,siteindex2,index,siteindex1)
                               - Prob->gMxElement(index-num,siteindex2,siteindex1,index);
                     Ctensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F0tensors[index][num][0]);
                     
                     alpha = - Prob->gMxElement(index-num,siteindex2,siteindex1,index); // Second line for Ctensors
                     Dtensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F1tensors[index][num][0]);
                  }
               }
            }
         }
      }
      
      //Qtensors : certain processes own certain Qtensors --- You don't want to locally parallellize when sending and receiving buffers!
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp single
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
      
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex = index+1+cnt2; //Corresponds to this site
         const int owner_q = MPIchemps2::owner_q( L, siteindex );
         #endif
         if ( index == 0 ){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( owner_q == MPIRANK )
            #endif
            {
               Qtensors[index][cnt2]->ClearStorage();
               Qtensors[index][cnt2]->AddTermSimple(MPS[index]);
            }
            
         } else {
         
            #ifdef CHEMPS2_MPI_COMPILATION
            const int owner_absigma = MPIchemps2::owner_absigma( index, siteindex );
            const int owner_cdf     = MPIchemps2::owner_cdf(  L, index, siteindex );
            if (( owner_q == owner_absigma ) && ( owner_q == owner_cdf ) && ( owner_q == MPIRANK )){ // No MPI needed
            #endif
         
               double * workmemBIS = new double[dimL*dimL];
               Qtensors[index][cnt2]->update(Qtensors[index-1][cnt2+1],MPS[index],workmem);
               Qtensors[index][cnt2]->AddTermSimple(MPS[index]);
               Qtensors[index][cnt2]->AddTermsL(Ltensors[index-1],MPS[index], workmemBIS, workmem);
               Qtensors[index][cnt2]->AddTermsAB(Atensors[index-1][cnt2+1][0], Btensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
               Qtensors[index][cnt2]->AddTermsCD(Ctensors[index-1][cnt2+1][0], Dtensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
               delete [] workmemBIS;
            
            #ifdef CHEMPS2_MPI_COMPILATION
            } else { //There's going to have to be some communication
            
               if (( owner_q == MPIRANK ) || ( owner_absigma == MPIRANK ) || ( owner_cdf == MPIRANK )){
               
                  TensorQ * tempQ = new TensorQ(index+1,denBK->gIrrep(siteindex),true,denBK,Prob,siteindex);
                  tempQ->ClearStorage();
                  
                  //Everyone creates his/her piece
                  double * workmemBIS = new double[dimL*dimL];
                  if ( owner_q == MPIRANK ){
                     tempQ->update(Qtensors[index-1][cnt2+1],MPS[index],workmem);
                     tempQ->AddTermSimple(MPS[index]);
                     tempQ->AddTermsL(Ltensors[index-1],MPS[index], workmemBIS, workmem);
                  }
                  if ( owner_absigma == MPIRANK ){
                     tempQ->AddTermsAB(Atensors[index-1][cnt2+1][0], Btensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
                  }
                  if ( owner_cdf == MPIRANK ){
                     tempQ->AddTermsCD(Ctensors[index-1][cnt2+1][0], Dtensors[index-1][cnt2+1][0], MPS[index], workmemBIS, workmem);
                  }
                  delete [] workmemBIS;
                  
                  //Add everything to owner_q's Qtensors[index][cnt2]: replace later with custom communication group?
                  int inc = 1;
                  int arraysize = tempQ->gKappa2index(tempQ->gNKappa());
                  double alpha = 1.0;
                  if ( owner_q == MPIRANK ){ dcopy_(&arraysize, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  if ( owner_q != owner_absigma ){
                     MPIchemps2::sendreceive_tensor( tempQ, owner_absigma, owner_q, 2*siteindex );
                     if ( owner_q == MPIRANK ){ daxpy_(&arraysize, &alpha, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  }
                  if (( owner_q != owner_cdf ) && ( owner_absigma != owner_cdf )){
                     MPIchemps2::sendreceive_tensor( tempQ, owner_cdf, owner_q, 2*siteindex+1 );
                     if ( owner_q == MPIRANK ){ daxpy_(&arraysize, &alpha, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  }
                  delete tempQ;
                  
               }
            }
            #endif
         }
      }
      
      delete [] workmem;
   
   }
   
   //Xtensors
   #ifdef CHEMPS2_MPI_COMPILATION
   const int owner_x = MPIchemps2::owner_x();
   #endif
   if ( index == 0 ){
   
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( owner_x == MPIRANK )
      #endif
      { Xtensors[index]->update(MPS[index]); }
      
   } else {
   
      #ifdef CHEMPS2_MPI_COMPILATION
      //Make sure that owner_x has all required tensors to construct X. Not as optimal as Q-tensor case, but easier hack.
      const int owner_q       = MPIchemps2::owner_q( L, index );
      const int owner_absigma = MPIchemps2::owner_absigma( index, index );
      const int owner_cdf     = MPIchemps2::owner_cdf( L, index, index );
      const int Idiff         = Irreps::directProd( denBK->gIrrep(index), denBK->gIrrep(index) ); // Will be num ^ num --> 0
      
      if ( owner_x != owner_q ){
         if ( owner_x == MPIRANK ){ Qtensors[index-1][0] = new TensorQ( index, denBK->gIrrep(index), true, denBK, Prob, index ); }
         if (( owner_x == MPIRANK ) || ( owner_q == MPIRANK )){ MPIchemps2::sendreceive_tensor( Qtensors[index-1][0], owner_q, owner_x, 3*L+3 ); }
      }
      
      if ( owner_x != owner_absigma ){
         if ( owner_x == MPIRANK ){ Atensors[index-1][0][0] = new TensorA( index, Idiff, true, denBK ); }
         if (( owner_x == MPIRANK ) || ( owner_absigma == MPIRANK )){ MPIchemps2::sendreceive_tensor( Atensors[index-1][0][0], owner_absigma, owner_x, 3*L+4 ); }
      }
      
      if ( owner_x != owner_cdf ){
         if ( owner_x == MPIRANK ){
            Ctensors[index-1][0][0] = new TensorC( index, Idiff, true, denBK );
            Dtensors[index-1][0][0] = new TensorD( index, Idiff, true, denBK );
         }
         if (( owner_x == MPIRANK ) || ( owner_cdf == MPIRANK )){
            MPIchemps2::sendreceive_tensor( Ctensors[index-1][0][0], owner_cdf, owner_x, 3*L+5 );
            MPIchemps2::sendreceive_tensor( Dtensors[index-1][0][0], owner_cdf, owner_x, 3*L+6 );
         }
      }
      
      if ( owner_x == MPIRANK ){
      #endif
      
      Xtensors[index]->update(MPS[index], Ltensors[index-1], Xtensors[index-1], Qtensors[index-1][0], Atensors[index-1][0][0], Ctensors[index-1][0][0], Dtensors[index-1][0][0]);
      
      #ifdef CHEMPS2_MPI_COMPILATION
         if ( owner_x != owner_q       ){ delete Qtensors[index-1][0];    }
         if ( owner_x != owner_absigma ){ delete Atensors[index-1][0][0]; }
         if ( owner_x != owner_cdf     ){ delete Ctensors[index-1][0][0];
                                          delete Dtensors[index-1][0][0]; }
      }
      #endif
      
   }
   
   //Otensors : certain processes own certain excitations
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
         #endif
         {
            if (index==0){
               Exc_Overlaps[state][index]->update(Exc_MPSs[state][index],MPS[index]);
            } else {
               Exc_Overlaps[state][index]->update(Exc_MPSs[state][index],MPS[index],Exc_Overlaps[state][index-1]);
            }
         }
      }
   }
   
   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_CALC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::updateMovingLeft(const int index){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const int dimL = denBK->gMaxDimAtBound(index+1);
   const int dimR = denBK->gMaxDimAtBound(index+2);
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   #pragma omp parallel
   {
   
      double * workmem = new double[dimL*dimR];

      //Ltensors : all processes own all Ltensors
      #pragma omp for schedule(static) nowait
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
         if (cnt2==0){
            Ltensors[index][cnt2]->makenew(MPS[index+1]);
         } else {
            Ltensors[index][cnt2]->update( Ltensors[index+1][cnt2-1] , MPS[index+1] , workmem );
         }
      }
      
      //Two-operator tensors : certain processes own certain two-operator tensors
      const int k1 = L-1-index;
      const int upperbound1 = k1*(k1+1)/2;
      //After this parallel region, WAIT because F0,F1,S0,S1[index][cnt2][cnt3==0] is required for the complementary operators
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static)
      #endif
      for (int glob=0; glob<upperbound1; glob++){
         const int cnt2 = trianglefunction(k1,glob);
         const int cnt3 = glob - (k1-1-cnt2)*(k1-cnt2)/2;
         //Operator[index][cnt2][cnt3] corresponds to site indices (index+1+cnt3) and (index+1+cnt2+cnt3)
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = index+1+cnt3;
         const int siteindex2 = index+1+cnt2+cnt3;
         #endif
         if (cnt3==0){ //Every MPI process owns the Operator[index][cnt2][cnt3==0]
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
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf(  L, siteindex1, siteindex2 ) == MPIRANK )
            #endif
            {
               F0tensors[index][cnt2][cnt3]->update(F0tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
               F1tensors[index][cnt2][cnt3]->update(F1tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( siteindex1, siteindex2 ) == MPIRANK )
            #endif
            {
               S0tensors[index][cnt2][cnt3]->update(S0tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem);
               if (cnt2>0){ S1tensors[index][cnt2][cnt3]->update(S1tensors[index+1][cnt2][cnt3-1],MPS[index+1],workmem); }
            }
         }
      }
         
      //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
      const int k2 = index+1;
      const int upperbound2 = k2*(k2+1)/2;
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp for schedule(dynamic)
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int glob=0; glob<upperbound2; glob++){
         const int cnt2 = trianglefunction(k2,glob);
         const int cnt3 = glob - (k2-1-cnt2)*(k2-cnt2)/2;
         //Operator[index][cnt2][cnt3] corresponds to site indices (index-cnt3) and (index-cnt2-cnt3)
         const int siteindex1  = index-cnt3-cnt2;
         const int siteindex2  = index-cnt3;
         const int irrep_prod  = Irreps::directProd(denBK->gIrrep(siteindex1),denBK->gIrrep(siteindex2));
         #ifdef CHEMPS2_MPI_COMPILATION
         const bool do_absigma = ( MPIchemps2::owner_absigma( siteindex1, siteindex2 ) == MPIRANK );
         const bool do_cdf     = ( MPIchemps2::owner_cdf(  L, siteindex1, siteindex2 ) == MPIRANK );
         #endif
         if (index==L-2){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_absigma )
            #endif
            {
               Atensors[index][cnt2][cnt3]->ClearStorage();
               if (cnt2>0){ Btensors[index][cnt2][cnt3]->ClearStorage(); }
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_cdf )
            #endif
            {
               Ctensors[index][cnt2][cnt3]->ClearStorage();
               Dtensors[index][cnt2][cnt3]->ClearStorage();
            }
         } else {
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_absigma )
            #endif
            {
               Atensors[index][cnt2][cnt3]->update(Atensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
               if (cnt2>0){ Btensors[index][cnt2][cnt3]->update(Btensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem); }
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( do_cdf )
            #endif
            {
               Ctensors[index][cnt2][cnt3]->update(Ctensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
               Dtensors[index][cnt2][cnt3]->update(Dtensors[index+1][cnt2][cnt3+1],MPS[index+1],workmem);
            }
         }
         for (int num=0; num<L-index-1; num++){
            if ( irrep_prod == S0tensors[index][num][0]->gIdiff() ){ //Then the matrix elements are not 0 due to symm.
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( do_absigma )
               #endif
               {
                  double alpha = Prob->gMxElement(siteindex1,siteindex2,index+1,index+1+num);
                  if ((cnt2==0) && (num==0)) alpha *= 0.5;
                  if ((cnt2>0) && (num>0)) alpha += Prob->gMxElement(siteindex1,siteindex2,index+1+num,index+1);
                  Atensors[index][cnt2][cnt3]->AddATerm(alpha,S0tensors[index][num][0]);
               
                  if ((num>0) && (cnt2>0)){
                     alpha = Prob->gMxElement(siteindex1,siteindex2,index+1,index+1+num)
                           - Prob->gMxElement(siteindex1,siteindex2,index+1+num,index+1);
                     Btensors[index][cnt2][cnt3]->AddATerm(alpha,S1tensors[index][num][0]);
                  }
               }
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( do_cdf )
               #endif
               {
                  double alpha = 2 * Prob->gMxElement(siteindex1,index+1,siteindex2,index+1+num)
                                   - Prob->gMxElement(siteindex1,index+1,index+1+num,siteindex2);
                  Ctensors[index][cnt2][cnt3]->AddATerm(alpha,F0tensors[index][num][0]);
                  
                  alpha = - Prob->gMxElement(siteindex1,index+1,index+1+num,siteindex2); // Second line for Ctensors
                  Dtensors[index][cnt2][cnt3]->AddATerm(alpha,F1tensors[index][num][0]);
                  
                  if (num>0){
                     alpha = 2 * Prob->gMxElement(siteindex1,index+1+num,siteindex2,index+1)
                               - Prob->gMxElement(siteindex1,index+1+num,index+1,siteindex2);
                     Ctensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F0tensors[index][num][0]);
                     
                     alpha = - Prob->gMxElement(siteindex1,index+1+num,index+1,siteindex2); // Second line for Ctensors
                     Dtensors[index][cnt2][cnt3]->AddATermTranspose(alpha,F1tensors[index][num][0]);
                  }
               }
            }
         }
      }
      
      //Qtensors : certain processes own certain Qtensors --- You don't want to locally parallellize when sending and receiving buffers!
      #ifdef CHEMPS2_MPI_COMPILATION
         #pragma omp single
      #else
         #pragma omp for schedule(static) nowait
      #endif
      for (int cnt2=0; cnt2<index+1; cnt2++){
      
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex = index-cnt2; //Corresponds to this site
         const int owner_q = MPIchemps2::owner_q( L, siteindex );
         #endif
         if ( index == L-2 ){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( owner_q == MPIRANK )
            #endif
            {
               Qtensors[index][cnt2]->ClearStorage();
               Qtensors[index][cnt2]->AddTermSimple(MPS[index+1]);
            }
         
         } else {
         
            #ifdef CHEMPS2_MPI_COMPILATION
            const int owner_absigma = MPIchemps2::owner_absigma( siteindex, index+1 );
            const int owner_cdf     = MPIchemps2::owner_cdf(  L, siteindex, index+1 );
            if (( owner_q == owner_absigma ) && ( owner_q == owner_cdf ) && ( owner_q == MPIRANK )){ // No MPI needed
            #endif
            
               double * workmemBIS = new double[dimR*dimR];
               Qtensors[index][cnt2]->update(Qtensors[index+1][cnt2+1],MPS[index+1],workmem);
               Qtensors[index][cnt2]->AddTermSimple(MPS[index+1]);
               Qtensors[index][cnt2]->AddTermsL(Ltensors[index+1],MPS[index+1], workmemBIS, workmem);
               Qtensors[index][cnt2]->AddTermsAB(Atensors[index+1][cnt2+1][0], Btensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
               Qtensors[index][cnt2]->AddTermsCD(Ctensors[index+1][cnt2+1][0], Dtensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
               delete [] workmemBIS;
            
            #ifdef CHEMPS2_MPI_COMPILATION
            } else { // There's going to have to be some communication
            
               if (( owner_q == MPIRANK ) || ( owner_absigma == MPIRANK ) || ( owner_cdf == MPIRANK )){
               
                  TensorQ * tempQ = new TensorQ(index+1,denBK->gIrrep(siteindex),false,denBK,Prob,siteindex);
                  tempQ->ClearStorage();
                  
                  //Everyone creates his/her piece
                  double * workmemBIS = new double[dimR*dimR];
                  if ( owner_q == MPIRANK ){
                     tempQ->update(Qtensors[index+1][cnt2+1],MPS[index+1],workmem);
                     tempQ->AddTermSimple(MPS[index+1]);
                     tempQ->AddTermsL(Ltensors[index+1],MPS[index+1], workmemBIS, workmem);
                  }
                  if ( owner_absigma == MPIRANK ){
                     tempQ->AddTermsAB(Atensors[index+1][cnt2+1][0], Btensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
                  }
                  if ( owner_cdf == MPIRANK ){
                     tempQ->AddTermsCD(Ctensors[index+1][cnt2+1][0], Dtensors[index+1][cnt2+1][0], MPS[index+1], workmemBIS, workmem);
                  }
                  delete [] workmemBIS;
                  
                  //Add everything to owner_q's Qtensors[index][cnt2]: replace later with custom communication group?
                  int inc = 1;
                  int arraysize = tempQ->gKappa2index(tempQ->gNKappa());
                  double alpha = 1.0;
                  if ( owner_q == MPIRANK ){ dcopy_(&arraysize, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  if ( owner_q != owner_absigma ){
                     MPIchemps2::sendreceive_tensor( tempQ, owner_absigma, owner_q, 2*siteindex );
                     if ( owner_q == MPIRANK ){ daxpy_(&arraysize, &alpha, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  }
                  if (( owner_q != owner_cdf ) && ( owner_absigma != owner_cdf )){
                     MPIchemps2::sendreceive_tensor( tempQ, owner_cdf, owner_q, 2*siteindex+1 );
                     if ( owner_q == MPIRANK ){ daxpy_(&arraysize, &alpha, tempQ->gStorage(), &inc, Qtensors[index][cnt2]->gStorage(), &inc); }
                  }
                  delete tempQ;
                  
               }
            }
            #endif
         }
      }
      
      delete [] workmem;
   
   }
   
   //Xtensors
   #ifdef CHEMPS2_MPI_COMPILATION
   const int owner_x = MPIchemps2::owner_x();
   #endif
   if ( index == L-2 ){
   
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( owner_x == MPIRANK )
      #endif
      { Xtensors[index]->update(MPS[index+1]); }
      
   } else {
   
      #ifdef CHEMPS2_MPI_COMPILATION
      //Make sure that owner_x has all required tensors to construct X. Not as optimal as Q-tensor case, but easier hack.
      const int owner_q       = MPIchemps2::owner_q( L, index+1 );
      const int owner_absigma = MPIchemps2::owner_absigma( index+1, index+1 );
      const int owner_cdf     = MPIchemps2::owner_cdf( L, index+1, index+1 );
      const int Idiff         = Irreps::directProd( denBK->gIrrep(index+1), denBK->gIrrep(index+1) ); // Will be num ^ num --> 0
      
      if ( owner_x != owner_q ){
         if ( owner_x == MPIRANK ){ Qtensors[index+1][0] = new TensorQ( index+2, denBK->gIrrep(index+1), false, denBK, Prob, index+1 ); }
         if (( owner_x == MPIRANK ) || ( owner_q == MPIRANK )){ MPIchemps2::sendreceive_tensor( Qtensors[index+1][0], owner_q, owner_x, 3*L+3 ); }
      }
      
      if ( owner_x != owner_absigma ){
         if ( owner_x == MPIRANK ){ Atensors[index+1][0][0] = new TensorA( index+2, Idiff, false, denBK ); }
         if (( owner_x == MPIRANK ) || ( owner_absigma == MPIRANK )){ MPIchemps2::sendreceive_tensor( Atensors[index+1][0][0], owner_absigma, owner_x, 3*L+4 ); }
      }
      
      if ( owner_x != owner_cdf ){
         if ( owner_x == MPIRANK ){
            Ctensors[index+1][0][0] = new TensorC( index+2, Idiff, false, denBK );
            Dtensors[index+1][0][0] = new TensorD( index+2, Idiff, false, denBK );
         }
         if (( owner_x == MPIRANK ) || ( owner_cdf == MPIRANK )){
            MPIchemps2::sendreceive_tensor( Ctensors[index+1][0][0], owner_cdf, owner_x, 3*L+5 );
            MPIchemps2::sendreceive_tensor( Dtensors[index+1][0][0], owner_cdf, owner_x, 3*L+6 );
         }
      }
      
      if ( owner_x == MPIRANK ){
      #endif
      
      Xtensors[index]->update(MPS[index+1], Ltensors[index+1], Xtensors[index+1], Qtensors[index+1][0], Atensors[index+1][0][0], Ctensors[index+1][0][0], Dtensors[index+1][0][0]);
      
      #ifdef CHEMPS2_MPI_COMPILATION
         if ( owner_x != owner_q       ){ delete Qtensors[index+1][0];    }
         if ( owner_x != owner_absigma ){ delete Atensors[index+1][0][0]; }
         if ( owner_x != owner_cdf     ){ delete Ctensors[index+1][0][0];
                                          delete Dtensors[index+1][0][0]; }
      }
      #endif
      
   }
   
   //Otensors
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
         #endif
         {
            if (index==L-2){
               Exc_Overlaps[state][index]->update(Exc_MPSs[state][index+1],MPS[index+1]);
            } else {
               Exc_Overlaps[state][index]->update(Exc_MPSs[state][index+1],MPS[index+1],Exc_Overlaps[state][index+1]);
            }
         }
      }
   }

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_CALC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::allocateTensors(const int index, const bool movingRight){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   if (movingRight){

      //Ltensors : all processes own all Ltensors
      //To right: Ltens[cnt][cnt2] = operator on site cnt-cnt2; at boundary cnt+1
      Ltensors[index] = new TensorL * [index+1];
      for (int cnt2=0; cnt2<(index+1) ; cnt2++){ Ltensors[index][cnt2] = new TensorL(index+1,denBK->gIrrep(index-cnt2),movingRight,denBK); }
   
      //Two-operator tensors : certain processes own certain two-operator tensors
      //To right: F0tens[cnt][cnt2][cnt3] = operators on sites cnt-cnt3-cnt2 and cnt-cnt3; at boundary cnt+1
      F0tensors[index] = new TensorF0 ** [index+1];
      F1tensors[index] = new TensorF1 ** [index+1];
      S0tensors[index] = new TensorS0 ** [index+1];
      S1tensors[index] = new TensorS1 ** [index+1];
      for (int cnt2=0; cnt2<(index+1); cnt2++){
         F0tensors[index][cnt2] = new TensorF0 * [index-cnt2+1];
         F1tensors[index][cnt2] = new TensorF1 * [index-cnt2+1];
         S0tensors[index][cnt2] = new TensorS0 * [index-cnt2+1];
         if (cnt2>0){ S1tensors[index][cnt2] = new TensorS1 * [index-cnt2+1]; }
         for (int cnt3=0; cnt3<(index-cnt2+1); cnt3++){
            const int Iprod = Irreps::directProd(denBK->gIrrep(index-cnt2-cnt3),denBK->gIrrep(index-cnt3));
            #ifdef CHEMPS2_MPI_COMPILATION
            if (( cnt3 == 0 ) || ( MPIchemps2::owner_cdf(L, index-cnt2-cnt3, index-cnt3) == MPIRANK )){
            #endif
               F0tensors[index][cnt2][cnt3] = new TensorF0(index+1,Iprod,movingRight,denBK);
               F1tensors[index][cnt2][cnt3] = new TensorF1(index+1,Iprod,movingRight,denBK);
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               F0tensors[index][cnt2][cnt3] = NULL;
               F1tensors[index][cnt2][cnt3] = NULL;
            }
            if (( cnt3 == 0 ) || ( MPIchemps2::owner_absigma(index-cnt2-cnt3, index-cnt3) == MPIRANK )){
            #endif
               S0tensors[index][cnt2][cnt3] = new TensorS0(index+1,Iprod,movingRight,denBK);
               if (cnt2>0){ S1tensors[index][cnt2][cnt3] = new TensorS1(index+1,Iprod,movingRight,denBK); }
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               S0tensors[index][cnt2][cnt3] = NULL;
               if (cnt2>0){ S1tensors[index][cnt2][cnt3] = NULL; }
            }
            #endif
         }
      }
   
      //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
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
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma(index+1+cnt3, index+1+cnt2+cnt3) == MPIRANK ){
            #endif
               Atensors[index][cnt2][cnt3] = new TensorA(index+1,Idiff,movingRight,denBK);
               if (cnt2>0){ Btensors[index][cnt2][cnt3] = new TensorB(index+1,Idiff,movingRight,denBK); }
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               Atensors[index][cnt2][cnt3] = NULL;
               if (cnt2>0){ Btensors[index][cnt2][cnt3] = NULL; }
            }
            if ( MPIchemps2::owner_cdf(L, index+1+cnt3, index+1+cnt2+cnt3) == MPIRANK ){
            #endif
               Ctensors[index][cnt2][cnt3] = new TensorC(index+1,Idiff,movingRight,denBK);
               Dtensors[index][cnt2][cnt3] = new TensorD(index+1,Idiff,movingRight,denBK);
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               Ctensors[index][cnt2][cnt3] = NULL;
               Dtensors[index][cnt2][cnt3] = NULL;
            }
            #endif
         }
      }
   
      //Qtensors
      //To right: Qtens[cnt][cnt2] = operator on site cnt+1+cnt2; at boundary cnt+1
      Qtensors[index] = new TensorQ * [L-1-index];
      for (int cnt2=0; cnt2<L-1-index; cnt2++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_q( L, index+1+cnt2 ) == MPIRANK ){
         #endif
            Qtensors[index][cnt2] = new TensorQ(index+1,denBK->gIrrep(index+1+cnt2),movingRight,denBK,Prob,index+1+cnt2);
         #ifdef CHEMPS2_MPI_COMPILATION
         } else { Qtensors[index][cnt2] = NULL; }
         #endif
      }
   
      //Xtensors : a certain process owns the Xtensors
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_x() == MPIRANK ){
      #endif
         Xtensors[index] = new TensorX(index+1,movingRight,denBK,Prob);
      #ifdef CHEMPS2_MPI_COMPILATION
      } else { Xtensors[index] = NULL; }
      #endif
      
      //Otensors : certain processes own certain excitations
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
            #endif
            { Exc_Overlaps[state][index] = new TensorO(index+1,movingRight,Exc_BKs[state],denBK,Prob); }
         }
      }
   
   } else {
   
      //Ltensors : all processes own all Ltensors
      //To left: Ltens[cnt][cnt2] = operator on site cnt+1+cnt2; at boundary cnt+1
      Ltensors[index] = new TensorL * [L-1-index];
      for (int cnt2=0; cnt2<L-1-index; cnt2++){ Ltensors[index][cnt2] = new TensorL(index+1,denBK->gIrrep(index+1+cnt2),movingRight,denBK); }
   
      //Two-operator tensors : certain processes own certain two-operator tensors
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
            #ifdef CHEMPS2_MPI_COMPILATION
            if (( cnt3 == 0 ) || ( MPIchemps2::owner_cdf(L, index+1+cnt3, index+1+cnt2+cnt3) == MPIRANK )){
            #endif
               F0tensors[index][cnt2][cnt3] = new TensorF0(index+1,Iprod,movingRight,denBK);
               F1tensors[index][cnt2][cnt3] = new TensorF1(index+1,Iprod,movingRight,denBK);
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               F0tensors[index][cnt2][cnt3] = NULL;
               F1tensors[index][cnt2][cnt3] = NULL;
            }
            if (( cnt3 == 0 ) || ( MPIchemps2::owner_absigma(index+1+cnt3, index+1+cnt2+cnt3) == MPIRANK )){
            #endif
               S0tensors[index][cnt2][cnt3] = new TensorS0(index+1,Iprod,movingRight,denBK);
               if (cnt2>0){ S1tensors[index][cnt2][cnt3] = new TensorS1(index+1,Iprod,movingRight,denBK); }
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               S0tensors[index][cnt2][cnt3] = NULL;
               if (cnt2>0){ S1tensors[index][cnt2][cnt3] = NULL; }
            }
            #endif
         }
      }
   
      //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
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
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma(index-cnt2-cnt3, index-cnt3) == MPIRANK ){
            #endif
               Atensors[index][cnt2][cnt3] = new TensorA(index+1,Idiff,movingRight,denBK);
               if (cnt2>0){ Btensors[index][cnt2][cnt3] = new TensorB(index+1,Idiff,movingRight,denBK); }
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               Atensors[index][cnt2][cnt3] = NULL;
               if (cnt2>0){ Btensors[index][cnt2][cnt3] = NULL; }
            }
            if ( MPIchemps2::owner_cdf(L, index-cnt2-cnt3, index-cnt3) == MPIRANK ){
            #endif
               Ctensors[index][cnt2][cnt3] = new TensorC(index+1,Idiff,movingRight,denBK);
               Dtensors[index][cnt2][cnt3] = new TensorD(index+1,Idiff,movingRight,denBK);
            #ifdef CHEMPS2_MPI_COMPILATION
            } else {
               Ctensors[index][cnt2][cnt3] = NULL;
               Dtensors[index][cnt2][cnt3] = NULL;
            }
            #endif
         }
      }
   
      //Qtensors : certain processes own certain Qtensors
      //To left: Qtens[cnt][cnt2] = operator on site cnt-cnt2; at boundary cnt+1
      Qtensors[index] = new TensorQ*[index+1];
      for (int cnt2=0; cnt2<index+1; cnt2++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_q(L, index-cnt2) == MPIRANK ){
         #endif
            Qtensors[index][cnt2] = new TensorQ(index+1,denBK->gIrrep(index-cnt2),movingRight,denBK,Prob,index-cnt2);
         #ifdef CHEMPS2_MPI_COMPILATION
         } else { Qtensors[index][cnt2] = NULL; }
         #endif
      }
   
      //Xtensors : a certain process owns the Xtensors
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_x() == MPIRANK ){
      #endif
         Xtensors[index] = new TensorX(index+1,movingRight,denBK,Prob);
      #ifdef CHEMPS2_MPI_COMPILATION
      } else { Xtensors[index] = NULL; }
      #endif
      
      //Otensors : certain processes own certain excitations
      if (Exc_activated){
         for (int state=0; state<nStates-1; state++){
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
            #endif
            { Exc_Overlaps[state][index] = new TensorO(index+1,movingRight,Exc_BKs[state],denBK,Prob); }
         }
      }
   
   }

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_ALLOC ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

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

void CheMPS2::DMRG::OperatorsOnDisk(const int index, const bool movingRight, const bool store){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const int Nbound = movingRight ? index+1 : L-1-index;
   const int Cbound = movingRight ? L-1-index : index+1;
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   std::stringstream thefilename;
   //The PID is different for each MPI process
   thefilename << tempfolder << "/" << CheMPS2::DMRG_OPERATOR_storage_prefix << thePID << "_index_" << index << ".h5";

   //The hdf5 file
   const hid_t file_id = ( store ) ? H5Fcreate( thefilename.str().c_str(), H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT )
                                   : H5Fopen(   thefilename.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );

   //Ltensors : all processes own all Ltensors
   for (int cnt2=0; cnt2<Nbound ; cnt2++){
      std::stringstream sstream;
      sstream << "/Ltensor_" << cnt2 ;
      if ( store ){ MY_HDF5_WRITE(file_id, sstream.str(), Ltensors[index][cnt2]); }
      else {        MY_HDF5_READ( file_id, sstream.str(), Ltensors[index][cnt2]); }
   }

   //Two-operator tensors : certain processes own certain two-operator tensors
   for (int cnt2=0; cnt2<Nbound ; cnt2++){
      for (int cnt3=0; cnt3<Nbound-cnt2 ; cnt3++){

         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = movingRight ? index - cnt2 - cnt3 : index + 1 + cnt3;
         const int siteindex2 = movingRight ? index - cnt3        : index + 1 + cnt2 + cnt3;
         if (( cnt3 == 0 ) || ( MPIchemps2::owner_cdf(L, siteindex1, siteindex2) == MPIRANK ))
         #endif
         {
            std::stringstream sstream1;
            sstream1 << "/F0tensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream1.str(), F0tensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream1.str(), F0tensors[index][cnt2][cnt3]); }
            
            std::stringstream sstream2;
            sstream2 << "/F1tensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream2.str(), F1tensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream2.str(), F1tensors[index][cnt2][cnt3]); }
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         if (( cnt3 == 0 ) || ( MPIchemps2::owner_absigma(siteindex1, siteindex2) == MPIRANK ))
         #endif
         {
            std::stringstream sstream3;
            sstream3 << "/S0tensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream3.str(), S0tensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream3.str(), S0tensors[index][cnt2][cnt3]); }
            
            if (cnt2>0){
               std::stringstream sstream4;
               sstream4 << "/S1tensor_" << cnt2 << "_" << cnt3 ;
               if ( store ){ MY_HDF5_WRITE(file_id, sstream4.str(), S1tensors[index][cnt2][cnt3]); }
               else {        MY_HDF5_READ( file_id, sstream4.str(), S1tensors[index][cnt2][cnt3]); }
            }
         }
      }
   }
   
   //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
   for (int cnt2=0; cnt2<Cbound ; cnt2++){
      for (int cnt3=0; cnt3<Cbound-cnt2 ; cnt3++){

         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = movingRight ? index + 1 + cnt3        : index - cnt2 - cnt3;
         const int siteindex2 = movingRight ? index + 1 + cnt2 + cnt3 : index - cnt3;
         if ( MPIchemps2::owner_absigma(siteindex1, siteindex2) == MPIRANK )
         #endif
         {
            std::stringstream sstream1;
            sstream1 << "/Atensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream1.str(), Atensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream1.str(), Atensors[index][cnt2][cnt3]); }

            if (cnt2>0){
               std::stringstream sstream2;
               sstream2 << "/Btensor_" << cnt2 << "_" << cnt3 ;
               if ( store ){ MY_HDF5_WRITE(file_id, sstream2.str(), Btensors[index][cnt2][cnt3]); }
               else {        MY_HDF5_READ( file_id, sstream2.str(), Btensors[index][cnt2][cnt3]); }
            }
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf(L, siteindex1, siteindex2) == MPIRANK )
         #endif
         {
            std::stringstream sstream3;
            sstream3 << "/Ctensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream3.str(), Ctensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream3.str(), Ctensors[index][cnt2][cnt3]); }

            std::stringstream sstream4;
            sstream4 << "/Dtensor_" << cnt2 << "_" << cnt3 ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream4.str(), Dtensors[index][cnt2][cnt3]); }
            else {        MY_HDF5_READ( file_id, sstream4.str(), Dtensors[index][cnt2][cnt3]); }
         }
      }
   }

   //Qtensors : certain processes own certain Qtensors
   for (int cnt2=0; cnt2<Cbound; cnt2++){

      #ifdef CHEMPS2_MPI_COMPILATION
      const int siteindex = movingRight ? index + 1 + cnt2 : index - cnt2;
      if ( MPIchemps2::owner_q(L, siteindex) == MPIRANK )
      #endif
      {
         std::stringstream sstream;
         sstream << "/Qtensor_" << cnt2 ;
         if ( store ){ MY_HDF5_WRITE(file_id, sstream.str(), Qtensors[index][cnt2]); }
         else {        MY_HDF5_READ( file_id, sstream.str(), Qtensors[index][cnt2]); }
      }
   }

   //Xtensors : a certain process owns the Xtensors
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::owner_x() == MPIRANK )
   #endif
   {
      std::string sPathX = "/Xtensor" ;
      if ( store ){ MY_HDF5_WRITE(file_id, sPathX, Xtensors[index]); }
      else {        MY_HDF5_READ( file_id, sPathX, Xtensors[index]); }
   }

   //Otensors : certain processes own certain excitations
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
         #endif
         {
            std::stringstream sstream;
            sstream << "/Otensor_" << state ;
            if ( store ){ MY_HDF5_WRITE(file_id, sstream.str(), Exc_Overlaps[state][index]); }
            else {        MY_HDF5_READ( file_id, sstream.str(), Exc_Overlaps[state][index]); }
         }
      }
   }

   H5Fclose(file_id);

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_DISKIO ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::deleteTensors(const int index, const bool movingRight){

   struct timeval start, end;
   gettimeofday(&start, NULL);

   const int Nbound = movingRight ? index+1 : L-1-index;
   const int Cbound = movingRight ? L-1-index : index+1;
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   //Ltensors : all processes own all Ltensors
   for (int cnt2=0; cnt2<Nbound; cnt2++){ delete Ltensors[index][cnt2]; }
   delete [] Ltensors[index];
   
   //Two-operator tensors : certain processes own certain two-operator tensors
   for (int cnt2=0; cnt2<Nbound; cnt2++){
      for (int cnt3=0; cnt3<Nbound-cnt2; cnt3++){
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = movingRight ? index - cnt2 - cnt3 : index + 1 + cnt3;
         const int siteindex2 = movingRight ? index - cnt3        : index + 1 + cnt2 + cnt3;
         if (( cnt3 == 0 ) || ( MPIchemps2::owner_cdf(L, siteindex1, siteindex2) == MPIRANK ))
         #endif
         {
            delete F0tensors[index][cnt2][cnt3];
            delete F1tensors[index][cnt2][cnt3];
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         if (( cnt3 == 0 ) || ( MPIchemps2::owner_absigma(siteindex1, siteindex2) == MPIRANK ))
         #endif
         {
            delete S0tensors[index][cnt2][cnt3];
            if (cnt2>0){ delete S1tensors[index][cnt2][cnt3]; }
         }
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
   
   //Complementary two-operator tensors : certain processes own certain complementary two-operator tensors
   for (int cnt2=0; cnt2<Cbound; cnt2++){
      for (int cnt3=0; cnt3<Cbound-cnt2; cnt3++){
         #ifdef CHEMPS2_MPI_COMPILATION
         const int siteindex1 = movingRight ? index + 1 + cnt3        : index - cnt2 - cnt3;
         const int siteindex2 = movingRight ? index + 1 + cnt2 + cnt3 : index - cnt3;
         if ( MPIchemps2::owner_absigma(siteindex1, siteindex2) == MPIRANK )
         #endif
         {
            delete Atensors[index][cnt2][cnt3];
            if (cnt2>0){ delete Btensors[index][cnt2][cnt3]; }
         }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf(L, siteindex1, siteindex2) == MPIRANK )
         #endif
         {
            delete Ctensors[index][cnt2][cnt3];
            delete Dtensors[index][cnt2][cnt3];
         }
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
   
   //Qtensors : certain processes own certain Qtensors
   for (int cnt2=0; cnt2<Cbound; cnt2++){
      #ifdef CHEMPS2_MPI_COMPILATION
      const int siteindex = movingRight ? index + 1 + cnt2 : index - cnt2;
      if ( MPIchemps2::owner_q(L, siteindex) == MPIRANK )
      #endif
      { delete Qtensors[index][cnt2]; }
   }
   delete [] Qtensors[index];
   
   //Xtensors
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::owner_x() == MPIRANK )
   #endif
   { delete Xtensors[index]; }
   
   //Otensors
   if (Exc_activated){
      for (int state=0; state<nStates-1; state++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIRANK )
         #endif
         { delete Exc_Overlaps[state][index]; }
      }
   }

   gettimeofday(&end, NULL);
   timings[ CHEMPS2_TIME_TENS_FREE ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

}

void CheMPS2::DMRG::deleteStoredOperators(){

   std::stringstream temp;
   temp << "rm " << tempfolder << "/" << CheMPS2::DMRG_OPERATOR_storage_prefix << thePID << "*.h5";
   int info = system(temp.str().c_str());
   std::cout << "Info on DMRG::operators rm call to system: " << info << std::endl;

}

