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

#include "TensorQ.h"
#include "Lapack.h"
#include "Gsl.h"

CheMPS2::TensorQ::TensorQ(const int indexIn, const int IdiffIn, const bool movingRightIn, const SyBookkeeper * denBKIn, const Problem * ProbIn, const int siteIn) : TensorSwap(indexIn, IdiffIn, movingRightIn, denBKIn){

   Prob = ProbIn;
   site = siteIn;

}

CheMPS2::TensorQ::~TensorQ(){

}

void CheMPS2::TensorQ::ClearStorage(){ Clear(); }

void CheMPS2::TensorQ::AddTermSimple(TensorT * denT){

   if (( movingRight) && (denBK->gIrrep(denT->gIndex()) == Idiff)){ AddTermSimpleRight(denT); }
   if ((!movingRight) && (denBK->gIrrep(denT->gIndex()) == Idiff)){ AddTermSimpleLeft( denT); }

}

void CheMPS2::TensorQ::AddTermSimpleRight(TensorT * denT){

   const double mxElement = Prob->gMxElement(index-1,index-1,index-1,site);
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimRU = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimRD = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
      int dimL  = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID);
      if (dimL>0){
      
         double * BlockTup = denT->gStorage(sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         double * BlockTdo = denT->gStorage(sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
         
         int fase = ((((sectorTwoSD[ikappa]+1-sectorTwoS1[ikappa])/2)%2)!=0)?-1:1;
         double alpha = fase * mxElement * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
         double beta = 1.0; //add
         char trans = 'T';
         char notr = 'N';
         
         dgemm_(&trans,&notr,&dimRU,&dimRD,&dimL,&alpha,BlockTup,&dimL,BlockTdo,&dimL,&beta,storage+kappa2index[ikappa],&dimRU);
      
      }
   }

}

void CheMPS2::TensorQ::AddTermSimpleLeft(TensorT * denT){

   const double mxElement = Prob->gMxElement(site,index,index,index);
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID = Irreps::directProd(Idiff,sectorI1[ikappa]);
      int dimLU = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimLD = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
      int dimR  = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
      if (dimR>0){
      
         double * BlockTup = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID,               sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         
         int fase = ((((sectorTwoS1[ikappa]+1-sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
         double alpha = fase * mxElement * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
         double beta = 1.0; //add
         char trans = 'T';
         char notr = 'N';
         
         dgemm_(&notr,&trans,&dimLU,&dimLD,&dimR,&alpha,BlockTup,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
      
      }
   }

}

void CheMPS2::TensorQ::AddTermsL(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2){

   if (movingRight){ AddTermsLRight(Ltensors, denT, workmem, workmem2); }
   else{ AddTermsLLeft( Ltensors, denT, workmem, workmem2); }

}

void CheMPS2::TensorQ::AddTermsLRight(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2){

   bool OneToAdd = false;
   for (int loca=0; loca<index-1; loca++){
      if (Ltensors[index-2-loca]->gIdiff() == Idiff){ OneToAdd = true; }
   }
   
   if (OneToAdd){
      for (int ikappa=0; ikappa<nKappa; ikappa++){
   
         const int ID = Irreps::directProd(sectorI1[ikappa],Idiff);
         int dimRU = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimRD = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
      
         //case 1
         int dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=0; loca<index-1; loca++){
               if (Ltensors[index-2-loca]->gIdiff() == Idiff){
                  double * BlockL = Ltensors[index-2-loca]->gStorage(sectorN1[ikappa]-1,sectorTwoSD[ikappa],ID,sectorN1[ikappa],sectorTwoS1[ikappa],sectorI1[ikappa]);
                  double alpha = Prob->gMxElement(loca,index-1,index-1,site);
                  int inc = 1;
                  daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            int fase = ((((sectorTwoSD[ikappa]+1-sectorTwoS1[ikappa])/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
            char totrans = 'T';
            // factor * Tup^T * L^T --> mem2
            dgemm_(&totrans,&totrans,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLD,&beta,workmem2,&dimRU);
            
            alpha = 1.0;
            beta = 1.0; //add
            totrans = 'N';
            // mem2 * Tdo --> storage
            dgemm_(&totrans,&totrans,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimRU);
         
         }
         
         //case 2
         dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         //dimLD same as case1
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=0; loca<index-1; loca++){
               if (Ltensors[index-2-loca]->gIdiff() == Idiff){
                  double * BlockL = Ltensors[index-2-loca]->gStorage(sectorN1[ikappa]-2,sectorTwoS1[ikappa],sectorI1[ikappa],sectorN1[ikappa]-1,sectorTwoSD[ikappa],ID);
                  double alpha = 2*Prob->gMxElement(loca,index-1,site,index-1) - Prob->gMxElement(loca,index-1,index-1,site);
                  int inc = 1;
                  daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            double alpha = 1.0; //factor = 1 in this case
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]-1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
            char trans = 'T';
            char notr = 'N';
            // factor * Tup^T * L --> mem2
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLU,&beta,workmem2,&dimRU);
            
            beta = 1.0; //add
            // mem2 * Tdo --> storage
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimRU);
         
         }
         
         //case 3
         for (int TwoSLU=sectorTwoS1[ikappa]-1; TwoSLU<=sectorTwoS1[ikappa]+1; TwoSLU+=2){
            for (int TwoSLD=sectorTwoSD[ikappa]-1; TwoSLD<=sectorTwoSD[ikappa]+1; TwoSLD+=2){
               if ((TwoSLD>=0) && (TwoSLU>=0) && (abs(TwoSLD-TwoSLU)<2)){
                  const int ILU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index-1));
                  const int ILD = Irreps::directProd(ID,              denBK->gIrrep(index-1));
                  dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, TwoSLU, ILU);
                  dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa]  , TwoSLD, ILD);
                  if ((dimLU>0) && (dimLD>0)){
                     int fase = ((((sectorTwoS1[ikappa]+TwoSLD)/2)%2)!=0)?-1:1;
                     double factor = fase * sqrt((TwoSLD+1)*(sectorTwoS1[ikappa]+1.0)) * gsl_sf_coupling_6j(sectorTwoS1[ikappa], sectorTwoSD[ikappa], 1, TwoSLD, TwoSLU, 1);
                  
                     int dimLUxLD = dimLU * dimLD;
                     for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
                     for (int loca=0; loca<index-1; loca++){
                        if (Ltensors[index-2-loca]->gIdiff() == Idiff){
                           double * BlockL = Ltensors[index-2-loca]->gStorage(sectorN1[ikappa]-1, TwoSLU, ILU, sectorN1[ikappa], TwoSLD, ILD);
                           double alpha = factor * Prob->gMxElement(loca,index-1,site,index-1);
                           if (TwoSLD==sectorTwoS1[ikappa]){ alpha += Prob->gMxElement(loca,index-1,index-1,site); }
                           int inc = 1;
                           daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
                        }
                     }

                     double alpha = 1.0;
                     double beta = 0.0; //set
         
                     double * BlockTup = denT->gStorage(sectorN1[ikappa]-1, TwoSLU, ILU, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
                     double * BlockTdo = denT->gStorage(sectorN1[ikappa]  , TwoSLD, ILD, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
                     char trans = 'T';
                     char notr = 'N';
                     // Tup^T * mem --> mem2
                     dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLU,&beta,workmem2,&dimRU);
            
                     beta = 1.0; //add
                     // mem2 * Tdo --> storage
                     dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimRU);
         
                  }
               }
            }
         }
   
      }
   }

}

void CheMPS2::TensorQ::AddTermsLLeft(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2){

   bool OneToAdd = false;
   for (int loca=index+1; loca<Prob->gL(); loca++){
      if (Ltensors[loca-index-1]->gIdiff() == Idiff){ OneToAdd = true; }
   }
   
   if (OneToAdd){
      for (int ikappa=0; ikappa<nKappa; ikappa++){
   
         const int ID = Irreps::directProd(sectorI1[ikappa],Idiff);
         int dimLU = denBK->gCurrentDim(index,   sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimLD = denBK->gCurrentDim(index,   sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
      
         //case 1
         int dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=index+1; loca<Prob->gL(); loca++){
               if (Ltensors[loca-index-1]->gIdiff() == Idiff){
                  double * BlockL = Ltensors[loca-index-1]->gStorage(sectorN1[ikappa]+1,sectorTwoSD[ikappa],ID,sectorN1[ikappa]+2,sectorTwoS1[ikappa],sectorI1[ikappa]);
                  double alpha = Prob->gMxElement(site,index,index,loca);
                  int inc = 1;
                  daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            int fase = ((((sectorTwoS1[ikappa]+1-sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0));
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
            char trans = 'T';
            char notr = 'N';
            // factor * Tup * L^T --> mem2
            dgemm_(&notr,&trans,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRD,&beta,workmem2,&dimLU);
            
            alpha = 1.0;
            beta = 1.0; //add
            // mem2 * Tdo^T --> storage
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimLU);
         
         }
         
         //case 2
         dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID);
         //dimRU same as case1
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=index+1; loca<Prob->gL(); loca++){
               if (Ltensors[loca-index-1]->gIdiff() == Idiff){
                  double * BlockL = Ltensors[loca-index-1]->gStorage(sectorN1[ikappa]+2,sectorTwoS1[ikappa],sectorI1[ikappa],sectorN1[ikappa]+3,sectorTwoSD[ikappa],ID);
                  double alpha = 2*Prob->gMxElement(site,index,loca,index) - Prob->gMxElement(site,index,index,loca);
                  int inc = 1;
                  daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            double alpha = 1.0; //factor = 1 in this case
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID);
            
            char notr = 'N';
            // factor * Tup * L --> mem2
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRU,&beta,workmem2,&dimLU);
            
            beta = 1.0; //add
            // mem2 * Tdo^T --> storage
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimLU);
         
         }
         
         //case 3
         for (int TwoSRU=sectorTwoS1[ikappa]-1; TwoSRU<=sectorTwoS1[ikappa]+1; TwoSRU+=2){
            for (int TwoSRD=sectorTwoSD[ikappa]-1; TwoSRD<=sectorTwoSD[ikappa]+1; TwoSRD+=2){
               if ((TwoSRD>=0) && (TwoSRU>=0) && (abs(TwoSRD-TwoSRU)<2)){
                  const int IRU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index));
                  const int IRD = Irreps::directProd(ID,              denBK->gIrrep(index));
                  dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa]+1, TwoSRU, IRU);
                  dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, TwoSRD, IRD);
                  if ((dimRU>0) && (dimRD>0)){
                     int fase = ((((sectorTwoSD[ikappa]+TwoSRU)/2)%2)!=0)?-1:1;
                     double factor1 = fase * sqrt((TwoSRU+1.0)/(sectorTwoSD[ikappa]+1.0)) * (TwoSRD+1) * gsl_sf_coupling_6j(sectorTwoS1[ikappa], sectorTwoSD[ikappa], 1, TwoSRD, TwoSRU, 1);
                     double factor2 = (TwoSRD+1.0)/(sectorTwoSD[ikappa]+1.0);
                  
                     int dimRUxRD = dimRU * dimRD;
                     for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
                     for (int loca=index+1; loca<Prob->gL(); loca++){
                        if (Ltensors[loca-index-1]->gIdiff() == Idiff){
                           double * BlockL = Ltensors[loca-index-1]->gStorage(sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+2, TwoSRD, IRD);
                           double alpha = factor1 * Prob->gMxElement(site,index,loca,index);
                           if (TwoSRU==sectorTwoSD[ikappa]){ alpha += factor2 * Prob->gMxElement(site,index,index,loca); }
                           int inc = 1;
                           daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
                        }
                     }

                     double alpha = 1.0;
                     double beta = 0.0; //set
         
                     double * BlockTup = denT->gStorage(sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+1, TwoSRU, IRU);
                     double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID,               sectorN1[ikappa]+2, TwoSRD, IRD);
            
                     char notr = 'N';
                     // Tup * mem --> mem2
                     dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRU,&beta,workmem2,&dimLU);
            
                     beta = 1.0; //add
                     // mem2 * Tdo^T --> storage
                     char trans = 'T';
                     dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimLU);
         
                  }
               }
            }
         }
   
      }
   }

}

void CheMPS2::TensorQ::AddTermsAB(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2){

   if (movingRight){ AddTermsABRight(denA, denB, denT, workmem, workmem2); }
   else{ AddTermsABLeft( denA, denB, denT, workmem, workmem2); }

}

void CheMPS2::TensorQ::AddTermsABRight(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimRU = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimRD = denBK->gCurrentDim(index, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
   
      //case 1
      const int ILU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index-1));
      for (int TwoSLU=sectorTwoS1[ikappa]-1; TwoSLU<=sectorTwoS1[ikappa]+1; TwoSLU+=2){
         int dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, TwoSLU,              ILU);
         int dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
         if ((dimLU>0) && (dimLD>0)){
         
            int fase = ((((TwoSLU + sectorTwoSD[ikappa] + 2)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(sectorTwoS1[ikappa]+1)) * gsl_sf_coupling_6j(1,2,1,sectorTwoSD[ikappa],sectorTwoS1[ikappa],TwoSLU);
            
            double alpha;
            double * mem;
            
            if (TwoSLU == sectorTwoSD[ikappa]){
            
               fase = ((((sectorTwoSD[ikappa]+1-sectorTwoS1[ikappa])/2)%2)!=0)?-1:1;
               double factorA = fase * sqrt( 0.5 * (sectorTwoS1[ikappa]+1.0) / (sectorTwoSD[ikappa]+1.0) );
            
               double * BlockA = denA->gStorage( sectorN1[ikappa]-1,TwoSLU,ILU,sectorN1[ikappa]+1,sectorTwoSD[ikappa],ID );
               double * BlockB = denB->gStorage( sectorN1[ikappa]-1,TwoSLU,ILU,sectorN1[ikappa]+1,sectorTwoSD[ikappa],ID );
               
               mem = workmem;
               for (int cnt=0; cnt<dimLU*dimLD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage(sectorN1[ikappa]-1,TwoSLU,ILU,sectorN1[ikappa]+1,sectorTwoSD[ikappa],ID);
            
            }
            
            double * BlockTup = denT->gStorage(sectorN1[ikappa]-1, TwoSLU,              ILU, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID,  sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
            char trans = 'T';
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,mem,&dimLU,&beta,workmem2,&dimRU);
            
            alpha = 1.0;
            beta = 1.0; //add
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimRU);
            
         }
      }

      //case 2
      const int ILD = Irreps::directProd(ID,denBK->gIrrep(index-1));
      for (int TwoSLD=sectorTwoSD[ikappa]-1; TwoSLD<=sectorTwoSD[ikappa]+1; TwoSLD+=2){
         int dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa],   TwoSLD,              ILD);
         if ((dimLU>0) && (dimLD>0)){
         
            int fase = ((((sectorTwoS1[ikappa] + sectorTwoSD[ikappa] + 1)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(TwoSLD+1)) * gsl_sf_coupling_6j(1,2,1,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSLD);
            
            double alpha;
            double * mem;
            
            if (TwoSLD == sectorTwoS1[ikappa]){

               double factorA = - sqrt(0.5);
            
               double * BlockA = denA->gStorage( sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], TwoSLD, ILD );
               double * BlockB = denB->gStorage( sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], TwoSLD, ILD );
               
               mem = workmem;
               for (int cnt=0; cnt<dimLU*dimLD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sectorN1[ikappa]-2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], TwoSLD, ILD );
            
            }
            
            double * BlockTup = denT->gStorage(sectorN1[ikappa]-2,sectorTwoS1[ikappa],sectorI1[ikappa],sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]  ,TwoSLD,             ILD,             sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
            
            char trans = 'T';
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,mem,&dimLU,&beta,workmem2,&dimRU);
            
            alpha = 1.0;
            beta = 1.0; //add
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimRU);
            
         }
      }
   }

}

void CheMPS2::TensorQ::AddTermsABLeft(TensorA * denA, TensorB * denB, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID  = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimLU = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimLD = denBK->gCurrentDim(index, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID);
   
      //case 1
      const int IRD = Irreps::directProd(ID, denBK->gIrrep(index));
      for (int TwoSRD=sectorTwoSD[ikappa]-1; TwoSRD<=sectorTwoSD[ikappa]+1; TwoSRD+=2){
         int dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, TwoSRD,              IRD);
         if ((dimRU>0) && (dimRD>0)){
         
            int fase = ((((TwoSRD + sectorTwoS1[ikappa] + 2)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0/(sectorTwoSD[ikappa]+1.0)) * (TwoSRD+1) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],TwoSRD,sectorTwoSD[ikappa]);
            
            double alpha;
            double * mem;
            
            if (TwoSRD == sectorTwoS1[ikappa]){
            
               fase = ((((sectorTwoS1[ikappa]+1-sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
               double factorA = fase * sqrt( 0.5 * (sectorTwoS1[ikappa]+1.0) / (sectorTwoSD[ikappa]+1.0) );
            
               double * BlockA = denA->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, TwoSRD, IRD );
               double * BlockB = denB->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, TwoSRD, IRD );
               
               mem = workmem;
               for (int cnt=0; cnt<dimRU*dimRD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, TwoSRD, IRD );
            
            }
            
            double * BlockTup = denT->gStorage(sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa],sectorN1[ikappa],  sectorTwoS1[ikappa], sectorI1[ikappa]);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1,sectorTwoSD[ikappa],ID,              sectorN1[ikappa]+2,TwoSRD,              IRD);
            
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,mem,&dimRU,&beta,workmem2,&dimLU);
            
            alpha = 1.0;
            beta = 1.0; //add
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
            
         }
      }

      //case 2
      const int IRU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index));
      for (int TwoSRU=sectorTwoS1[ikappa]-1; TwoSRU<=sectorTwoS1[ikappa]+1; TwoSRU+=2){
         int dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa]+1, TwoSRU,              IRU);
         int dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID);
         if ((dimRU>0) && (dimRD>0)){
         
            int fase = ((((sectorTwoS1[ikappa] + sectorTwoSD[ikappa] + 1)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(TwoSRU+1)) * gsl_sf_coupling_6j(1,1,2,TwoSRU,sectorTwoSD[ikappa],sectorTwoS1[ikappa]);
            
            double alpha;
            double * mem;
            
            if (TwoSRU == sectorTwoSD[ikappa]){

               double factorA = - sqrt(0.5);
            
               double * BlockA = denA->gStorage( sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID );
               double * BlockB = denB->gStorage( sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID );
               
               mem = workmem;
               for (int cnt=0; cnt<dimRU*dimRD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID );
            
            }
            
            double * BlockTup = denT->gStorage(sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+1, TwoSRU, IRU);
            double * BlockTdo = denT->gStorage(sectorN1[ikappa]+1, sectorTwoSD[ikappa], ID, sectorN1[ikappa]+3, sectorTwoSD[ikappa], ID);
            
            char notr = 'N';
            double beta = 0.0; //set
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,mem,&dimRU,&beta,workmem2,&dimLU);
            
            alpha = 1.0;
            beta = 1.0; //add
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
            
         }
      }
   }

}

void CheMPS2::TensorQ::AddTermsCD(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2){

   if (movingRight){ AddTermsCDRight(denC, denD, denT, workmem, workmem2); }
   else{ AddTermsCDLeft(denC, denD, denT, workmem, workmem2); }
   
}

void CheMPS2::TensorQ::AddTermsCDRight(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IRD  = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimRU = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimRD = denBK->gCurrentDim(index, sectorN1[ikappa]+1, sectorTwoSD[ikappa], IRD);
   
      //case 1
      const int ILD = Irreps::directProd(IRD,denBK->gIrrep(index-1));
      for (int TwoSLD=sectorTwoSD[ikappa]-1; TwoSLD<=sectorTwoSD[ikappa]+1; TwoSLD+=2){
         int dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa], TwoSLD,              ILD);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            
            //first set to D
            int fase = ((((sectorTwoS1[ikappa]+sectorTwoSD[ikappa]+1)/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(TwoSLD+1)) * gsl_sf_coupling_6j(1,2,1,sectorTwoS1[ikappa],sectorTwoSD[ikappa],TwoSLD);
            double * block = denD->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], TwoSLD, ILD );
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSLD==sectorTwoS1[ikappa]){
               factor = sqrt(0.5);
               block = denC->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], TwoSLD, ILD );
               int inc = 1;
               daxpy_(&dimLUxLD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa] );
            double * BlockTdo = denT->gStorage( sectorN1[ikappa], TwoSLD, ILD, sectorN1[ikappa]+1, sectorTwoSD[ikappa], IRD );
            
            char trans = 'T';
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLU,&beta,workmem2,&dimRU);
            beta = 1.0;
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimRU);
            
         }
      }

      //case 2
      const int ILU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index-1));
      for (int TwoSLU=sectorTwoS1[ikappa]-1; TwoSLU<=sectorTwoS1[ikappa]+1; TwoSLU+=2){
         int dimLU = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, TwoSLU,              ILU);
         int dimLD = denBK->gCurrentDim(index-1, sectorN1[ikappa]-1, sectorTwoSD[ikappa], IRD);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            
            //first set to D
            int fase = ((((TwoSLU + sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(sectorTwoS1[ikappa]+1)) * gsl_sf_coupling_6j(1,2,1,sectorTwoSD[ikappa],sectorTwoS1[ikappa],TwoSLU);
            double * block = denD->gStorage( sectorN1[ikappa]-1, TwoSLU, ILU, sectorN1[ikappa]-1, sectorTwoSD[ikappa], IRD );
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSLU==sectorTwoSD[ikappa]){
               fase = ((((sectorTwoSD[ikappa]+1-sectorTwoS1[ikappa])/2)%2)!=0)?-1:1;
               factor = fase * sqrt(0.5 * ( sectorTwoS1[ikappa]+1.0 ) / ( sectorTwoSD[ikappa] + 1.0 ) );
               block = denC->gStorage( sectorN1[ikappa]-1, TwoSLU, ILU, sectorN1[ikappa]-1, sectorTwoSD[ikappa], IRD );
               int inc = 1;
               daxpy_(&dimLUxLD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sectorN1[ikappa]-1, TwoSLU, ILU, sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa] );
            double * BlockTdo = denT->gStorage( sectorN1[ikappa]-1, sectorTwoSD[ikappa], IRD, sectorN1[ikappa]+1, sectorTwoSD[ikappa], IRD );
            
            char trans = 'T';
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLU,&beta,workmem2,&dimRU);
            beta = 1.0;
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimRU);
            
         }
      }
   }

}

void CheMPS2::TensorQ::AddTermsCDLeft(TensorC * denC, TensorD * denD, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ILD  = Irreps::directProd(sectorI1[ikappa],Idiff);
      int dimLU = denBK->gCurrentDim(index, sectorN1[ikappa],   sectorTwoS1[ikappa], sectorI1[ikappa]);
      int dimLD = denBK->gCurrentDim(index, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD);
   
      //case 1
      const int IRU = Irreps::directProd(sectorI1[ikappa],denBK->gIrrep(index));
      for (int TwoSRU=sectorTwoS1[ikappa]-1; TwoSRU<=sectorTwoS1[ikappa]+1; TwoSRU+=2){
         int dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa]+1, TwoSRU, IRU);
         int dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            
            //first set to D
            int fase = ((((sectorTwoS1[ikappa]+TwoSRU+3)/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0/(sectorTwoSD[ikappa]+1.0)) * (TwoSRU+1) * gsl_sf_coupling_6j(1,1,2,TwoSRU,sectorTwoSD[ikappa],sectorTwoS1[ikappa]);
            double * block = denD->gStorage( sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD );
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSRU==sectorTwoSD[ikappa]){
               factor = sqrt(0.5);
               block = denC->gStorage( sectorN1[ikappa]+1, TwoSRU, IRU, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD );
               int inc = 1;
               daxpy_(&dimRUxRD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sectorN1[ikappa], sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+1, TwoSRU, IRU );
            double * BlockTdo = denT->gStorage( sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD, sectorN1[ikappa]+1, sectorTwoSD[ikappa], ILD );
            
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRU,&beta,workmem2,&dimLU);
            beta = 1.0;
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
            
         }
      }

      //case 2
      const int IRD = Irreps::directProd(ILD,denBK->gIrrep(index));
      for (int TwoSRD=sectorTwoSD[ikappa]-1; TwoSRD<=sectorTwoSD[ikappa]+1; TwoSRD+=2){
         int dimRU = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa]);
         int dimRD = denBK->gCurrentDim(index+1, sectorN1[ikappa]+2, TwoSRD,              IRD);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            
            //first set to D
            int fase = (((TwoSRD+1)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(TwoSRD+1.0)*(sectorTwoS1[ikappa]+1.0)/(sectorTwoSD[ikappa]+1.0)) * gsl_sf_coupling_6j(1,1,2,sectorTwoS1[ikappa],TwoSRD,sectorTwoSD[ikappa]);
            double * block = denD->gStorage( sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, TwoSRD, IRD );
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSRD==sectorTwoS1[ikappa]){
               fase = ((((sectorTwoS1[ikappa]+1-sectorTwoSD[ikappa])/2)%2)!=0)?-1:1;
               factor = fase * sqrt(0.5 * ( sectorTwoS1[ikappa]+1.0 ) / ( sectorTwoSD[ikappa] + 1.0 ) );
               block = denC->gStorage( sectorN1[ikappa]+2, sectorTwoS1[ikappa], sectorI1[ikappa], sectorN1[ikappa]+2, TwoSRD, IRD );
               int inc = 1;
               daxpy_(&dimRUxRD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sectorN1[ikappa],  sectorTwoS1[ikappa],sectorI1[ikappa],sectorN1[ikappa]+2,sectorTwoS1[ikappa],sectorI1[ikappa] );
            double * BlockTdo = denT->gStorage( sectorN1[ikappa]+1,sectorTwoSD[ikappa],ILD,             sectorN1[ikappa]+2,TwoSRD,             IRD);
            
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRU,&beta,workmem2,&dimLU);
            beta = 1.0;
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
            
         }
      }
   }

}



