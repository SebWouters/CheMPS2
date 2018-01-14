/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

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
#include "Wigner.h"

CheMPS2::TensorQ::TensorQ(const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * denBK, const Problem * Prob, const int site) :
TensorOperator(boundary_index,
               1, //two_j
               1, //n_elec
               Idiff,
               moving_right,
               true, //prime_last
               true, //jw_phase (three 2nd quantized operators)
               denBK,
               denBK){

   this->Prob = Prob;
   this->site = site;

}

CheMPS2::TensorQ::~TensorQ(){ }

void CheMPS2::TensorQ::AddTermSimple(TensorT * denT){

   if (( moving_right) && (bk_up->gIrrep(denT->gIndex()) == n_irrep)){ AddTermSimpleRight(denT); }
   if ((!moving_right) && (bk_up->gIrrep(denT->gIndex()) == n_irrep)){ AddTermSimpleLeft( denT); }

}

void CheMPS2::TensorQ::AddTermSimpleRight(TensorT * denT){

   const double mxElement = Prob->gMxElement(index-1,index-1,index-1,site);
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID = Irreps::directProd( n_irrep, sector_irrep_up[ikappa] );
      int dimRU = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimRD = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
      int dimL  = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID);
      if (dimL>0){
      
         double * BlockTup = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
         
         int fase = ((((sector_spin_down[ikappa]+1-sector_spin_up[ikappa])/2)%2)!=0)?-1:1;
         double alpha = fase * mxElement * sqrt((sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1.0));
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
      const int ID = Irreps::directProd( n_irrep, sector_irrep_up[ikappa] );
      int dimLU = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimLD = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
      int dimR  = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      if (dimR>0){
      
         double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         
         int fase = ((((sector_spin_up[ikappa]+1-sector_spin_down[ikappa])/2)%2)!=0)?-1:1;
         double alpha = fase * mxElement * sqrt((sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1.0));
         double beta = 1.0; //add
         char trans = 'T';
         char notr = 'N';
         
         dgemm_(&notr,&trans,&dimLU,&dimLD,&dimR,&alpha,BlockTup,&dimLU,BlockTdo,&dimLD,&beta,storage+kappa2index[ikappa],&dimLU);
      
      }
   }

}

void CheMPS2::TensorQ::AddTermsL(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2){

   if (moving_right){ AddTermsLRight(Ltensors, denT, workmem, workmem2); }
   else{ AddTermsLLeft( Ltensors, denT, workmem, workmem2); }

}

void CheMPS2::TensorQ::AddTermsLRight(TensorL ** Ltensors, TensorT * denT, double * workmem, double * workmem2){

   bool OneToAdd = false;
   for (int loca=0; loca<index-1; loca++){
      if (Ltensors[index-2-loca]->get_irrep() == n_irrep){ OneToAdd = true; }
   }
   
   if (OneToAdd){
      for (int ikappa=0; ikappa<nKappa; ikappa++){
   
         const int ID = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
         int dimRU = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimRD = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
      
         //case 1
         int dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=0; loca<index-1; loca++){
               if (Ltensors[index-2-loca]->get_irrep() == n_irrep){
                  double * BlockL = Ltensors[index-2-loca]->gStorage(sector_nelec_up[ikappa]-1,sector_spin_down[ikappa],ID,sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa]);
                  double alpha = Prob->gMxElement(loca,site,index-1,index-1);
                  int inc = 1;
                  daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            int fase = ((((sector_spin_down[ikappa]+1-sector_spin_up[ikappa])/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1.0));
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
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
         dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         //dimLD same as case1
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=0; loca<index-1; loca++){
               if (Ltensors[index-2-loca]->get_irrep() == n_irrep){
                  double * BlockL = Ltensors[index-2-loca]->gStorage(sector_nelec_up[ikappa]-2,sector_spin_up[ikappa],sector_irrep_up[ikappa],sector_nelec_up[ikappa]-1,sector_spin_down[ikappa],ID);
                  double alpha = 2*Prob->gMxElement(loca,index-1,site,index-1) - Prob->gMxElement(loca,index-1,index-1,site);
                  int inc = 1;
                  daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            double alpha = 1.0; //factor = 1 in this case
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa]-2, sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
            char trans = 'T';
            char notr = 'N';
            // factor * Tup^T * L --> mem2
            dgemm_(&trans,&notr,&dimRU,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,workmem,&dimLU,&beta,workmem2,&dimRU);
            
            beta = 1.0; //add
            // mem2 * Tdo --> storage
            dgemm_(&notr,&notr,&dimRU,&dimRD,&dimLD,&alpha,workmem2,&dimRU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimRU);
         
         }
         
         //case 3
         for (int TwoSLU=sector_spin_up[ikappa]-1; TwoSLU<=sector_spin_up[ikappa]+1; TwoSLU+=2){
            for (int TwoSLD=sector_spin_down[ikappa]-1; TwoSLD<=sector_spin_down[ikappa]+1; TwoSLD+=2){
               if ((TwoSLD>=0) && (TwoSLU>=0) && (abs(TwoSLD-TwoSLU)<2)){
                  const int ILU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index-1));
                  const int ILD = Irreps::directProd(ID,              bk_up->gIrrep(index-1));
                  dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, TwoSLU, ILU);
                  dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]  , TwoSLD, ILD);
                  if ((dimLU>0) && (dimLD>0)){
                     int fase = ((((sector_spin_up[ikappa]+TwoSLD)/2)%2)!=0)?-1:1;
                     double factor = fase * sqrt((TwoSLD+1)*(sector_spin_up[ikappa]+1.0))
                                   * Wigner::wigner6j( sector_spin_up[ikappa], sector_spin_down[ikappa], 1, TwoSLD, TwoSLU, 1 );
                  
                     int dimLUxLD = dimLU * dimLD;
                     for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = 0.0; }
         
                     for (int loca=0; loca<index-1; loca++){
                        if (Ltensors[index-2-loca]->get_irrep() == n_irrep){
                           double * BlockL = Ltensors[index-2-loca]->gStorage(sector_nelec_up[ikappa]-1, TwoSLU, ILU, sector_nelec_up[ikappa], TwoSLD, ILD);
                           double alpha = factor * Prob->gMxElement(loca,index-1,site,index-1);
                           if (TwoSLD==sector_spin_up[ikappa]){ alpha += Prob->gMxElement(loca,index-1,index-1,site); }
                           int inc = 1;
                           daxpy_(&dimLUxLD, &alpha, BlockL, &inc, workmem, &inc);
                        }
                     }

                     double alpha = 1.0;
                     double beta = 0.0; //set
         
                     double * BlockTup = denT->gStorage(sector_nelec_up[ikappa]-1, TwoSLU, ILU, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
                     double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]  , TwoSLD, ILD, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
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
      if (Ltensors[loca-index-1]->get_irrep() == n_irrep){ OneToAdd = true; }
   }
   
   if (OneToAdd){
      for (int ikappa=0; ikappa<nKappa; ikappa++){
   
         const int ID = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
         int dimLU = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimLD = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
      
         //case 1
         int dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=index+1; loca<Prob->gL(); loca++){
               if (Ltensors[loca-index-1]->get_irrep() == n_irrep){
                  double * BlockL = Ltensors[loca-index-1]->gStorage(sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ID,sector_nelec_up[ikappa]+2,sector_spin_up[ikappa],sector_irrep_up[ikappa]);
                  double alpha = Prob->gMxElement(site,loca,index,index);
                  int inc = 1;
                  daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            int fase = ((((sector_spin_up[ikappa]+1-sector_spin_down[ikappa])/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1.0));
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
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
         dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID);
         //dimRU same as case1
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
            for (int loca=index+1; loca<Prob->gL(); loca++){
               if (Ltensors[loca-index-1]->get_irrep() == n_irrep){
                  double * BlockL = Ltensors[loca-index-1]->gStorage(sector_nelec_up[ikappa]+2,sector_spin_up[ikappa],sector_irrep_up[ikappa],sector_nelec_up[ikappa]+3,sector_spin_down[ikappa],ID);
                  double alpha = 2*Prob->gMxElement(site,index,loca,index) - Prob->gMxElement(site,index,index,loca);
                  int inc = 1;
                  daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
               }
            }

            double alpha = 1.0; //factor = 1 in this case
            double beta = 0.0; //set
         
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID);
            
            char notr = 'N';
            // factor * Tup * L --> mem2
            dgemm_(&notr,&notr,&dimLU,&dimRD,&dimRU,&alpha,BlockTup,&dimLU,workmem,&dimRU,&beta,workmem2,&dimLU);
            
            beta = 1.0; //add
            // mem2 * Tdo^T --> storage
            char trans = 'T';
            dgemm_(&notr,&trans,&dimLU,&dimLD,&dimRD,&alpha,workmem2,&dimLU, BlockTdo, &dimLD, &beta, storage+kappa2index[ikappa], &dimLU);
         
         }
         
         //case 3
         for (int TwoSRU=sector_spin_up[ikappa]-1; TwoSRU<=sector_spin_up[ikappa]+1; TwoSRU+=2){
            for (int TwoSRD=sector_spin_down[ikappa]-1; TwoSRD<=sector_spin_down[ikappa]+1; TwoSRD+=2){
               if ((TwoSRD>=0) && (TwoSRU>=0) && (abs(TwoSRD-TwoSRU)<2)){
                  const int IRU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index));
                  const int IRD = Irreps::directProd(ID,              bk_up->gIrrep(index));
                  dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+1, TwoSRU, IRU);
                  dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, TwoSRD, IRD);
                  if ((dimRU>0) && (dimRD>0)){
                     int fase = ((((sector_spin_down[ikappa]+TwoSRU)/2)%2)!=0)?-1:1;
                     double factor1 = fase * sqrt((TwoSRU+1.0)/(sector_spin_down[ikappa]+1.0)) * (TwoSRD+1)
                                    * Wigner::wigner6j( sector_spin_up[ikappa], sector_spin_down[ikappa], 1, TwoSRD, TwoSRU, 1 );
                     double factor2 = (TwoSRD+1.0)/(sector_spin_down[ikappa]+1.0);
                  
                     int dimRUxRD = dimRU * dimRD;
                     for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = 0.0; }
         
                     for (int loca=index+1; loca<Prob->gL(); loca++){
                        if (Ltensors[loca-index-1]->get_irrep() == n_irrep){
                           double * BlockL = Ltensors[loca-index-1]->gStorage(sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+2, TwoSRD, IRD);
                           double alpha = factor1 * Prob->gMxElement(site,index,loca,index);
                           if (TwoSRU==sector_spin_down[ikappa]){ alpha += factor2 * Prob->gMxElement(site,index,index,loca); }
                           int inc = 1;
                           daxpy_(&dimRUxRD, &alpha, BlockL, &inc, workmem, &inc);
                        }
                     }

                     double alpha = 1.0;
                     double beta = 0.0; //set
         
                     double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+1, TwoSRU, IRU);
                     double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+2, TwoSRD, IRD);
            
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

void CheMPS2::TensorQ::AddTermsAB(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2){

   if (moving_right){ AddTermsABRight(denA, denB, denT, workmem, workmem2); }
   else{ AddTermsABLeft( denA, denB, denT, workmem, workmem2); }

}

void CheMPS2::TensorQ::AddTermsABRight(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
      int dimRU = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimRD = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
   
      //case 1
      const int ILU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index-1));
      for (int TwoSLU=sector_spin_up[ikappa]-1; TwoSLU<=sector_spin_up[ikappa]+1; TwoSLU+=2){
         int dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, TwoSLU,                 ILU);
         int dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
         if ((dimLU>0) && (dimLD>0)){
         
            int fase = ((((TwoSLU + sector_spin_down[ikappa] + 2)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(sector_spin_up[ikappa]+1))
                           * Wigner::wigner6j( 1, 2, 1, sector_spin_down[ikappa], sector_spin_up[ikappa], TwoSLU );
            
            double alpha;
            double * mem;
            
            if (TwoSLU == sector_spin_down[ikappa]){
            
               fase = ((((sector_spin_down[ikappa]+1-sector_spin_up[ikappa])/2)%2)!=0)?-1:1;
               double factorA = fase * sqrt( 0.5 * (sector_spin_up[ikappa]+1.0) / (sector_spin_down[ikappa]+1.0) );
            
               double * BlockA = denA->gStorage( sector_nelec_up[ikappa]-1,TwoSLU,ILU,sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ID );
               double * BlockB = denB->gStorage( sector_nelec_up[ikappa]-1,TwoSLU,ILU,sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ID );
               
               mem = workmem;
               for (int cnt=0; cnt<dimLU*dimLD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage(sector_nelec_up[ikappa]-1,TwoSLU,ILU,sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ID);
            
            }
            
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa]-1, TwoSLU,                 ILU, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,  sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
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
      const int ILD = Irreps::directProd(ID,bk_up->gIrrep(index-1));
      for (int TwoSLD=sector_spin_down[ikappa]-1; TwoSLD<=sector_spin_down[ikappa]+1; TwoSLD+=2){
         int dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         int dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa],   TwoSLD,              ILD);
         if ((dimLU>0) && (dimLD>0)){
         
            int fase = ((((sector_spin_up[ikappa] + sector_spin_down[ikappa] + 1)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(TwoSLD+1)) * Wigner::wigner6j( 1, 2, 1, sector_spin_up[ikappa], sector_spin_down[ikappa], TwoSLD );
            
            double alpha;
            double * mem;
            
            if (TwoSLD == sector_spin_up[ikappa]){

               double factorA = - sqrt(0.5);
            
               double * BlockA = denA->gStorage( sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], TwoSLD, ILD );
               double * BlockB = denB->gStorage( sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], TwoSLD, ILD );
               
               mem = workmem;
               for (int cnt=0; cnt<dimLU*dimLD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], TwoSLD, ILD );
            
            }
            
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa]-2,sector_spin_up[ikappa],sector_irrep_up[ikappa],sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]  ,TwoSLD,             ILD,             sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
            
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

void CheMPS2::TensorQ::AddTermsABLeft(TensorOperator * denA, TensorOperator * denB, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ID  = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
      int dimLU = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimLD = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID);
   
      //case 1
      const int IRD = Irreps::directProd(ID, bk_up->gIrrep(index));
      for (int TwoSRD=sector_spin_down[ikappa]-1; TwoSRD<=sector_spin_down[ikappa]+1; TwoSRD+=2){
         int dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         int dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, TwoSRD,              IRD);
         if ((dimRU>0) && (dimRD>0)){
         
            int fase = ((((TwoSRD + sector_spin_up[ikappa] + 2)/2)%2)!=0)?-1:1;
            const double factorB = fase * sqrt(3.0/(sector_spin_down[ikappa]+1.0)) * (TwoSRD+1)
                                 * Wigner::wigner6j( 1, 1, 2, sector_spin_up[ikappa], TwoSRD, sector_spin_down[ikappa] );
            
            double alpha;
            double * mem;
            
            if (TwoSRD == sector_spin_up[ikappa]){
            
               fase = ((((sector_spin_up[ikappa]+1-sector_spin_down[ikappa])/2)%2)!=0)?-1:1;
               double factorA = fase * sqrt( 0.5 * (sector_spin_up[ikappa]+1.0) / (sector_spin_down[ikappa]+1.0) );
            
               double * BlockA = denA->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, TwoSRD, IRD );
               double * BlockB = denB->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, TwoSRD, IRD );
               
               mem = workmem;
               for (int cnt=0; cnt<dimRU*dimRD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, TwoSRD, IRD );
            
            }
            
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],  sector_spin_up[ikappa],   sector_irrep_up[ikappa],sector_nelec_up[ikappa],  sector_spin_up[ikappa], sector_irrep_up[ikappa]);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ID,              sector_nelec_up[ikappa]+2,TwoSRD,              IRD);
            
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
      const int IRU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index));
      for (int TwoSRU=sector_spin_up[ikappa]-1; TwoSRU<=sector_spin_up[ikappa]+1; TwoSRU+=2){
         int dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+1, TwoSRU,                 IRU);
         int dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID);
         if ((dimRU>0) && (dimRD>0)){
         
            int fase = ((((sector_spin_up[ikappa] + sector_spin_down[ikappa] + 1)/2)%2)!=0)?-1:1;
            double factorB = fase * sqrt(3.0*(TwoSRU+1)) * Wigner::wigner6j( 1, 1, 2, TwoSRU, sector_spin_down[ikappa], sector_spin_up[ikappa] );
            
            double alpha;
            double * mem;
            
            if (TwoSRU == sector_spin_down[ikappa]){

               double factorA = - sqrt(0.5);
            
               double * BlockA = denA->gStorage( sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID );
               double * BlockB = denB->gStorage( sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID );
               
               mem = workmem;
               for (int cnt=0; cnt<dimRU*dimRD; cnt++){ mem[cnt] = factorA * BlockA[cnt] + factorB * BlockB[cnt]; }
               alpha = 1.0;
               
            } else {
            
               alpha = factorB;
               mem = denB->gStorage( sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID );
            
            }
            
            double * BlockTup = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+1, TwoSRU,                 IRU);
            double * BlockTdo = denT->gStorage(sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ID,               sector_nelec_up[ikappa]+3, sector_spin_down[ikappa], ID);
            
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

void CheMPS2::TensorQ::AddTermsCD(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2){

   if (moving_right){ AddTermsCDRight(denC, denD, denT, workmem, workmem2); }
   else{ AddTermsCDLeft(denC, denD, denT, workmem, workmem2); }
   
}

void CheMPS2::TensorQ::AddTermsCDRight(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IRD  = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
      int dimRU = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimRD = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IRD);
   
      //case 1
      const int ILD = Irreps::directProd(IRD,bk_up->gIrrep(index-1));
      for (int TwoSLD=sector_spin_down[ikappa]-1; TwoSLD<=sector_spin_down[ikappa]+1; TwoSLD+=2){
         int dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         int dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa], TwoSLD,              ILD);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            
            //first set to D
            int fase = ((((sector_spin_up[ikappa]+sector_spin_down[ikappa]+1)/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(TwoSLD+1)) * Wigner::wigner6j( 1, 2, 1, sector_spin_up[ikappa], sector_spin_down[ikappa], TwoSLD );
            double * block = denD->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], TwoSLD, ILD );
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSLD==sector_spin_up[ikappa]){
               factor = sqrt(0.5);
               block = denC->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], TwoSLD, ILD );
               int inc = 1;
               daxpy_(&dimLUxLD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa] );
            double * BlockTdo = denT->gStorage( sector_nelec_up[ikappa], TwoSLD, ILD, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IRD );
            
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
      const int ILU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index-1));
      for (int TwoSLU=sector_spin_up[ikappa]-1; TwoSLU<=sector_spin_up[ikappa]+1; TwoSLU+=2){
         int dimLU = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, TwoSLU,                 ILU);
         int dimLD = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IRD);
         
         if ((dimLU>0) && (dimLD>0)){
         
            int dimLUxLD = dimLU * dimLD;
            
            //first set to D
            int fase = ((((TwoSLU + sector_spin_down[ikappa])/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(sector_spin_up[ikappa]+1)) * Wigner::wigner6j( 1, 2, 1, sector_spin_down[ikappa], sector_spin_up[ikappa], TwoSLU );
            double * block = denD->gStorage( sector_nelec_up[ikappa]-1, TwoSLU, ILU, sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IRD );
            for (int cnt=0; cnt<dimLUxLD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSLU==sector_spin_down[ikappa]){
               fase = ((((sector_spin_down[ikappa]+1-sector_spin_up[ikappa])/2)%2)!=0)?-1:1;
               factor = fase * sqrt(0.5 * ( sector_spin_up[ikappa]+1.0 ) / ( sector_spin_down[ikappa] + 1.0 ) );
               block = denC->gStorage( sector_nelec_up[ikappa]-1, TwoSLU, ILU, sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IRD );
               int inc = 1;
               daxpy_(&dimLUxLD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sector_nelec_up[ikappa]-1, TwoSLU,                 ILU, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa] );
            double * BlockTdo = denT->gStorage( sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IRD, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IRD );
            
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

void CheMPS2::TensorQ::AddTermsCDLeft(TensorOperator * denC, TensorOperator * denD, TensorT * denT, double * workmem, double * workmem2){

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int ILD  = Irreps::directProd( sector_irrep_up[ikappa], n_irrep );
      int dimLU = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimLD = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD);
   
      //case 1
      const int IRU = Irreps::directProd(sector_irrep_up[ikappa],bk_up->gIrrep(index));
      for (int TwoSRU=sector_spin_up[ikappa]-1; TwoSRU<=sector_spin_up[ikappa]+1; TwoSRU+=2){
         int dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+1, TwoSRU,                 IRU);
         int dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            
            //first set to D
            int fase = ((((sector_spin_up[ikappa]+TwoSRU+3)/2)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0/(sector_spin_down[ikappa]+1.0)) * ( TwoSRU + 1 )
                          * Wigner::wigner6j( 1, 1, 2, TwoSRU, sector_spin_down[ikappa], sector_spin_up[ikappa] );
            double * block = denD->gStorage( sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD );
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSRU==sector_spin_down[ikappa]){
               factor = sqrt(0.5);
               block = denC->gStorage( sector_nelec_up[ikappa]+1, TwoSRU, IRU, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD );
               int inc = 1;
               daxpy_(&dimRUxRD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], sector_nelec_up[ikappa]+1, TwoSRU,                 IRU );
            double * BlockTdo = denT->gStorage( sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD,              sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], ILD );
            
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
      const int IRD = Irreps::directProd(ILD,bk_up->gIrrep(index));
      for (int TwoSRD=sector_spin_down[ikappa]-1; TwoSRD<=sector_spin_down[ikappa]+1; TwoSRD+=2){
         int dimRU = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         int dimRD = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, TwoSRD,              IRD);
         
         if ((dimRU>0) && (dimRD>0)){
         
            int dimRUxRD = dimRU * dimRD;
            
            //first set to D
            int fase = (((TwoSRD+1)%2)!=0)?-1:1;
            double factor = fase * sqrt(3.0*(TwoSRD+1.0)*(sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1.0))
                          * Wigner::wigner6j( 1, 1, 2, sector_spin_up[ikappa], TwoSRD, sector_spin_down[ikappa] );
            double * block = denD->gStorage( sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, TwoSRD, IRD );
            for (int cnt=0; cnt<dimRUxRD; cnt++){ workmem[cnt] = factor * block[cnt]; }
            
            //add C
            if (TwoSRD==sector_spin_up[ikappa]){
               fase = ((((sector_spin_up[ikappa]+1-sector_spin_down[ikappa])/2)%2)!=0)?-1:1;
               factor = fase * sqrt(0.5 * ( sector_spin_up[ikappa]+1.0 ) / ( sector_spin_down[ikappa] + 1.0 ) );
               block = denC->gStorage( sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, TwoSRD, IRD );
               int inc = 1;
               daxpy_(&dimRUxRD, &factor, block, &inc, workmem, &inc);
            }
            
            double * BlockTup = denT->gStorage( sector_nelec_up[ikappa],  sector_spin_up[ikappa],   sector_irrep_up[ikappa],sector_nelec_up[ikappa]+2,sector_spin_up[ikappa],sector_irrep_up[ikappa] );
            double * BlockTdo = denT->gStorage( sector_nelec_up[ikappa]+1,sector_spin_down[ikappa],ILD,             sector_nelec_up[ikappa]+2,TwoSRD,             IRD);
            
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



