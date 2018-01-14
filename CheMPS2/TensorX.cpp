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

#include "TensorX.h"
#include "Lapack.h"
#include "Wigner.h"

CheMPS2::TensorX::TensorX(const int boundary_index, const bool moving_right, const SyBookkeeper * denBK, const Problem * Prob) :
TensorOperator(boundary_index,
               0, //two_j
               0, //n_elec
               0, //n_irrep
               moving_right,
               true,  //prime_last (doesn't matter for spin-0 tensors)
               false, //jw_phase (four 2nd quantized operators)
               denBK,
               denBK){

   this->Prob = Prob;

}

CheMPS2::TensorX::~TensorX(){ }

void CheMPS2::TensorX::update(TensorT * denT){

   if (moving_right){
      //PARALLEL
      #pragma omp parallel for schedule(dynamic)
      for (int ikappa=0; ikappa<nKappa; ikappa++){ makenewRight(ikappa, denT); }
   } else {
      //PARALLEL
      #pragma omp parallel for schedule(dynamic)
      for (int ikappa=0; ikappa<nKappa; ikappa++){ makenewLeft(ikappa, denT); }
   }

}

void CheMPS2::TensorX::update(TensorT * denT, TensorL ** Ltensors, TensorX * Xtensor, TensorQ * Qtensor, TensorOperator * Atensor, TensorOperator * Ctensor, TensorOperator * Dtensor){

   if (moving_right){
      //PARALLEL
      #pragma omp parallel
      {
      
         const bool doOtherThings = (index>1) ? true : false ;
         const int dimL     = (doOtherThings) ? bk_up->gMaxDimAtBound(index-1) : 0 ;
         const int dimR     = (doOtherThings) ? bk_up->gMaxDimAtBound(index)   : 0 ;
         double * workmemLL = (doOtherThings) ? new double[dimL*dimL] : NULL ;
         double * workmemLR = (doOtherThings) ? new double[dimL*dimR] : NULL ;
         double * workmemRR = (doOtherThings) ? new double[dimR*dimR] : NULL ;
      
         #pragma omp for schedule(dynamic)
         for (int ikappa=0; ikappa<nKappa; ikappa++){
            makenewRight(ikappa, denT);
            if (doOtherThings){
               update_moving_right(ikappa, Xtensor, denT, denT, workmemLR);
               addTermQLRight(ikappa, denT, Ltensors, Qtensor, workmemRR, workmemLR, workmemLL);
               addTermARight(ikappa, denT, Atensor, workmemRR, workmemLR);
               addTermCRight(ikappa, denT, Ctensor, workmemLR);
               addTermDRight(ikappa, denT, Dtensor, workmemLR);
            }
         }
         
         if (doOtherThings){
            delete [] workmemLL;
            delete [] workmemLR;
            delete [] workmemRR;
         }
      
      }
   } else {
      //PARALLEL
      #pragma omp parallel
      {
      
         const bool doOtherThings = (index<Prob->gL()-1) ? true : false ;
         const int dimL     = (doOtherThings) ? bk_up->gMaxDimAtBound(index)   : 0 ;
         const int dimR     = (doOtherThings) ? bk_up->gMaxDimAtBound(index+1) : 0 ;
         double * workmemLL = (doOtherThings) ? new double[dimL*dimL] : NULL ;
         double * workmemLR = (doOtherThings) ? new double[dimL*dimR] : NULL ;
         double * workmemRR = (doOtherThings) ? new double[dimR*dimR] : NULL ;
      
         #pragma omp for schedule(dynamic)
         for (int ikappa=0; ikappa<nKappa; ikappa++){
            makenewLeft(ikappa, denT);
            if (doOtherThings){
               update_moving_left(ikappa, Xtensor, denT, denT, workmemLR);
               addTermQLLeft(ikappa, denT, Ltensors, Qtensor, workmemLL, workmemLR, workmemRR);
               addTermALeft(ikappa, denT, Atensor, workmemLR, workmemLL);
               addTermCLeft(ikappa, denT, Ctensor, workmemLR);
               addTermDLeft(ikappa, denT, Dtensor, workmemLR);
            }
         }
         
         if (doOtherThings){
            delete [] workmemLL;
            delete [] workmemLR;
            delete [] workmemRR;
         }
      
      }
   }

}

void CheMPS2::TensorX::makenewRight(const int ikappa, TensorT * denT){
   
   int dimR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimL = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   double alpha = Prob->gMxElement(index-1,index-1,index-1,index-1);
      
   if ((dimL>0) && (fabs(alpha)>0.0)){
      
      double * BlockT = denT->gStorage(sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      char trans = 'T';
      char notr = 'N';
      double beta = 0.0; //because there's only 1 term contributing per kappa, we might as well set it i.o. adding
      dgemm_(&trans,&notr,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
      
   } else {
      for (int cnt=kappa2index[ikappa]; cnt<kappa2index[ikappa+1]; cnt++){ storage[cnt] = 0.0; }
   }

}

void CheMPS2::TensorX::makenewLeft(const int ikappa, TensorT * denT){
   
   int dimL = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimR = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   double alpha = Prob->gMxElement(index,index,index,index);

   if ((dimR>0) && (fabs(alpha)>0.0)){
   
      double * BlockT = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      char trans = 'T';
      char notr = 'N';
      double beta = 0.0; //set, not add (only 1 term)
      dgemm_(&notr,&trans,&dimL,&dimL,&dimR,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimL);
      
   } else {
      for (int cnt=kappa2index[ikappa]; cnt<kappa2index[ikappa+1]; cnt++){ storage[cnt] = 0.0; }
   }

}


void CheMPS2::TensorX::addTermQLRight(const int ikappa, TensorT * denT, TensorL ** Lprev, TensorQ * Qprev, double * workmemRR, double * workmemLR, double * workmemLL){

   int dimR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimTot = dimR * dimR;
   for (int cnt=0; cnt<dimTot; cnt++){ workmemRR[cnt] = 0.0; }
   
   for (int geval=0; geval<4; geval++){
      int NLup,TwoSLup,ILup,NLdown,TwoSLdown,ILdown;
      switch(geval){
         case 0:
            NLup = sector_nelec_up[ikappa]-1;
            TwoSLup = sector_spin_up[ikappa]-1;
            ILup = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
            NLdown = sector_nelec_up[ikappa];
            TwoSLdown = sector_spin_up[ikappa];
            ILdown = sector_irrep_up[ikappa];
            break;
         case 1:
            NLup = sector_nelec_up[ikappa]-1;
            TwoSLup = sector_spin_up[ikappa]+1;
            ILup = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
            NLdown = sector_nelec_up[ikappa];
            TwoSLdown = sector_spin_up[ikappa];
            ILdown = sector_irrep_up[ikappa];
            break;
         case 2:
            NLup = sector_nelec_up[ikappa]-2;
            TwoSLup = sector_spin_up[ikappa];
            ILup = sector_irrep_up[ikappa];
            NLdown = sector_nelec_up[ikappa]-1;
            TwoSLdown = sector_spin_up[ikappa]-1;
            ILdown = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
            break;
         case 3:
            NLup = sector_nelec_up[ikappa]-2;
            TwoSLup = sector_spin_up[ikappa];
            ILup = sector_irrep_up[ikappa];
            NLdown = sector_nelec_up[ikappa]-1;
            TwoSLdown = sector_spin_up[ikappa]+1;
            ILdown = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
            break;
      }
      int dimLup   = bk_up->gCurrentDim(index-1, NLup,   TwoSLup,   ILup);
      int dimLdown = bk_up->gCurrentDim(index-1, NLdown, TwoSLdown, ILdown);

      if ((dimLup>0) && (dimLdown>0)){
         double * BlockTup   = denT->gStorage(NLup,   TwoSLup,   ILup,   sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         double * BlockTdown = denT->gStorage(NLdown, TwoSLdown, ILdown, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         double * BlockQ    = Qprev->gStorage(NLup,   TwoSLup,   ILup,   NLdown,           TwoSLdown,           ILdown);
         
         double factor;
         double * ptr;
         if (geval<2){
         
            factor = 1.0;
            ptr = BlockQ;
            
         } else {
         
            int fase = ((((sector_spin_up[ikappa]+1-TwoSLdown)/2)%2)!=0)?-1:1;
            factor = fase * sqrt((TwoSLdown + 1.0)/(sector_spin_up[ikappa] + 1.0));
            
            int dimLupdown = dimLup * dimLdown;
            int inc = 1;
            ptr = workmemLL;
            dcopy_(&dimLupdown,BlockQ,&inc,ptr,&inc);
            
            for (int loca=0; loca<index-1; loca++){
               if (bk_up->gIrrep(index-1) == bk_up->gIrrep(loca)){
                  double alpha = Prob->gMxElement(loca, index-1, index-1, index-1);
                  double * BlockL = Lprev[index-2-loca]->gStorage(NLup, TwoSLup, ILup, NLdown, TwoSLdown, ILdown);
                  daxpy_(&dimLupdown,&alpha,BlockL,&inc,ptr,&inc);
               }
            }
            
         }

         //factor * Tup^T * L --> mem2 //set
         char trans = 'T';
         char notr = 'N';
         double beta = 0.0;
         dgemm_(&trans,&notr,&dimR,&dimLdown,&dimLup,&factor,BlockTup,&dimLup,ptr,&dimLup,&beta,workmemLR,&dimR);

         //mem2 * Tdown --> mem //add
         factor = 1.0;
         beta = 1.0;
         dgemm_(&notr,&notr,&dimR,&dimR,&dimLdown,&factor,workmemLR,&dimR,BlockTdown,&dimLdown,&beta,workmemRR,&dimR);

      }
   }
   //mem + mem^T --> storage
   for (int irow = 0; irow<dimR; irow++){
      for (int icol = irow; icol<dimR; icol++){
         workmemRR[irow + dimR * icol] += workmemRR[icol + dimR * irow];
         workmemRR[icol + dimR * irow] = workmemRR[irow + dimR * icol];
      }
   }
   int inc = 1;
   double alpha = 1.0;
   daxpy_(&dimTot,&alpha,workmemRR,&inc,storage + kappa2index[ikappa],&inc);

}

void CheMPS2::TensorX::addTermQLLeft(const int ikappa, TensorT * denT, TensorL ** Lprev, TensorQ * Qprev, double * workmemLL, double * workmemLR, double * workmemRR){

   int dimL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimTot = dimL * dimL;
   for (int cnt=0; cnt<dimTot; cnt++){ workmemLL[cnt] = 0.0; }
   
   for (int geval=0; geval<4; geval++){
      int NRup,TwoSRup,IRup,NRdown,TwoSRdown,IRdown;
      switch(geval){
         case 0:
            NRup = sector_nelec_up[ikappa];
            TwoSRup = sector_spin_up[ikappa];
            IRup = sector_irrep_up[ikappa];
            NRdown = sector_nelec_up[ikappa]+1;
            TwoSRdown = sector_spin_up[ikappa]-1;
            IRdown = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
            break;
         case 1:
            NRup = sector_nelec_up[ikappa];
            TwoSRup = sector_spin_up[ikappa];
            IRup = sector_irrep_up[ikappa];
            NRdown = sector_nelec_up[ikappa]+1;
            TwoSRdown = sector_spin_up[ikappa]+1;
            IRdown = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
            break;
         case 2:
            NRup = sector_nelec_up[ikappa]+1;
            TwoSRup = sector_spin_up[ikappa]-1;
            IRup = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
            NRdown = sector_nelec_up[ikappa]+2;
            TwoSRdown = sector_spin_up[ikappa];
            IRdown = sector_irrep_up[ikappa];
            break;
         case 3:
            NRup = sector_nelec_up[ikappa]+1;
            TwoSRup = sector_spin_up[ikappa]+1;
            IRup = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
            NRdown = sector_nelec_up[ikappa]+2;
            TwoSRdown = sector_spin_up[ikappa];
            IRdown = sector_irrep_up[ikappa];
            break;
      }
      int dimRup   = bk_up->gCurrentDim(index+1, NRup,   TwoSRup,   IRup);
      int dimRdown = bk_up->gCurrentDim(index+1, NRdown, TwoSRdown, IRdown);

      if ((dimRup>0) && (dimRdown>0)){
         double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], NRup,   TwoSRup,   IRup);
         double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], NRdown, TwoSRdown, IRdown);
         double * BlockQ    = Qprev->gStorage(NRup,             TwoSRup,             IRup,             NRdown, TwoSRdown, IRdown);
         
         double factor;
         double * ptr;
         if (geval<2){
         
            factor = (TwoSRdown + 1.0)/(sector_spin_up[ikappa] + 1.0);
            ptr = BlockQ;
            
         } else {
         
            int fase = ((((sector_spin_up[ikappa]+1-TwoSRup)/2)%2)!=0)?-1:1;
            factor = fase * sqrt((TwoSRup + 1.0)/(sector_spin_up[ikappa] + 1.0));
         
            int dimRupdown = dimRup * dimRdown;
            ptr = workmemRR;
            int inc = 1;
            dcopy_(&dimRupdown,BlockQ,&inc,ptr,&inc);
            
            for (int loca=index+1; loca<Prob->gL(); loca++){
               if (bk_up->gIrrep(index) == bk_up->gIrrep(loca)){
                  double alpha = Prob->gMxElement(index,index,index,loca);
                  double * BlockL = Lprev[loca-index-1]->gStorage(NRup, TwoSRup, IRup, NRdown, TwoSRdown, IRdown);
                  daxpy_(&dimRupdown,&alpha,BlockL,&inc,ptr,&inc);
               }
            }
         
         }

         //factor * Tup * L --> mem2 //set
         char notr = 'N';
         double beta = 0.0;//set
         dgemm_(&notr,&notr,&dimL,&dimRdown,&dimRup,&factor,BlockTup,&dimL,ptr,&dimRup,&beta,workmemLR,&dimL);
            
         //mem2 * Tdown^T --> mem //add
         char trans = 'T';
         factor = 1.0;
         beta = 1.0;
         dgemm_(&notr,&trans,&dimL,&dimL,&dimRdown,&factor,workmemLR,&dimL,BlockTdown,&dimL,&beta,workmemLL,&dimL);

      }
   }
   //mem + mem^T --> storage
   for (int irow = 0; irow<dimL; irow++){
      for (int icol = irow; icol<dimL; icol++){
         workmemLL[irow + dimL * icol] += workmemLL[icol + dimL * irow];
         workmemLL[icol + dimL * irow] = workmemLL[irow + dimL * icol];
      }
   }
   int inc = 1;
   double alpha = 1.0;
   daxpy_(&dimTot,&alpha,workmemLL,&inc,storage + kappa2index[ikappa],&inc);

}

void CheMPS2::TensorX::addTermARight(const int ikappa, TensorT * denT, TensorOperator * Aprev, double * workmemRR, double * workmemLR){

   int dimR     = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimLup   = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimLdown = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);

   if ((dimLup>0) && (dimLdown>0)){
      
      double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      double * BlockA    = Aprev->gStorage(sector_nelec_up[ikappa]-2, sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);

      //factor * Tup^T * A --> mem2 //set
      char trans = 'T';
      char notr = 'N';
      double factor = sqrt(2.0);
      double beta = 0.0; //set
      dgemm_(&trans,&notr,&dimR,&dimLdown,&dimLup,&factor,BlockTup,&dimLup,BlockA,&dimLup,&beta,workmemLR,&dimR);

      //mem2 * Tdown --> mem //set
      factor = 1.0;
      dgemm_(&notr,&notr,&dimR,&dimR,&dimLdown,&factor,workmemLR,&dimR,BlockTdown,&dimLdown,&beta,workmemRR,&dimR);
      
      //mem + mem^T --> storage
      for (int irow = 0; irow<dimR; irow++){
         for (int icol = irow; icol<dimR; icol++){
            workmemRR[irow + dimR * icol] += workmemRR[icol + dimR * irow];
            workmemRR[icol + dimR * irow] = workmemRR[irow + dimR * icol];
         }
      }
      int dimTot = dimR * dimR;
      int inc = 1;
      double alpha = 1.0;
      daxpy_(&dimTot,&alpha,workmemRR,&inc,storage + kappa2index[ikappa],&inc);

   }

}

void CheMPS2::TensorX::addTermALeft(const int ikappa, TensorT * denT, TensorOperator * Aprev, double * workmemLR, double * workmemLL){

   int dimL     = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimRup   = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   int dimRdown = bk_up->gCurrentDim(index+1, sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);

   if ((dimRup>0) && (dimRdown>0)){
      
      double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa],   sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      double * BlockA    = Aprev->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+2, sector_spin_up[ikappa], sector_irrep_up[ikappa]);

      //factor * Tup * A --> mem2 //set
      char notr = 'N';
      double factor = sqrt(2.0);
      double beta = 0.0; //set
      dgemm_(&notr,&notr,&dimL,&dimRdown,&dimRup,&factor,BlockTup,&dimL,BlockA,&dimRup,&beta,workmemLR,&dimL);

      //mem2 * Tdown^T --> mem //set
      char trans = 'T';
      factor = 1.0;
      dgemm_(&notr,&trans,&dimL,&dimL,&dimRdown,&factor,workmemLR,&dimL,BlockTdown,&dimL,&beta,workmemLL,&dimL);
      
      //mem + mem^T --> storage
      for (int irow = 0; irow<dimL; irow++){
         for (int icol = irow; icol<dimL; icol++){
            workmemLL[irow + dimL * icol] += workmemLL[icol + dimL * irow];
            workmemLL[icol + dimL * irow] = workmemLL[irow + dimL * icol];
         }
      }
      int dimTot = dimL * dimL;
      int inc = 1;
      double alpha = 1.0;
      daxpy_(&dimTot,&alpha,workmemLL,&inc,storage + kappa2index[ikappa],&inc);

   }

}

void CheMPS2::TensorX::addTermCRight(const int ikappa, TensorT * denT, TensorOperator * denC, double * workmemLR){

   int dimR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   for (int geval=0; geval<3; geval++){
      int NL, TwoSL, IL;
      switch(geval){
         case 0:
            NL = sector_nelec_up[ikappa]-1;
            TwoSL = sector_spin_up[ikappa]-1;
            IL = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index-1) );
            break;
         case 1:
            NL = sector_nelec_up[ikappa]-1;
            TwoSL = sector_spin_up[ikappa]+1;
            IL = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index-1) );
            break;
         case 2:
            NL = sector_nelec_up[ikappa]-2;
            TwoSL = sector_spin_up[ikappa];
            IL = sector_irrep_up[ikappa];
            break;
      }
      int dimL = bk_up->gCurrentDim(index-1,NL,TwoSL,IL);
      if (dimL>0){
      
         double * BlockC = denC->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
         double * BlockT = denT->gStorage(NL,TwoSL,IL,sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa]);

         double factor = (geval<2)?sqrt(0.5):sqrt(2.0);
         double beta = 0.0; //set
         char totrans = 'T';
         dgemm_(&totrans, &totrans, &dimR, &dimL, &dimL, &factor, BlockT, &dimL, BlockC, &dimL, &beta, workmemLR, &dimR);
         
         totrans = 'N';
         factor = 1.0;
         beta = 1.0; //add
         dgemm_(&totrans, &totrans, &dimR, &dimR, &dimL, &factor, workmemLR, &dimR, BlockT, &dimL, &beta, storage+kappa2index[ikappa], &dimR);

      }
   }

}

void CheMPS2::TensorX::addTermCLeft(const int ikappa, TensorT * denT, TensorOperator * denC, double * workmemLR){

   int dimL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   for (int geval=0; geval<3; geval++){
      int NR, TwoSR, IR;
      switch(geval){
         case 0:
            NR = sector_nelec_up[ikappa]+1;
            TwoSR = sector_spin_up[ikappa]-1;
            IR = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index) );
            break;
         case 1:
            NR = sector_nelec_up[ikappa]+1;
            TwoSR = sector_spin_up[ikappa]+1;
            IR = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index) );
            break;
         case 2:
            NR = sector_nelec_up[ikappa]+2;
            TwoSR = sector_spin_up[ikappa];
            IR = sector_irrep_up[ikappa];
            break;
      }
      int dimR = bk_up->gCurrentDim(index+1,NR,TwoSR,IR);
      if (dimR>0){
      
         double * BlockC = denC->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
         double * BlockT = denT->gStorage(sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa],NR,TwoSR,IR);

         double factor = (geval<2)?(sqrt(0.5)*(TwoSR+1.0)/(sector_spin_up[ikappa]+1.0)):sqrt(2.0);
         double beta = 0.0; //set
         char trans = 'T';
         char notr = 'N';
         dgemm_(&notr, &trans, &dimL, &dimR, &dimR, &factor, BlockT, &dimL, BlockC, &dimR, &beta, workmemLR, &dimL);
         
         factor = 1.0;
         beta = 1.0; //add
         dgemm_(&notr, &trans, &dimL, &dimL, &dimR, &factor, workmemLR, &dimL, BlockT, &dimL, &beta, storage+kappa2index[ikappa], &dimL);

      }
   }

}

void CheMPS2::TensorX::addTermDRight(const int ikappa, TensorT * denT, TensorOperator * denD, double * workmemLR){

   int dimR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   
   const int IL = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index-1) );
   const int NL = sector_nelec_up[ikappa]-1;
   
   for (int geval=0; geval<4; geval++){
      int TwoSLup, TwoSLdown;
      switch(geval){
         case 0:
            TwoSLup   = sector_spin_up[ikappa]-1;
            TwoSLdown = sector_spin_up[ikappa]-1;
            break;
         case 1:
            TwoSLup   = sector_spin_up[ikappa]+1;
            TwoSLdown = sector_spin_up[ikappa]-1;
            break;
         case 2:
            TwoSLup   = sector_spin_up[ikappa]-1;
            TwoSLdown = sector_spin_up[ikappa]+1;
            break;
         case 3:
            TwoSLup   = sector_spin_up[ikappa]+1;
            TwoSLdown = sector_spin_up[ikappa]+1;
            break;
      }
      
      int dimLup   = bk_up->gCurrentDim(index-1,NL,TwoSLup,  IL);
      int dimLdown = bk_up->gCurrentDim(index-1,NL,TwoSLdown,IL);
      
      if ((dimLup>0) && (dimLdown>0)){
      
         double * BlockD = denD->gStorage(NL,TwoSLdown,IL,NL,TwoSLup,IL);
         double * BlockTup   = denT->gStorage(NL,TwoSLup,  IL,sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa]);
         double * BlockTdown = (TwoSLup==TwoSLdown)? BlockTup : denT->gStorage(NL,TwoSLdown,IL,sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa]);
         
         int fase = ((((TwoSLdown + sector_spin_up[ikappa] + 1)/2)%2)!=0)?-1:1;
         double factor = fase * sqrt(3.0 * (TwoSLup+1))
                       * Wigner::wigner6j( 1, 1, 2, TwoSLup, TwoSLdown, sector_spin_up[ikappa] );
         double beta = 0.0; //set
         char totrans = 'T';
         dgemm_(&totrans, &totrans, &dimR, &dimLdown, &dimLup, &factor, BlockTup, &dimLup, BlockD, &dimLdown, &beta, workmemLR, &dimR);
         
         totrans = 'N';
         factor = 1.0;
         beta = 1.0; //add
         dgemm_(&totrans, &totrans, &dimR, &dimR, &dimLdown, &factor, workmemLR, &dimR, BlockTdown, &dimLdown, &beta, storage+kappa2index[ikappa], &dimR);
         
      }
   }

}

void CheMPS2::TensorX::addTermDLeft(const int ikappa, TensorT * denT, TensorOperator * denD, double * workmemLR){

   int dimL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
   
   const int NR = sector_nelec_up[ikappa]+1;
   const int IR = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index) );
   
   for (int geval=0; geval<4; geval++){
      int TwoSRup, TwoSRdown;
      switch(geval){
         case 0:
            TwoSRup   = sector_spin_up[ikappa] - 1;
            TwoSRdown = sector_spin_up[ikappa] - 1;
            break;
         case 1:
            TwoSRup   = sector_spin_up[ikappa] + 1;
            TwoSRdown = sector_spin_up[ikappa] - 1;
            break;
         case 2:
            TwoSRup   = sector_spin_up[ikappa] - 1;
            TwoSRdown = sector_spin_up[ikappa] + 1;
            break;
         case 3:
            TwoSRup   = sector_spin_up[ikappa] + 1;
            TwoSRdown = sector_spin_up[ikappa] + 1;
            break;
      }
      
      int dimRup   = bk_up->gCurrentDim(index+1,NR,TwoSRup,  IR);
      int dimRdown = bk_up->gCurrentDim(index+1,NR,TwoSRdown,IR);
      
      if ((dimRup>0) && (dimRdown>0)){
      
         double * BlockD = denD->gStorage(NR,TwoSRdown,IR,NR,TwoSRup,IR);
         double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa],NR,TwoSRup,IR);
         double * BlockTdown = (TwoSRup == TwoSRdown)? BlockTup : denT->gStorage(sector_nelec_up[ikappa],sector_spin_up[ikappa],sector_irrep_up[ikappa],NR,TwoSRdown,IR);
         
         int fase = ((((sector_spin_up[ikappa] + TwoSRdown + 3)/2)%2)!=0)?-1:1;
         double factor = fase*sqrt(3.0 *(TwoSRup+1))*((TwoSRdown + 1.0)/(sector_spin_up[ikappa]+1.0))
                       * Wigner::wigner6j( 1, 1, 2, TwoSRup, TwoSRdown, sector_spin_up[ikappa] );
         double beta = 0.0; //set
         char trans = 'T';
         char notr = 'N';
         dgemm_(&notr, &trans, &dimL, &dimRdown, &dimRup, &factor, BlockTup, &dimL, BlockD, &dimRdown, &beta, workmemLR, &dimL);
         
         factor = 1.0;
         beta = 1.0; //add
         dgemm_(&notr, &trans, &dimL, &dimL, &dimRdown, &factor, workmemLR, &dimL, BlockTdown, &dimL, &beta, storage+kappa2index[ikappa], &dimL);

      }
   }

}


