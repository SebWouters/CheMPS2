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

#include "TensorF0.h"
#include "Lapack.h"

CheMPS2::TensorF0::TensorF0( const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * denBK ) :
TensorOperator(boundary_index,
               0, // two_j
               0, // n_elec
               Idiff,
               moving_right,
               true,  // prime_last (doesn't matter for spin-0)
               false, // jw_phase (two 2nd quantized operators)
               denBK,
               denBK){ }

CheMPS2::TensorF0::~TensorF0(){ }

void CheMPS2::TensorF0::makenew(TensorT * denT){

   if (moving_right){ makenewRight(denT); }
   else{ makenewLeft( denT); }

}

void CheMPS2::TensorF0::makenew(TensorL * denL, TensorT * denT, double * workmem){

   if (moving_right){ makenewRight(denL, denT, workmem); }
   else{ makenewLeft( denL, denT, workmem); }

}

void CheMPS2::TensorF0::makenewRight(TensorT * denT){ //Idiff = Itrivial

   clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      const double sqrt_of_2 = sqrt(2.0);
      for (int geval=0; geval<3; geval++){
         int TwoSL, NL, IL;
         switch(geval){
            case 0:
               NL = sector_nelec_up[ikappa]-2;
               TwoSL = sector_spin_up[ikappa];
               IL = sector_irrep_up[ikappa];
               break;
            case 1:
               NL = sector_nelec_up[ikappa]-1;
               TwoSL = sector_spin_up[ikappa]-1;
               IL = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
               break;
            case 2:
               NL = sector_nelec_up[ikappa]-1;
               TwoSL = sector_spin_up[ikappa]+1;
               IL = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
               break;
         }
         int dimL  = bk_up->gCurrentDim(index-1, NL, TwoSL, IL);
         if (dimL>0){
         
            double * BlockT   = denT->gStorage(NL, TwoSL, IL, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);

            char trans = 'T';
            char notrans = 'N';
            double alpha = (geval==0)?sqrt_of_2:(0.5*sqrt_of_2);
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
         
         }
      }
   }

}

void CheMPS2::TensorF0::makenewLeft(TensorT * denT){ //Idiff = Itrivial

   clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      const double sqrt_of_2 = sqrt(2.0);
      for (int geval=0; geval<3; geval++){
         int TwoSR, NR, IR;
         switch(geval){
            case 0:
               NR = sector_nelec_up[ikappa]+2;
               TwoSR = sector_spin_up[ikappa];
               IR = sector_irrep_up[ikappa];
               break;
            case 1:
               NR = sector_nelec_up[ikappa]+1;
               TwoSR = sector_spin_up[ikappa]-1;
               IR = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               break;
            case 2:
               NR = sector_nelec_up[ikappa]+1;
               TwoSR = sector_spin_up[ikappa]+1;
               IR = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               break;
         }
         int dimR  = bk_up->gCurrentDim(index+1, NR, TwoSR, IR);
         if (dimR>0){
         
            double * BlockT   = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], NR, TwoSR, IR);

            char trans = 'T';
            char notrans = 'N';
            double alpha = sqrt_of_2;
            if (geval>=1){ alpha *= 0.5 * (TwoSR + 1.0) / (sector_spin_up[ikappa] + 1.0); }
            double beta = 1.0; //add
            dgemm_(&notrans,&trans,&dimL,&dimL,&dimR,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimL);
         
         }
      }
   }

}

void CheMPS2::TensorF0::makenewRight(TensorL * denL, TensorT * denT, double * workmem){

   clear();

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDR = Irreps::directProd(n_irrep,sector_irrep_up[ikappa]);
      int dimUR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      int dimDR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], IDR             );
      
      for (int geval=0; geval<4; geval++){
         int NLU,TwoSLU,ILU,TwoSLD,ILD; //NLD = NLU+1
         switch(geval){
            case 0:
               NLU = sector_nelec_up[ikappa]-1;
               TwoSLU = sector_spin_up[ikappa]-1;
               ILU = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index-1) );
               TwoSLD = sector_spin_up[ikappa];
               ILD = IDR;
               break;
            case 1:
               NLU = sector_nelec_up[ikappa]-1;
               TwoSLU = sector_spin_up[ikappa]+1;
               ILU = Irreps::directProd( sector_irrep_up[ikappa], bk_up->gIrrep(index-1) );
               TwoSLD = sector_spin_up[ikappa];
               ILD = IDR;
               break;
            case 2:
               NLU = sector_nelec_up[ikappa]-2;
               TwoSLU = sector_spin_up[ikappa];
               ILU = sector_irrep_up[ikappa];
               TwoSLD = sector_spin_up[ikappa]-1;
               ILD = Irreps::directProd( ILU, denL->get_irrep() );
               break;
            case 3:
               NLU = sector_nelec_up[ikappa]-2;
               TwoSLU = sector_spin_up[ikappa];
               ILU = sector_irrep_up[ikappa];
               TwoSLD = sector_spin_up[ikappa]+1;
               ILD = Irreps::directProd( ILU, denL->get_irrep() );
               break;
         }
         int dimLU = bk_up->gCurrentDim(index-1, NLU,   TwoSLU, ILU);
         int dimLD = bk_up->gCurrentDim(index-1, NLU+1, TwoSLD, ILD);
         if ((dimLU>0) && (dimLD>0)){
         
            double * BlockTup   = denT->gStorage(NLU,   TwoSLU, ILU, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
            double * BlockTdown = denT->gStorage(NLU+1, TwoSLD, ILD, sector_nelec_up[ikappa], sector_spin_up[ikappa], IDR);
            double * BlockL     = denL->gStorage(NLU,   TwoSLU, ILU, NLU+1,            TwoSLD,              ILD);
            
            //factor * Tup^T * L -> mem
            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               alpha = sqrt(0.5);
            } else {
               int fase = ((((sector_spin_up[ikappa] + 1 - TwoSLD)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt( 0.5 * ( TwoSLD + 1.0 ) / ( sector_spin_up[ikappa] + 1.0 ) );
            }
            double beta = 0.0; //set
            dgemm_(&trans,&notrans,&dimUR,&dimLD,&dimLU,&alpha,BlockTup,&dimLU,BlockL,&dimLU,&beta,workmem,&dimUR);
            
            //mem * Tdown -> storage
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&notrans,&dimUR,&dimDR,&dimLD,&alpha,workmem,&dimUR,BlockTdown,&dimLD,&beta,storage+kappa2index[ikappa],&dimUR);
         
         }
      }
   }

}

void CheMPS2::TensorF0::makenewLeft(TensorL * denL, TensorT * denT, double * workmem){

   clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDL = Irreps::directProd(n_irrep,sector_irrep_up[ikappa]);
      int dimUL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      int dimDL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa], sector_spin_up[ikappa], IDL             );
      
      for (int geval=0; geval<4; geval++){
         int NRU,TwoSRU,IRU,TwoSRD,IRD; //NRD = NRU+1
         switch(geval){
            case 0:
               NRU = sector_nelec_up[ikappa];
               TwoSRU = sector_spin_up[ikappa];
               IRU = sector_irrep_up[ikappa];
               TwoSRD = sector_spin_up[ikappa] - 1;
               IRD = Irreps::directProd( IRU , denL->get_irrep() );
               break;
            case 1:
               NRU = sector_nelec_up[ikappa];
               TwoSRU = sector_spin_up[ikappa];
               IRU = sector_irrep_up[ikappa];
               TwoSRD = sector_spin_up[ikappa] + 1;
               IRD = Irreps::directProd( IRU , denL->get_irrep() );
               break;
            case 2:
               NRU = sector_nelec_up[ikappa] + 1;
               TwoSRU = sector_spin_up[ikappa] - 1;
               IRU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               TwoSRD = sector_spin_up[ikappa];
               IRD = IDL;
               break;
            case 3:
               NRU = sector_nelec_up[ikappa] + 1;
               TwoSRU = sector_spin_up[ikappa] + 1;
               IRU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               TwoSRD = sector_spin_up[ikappa];
               IRD = IDL;
               break;
         }
         int dimRU = bk_up->gCurrentDim(index+1, NRU,   TwoSRU, IRU);
         int dimRD = bk_up->gCurrentDim(index+1, NRU+1, TwoSRD, IRD);
         if ((dimRU>0) && (dimRD>0)){
         
            double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], NRU,   TwoSRU, IRU);
            double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], IDL,              NRU+1, TwoSRD, IRD);
            double * BlockL     = denL->gStorage(NRU,              TwoSRU,              IRU,              NRU+1, TwoSRD, IRD);
            
            //factor * Tup * L -> mem
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               alpha = sqrt(0.5) * (TwoSRD + 1.0) / (sector_spin_up[ikappa] + 1.0);
            } else {
               int fase = ((((sector_spin_up[ikappa] - TwoSRU + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt( 0.5 * (TwoSRU + 1.0) / (sector_spin_up[ikappa] + 1.0) ) ;
            }
            double beta = 0.0; //set
            dgemm_(&notrans,&notrans,&dimUL,&dimRD,&dimRU,&alpha,BlockTup,&dimUL,BlockL,&dimRU,&beta,workmem,&dimUL);
            
            //mem * Tdown^T -> storage
            char trans = 'T';
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimRD,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
         }
      }
   }

}

