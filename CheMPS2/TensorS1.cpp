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

#include "TensorS1.h"
#include "Lapack.h"
#include "Wigner.h"

CheMPS2::TensorS1::TensorS1(const int boundary_index, const int Idiff, const bool moving_right, const SyBookkeeper * denBK) :
TensorOperator(boundary_index,
               2, // two_j
               2, // n_elec
               Idiff,
               moving_right,
               true,  // prime_last
               false, // jw_phase (two 2nd quantized operators)
               denBK,
               denBK){ }

CheMPS2::TensorS1::~TensorS1(){ }

void CheMPS2::TensorS1::makenew(TensorL * denL, TensorT * denT, double * workmem){

   if (moving_right){ makenewRight(denL, denT, workmem); }
   else{ makenewLeft( denL, denT, workmem); }

}

void CheMPS2::TensorS1::makenewRight(TensorL * denL, TensorT * denT, double * workmem){

   clear();

   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDR = Irreps::directProd(n_irrep,sector_irrep_up[ikappa]);
      int dimUR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimDR = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+2, sector_spin_down[ikappa], IDR             );
      
      for (int geval=0; geval<4; geval++){
         int NLU,TwoSLU,ILU,TwoSLD,ILD; //NLD = NLU+1
         switch(geval){
            case 0:
               NLU = sector_nelec_up[ikappa];
               TwoSLU = sector_spin_up[ikappa];
               ILU = sector_irrep_up[ikappa];
               TwoSLD = sector_spin_down[ikappa]-1;
               ILD = Irreps::directProd( ILU, denL->get_irrep() );
               break;
            case 1:
               NLU = sector_nelec_up[ikappa];
               TwoSLU = sector_spin_up[ikappa];
               ILU = sector_irrep_up[ikappa];
               TwoSLD = sector_spin_down[ikappa]+1;
               ILD = Irreps::directProd( ILU, denL->get_irrep() );
               break;
            case 2:
               NLU = sector_nelec_up[ikappa]-1;
               TwoSLU = sector_spin_up[ikappa]-1;
               ILU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
               TwoSLD = sector_spin_down[ikappa];
               ILD = IDR;
               break;
            case 3:
               NLU = sector_nelec_up[ikappa]-1;
               TwoSLU = sector_spin_up[ikappa]+1;
               ILU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
               TwoSLD = sector_spin_down[ikappa];
               ILD = IDR;
               break;
         }
         int dimLU = bk_up->gCurrentDim(index-1, NLU,   TwoSLU, ILU);
         int dimLD = bk_up->gCurrentDim(index-1, NLU+1, TwoSLD, ILD);
         if ((dimLU>0) && (dimLD>0) && (abs(TwoSLU-TwoSLD)<2)){
         
            double * BlockTup   = denT->gStorage(NLU,   TwoSLU, ILU, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
            double * BlockTdown = denT->gStorage(NLU+1, TwoSLD, ILD, sector_nelec_up[ikappa]+2, sector_spin_down[ikappa], IDR             );
            double * BlockL     = denL->gStorage(NLU,   TwoSLU, ILU, NLU+1,              TwoSLD,                 ILD             );
            
            //factor * Tup^T * L -> mem
            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               int fase = ((((sector_spin_up[ikappa] + sector_spin_down[ikappa] + 2)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0*(TwoSLD+1)) * Wigner::wigner6j(1,1,2,sector_spin_up[ikappa],sector_spin_down[ikappa],TwoSLD);
            } else {
               int fase = ((((TwoSLU + sector_spin_down[ikappa] + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0*(sector_spin_up[ikappa]+1)) * Wigner::wigner6j(1,1,2,sector_spin_up[ikappa],sector_spin_down[ikappa],TwoSLU);
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

void CheMPS2::TensorS1::makenewLeft(TensorL * denL, TensorT * denT, double * workmem){

   clear();
   
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      const int IDL = Irreps::directProd(n_irrep,sector_irrep_up[ikappa]);
      int dimUL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
      int dimDL = bk_up->gCurrentDim(index, sector_nelec_up[ikappa]+2, sector_spin_down[ikappa], IDL             );
      
      for (int geval=0; geval<4; geval++){
         int NRU,TwoSRU,IRU,TwoSRD,IRD; //NRD = NRU+1
         switch(geval){
            case 0:
               NRU = sector_nelec_up[ikappa]+1;
               TwoSRU = sector_spin_up[ikappa]-1;
               IRU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               TwoSRD = sector_spin_down[ikappa];
               IRD = IDL;
               break;
            case 1:
               NRU = sector_nelec_up[ikappa]+1;
               TwoSRU = sector_spin_up[ikappa]+1;
               IRU = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index) );
               TwoSRD = sector_spin_down[ikappa];
               IRD = IDL;
               break;
            case 2:
               NRU = sector_nelec_up[ikappa]+2;
               TwoSRU = sector_spin_up[ikappa];
               IRU = sector_irrep_up[ikappa];
               TwoSRD = sector_spin_down[ikappa]-1;
               IRD = Irreps::directProd( sector_irrep_up[ikappa] , denL->get_irrep() );
               break;
            case 3:
               NRU = sector_nelec_up[ikappa]+2;
               TwoSRU = sector_spin_up[ikappa];
               IRU = sector_irrep_up[ikappa];
               TwoSRD = sector_spin_down[ikappa]+1;
               IRD = Irreps::directProd( sector_irrep_up[ikappa] , denL->get_irrep() );
               break;
         }
         int dimRU = bk_up->gCurrentDim(index+1, NRU,   TwoSRU, IRU);
         int dimRD = bk_up->gCurrentDim(index+1, NRU+1, TwoSRD, IRD);
         if ((dimRU>0) && (dimRD>0) && (abs(TwoSRD-TwoSRU)<2)){
         
            double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa], NRU,   TwoSRU, IRU);
            double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa]+2, sector_spin_down[ikappa], IDL,              NRU+1, TwoSRD, IRD);
            double * BlockL     = denL->gStorage(NRU,                TwoSRU,                 IRU,              NRU+1, TwoSRD, IRD);
            
            //factor * Tup * L -> mem
            char notrans = 'N';
            double alpha = 1.0;
            if (geval<=1){
               int fase = ((((sector_spin_up[ikappa] + sector_spin_down[ikappa] + 2)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0 * (TwoSRU+1)) * Wigner::wigner6j( 1, 1, 2, sector_spin_up[ikappa], sector_spin_down[ikappa], TwoSRU );
            } else {
               int fase = ((((sector_spin_up[ikappa] + TwoSRD + 1)/2)%2)!=0)?-1:1;
               alpha = fase * sqrt(3.0 / (sector_spin_down[ikappa] + 1.0)) * (TwoSRD + 1)
                     * Wigner::wigner6j( 1, 1, 2, sector_spin_up[ikappa], sector_spin_down[ikappa], TwoSRD );
            }
            double beta = 0.0; //set
            dgemm_(&notrans,&notrans,&dimUL,&dimRD,&dimRU,&alpha,BlockTup,&dimUL,BlockL,&dimRU,&beta,workmem,&dimUL);
            
            //mem * Tdown -> storage
            char trans = 'T';
            alpha = 1.0;
            beta = 1.0; // add
            dgemm_(&notrans,&trans,&dimUL,&dimDL,&dimRD,&alpha,workmem,&dimUL,BlockTdown,&dimDL,&beta,storage+kappa2index[ikappa],&dimUL);
         
         }
      }
   }

}

