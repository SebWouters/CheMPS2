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

#include <math.h>

#include "TensorKM.h"
#include "Lapack.h"

CheMPS2::TensorKM::TensorKM( const int boundary_index, const char identity, const int Idiff, const SyBookkeeper * denBK ) :
TensorOperator( boundary_index,
                1, // two_j
                1, // n_elec
                Idiff,
                true, // TensorKM only exists moving left to right
                true, // prime_last
                false, // No jw_phase when updating (two-orbital mutual information!)
                denBK,
                denBK ){

   this->identity = identity;

}

CheMPS2::TensorKM::~TensorKM(){ }

void CheMPS2::TensorKM::construct(TensorT * denT){

   clear();
   
   if ( identity == 'K' ){

      for (int ikappa=0; ikappa<nKappa; ikappa++){

         const int IDR = Irreps::directProd( n_irrep, sector_irrep_up[ikappa] );
         int dimUR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimDR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IDR             );
         int dimL = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);

         if (dimL>0){

            double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa],   sector_spin_up[ikappa],   sector_irrep_up[ikappa]);
            double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa], sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IDR);

            char trans = 'T';
            char notrans = 'N';
            double alpha = 1.0;
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);

         }
      }
   }
   
   if ( identity == 'M' ){
   
      for (int ikappa=0; ikappa<nKappa; ikappa++){

         const int IDR = Irreps::directProd( n_irrep, sector_irrep_up[ikappa] );
         int dimUR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa],   sector_spin_up[ikappa],    sector_irrep_up[ikappa]);
         int dimDR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IDR             );
         int dimL = bk_up->gCurrentDim(index-1, sector_nelec_up[ikappa]-1,  sector_spin_down[ikappa], IDR);

         if (dimL>0){

            double * BlockTup   = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IDR, sector_nelec_up[ikappa],   sector_spin_up[ikappa],   sector_irrep_up[ikappa]);
            double * BlockTdown = denT->gStorage(sector_nelec_up[ikappa]-1, sector_spin_down[ikappa], IDR, sector_nelec_up[ikappa]+1, sector_spin_down[ikappa], IDR);

            char trans = 'T';
            char notrans = 'N';
            int fase = ((((sector_spin_down[ikappa] - sector_spin_up[ikappa] + 1)/2)%2)!=0)?-1:1;
            double alpha = fase * sqrt((sector_spin_up[ikappa]+1.0)/(sector_spin_down[ikappa]+1));
            double beta = 1.0; //add
            dgemm_(&trans,&notrans,&dimUR,&dimDR,&dimL,&alpha,BlockTup,&dimL,BlockTdown,&dimL,&beta,storage+kappa2index[ikappa],&dimUR);

         }
      }
   }

}


