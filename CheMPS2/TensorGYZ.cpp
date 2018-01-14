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

#include "TensorGYZ.h"
#include "Lapack.h"

CheMPS2::TensorGYZ::TensorGYZ(const int boundary_index, const char identity, const SyBookkeeper * denBK) :
TensorOperator(boundary_index,
               0,     // two_j
               0,     // n_elec = 0
               0,     // n_irrep = I_trivial = 0
               true,  // TensorGYZ only exists moving left to right
               true,  // prime_last (doesn't matter for spin-0 tensors)
               false, // No jw_phase when updating (two-orbital mutual information!)
               denBK,
               denBK){

   this->identity = identity;

}

CheMPS2::TensorGYZ::~TensorGYZ(){ }

void CheMPS2::TensorGYZ::construct(TensorT * denT){

   for (int ikappa=0; ikappa<nKappa; ikappa++){

      int NL       = -1;
      int IL       = -1;
      int TwoSL    = -1;
      double alpha = 1.0;

      if (identity=='Y'){
         NL    = sector_nelec_up[ikappa];
         TwoSL = sector_spin_up[ikappa];
         IL    = sector_irrep_up[ikappa];
      }

      if (identity=='Z'){
         NL    = sector_nelec_up[ikappa]-2;
         TwoSL = sector_spin_up[ikappa];
         IL    = sector_irrep_up[ikappa];
      }

      if (identity=='G'){
         NL    = sector_nelec_up[ikappa]-1;
         TwoSL = sector_spin_up[ikappa]-1;
         IL    = Irreps::directProd( sector_irrep_up[ikappa] , bk_up->gIrrep(index-1) );
         alpha = sqrt(0.5);
      }

      int dimR = bk_up->gCurrentDim(index,   sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
      int dimL = bk_up->gCurrentDim(index-1, NL,               TwoSL,               IL);

      if (dimL>0){
         double * BlockT = denT->gStorage(NL, TwoSL, IL, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
         char trans = 'T';
         char notr = 'N';
         double beta = 0.0;
         dgemm_(&trans,&notr,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
      } else {
         for (int cnt=kappa2index[ikappa]; cnt<kappa2index[ikappa+1]; cnt++){ storage[cnt] = 0.0; }
      }

      if (identity=='G'){
         TwoSL = sector_spin_up[ikappa]+1;
         dimL  = bk_up->gCurrentDim(index-1, NL, TwoSL, IL);

         if (dimL>0){
            double * BlockT = denT->gStorage(NL, TwoSL, IL, sector_nelec_up[ikappa], sector_spin_up[ikappa], sector_irrep_up[ikappa]);
            char trans = 'T';
            char notr = 'N';
            double beta = 1.0; //ADD NOW!!!
            dgemm_(&trans,&notr,&dimR,&dimR,&dimL,&alpha,BlockT,&dimL,BlockT,&dimL,&beta,storage+kappa2index[ikappa],&dimR);
         }
      }

   }

}


