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

#include "Cumulant.h"

CheMPS2::Cumulant::Cumulant(){ }

CheMPS2::Cumulant::~Cumulant(){ }

double lambda2_dmrg(Problem * prob, TwoDM * the2DM, const int dmrg_i, const int dmrg_j, const int dmrg_p, const int dmrg_q){

   const int irrep_i = prob->gIrrep( dmrg_i );
   const int irrep_j = prob->gIrrep( dmrg_j );
   const int irrep_p = prob->gIrrep( dmrg_p );
   const int irrep_q = prob->gIrrep( dmrg_q );
   
   if ( Irreps::directProd( irrep_i, irrep_j ) == Irreps::directProd( irrep_p, irrep_q ) ){
      const double value = the2DM->getTwoDMA_DMRG( dmrg_i, dmrg_j, dmrg_p, dmrg_q )
                         - the2DM->get1RDM_DMRG( dmrg_i, dmrg_p ) * the2DM->get1RDM_DMRG( dmrg_j, dmrg_q )
                         + the2DM->get1RDM_DMRG( dmrg_i, dmrg_q ) * the2DM->get1RDM_DMRG( dmrg_j, dmrg_p ) * 0.5;
      return value;
   }
   
   return 0.0;

}

double gamma4_ham(Problem * prob, ThreeDM * the3DM, TwoDM * the2DM, const int ham_i, const int ham_j, const int ham_k, const int ham_l,
                                                                    const int ham_p, const int ham_q, const int ham_r, const int ham_s){
   const double value = 0.0;
   return value;

}




