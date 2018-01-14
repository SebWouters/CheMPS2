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

#include "DMRGSCFwtilde.h"

CheMPS2::DMRGSCFwtilde::DMRGSCFwtilde(DMRGSCFindices * iHandler_in){

   iHandler = iHandler_in;
   
   Nocc_dmrg = new int[ iHandler->getNirreps() ];
   for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
      Nocc_dmrg[ irrep ] = iHandler->getNOCC( irrep ) + iHandler->getNDMRG( irrep );
   }
   
   wmattilde = new double***[ iHandler->getNirreps() ];
   for (int irrep_pq = 0; irrep_pq < iHandler->getNirreps(); irrep_pq++){
      wmattilde[ irrep_pq ] = new double**[ iHandler->getNirreps() ];
      for (int irrep_rs = 0; irrep_rs < iHandler->getNirreps(); irrep_rs++){
         const unsigned int sizeblock_pr = Nocc_dmrg[ irrep_pq ] * Nocc_dmrg[ irrep_rs ];
         const unsigned int sizeblock_qs = iHandler->getNORB( irrep_pq ) * iHandler->getNORB( irrep_rs );
         wmattilde[ irrep_pq ][ irrep_rs ] = new double*[ sizeblock_pr ];
         for (unsigned int combined_pr = 0; combined_pr < sizeblock_pr; combined_pr++){
            wmattilde[ irrep_pq ][ irrep_rs ][ combined_pr ] = new double[ sizeblock_qs ];
         }
      }
   }
   
}

CheMPS2::DMRGSCFwtilde::~DMRGSCFwtilde(){

   for (int irrep_pq = 0; irrep_pq < iHandler->getNirreps(); irrep_pq++){
      for (int irrep_rs = 0; irrep_rs < iHandler->getNirreps(); irrep_rs++){
         const unsigned int sizeblock_pr = Nocc_dmrg[ irrep_pq ] * Nocc_dmrg[ irrep_rs ];
         for (unsigned int combined_pr = 0; combined_pr < sizeblock_pr; combined_pr++){
            delete [] wmattilde[ irrep_pq ][ irrep_rs ][ combined_pr ];
         }
         delete [] wmattilde[ irrep_pq ][ irrep_rs ];
      }
      delete [] wmattilde[ irrep_pq ];
   }
   delete [] wmattilde;

   delete [] Nocc_dmrg;
   
}

void CheMPS2::DMRGSCFwtilde::clear(){

   for (int irrep_pq = 0; irrep_pq < iHandler->getNirreps(); irrep_pq++){
      for (int irrep_rs = 0; irrep_rs < iHandler->getNirreps(); irrep_rs++){
         const unsigned int sizeblock_pr = Nocc_dmrg[ irrep_pq ] * Nocc_dmrg[ irrep_rs ];
         const unsigned int sizeblock_qs = iHandler->getNORB( irrep_pq ) * iHandler->getNORB( irrep_rs );
         for (unsigned int combined_pr = 0; combined_pr < sizeblock_pr; combined_pr++){
            for (unsigned int combined_qs = 0; combined_qs < sizeblock_qs; combined_qs++){
               wmattilde[ irrep_pq ][ irrep_rs ][ combined_pr ][ combined_qs ] = 0.0;
            }
         }
      }
   }
   
}

void CheMPS2::DMRGSCFwtilde::set(const int irrep_pq, const int irrep_rs, const int p, const int q, const int r, const int s, const double val){

   wmattilde[ irrep_pq ][ irrep_rs ][ p + Nocc_dmrg[ irrep_pq ] * r ][ q + iHandler->getNORB(irrep_pq) * s ] = val;

}

double CheMPS2::DMRGSCFwtilde::get(const int irrep_pq, const int irrep_rs, const int p, const int q, const int r, const int s) const{

   return wmattilde[ irrep_pq ][ irrep_rs ][ p + Nocc_dmrg[ irrep_pq ] * r ][ q + iHandler->getNORB(irrep_pq) * s ];

}

double * CheMPS2::DMRGSCFwtilde::getBlock(const int irrep_pq, const int irrep_rs, const int p, const int r){

   return wmattilde[ irrep_pq ][ irrep_rs ][ p + Nocc_dmrg[ irrep_pq ] * r ];

}

