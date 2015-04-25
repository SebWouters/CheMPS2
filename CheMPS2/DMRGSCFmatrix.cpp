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

#include "DMRGSCFmatrix.h"

CheMPS2::DMRGSCFmatrix::DMRGSCFmatrix(DMRGSCFindices * iHandler_in){

   iHandler = iHandler_in;
   
   entries = new double*[ iHandler->getNirreps() ];
   for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
      entries[ irrep ] = new double[ iHandler->getNORB( irrep ) * iHandler->getNORB( irrep ) ];
   }
   
}

CheMPS2::DMRGSCFmatrix::~DMRGSCFmatrix(){

   for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
      delete [] entries[ irrep ];
   }
   delete [] entries;
   
}

void CheMPS2::DMRGSCFmatrix::set(const int irrep, const int p, const int q, const double val){

   entries[ irrep ][ p + iHandler->getNORB( irrep ) * q ] = val;

}

void CheMPS2::DMRGSCFmatrix::clear(){

   for (int irrep = 0; irrep < iHandler->getNirreps(); irrep++){
      for (int counter = 0; counter < iHandler->getNORB( irrep ) * iHandler->getNORB( irrep ); counter++){
         entries[ irrep ][ counter ] = 0.0;
      }
   }

}

double CheMPS2::DMRGSCFmatrix::get(const int irrep, const int p, const int q) const{

   return entries[ irrep ][ p + iHandler->getNORB( irrep ) * q ];

}

double * CheMPS2::DMRGSCFmatrix::getBlock(const int irrep){

   return entries[ irrep ];

}

