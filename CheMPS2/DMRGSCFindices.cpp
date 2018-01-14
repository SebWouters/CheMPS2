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
#include <assert.h>
#include <iostream>
#include <algorithm>

#include "DMRGSCFindices.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::DMRGSCFindices::DMRGSCFindices(const int L, const int Group, int * NOCCin, int * NDMRGin, int * NVIRTin){

   this->L = L;
   SymmInfo.setGroup(Group);
   this->Nirreps = SymmInfo.getNumberOfIrreps();
   
   NORB  = new int[Nirreps];
   NOCC  = new int[Nirreps];
   NDMRG = new int[Nirreps];
   NVIRT = new int[Nirreps];
   NORBcumulative  = new int[Nirreps+1];
   NDMRGcumulative = new int[Nirreps+1];
   
   int totalNumOrbs = 0;
   NORBcumulative[0]  = 0;
   NDMRGcumulative[0] = 0;
   for (int irrep=0; irrep<Nirreps; irrep++){
   
      assert( NOCCin [irrep]>=0 );
      assert( NDMRGin[irrep]>=0 );
      assert( NVIRTin[irrep]>=0 );
      
      NORB[ irrep] = NOCCin[ irrep] + NDMRGin[irrep] + NVIRTin[irrep];
      NOCC[ irrep] = NOCCin[ irrep];
      NDMRG[irrep] = NDMRGin[irrep];
      NVIRT[irrep] = NVIRTin[irrep];
      
      totalNumOrbs += NORB[irrep];
      
      NORBcumulative[ irrep+1] = NORBcumulative[ irrep] + NORB[irrep];
      NDMRGcumulative[irrep+1] = NDMRGcumulative[irrep] + NDMRG[irrep];
      
   }
   assert( totalNumOrbs==L );
   
   irrepOfEachDMRGorbital = new int[NDMRGcumulative[Nirreps]];
   irrepOfEachOrbital = new int[L];
   for (int irrep=0; irrep<Nirreps; irrep++){
      for (int cnt=0; cnt<NDMRG[irrep]; cnt++){
         irrepOfEachDMRGorbital[ NDMRGcumulative[irrep] + cnt ] = irrep;
      }
      for (int cnt=0; cnt<NORB[irrep]; cnt++){
         irrepOfEachOrbital[ NORBcumulative[irrep] + cnt ] = irrep;
      }
   }

}

CheMPS2::DMRGSCFindices::~DMRGSCFindices(){

   delete [] NORB;
   delete [] NOCC;
   delete [] NDMRG;
   delete [] NVIRT;
   delete [] NORBcumulative;
   delete [] NDMRGcumulative;
   delete [] irrepOfEachDMRGorbital;
   delete [] irrepOfEachOrbital;

}

int CheMPS2::DMRGSCFindices::getL() const{ return L; }

int CheMPS2::DMRGSCFindices::getGroupNumber() const{ return SymmInfo.getGroupNumber(); }

int CheMPS2::DMRGSCFindices::getNirreps() const{ return Nirreps; }

int CheMPS2::DMRGSCFindices::getNORB(const int irrep) const{ return NORB[irrep]; }

int CheMPS2::DMRGSCFindices::getNOCC(const int irrep) const{ return NOCC[irrep]; }

int CheMPS2::DMRGSCFindices::getNDMRG(const int irrep) const{ return NDMRG[irrep]; }

int CheMPS2::DMRGSCFindices::getNVIRT(const int irrep) const{ return NVIRT[irrep]; }

int CheMPS2::DMRGSCFindices::getDMRGcumulative(const int irrep) const{ return NDMRGcumulative[irrep]; }

int CheMPS2::DMRGSCFindices::getOrigNOCCstart(const int irrep) const{ return NORBcumulative[irrep]; }

int CheMPS2::DMRGSCFindices::getOrigNDMRGstart(const int irrep) const{ return NORBcumulative[irrep] + NOCC[irrep]; }

int CheMPS2::DMRGSCFindices::getOrigNVIRTstart(const int irrep) const{ return NORBcumulative[irrep+1] - NVIRT[irrep]; }

int * CheMPS2::DMRGSCFindices::getIrrepOfEachDMRGorbital(){ return irrepOfEachDMRGorbital; }

int CheMPS2::DMRGSCFindices::getOrbitalIrrep(const int index) const{ return irrepOfEachOrbital[index]; }

int CheMPS2::DMRGSCFindices::getNOCCsum() const{

   int total = 0;
   for ( int irrep = 0; irrep < getNirreps(); irrep++ ){ total += getNOCC( irrep ); }
   return total;

}

int CheMPS2::DMRGSCFindices::getNORBmax() const{

   int the_max = 0;
   for ( int irrep = 0; irrep < getNirreps(); irrep++ ){ the_max = max( the_max, getNORB( irrep ) ); }
   return the_max;

}

int CheMPS2::DMRGSCFindices::getROTparamsize() const{

   int paramsize = 0;
   for ( int irrep = 0; irrep < getNirreps(); irrep++ ){ paramsize += ( getNORB( irrep ) * ( getNORB( irrep ) - 1 ) ) / 2; }
   return paramsize;

}

void CheMPS2::DMRGSCFindices::Print() const{

   cout << "NORB  = [ ";
   for (int irrep=0; irrep<Nirreps-1; irrep++){ cout << NORB[irrep] << " , "; }
   cout << NORB[Nirreps-1] << " ]" << endl;

   cout << "NOCC  = [ ";
   for (int irrep=0; irrep<Nirreps-1; irrep++){ cout << NOCC[irrep] << " , "; }
   cout << NOCC[Nirreps-1] << " ]" << endl;
   
   cout << "NDMRG = [ ";
   for (int irrep=0; irrep<Nirreps-1; irrep++){ cout << NDMRG[irrep] << " , "; }
   cout << NDMRG[Nirreps-1] << " ]" << endl;
   
   cout << "NVIRT = [ ";
   for (int irrep=0; irrep<Nirreps-1; irrep++){ cout << NVIRT[irrep] << " , "; }
   cout << NVIRT[Nirreps-1] << " ]" << endl;

}


