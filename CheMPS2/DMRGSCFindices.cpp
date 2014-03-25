/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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
#include <iostream>

#include "DMRGSCFindices.h"

using std::cerr;
using std::cout;
using std::endl;

CheMPS2::DMRGSCFindices::DMRGSCFindices(const int L, const int Group, int * NOCCin, int * NDMRGin, int * NVIRTin){

   this->L = L;
   this->Group = Group;
   SymmInfo.setGroup(Group);
   this->Nirreps = SymmInfo.getNumberOfIrreps();
   
   NORB  = new int[Nirreps];
   NOCC  = new int[Nirreps];
   NDMRG = new int[Nirreps];
   NVIRT = new int[Nirreps];
   NORBcumulative  = new int[Nirreps+1];
   NDMRGcumulative = new int[Nirreps+1];
   
   int sum_check = 0;
   NORBcumulative[0]  = 0;
   NDMRGcumulative[0] = 0;
   for (int irrep=0; irrep<Nirreps; irrep++){
   
      if (NOCCin[irrep]  < 0){ cerr << "DMRGSCFindices::DMRGSCFindices : NOCC["  << irrep << "] = " << NOCCin[ irrep] << endl; }
      if (NDMRGin[irrep] < 0){ cerr << "DMRGSCFindices::DMRGSCFindices : NDMRG[" << irrep << "] = " << NDMRGin[irrep] << endl; }
      if (NVIRTin[irrep] < 0){ cerr << "DMRGSCFindices::DMRGSCFindices : NVIRT[" << irrep << "] = " << NVIRTin[irrep] << endl; }
      
      NORB[ irrep] = NOCCin[ irrep] + NDMRGin[irrep] + NVIRTin[irrep];
      NOCC[ irrep] = NOCCin[ irrep];
      NDMRG[irrep] = NDMRGin[irrep];
      NVIRT[irrep] = NVIRTin[irrep];
      
      sum_check += NORB[irrep];
      
      NORBcumulative[ irrep+1] = NORBcumulative[ irrep] + NORB[irrep];
      NDMRGcumulative[irrep+1] = NDMRGcumulative[irrep] + NDMRG[irrep];
      
   }
   if (sum_check != L){ cerr << "DMRGSCFindices::DMRGSCFindices : Sum over all OCC, DMRG and VIRT orbitals is not L." << endl; }
   
   irrepOfEachDMRGorbital = new int[NDMRGcumulative[Nirreps]];
   for (int irrep=0; irrep<Nirreps; irrep++){
      for (int cnt=0; cnt<NDMRG[irrep]; cnt++){
         irrepOfEachDMRGorbital[ NDMRGcumulative[irrep] + cnt ] = irrep;
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

}

int CheMPS2::DMRGSCFindices::getL() const{ return L; }

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


