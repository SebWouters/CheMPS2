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

#include "DMRGSCFoptions.h"

CheMPS2::DMRGSCFoptions::DMRGSCFoptions(){

   DoDIIS             = CheMPS2::DMRGSCF_doDIIS;
   DIISGradientBranch = CheMPS2::DMRGSCF_DIISgradientBranch;
   NumDIISVecs        = CheMPS2::DMRGSCF_numDIISvecs;
   StoreDIIS          = CheMPS2::DMRGSCF_storeDIIS;
   DIISStorageName    = CheMPS2::DMRGSCF_diis_storage_name;
   
   MaxIterations      = CheMPS2::DMRGSCF_maxIterations;
   GradientThreshold  = CheMPS2::DMRGSCF_gradientNormThreshold;
   StoreUnitary       = CheMPS2::DMRGSCF_storeUnitary;
   UnitaryStorageName = CheMPS2::DMRGSCF_unitary_storage_name;
   StateAveraging     = CheMPS2::DMRGSCF_stateAveraged;
   
   WhichActiveSpace   = CheMPS2::DMRGSCF_whichActiveSpace;
   DumpCorrelations   = CheMPS2::DMRGSCF_dumpCorrelations;
   StartLocRandom     = CheMPS2::DMRGSCF_startLocRandom;

}

CheMPS2::DMRGSCFoptions::~DMRGSCFoptions(){ }

bool   CheMPS2::DMRGSCFoptions::getDoDIIS() const{             return DoDIIS;             }
double CheMPS2::DMRGSCFoptions::getDIISGradientBranch() const{ return DIISGradientBranch; }
int    CheMPS2::DMRGSCFoptions::getNumDIISVecs() const{        return NumDIISVecs;        }
bool   CheMPS2::DMRGSCFoptions::getStoreDIIS() const{          return StoreDIIS;          }
string CheMPS2::DMRGSCFoptions::getDIISStorageName() const{    return DIISStorageName;    }
int    CheMPS2::DMRGSCFoptions::getMaxIterations() const{      return MaxIterations;      }
double CheMPS2::DMRGSCFoptions::getGradientThreshold() const{  return GradientThreshold;  }
bool   CheMPS2::DMRGSCFoptions::getStoreUnitary() const{       return StoreUnitary;       }
string CheMPS2::DMRGSCFoptions::getUnitaryStorageName() const{ return UnitaryStorageName; }
bool   CheMPS2::DMRGSCFoptions::getStateAveraging() const{     return StateAveraging;     }
int    CheMPS2::DMRGSCFoptions::getWhichActiveSpace() const{   return WhichActiveSpace;   }
bool   CheMPS2::DMRGSCFoptions::getDumpCorrelations() const{   return DumpCorrelations;   }
bool   CheMPS2::DMRGSCFoptions::getStartLocRandom() const{     return StartLocRandom;     }

void CheMPS2::DMRGSCFoptions::setDoDIIS(const bool DoDIIS_in){                           DoDIIS             = DoDIIS_in;             }
void CheMPS2::DMRGSCFoptions::setDIISGradientBranch(const double DIISGradientBranch_in){ DIISGradientBranch = DIISGradientBranch_in; }
void CheMPS2::DMRGSCFoptions::setNumDIISVecs(const int NumDIISVecs_in){                  NumDIISVecs        = NumDIISVecs_in;        }
void CheMPS2::DMRGSCFoptions::setStoreDIIS(const bool StoreDIIS_in){                     StoreDIIS          = StoreDIIS_in;          }
void CheMPS2::DMRGSCFoptions::setDIISStorageName(const string DIISStorageName_in){       DIISStorageName    = DIISStorageName_in;    }
void CheMPS2::DMRGSCFoptions::setMaxIterations(const int MaxIterations_in){              MaxIterations      = MaxIterations_in;      }
void CheMPS2::DMRGSCFoptions::setGradientThreshold(const double GradientThreshold_in){   GradientThreshold  = GradientThreshold_in;  }
void CheMPS2::DMRGSCFoptions::setStoreUnitary(const bool StoreUnitary_in){               StoreUnitary       = StoreUnitary_in;       }
void CheMPS2::DMRGSCFoptions::setUnitaryStorageName(const string UnitaryStorageName_in){ UnitaryStorageName = UnitaryStorageName_in; }
void CheMPS2::DMRGSCFoptions::setStateAveraging(const bool StateAveraging_in){           StateAveraging     = StateAveraging_in;     }
void CheMPS2::DMRGSCFoptions::setWhichActiveSpace(const int WhichActiveSpace_in){        WhichActiveSpace   = WhichActiveSpace_in;   }
void CheMPS2::DMRGSCFoptions::setDumpCorrelations(const bool DumpCorrelations_in){       DumpCorrelations   = DumpCorrelations_in;   }
void CheMPS2::DMRGSCFoptions::setStartLocRandom(const bool StartLocRandom_in){           StartLocRandom     = StartLocRandom_in;     }



