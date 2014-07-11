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

#include <stdlib.h> /*srand, rand*/
#include <iostream>
#include <time.h> /*time*/
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "CASSCF.h"

using namespace std;

int main(void){

   cout.precision(15);
   srand(time(NULL));
   
   //The path to the matrix elements
   string matrixelements = "../../tests/matrixelements/O2_CCPVDZ.dat";
   struct stat stFileInfo;
   int intStat = stat(matrixelements.c_str(),&stFileInfo);
   if (intStat != 0){
      cout << "Please set the correct relative path to tests/matrixelements/O2_CCPVDZ.dat in tests/test6.cpp for the compiled binary test6 to work." << endl;
      return 628788;
   }
   
   //Setup CASSCF
   CheMPS2::CASSCF koekoek(matrixelements);
   
   int * Nocc  = new int[koekoek.getNumberOfIrreps()];
   int * NDMRG = new int[koekoek.getNumberOfIrreps()];
   int * Nvirt = new int[koekoek.getNumberOfIrreps()];
   
   Nocc[0] = Nocc[5] = 1;
   Nocc[1] = Nocc[2] = Nocc[3] = Nocc[4] = Nocc[6] = Nocc[7] = 0;
   
   NDMRG[0] = NDMRG[5] = 2;
   NDMRG[1] = NDMRG[4] = 0;
   NDMRG[2] = NDMRG[3] = NDMRG[6] = NDMRG[7] = 2;
   
   Nvirt[0] = Nvirt[5] = 4;
   Nvirt[1] = Nvirt[2] = Nvirt[3] = Nvirt[4] = Nvirt[6] = Nvirt[7] = 1;
   
   koekoek.setupStart(Nocc,NDMRG,Nvirt);
   
   delete [] Nocc;
   delete [] NDMRG;
   delete [] Nvirt;
   
   //Setup symmetry sector
   int N = 16;
   int TwoS = 0;
   int Irrep = 0;
   
   //Setup convergence scheme
   CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme(1);
   int D = 1000;
   double Econv = 1e-8;
   int maxSweeps = 20;
   double noisePrefactor = 0.0;
   OptScheme->setInstruction(0,D,Econv,maxSweeps,noisePrefactor);

   //Run CASSCF
   int rootNum = 1; //Ground state only
   bool doDIIS = true;
   double Energy = koekoek.doCASSCFnewtonraphson(N, TwoS, Irrep, OptScheme, rootNum, doDIIS);
   
   //Clean up
   if (CheMPS2::DMRGSCF_storeUnitary){ koekoek.deleteStoredUnitary(); }
   if (CheMPS2::DMRGSCF_storeDIIS){ koekoek.deleteStoredDIIS(); }
   delete OptScheme;
   
   //Check succes
   bool success = (fabs(Energy + 149.690485081255) < 1e-10) ? true : false;
   cout << "================> Did test 6 succeed : ";
   if (success){ cout << "yes" << endl; }
   else { cout << "no" << endl; }

   return 0;

}


