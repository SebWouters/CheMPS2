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

#include <iostream>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "EdmistonRuedenberg.h"

using namespace std;

int main(void){

   cout.precision(15);
   
   //The path to the matrix elements
   string matrixelements = "../../tests/matrixelements/N2_CCPVDZ.dat";
   struct stat stFileInfo;
   int intStat = stat(matrixelements.c_str(),&stFileInfo);
   if (intStat != 0){
      cout << "Please set the correct relative path to tests/matrixelements/N2_CCPVDZ.dat in tests/test9.cpp for the compiled binary test9 to work." << endl;
      return 628788;
   }
   
   //The Hamiltonian
   CheMPS2::Hamiltonian * theHam = new CheMPS2::Hamiltonian(matrixelements);
   CheMPS2::Irreps SymmInfo(theHam->getNGroup());
   
   //The maximum linear size for rotations
   int * Isizes = new int[SymmInfo.getNumberOfIrreps()];
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){ Isizes[irrep] = 0; }
   for (int orb=0; orb<theHam->getL(); orb++){ Isizes[theHam->getOrbitalIrrep(orb)]++; }
   int maxlinsize = 0;
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      if (Isizes[irrep] > maxlinsize){ maxlinsize = Isizes[irrep]; }
   }
   
   //Work memory for the Edmiston-Ruedenberg localization procedure
   double * workmem1 = new double[maxlinsize*maxlinsize*maxlinsize*maxlinsize];
   double * workmem2 = new double[maxlinsize*maxlinsize*maxlinsize*maxlinsize];

   //The localizer
   const int printLevel = 2; //Print a little extra information
   CheMPS2::EdmistonRuedenberg theLocalizer(theHam, printLevel);
   const double Cost = theLocalizer.Optimize(workmem1, workmem2);
   
   //Deallocate stuff
   delete theHam;
   delete [] Isizes;
   delete [] workmem1;
   delete [] workmem2;
   
   //Check succes
   bool success = (fabs(Cost - 20.6290770115892) < 1e-10) ? true : false;
   cout << "================> Did test 9 succeed : ";
   if (success){ cout << "yes" << endl; }
   else { cout << "no" << endl; }

   return 0;

}


