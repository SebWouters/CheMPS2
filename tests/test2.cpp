/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013 Sebastian Wouters

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

#include "DMRG.h"

using namespace std;

int main(void){

   cout.precision(15);
   srand(time(NULL));
   
   //The path to the matrix elements
   string matrixelements = "../../tests/matrixelements/H6_N6_S0_d2h_I0.dat";
   struct stat stFileInfo;
   int intStat = stat(matrixelements.c_str(),&stFileInfo);
   if (intStat != 0){
      cout << "Please set the correct relative path to tests/matrixelements/H6_N6_S0_d2h_I0.dat in tests/test2.cpp for the compiled binary test2 to work." << endl;
      return 628788;
   }
   
   //The Hamiltonian
   CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(matrixelements);
   cout << "The group was found to be " << CheMPS2::Irreps::getGroupName(Ham->getNGroup()) << endl;
   
   //The targeted state
   int TwoS = 0;
   int N = 6;
   int Irrep = 0;
   CheMPS2::Problem * Prob = new CheMPS2::Problem(Ham, TwoS, N, Irrep);
   Prob->SetupReorderD2h();
   
   //The convergence scheme
   CheMPS2::ConvergenceScheme * OptScheme = new CheMPS2::ConvergenceScheme(2);
   int D = 30;
   double Econv = 1e-10;
   int maxSweeps = 3;
   double noisePrefactor = 0.1;
   OptScheme->setInstruction(0,D,Econv,maxSweeps,noisePrefactor);
   D = 1000;
   maxSweeps = 10;
   noisePrefactor = 0.0;
   OptScheme->setInstruction(1,D,Econv,maxSweeps,noisePrefactor);
   
   //Run ground state calculation
   CheMPS2::DMRG * theDMRG = new CheMPS2::DMRG(Prob,OptScheme);
   double Energy = theDMRG->Solve();
   theDMRG->calc2DM();
   
   //Clean up
   if (CheMPS2::DMRG_storeMpsOnDisk){ theDMRG->deleteStoredMPS(); }
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
   delete theDMRG;
   delete OptScheme;
   delete Prob;
   delete Ham;
   
   //Check succes
   bool success = (fabs(Energy + 3.33351730146068) < 1e-10) ? true : false;
   cout << "================> Did test 2 succeed : ";
   if (success){ cout << "yes" << endl; }
   else { cout << "no" << endl; }

   return 0;

}


