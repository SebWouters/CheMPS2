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
#include <sstream>

#include "Initialize.h"
#include "Hamiltonian.h"

using namespace std;

int main(void){

   CheMPS2::Initialize::Init();
   
   //The path to the matrix elements
   string matrixelements = "../../tests/matrixelements/O2_CCPVDZ.dat";
   struct stat stFileInfo;
   int intStat = stat(matrixelements.c_str(),&stFileInfo);
   if (intStat != 0){
      cout << "Please set the correct relative path to tests/matrixelements/O2_CCPVDZ.dat in tests/test7.cpp for the compiled binary test7 to work." << endl;
      return 628788;
   }
   
   //The Hamiltonian
   CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(matrixelements);
   Ham->save();
   int L = Ham->getL();
   int GroupNumber = Ham->getNGroup();
   int * OrbitalIrreps = new int[L];
   for (int cnt=0; cnt<L; cnt++){
      OrbitalIrreps[cnt] = Ham->getOrbitalIrrep(cnt);
   }
   
   //The loaded Hamiltonian
   CheMPS2::Hamiltonian * HamLoaded = new CheMPS2::Hamiltonian(L, GroupNumber, OrbitalIrreps);
   HamLoaded->read();
   
   //Test whether it's the same thing
   double RMS = 0.0;
   double temp = 0.0;
   
   temp = Ham->getEconst() - HamLoaded->getEconst();
   RMS += temp * temp;
   for (int i1=0; i1<L; i1++){
      for (int i2=0; i2<L; i2++){
         temp = Ham->getTmat(i1,i2) - HamLoaded->getTmat(i1,i2);
         RMS += temp * temp;
         for (int i3=0; i3<L; i3++){
            for (int i4=0; i4<L; i4++){
               temp = Ham->getVmat(i1,i2,i3,i4) - HamLoaded->getVmat(i1,i2,i3,i4);
               RMS += temp * temp;
            }
         }
      }
   }
   RMS = sqrt(RMS);
   
   cout << "The RMS difference between Ham and HamLoaded is " << RMS << endl;
   
   //Clean up
   delete [] OrbitalIrreps;
   delete Ham;
   delete HamLoaded;
   stringstream temp1;
   temp1 << "rm " << CheMPS2::HAMILTONIAN_TmatStorageName;
   int info = system(temp1.str().c_str());
   stringstream temp2;
   temp2 << "rm " << CheMPS2::HAMILTONIAN_VmatStorageName;
   info = system(temp2.str().c_str());
   stringstream temp3;
   temp3 << "rm " << CheMPS2::HAMILTONIAN_ParentStorageName;
   info = system(temp3.str().c_str());
   
   //Check succes
   bool success = (RMS < 1e-10) ? true : false;
   cout << "================> Did test 7 succeed : ";
   if (success){ cout << "yes" << endl; }
   else { cout << "no" << endl; }

   return 0;

}


