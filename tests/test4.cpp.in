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

#include "Initialize.h"
#include "DMRG.h"

using namespace std;

int main(void){

   CheMPS2::Initialize::Init();
   
   //The Hamiltonian: 1D Hubbard model
   int L = 10;
   int Group = 0;
   double U = 2.0;
   double T = -1.0;
   int * irreps = new int[L];
   for (int cnt=0; cnt<L; cnt++){ irreps[cnt] = 0; }
   CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(L,Group,irreps);
   Ham->setEconst(0.0);
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         Ham->setTmat(cnt1,cnt2,0.0);
         for (int cnt3=0; cnt3<L; cnt3++){
            for (int cnt4=0; cnt4<L; cnt4++){
               Ham->setVmat(cnt1,cnt2,cnt3,cnt4,0.0);
            }
         }
      }
   }
   for (int cnt=0; cnt<L; cnt++){ Ham->setVmat(cnt,cnt,cnt,cnt,U); }
   for (int cnt=0; cnt<L-1; cnt++){ Ham->setTmat(cnt,cnt+1,T); }
   
   //The targeted state
   int TwoS = 5;
   int N = 9;
   int Irrep = 0;
   CheMPS2::Problem * Prob = new CheMPS2::Problem(Ham, TwoS, N, Irrep);
   
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
   theDMRG->calc2DMandCorrelations();
   theDMRG->getCorrelations()->Print();
   
   //Clean up
   if (CheMPS2::DMRG_storeMpsOnDisk){ theDMRG->deleteStoredMPS(); }
   if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
   delete theDMRG;
   delete OptScheme;
   delete Prob;
   delete Ham;
   delete [] irreps;
   
   //Check succes
   bool success = (fabs(Energy + 6.24546701457254) < 1e-10) ? true : false;
   cout << "================> Did test 4 succeed : ";
   if (success){ cout << "yes" << endl; }
   else { cout << "no" << endl; }

   return 0;

}


