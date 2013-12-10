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

/*! \mainpage
 * \verbinclude ../README.md
 */

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include <sys/stat.h>

#include "DMRG.h"

using std::cout;
using std::endl;

CheMPS2::DMRG::DMRG(Problem * Probin, ConvergenceScheme * OptSchemeIn){

   PrintLicense();

   if (Probin->checkConsistency()){
      Prob = Probin;
   } else {
      cout << "DMRG::DMRG : The Problem consistency check failed." << endl;
      Prob = NULL;
   }

   OptScheme = OptSchemeIn;
   RNstorage = rand();
   nStates = 1;

   Ltensors = new TensorL ** [Prob->gL()-1];
   F0tensors = new TensorF0 *** [Prob->gL()-1];
   F1tensors = new TensorF1 *** [Prob->gL()-1];
   S0tensors = new TensorS0 *** [Prob->gL()-1];
   S1tensors = new TensorS1 *** [Prob->gL()-1];
   Atensors = new TensorA *** [Prob->gL()-1];
   Btensors = new TensorB *** [Prob->gL()-1];
   Ctensors = new TensorC *** [Prob->gL()-1];
   Dtensors = new TensorD *** [Prob->gL()-1];
   Qtensors = new TensorQ ** [Prob->gL()-1];
   Xtensors = new TensorX * [Prob->gL()-1];
   isAllocated = new int[Prob->gL()-1]; //0 not allocated, 1 allocated with movingRight true, 2 allocated with movingRight false
   
   for (int cnt=0; cnt<Prob->gL()-1; cnt++){ isAllocated[cnt] = 0; }
   
   the2DMallocated = false;
   Exc_activated = false;
   
   setupBookkeeperAndMPS();
   PreSolve();

}

void CheMPS2::DMRG::setupBookkeeperAndMPS(){

   std::stringstream sstream;
   sstream << "CheMPS2_MPS" << nStates-1 << ".h5";
   MPSstoragename.assign( sstream.str() );
   
   denBK = new SyBookkeeper(Prob,OptScheme->getD(0));
   if (!(denBK->IsPossible())){
      delete denBK;
      cout << "DMRG::DMRG : The desired symmetry is not possible." << endl;
      denBK = NULL; //Now all the rest will fail too.
   }
   
   struct stat stFileInfo;
   int intStat = stat(MPSstoragename.c_str(),&stFileInfo);
   loadedMPS = ((CheMPS2::DMRG_storeMpsOnDisk) && (intStat==0))? true : false ;
   
   if (loadedMPS){ loadDIM(MPSstoragename,denBK); }
   
   MPS = new TensorT * [Prob->gL()];
   for (int cnt=0; cnt<Prob->gL(); cnt++){ MPS[cnt] = new TensorT(cnt,denBK->gIrrep(cnt),denBK); }
   
   //Left normalize them all and calculate all diagrams up to the end
   if (loadedMPS){
      bool isConverged;
      loadMPS(MPSstoragename, MPS, &isConverged);
      cout << "Loaded MPS " << MPSstoragename << " converged y/n? : " << isConverged << endl;
   } else {
      for (int cnt=0; cnt<Prob->gL(); cnt++){
         TensorDiag * Dstor = new TensorDiag(cnt+1,denBK);
         MPS[cnt]->random();
         MPS[cnt]->QR(Dstor);
         delete Dstor;
      }
   }

}

CheMPS2::DMRG::~DMRG(){

   if (denBK!=NULL) delete denBK;
   
   deleteAllBoundaryOperators();
   
   delete [] Ltensors;
   delete [] F0tensors;
   delete [] F1tensors;
   delete [] S0tensors;
   delete [] S1tensors;
   delete [] Atensors;
   delete [] Btensors;
   delete [] Ctensors;
   delete [] Dtensors;
   delete [] Qtensors;
   delete [] Xtensors;
   delete [] isAllocated;
   
   for (int cnt=0; cnt<Prob->gL(); cnt++){ delete MPS[cnt]; }
   delete [] MPS;
   
   if (Exc_activated){
      delete [] Exc_Eshifts;
      for (int cnt=0; cnt<nStates-1; cnt++){
         for (int cnt2=0; cnt2<Prob->gL(); cnt2++){ delete Exc_MPSs[cnt][cnt2]; }
         delete [] Exc_MPSs[cnt];
      }
      delete [] Exc_MPSs;
      for (int cnt=0; cnt<nStates-1; cnt++){ delete Exc_BKs[cnt]; }
      delete [] Exc_BKs;
      for (int cnt=0; cnt<nStates-1; cnt++){ delete [] Exc_Overlaps[cnt]; } //The rest is allocated and deleted at DMRGoperators.cpp
      delete [] Exc_Overlaps;
   }
   
   if (the2DMallocated){ delete the2DM; }

}

void CheMPS2::DMRG::PreSolve(){
   
   for (int cnt=0; cnt<Prob->gL()-2; cnt++){ updateMovingRightSafeFirstTime(cnt); }
   
   MinEnergy = 1e8;
   MaxDiscWeightLastSweep = 0.0;

}

double CheMPS2::DMRG::Solve(){

   bool change = (MinEnergy<1e8) ? true : false; //1 sweep from right to left: fixed virtual dimensions
   
   for (int instruction=0; instruction < OptScheme->getNInstructions(); instruction++){
   
      double Energy = 0.0;
      double EnergyPrevious = 1.0;
   
      int nIterations = 0;
      
      while ( (fabs(Energy-EnergyPrevious) > OptScheme->getEconv(instruction) ) && ( nIterations < OptScheme->getMaxSweeps(instruction) )){
      
         EnergyPrevious = Energy;
         Energy = sweepleft(change, instruction);
         cout << "***  The max. disc. weight at last sweep is " << MaxDiscWeightLastSweep << endl;
         if (!change) change = true; //rest of sweeps: variable virtual dimensions
         Energy = sweepright(change, instruction);
         cout << "***  The max. disc. weight at last sweep is " << MaxDiscWeightLastSweep << endl;
         if (CheMPS2::DMRG_storeMpsOnDisk){ saveMPS(MPSstoragename, MPS, denBK, false); }
         
         nIterations++;
         
         cout << "*** Number of leftright sweep iterations is " << nIterations << endl; 
         cout << "***                The energy difference is " << fabs(Energy-EnergyPrevious) << endl;
         cout << "***                           The energy is " << Energy << endl;
         if (Exc_activated){ calcOverlapsWithLowerStates(); }
      
      }
      
      cout <<    "****************************************************************************" << endl;
      cout <<    "***    Performed instruction " << instruction << endl;
      cout <<    "***    Number of reduced DMRG basis states D = " << OptScheme->getD(instruction) << endl;
      cout <<    "***    The min. energy during the sweeps is " << MinEnergy << endl;
      cout <<    "***    The max. discarded weight during the last sweep is " << MaxDiscWeightLastSweep << endl;
      cout <<    "****************************************************************************" << endl;
   
   }
   
   return MinEnergy;

}

double CheMPS2::DMRG::sweepleft(const bool change, const int instruction){

   double Energy = 0.0;
   double NoiseLevel = OptScheme->getNoisePrefactor(instruction) * MaxDiscWeightLastSweep;
   MaxDiscWeightLastSweep = 0.0;

   for (int index = Prob->gL()-2; index>0; index--){
      //Construct S
      Sobject * denS = new Sobject(index,denBK->gIrrep(index),denBK->gIrrep(index+1),denBK);
      denS->Join(MPS[index],MPS[index+1]);
      
      //Feed everything to the solver
      Heff Solver(denBK, Prob);
      double ** VeffTilde = NULL;
      if (Exc_activated){
         VeffTilde = new double*[nStates-1];
         for (int cnt=0; cnt<nStates-1; cnt++){
            VeffTilde[cnt] = new double[denS->gKappa2index(denS->gNKappa())];
            calcVeffTilde(VeffTilde[cnt], denS, cnt);
         }
      }
      Energy = Solver.SolveDAVIDSON(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nStates-1, VeffTilde);
      if (Exc_activated){
         //calcOverlapsWithLowerStatesDuringSweeps_debug(VeffTilde, denS);
         for (int cnt=0; cnt<nStates-1; cnt++){ delete [] VeffTilde[cnt]; }
         delete [] VeffTilde;
      }
      Energy += Prob->gEconst();
      if (Energy<MinEnergy){ MinEnergy = Energy; }
      
      //Decompose the S-object
      if (NoiseLevel>0.0){ denS->addNoise(NoiseLevel); }
      double discWeight = denS->Split(MPS[index],MPS[index+1],OptScheme->getD(instruction),false,change);
      delete denS;
      if (discWeight > MaxDiscWeightLastSweep){ MaxDiscWeightLastSweep = discWeight; }
      
      //Print info
      cout << "Energy at sites (" << index << ", " << (index+1) << ") is " << Energy << endl;
      if (CheMPS2::DMRG_printDiscardedWeight && change){ cout << "   Info(DMRG) : Discarded weight in SVD decomp. (non-reduced) = " << discWeight << endl; }
      
      //Prepare for next step
      updateMovingLeftSafe(index);

   }
   
   return Energy;

}

double CheMPS2::DMRG::sweepright(const bool change, const int instruction){

   double Energy=0.0;
   double NoiseLevel = OptScheme->getNoisePrefactor(instruction) * MaxDiscWeightLastSweep;
   MaxDiscWeightLastSweep = 0.0;

   for (int index = 0; index<Prob->gL()-2; index++){
      //Construct S
      Sobject * denS = new Sobject(index,denBK->gIrrep(index),denBK->gIrrep(index+1),denBK);
      denS->Join(MPS[index],MPS[index+1]);
      
      //Feed everything to solver
      Heff Solver(denBK, Prob);
      double ** VeffTilde = NULL;
      if (Exc_activated){
         VeffTilde = new double*[nStates-1];
         for (int cnt=0; cnt<nStates-1; cnt++){
            VeffTilde[cnt] = new double[denS->gKappa2index(denS->gNKappa())];
            calcVeffTilde(VeffTilde[cnt], denS, cnt);
         }
      }
      Energy = Solver.SolveDAVIDSON(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nStates-1, VeffTilde);
      if (Exc_activated){
         //calcOverlapsWithLowerStatesDuringSweeps_debug(VeffTilde, denS);
         for (int cnt=0; cnt<nStates-1; cnt++){ delete [] VeffTilde[cnt]; }
         delete [] VeffTilde;
      }
      Energy += Prob->gEconst();
      if (Energy<MinEnergy){ MinEnergy = Energy; }
      
      //Decompose the S-object
      if (NoiseLevel>0.0){ denS->addNoise(NoiseLevel); }
      double discWeight = denS->Split(MPS[index],MPS[index+1],OptScheme->getD(instruction),true,change);
      delete denS;
      if (discWeight > MaxDiscWeightLastSweep){ MaxDiscWeightLastSweep = discWeight; }
      
      //Print info
      cout << "Energy at sites (" << index << ", " << (index+1) << ") is " << Energy << endl;
      if (CheMPS2::DMRG_printDiscardedWeight && change){ cout << "   Info(DMRG) : Discarded weight in SVD decomp. (non-reduced) = " << discWeight << endl; }
      
      //Prepare for next step
      updateMovingRightSafe(index);

   }
   
   return Energy;

}

void CheMPS2::DMRG::activateExcitations(const int maxExcIn){

   Exc_activated = true;
   maxExc = maxExcIn;
   Exc_Eshifts = new double[maxExc];
   Exc_MPSs = new TensorT ** [maxExc];
   Exc_BKs = new SyBookkeeper * [maxExc];
   Exc_Overlaps = new TensorO ** [maxExc];

}

void CheMPS2::DMRG::newExcitation(const double EshiftIn){

   if (!Exc_activated){
   
      cout << "DMRG::activateExcitations has not been called yet!" << endl;
      return;
      
   }

   if (nStates-1<maxExc){
   
      Exc_Eshifts[nStates-1] = EshiftIn;
      Exc_MPSs[nStates-1] = MPS;
      Exc_BKs[nStates-1] = denBK;
      Exc_Overlaps[nStates-1] = new TensorO * [Prob->gL()-1];
      
      deleteAllBoundaryOperators();
      
      nStates++;
      
      setupBookkeeperAndMPS();
      PreSolve();
   
   } else {
   
      cout << "Max. number of excitations has been reached!" << endl;
      
   }

}

