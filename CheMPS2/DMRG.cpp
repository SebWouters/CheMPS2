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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <assert.h>
#include <unistd.h>

#include "DMRG.h"

using std::cout;
using std::cerr;
using std::endl;

CheMPS2::DMRG::DMRG(Problem * ProbIn, ConvergenceScheme * OptSchemeIn, const bool makechkpt, const string tmpfolder){

   PrintLicense();
   
   assert( ProbIn->checkConsistency() );
   Prob = ProbIn;
   L = Prob->gL();
   OptScheme = OptSchemeIn;
   thePID = getpid();
   nStates = 1;

   Ltensors = new TensorL ** [L-1];
   F0tensors = new TensorF0 *** [L-1];
   F1tensors = new TensorF1 *** [L-1];
   S0tensors = new TensorS0 *** [L-1];
   S1tensors = new TensorS1 *** [L-1];
   Atensors = new TensorA *** [L-1];
   Btensors = new TensorB *** [L-1];
   Ctensors = new TensorC *** [L-1];
   Dtensors = new TensorD *** [L-1];
   Qtensors = new TensorQ ** [L-1];
   Xtensors = new TensorX * [L-1];
   isAllocated = new int[L-1]; //0 not allocated, 1 allocated with movingRight true, 2 allocated with movingRight false
   
   for (int cnt=0; cnt<L-1; cnt++){ isAllocated[cnt] = 0; }
   
   the2DMallocated = false;
   the2DM = NULL;
   theCorrAllocated = false;
   theCorr = NULL;
   Exc_activated = false;
   makecheckpoints = makechkpt;
   tempfolder = tmpfolder;
   
   setupBookkeeperAndMPS();
   PreSolve();

}

void CheMPS2::DMRG::setupBookkeeperAndMPS(){

   std::stringstream sstream;
   sstream << CheMPS2::DMRG_MPS_storage_prefix << nStates-1 << ".h5";
   MPSstoragename.assign( sstream.str() );
   
   denBK = new SyBookkeeper(Prob, OptScheme->getD(0));
   assert( denBK->IsPossible() );
   
   struct stat stFileInfo;
   int intStat = stat(MPSstoragename.c_str(),&stFileInfo);
   loadedMPS = ((makecheckpoints) && (intStat==0))? true : false ;
   
   if (loadedMPS){ loadDIM(MPSstoragename,denBK); }
   
   MPS = new TensorT * [L];
   for (int cnt=0; cnt<L; cnt++){ MPS[cnt] = new TensorT(cnt,denBK->gIrrep(cnt),denBK); }
   
   //Left normalize them all and calculate all diagrams up to the end
   if (loadedMPS){
      bool isConverged;
      loadMPS(MPSstoragename, MPS, &isConverged);
      cout << "Loaded MPS " << MPSstoragename << " converged y/n? : " << isConverged << endl;
   } else {
      for (int cnt=0; cnt<L; cnt++){
         TensorDiag * Dstor = new TensorDiag(cnt+1,denBK);
         MPS[cnt]->random();
         MPS[cnt]->QR(Dstor);
         delete Dstor;
      }
   }

}

CheMPS2::DMRG::~DMRG(){

   delete denBK;
   
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
   
   for (int cnt = 0; cnt < L; cnt++){ delete MPS[cnt]; }
   delete [] MPS;
   
   if (Exc_activated){
      delete [] Exc_Eshifts;
      for (int state = 0; state < nStates-1; state++){
         for (int orb = 0; orb < L; orb++){ delete Exc_MPSs[ state ][ orb ]; }
         delete [] Exc_MPSs[ state ];
      }
      delete [] Exc_MPSs;
      for (int cnt=0; cnt<nStates-1; cnt++){ delete Exc_BKs[cnt]; }
      delete [] Exc_BKs;
      for (int cnt=0; cnt<nStates-1; cnt++){ delete [] Exc_Overlaps[cnt]; } //The rest is allocated and deleted at DMRGoperators.cpp
      delete [] Exc_Overlaps;
   }
   
   if (the2DMallocated){ delete the2DM; }
   if (theCorrAllocated){ delete theCorr; }

}

void CheMPS2::DMRG::PreSolve(){
   
   for (int cnt=0; cnt<L-2; cnt++){ updateMovingRightSafeFirstTime(cnt); }
   
   TotalMinEnergy = 1e8;
   MaxDiscWeightLastSweep = 0.0;

}

double CheMPS2::DMRG::Solve(){

   bool change = (TotalMinEnergy<1e8) ? true : false; //1 sweep from right to left: fixed virtual dimensions
   
   double Energy = 0.0;
   
   for (int instruction=0; instruction < OptScheme->getNInstructions(); instruction++){
   
      int nIterations = 0;
      double EnergyPrevious = Energy + 10 * OptScheme->getEconv(instruction); //Guarantees that there's always at least 1 left-right sweep
      
      while ( (fabs(Energy-EnergyPrevious) > OptScheme->getEconv(instruction) ) && ( nIterations < OptScheme->getMaxSweeps(instruction) )){
      
         struct timeval start, end;
         EnergyPrevious = Energy;
         gettimeofday(&start, NULL);
         Energy = sweepleft(change, instruction);
         gettimeofday(&end, NULL);
         double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
         cout << "***  Information on left sweep " << nIterations << " of instruction " << instruction << ":" << endl;
         cout << "***     Elapsed wall time        = " << elapsed << " seconds" << endl;
         cout << "***     Minimum energy           = " << LastMinEnergy << endl;
         cout << "***     Maximum discarded weight = " << MaxDiscWeightLastSweep << endl;
         if (!change) change = true; //rest of sweeps: variable virtual dimensions
         gettimeofday(&start, NULL);
         Energy = sweepright(change, instruction);
         gettimeofday(&end, NULL);
         elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
         cout << "***  Information on right sweep " << nIterations << " of instruction " << instruction << ":" << endl;
         cout << "***     Elapsed wall time        = " << elapsed << " seconds" << endl;
         cout << "***     Minimum energy           = " << LastMinEnergy << endl;
         cout << "***     Maximum discarded weight = " << MaxDiscWeightLastSweep << endl;
         if (makecheckpoints){ saveMPS(MPSstoragename, MPS, denBK, false); }
         
         nIterations++;
         
         cout << "***  Energy difference with respect to previous leftright sweep = " << fabs(Energy-EnergyPrevious) << endl;
         if (Exc_activated){ calcOverlapsWithLowerStates(); }
      
      }
      
      cout <<    "****************************************************************************" << endl;
      cout <<    "***  Information on completed instruction " << instruction << ":" << endl;
      cout <<    "***     The reduced virtual dimension DSU(2)               = " << OptScheme->getD(instruction) << endl;
      cout <<    "***     Minimum energy encountered during all instructions = " << TotalMinEnergy << endl;
      cout <<    "***     Minimum energy encountered during the last sweep   = " << LastMinEnergy << endl;
      cout <<    "***     Maximum discarded weight during the last sweep     = " << MaxDiscWeightLastSweep << endl;
      cout <<    "****************************************************************************" << endl;
   
   }
   
   return TotalMinEnergy;

}

double CheMPS2::DMRG::sweepleft(const bool change, const int instruction){

   double Energy = 0.0;
   double NoiseLevel = OptScheme->getNoisePrefactor(instruction) * MaxDiscWeightLastSweep;
   MaxDiscWeightLastSweep = 0.0;
   LastMinEnergy = 1e8;

   for (int index = L-2; index>0; index--){
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
      if (Energy<TotalMinEnergy){ TotalMinEnergy = Energy; }
      if (Energy<LastMinEnergy){  LastMinEnergy  = Energy; }
      
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
   LastMinEnergy = 1e8;

   for (int index = 0; index<L-2; index++){
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
      if (Energy<TotalMinEnergy){ TotalMinEnergy = Energy; }
      if (Energy<LastMinEnergy){  LastMinEnergy  = Energy; }
      
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
   
      cerr << "DMRG::newExcitation : DMRG::activateExcitations has not been called yet!" << endl;
      return;
      
   }

   if (nStates-1<maxExc){
   
      Exc_Eshifts[nStates-1] = EshiftIn;
      Exc_MPSs[nStates-1] = MPS;
      Exc_BKs[nStates-1] = denBK;
      Exc_Overlaps[nStates-1] = new TensorO * [L-1];
      
      if (the2DMallocated){
         delete the2DM;
         the2DMallocated = false;
      }
      if (theCorrAllocated){
         delete theCorr;
         theCorrAllocated = false;
      }
      deleteAllBoundaryOperators();
      
      nStates++;
      
      setupBookkeeperAndMPS();
      PreSolve();
   
   } else {
   
      cerr << "DMRG::newExcitation : Maximum number of excitations has been reached!" << endl;
      
   }

}

