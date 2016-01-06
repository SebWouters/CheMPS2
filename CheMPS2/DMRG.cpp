/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

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
#include "MPIchemps2.h"

using std::cout;
using std::cerr;
using std::endl;

CheMPS2::DMRG::DMRG(Problem * ProbIn, ConvergenceScheme * OptSchemeIn, const bool makechkpt, const string tmpfolder){

   #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER ){ PrintLicense(); }
   #else
      PrintLicense();
   #endif
   
   assert( ProbIn->checkConsistency() );
   Prob = ProbIn;
   L = Prob->gL();
   Prob->construct_mxelem();
   OptScheme = OptSchemeIn;
   thePID = getpid(); //PID is unique for each MPI process
   nStates = 1;

   Ltensors = new TensorL ** [L-1];
   F0tensors = new TensorF0 *** [L-1];
   F1tensors = new TensorF1 *** [L-1];
   S0tensors = new TensorS0 *** [L-1];
   S1tensors = new TensorS1 *** [L-1];
   Atensors = new TensorOperator *** [L-1];
   Btensors = new TensorOperator *** [L-1];
   Ctensors = new TensorOperator *** [L-1];
   Dtensors = new TensorOperator *** [L-1];
   Qtensors = new TensorQ ** [L-1];
   Xtensors = new TensorX * [L-1];
   isAllocated = new int[L-1]; //0 not allocated, 1 allocated with movingRight true, 2 allocated with movingRight false
   
   tensor_3rdm_a_J0_doublet = NULL;
   tensor_3rdm_a_J1_doublet = NULL;
   tensor_3rdm_a_J1_quartet = NULL;
   tensor_3rdm_b_J0_doublet = NULL;
   tensor_3rdm_b_J1_doublet = NULL;
   tensor_3rdm_b_J1_quartet = NULL;
   tensor_3rdm_c_J0_doublet = NULL;
   tensor_3rdm_c_J1_doublet = NULL;
   tensor_3rdm_c_J1_quartet = NULL;
   tensor_3rdm_d_J0_doublet = NULL;
   tensor_3rdm_d_J1_doublet = NULL;
   tensor_3rdm_d_J1_quartet = NULL;
   
   Gtensors = NULL;
   Ytensors = NULL;
   Ztensors = NULL;
   Ktensors = NULL;
   Mtensors = NULL;
   
   for (int cnt=0; cnt<L-1; cnt++){ isAllocated[cnt] = 0; }
   for (int timecnt=0; timecnt<CHEMPS2_TIME_VECLENGTH; timecnt++){ timings[timecnt]=0.0; } // Clear here so that valgrind can never complain :-)
   num_double_write_disk = 0;
   num_double_read_disk  = 0;
   
   the2DM  = NULL;
   the3DM  = NULL;
   theCorr = NULL;
   Exc_activated = false;
   makecheckpoints = makechkpt;
   tempfolder = tmpfolder;
   
   setupBookkeeperAndMPS();
   PreSolve();

}

void CheMPS2::DMRG::setupBookkeeperAndMPS(){
   
   denBK = new SyBookkeeper(Prob, OptScheme->getD(0));
   assert( denBK->IsPossible() );
   
   std::stringstream sstream;
   sstream << CheMPS2::DMRG_MPS_storage_prefix << nStates-1 << ".h5";
   MPSstoragename.assign( sstream.str() );
   struct stat stFileInfo;
   int intStat = stat(MPSstoragename.c_str(),&stFileInfo);
   loadedMPS = ((makecheckpoints) && (intStat==0))? true : false ;
   #ifdef CHEMPS2_MPI_COMPILATION
   assert( MPIchemps2::all_booleans_equal( loadedMPS ) );
   #endif
   
   if (loadedMPS){ loadDIM(MPSstoragename,denBK); }
   
   MPS = new TensorT * [L];
   for (int cnt=0; cnt<L; cnt++){ MPS[cnt] = new TensorT(cnt,denBK->gIrrep(cnt),denBK); }
   
   if (loadedMPS){
      bool isConverged;
      loadMPS(MPSstoragename, MPS, &isConverged);
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER )
      #endif
      { cout << "Loaded MPS " << MPSstoragename << " converged y/n? : " << isConverged << endl; }
   } else {
      for (int cnt=0; cnt<L; cnt++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER ){
         #endif
            TensorOperator * diag = new TensorOperator(cnt+1, 0, 0, 0, true, true, false, denBK); // (J,N,I) = (0,0,0) and (moving_right, prime_last, jw_phase) = (true, true, false)
            MPS[cnt]->random();
            MPS[cnt]->QR(diag);
            delete diag;
         #ifdef CHEMPS2_MPI_COMPILATION
         }
         MPIchemps2::broadcast_tensor(MPS[cnt], MPI_CHEMPS2_MASTER);
         #endif
      }
   }

}

CheMPS2::DMRG::~DMRG(){

   if ( the2DM  != NULL ){ delete the2DM;  }
   if ( the3DM  != NULL ){ delete the3DM;  }
   if ( theCorr != NULL ){ delete theCorr; }
   
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
   
   for (int site=0; site<L; site++){ delete MPS[site]; }
   delete [] MPS;
   
   if (Exc_activated){
      delete [] Exc_Eshifts;
      for (int state = 0; state < nStates-1; state++){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIchemps2::mpi_rank() )
         #endif
         {
            for (int orb = 0; orb < L; orb++){ delete Exc_MPSs[ state ][ orb ]; }
            delete [] Exc_MPSs[ state ];
            delete Exc_BKs[ state ];
            delete [] Exc_Overlaps[ state ]; //The rest is allocated and deleted at DMRGoperators.cpp
         }
      }
      delete [] Exc_MPSs;
      delete [] Exc_BKs;
      delete [] Exc_Overlaps;
   }
   
   delete denBK;

}

void CheMPS2::DMRG::PreSolve(){
   
   for (int cnt=0; cnt<L-2; cnt++){ updateMovingRightSafeFirstTime(cnt); }
   
   TotalMinEnergy = 1e8;
   MaxDiscWeightLastSweep = 0.0;

}

double CheMPS2::DMRG::Solve(){

   bool change = (TotalMinEnergy<1e8) ? true : false; //1 sweep from right to left: fixed virtual dimensions
   
   double Energy = 0.0;
   
   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif
   
   for (int instruction=0; instruction < OptScheme->getNInstructions(); instruction++){
   
      int nIterations = 0;
      double EnergyPrevious = Energy + 10 * OptScheme->getEconv(instruction); //Guarantees that there's always at least 1 left-right sweep
      
      while ( (fabs(Energy-EnergyPrevious) > OptScheme->getEconv(instruction) ) && ( nIterations < OptScheme->getMaxSweeps(instruction) )){
      
         for (int timecnt=0; timecnt<CHEMPS2_TIME_VECLENGTH; timecnt++){ timings[timecnt]=0.0; }
         num_double_write_disk = 0;
         num_double_read_disk  = 0;
         struct timeval start, end;
         EnergyPrevious = Energy;
         gettimeofday(&start, NULL);
         Energy = sweepleft(change, instruction, am_i_master); // Only relevant call in this block of code
         gettimeofday(&end, NULL);
         double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
         if ( am_i_master ){
            cout << "******************************************************************" << endl;
            cout << "***  Information on left sweep " << nIterations << " of instruction " << instruction << ":" << endl;
            cout << "***     Elapsed wall time        = " << elapsed << " seconds" << endl;
            cout << "***       |--> S.join            = " << timings[ CHEMPS2_TIME_S_JOIN  ] << " seconds" << endl;
            cout << "***       |--> S.solve           = " << timings[ CHEMPS2_TIME_S_SOLVE ] << " seconds" << endl;
            cout << "***       |--> S.split           = " << timings[ CHEMPS2_TIME_S_SPLIT ] << " seconds" << endl;
            print_tensor_update_performance();
            cout << "***     Minimum energy           = " << LastMinEnergy << endl;
            cout << "***     Maximum discarded weight = " << MaxDiscWeightLastSweep << endl;
            cout << "******************************************************************" << endl;
         }
         if (!change) change = true; //rest of sweeps: variable virtual dimensions
         for (int timecnt=0; timecnt<CHEMPS2_TIME_VECLENGTH; timecnt++){ timings[timecnt]=0.0; }
         num_double_write_disk = 0;
         num_double_read_disk  = 0;
         gettimeofday(&start, NULL);
         Energy = sweepright(change, instruction, am_i_master); // Only relevant call in this block of code
         gettimeofday(&end, NULL);
         elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
         if ( am_i_master ){
            cout << "******************************************************************" << endl;
            cout << "***  Information on right sweep " << nIterations << " of instruction " << instruction << ":" << endl;
            cout << "***     Elapsed wall time        = " << elapsed << " seconds" << endl;
            cout << "***       |--> S.join            = " << timings[ CHEMPS2_TIME_S_JOIN  ] << " seconds" << endl;
            cout << "***       |--> S.solve           = " << timings[ CHEMPS2_TIME_S_SOLVE ] << " seconds" << endl;
            cout << "***       |--> S.split           = " << timings[ CHEMPS2_TIME_S_SPLIT ] << " seconds" << endl;
            print_tensor_update_performance();
            cout << "***     Minimum energy           = " << LastMinEnergy << endl;
            cout << "***     Maximum discarded weight = " << MaxDiscWeightLastSweep << endl;
            cout << "***     Energy difference with respect to previous leftright sweep = " << fabs(Energy-EnergyPrevious) << endl;
            cout << "******************************************************************" << endl;
            if ( makecheckpoints ){ saveMPS(MPSstoragename, MPS, denBK, false); } // Only the master proc makes MPS checkpoints !!
         }
         
         nIterations++;
         if (Exc_activated){ calcOverlapsWithLowerStates(); }
      
      }
      
      if ( am_i_master ){
         cout << "***  Information on completed instruction " << instruction << ":" << endl;
         cout << "***     The reduced virtual dimension DSU(2)               = " << OptScheme->getD(instruction) << endl;
         cout << "***     Minimum energy encountered during all instructions = " << TotalMinEnergy << endl;
         cout << "***     Minimum energy encountered during the last sweep   = " << LastMinEnergy << endl;
         cout << "***     Maximum discarded weight during the last sweep     = " << MaxDiscWeightLastSweep << endl;
         cout << "******************************************************************" << endl;
      }
   
   }
   
   return TotalMinEnergy;

}

double CheMPS2::DMRG::sweepleft(const bool change, const int instruction, const bool am_i_master){

   double Energy = 0.0;
   double NoiseLevel = OptScheme->getNoisePrefactor(instruction) * MaxDiscWeightLastSweep;
   MaxDiscWeightLastSweep = 0.0;
   LastMinEnergy = 1e8;
   struct timeval start, end;
   
   for (int index = L-2; index>0; index--){
      //Construct S
      gettimeofday(&start, NULL);
      Sobject * denS = new Sobject(index,denBK->gIrrep(index),denBK->gIrrep(index+1),denBK);
      //Each MPI process joins the MPS tensors. Before a matrix-vector multiplication the vector is broadcasted anyway.
      denS->Join(MPS[index],MPS[index+1]);
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_JOIN ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Feed everything to the solver
      gettimeofday(&start, NULL);
      Heff Solver(denBK, Prob);
      double ** VeffTilde = NULL;
      if (Exc_activated){ VeffTilde = prepare_excitations(denS); }
      //Each MPI process returns the correct energy. Only MPI_CHEMPS2_MASTER has the correct denS solution.
      Energy = Solver.SolveDAVIDSON(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nStates-1, VeffTilde);
      if (Exc_activated){ cleanup_excitations(VeffTilde); }
      Energy += Prob->gEconst();
      if (Energy<TotalMinEnergy){ TotalMinEnergy = Energy; }
      if (Energy<LastMinEnergy){  LastMinEnergy  = Energy; }
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_SOLVE ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Decompose the S-object
      gettimeofday(&start, NULL);
      if (( NoiseLevel>0.0 ) && ( am_i_master )){ denS->addNoise(NoiseLevel); }
      //MPI_CHEMPS2_MASTER decomposes denS. Each MPI process returns the correct discWeight and now has the new MPS tensors set.
      double discWeight = denS->Split(MPS[index],MPS[index+1],OptScheme->getD(instruction),false,change);
      delete denS;
      if (discWeight > MaxDiscWeightLastSweep){ MaxDiscWeightLastSweep = discWeight; }
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_SPLIT ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Print info
      if ( am_i_master ){
         cout << "Energy at sites (" << index << ", " << (index+1) << ") is " << Energy << endl;
         if (CheMPS2::DMRG_printDiscardedWeight && change){ cout << "   Info(DMRG) : Discarded weight in SVD decomp. (non-reduced) = " << discWeight << endl; }
      }
      
      //Prepare for next step
      gettimeofday(&start, NULL);
      updateMovingLeftSafe(index);
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_TENS_TOTAL ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

   }
   
   return Energy;

}

double CheMPS2::DMRG::sweepright(const bool change, const int instruction, const bool am_i_master){

   double Energy=0.0;
   double NoiseLevel = OptScheme->getNoisePrefactor(instruction) * MaxDiscWeightLastSweep;
   MaxDiscWeightLastSweep = 0.0;
   LastMinEnergy = 1e8;
   struct timeval start, end;
   
   for (int index = 0; index<L-2; index++){
      //Construct S
      gettimeofday(&start, NULL);
      Sobject * denS = new Sobject(index,denBK->gIrrep(index),denBK->gIrrep(index+1),denBK);
      //Each MPI process joins the MPS tensors. Before a matrix-vector multiplication the vector is broadcasted anyway.
      denS->Join(MPS[index],MPS[index+1]);
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_JOIN ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Feed everything to solver
      gettimeofday(&start, NULL);
      Heff Solver(denBK, Prob);
      double ** VeffTilde = NULL;
      if (Exc_activated){ VeffTilde = prepare_excitations(denS); }
      //Each MPI process returns the correct energy. Only MPI_CHEMPS2_MASTER has the correct denS solution.
      Energy = Solver.SolveDAVIDSON(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nStates-1, VeffTilde);
      if (Exc_activated){ cleanup_excitations(VeffTilde); }
      Energy += Prob->gEconst();
      if (Energy<TotalMinEnergy){ TotalMinEnergy = Energy; }
      if (Energy<LastMinEnergy){  LastMinEnergy  = Energy; }
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_SOLVE ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Decompose the S-object
      gettimeofday(&start, NULL);
      if (( NoiseLevel>0.0 ) && ( am_i_master )){ denS->addNoise(NoiseLevel); }
      //MPI_CHEMPS2_MASTER decomposes denS. Each MPI process returns the correct discWeight and now has the new MPS tensors set.
      double discWeight = denS->Split(MPS[index],MPS[index+1],OptScheme->getD(instruction),true,change);
      delete denS;
      if (discWeight > MaxDiscWeightLastSweep){ MaxDiscWeightLastSweep = discWeight; }
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_S_SPLIT ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
      
      //Print info
      if ( am_i_master ){
         cout << "Energy at sites (" << index << ", " << (index+1) << ") is " << Energy << endl;
         if (CheMPS2::DMRG_printDiscardedWeight && change){ cout << "   Info(DMRG) : Discarded weight in SVD decomp. (non-reduced) = " << discWeight << endl; }
      }
      
      //Prepare for next step
      gettimeofday(&start, NULL);
      updateMovingRightSafe(index);
      gettimeofday(&end, NULL);
      timings[ CHEMPS2_TIME_TENS_TOTAL ] += (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);

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

   assert( Exc_activated );
   assert( nStates-1 < maxExc );
   
   if ( the2DM  != NULL ){ delete the2DM;  the2DM  = NULL; }
   if ( the3DM  != NULL ){ delete the3DM;  the3DM  = NULL; }
   if ( theCorr != NULL ){ delete theCorr; theCorr = NULL; }
   deleteAllBoundaryOperators();

   Exc_Eshifts[nStates-1] = EshiftIn;
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::owner_specific_excitation( L, nStates-1 ) == MPIchemps2::mpi_rank() ){
   #endif
      Exc_MPSs[nStates-1] = MPS;
      Exc_BKs[nStates-1] = denBK;
      Exc_Overlaps[nStates-1] = new TensorO*[L-1];
   #ifdef CHEMPS2_MPI_COMPILATION
   } else {
      for (int site=0; site<L; site++){ delete MPS[site]; }
      delete [] MPS;
      delete denBK;
   }
   #endif
   
   nStates++;
   
   setupBookkeeperAndMPS();
   PreSolve();

}

