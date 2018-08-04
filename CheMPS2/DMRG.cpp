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
using std::endl;

CheMPS2::DMRG::DMRG( Problem * ProbIn, ConvergenceScheme * OptSchemeIn, const bool makechkpt, const string tmpfolder, int * occupancies ){

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

   Ltensors  = new TensorL ** [ L - 1 ];
   F0tensors = new TensorF0 *** [ L - 1 ];
   F1tensors = new TensorF1 *** [ L - 1 ];
   S0tensors = new TensorS0 *** [ L - 1 ];
   S1tensors = new TensorS1 *** [ L - 1 ];
   Atensors  = new TensorOperator *** [ L - 1 ];
   Btensors  = new TensorOperator *** [ L - 1 ];
   Ctensors  = new TensorOperator *** [ L - 1 ];
   Dtensors  = new TensorOperator *** [ L - 1 ];
   Qtensors  = new TensorQ ** [ L - 1 ];
   Xtensors  = new TensorX * [ L - 1 ];
   isAllocated = new int[ L - 1 ]; // 0 not allocated; 1 moving right; 2 moving left

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

   for ( int cnt = 0; cnt < L - 1; cnt++ ){ isAllocated[ cnt ] = 0; }
   for ( int timecnt = 0; timecnt < CHEMPS2_TIME_VECLENGTH; timecnt++ ){ timings[ timecnt ] = 0.0; }
   num_double_write_disk = 0;
   num_double_read_disk  = 0;
   
   the2DM  = NULL;
   the3DM  = NULL;
   theCorr = NULL;
   Exc_activated = false;
   makecheckpoints = makechkpt;
   tempfolder = tmpfolder;
   
   setupBookkeeperAndMPS( occupancies );
   PreSolve();

}

void CheMPS2::DMRG::setupBookkeeperAndMPS( int * occupancies ){

   denBK = new SyBookkeeper( Prob, OptScheme->get_D( 0 ) );
   assert( denBK->IsPossible() );

   std::stringstream sstream;
   sstream << CheMPS2::DMRG_MPS_storage_prefix << nStates-1 << ".h5";
   MPSstoragename.assign( sstream.str() );
   struct stat stFileInfo;
   int intStat = stat( MPSstoragename.c_str(), &stFileInfo );
   loadedMPS = (( makecheckpoints ) && ( intStat==0 )) ? true : false;
   #ifdef CHEMPS2_MPI_COMPILATION
   assert( MPIchemps2::all_booleans_equal( loadedMPS ) );
   #endif

   if ( loadedMPS ){ loadDIM( MPSstoragename, denBK ); }

   // Convert occupancies from HAM to DMRG orbitals
   if (( occupancies != NULL ) && ( Prob->gReorder() )){
      int * tmp_cpy_occ = new int[ L ];
      for ( int cnt = 0; cnt < L; cnt++ ){ tmp_cpy_occ[ cnt ] = occupancies[ cnt ]; }
      for ( int cnt = 0; cnt < L; cnt++ ){ occupancies[ cnt ] = tmp_cpy_occ[ Prob->gf2( cnt ) ]; }
      delete [] tmp_cpy_occ;
   }

   // Set to ROHF dimensions
   /*if (( !loadedMPS ) && ( occupancies != NULL )){
      int left_n  = 0;
      int left_i  = 0;
      int left_2s = 0;
      for ( int site = 0; site < denBK->gL(); site++ ){
         for ( int N = denBK->gNmin( site ); N <= denBK->gNmax( site ); N++ ){
            for ( int TwoS = denBK->gTwoSmin( site, N ); TwoS <= denBK->gTwoSmax( site, N ); TwoS+=2 ){
               for ( int Irrep = 0; Irrep < denBK->getNumberOfIrreps(); Irrep++ ){
                  denBK->SetDim( site, N, TwoS, Irrep, 0 );
               }
            }
         }
         denBK->SetDim( site, left_n, left_2s, left_i, 1 );
         left_n  = left_n + occupancies[ site ];
         left_i  = (( occupancies[ site ] == 1 ) ? Irreps::directProd( left_i, Prob->gIrrep( site ) ) : left_i );
         left_2s = (( occupancies[ site ] == 1 ) ? ( left_2s + 1 ) : left_2s );
      }
      assert( left_n  == Prob->gN() );
      assert( left_i  == Prob->gIrrep() );
      assert( left_2s == Prob->gTwoS() );
   }*/

   MPS = new TensorT * [ L ];
   for ( int cnt = 0; cnt < L; cnt++ ){ MPS[ cnt ] = new TensorT( cnt, denBK ); }

   if ( loadedMPS ){
      bool isConverged;
      loadMPS( MPSstoragename, MPS, &isConverged );
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER )
      #endif
      { cout << "Loaded MPS " << MPSstoragename << " converged y/n? : " << isConverged << endl; }
   } else {
      #ifdef CHEMPS2_MPI_COMPILATION
         const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
      #else
         const bool am_i_master = true;
      #endif
      if ( occupancies == NULL ){
         for ( int site = 0; site < L; site++ ){
            if ( am_i_master ){ MPS[ site ]->random(); }
            left_normalize( MPS[ site ], NULL );
         }
      } else {
         assert( Prob->check_rohf_occ( occupancies ) ); // Check compatibility
         int left_n  = 0;
         int left_i  = 0;
         int left_2s = 0;
         for ( int site = 0; site < L; site++ ){
            const int right_n  = left_n + occupancies[ site ];
            const int right_i  = (( occupancies[ site ] == 1 ) ? Irreps::directProd( left_i, Prob->gIrrep( site ) ) : left_i );
            const int right_2s = (( occupancies[ site ] == 1 ) ? ( left_2s + 1 ) : left_2s );
            const int dimL = denBK->gCurrentDim( site,      left_n,  left_2s,  left_i );
            const int dimR = denBK->gCurrentDim( site + 1, right_n, right_2s, right_i );
            assert( dimL > 0 );
            assert( dimR > 0 );
            if ( am_i_master ){
               MPS[ site ]->random();
               for ( int NL = right_n - 2; NL <= right_n; NL++ ){
                  const int DS = (( right_n == NL + 1 ) ? 1 : 0 );
                  const int IL = (( right_n == NL + 1 ) ? Irreps::directProd( right_i, Prob->gIrrep( site ) ) : right_i );
                  for ( int TwoSL = right_2s - DS; TwoSL <= right_2s + DS; TwoSL+=2 ){
                     const int dimL2 = denBK->gCurrentDim( site, NL, TwoSL, IL );
                     if ( dimL2 > 0 ){
                        double * space = MPS[ site ]->gStorage( NL, TwoSL, IL, right_n, right_2s, right_i );
                        for ( int row = 0; row < dimL2; row++ ){ space[ row + dimL2 * 0 ] = 0.0; }
                        if (( NL == left_n ) && ( TwoSL == left_2s ) && ( IL == left_i )){ space[ 0 + dimL * 0  ] = 42; }
                     }
                  }
               }
            }
            left_normalize( MPS[ site ], NULL );
            left_n  = right_n;
            left_i  = right_i;
            left_2s = right_2s;
         }
         assert( left_n  == Prob->gN() );
         assert( left_i  == Prob->gIrrep() );
         assert( left_2s == Prob->gTwoS() );
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

   for ( int site = 0; site < L; site++ ){ delete MPS[ site ]; }
   delete [] MPS;

   if ( Exc_activated ){
      delete [] Exc_Eshifts;
      for ( int state = 0; state < nStates - 1; state++ ){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_specific_excitation( L, state ) == MPIchemps2::mpi_rank() )
         #endif
         {
            for ( int orb = 0; orb < L; orb++ ){ delete Exc_MPSs[ state ][ orb ]; }
            delete [] Exc_MPSs[ state ];
            delete Exc_BKs[ state ];
            delete [] Exc_Overlaps[ state ]; // The rest is allocated and deleted at DMRGoperators.cpp
         }
      }
      delete [] Exc_MPSs;
      delete [] Exc_BKs;
      delete [] Exc_Overlaps;
   }

   delete denBK;

}

void CheMPS2::DMRG::PreSolve(){

   deleteAllBoundaryOperators();

   for ( int cnt = 0; cnt < L - 2; cnt++ ){ updateMovingRightSafeFirstTime( cnt ); }

   TotalMinEnergy = 1e8;
   MaxDiscWeightLastSweep = 0.0;

}

double CheMPS2::DMRG::Solve(){

   bool change = ( TotalMinEnergy < 1e8 ) ? true : false; // 1 sweep from right to left: fixed virtual dimensions

   double Energy = 0.0;

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   for ( int instruction = 0; instruction < OptScheme->get_number(); instruction++ ){

      int nIterations = 0;
      double EnergyPrevious = Energy + 10 * OptScheme->get_energy_conv( instruction ); // Guarantees that there's always at least 1 left-right sweep

      while (( fabs( Energy - EnergyPrevious ) > OptScheme->get_energy_conv( instruction ) ) && ( nIterations < OptScheme->get_max_sweeps( instruction ) )){

         for ( int timecnt = 0; timecnt < CHEMPS2_TIME_VECLENGTH; timecnt++ ){ timings[ timecnt ] = 0.0; }
         num_double_write_disk = 0;
         num_double_read_disk  = 0;
         struct timeval start, end;
         EnergyPrevious = Energy;
         gettimeofday( &start, NULL );
         Energy = sweepleft( change, instruction, am_i_master ); // Only relevant call in this block of code
         gettimeofday( &end, NULL );
         double elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
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
         }
         if ( Exc_activated ){ calc_overlaps( false ); }
         if ( am_i_master ){
            cout << "******************************************************************" << endl;
         }
         change = true; //rest of sweeps: variable virtual dimensions
         for ( int timecnt = 0; timecnt < CHEMPS2_TIME_VECLENGTH; timecnt++ ){ timings[ timecnt ] = 0.0; }
         num_double_write_disk = 0;
         num_double_read_disk  = 0;
         gettimeofday( &start, NULL );
         Energy = sweepright( change, instruction, am_i_master ); // Only relevant call in this block of code
         gettimeofday( &end, NULL );
         elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
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
         }
         if ( Exc_activated ){ calc_overlaps( true ); }
         if ( am_i_master ){
            cout << "******************************************************************" << endl;
            if ( makecheckpoints ){ saveMPS( MPSstoragename, MPS, denBK, false ); } // Only the master proc makes MPS checkpoints !!
         }

         nIterations++;

      }

      if ( am_i_master ){
         cout << "***  Information on completed instruction " << instruction << ":" << endl;
         cout << "***     The reduced virtual dimension DSU(2)               = " << OptScheme->get_D(instruction) << endl;
         cout << "***     The total number of reduced MPS variables          = " << get_num_mps_var() << endl;
         cout << "***     Minimum energy encountered during all instructions = " << TotalMinEnergy << endl;
         cout << "***     Minimum energy encountered during the last sweep   = " << LastMinEnergy << endl;
         cout << "***     Maximum discarded weight during the last sweep     = " << MaxDiscWeightLastSweep << endl;
         cout << "******************************************************************" << endl;
      }

   }

   return TotalMinEnergy;

}

double CheMPS2::DMRG::sweepleft( const bool change, const int instruction, const bool am_i_master ){

   double Energy = 0.0;
   const double noise_level = fabs( OptScheme->get_noise_prefactor( instruction ) ) * MaxDiscWeightLastSweep;
   const double dvdson_rtol = OptScheme->get_dvdson_rtol( instruction );
   const int vir_dimension  = OptScheme->get_D( instruction );
   MaxDiscWeightLastSweep = 0.0;
   LastMinEnergy = 1e8;

   for ( int index = L - 2; index > 0; index-- ){

      Energy = solve_site( index, dvdson_rtol, noise_level, vir_dimension, am_i_master, false, change );
      if ( Energy < TotalMinEnergy ){ TotalMinEnergy = Energy; }
      if ( Energy < LastMinEnergy  ){  LastMinEnergy = Energy; }
      if ( am_i_master ){
         cout << "Energy at sites (" << index << ", " << index + 1 << ") is " << Energy << endl;
      }

      // Prepare for next step
      struct timeval start, end;
      gettimeofday( &start, NULL );
      updateMovingLeftSafe( index );
      gettimeofday( &end, NULL );
      timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   }

   return Energy;

}

double CheMPS2::DMRG::sweepright( const bool change, const int instruction, const bool am_i_master ){

   double Energy = 0.0;
   const double noise_level = fabs( OptScheme->get_noise_prefactor( instruction ) ) * MaxDiscWeightLastSweep;
   const double dvdson_rtol = OptScheme->get_dvdson_rtol( instruction );
   const int vir_dimension  = OptScheme->get_D( instruction );
   MaxDiscWeightLastSweep = 0.0;
   LastMinEnergy = 1e8;

   for ( int index = 0; index < L - 2; index++ ){

      Energy = solve_site( index, dvdson_rtol, noise_level, vir_dimension, am_i_master, true, change );
      if ( Energy < TotalMinEnergy ){ TotalMinEnergy = Energy; }
      if ( Energy < LastMinEnergy  ){  LastMinEnergy = Energy; }
      if ( am_i_master ){
         cout << "Energy at sites (" << index << ", " << index + 1 << ") is " << Energy << endl;
      }

      // Prepare for next step
      struct timeval start, end;
      gettimeofday( &start, NULL );
      updateMovingRightSafe( index );
      gettimeofday( &end, NULL );
      timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   }

   return Energy;

}

double CheMPS2::DMRG::solve_site( const int index, const double dvdson_rtol, const double noise_level, const int virtual_dimension, const bool am_i_master, const bool moving_right, const bool change ){

   struct timeval start, end;

   // Construct two-site object S. Each MPI process joins the MPS tensors. Before a matrix-vector multiplication the vector is broadcasted anyway.
   gettimeofday( &start, NULL );
   Sobject * denS = new Sobject( index, denBK );
   denS->Join( MPS[ index ], MPS[ index + 1 ] );
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_TIME_S_JOIN ] += ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // Feed everything to the solver. Each MPI process returns the correct energy. Only MPI_CHEMPS2_MASTER has the correct denS solution.
   gettimeofday( &start, NULL );
   Heff Solver( denBK, Prob, dvdson_rtol );
   double ** VeffTilde = NULL;
   if ( Exc_activated ){ VeffTilde = prepare_excitations( denS ); }
   double Energy = Solver.SolveDAVIDSON( denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nStates - 1, VeffTilde );
   Energy += Prob->gEconst();
   if ( Exc_activated ){ cleanup_excitations( VeffTilde ); }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_TIME_S_SOLVE ] += ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // Decompose the S-object. MPI_CHEMPS2_MASTER decomposes denS. Each MPI process returns the correct discWeight. Each MPI process has the new MPS tensors set.
   gettimeofday( &start, NULL );
   if (( noise_level > 0.0 ) && ( am_i_master )){ denS->addNoise( noise_level ); }
   const double discWeight = denS->Split( MPS[ index ], MPS[ index + 1 ], virtual_dimension, moving_right, change );
   delete denS;
   if ( discWeight > MaxDiscWeightLastSweep ){ MaxDiscWeightLastSweep = discWeight; }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_TIME_S_SPLIT ] += ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   return Energy;

}

int CheMPS2::DMRG::get_num_mps_var() const{

   int num_var = 0;
   for ( int site = 0; site < L; site++ ){
      num_var += MPS[ site ]->gKappa2index( MPS[ site ]->gNKappa() );
   }
   return num_var;

}

void CheMPS2::DMRG::activateExcitations( const int maxExcIn ){

   Exc_activated = true;
   maxExc = maxExcIn;
   Exc_Eshifts = new double[ maxExc ];
   Exc_MPSs = new TensorT ** [ maxExc ];
   Exc_BKs = new SyBookkeeper * [ maxExc ];
   Exc_Overlaps = new TensorO ** [ maxExc ];

}

void CheMPS2::DMRG::newExcitation( const double EshiftIn ){

   assert( Exc_activated );
   assert( nStates - 1 < maxExc );

   if ( the2DM  != NULL ){ delete the2DM;  the2DM  = NULL; }
   if ( the3DM  != NULL ){ delete the3DM;  the3DM  = NULL; }
   if ( theCorr != NULL ){ delete theCorr; theCorr = NULL; }
   deleteAllBoundaryOperators();

   Exc_Eshifts[ nStates - 1 ] = EshiftIn;
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::owner_specific_excitation( L, nStates - 1 ) == MPIchemps2::mpi_rank() ){
   #endif
      Exc_MPSs[ nStates - 1 ] = MPS;
      Exc_BKs[ nStates - 1 ] = denBK;
      Exc_Overlaps[ nStates - 1 ] = new TensorO*[ L - 1 ];
   #ifdef CHEMPS2_MPI_COMPILATION
   } else {
      for ( int site = 0; site < L; site++ ){ delete MPS[ site ]; }
      delete [] MPS;
      delete denBK;
   }
   #endif

   nStates++;

   setupBookkeeperAndMPS();
   PreSolve();

}

