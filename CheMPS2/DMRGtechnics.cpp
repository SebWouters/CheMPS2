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

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>

#include "DMRG.h"
#include "Lapack.h"
#include "Heff.h"
#include "MPIchemps2.h"
#include "Wigner.h"
#include "Special.h"

using std::cout;
using std::endl;

void CheMPS2::DMRG::calc_rdms_and_correlations( const bool do_3rdm, const bool disk_3rdm ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   // Reset timings
   for ( int timecnt = 0; timecnt < CHEMPS2_TIME_VECLENGTH; timecnt++ ){ timings[ timecnt ] = 0.0; }
   num_double_write_disk = 0;
   num_double_read_disk  = 0;
   struct timeval start_global, end_global, start_part, end_part;
   gettimeofday( &start_global, NULL );

   // Get the whole MPS into left-canonical form
   gettimeofday( &start_part, NULL );
   left_normalize( MPS[ L - 2 ], MPS[ L - 1 ] );
   left_normalize( MPS[ L - 1 ], NULL         );
   gettimeofday( &end_part, NULL );
   timings[ CHEMPS2_TIME_S_SPLIT ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

   if ( am_i_master ){
      if ( do_3rdm ){
         cout << "****************************************************" << endl;
         cout << "***  2-RDM, 3-RDM, and Correlations calculation  ***" << endl;
         cout << "****************************************************" << endl;
      } else {
         cout << "********************************************" << endl;
         cout << "***  2-RDM and Correlations calculation  ***" << endl;
         cout << "********************************************" << endl;
      }
   }

   // Make the renormalized operators one site further ( one-dot )
   gettimeofday( &start_part, NULL );
   updateMovingRightSafe( L - 2 );
   gettimeofday( &end_part, NULL );
   timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

   // Calculate the 2DM
   if ( the2DM != NULL ){ delete the2DM; the2DM = NULL; }
   the2DM = new TwoDM( denBK, Prob );

   for ( int siteindex = L - 1; siteindex >= 0; siteindex-- ){

      gettimeofday( &start_part, NULL );
      // Specific 2-RDM entries are internally added per MPI processes; after which an allreduce is called
      the2DM->FillSite( MPS[ siteindex ], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors );
      gettimeofday( &end_part, NULL );
      timings[ CHEMPS2_TIME_S_SOLVE ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

      if ( siteindex > 0 ){

         gettimeofday( &start_part, NULL );
         right_normalize( MPS[ siteindex - 1 ], MPS[ siteindex ] );
         gettimeofday( &end_part, NULL );
         timings[ CHEMPS2_TIME_S_SPLIT ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

         gettimeofday( &start_part, NULL );
         updateMovingLeftSafe2DM( siteindex - 1 );
         gettimeofday( &end_part, NULL );
         timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );
      }
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   gettimeofday( &start_part, NULL );
   the2DM->mpi_allreduce();
   gettimeofday( &end_part, NULL );
   timings[ CHEMPS2_TIME_S_SOLVE ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );
   #endif

   the2DM->correct_higher_multiplicities();

   // Trace, energy, and NOON
   if ( am_i_master ){
      cout << "   N(N-1)                     = " << denBK->gN() * ( denBK->gN() - 1 ) << endl;
      cout << "   Double trace of DMRG 2-RDM = " << the2DM->trace() << endl;
      cout << "   Econst + 0.5 * trace(2DM-A * Ham) = " << the2DM->energy() << endl;
      the2DM->print_noon();
   }

   // Calculate the 3DM and Correlations
   if ( the3DM  != NULL ){ delete the3DM;  the3DM  = NULL; }
   if ( theCorr != NULL ){ delete theCorr; theCorr = NULL; }
   if ( do_3rdm ){ the3DM = new ThreeDM( denBK, Prob, disk_3rdm ); }
   theCorr = new Correlations( denBK, Prob, the2DM );
   if ( am_i_master ){
      Gtensors = new TensorGYZ*[ L - 1 ];
      Ytensors = new TensorGYZ*[ L - 1 ];
      Ztensors = new TensorGYZ*[ L - 1 ];
      Ktensors = new TensorKM *[ L - 1 ];
      Mtensors = new TensorKM *[ L - 1 ];
   }
   if ( do_3rdm ){
      tensor_3rdm_a_J0_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_a_J1_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_a_J1_quartet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_b_J0_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_b_J1_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_b_J1_quartet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_c_J0_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_c_J1_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_c_J1_quartet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_d_J0_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_d_J1_doublet = new Tensor3RDM****[ L - 1 ];
      tensor_3rdm_d_J1_quartet = new Tensor3RDM****[ L - 1 ];

      // Calculate the leftmost site contribution to the 3-RDM
      gettimeofday( &start_part, NULL );
      the3DM->fill_site( MPS[ 0 ], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors,
                         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );
      gettimeofday( &end_part, NULL );
      timings[ CHEMPS2_TIME_S_SOLVE ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );
   }

   for ( int siteindex = 1; siteindex < L; siteindex++ ){

      // Change MPS gauge
      gettimeofday( &start_part, NULL );
      left_normalize( MPS[ siteindex - 1 ], MPS[ siteindex ] );
      gettimeofday( &end_part, NULL );
      timings[ CHEMPS2_TIME_S_SPLIT ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

      // Update 2-RDM, 3-RDM, and Correlations tensors
      gettimeofday( &start_part, NULL );
      if ( do_3rdm ){ update_safe_3rdm_operators( siteindex ); }
      updateMovingRightSafe2DM( siteindex - 1 );
      if ( am_i_master ){ update_correlations_tensors( siteindex ); }
      gettimeofday( &end_part, NULL );
      timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

      // Calculate Correlation and 3-RDM diagrams. Specific contributions per MPI process. Afterwards an MPI allreduce/bcast is required.
      gettimeofday( &start_part, NULL );
      if ( am_i_master ){ theCorr->FillSite( MPS[ siteindex ], Gtensors, Ytensors, Ztensors, Ktensors, Mtensors ); }
      if ( do_3rdm ){ the3DM->fill_site( MPS[ siteindex ], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors,
                                         tensor_3rdm_a_J0_doublet[ siteindex - 1 ], tensor_3rdm_a_J1_doublet[ siteindex - 1 ], tensor_3rdm_a_J1_quartet[ siteindex - 1 ],
                                         tensor_3rdm_b_J0_doublet[ siteindex - 1 ], tensor_3rdm_b_J1_doublet[ siteindex - 1 ], tensor_3rdm_b_J1_quartet[ siteindex - 1 ],
                                         tensor_3rdm_c_J0_doublet[ siteindex - 1 ], tensor_3rdm_c_J1_doublet[ siteindex - 1 ], tensor_3rdm_c_J1_quartet[ siteindex - 1 ],
                                         tensor_3rdm_d_J0_doublet[ siteindex - 1 ], tensor_3rdm_d_J1_doublet[ siteindex - 1 ], tensor_3rdm_d_J1_quartet[ siteindex - 1 ] ); }
      gettimeofday( &end_part, NULL );
      timings[ CHEMPS2_TIME_S_SOLVE ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );
   }

   // Delete the renormalized operators from boundary L-2 and load the ones from boundary L-3
   gettimeofday( &start_part, NULL );
   assert( isAllocated[ L - 2 ] == 1 );                      // Renormalized operators exist on the last boundary (L-2) and are moving to the right.
   assert( isAllocated[ L - 3 ] == 0 );                      // Renormalized operators do not exist on boundary L-3.
     deleteTensors( L - 2, true ); isAllocated[ L - 2 ] = 0; // Delete the renormalized operators on the last boundary (L-2).
   allocateTensors( L - 3, true ); isAllocated[ L - 3 ] = 1; // Create the renormalized operators on boundary L-3.
   OperatorsOnDisk( L - 3, true, false );                    // Load the renormalized operators on boundary L-3.
   gettimeofday( &end_part, NULL );
   timings[ CHEMPS2_TIME_TENS_TOTAL ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );

   #ifdef CHEMPS2_MPI_COMPILATION
   gettimeofday( &start_part, NULL );
   theCorr->mpi_broadcast();
   if (do_3rdm){ the3DM->mpi_allreduce(); }
   gettimeofday( &end_part, NULL );
   timings[ CHEMPS2_TIME_S_SOLVE ] += ( end_part.tv_sec - start_part.tv_sec ) + 1e-6 * ( end_part.tv_usec - start_part.tv_usec );
   #endif

   if ( do_3rdm ){
      the3DM->correct_higher_multiplicities();
      delete_3rdm_operators( L - 1 );
      delete [] tensor_3rdm_a_J0_doublet;
      delete [] tensor_3rdm_a_J1_doublet;
      delete [] tensor_3rdm_a_J1_quartet;
      delete [] tensor_3rdm_b_J0_doublet;
      delete [] tensor_3rdm_b_J1_doublet;
      delete [] tensor_3rdm_b_J1_quartet;
      delete [] tensor_3rdm_c_J0_doublet;
      delete [] tensor_3rdm_c_J1_doublet;
      delete [] tensor_3rdm_c_J1_quartet;
      delete [] tensor_3rdm_d_J0_doublet;
      delete [] tensor_3rdm_d_J1_doublet;
      delete [] tensor_3rdm_d_J1_quartet;
   }

   if ( am_i_master ){
      for ( int previousindex = 0; previousindex < L - 1; previousindex++ ){
         delete Gtensors[ previousindex ];
         delete Ytensors[ previousindex ];
         delete Ztensors[ previousindex ];
         delete Ktensors[ previousindex ];
         delete Mtensors[ previousindex ];
      }
      delete [] Gtensors;
      delete [] Ytensors;
      delete [] Ztensors;
      delete [] Ktensors;
      delete [] Mtensors;
   }

   gettimeofday( &end_global, NULL );
   const double elapsed_global = ( end_global.tv_sec - start_global.tv_sec ) + 1e-6 * ( end_global.tv_usec - start_global.tv_usec );

   if ( am_i_master ){
      cout << "   Single-orbital entropies (Hamiltonian index order is used!) = [ ";
      for ( int index = 0; index < L - 1; index++ ){ cout << theCorr->SingleOrbitalEntropy_HAM( index ) << " , "; }
      cout << theCorr->SingleOrbitalEntropy_HAM( L - 1 ) << " ]." << endl;
      for ( int power = 0; power <= 2; power++ ){
         cout << "   Idistance(" << power << ") = " << theCorr->MutualInformationDistance( (double) power ) << endl;
      }
      if ( do_3rdm ){ cout << "   N(N-1)(N-2)                = " << denBK->gN() * ( denBK->gN() - 1 ) * ( denBK->gN() - 2 ) << endl;
                      cout << "   Triple trace of DMRG 3-RDM = " << the3DM->trace() << endl;
                      cout << "***********************************************************" << endl;
                      cout << "***  Timing information 2-RDM, 3-RDM, and Correlations  ***" << endl;
                      cout << "***********************************************************" << endl; }
               else { cout << "***************************************************" << endl;
                      cout << "***  Timing information 2-RDM and Correlations  ***" << endl;
                      cout << "***************************************************" << endl; }
                      cout << "***     Elapsed wall time        = " << elapsed_global << " seconds" << endl;
                      cout << "***       |--> MPS gauge change  = " << timings[ CHEMPS2_TIME_S_SPLIT ] << " seconds" << endl;
                      cout << "***       |--> Diagram calc      = " << timings[ CHEMPS2_TIME_S_SOLVE ] << " seconds" << endl;
                      print_tensor_update_performance();
      if ( do_3rdm ){ cout << "***********************************************************" << endl; }
               else { cout << "***************************************************" << endl; }
   }

}

void CheMPS2::DMRG::print_tensor_update_performance() const{

    cout << "***       |--> Tensor update     = " << timings[ CHEMPS2_TIME_TENS_TOTAL ] << " seconds" << endl;
    cout << "***              |--> create     = " << timings[ CHEMPS2_TIME_TENS_ALLOC ] << " seconds" << endl;
    cout << "***              |--> destroy    = " << timings[ CHEMPS2_TIME_TENS_FREE  ] << " seconds" << endl;
    cout << "***              |--> disk write = " << timings[ CHEMPS2_TIME_DISK_WRITE ] << " seconds" << endl;
    cout << "***              |--> disk read  = " << timings[ CHEMPS2_TIME_DISK_READ  ] << " seconds" << endl;
    cout << "***              |--> calc       = " << timings[ CHEMPS2_TIME_TENS_CALC  ] << " seconds" << endl;
    cout << "***     Disk write bandwidth     = " << num_double_write_disk * sizeof(double) / ( timings[ CHEMPS2_TIME_DISK_WRITE ] * 1048576 ) << " MB/s" << endl;
    cout << "***     Disk read  bandwidth     = " << num_double_read_disk  * sizeof(double) / ( timings[ CHEMPS2_TIME_DISK_READ  ] * 1048576 ) << " MB/s" << endl;

}

void CheMPS2::DMRG::left_normalize( TensorT * left_mps, TensorT * right_mps ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      const int siteindex = left_mps->gIndex();
      const SyBookkeeper * theBK = left_mps->gBK();
      // (J,N,I) = (0,0,0) and (moving_right, prime_last, jw_phase) = (true, true, false)
      TensorOperator * temp = new TensorOperator( siteindex + 1, 0, 0, 0, true, true, false, theBK, theBK );
      left_mps->QR( temp );
      if ( right_mps != NULL ){ right_mps->LeftMultiply( temp ); }
      delete temp;
   }
   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_tensor( left_mps, MPI_CHEMPS2_MASTER );
   if ( right_mps != NULL ){ MPIchemps2::broadcast_tensor( right_mps, MPI_CHEMPS2_MASTER ); }
   #endif

}

void CheMPS2::DMRG::right_normalize( TensorT * left_mps, TensorT * right_mps ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      const int siteindex = right_mps->gIndex();
      const SyBookkeeper * theBK = right_mps->gBK();
      // (J,N,I) = (0,0,0) and (moving_right, prime_last, jw_phase) = (true, true, false)
      TensorOperator * temp = new TensorOperator( siteindex, 0, 0, 0, true, true, false, theBK, theBK );
      right_mps->LQ( temp );
      if ( left_mps != NULL ){ left_mps->RightMultiply( temp ); }
      delete temp;
   }
   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_tensor( right_mps, MPI_CHEMPS2_MASTER );
   if ( left_mps != NULL ){ MPIchemps2::broadcast_tensor( left_mps, MPI_CHEMPS2_MASTER ); }
   #endif

}

double CheMPS2::DMRG::getSpecificCoefficient(int * coeff) const{
   
   int * alpha = new int[ L ];
   int * beta  = new int[ L ];
   for ( int orb=0; orb<L; orb++ ){
      assert( ( coeff[orb] >= 0 ) && ( coeff[orb] <= 2 ) );
      if ( coeff[orb] == 0 ){ alpha[orb] = 0; beta[orb] = 0; }
      if ( coeff[orb] == 1 ){ alpha[orb] = 1; beta[orb] = 0; }
      if ( coeff[orb] == 2 ){ alpha[orb] = 1; beta[orb] = 1; }
   }
   const double FCIcoeff = getFCIcoefficient( alpha, beta );
   delete [] alpha;
   delete [] beta;
   return FCIcoeff;

}

double CheMPS2::DMRG::getFCIcoefficient(int * alpha, int * beta, const bool mpi_chemps2_master_only) const{

   //DMRGcoeff = alpha/beta[Hamindex = Prob->gf2(DMRGindex)]

   //Check if it's possible
   {
      int nTot = 0;
      int twoSz = 0;
      int iTot = 0;
      for (int DMRGindex=0; DMRGindex<L; DMRGindex++){
         const int HamIndex = (Prob->gReorder()) ? Prob->gf2(DMRGindex) : DMRGindex;
         assert( ( alpha[HamIndex] == 0 ) || ( alpha[HamIndex] == 1 ) );
         assert( (  beta[HamIndex] == 0 ) || (  beta[HamIndex] == 1 ) );
         nTot  += alpha[HamIndex] + beta[HamIndex];
         twoSz += alpha[HamIndex] - beta[HamIndex];
         if ((alpha[HamIndex]+beta[HamIndex])==1){ iTot = Irreps::directProd(iTot,denBK->gIrrep(DMRGindex)); }
      }
      if ( Prob->gN() != nTot ){
         cout << "DMRG::getFCIcoefficient : Ndesired = " << Prob->gN() << " and Ntotal in alpha and beta strings = " << nTot << endl;
         return 0.0;
      }
      // 2Sz can be -Prob->2S() ; -Prob->2S()+2 ; -Prob->2S()+4 ; ... ; Prob->2S()
      if ( ( Prob->gTwoS() < twoSz ) || ( twoSz < -Prob->gTwoS() ) || ( ( Prob->gTwoS() - twoSz ) % 2 != 0 ) ){
         cout << "DMRG::getFCIcoefficient : 2Sdesired = " << Prob->gTwoS() << " and 2Sz in alpha and beta strings = " << twoSz << endl;
         return 0.0;
      }
      if ( Prob->gIrrep() != iTot ){
         cout << "DMRG::getFCIcoefficient : Idesired = " << Prob->gIrrep() << " and Irrep of alpha and beta strings = " << iTot << endl;
         return 0.0;
      }
   }
   
   double theCoeff = 2.0; // A FCI coefficient always lies in between -1.0 and 1.0
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER ) || ( mpi_chemps2_master_only == false ))
   #endif
   {
   
      //Construct necessary arrays
      int Dmax = 1;
      for (int DMRGindex=1; DMRGindex<L; DMRGindex++){
         const int DtotBound = denBK->gTotDimAtBound(DMRGindex);
         if (DtotBound>Dmax){ Dmax = DtotBound; }
      }
      double * arrayL = new double[Dmax];
      double * arrayR = new double[Dmax];
      int * twoSL = new int[L];
      int * twoSR = new int[L];
      int * jumpL = new int[L+1];
      int * jumpR = new int[L+1];
      
      //Start the iterator
      int num_SL = 0;
      jumpL[num_SL] = 0;
      int dimFirst = 1;
      jumpL[num_SL+1] = jumpL[num_SL] + dimFirst;
      twoSL[num_SL] = 0;
      num_SL++;
      arrayL[0] = 1.0;
      int NL = 0;
      int IL = 0;
      int twoSLz = 0;
      
      for (int DMRGindex=0; DMRGindex<L; DMRGindex++){
      
         //Clear the right array
         for (int count = 0; count < Dmax; count++){ arrayR[count] = 0.0; }
         
         //The local occupation
         const int HamIndex = (Prob->gReorder()) ? Prob->gf2(DMRGindex) : DMRGindex;
         const int Nlocal   = alpha[HamIndex] + beta[HamIndex];
         const int twoSzloc = alpha[HamIndex] - beta[HamIndex];
         
         //The right symmetry sectors
         const int NR     = NL + Nlocal;
         const int twoSRz = twoSLz + twoSzloc;
         const int IR     = (( Nlocal == 1 ) ? (Irreps::directProd(IL,denBK->gIrrep(DMRGindex))) : IL);
         
         int num_SR = 0;
         jumpR[num_SR] = 0;
         const int spread = ( ( Nlocal == 1 ) ? 1 : 0 );
         for ( int cntSL = 0; cntSL < num_SL; cntSL++ ){
            for ( int TwoSRattempt = twoSL[cntSL] - spread; TwoSRattempt <= twoSL[cntSL] + spread; TwoSRattempt+=2 ){
               bool encountered = false;
               for ( int cntSR = 0; cntSR < num_SR; cntSR++ ){
                  if ( twoSR[cntSR] == TwoSRattempt ){
                     encountered = true;
                  }
               }
               if ( encountered == false ){
                  const int dimR = denBK->gCurrentDim(DMRGindex+1,NR,TwoSRattempt,IR);
                  if ( dimR > 0 ){
                     jumpR[num_SR+1] = jumpR[num_SR] + dimR;
                     twoSR[num_SR] = TwoSRattempt;
                     num_SR++;
                  }
               }
            }
         }
         assert( jumpR[num_SR] <= Dmax );
         
         for ( int cntSR = 0; cntSR < num_SR; cntSR++ ){
            int TwoSRvalue = twoSR[ cntSR ];
            int dimR = jumpR[ cntSR+1 ] - jumpR[ cntSR ];
            for ( int TwoSLvalue = TwoSRvalue - spread; TwoSLvalue <= TwoSRvalue + spread; TwoSLvalue += 2 ){
            
               int indexSL = -1;
               for ( int cntSL = 0; cntSL < num_SL; cntSL++ ){
                  if ( twoSL[cntSL] == TwoSLvalue ){
                     indexSL = cntSL;
                     cntSL = num_SL; //exit loop
                  }
               }
               if ( indexSL != -1 ){
                  int dimL = jumpL[ indexSL+1 ] - jumpL[ indexSL ];
                  double * Tblock = MPS[DMRGindex]->gStorage(NL,TwoSLvalue,IL,NR,TwoSRvalue,IR);
                  double prefactor = sqrt( TwoSRvalue + 1 )
                                   * Wigner::wigner3j(TwoSLvalue, spread, TwoSRvalue, twoSLz, twoSzloc, -twoSRz)
                                   * Special::phase( -TwoSLvalue + spread - twoSRz );
                  double add2array = 1.0;
                  char notrans = 'N';
                  dgemm_( &notrans, &notrans, &dimFirst, &dimR, &dimL, &prefactor, arrayL + jumpL[indexSL], &dimFirst, Tblock, &dimL, &add2array, arrayR + jumpR[cntSR], &dimFirst);
               }
            }
         }
         
         //Swap L <--> R
         {
            double * temp = arrayR;
            arrayR = arrayL;
            arrayL = temp;
            int * temp2 = twoSR;
            twoSR = twoSL;
            twoSL = temp2;
            temp2 = jumpR;
            jumpR = jumpL;
            jumpL = temp2;
            num_SL = num_SR;
            NL = NR;
            IL = IR;
            twoSLz = twoSRz;
         }
      }
      
      theCoeff = arrayL[0];
      
      assert(   num_SL == 1              );
      assert( jumpL[1] == 1              );
      assert( twoSL[0] == Prob->gTwoS()  );
      assert(       NL == Prob->gN()     );
      assert(       IL == Prob->gIrrep() );
      
      delete [] arrayL;
      delete [] arrayR;
      delete [] twoSL;
      delete [] twoSR;
      delete [] jumpL;
      delete [] jumpR;
   
   }
   
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( mpi_chemps2_master_only ){ MPIchemps2::broadcast_array_double( &theCoeff, 1, MPI_CHEMPS2_MASTER ); }
   #endif
   return theCoeff;

}

double ** CheMPS2::DMRG::prepare_excitations(Sobject * denS){

   double ** VeffTilde = new double*[nStates-1];
   for (int state=0; state<nStates-1; state++){
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_specific_excitation(L,state) == MPIchemps2::mpi_rank() ){
      #endif
         VeffTilde[state] = new double[denS->gKappa2index(denS->gNKappa())];
         calcVeffTilde(VeffTilde[state], denS, state);
      #ifdef CHEMPS2_MPI_COMPILATION
      } else { VeffTilde[state] = NULL; }
      #endif
   }
   return VeffTilde;

}

void CheMPS2::DMRG::cleanup_excitations(double ** VeffTilde) const{

   for (int state=0; state<nStates-1; state++){
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_specific_excitation(L,state) == MPIchemps2::mpi_rank() )
      #endif
      {
         delete [] VeffTilde[state];
      }
   }
   delete [] VeffTilde;

}

void CheMPS2::DMRG::calcVeffTilde(double * result, Sobject * currentS, int state_number){

   int dimTot = currentS->gKappa2index(currentS->gNKappa());
   for (int cnt=0; cnt<dimTot; cnt++){ result[cnt] = 0.0; }
   int index = currentS->gIndex();
   
   const int dimL = std::max(denBK->gMaxDimAtBound(index),   Exc_BKs[state_number]->gMaxDimAtBound(index)   );
   const int dimR = std::max(denBK->gMaxDimAtBound(index+2), Exc_BKs[state_number]->gMaxDimAtBound(index+2) );
   double * workmem = new double[dimL * dimR];
   
   //Construct Sup
   Sobject * Sup = new Sobject( index, Exc_BKs[ state_number ] );
   Sup->Join(Exc_MPSs[state_number][index],Exc_MPSs[state_number][index+1]);
   
   //Construct VeffTilde
   const double prefactor = sqrt(Exc_Eshifts[state_number]) / (Prob->gTwoS() + 1.0);
   for (int ikappa=0; ikappa<currentS->gNKappa(); ikappa++){
      int NL    = currentS->gNL(ikappa);
      int TwoSL = currentS->gTwoSL(ikappa);
      int IL    = currentS->gIL(ikappa);
      int N1    = currentS->gN1(ikappa);
      int N2    = currentS->gN2(ikappa);
      int TwoJ  = currentS->gTwoJ(ikappa);
      int NR    = currentS->gNR(ikappa);
      int TwoSR = currentS->gTwoSR(ikappa);
      int IR    = currentS->gIR(ikappa);
      
      //Check if block also exists for other MPS
      int kappaSup = Sup->gKappa(NL, TwoSL, IL, N1, N2, TwoJ, NR, TwoSR, IR);
      if (kappaSup!=-1){
      
         int dimLdown =                 denBK->gCurrentDim(index,  NL,TwoSL,IL);
         int dimLup   = Exc_BKs[state_number]->gCurrentDim(index,  NL,TwoSL,IL);
         int dimRdown =                 denBK->gCurrentDim(index+2,NR,TwoSR,IR);
         int dimRup   = Exc_BKs[state_number]->gCurrentDim(index+2,NR,TwoSR,IR);
         
         //Do sqrt( (TwoJR+1) * Eshift ) / (TwoStarget+1) times (OL * Sup)_{block} --> workmem
         double * SupPart = Sup->gStorage() + Sup->gKappa2index(kappaSup);
         double alpha = prefactor * sqrt(TwoSR+1.0);
         if (index==0){

            int dimBlock = dimLup * dimRup;
            int inc = 1;
            dcopy_(&dimBlock,SupPart,&inc,workmem,&inc);
            dscal_(&dimBlock,&alpha,workmem,&inc);
            
         } else {
            
            char notrans = 'N';
            double beta = 0.0;
            double * Opart = Exc_Overlaps[state_number][index-1]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
            dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimLup,&alpha,Opart,&dimLdown,SupPart,&dimLup,&beta,workmem,&dimLdown);
            
         }
         
         //Do (workmem * OR)_{block} --> result + jumpCurrentS
         int jumpCurrentS = currentS->gKappa2index(ikappa);
         if (index==L-2){
         
            int dimBlock = dimLdown * dimRdown;
            int inc = 1;
            dcopy_(&dimBlock, workmem, &inc, result + jumpCurrentS, &inc);
         
         } else {
         
            char trans = 'T';
            char notrans = 'N';
            alpha = 1.0;
            double beta = 0.0; //set
            double * Opart = Exc_Overlaps[state_number][index+1]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
            dgemm_(&notrans,&trans,&dimLdown,&dimRdown,&dimRup,&alpha,workmem,&dimLdown,Opart,&dimRdown,&beta,result+jumpCurrentS,&dimLdown);
         
         }
      }
   }
   
   //Deallocate everything
   delete Sup;
   delete [] workmem;

}

void CheMPS2::DMRG::calc_overlaps( const bool moving_right ){

   for ( int state = 0; state < nStates - 1; state++ ){

      double overlap;
      #ifdef CHEMPS2_MPI_COMPILATION
      const int OWNER = MPIchemps2::owner_specific_excitation( L, state );
      if ( OWNER == MPIchemps2::mpi_rank() ){
      #endif
         if ( moving_right ){
            TensorO * first  = new TensorO( L - 1, moving_right, denBK, Exc_BKs[ state ] );
            TensorO * second = new TensorO( L,     moving_right, denBK, Exc_BKs[ state ] );
             first->update_ownmem( MPS[ L - 2 ], Exc_MPSs[ state ][ L - 2 ], Exc_Overlaps[ state ][ L - 3 ] );
            second->update_ownmem( MPS[ L - 1 ], Exc_MPSs[ state ][ L - 1 ], first );
            overlap = second->gStorage()[ 0 ];
            delete second;
            delete first;
         } else {
            TensorO * first  = new TensorO( 1, moving_right, denBK, Exc_BKs[ state ] );
            TensorO * second = new TensorO( 0, moving_right, denBK, Exc_BKs[ state ] );
             first->update_ownmem( MPS[ 1 ], Exc_MPSs[ state ][ 1 ], Exc_Overlaps[ state ][ 1 ] );
            second->update_ownmem( MPS[ 0 ], Exc_MPSs[ state ][ 0 ], first );
            overlap = second->gStorage()[ 0 ] / ( Prob->gTwoS() + 1 );
            delete second;
            delete first;
         }
      #ifdef CHEMPS2_MPI_COMPILATION
      }
      MPIchemps2::broadcast_array_double( &overlap, 1, OWNER );

      if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER )
      #endif
      { cout << "***     The overlap < " << nStates-1 << " | " << state << " > = " << overlap << endl; }

   }

}

