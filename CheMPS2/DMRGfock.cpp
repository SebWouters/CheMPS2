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
#include "Excitation.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

void CheMPS2::DMRG::Symm4RDM( double * output, const int Y, const int Z, const bool last_case ){

   struct timeval start, end;
   gettimeofday( &start, NULL );

   assert( the3DM != NULL );

   symm_4rdm_helper( output, Y, Z, 1.0, 1.0, false, 0.5 ); // output = 0.5 *   3rdm[ ( 1 + E_{YZ} + E_{ZY} ) | 0 > ]
   symm_4rdm_helper( output, Y, Z, 1.0, 0.0, true, -0.5 ); // output = 0.5 * ( 3rdm[ ( 1 + E_{YZ} + E_{ZY} ) | 0 > ] - 3rdm[ E_{YZ} + E_{ZY} | 0 > ] )

   for ( int r = 0; r < L; r++ ){
      for ( int q = 0; q < L; q++ ){
         for ( int p = 0; p < L; p++ ){
            for ( int k = 0; k < L; k++ ){
               for ( int j = 0; j < L; j++ ){
                  for ( int i = 0; i < L; i++ ){
                     output[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( i, j, k, p, q, r );
                  }
                  output[ Y + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( Z, j, k, p, q, r );
                  output[ Z + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( Y, j, k, p, q, r );
                  output[ j + L * ( Y + L * ( k + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, Z, k, p, q, r );
                  output[ j + L * ( Z + L * ( k + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, Y, k, p, q, r );
                  output[ j + L * ( k + L * ( Y + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, k, Z, p, q, r );
                  output[ j + L * ( k + L * ( Z + L * ( p + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, k, Y, p, q, r );
                  output[ j + L * ( k + L * ( p + L * ( Y + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, k, p, Z, q, r );
                  output[ j + L * ( k + L * ( p + L * ( Z + L * ( q + L * r )))) ] -= 0.5 * the3DM->get_ham_index( j, k, p, Y, q, r );
                  output[ j + L * ( k + L * ( p + L * ( q + L * ( Y + L * r )))) ] -= 0.5 * the3DM->get_ham_index( k, j, p, Z, q, r );
                  output[ j + L * ( k + L * ( p + L * ( q + L * ( Z + L * r )))) ] -= 0.5 * the3DM->get_ham_index( k, j, p, Y, q, r );
                  output[ j + L * ( k + L * ( p + L * ( q + L * ( r + L * Y )))) ] -= 0.5 * the3DM->get_ham_index( p, j, k, Z, q, r );
                  output[ j + L * ( k + L * ( p + L * ( q + L * ( r + L * Z )))) ] -= 0.5 * the3DM->get_ham_index( p, j, k, Y, q, r );
               }
            }
         }
      }
   }

   if ( last_case ){ PreSolve(); } // Need to set up the renormalized operators again to continue sweeping

   gettimeofday( &end, NULL );
   const double elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER )
   #endif
   { cout << "CheMPS2::DMRG::Symm4RDM( " << Y << " , " << Z << " ) : Elapsed wall time = " << elapsed << " seconds." << endl; }

}

void CheMPS2::DMRG::symm_4rdm_helper( double * output, const int ham_orb1, const int ham_orb2, const double alpha, const double beta, const bool add, const double factor ){

   // Figure out the DMRG orbitals, in order
   assert( ham_orb1 >= 0 );
   assert( ham_orb2 >= 0 );
   assert( ham_orb1 <  L );
   assert( ham_orb2 <  L );
   const int  first_dmrg_orb = (( Prob->gReorder() ) ? Prob->gf1( ham_orb1 ) : ham_orb1 );
   const int second_dmrg_orb = (( Prob->gReorder() ) ? Prob->gf1( ham_orb2 ) : ham_orb2 );
   const int dmrg_orb1 = (( first_dmrg_orb <= second_dmrg_orb ) ?  first_dmrg_orb : second_dmrg_orb );
   const int dmrg_orb2 = (( first_dmrg_orb <= second_dmrg_orb ) ? second_dmrg_orb :  first_dmrg_orb );

   // Make a back-up of the entirely left-normalized MPS and the corresponding bookkeeper
   SyBookkeeper * oldBK = denBK;
   if ( dmrg_orb1 != dmrg_orb2 ){ denBK = new SyBookkeeper( *oldBK ); }
   TensorT ** backup_mps = new TensorT * [ L ];
   for ( int orbital = 0; orbital < L; orbital++ ){
      backup_mps[ orbital ] = MPS[ orbital ];
      MPS[ orbital ] = new TensorT( orbital, denBK ); // denBK is now a DIFFERENT pointer than backup_mps[ orbital ]->gBK()
      int totalsize = MPS[ orbital ]->gKappa2index( MPS[ orbital ]->gNKappa() );
      int inc1 = 1;
      dcopy_( &totalsize, backup_mps[ orbital ]->gStorage(), &inc1, MPS[ orbital ]->gStorage(), &inc1 );
   }
   deleteAllBoundaryOperators();

   // Change the gauge so that the non-orthonormal MPS tensor is on site dmrg_orb2
   for ( int siteindex = L - 1; siteindex > dmrg_orb2; siteindex-- ){
      right_normalize( MPS[ siteindex - 1 ], MPS[ siteindex ] );
      updateMovingLeftSafeFirstTime( siteindex - 1 );
   }

   // Solve
   solve_fock( dmrg_orb1, dmrg_orb2, alpha, beta );

   // Further right normalize the wavefunction except for the first MPS tensor ( contains the norm )
   for ( int siteindex = dmrg_orb2; siteindex > 0; siteindex-- ){
      right_normalize( MPS[ siteindex - 1 ], MPS[ siteindex ] );
      updateMovingLeftSafeFirstTime( siteindex - 1 );
   }

   ThreeDM * helper3rdm = new ThreeDM( denBK, Prob, false );
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

   // Leftmost contribution to the helper3rdm
   helper3rdm->fill_site( MPS[ 0 ], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );

   // Other contributions to the helper3rdm
   for ( int siteindex = 1; siteindex < L; siteindex++ ){

      /* Change the MPS gauge */
      left_normalize( MPS[ siteindex - 1 ], MPS[ siteindex ] );

      /* Update the required renormalized operators */
      update_safe_3rdm_operators( siteindex );
      updateMovingRightSafe2DM( siteindex - 1 );

      /* Current contribution to helper3rdm */
      helper3rdm->fill_site( MPS[ siteindex ], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors,
                             tensor_3rdm_a_J0_doublet[ siteindex - 1 ], tensor_3rdm_a_J1_doublet[ siteindex - 1 ], tensor_3rdm_a_J1_quartet[ siteindex - 1 ],
                             tensor_3rdm_b_J0_doublet[ siteindex - 1 ], tensor_3rdm_b_J1_doublet[ siteindex - 1 ], tensor_3rdm_b_J1_quartet[ siteindex - 1 ],
                             tensor_3rdm_c_J0_doublet[ siteindex - 1 ], tensor_3rdm_c_J1_doublet[ siteindex - 1 ], tensor_3rdm_c_J1_quartet[ siteindex - 1 ],
                             tensor_3rdm_d_J0_doublet[ siteindex - 1 ], tensor_3rdm_d_J1_doublet[ siteindex - 1 ], tensor_3rdm_d_J1_quartet[ siteindex - 1 ] );

   }

   // Collect all data
   #ifdef CHEMPS2_MPI_COMPILATION
   helper3rdm->mpi_allreduce();
   #endif

   // Copy the contributions
   helper3rdm->correct_higher_multiplicities();
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
   helper3rdm->fill_ham_index( factor, add, output, 0, L );

   // Throw out the changed MPS and place back the original left-normalized MPS
   for ( int orbital = 0; orbital < L; orbital++ ){
      delete MPS[ orbital ];
      MPS[ orbital ] = backup_mps[ orbital ];
   }
   delete [] backup_mps;
   if ( dmrg_orb1 != dmrg_orb2 ){
      delete denBK;
      denBK = oldBK;
   }
   delete helper3rdm;
   deleteAllBoundaryOperators();

}

void CheMPS2::DMRG::solve_fock( const int dmrg_orb1, const int dmrg_orb2, const double alpha, const double beta ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( dmrg_orb1 == dmrg_orb2 ){
      MPS[ dmrg_orb1 ]->number_operator( 2 * alpha, beta ); // alpha * ( E_zz + E_zz ) + beta * 1
      return;
   }

   double sweep_inproduct = 0.0;

   if ( dmrg_orb1 + 1 == dmrg_orb2 ){
      Sobject * newS = new Sobject( dmrg_orb1, denBK );
      if ( am_i_master ){
         Sobject * oldS = new Sobject( dmrg_orb1, denBK );
         oldS->Join( MPS[ dmrg_orb1 ], MPS[ dmrg_orb2 ] );
         sweep_inproduct = Excitation::matvec( denBK, denBK, dmrg_orb1, dmrg_orb2, alpha, alpha, beta, newS, oldS, NULL, NULL, NULL );
         delete oldS;
      }
      // MPI_CHEMPS2_MASTER decomposes newS. Each MPI process returns the correct discarded_weight. Each MPI process has the new MPS tensors set.
      const double discarded_weight = newS->Split( MPS[ dmrg_orb1 ], MPS[ dmrg_orb2 ], OptScheme->get_D( OptScheme->get_number() - 1 ), true, true );
      delete newS;
   }

   if ( dmrg_orb1 + 1 < dmrg_orb2 ){

      SyBookkeeper * newBK = denBK;
      denBK = new SyBookkeeper( *newBK );
      newBK->restart( dmrg_orb1 + 1, dmrg_orb2, OptScheme->get_D( OptScheme->get_number() - 1 ) );
      TensorT ** old_mps = new TensorT * [ L ];
      for ( int orbital = dmrg_orb1; orbital <= dmrg_orb2; orbital++ ){
         old_mps[ orbital ] = MPS[ orbital ];
         old_mps[ orbital ]->sBK( denBK );
         MPS[ orbital ] = new TensorT( orbital, newBK );
         MPS[ orbital ]->random();
         left_normalize( MPS[ orbital ], NULL ); // MPI_CHEMPS2_MASTER broadcasts MPS[ orbital ] ( left-normalized ).
      }

      TensorO ** overlaps = NULL;
      TensorL ** regular  = NULL;
      TensorL ** trans    = NULL;
      if ( am_i_master ){
         overlaps = new TensorO*[ L - 1 ];
         regular  = new TensorL*[ L - 1 ];
         trans    = new TensorL*[ L - 1 ];
         for ( int cnt = 0; cnt < L - 1; cnt++ ){
            overlaps[ cnt ] = NULL;
             regular[ cnt ] = NULL;
               trans[ cnt ] = NULL;
         }
         for ( int orbital = dmrg_orb1; orbital < dmrg_orb2 - 1; orbital++ ){
            solve_fock_update_helper( orbital, dmrg_orb1, dmrg_orb2, true, MPS, old_mps, newBK, denBK, overlaps, regular, trans );
         }
      }

      // Sweeps
      bool change = false;
      for ( int instruction = 0; instruction < OptScheme->get_number(); instruction++ ){
         int num_iterations = 0;
         double previous_inproduct = sweep_inproduct + 10 * OptScheme->get_energy_conv( instruction );
         while (( fabs( sweep_inproduct - previous_inproduct ) > OptScheme->get_energy_conv( instruction ) ) && ( num_iterations < OptScheme->get_max_sweeps( instruction ) )){
            {
               const double noise_level = fabs( OptScheme->get_noise_prefactor( instruction ) ) * MaxDiscWeightLastSweep;
               MaxDiscWeightLastSweep   = 0.0;
               for ( int index = dmrg_orb2 - 1; index > dmrg_orb1; index-- ){
                  Sobject * newS = new Sobject( index, newBK );
                  if ( am_i_master ){
                     Sobject * oldS = new Sobject( index, denBK );
                     oldS->Join( old_mps[ index ], old_mps[ index + 1 ] );
                     sweep_inproduct = Excitation::matvec( newBK, denBK, dmrg_orb1, dmrg_orb2, alpha, alpha, beta, newS, oldS, overlaps, regular, trans );
                     delete oldS;
                     if ( noise_level > 0.0 ){ newS->addNoise( noise_level ); }
                  }
                  // MPI_CHEMPS2_MASTER decomposes newS. Each MPI process returns the correct discarded_weight. Each MPI process has the new MPS tensors set.
                  const double discarded_weight = newS->Split( MPS[ index ], MPS[ index + 1 ], OptScheme->get_D( instruction ), false, change );
                  if ( discarded_weight > MaxDiscWeightLastSweep ){ MaxDiscWeightLastSweep = discarded_weight; }
                  delete newS;
                  if ( am_i_master ){
                     solve_fock_update_helper( index, dmrg_orb1, dmrg_orb2, false, MPS, old_mps, newBK, denBK, overlaps, regular, trans );
                  }
               }
            }
            change = true;
            {
               const double noise_level = fabs( OptScheme->get_noise_prefactor( instruction ) ) * MaxDiscWeightLastSweep;
               MaxDiscWeightLastSweep   = 0.0;
               for ( int index = dmrg_orb1; index < dmrg_orb2 - 1; index++ ){
                  Sobject * newS = new Sobject( index, newBK );
                  if ( am_i_master ){
                     Sobject * oldS = new Sobject( index, denBK );
                     oldS->Join( old_mps[ index ], old_mps[ index + 1 ] );
                     sweep_inproduct = Excitation::matvec( newBK, denBK, dmrg_orb1, dmrg_orb2, alpha, alpha, beta, newS, oldS, overlaps, regular, trans );
                     delete oldS;
                     if ( noise_level > 0.0 ){ newS->addNoise( noise_level ); }
                  }
                  // MPI_CHEMPS2_MASTER decomposes newS. Each MPI process returns the correct discarded_weight. Each MPI process has the new MPS tensors set.
                  const double discarded_weight = newS->Split( MPS[ index ], MPS[ index + 1 ], OptScheme->get_D( instruction ), true, change );
                  if ( discarded_weight > MaxDiscWeightLastSweep ){ MaxDiscWeightLastSweep = discarded_weight; }
                  delete newS;
                  if ( am_i_master ){
                     solve_fock_update_helper( index, dmrg_orb1, dmrg_orb2, true, MPS, old_mps, newBK, denBK, overlaps, regular, trans );
                  }
               }
            }
            #ifdef CHEMPS2_MPI_COMPILATION
            CheMPS2::MPIchemps2::broadcast_array_double( &sweep_inproduct, 1, MPI_CHEMPS2_MASTER );
            #endif
            previous_inproduct = sweep_inproduct;
            num_iterations++;
         }
      }

      if ( am_i_master ){
         for ( int index = 0; index < L - 1; index++ ){
            if ( overlaps[ index ] != NULL ){ delete overlaps[ index ]; }
            if (  regular[ index ] != NULL ){ delete  regular[ index ]; }
            if (    trans[ index ] != NULL ){ delete    trans[ index ]; }
         }
         delete [] overlaps;
         delete [] regular;
         delete [] trans;
      }

      for ( int orbital = dmrg_orb1; orbital <= dmrg_orb2; orbital++ ){ delete old_mps[ orbital ]; }
      delete [] old_mps;
      delete denBK;
      denBK = newBK;

   }

   if ( am_i_master ){
      const double rdm_inproduct = 2 * alpha * the2DM->get1RDM_DMRG( dmrg_orb1, dmrg_orb2 ) + beta;
      cout << "DMRG::solve_fock : Accuracy = " << fabs( sweep_inproduct / ( Prob->gTwoS() + 1 ) - rdm_inproduct ) << endl;
   }

}

void CheMPS2::DMRG::solve_fock_update_helper( const int index, const int dmrg_orb1, const int dmrg_orb2, const bool moving_right, TensorT ** new_mps, TensorT ** old_mps, SyBookkeeper * new_bk, SyBookkeeper * old_bk, TensorO ** overlaps, TensorL ** regular, TensorL ** trans ){

   if ( overlaps[ index ] != NULL ){ delete overlaps[ index ]; }
   if (  regular[ index ] != NULL ){ delete  regular[ index ]; }
   if (    trans[ index ] != NULL ){ delete    trans[ index ]; }

   const int Idiff = new_bk->gIrrep( dmrg_orb1 );
   assert( Idiff == new_bk->gIrrep( dmrg_orb2 ) );
   overlaps[ index ] = new TensorO( index + 1,        moving_right, new_bk, old_bk );
    regular[ index ] = new TensorL( index + 1, Idiff, moving_right, new_bk, old_bk );
      trans[ index ] = new TensorL( index + 1, Idiff, moving_right, old_bk, new_bk );

   if ( moving_right ){
      if ( index == dmrg_orb1 ){
         overlaps[ index ]->create( new_mps[ index ], old_mps[ index ] );
          regular[ index ]->create( new_mps[ index ], old_mps[ index ], NULL, NULL );
            trans[ index ]->create( old_mps[ index ], new_mps[ index ], NULL, NULL );
      } else {
         const int dimL = std::max( new_bk->gMaxDimAtBound( index     ), old_bk->gMaxDimAtBound( index     ) );
         const int dimR = std::max( new_bk->gMaxDimAtBound( index + 1 ), old_bk->gMaxDimAtBound( index + 1 ) );
         double * workmem = new double[ dimL * dimR ];
         overlaps[ index ]->update( overlaps[ index - 1 ], new_mps[ index ], old_mps[ index ], workmem );
          regular[ index ]->update(  regular[ index - 1 ], new_mps[ index ], old_mps[ index ], workmem );
            trans[ index ]->update(    trans[ index - 1 ], old_mps[ index ], new_mps[ index ], workmem );
         delete [] workmem;
      }
   } else {
      if ( index + 1 == dmrg_orb2 ){
         overlaps[ index ]->create( new_mps[ index + 1 ], old_mps[ index + 1 ] );
          regular[ index ]->create( new_mps[ index + 1 ], old_mps[ index + 1 ], NULL, NULL );
            trans[ index ]->create( old_mps[ index + 1 ], new_mps[ index + 1 ], NULL, NULL );
      } else {
         const int dimL = std::max( new_bk->gMaxDimAtBound( index + 1 ), old_bk->gMaxDimAtBound( index + 1 ) );
         const int dimR = std::max( new_bk->gMaxDimAtBound( index + 2 ), old_bk->gMaxDimAtBound( index + 2 ) );
         double * workmem = new double[ dimL * dimR ];
         overlaps[ index ]->update( overlaps[ index + 1 ], new_mps[ index + 1 ], old_mps[ index + 1 ], workmem );
          regular[ index ]->update(  regular[ index + 1 ], new_mps[ index + 1 ], old_mps[ index + 1 ], workmem );
            trans[ index ]->update(    trans[ index + 1 ], old_mps[ index + 1 ], new_mps[ index + 1 ], workmem );
         delete [] workmem;
      }
   }

}


