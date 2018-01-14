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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "CASSCF.h"
#include "Lapack.h"
#include "Special.h"
#include "MPIchemps2.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

CheMPS2::CASSCF::CASSCF( Hamiltonian * ham_in, int * docc, int * socc, int * nocc, int * ndmrg, int * nvirt, const string new_tmp_folder ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   NUCL_ORIG = ham_in->getEconst();
   TMAT_ORIG = ham_in->getTmat();
   VMAT_ORIG = ham_in->getVmat();

   L = ham_in->getL();
   SymmInfo.setGroup( ham_in->getNGroup() );
   num_irreps = SymmInfo.getNumberOfIrreps();
   successful_solve = false;

   if (( am_i_master ) && ( docc != NULL ) && ( socc != NULL )){
      cout << "DOCC = [ ";
      for ( int irrep = 0; irrep < num_irreps-1; irrep++ ){ cout << docc[ irrep ] << " , "; }
      cout << docc[ num_irreps - 1 ] << " ]" << endl;
      cout << "SOCC = [ ";
      for ( int irrep = 0; irrep < num_irreps-1; irrep++ ){ cout << socc[ irrep ] << " , "; }
      cout << socc[ num_irreps - 1 ] << " ]" << endl;
   }

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int norb_in  = nocc[ irrep ] + ndmrg[ irrep ] + nvirt[ irrep ];
      const int norb_ham = VMAT_ORIG->get_irrep_size( irrep );
      if (( norb_ham != norb_in ) && ( am_i_master )){
         cout << "CASSCF::CASSCF : nocc[" << irrep << "] + ndmrg[" << irrep << "] + nvirt[" << irrep << "] = " << norb_in
              << " and in the Hamiltonian norb[" << irrep << "] = " << norb_ham << "." << endl;
      }
      assert( norb_ham == norb_in );
   }

   iHandler = new DMRGSCFindices( L, SymmInfo.getGroupNumber(), nocc, ndmrg, nvirt );
   unitary  = new DMRGSCFunitary( iHandler );

   // Allocate space for the DMRG 1DM and 2DM
   nOrbDMRG = iHandler->getDMRGcumulative( num_irreps );
   DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
   DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];

   // To store the F-matrix and Q-matrix(occ,act)
   theFmatrix = new DMRGSCFmatrix( iHandler );  theFmatrix->clear();
   theQmatOCC = new DMRGSCFmatrix( iHandler );  theQmatOCC->clear();
   theQmatACT = new DMRGSCFmatrix( iHandler );  theQmatACT->clear();
   theQmatWORK= new DMRGSCFmatrix( iHandler ); theQmatWORK->clear();
   theTmatrix = new DMRGSCFmatrix( iHandler );  theTmatrix->clear();

   if ( am_i_master ){
      if (( docc != NULL ) && ( socc != NULL )){ checkHF( docc, socc ); } // Print the MO info. This requires the iHandler to be created...
      iHandler->Print();
      cout << "DMRGSCF::setupStart : Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;
   }

   this->tmp_folder = new_tmp_folder;

}

CheMPS2::CASSCF::~CASSCF(){

   delete [] DMRG1DM;
   delete [] DMRG2DM;

   //The following objects depend on iHandler: delete them first
   delete theFmatrix;
   delete theQmatOCC;
   delete theQmatACT;
   delete theQmatWORK;
   delete theTmatrix;
   delete unitary;

   delete iHandler;

}

int CheMPS2::CASSCF::get_num_irreps(){ return num_irreps; }

void CheMPS2::CASSCF::copy2DMover( TwoDM * theDMRG2DM, const int LAS, double * two_dm ){

   for ( int i1 = 0; i1 < LAS; i1++ ){
      for ( int i2 = 0; i2 < LAS; i2++ ){
         for ( int i3 = 0; i3 < LAS; i3++ ){
            for ( int i4 = 0; i4 < LAS; i4++ ){
               // The assignment has been changed to an addition for state-averaged calculations!
               two_dm[ i1 + LAS * ( i2 + LAS * ( i3 + LAS * i4 ) ) ] += theDMRG2DM->getTwoDMA_HAM( i1, i2, i3, i4 );
            }
         }
      }
   }

}

void CheMPS2::CASSCF::setDMRG1DM( const int num_elec, const int LAS, double * one_dm, double * two_dm ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      const double prefactor = 1.0 / ( num_elec - 1 );
      for ( int cnt1 = 0; cnt1 < LAS; cnt1++ ){
         for ( int cnt2 = cnt1; cnt2 < LAS; cnt2++ ){
            double value = 0.0;
            for ( int sum = 0; sum < LAS; sum++ ){ value += two_dm[ cnt1 + LAS * ( sum + LAS * ( cnt2 + LAS * sum ) ) ]; }
            one_dm[ cnt1 + LAS * cnt2 ] = prefactor * value;
            one_dm[ cnt2 + LAS * cnt1 ] = one_dm[ cnt1 + LAS * cnt2 ];
         }
      }
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   MPIchemps2::broadcast_array_double( one_dm, LAS * LAS, MPI_CHEMPS2_MASTER );
   #endif

}

void CheMPS2::CASSCF::fillLocalizedOrbitalRotations( DMRGSCFunitary * umat, DMRGSCFindices * idx, double * eigenvecs ){

   const int n_irreps = idx->getNirreps();
   const int tot_dmrg = idx->getDMRGcumulative( n_irreps );
   const int size = tot_dmrg * tot_dmrg;

   for ( int cnt = 0; cnt < size; cnt++ ){ eigenvecs[ cnt ] = 0.0; }
   int passed = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NDMRG = idx->getNDMRG( irrep );
      if ( NDMRG > 0 ){

         double * Ublock = umat->getBlock( irrep );
         double * Eblock = eigenvecs + passed * ( 1 + tot_dmrg );

         for ( int row = 0; row < NDMRG; row++ ){
            for ( int col = 0; col < NDMRG; col++ ){
               Eblock[ row + tot_dmrg * col ] = Ublock[ col + NDMRG * row ]; //Eigs = Unit^T
            }
         }

      }

      passed += NDMRG;

   }

}

void CheMPS2::CASSCF::rotateOldToNew( DMRGSCFmatrix * myMatrix ){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      int NORB = iHandler->getNORB( irrep );
      if ( NORB > 0 ){
         double * Umat = unitary->getBlock( irrep );
         double * work = theQmatWORK->getBlock( irrep );
         double * block = myMatrix->getBlock( irrep );
         double one = 1.0;
         double set = 0.0;
         char trans   = 'T';
         char notrans = 'N';
         dgemm_( &notrans, &notrans, &NORB, &NORB, &NORB, &one, Umat, &NORB, block, &NORB, &set, work,  &NORB );
         dgemm_( &notrans, &trans,   &NORB, &NORB, &NORB, &one, work, &NORB, Umat,  &NORB, &set, block, &NORB );
      }
   }

}

void CheMPS2::CASSCF::constructCoulombAndExchangeMatrixInOrigIndices( DMRGSCFmatrix * density, DMRGSCFmatrix * result ){

  for ( int irrepQ = 0; irrepQ < num_irreps; irrepQ++ ){

      const int linearsizeQ = iHandler->getNORB( irrepQ );
      const int triangsizeQ = ( linearsizeQ * ( linearsizeQ + 1 ) ) / 2;
      int myindices[ 2 ];

      #pragma omp parallel for schedule(static)
      for ( int combinedindex = 0; combinedindex < triangsizeQ; combinedindex++ ){

         Special::invert_triangle_two( combinedindex, myindices );
         const int rowQ = myindices[ 0 ];
         const int colQ = myindices[ 1 ];

         double value = 0.0;

         for ( int irrepN = 0; irrepN < num_irreps; irrepN++ ){
            const int linearsizeN = iHandler->getNORB( irrepN );
            for ( int rowN = 0; rowN < linearsizeN; rowN++ ){

               value += density->get( irrepN, rowN, rowN ) * ( VMAT_ORIG->get( irrepQ, irrepN, irrepQ, irrepN, rowQ, rowN, colQ, rowN )
                                                       - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, rowN, rowN ) );

               for ( int colN = rowN + 1; colN < linearsizeN; colN++ ){

                  value += density->get( irrepN, rowN, colN ) * ( 2 * VMAT_ORIG->get( irrepQ, irrepN, irrepQ, irrepN, rowQ, rowN, colQ, colN )
                                                              - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, rowN, colN )
                                                              - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, colN, rowN ) );

               }
            }
         }

         result->set( irrepQ, rowQ, colQ, value );
         result->set( irrepQ, colQ, rowQ, value );

      }
   }

}

void CheMPS2::CASSCF::buildQmatOCC(){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){

         int NORB = iHandler->getNORB( irrep );
         if ( NORB > 0 ){
            int NOCC = iHandler->getNOCC( irrep );
            double two = 2.0;
            double set = 0.0;
            char trans   = 'T';
            char notrans = 'N';
            double * Umat = unitary->getBlock( irrep );
            double * work = theQmatWORK->getBlock( irrep );
            dgemm_( &trans, &notrans, &NORB, &NORB, &NOCC, &two, Umat, &NORB, Umat, &NORB, &set, work, &NORB );
         }
      }
      constructCoulombAndExchangeMatrixInOrigIndices( theQmatWORK, theQmatOCC );
      rotateOldToNew( theQmatOCC );
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   theQmatOCC->broadcast( MPI_CHEMPS2_MASTER );
   #endif

}

void CheMPS2::CASSCF::buildQmatACT(){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){

         int NORB = iHandler->getNORB( irrep );
         if ( NORB > 0 ){
            int NACT = iHandler->getNDMRG( irrep );
            double one = 1.0;
            double set = 0.0;
            char trans   = 'T';
            char notrans = 'N';
            double * Umat  =     unitary->getBlock( irrep ) + iHandler->getNOCC( irrep );
            double * work  = theQmatWORK->getBlock( irrep );
            double * work2 =  theQmatACT->getBlock( irrep );
            double * RDM = DMRG1DM + iHandler->getDMRGcumulative( irrep ) * ( 1 + nOrbDMRG );
            dgemm_( &trans,   &notrans, &NORB, &NACT, &NACT, &one, Umat,  &NORB, RDM,  &nOrbDMRG, &set, work2, &NORB );
            dgemm_( &notrans, &notrans, &NORB, &NORB, &NACT, &one, work2, &NORB, Umat, &NORB,     &set, work,  &NORB );
         }
      }
      constructCoulombAndExchangeMatrixInOrigIndices( theQmatWORK, theQmatACT );
      rotateOldToNew( theQmatACT );
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   theQmatACT->broadcast( MPI_CHEMPS2_MASTER );
   #endif

}

void CheMPS2::CASSCF::buildTmatrix(){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   if ( am_i_master ){
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int NumORB = iHandler->getNORB( irrep );
         for ( int row = 0; row < NumORB; row++ ){
            for ( int col = 0; col < NumORB; col++ ){
               theTmatrix->set( irrep, row, col, TMAT_ORIG->get( irrep, row, col ) );
            }
         }
      }
      rotateOldToNew( theTmatrix );
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   theTmatrix->broadcast( MPI_CHEMPS2_MASTER );
   #endif

}

void CheMPS2::CASSCF::fillConstAndTmatDMRG( Hamiltonian * HamDMRG ) const{

   //Constant part of the energy
   double constant = NUCL_ORIG;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = iHandler->getNOCC( irrep );
      for ( int occ = 0; occ < NOCC; occ++ ){
         constant += ( 2 * theTmatrix->get( irrep, occ, occ )
                         + theQmatOCC->get( irrep, occ, occ ) );
      }
   }
   HamDMRG->setEconst( constant );

   //One-body terms: diagonal in the irreps
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int JUMP = iHandler->getDMRGcumulative( irrep );
      const int NACT = iHandler->getNDMRG( irrep );
      const int NOCC = iHandler->getNOCC( irrep );
      for ( int cnt1 = 0; cnt1 < NACT; cnt1++ ){
         for ( int cnt2 = cnt1; cnt2 < NACT; cnt2++ ){
            HamDMRG->setTmat( JUMP + cnt1, JUMP + cnt2, ( theTmatrix->get( irrep, NOCC + cnt1, NOCC + cnt2 )
                                                        + theQmatOCC->get( irrep, NOCC + cnt1, NOCC + cnt2 ) ) );
         }
      }
   }

}

void CheMPS2::CASSCF::copy_active( double * origin, DMRGSCFmatrix * result, const DMRGSCFindices * idx, const bool one_rdm ){

   result->clear();

   const int n_irreps = idx->getNirreps();
   const int tot_dmrg = idx->getDMRGcumulative( n_irreps );

   int passed = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NOCC = idx->getNOCC( irrep );
      const int NACT = idx->getNDMRG( irrep );

      if ( one_rdm ){
         for ( int orb = 0; orb < NOCC; orb++ ){
            result->set( irrep, orb, orb, 2.0 );
         }
      }

      for ( int row = 0; row < NACT; row++ ){
         for ( int col = 0; col < NACT; col++ ){
            result->set( irrep, NOCC + row, NOCC + col, origin[ passed + row + tot_dmrg * ( passed + col ) ] );
         }
      }

      passed += NACT;
   }

   assert( passed == tot_dmrg );

}

void CheMPS2::CASSCF::copy_active( const DMRGSCFmatrix * origin, double * result, const DMRGSCFindices * idx ){

   const int n_irreps = idx->getNirreps();
   const int tot_dmrg = idx->getDMRGcumulative( n_irreps );

   for ( int cnt = 0; cnt < tot_dmrg * tot_dmrg; cnt++ ){ result[ cnt ] = 0.0; }

   int passed = 0;
   for ( int irrep = 0; irrep < n_irreps; irrep++ ){

      const int NOCC = idx->getNOCC( irrep );
      const int NACT = idx->getNDMRG( irrep );

      for ( int row = 0; row < NACT; row++ ){
         for ( int col = 0; col < NACT; col++ ){
            result[ passed + row + tot_dmrg * ( passed + col ) ] = origin->get( irrep, NOCC + row, NOCC + col );
         }
      }

      passed += NACT;
   }

   assert( passed == tot_dmrg );

}

void CheMPS2::CASSCF::block_diagonalize( const char space, const DMRGSCFmatrix * Mat, DMRGSCFunitary * Umat, double * work1, double * work2, const DMRGSCFindices * idx, const bool invert, double * two_dm, double * three_dm, double * contract ){

   #ifdef CHEMPS2_MPI_COMPILATION
      const bool am_i_master = ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER );
   #else
      const bool am_i_master = true;
   #endif

   const int n_irreps = idx->getNirreps();
   const int tot_dmrg = idx->getDMRGcumulative( n_irreps );

   if ( am_i_master ){

      for ( int irrep = 0; irrep < n_irreps; irrep++ ){

         int NTOTAL  = idx->getNORB( irrep );
         int NROTATE = (( space == 'O' ) ? idx->getNOCC( irrep ) : (( space == 'A' ) ? idx->getNDMRG( irrep ) : idx->getNVIRT( irrep )));
         int NSHIFT  = (( space == 'O' ) ? 0                     : (( space == 'A' ) ? idx->getNOCC( irrep )  : NTOTAL - NROTATE      ));
         int NJUMP   = idx->getDMRGcumulative( irrep );
         if ( NROTATE > 1 ){

            // Diagonalize the relevant block of Mat
            {
               for ( int row = 0; row < NROTATE; row++ ){
                  for ( int col = 0; col < NROTATE; col++ ){
                     work1[ row + NROTATE * col ] = Mat->get( irrep, NSHIFT + row, NSHIFT + col );
                  }
               }
               char jobz = 'V';
               char uplo = 'U';
               int info;
               int size = max( 3 * NROTATE - 1, NROTATE * NROTATE );
               dsyev_( &jobz, &uplo, &NROTATE, work1, &NROTATE, work2 + size, work2, &size, &info );
            }

            // Invert the order (large to small)
            if ( invert ){
               for ( int col = 0; col < NROTATE / 2; col++ ){
                  for ( int row = 0; row < NROTATE; row++ ){
                     const double temp = work1[ row + NROTATE * ( NROTATE - 1 - col ) ];
                     work1[ row + NROTATE * ( NROTATE - 1 - col ) ] = work1[ row + NROTATE * col ];
                     work1[ row + NROTATE * col ] = temp;
                  }
               }
            }

            // Adjust the u-matrix accordingly
            double * umatrix = Umat->getBlock( irrep ) + NSHIFT;
            for ( int row = 0; row < NROTATE; row++ ){
               for ( int col = 0; col < NTOTAL; col++ ){
                  work2[ row + NROTATE * col ] = umatrix[ row + NTOTAL * col ];
               }
            }
            char trans   = 'T';
            char notrans = 'N';
            double one = 1.0;
            double set = 0.0;
            dgemm_( &trans, &notrans, &NROTATE, &NTOTAL, &NROTATE, &one, work1, &NROTATE, work2, &NROTATE, &set, umatrix, &NTOTAL );

            // Adjust the two_dm, three_dm, and contract objects accordingly
            if ( space == 'A' ){
               if (   two_dm != NULL ){ rotate_active_space_object( 4,   two_dm, work2, work1, tot_dmrg, NJUMP, NROTATE ); }
               if ( three_dm != NULL ){ rotate_active_space_object( 6, three_dm, work2, work1, tot_dmrg, NJUMP, NROTATE ); }
               if ( contract != NULL ){ rotate_active_space_object( 6, contract, work2, work1, tot_dmrg, NJUMP, NROTATE ); }
            }
         }
      }
   }

   #ifdef CHEMPS2_MPI_COMPILATION
   Umat->broadcast( MPI_CHEMPS2_MASTER );
   if ( space == 'A' ){
      if (   two_dm != NULL ){ MPIchemps2::broadcast_array_double(   two_dm, tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg,                       MPI_CHEMPS2_MASTER ); }
      if ( three_dm != NULL ){ MPIchemps2::broadcast_array_double( three_dm, tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg, MPI_CHEMPS2_MASTER ); }
      if ( contract != NULL ){ MPIchemps2::broadcast_array_double( contract, tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg, MPI_CHEMPS2_MASTER ); }
   }
   #endif

}

void CheMPS2::CASSCF::rotate_active_space_object( const int num_indices, double * object, double * work, double * rotation, const int LAS, const int NJUMP, const int NROTATE ){

   assert( num_indices >= 2 );
   assert( num_indices <= 6 );

   int power[] = { 1,
                   LAS,
                   LAS * LAS,
                   LAS * LAS * LAS,
                   LAS * LAS * LAS * LAS,
                   LAS * LAS * LAS * LAS * LAS,
                   LAS * LAS * LAS * LAS * LAS * LAS };

   for ( int rot_index = num_indices - 1; rot_index >= 0; rot_index-- ){
      for ( int block = 0; block < power[ num_indices - 1 - rot_index ]; block++ ){
         double * mat = object + power[ rot_index ] * NJUMP + power[ rot_index + 1 ] * block;
         int ROTDIM = NROTATE;
         char notrans = 'N';
         double one = 1.0;
         double set = 0.0;
         dgemm_( &notrans, &notrans, power + rot_index, &ROTDIM, &ROTDIM, &one, mat, power + rot_index, rotation, &ROTDIM, &set, work, power + rot_index );
         int size = power[ rot_index ] * NROTATE;
         int inc1 = 1;
         dcopy_( &size, work, &inc1, mat, &inc1 );
      }
   }

}


