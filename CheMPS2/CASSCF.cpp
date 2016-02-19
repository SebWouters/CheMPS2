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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>
#include <assert.h>

#include "CASSCF.h"
#include "Lapack.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

CheMPS2::CASSCF::CASSCF(Hamiltonian * ham_in, int * docc, int * socc, int * nocc_in, int * ndmrg_in, int * nvirt_in){

   NUCL_ORIG = ham_in->getEconst();
   TMAT_ORIG = ham_in->getTmat();
   VMAT_ORIG = ham_in->getVmat();

   L = ham_in->getL();
   SymmInfo.setGroup( ham_in->getNGroup() );
   num_irreps = SymmInfo.getNumberOfIrreps();

   cout << "DOCC = [ ";
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ cout << docc[ irrep ] << " , "; }
   cout << docc[ num_irreps - 1 ] << " ]" << endl;
   cout << "SOCC = [ ";
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){ cout << socc[ irrep ] << " , "; }
   cout << socc[ num_irreps - 1 ] << " ]" << endl;

   iHandler = new DMRGSCFindices( L, SymmInfo.getGroupNumber(), nocc_in, ndmrg_in, nvirt_in );
   unitary  = new DMRGSCFunitary( iHandler );
   theDIIS = NULL;
   theRotatedTEI = new DMRGSCFintegrals( iHandler );

   //Allocate space for the DMRG 1DM and 2DM
   nOrbDMRG = iHandler->getDMRGcumulative( num_irreps );
   DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
   DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];

   //To calculate the F-matrix and Q-matrix(occ,act) elements only once, and to store them for future access
   theFmatrix = new DMRGSCFmatrix( iHandler );  theFmatrix->clear();
   theQmatOCC = new DMRGSCFmatrix( iHandler );  theQmatOCC->clear();
   theQmatACT = new DMRGSCFmatrix( iHandler );  theQmatACT->clear();
   theQmatWORK= new DMRGSCFmatrix( iHandler ); theQmatWORK->clear();
   theTmatrix = new DMRGSCFmatrix( iHandler );  theTmatrix->clear();

   //To calculate the w_tilde elements only once, and store them for future access
   wmattilde = new DMRGSCFwtilde( iHandler );

   //Print the MO info. This requires the indexHandler to be created...
   checkHF( docc, socc );

   //Print what we have just set up.
   iHandler->Print();

   cout << "DMRGSCF::setupStart : Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;

}

CheMPS2::CASSCF::~CASSCF(){
   
   delete theRotatedTEI;

   delete [] DMRG1DM;
   delete [] DMRG2DM;
   
   //The following objects depend on iHandler: delete them first
   delete theFmatrix;
   delete theQmatOCC;
   delete theQmatACT;
   delete theQmatWORK;
   delete theTmatrix;
   delete wmattilde;
   delete unitary;
   
   delete iHandler;
   if (theDIIS!=NULL){ delete theDIIS; }

}

int CheMPS2::CASSCF::get_num_irreps(){ return num_irreps; }

void CheMPS2::CASSCF::copy2DMover(TwoDM * theDMRG2DM, const int totOrbDMRG, double * localDMRG2DM){

   for (int i1=0; i1<totOrbDMRG; i1++){
      for (int i2=0; i2<totOrbDMRG; i2++){
         for (int i3=0; i3<totOrbDMRG; i3++){
            for (int i4=0; i4<totOrbDMRG; i4++){
               // The assignment has been changed to an addition for state-averaged calculations!
               localDMRG2DM[i1 + totOrbDMRG * ( i2 + totOrbDMRG * (i3 + totOrbDMRG * i4 ) ) ] += theDMRG2DM->getTwoDMA_HAM(i1, i2, i3, i4);
            }
         }
      }
   }

}

void CheMPS2::CASSCF::copy3DMover(ThreeDM * theDMRG3DM, const int numL, double * three_dm){

   for ( int i = 0; i < numL; i++ ){
      for ( int j = 0; j < numL; j++ ){
         for ( int k = 0; k < numL; k++ ){
            for ( int l = 0; l < numL; l++ ){
               for ( int m = 0; m < numL; m++ ){
                  for ( int n = 0; n < numL; n++ ){
                     three_dm[ i + numL * ( j + numL * ( k + numL * ( l + numL * ( m + numL * n ) ) ) ) ] = theDMRG3DM->get_ham_index(i, j, k, l, m, n);
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::CASSCF::setDMRG1DM(const int num_elec, const int numL, double * localDMRG1DM, double * localDMRG2DM){

   const double prefactor = 1.0/( num_elec - 1.0 );

   for ( int cnt1 = 0; cnt1 < numL; cnt1++ ){
      for ( int cnt2 = cnt1; cnt2 < numL; cnt2++ ){
         double value = 0.0;
         for ( int sum = 0; sum < numL; sum++ ){ value += localDMRG2DM[ cnt1 + numL * ( sum + numL * ( cnt2 + numL * sum ) ) ]; }
         localDMRG1DM[ cnt1 + numL * cnt2 ] = prefactor * value;
         localDMRG1DM[ cnt2 + numL * cnt1 ] = localDMRG1DM[ cnt1 + numL * cnt2 ];
      }
   }

}

void CheMPS2::CASSCF::fillLocalizedOrbitalRotations(CheMPS2::DMRGSCFunitary * unitary, CheMPS2::DMRGSCFindices * localIdx, double * eigenvecs){

   const int numIrreps = localIdx->getNirreps();
   const int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   const int size = totOrbDMRG * totOrbDMRG;
   for (int cnt=0; cnt<size; cnt++){ eigenvecs[cnt] = 0.0; }
   int passed = 0;
   for (int irrep=0; irrep<numIrreps; irrep++){

      const int NDMRG = localIdx->getNDMRG(irrep);
      if (NDMRG>0){

         double * blockUnit = unitary->getBlock(irrep);
         double * blockEigs = eigenvecs + passed * ( 1 + totOrbDMRG );

         for (int row=0; row<NDMRG; row++){
            for (int col=0; col<NDMRG; col++){
               blockEigs[row + totOrbDMRG * col] = blockUnit[col + NDMRG * row]; //Eigs = Unit^T
            }
         }

      }

      passed += NDMRG;

   }

}

void CheMPS2::CASSCF::rotateOldToNew(DMRGSCFmatrix * myMatrix){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
   
      int linsize = iHandler->getNORB(irrep);
      double * Umat = unitary->getBlock(irrep);
      double * work = theQmatWORK->getBlock(irrep);
      double * block = myMatrix->getBlock(irrep);
      double alpha = 1.0;
      double beta  = 0.0;
      char trans   = 'T';
      char notrans = 'N';
      dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, Umat, &linsize, block, &linsize, &beta, work,  &linsize);
      dgemm_(&notrans, &trans,   &linsize, &linsize, &linsize, &alpha, work, &linsize, Umat,  &linsize, &beta, block, &linsize);
      
   }

}

void CheMPS2::CASSCF::constructCoulombAndExchangeMatrixInOrigIndices(DMRGSCFmatrix * densityMatrix, DMRGSCFmatrix * resultMatrix){

  for ( int irrepQ = 0; irrepQ < num_irreps; irrepQ++ ){

      const int linearsizeQ = iHandler->getNORB(irrepQ);
      const int numberOfUniqueIndices = (linearsizeQ * (linearsizeQ + 1))/2;

      #pragma omp parallel for schedule(static)
      for (int combinedindex = 0; combinedindex < numberOfUniqueIndices; combinedindex++){

         int colQ = 1;
         while ( (colQ*(colQ+1))/2 <= combinedindex ){ colQ++; }
         colQ -= 1;
         int rowQ = combinedindex - (colQ*(colQ+1))/2;

         double theValue = 0.0;

         for ( int irrepN = 0; irrepN < num_irreps; irrepN++ ){
            const int linearsizeN = iHandler->getNORB( irrepN );
            for (int rowN = 0; rowN < linearsizeN; rowN++){

               theValue += densityMatrix->get(irrepN, rowN, rowN) * ( VMAT_ORIG->get( irrepQ, irrepN, irrepQ, irrepN, rowQ, rowN, colQ, rowN )
                                                              - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, rowN, rowN ) );

               for (int colN = rowN+1; colN < linearsizeN; colN++){

                  theValue += densityMatrix->get(irrepN, rowN, colN) * ( 2 * VMAT_ORIG->get( irrepQ, irrepN, irrepQ, irrepN, rowQ, rowN, colQ, colN )
                                                                     - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, rowN, colN )
                                                                     - 0.5 * VMAT_ORIG->get( irrepQ, irrepQ, irrepN, irrepN, rowQ, colQ, colN, rowN ) );

               }
            }
         }

         resultMatrix->set( irrepQ, rowQ, colQ, theValue );
         resultMatrix->set( irrepQ, colQ, rowQ, theValue );

      }
   }

}

void CheMPS2::CASSCF::buildQmatOCC(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
   
      int linsize = iHandler->getNORB(irrep);
      int NOCC    = iHandler->getNOCC(irrep);
      double alpha = 2.0;
      double beta  = 0.0;
      char trans   = 'T';
      char notrans = 'N';
      double * Umat = unitary->getBlock(irrep);
      double * work = theQmatWORK->getBlock(irrep);
      dgemm_(&trans, &notrans, &linsize, &linsize, &NOCC, &alpha, Umat, &linsize, Umat, &linsize, &beta, work, &linsize);
      
   }
   
   constructCoulombAndExchangeMatrixInOrigIndices( theQmatWORK, theQmatOCC );
   rotateOldToNew( theQmatOCC );

}

void CheMPS2::CASSCF::buildQmatACT(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
   
      int linsize = iHandler->getNORB(irrep);
      int NDMRG   = iHandler->getNDMRG(irrep);
      double alpha = 1.0;
      double beta  = 0.0;
      char trans   = 'T';
      char notrans = 'N';
      double * Umat  =     unitary->getBlock(irrep) + iHandler->getNOCC(irrep);
      double * work  = theQmatWORK->getBlock(irrep);
      double * work2 =  theQmatACT->getBlock(irrep);
      double * RDM = DMRG1DM + iHandler->getDMRGcumulative(irrep) * ( 1 + nOrbDMRG );
      dgemm_(&trans,   &notrans, &linsize, &NDMRG,   &NDMRG, &alpha, Umat,  &linsize, RDM,  &nOrbDMRG, &beta, work2, &linsize);
      dgemm_(&notrans, &notrans, &linsize, &linsize, &NDMRG, &alpha, work2, &linsize, Umat, &linsize,  &beta, work,  &linsize);
      
   }
   
   constructCoulombAndExchangeMatrixInOrigIndices( theQmatWORK, theQmatACT );
   rotateOldToNew( theQmatACT );

}

void CheMPS2::CASSCF::buildTmatrix(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NumORB = iHandler->getNORB(irrep);
      for (int row = 0; row < NumORB; row++){
         for (int col = 0; col < NumORB; col++){
            theTmatrix->set( irrep, row, col, TMAT_ORIG->get( irrep, row, col ) );
         }
      }
   }

   rotateOldToNew( theTmatrix );

}

void CheMPS2::CASSCF::fillConstAndTmatDMRG(Hamiltonian * HamDMRG) const{

   //Constant part of the energy
   double value = NUCL_ORIG;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      for (int orb = 0; orb < iHandler->getNOCC(irrep); orb++){
         value += 2 * theTmatrix->get(irrep, orb, orb) + theQmatOCC->get(irrep, orb, orb);
      }
   }
   HamDMRG->setEconst(value);
   
   //One-body terms: diagonal in the irreps
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int passedDMRG  = iHandler->getDMRGcumulative(irrep);
      const int linsizeDMRG = iHandler->getNDMRG(irrep);
      const int NumOCC      = iHandler->getNOCC(irrep);
      for (int cnt1=0; cnt1<linsizeDMRG; cnt1++){
         for (int cnt2=cnt1; cnt2<linsizeDMRG; cnt2++){
            HamDMRG->setTmat( passedDMRG+cnt1, passedDMRG+cnt2, theTmatrix->get(irrep, NumOCC+cnt1, NumOCC+cnt2) + theQmatOCC->get(irrep, NumOCC+cnt1, NumOCC+cnt2) );
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

void CheMPS2::CASSCF::block_diagonalize( const char space, const DMRGSCFmatrix * Mat, DMRGSCFunitary * Umat, double * work1, double * work2, const DMRGSCFindices * idx, const bool invert, double * localDMRG2RDM ){

   const int n_irreps = idx->getNirreps();
   const int tot_dmrg = idx->getDMRGcumulative( n_irreps );

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
            for ( int col = 0; col < NROTATE/2; col++ ){
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

         // Adjust the 2-RDM accordingly
         if (( localDMRG2RDM != NULL ) && ( space == 'A' )){
            int power[] = { 1,
                            tot_dmrg,
                            tot_dmrg * tot_dmrg,
                            tot_dmrg * tot_dmrg * tot_dmrg,
                            tot_dmrg * tot_dmrg * tot_dmrg * tot_dmrg };
            for ( int index = 3; index >= 0; index-- ){
               for ( int cnt = 0; cnt < power[ 3 - index ]; cnt++ ){
                  double * two_dm_ptr = localDMRG2RDM + power[ index ] * NJUMP + power[ index + 1 ] * cnt;
                  dgemm_( &notrans, &notrans, power + index, &NROTATE, &NROTATE, &one, two_dm_ptr, power + index, work1, &NROTATE, &set, work2, power + index );
                  int size = power[ index ] * NROTATE;
                  int inc1 = 1;
                  dcopy_( &size, work2, &inc1, two_dm_ptr, &inc1 );
               }
            }
         }
      }
   }

}


