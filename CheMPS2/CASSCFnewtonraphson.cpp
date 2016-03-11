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
#include <string>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <sstream>
#include <assert.h>

#include "CASSCF.h"
#include "Lapack.h"
#include "DMRGSCFVmatRotations.h"
#include "EdmistonRuedenberg.h"
#include "Davidson.h"
#include "FCI.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

void CheMPS2::CASSCF::delete_file( const string filename ){

   std::stringstream temp;
   temp << "rm " << filename;
   int info = system( temp.str().c_str() );
   cout << "Info on system( " << temp.str() << " ) = " << info << endl;

}

double CheMPS2::CASSCF::solve( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * scf_options ){

   const int num_elec = Nelectrons - 2 * iHandler->getNOCCsum();
   assert( num_elec >= 0 );
   assert(( OptScheme != NULL ) || (( OptScheme == NULL ) && ( rootNum == 1 )));

   // Convergence variables
   double gradNorm = 1.0;
   double updateNorm = 1.0;
   double * gradient = new double[ unitary->getNumVariablesX() ];
   for ( int cnt = 0; cnt < unitary->getNumVariablesX(); cnt++ ){ gradient[ cnt ] = 0.0; }
   double * diis_vec = NULL;
   double Energy = 1e8;

   // The CheMPS2::Problem for the inner DMRG calculation
   Hamiltonian * HamDMRG = new Hamiltonian(nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital());
   Problem * Prob = new Problem(HamDMRG, TwoS, num_elec, Irrep);
   Prob->SetupReorderD2h(); //Doesn't matter if the group isn't D2h, Prob checks it.

   // Determine the maximum NORB(irrep) and the max_block_size for the ERI orbital rotation
   const int  maxlinsize     = iHandler->getNORBmax();
   const long long fullsize  = ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize ) * ((long long) maxlinsize );
   const bool use_disk_tfo   = (( fullsize > CheMPS2::DMRGSCF_max_mem_eri_tfo ) ? true : false );
   const string tmp_filename = CheMPS2::defaultTMPpath + "/" + CheMPS2::DMRGSCF_eri_storage_name;

   // Allocate 2-body rotation memory
   const int linsize_power4  = (( use_disk_tfo ) ? CheMPS2::DMRGSCF_max_mem_eri_tfo : fullsize );
   const int dmrgsize_power4 = nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG;
   const int work_mem_size   = max( max( linsize_power4 , maxlinsize*maxlinsize*4 ) , dmrgsize_power4 ); //For (ERI rotation, update unitary, block diagonalize, orbital localization)
   double * mem1 = new double[ work_mem_size ];
   double * mem2 = new double[ work_mem_size ];

   //The two-body rotator and Edmiston-Ruedenberg active space localizer
   EdmistonRuedenberg * theLocalizer = NULL;
   if ( scf_options->getWhichActiveSpace() == 2 ){ theLocalizer = new EdmistonRuedenberg( HamDMRG->getVmat(), iHandler->getGroupNumber() ); }

   //Load unitary from disk
   if ( scf_options->getStoreUnitary() ){
      struct stat file_info;
      int master_stat = stat( (scf_options->getUnitaryStorageName()).c_str(), &file_info );
      if ( master_stat == 0 ){ unitary->loadU( scf_options->getUnitaryStorageName() ); }
   }

   //Load DIIS from disk
   DIIS * diis = NULL;
   if (( scf_options->getDoDIIS() ) && ( scf_options->getStoreDIIS() )){
      struct stat file_info;
      int master_stat = stat( (scf_options->getDIISStorageName()).c_str(), &file_info );
      if ( master_stat == 0 ){
         const int diis_vec_size = iHandler->getROTparamsize();
         diis = new DIIS( diis_vec_size, unitary->getNumVariablesX(), scf_options->getNumDIISVecs() );
         diis->loadDIIS( scf_options->getDIISStorageName() );
         diis_vec = new double[ diis_vec_size ];
      }
   }

   int nIterations = 0;

   /*******************************
   ***   Actual DMRGSCF loops   ***
   *******************************/
   while (( gradNorm > scf_options->getGradientThreshold() ) && ( nIterations < scf_options->getMaxIterations() )){

      nIterations++;

      //Update the unitary transformation
      if ( unitary->getNumVariablesX() > 0 ){

         unitary->updateUnitary( mem1, mem2, gradient, true, true ); //multiply = compact = true

         if (( scf_options->getDoDIIS() ) && ( updateNorm <= scf_options->getDIISGradientBranch() )){
            if ( scf_options->getWhichActiveSpace() == 1 ){
               cout << "DMRGSCF::solve : DIIS has started. Active space not rotated to NOs anymore!" << endl;
            }
            if ( scf_options->getWhichActiveSpace() == 2 ){
               cout << "DMRGSCF::solve : DIIS has started. Active space not rotated to localized orbitals anymore!" << endl;
            }
            if ( diis == NULL ){
               const int diis_vec_size = iHandler->getROTparamsize();
               diis = new DIIS( diis_vec_size, unitary->getNumVariablesX(), scf_options->getNumDIISVecs() );
               diis_vec = new double[ diis_vec_size ];
               unitary->makeSureAllBlocksDetOne( mem1, mem2 );
            }
            unitary->getLog( diis_vec, mem1, mem2 );
            diis->appendNew( gradient, diis_vec );
            diis->calculateParam( diis_vec );
            unitary->updateUnitary( mem1, mem2, diis_vec, false, false ); //multiply = compact = false
         }
      }
      if (( scf_options->getStoreUnitary() ) && ( gradNorm != 1.0 )){ unitary->saveU( scf_options->getUnitaryStorageName() ); }
      if (( scf_options->getStoreDIIS() ) && ( updateNorm != 1.0 ) && ( diis != NULL )){ diis->saveDIIS( scf_options->getDIISStorageName() ); }
      int master_diis = (( diis != NULL ) ? 1 : 0 );

      //Fill HamDMRG
      buildQmatOCC();
      buildTmatrix();
      fillConstAndTmatDMRG( HamDMRG );
      if ( use_disk_tfo ){ DMRGSCFVmatRotations::blockwise_disk( VMAT_ORIG, HamDMRG->getVmat(), 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename ); }
      else {               DMRGSCFVmatRotations::full( VMAT_ORIG, HamDMRG->getVmat(), 'A', iHandler, unitary, mem1, mem2 ); }

      //Localize the active space and reorder the orbitals within each irrep based on the exchange matrix
      if (( scf_options->getWhichActiveSpace() == 2 ) && ( master_diis == 0 )){ //When the DIIS has started: stop
         theLocalizer->Optimize(mem1, mem2, scf_options->getStartLocRandom()); //Default EDMISTONRUED_gradThreshold and EDMISTONRUED_maxIter used
         theLocalizer->FiedlerExchange(maxlinsize, mem1, mem2);
         fillLocalizedOrbitalRotations(theLocalizer->getUnitary(), iHandler, mem1);
         unitary->rotateActiveSpaceVectors(mem1, mem2);
         buildQmatOCC(); //With an updated unitary, the Qocc, Tmat, and HamDMRG objects need to be updated as well.
         buildTmatrix();
         fillConstAndTmatDMRG( HamDMRG );
         if ( use_disk_tfo ){ DMRGSCFVmatRotations::blockwise_disk( VMAT_ORIG, HamDMRG->getVmat(), 'A', iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename ); }
         else {               DMRGSCFVmatRotations::full( VMAT_ORIG, HamDMRG->getVmat(), 'A', iHandler, unitary, mem1, mem2 ); }
         cout << "DMRGSCF::solve : Rotated the active space to localized orbitals, sorted according to the exchange matrix." << endl;
      }

      if (( OptScheme == NULL ) && ( rootNum == 1 )){ // Do FCI, and calculate the 2DM

         const int nalpha = ( num_elec + TwoS ) / 2;
         const int nbeta  = ( num_elec - TwoS ) / 2;
         const double workmem = 1000.0; // 1GB
         const int verbose = 2;
         CheMPS2::FCI * theFCI = new CheMPS2::FCI( HamDMRG, nalpha, nbeta, Irrep, workmem, verbose );
         double * inoutput = new double[ theFCI->getVecLength(0) ];
         theFCI->ClearVector( theFCI->getVecLength(0), inoutput );
         inoutput[ theFCI->LowestEnergyDeterminant() ] = 1.0;
         Energy = theFCI->GSDavidson( inoutput );
         theFCI->Fill2RDM( inoutput, DMRG2DM );
         delete theFCI;
         delete [] inoutput;

      } else { //Do the DMRG sweeps, and calculate the 2DM
      
         for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] = 0.0; } //Clear the 2-RDM (to allow for state-averaged calculations)
         DMRG * theDMRG = new DMRG(Prob, OptScheme);
         for (int state = 0; state < rootNum; state++){
            if (state > 0){ theDMRG->newExcitation( fabs( Energy ) ); }
            Energy = theDMRG->Solve();
            if ( scf_options->getStateAveraging() ){ // When SA-DMRGSCF: 2DM += current 2DM
               theDMRG->calc2DMandCorrelations();
               copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
            }
            if ((state == 0) && (rootNum > 1)){ theDMRG->activateExcitations( rootNum-1 ); }
         }
         if ( !( scf_options->getStateAveraging() )){ // When SS-DMRGSCF: 2DM += last 2DM
            theDMRG->calc2DMandCorrelations();
            copy2DMover( theDMRG->get2DM(), nOrbDMRG, DMRG2DM );
         }
         if (scf_options->getDumpCorrelations()){ theDMRG->getCorrelations()->Print(); } // Correlations have been calculated in the loop (SA) or outside of the loop (SS)
         if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
         if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
         delete theDMRG;
         if ((scf_options->getStateAveraging()) && (rootNum > 1)){
            const double averagingfactor = 1.0 / rootNum;
            for ( int cnt = 0; cnt < dmrgsize_power4; cnt++ ){ DMRG2DM[ cnt ] *= averagingfactor; }
         }
         
      }
      setDMRG1DM(num_elec, nOrbDMRG, DMRG1DM, DMRG2DM);

      //Possibly rotate the active space to the natural orbitals
      if (( scf_options->getWhichActiveSpace() == 1 ) && ( master_diis == 0 )){ //When the DIIS has started: stop
         copy_active( DMRG1DM, theQmatWORK, iHandler, true );
         block_diagonalize( 'A', theQmatWORK, unitary, mem1, mem2, iHandler, true, DMRG2DM ); // Unitary is updated and DMRG2DM rotated
         setDMRG1DM( num_elec, nOrbDMRG, DMRG1DM, DMRG2DM );
         buildQmatOCC(); //With an updated unitary, the Qocc and Tmat matrices need to be updated as well.
         buildTmatrix();
         cout << "DMRGSCF::solve : Rotated the active space to natural orbitals, sorted according to the NOON." << endl;
      }

      //Calculate the matrix elements needed to calculate the gradient and hessian
      buildQmatACT();
      if ( use_disk_tfo ){ DMRGSCFVmatRotations::blockwise_disk( VMAT_ORIG, theRotatedTEI, iHandler, unitary, mem1, mem2, work_mem_size, tmp_filename ); }
      else {               DMRGSCFVmatRotations::full( VMAT_ORIG, theRotatedTEI, iHandler, unitary, mem1, mem2 ); }
      buildFmat( theFmatrix, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);
      buildWtilde(wmattilde, theTmatrix, theQmatOCC, theQmatACT, iHandler, theRotatedTEI, DMRG2DM, DMRG1DM);

      //Calculate the gradient, hessian and corresponding update. On return, gradient contains the rescaled gradient == the update.
      augmentedHessianNR(theFmatrix, wmattilde, iHandler, unitary, gradient, &updateNorm, &gradNorm);

   }

   delete [] mem1;
   delete [] mem2;

   delete Prob;
   delete HamDMRG;
   delete [] gradient;
   if ( diis_vec != NULL ){ delete [] diis_vec; }
   if ( diis != NULL ){ delete diis; }
   if ( theLocalizer != NULL ){ delete theLocalizer; }

   return Energy;

}

void CheMPS2::CASSCF::augmentedHessianNR(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * theupdate, double * updateNorm, double * gradNorm){

   /* A good read to understand
            (1) how the augmented Hessian arises from a rational function optimization
            (2) where the parameter lambda in Eq. (22) of Yanai, IJQC 109, 2178-2190 (2009) comes from
            (3) why the smallest algebraic eigenvalue + corresponding eigenvector should be retained for minimizations
      Banerjee, Adams, Simons, Shepard, "Search for stationary points on surfaces",
      J. Phys. Chem. 1985, volume 89, pages 52-57, doi:10.1021/j100247a015  */
   
   //Calculate the gradient
   const int x_linearlength = localUmat->getNumVariablesX();
   gradNorm[0] = calcGradient(localFmat, localIdx, localUmat, theupdate);
   
   //Calculate the Hessian
   int dim = x_linearlength + 1;
   int size = dim * dim;
   double * hessian = new double[size];
   calcHessian(localFmat, localwtilde, localIdx, localUmat, hessian, dim);
   
   //Augment the gradient into the Hessian matrix
   for (int cnt=0; cnt<x_linearlength; cnt++){
      hessian[cnt + x_linearlength*dim] = theupdate[cnt];
      hessian[x_linearlength + dim*cnt] = theupdate[cnt];
   }
   hessian[x_linearlength + dim*x_linearlength] = 0.0;

   //Find the lowest eigenvalue and corresponding eigenvector of the augmented hessian
   {
      Davidson deBoskabouter( dim, CheMPS2::DAVIDSON_NUM_VEC,
                                   CheMPS2::DAVIDSON_NUM_VEC_KEEP,
                                   CheMPS2::DAVIDSON_FCI_RTOL,
                                   CheMPS2::DAVIDSON_PRECOND_CUTOFF, false ); // No debug printing
      double ** whichpointers = new double*[2];

      char instruction = deBoskabouter.FetchInstruction( whichpointers );
      assert( instruction == 'A' );
      for (int cnt = 0; cnt < dim; cnt++){ whichpointers[1][cnt] = hessian[cnt*(1+dim)]; } // Preconditioner = diagonal elements of the augmented Hessian
      for (int cnt = 0; cnt < x_linearlength; cnt++){ // Initial guess = [ -gradient / diag(hessian) , 1 ]
         const double denom = ( whichpointers[1][cnt] > CheMPS2::DAVIDSON_PRECOND_CUTOFF ) ? whichpointers[1][cnt] : CheMPS2::DAVIDSON_PRECOND_CUTOFF;
         whichpointers[0][cnt] = - theupdate[cnt] / denom;
      }
      whichpointers[0][x_linearlength] = 1.0;

      instruction = deBoskabouter.FetchInstruction( whichpointers );
      while ( instruction == 'B' ){
         char notrans = 'N';
         int one = 1;
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_(&notrans, &notrans, &dim, &one, &dim, &alpha, hessian, &dim, whichpointers[0], &dim, &beta, whichpointers[1], &dim);
         instruction = deBoskabouter.FetchInstruction( whichpointers );
      }

      assert( instruction == 'C' );
      double scalar = 1.0 / whichpointers[0][x_linearlength];
      cout << "DMRGSCF::augmentedHessianNR : Augmented Hessian update found with " << deBoskabouter.GetNumMultiplications() << " Davidson iterations." << endl;
      if (CheMPS2::DMRGSCF_debugPrint){
         cout << "DMRGSCF::augmentedHessianNR : Lowest eigenvalue = " << whichpointers[1][0] << endl;
         cout << "DMRGSCF::augmentedHessianNR : The last number of the eigenvector (which will be rescaled to one) = " << scalar << endl;
      }
      for (int cnt = 0; cnt < x_linearlength; cnt++){ theupdate[cnt] = scalar * whichpointers[0][cnt]; }
      delete [] whichpointers;
   }
   
   //Calculate the update norm
   updateNorm[0] = 0.0;
   for (int cnt = 0; cnt < x_linearlength; cnt++){ updateNorm[0] += theupdate[cnt] * theupdate[cnt]; }
   updateNorm[0] = sqrt(updateNorm[0]);
   cout << "DMRGSCF::augmentedHessianNR : Norm of the update = " << updateNorm[0] << endl;
   
   delete [] hessian;

}

double CheMPS2::CASSCF::calcGradient(const DMRGSCFmatrix * localFmat, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * gradient){

   for (int cnt=0; cnt<localUmat->getNumVariablesX(); cnt++){
      const int index1 = localUmat->getFirstIndex(cnt);
      const int index2 = localUmat->getSecondIndex(cnt);
      // irrep1 == irrep2 due to construction DMRGSCFunitary
      const int irrep  = localIdx->getOrbitalIrrep( index1 );
      const int shift  = localIdx->getOrigNOCCstart( irrep );
      const int relIndex1 = index1 - shift;
      const int relIndex2 = index2 - shift;
      gradient[cnt] = 2 * ( localFmat->get( irrep, relIndex1, relIndex2 ) - localFmat->get( irrep, relIndex2, relIndex1 ) );
   }
   
   double gradNorm = 0.0;
   for (int cnt=0; cnt<localUmat->getNumVariablesX(); cnt++){ gradNorm += gradient[cnt] * gradient[cnt]; }
   gradNorm = sqrt(gradNorm);
   cout << "DMRGSCF::calcGradient : Norm of the gradient = " << gradNorm << endl;
   return gradNorm;

}

void CheMPS2::CASSCF::calcHessian(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const DMRGSCFunitary * localUmat, double * hessian, const int rowjump){

   const int lindim = ( localUmat->getNumVariablesX() * ( localUmat->getNumVariablesX() + 1 ))/2;
   
   #pragma omp parallel for schedule(static)
   for (int count=0; count<lindim; count++){
   
      int col = 1;
      while ( (col*(col+1))/2 <= count ) col++;
      col -= 1;
      int row = count - (col*(col+1))/2;
      
      const int p_index = localUmat->getFirstIndex(row);
      const int q_index = localUmat->getSecondIndex(row);
      // irrep_p == irrep_q due to construction DMRGSCFunitary
      const int irrep_pq = localIdx->getOrbitalIrrep( p_index );
      
      const int r_index = localUmat->getFirstIndex(col);
      const int s_index = localUmat->getSecondIndex(col);
      // irrep_r == irrep_s due to construction DMRGSCFunitary
      const int irrep_rs = localIdx->getOrbitalIrrep( r_index );
      
      const int rel_p_index = p_index - localIdx->getOrigNOCCstart( irrep_pq );
      const int rel_q_index = q_index - localIdx->getOrigNOCCstart( irrep_pq );
      const int rel_r_index = r_index - localIdx->getOrigNOCCstart( irrep_rs );
      const int rel_s_index = s_index - localIdx->getOrigNOCCstart( irrep_rs );
      
      hessian[row + rowjump * col] = Wmat(localFmat, localwtilde, localIdx, irrep_pq, irrep_rs, rel_p_index, rel_q_index, rel_r_index, rel_s_index)
                                   - Wmat(localFmat, localwtilde, localIdx, irrep_pq, irrep_rs, rel_q_index, rel_p_index, rel_r_index, rel_s_index)
                                   - Wmat(localFmat, localwtilde, localIdx, irrep_pq, irrep_rs, rel_p_index, rel_q_index, rel_s_index, rel_r_index)
                                   + Wmat(localFmat, localwtilde, localIdx, irrep_pq, irrep_rs, rel_q_index, rel_p_index, rel_s_index, rel_r_index);
      hessian[col + rowjump * row] = hessian[row + rowjump * col];
      
   }

}

double CheMPS2::CASSCF::Wmat(const DMRGSCFmatrix * localFmat, const DMRGSCFwtilde * localwtilde, const DMRGSCFindices * localIdx, const int irrep_pq, const int irrep_rs, const int relindexP, const int relindexQ, const int relindexR, const int relindexS){

   double value = 0.0;
   
   if ( ( irrep_pq == irrep_rs ) && ( relindexQ == relindexR ) ) {
      value = localFmat->get(irrep_pq, relindexP, relindexS) + localFmat->get(irrep_pq, relindexS, relindexP);
   }
   
   if (relindexP >= localIdx->getNOCC(irrep_pq) + localIdx->getNDMRG(irrep_pq)){ return value; }
   if (relindexR >= localIdx->getNOCC(irrep_rs) + localIdx->getNDMRG(irrep_rs)){ return value; } //index1 and index3 are now certainly not virtual!
   
   value += localwtilde->get( irrep_pq, irrep_rs, relindexP, relindexQ, relindexR, relindexS );
   return value;

}

void CheMPS2::CASSCF::buildWtilde(DMRGSCFwtilde * localwtilde, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM){

   localwtilde->clear();
   const int numIrreps  = localIdx->getNirreps();
   const int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   for (int irrep_pq = 0; irrep_pq < numIrreps; irrep_pq++){
   
      const int NumOCCpq  = localIdx->getNOCC(  irrep_pq );
      const int NumDMRGpq = localIdx->getNDMRG( irrep_pq );
      const int NumORBpq  = localIdx->getNORB(  irrep_pq );
      const int NumCOREpq = NumOCCpq + NumDMRGpq;
      
      //If irrep_pq == irrep_rs and P == R occupied --> QS only active or virtual
      #pragma omp parallel for schedule(static)
      for (int relindexP = 0; relindexP < NumOCCpq; relindexP++){
         double * subblock = localwtilde->getBlock( irrep_pq, irrep_pq, relindexP, relindexP);
         for (int relindexS = NumOCCpq; relindexS < NumORBpq; relindexS++){
            for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 4 * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                   + localJKocc->get(irrep_pq, relindexQ, relindexS)
                                                                   + localJKact->get(irrep_pq, relindexQ, relindexS) );
            }
         }
      }
      
      //If irrep_pq == irrep_rs and P,R active --> QS only occupied or virtual
      #pragma omp parallel for schedule(static)
      for (int combined = 0; combined < NumDMRGpq*NumDMRGpq; combined++){
         
         const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
         const int relindexR  = NumOCCpq + ( combined / NumDMRGpq );
         double * subblock    = localwtilde->getBlock( irrep_pq, irrep_pq, relindexP, relindexR );
         const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
         const int DMRGindexR = relindexR - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
         const double OneDMvalue = local1DM[ DMRGindexP + totOrbDMRG * DMRGindexR ];
         
         for (int relindexS = 0; relindexS < NumOCCpq; relindexS++){
            for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
            for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
         }
         for (int relindexS = NumCOREpq; relindexS < NumORBpq; relindexS++){
            for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
            for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
               subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue * ( localTmat->get( irrep_pq, relindexQ, relindexS)
                                                                                + localJKocc->get(irrep_pq, relindexQ, relindexS) );
            }
         }
      }
   
      for (int irrep_rs = 0; irrep_rs < numIrreps; irrep_rs++){
      
         const int NumOCCrs  = localIdx->getNOCC(  irrep_rs );
         const int NumDMRGrs = localIdx->getNDMRG( irrep_rs );
         const int NumORBrs  = localIdx->getNORB(  irrep_rs );
         const int NumCORErs = NumOCCrs + NumDMRGrs;
         const int productirrep = Irreps::directProd( irrep_pq, irrep_rs );
      
         // P and R occupied --> QS only active or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumOCCpq*NumOCCrs; combined++){
         
            const int relindexP = combined % NumOCCpq;
            const int relindexR = combined / NumOCCpq;
            double * subblock   = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR);
            
            for (int relindexS = NumOCCrs; relindexS < NumORBrs; relindexS++){
               for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
                  subblock[ relindexQ + NumORBpq * relindexS ] +=  
                      4 * ( 4 * theInts->FourIndexAPI(irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relindexR)
                              - theInts->FourIndexAPI(irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relindexR)
                              - theInts->FourIndexAPI(irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relindexP) );
               }
            }
         } // End combined P and R occupied
         
         // P and R active --> QS only occupied or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumDMRGpq*NumDMRGrs; combined++){
         
            const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
            const int relindexR  = NumOCCrs + ( combined / NumDMRGpq );
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
            const int DMRGindexR = relindexR - NumOCCrs + localIdx->getDMRGcumulative( irrep_rs );
            
            for (int irrep_alpha = 0; irrep_alpha < numIrreps; irrep_alpha++){
               
               const int irrep_beta   = Irreps::directProd( irrep_alpha, productirrep );
               const int NumDMRGalpha = localIdx->getNDMRG( irrep_alpha );
               const int NumDMRGbeta  = localIdx->getNDMRG( irrep_beta  );
               
               for (int alpha = 0; alpha < NumDMRGalpha; alpha++){
               
                  const int DMRGalpha = localIdx->getDMRGcumulative( irrep_alpha ) + alpha;
                  const int relalpha  = localIdx->getNOCC( irrep_alpha ) + alpha;
                  
                  for (int beta = 0; beta < NumDMRGbeta; beta++){
                  
                     const int DMRGbeta = localIdx->getDMRGcumulative( irrep_beta ) + beta;
                     const int relbeta  = localIdx->getNOCC( irrep_beta ) + beta;
                     
                     const double TwoDMvalue1  = local2DM[ DMRGindexR + totOrbDMRG * ( DMRGalpha + totOrbDMRG * ( DMRGindexP + totOrbDMRG * DMRGbeta ) ) ];
                     const double TwoDMvalue23 = local2DM[ DMRGindexR + totOrbDMRG * ( DMRGalpha + totOrbDMRG * ( DMRGbeta + totOrbDMRG * DMRGindexP ) ) ]
                                               + local2DM[ DMRGindexR + totOrbDMRG * ( DMRGindexP + totOrbDMRG * ( DMRGbeta + totOrbDMRG * DMRGalpha ) ) ];
                     
                     for (int relindexS = 0; relindexS < NumOCCrs; relindexS++){
                        for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                        for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                     }
                     for (int relindexS = NumCORErs; relindexS < NumORBrs; relindexS++){
                        for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                        for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                           subblock[ relindexQ + NumORBpq * relindexS ] +=
                              2 * ( TwoDMvalue1  * theInts->FourIndexAPI( irrep_pq, irrep_alpha, irrep_rs, irrep_beta, relindexQ, relalpha, relindexS, relbeta)
                                  + TwoDMvalue23 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_alpha, irrep_beta, relindexQ, relindexS, relalpha, relbeta) );
                        }
                     }
                  }
               }
            }
         } // End combined P and R active
         
         // P active and R occupied  -->  Q occupied or virtual  //  S active or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumDMRGpq*NumOCCrs; combined++){
            
            const int relindexP  = NumOCCpq + ( combined % NumDMRGpq );
            const int relindexR  = combined / NumDMRGpq;
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexP = relindexP - NumOCCpq + localIdx->getDMRGcumulative( irrep_pq );
            
            for (int alpha = 0; alpha < NumDMRGpq; alpha++){
            
               const int DMRGalpha     = localIdx->getDMRGcumulative( irrep_pq ) + alpha;
               const int relalpha      = localIdx->getNOCC( irrep_pq ) + alpha;
               const double OneDMvalue = local1DM[ DMRGalpha + totOrbDMRG * DMRGindexP ];
               
               for (int relindexS = NumOCCrs; relindexS < NumORBrs; relindexS++){
                  for (int relindexQ = 0; relindexQ < NumOCCpq; relindexQ++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relalpha, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relalpha, relindexS, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relalpha) );
                  }
                  for (int relindexQ = NumCOREpq; relindexQ < NumORBpq; relindexQ++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relalpha, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relalpha, relindexS, relindexR)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relindexR, relalpha) );
                  }
               }
            }
         } // End combined P active and R occupied
         
         // P occupied and R active  -->  Q active or virtual  //  S occupied or virtual
         #pragma omp parallel for schedule(static)
         for (int combined = 0; combined < NumOCCpq*NumDMRGrs; combined++){
            
            const int relindexP  = combined % NumOCCpq;
            const int relindexR  = NumOCCrs + ( combined / NumOCCpq );
            double * subblock    = localwtilde->getBlock( irrep_pq, irrep_rs, relindexP, relindexR );
            const int DMRGindexR = relindexR - NumOCCrs + localIdx->getDMRGcumulative( irrep_rs );
            
            for (int beta = 0; beta < NumDMRGrs; beta++){
            
               const int DMRGbeta      = localIdx->getDMRGcumulative( irrep_rs ) + beta;
               const int relbeta       = localIdx->getNOCC( irrep_rs ) + beta;
               const double OneDMvalue = local1DM[ DMRGindexR + totOrbDMRG * DMRGbeta ];
               
               for (int relindexQ = NumOCCpq; relindexQ < NumORBpq; relindexQ++){
                  for (int relindexS = 0; relindexS < NumOCCrs; relindexS++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relbeta, relindexP) );
                  }
                  for (int relindexS = NumCORErs; relindexS < NumORBrs; relindexS++){
                     subblock[ relindexQ + NumORBpq * relindexS ] += 2 * OneDMvalue *
                         ( 4 * theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_pq, irrep_rs, relindexQ, relindexS, relindexP, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_pq, irrep_rs, irrep_rs, relindexQ, relindexP, relindexS, relbeta)
                             - theInts->FourIndexAPI( irrep_pq, irrep_rs, irrep_rs, irrep_pq, relindexQ, relindexS, relbeta, relindexP) );
                  }
               }
            }
         } // End combined P occupied and R active
         
      }
   }

}

void CheMPS2::CASSCF::buildFmat(DMRGSCFmatrix * localFmat, const DMRGSCFmatrix * localTmat, const DMRGSCFmatrix * localJKocc, const DMRGSCFmatrix * localJKact, const DMRGSCFindices * localIdx, const DMRGSCFintegrals * theInts, double * local2DM, double * local1DM){
   
   localFmat->clear();
   const int numIrreps  = localIdx->getNirreps();
   const int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   for (int irrep_pq = 0; irrep_pq < numIrreps; irrep_pq++){
   
      const int NumORB  = localIdx->getNORB(  irrep_pq );
      const int NumOCC  = localIdx->getNOCC(  irrep_pq );
      const int NumDMRG = localIdx->getNDMRG( irrep_pq );
      const int NumOCCDMRG = NumOCC + NumDMRG;
      
      #pragma omp parallel for schedule(static)
      for (int p = 0; p < NumOCC; p++){
         for (int q = 0; q < NumORB; q++){
            localFmat->set( irrep_pq, p, q, 2 * ( localTmat->get(  irrep_pq, q, p )
                                                + localJKocc->get( irrep_pq, q, p )
                                                + localJKact->get( irrep_pq, q, p ) ) );
         }
      }
      
      #pragma omp parallel for schedule(static)
      for (int p = NumOCC; p < NumOCCDMRG; p++){
         const int DMRGindex_p = p - NumOCC + localIdx->getDMRGcumulative( irrep_pq );
         
         //One-body terms --> matrix multiplication?
         for (int r = NumOCC; r < NumOCCDMRG; r++){
            const double OneDMvalue = local1DM[ DMRGindex_p + totOrbDMRG * ( DMRGindex_p + r - p ) ];
            for (int q = 0; q < NumORB; q++){
               localFmat->getBlock(irrep_pq)[ p + NumORB * q ] += OneDMvalue * ( localTmat->get( irrep_pq, q, r ) + localJKocc->get( irrep_pq, q, r) );
            }
         }
         
         //Two-body terms --> matrix multiplication possible?
         for (int irrep_r = 0; irrep_r < numIrreps; irrep_r++){
            const int irrep_product = Irreps::directProd(irrep_pq, irrep_r);
            for (int irrep_s = 0; irrep_s < numIrreps; irrep_s++){
               const int irrep_t = Irreps::directProd(irrep_product, irrep_s);
               for (int r = localIdx->getNOCC(irrep_r); r < localIdx->getNOCC(irrep_r) + localIdx->getNDMRG(irrep_r); r++){
                  const int DMRGindex_r = r - localIdx->getNOCC( irrep_r ) + localIdx->getDMRGcumulative( irrep_r );
                  for (int s = localIdx->getNOCC(irrep_s); s < localIdx->getNOCC(irrep_s) + localIdx->getNDMRG(irrep_s); s++){
                     const int DMRGindex_s = s - localIdx->getNOCC( irrep_s ) + localIdx->getDMRGcumulative( irrep_s );
                     for (int t = localIdx->getNOCC(irrep_t); t < localIdx->getNOCC(irrep_t) + localIdx->getNDMRG(irrep_t); t++){
                        const int DMRGindex_t = t - localIdx->getNOCC( irrep_t ) + localIdx->getDMRGcumulative( irrep_t );
                        const double TwoDMvalue = local2DM[ DMRGindex_p + totOrbDMRG * ( DMRGindex_r + totOrbDMRG * ( DMRGindex_s + totOrbDMRG * DMRGindex_t ) ) ];
                        for (int q = 0; q < NumORB; q++){
                           localFmat->getBlock(irrep_pq)[ p + NumORB * q ] += TwoDMvalue * theInts->FourIndexAPI(irrep_pq, irrep_r, irrep_s, irrep_t, q, r, s, t);
                        }
                     }
                  }
               }
            }
         }
      }
   }

}


