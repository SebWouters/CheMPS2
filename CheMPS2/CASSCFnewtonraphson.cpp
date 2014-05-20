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

#include <stdlib.h>
#include <iostream>
#include <string>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>

#include "CASSCF.h"
#include "Lapack.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;
using std::max;

double CheMPS2::CASSCF::doCASSCFnewtonraphson(const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum){

   //Convergence variables
   double gradNorm = 1.0;
   double Energy;
   
   //The CheMPS2::Problem for the inner DMRG calculation
   Hamiltonian * HamDMRG = new Hamiltonian(nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital());
   int N = Nelectrons;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){ N -= 2*iHandler->getNOCC(irrep); }
   Problem * Prob = new Problem(HamDMRG, TwoS, N, Irrep);
   Prob->SetupReorderD2h(); //Doesn't matter if the group isn't D2h, Prob checks it.
   
   //Determine the maximum NORB(irrep); and the maximum NORB(irrep) which is OK according to the cutoff.
   int maxlinsize   = 0;
   int maxlinsizeOK = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      const int linsize_irrep = iHandler->getNORB(irrep);
      if  (linsize_irrep > maxlinsize  )                                                        {  maxlinsize  = linsize_irrep; }
      if ((linsize_irrep > maxlinsizeOK) && (linsize_irrep <= CheMPS2::CASSCF_maxlinsizeCutoff)){ maxlinsizeOK = linsize_irrep; }
   }
   const bool doBlockWise = (maxlinsize <= CheMPS2::CASSCF_maxlinsizeCutoff) ? false : true; //Only if bigger, do we want to work blockwise
   
   //Determine the blocksize for the 2-body transformation
   int maxBlockSize = maxlinsize;
   if (doBlockWise){
      int factor   = (int) (ceil( (1.0 * maxlinsize) / CheMPS2::CASSCF_maxlinsizeCutoff ) + 0.01);
      maxBlockSize = max( (int) (ceil( (1.0 * maxlinsize) / factor ) + 0.01) , maxlinsizeOK ); //If a particular index can be rotated at once....
      cout << "DMRGSCF info: The max. # orb per irrep           = " << maxlinsize << endl;
      cout << "              The size cutoff for 2-body transfo = " << CheMPS2::CASSCF_maxlinsizeCutoff << endl;
      cout << "              The max. # orb per irrep <= cutoff = " << maxlinsizeOK << endl;
      cout << "              The blocksize for piecewise tfo    = " << maxBlockSize << endl;
   }
   
   //Allocate 2-body rotation memory: One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
   const int maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!
   const int sizeWorkmem1 = max( max( maxBSpower4 , maxlinsize*maxlinsize*3 ) , nOrbDMRG * nOrbDMRG ); //For (2-body tfo , updateUnitary, calcNOON)
   const int sizeWorkmem2 = max( max( maxBSpower4 , maxlinsize*maxlinsize*2 ) , (CheMPS2::CASSCF_rotate2DMtoNO) ? nOrbDMRG*nOrbDMRG*nOrbDMRG*nOrbDMRG : nOrbDMRG*(nOrbDMRG + 1) ); //For (2-body tfo, updateUnitary and rotateUnitaryNOeigenvecs, rotate2DMand1DM or calcNOON)
   double * mem1 = new double[sizeWorkmem1];
   double * mem2 = new double[sizeWorkmem2];
   double * mem3 = NULL;
   if (doBlockWise){ mem3 = new double[maxBSpower4]; }
   
   //Load unitary from disk
   if (CheMPS2::CASSCF_storeUnitary){
      struct stat stFileInfo;
      int intStat = stat(CheMPS2::CASSCF_unitaryStorageName.c_str(),&stFileInfo);
      if (intStat==0){ unitary->loadU(); }
   }

   /*******************************
   ***   Actual DMRGSCF loops   ***
   *******************************/
   while (gradNorm > CheMPS2::CASSCF_gradientNormThreshold){
   
      //Update the unitary transformation based on the previous unitary transformation and the xmatrix
      if (unitary->getNumVariablesX() > 0){ unitary->updateUnitary(mem1, mem2); }
      if ((CheMPS2::CASSCF_storeUnitary) && (gradNorm!=1.0)){ unitary->saveU(); }
   
      //Fill HamDMRG
      buildQmatrixOCC();
      buildOneBodyMatrixElements();
      if (doBlockWise){ fillHamDMRGBlockWise(HamDMRG, mem1, mem2, mem3, maxBlockSize); }
      else {            fillHamDMRG(HamDMRG, mem1, mem2); }
      
      //Do the DMRG sweeps, and calculate the 2DM
      DMRG * theDMRG = new DMRG(Prob,OptScheme);
      Energy = theDMRG->Solve();
      if (rootNum>1){
         theDMRG->activateExcitations(rootNum-1);
         for (int exc=0; exc<rootNum-1; exc++){
            theDMRG->newExcitation(fabs(Energy));
            Energy = theDMRG->Solve();
         }
      }
      theDMRG->calc2DM();
      copy2DMover(theDMRG->get2DM());
      setDMRG1DM(N);
      
      //Calculate the NOON and possibly rotate the active space to the natural orbitals
      calcNOON(mem1, mem2);
      if (CheMPS2::CASSCF_rotate2DMtoNO){
         rotate2DMand1DM(N, mem1, mem2);
         unitary->rotateUnitaryNOeigenvecs(mem1, mem2);
         buildQmatrixOCC(); //With an updated unitary, the Qocc and Tmat matrices need to be updated as well.
         buildOneBodyMatrixElements();
      }
      
      //Calculate the matrix elements needed to calculate the gradient and hessian
      buildQmatrixACT();
      if (doBlockWise){ fillVmatRotatedBlockWise(mem1, mem2, mem3, maxBlockSize, true); }
      else{             fillVmatRotated(mem1, mem2); }
      buildFmat();

      //Calculate the gradient, hessian and corresponding update
      gradNorm = updateXmatrixAugmentedHessianNR();
      
      //Print the coefficients to discern the 1Ag and 3B1u states in the carbon dimer
      /*PrintCoeff_C2(theDMRG);*/
      
      if (CheMPS2::DMRG_storeMpsOnDisk){        theDMRG->deleteStoredMPS();       }
      if (CheMPS2::DMRG_storeRenormOptrOnDisk){ theDMRG->deleteStoredOperators(); }
      delete theDMRG;
   
   }
   
   delete [] mem1;
   delete [] mem2;
   if (doBlockWise){ delete [] mem3; }
   
   delete Prob;
   delete HamDMRG;
   
   return Energy;

}

double CheMPS2::CASSCF::updateXmatrixAugmentedHessianNR(){

   /* A good read to understand
            (1) how the augmented Hessian arises from a rational function optimization
            (2) where the parameter lambda in Eq. (22) of Yanai, IJQC 109, 2178-2190 (2009) comes from
            (3) why the smallest algebraic eigenvalue + corresponding eigenvector should be retained for minimizations
      Banerjee, Adams, Simons, Shepard, "Search for stationary points on surfaces",
      J. Phys. Chem. 1985, volume 89, pages 52-57, doi:10.1021/j100247a015  */
   
   //Calculate the gradient
   int x_linearlength = unitary->getNumVariablesX();
   double * gradient = new double[x_linearlength];
   double gradNorm = calcGradient(gradient);
   
   //Calculate the Hessian
   int aug_linlength = x_linearlength+1;
   int size = aug_linlength * aug_linlength;
   double * hessian = new double[size];
   calcHessian(hessian, x_linearlength+1);
   
   //Augment the gradient into the Hessian matrix
   for (int cnt=0; cnt<x_linearlength; cnt++){
      hessian[cnt + x_linearlength*aug_linlength] = gradient[cnt];
      hessian[x_linearlength + aug_linlength*cnt] = gradient[cnt];
   }
   hessian[x_linearlength + aug_linlength*x_linearlength] = 0.0;
   
   //Find its lowest eigenvalue and vector
   double * work = new double[size];
   double * eigen = new double[aug_linlength];
   char jobz = 'V';
   char uplo = 'U';
   int info;
   dsyev_(&jobz,&uplo,&aug_linlength,hessian,&aug_linlength,eigen,work,&size,&info);
   
   if (CheMPS2::CASSCF_debugPrint){
      cout << "Lowest eigenvalue = " << eigen[0] << endl;
      cout << "The last number of the eigenvector (which will be rescaled to one) = " << hessian[x_linearlength] << endl;
   }
   double scalar = 1.0/hessian[x_linearlength];
   int inc = 1;
   dscal_(&x_linearlength,&scalar,hessian,&inc);
   
   //Copy the new x back to xmatrix.
   unitary->copyXsolutionBack(hessian);
   
   delete [] hessian;
   delete [] gradient;
   delete [] eigen;
   delete [] work;
   
   return gradNorm;

}

double CheMPS2::CASSCF::updateXmatrixNewtonRaphson(){
   
   //Calculate the gradient
   int x_linearlength = unitary->getNumVariablesX();
   double * gradient = new double[x_linearlength];
   double gradNorm = calcGradient(gradient);
   
   //Calculate the Hessian
   int size = x_linearlength * x_linearlength;
   double * hessian = new double[size];
   calcHessian(hessian, x_linearlength);
   
   //Invert the Hessian
   double * inverse = new double[size];
   double * vector = new double[x_linearlength];
   char jobz = 'V';
   char uplo = 'U';
   int info;
   dsyev_(&jobz,&uplo,&x_linearlength,hessian,&x_linearlength,vector,inverse,&size,&info);
   
   if (vector[0]<=0.0){
      cout << "DMRGSCF :: Eigenvalues of the Hessian = [ ";
      for (int cnt=0; cnt<x_linearlength-1; cnt++){ cout << vector[cnt] << " , "; }
      cout << vector[x_linearlength-1] << " ]" << endl;
   }
   
   for (int cnt=0; cnt<x_linearlength; cnt++){
      double value = 1.0/sqrt(vector[cnt]);
      for (int cnt2=0; cnt2<x_linearlength; cnt2++){
         hessian[cnt2 + x_linearlength*cnt] *= value;
      }
   }
   char notr = 'N';
   char tran = 'T';
   double afac = 1.0;
   double bfac = 0.0; //set
   dgemm_(&notr,&tran,&x_linearlength,&x_linearlength,&x_linearlength,&afac,hessian,&x_linearlength,hessian,&x_linearlength,&bfac,inverse,&x_linearlength);
   
   //Calculate the new x --> Eq. (6c) of the Siegbahn paper
   afac = -1.0;
   int one = 1;
   dgemm_(&notr,&notr,&x_linearlength,&one,&x_linearlength,&afac,inverse,&x_linearlength,gradient,&x_linearlength,&bfac,vector,&x_linearlength);
   
   //Copy the new x back to xmatrix.
   unitary->copyXsolutionBack(vector);
   
   delete [] hessian;
   delete [] gradient;
   delete [] vector;
   delete [] inverse;
   
   return gradNorm;

}

double CheMPS2::CASSCF::calcGradient(double * gradient){

   #pragma omp parallel for schedule(static)
   for (int cnt=0; cnt<unitary->getNumVariablesX(); cnt++){
      gradient[cnt] = 2*( Fmat( unitary->getFirstIndex(cnt), unitary->getSecondIndex(cnt) ) - Fmat( unitary->getSecondIndex(cnt), unitary->getFirstIndex(cnt) ));
   }
   
   double gradNorm = 0.0;
   for (int cnt=0; cnt<unitary->getNumVariablesX(); cnt++){ gradNorm += gradient[cnt] * gradient[cnt]; }
   gradNorm = sqrt(gradNorm);
   cout << "DMRGSCF :: Norm of the gradient = " << gradNorm << endl;
   return gradNorm;

}

void CheMPS2::CASSCF::calcHessian(double * hessian, const int rowjump){

   const int lindim = ( unitary->getNumVariablesX() * ( unitary->getNumVariablesX() + 1 ))/2;
   
   #pragma omp parallel for schedule(static)
   for (int count=0; count<lindim; count++){
   
      int col = 1;
      while ( (col*(col+1))/2 <= count ) col++;
      col -= 1;
      int row = count - (col*(col+1))/2;
      
      int p_index = unitary->getFirstIndex(row);
      int q_index = unitary->getSecondIndex(row);
      int r_index = unitary->getFirstIndex(col);
      int s_index = unitary->getSecondIndex(col);
      hessian[row + rowjump * col] = Wmat(p_index,q_index,r_index,s_index)
                                   - Wmat(q_index,p_index,r_index,s_index)
                                   - Wmat(p_index,q_index,s_index,r_index)
                                   + Wmat(q_index,p_index,s_index,r_index);
      hessian[col + rowjump * row] = hessian[row + rowjump * col];
      
   }

}

double CheMPS2::CASSCF::Wmat(const int index1, const int index2, const int index3, const int index4) const{

   const int irrep1 = HamOrig->getOrbitalIrrep(index1);
   const int irrep2 = HamOrig->getOrbitalIrrep(index2);
   if (irrep1 != irrep2){ return 0.0; } //From now on: irrep1 == irrep2
   
   const int irrep3 = HamOrig->getOrbitalIrrep(index3);
   const int irrep4 = HamOrig->getOrbitalIrrep(index4);
   if (irrep3 != irrep4){ return 0.0; } //From now on: irrep3 == irrep4
   
   double value = 0.0;
   
   if ( (irrep1 == irrep3) && (index2 == index3) ) { //irrep1 == irrep2 == irrep3 == irrep4
      value = Fmat(index1,index4) + Fmat(index4,index1);
   }
   
   if (index1 >= iHandler->getOrigNVIRTstart(irrep1)){ return value; }
   if (index3 >= iHandler->getOrigNVIRTstart(irrep3)){ return value; } //index1 and index3 are now certainly not virtual!
   
   const int relIndex1 = index1 - iHandler->getOrigNOCCstart(irrep1);
   const int relIndex2 = index2 - iHandler->getOrigNOCCstart(irrep2);
   const int relIndex3 = index3 - iHandler->getOrigNOCCstart(irrep3);
   const int relIndex4 = index4 - iHandler->getOrigNOCCstart(irrep4);
   
   if (index1 < iHandler->getOrigNDMRGstart(irrep1)){
      if (index3 < iHandler->getOrigNDMRGstart(irrep3)){
      
         // (index1,index3) (occupied,occupied) --> (alpha,beta) can be (occupied,occupied) or (active,active)
         if (index1==index3){
         
            // Part of one-body matrix elements if 1-RDM occ indices
            value += 4 * ( TmatRotated(index2,index4) + QmatOCC(index2,index4) + QmatACT(index2, index4) );
            value += 12 * VmatRotated->get(irrep2, irrep4, irrep1, irrep1, relIndex2, relIndex4, relIndex1, relIndex1)
                   - 4  * VmatRotated->get(irrep2, irrep1, irrep4, irrep1, relIndex2, relIndex1, relIndex4, relIndex1);
            return value;
            
         } else { //index1 != index3
         
            value += 16 * VmatRotated->get(irrep2, irrep4, irrep1, irrep3, relIndex2, relIndex4, relIndex1, relIndex3)
                   - 4  * VmatRotated->get(irrep2, irrep1, irrep4, irrep3, relIndex2, relIndex1, relIndex4, relIndex3)
                   - 4  * VmatRotated->get(irrep2, irrep4, irrep3, irrep1, relIndex2, relIndex4, relIndex3, relIndex1);
            return value;
            
         }
      
      } else {
      
         const int DMRGindex3 = iHandler->getDMRGcumulative(irrep3) + index3 - iHandler->getOrigNDMRGstart(irrep3);
      
         // (index1,index3) (occupied,active) --> (alpha,beta) can be (active,occupied) or (occupied,active)
         for (int alpha_index=iHandler->getOrigNDMRGstart(irrep3); alpha_index<iHandler->getOrigNVIRTstart(irrep3); alpha_index++){
            const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep3) + alpha_index - iHandler->getOrigNDMRGstart(irrep3);
            const int relALPHA = alpha_index - iHandler->getOrigNOCCstart(irrep3);
            value += 2 * DMRG1DM[ DMRGindex3 + nOrbDMRG * DMRGindexALPHA ]
                          * ( 4 * VmatRotated->get(irrep2, irrep4, irrep1, irrep3, relIndex2, relIndex4, relIndex1, relALPHA)
                                - VmatRotated->get(irrep2, irrep4, irrep3, irrep1, relIndex2, relIndex4, relALPHA,  relIndex1)
                                - VmatRotated->get(irrep2, irrep1, irrep4, irrep3, relIndex2, relIndex1, relIndex4, relALPHA) );
         }
         return value;
         
      }
   } else {
      
      const int DMRGindex1 = iHandler->getDMRGcumulative(irrep1) + index1 - iHandler->getOrigNDMRGstart(irrep1);
   
      if (index3 < iHandler->getOrigNDMRGstart(irrep3)){
      
         // (index1,index3) (active,occupied) --> (alpha,beta) can be (active,occupied) or (occupied,active)
         for (int alpha_index=iHandler->getOrigNDMRGstart(irrep1); alpha_index<iHandler->getOrigNVIRTstart(irrep1); alpha_index++){
            const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep1) + alpha_index - iHandler->getOrigNDMRGstart(irrep1);
            const int relALPHA = alpha_index - iHandler->getOrigNOCCstart(irrep1);
            value += 2 * DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindexALPHA ]
                          * ( 4 * VmatRotated->get(irrep2, irrep4, irrep1, irrep3, relIndex2, relIndex4, relALPHA,  relIndex3)
                                - VmatRotated->get(irrep2, irrep1, irrep4, irrep3, relIndex2, relALPHA,  relIndex4, relIndex3)
                                - VmatRotated->get(irrep2, irrep4, irrep3, irrep1, relIndex2, relIndex4, relIndex3, relALPHA) );
         }
         return value;
         
      } else {
      
         const int DMRGindex3 = iHandler->getDMRGcumulative(irrep3) + index3 - iHandler->getOrigNDMRGstart(irrep3);
      
         // (index1,index3) (active,active) --> (alpha,beta) can be (occupied,occupied) or (active,active)
         // Case1: (alpha,beta)==(occ,occ) --> alpha == beta
         // Part of one-body matrix elements if 1-RDM active indices
         value += 2 * DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindex3 ] * ( TmatRotated(index2,index4) + QmatOCC(index2,index4) );
         
         // Case2: (alpha,beta)==(act,act)
         const int productIrrep = SymmInfo.directProd(irrep1,irrep3);
         for (int irrep_alpha=0; irrep_alpha<numberOfIrreps; irrep_alpha++){
            int irrep_beta = SymmInfo.directProd(productIrrep,irrep_alpha);
            for (int alpha_index=iHandler->getOrigNDMRGstart(irrep_alpha); alpha_index<iHandler->getOrigNVIRTstart(irrep_alpha); alpha_index++){
               const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep_alpha) + alpha_index - iHandler->getOrigNDMRGstart(irrep_alpha);
               const int relALPHA = alpha_index - iHandler->getOrigNOCCstart(irrep_alpha);
               for (int beta_index=iHandler->getOrigNDMRGstart(irrep_beta); beta_index<iHandler->getOrigNVIRTstart(irrep_beta); beta_index++){
                  const int DMRGindexBETA = iHandler->getDMRGcumulative(irrep_beta) + beta_index - iHandler->getOrigNDMRGstart(irrep_beta);
                  const int relBETA = beta_index - iHandler->getOrigNOCCstart(irrep_beta);
                  const double TwoDMval1 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindexALPHA + nOrbDMRG * ( DMRGindex1 + nOrbDMRG * DMRGindexBETA ) ) ];
                  const double TwoDMval2 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindexALPHA + nOrbDMRG * ( DMRGindexBETA + nOrbDMRG * DMRGindex1 ) ) ];
                  const double TwoDMval3 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindex1 + nOrbDMRG * ( DMRGindexBETA + nOrbDMRG * DMRGindexALPHA ) ) ];
                  value += 2 * (   TwoDMval1 * VmatRotated->get(irrep2, irrep_alpha, irrep4, irrep_beta, relIndex2, relALPHA, relIndex4, relBETA)
                     + (TwoDMval2+TwoDMval3) * VmatRotated->get(irrep2, irrep4, irrep_alpha, irrep_beta, relIndex2, relIndex4, relALPHA, relBETA) );
               }
            }
         }
         return value;
         
      }
   }
   
   return value;

}

void CheMPS2::CASSCF::buildFmat(){

   for (int cnt=0; cnt<numberOfIrreps; cnt++){
      const int nOrbitals = iHandler->getNORB(cnt);
      for (int cnt2=0; cnt2<nOrbitals*nOrbitals; cnt2++){
         int row = cnt2 % nOrbitals;
         int col = cnt2 / nOrbitals;
         Fmatrix[cnt][cnt2] = FmatHelper(iHandler->getOrigNOCCstart(cnt) + row, iHandler->getOrigNOCCstart(cnt) + col);
      }
   }

}

double CheMPS2::CASSCF::Fmat(const int index1, const int index2) const{

   const int irrep1 = HamOrig->getOrbitalIrrep(index1);
   const int irrep2 = HamOrig->getOrbitalIrrep(index2);
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   return Fmatrix[irrep1][ index1 - iHandler->getOrigNOCCstart(irrep1) + iHandler->getNORB(irrep1) * ( index2 - iHandler->getOrigNOCCstart(irrep1) ) ];

}

double CheMPS2::CASSCF::FmatHelper(const int index1, const int index2) const{
   
   const int irrep1 = HamOrig->getOrbitalIrrep(index1);
   const int irrep2 = HamOrig->getOrbitalIrrep(index2);
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   if (index1 >= iHandler->getOrigNVIRTstart(irrep1)){ return 0.0; } //index1 virtual: 1/2DM returns 0.0
   
   double value = 0.0;
   
   if (index1 < iHandler->getOrigNDMRGstart(irrep1)){ //index1 occupied
   
      value += 2 * ( TmatRotated(index2,index1) + QmatOCC(index2,index1) + QmatACT(index2,index1) );
   
   } else { //index1 active
   
      const int DMRGindex1 = iHandler->getDMRGcumulative(irrep1) + index1 - iHandler->getOrigNDMRGstart(irrep1);
   
      //1DM will only return non-zero if r_index also active, and corresponds to the same irrep
      for (int r_index=iHandler->getOrigNDMRGstart(irrep1); r_index<iHandler->getOrigNVIRTstart(irrep1); r_index++){
         const int DMRGindexR = iHandler->getDMRGcumulative(irrep1) + r_index - iHandler->getOrigNDMRGstart(irrep1);
         value += DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindexR ] * ( TmatRotated(index2,r_index) + QmatOCC(index2,r_index) );
      }
      
      //All the summation indices are active
      const int relIndex2 = index2 - iHandler->getOrigNOCCstart(irrep2);
      for (int irrep_r=0; irrep_r<numberOfIrreps; irrep_r++){
         const int productIrrep = SymmInfo.directProd(irrep1,irrep_r);
         for (int irrep_s=0; irrep_s<numberOfIrreps; irrep_s++){
            int irrep_t = SymmInfo.directProd(productIrrep,irrep_s);
            for (int r_index=iHandler->getOrigNDMRGstart(irrep_r); r_index<iHandler->getOrigNVIRTstart(irrep_r); r_index++){
               const int DMRGindexR = iHandler->getDMRGcumulative(irrep_r) + r_index - iHandler->getOrigNDMRGstart(irrep_r);
               const int relR = r_index - iHandler->getOrigNOCCstart(irrep_r);
               for (int s_index=iHandler->getOrigNDMRGstart(irrep_s); s_index<iHandler->getOrigNVIRTstart(irrep_s); s_index++){
                  const int DMRGindexS = iHandler->getDMRGcumulative(irrep_s) + s_index - iHandler->getOrigNDMRGstart(irrep_s);
                  const int relS = s_index - iHandler->getOrigNOCCstart(irrep_s);
                  for (int t_index=iHandler->getOrigNDMRGstart(irrep_t); t_index<iHandler->getOrigNVIRTstart(irrep_t); t_index++){
                     const int DMRGindexT = iHandler->getDMRGcumulative(irrep_t) + t_index - iHandler->getOrigNDMRGstart(irrep_t);
                     const int relT = t_index - iHandler->getOrigNOCCstart(irrep_t);
                     const double TwoDMvalue = DMRG2DM[ DMRGindex1 + nOrbDMRG * ( DMRGindexR + nOrbDMRG * ( DMRGindexS + nOrbDMRG * DMRGindexT ) ) ];
                     value += TwoDMvalue * VmatRotated->get(irrep2, irrep_r, irrep_s, irrep_t, relIndex2, relR, relS, relT);
                  }
               }
            }
         }
      }
   
   }
   
   return value;

}


