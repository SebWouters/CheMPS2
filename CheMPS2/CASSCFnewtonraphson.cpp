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

   double gradNorm = 1.0;
   double Energy;
   
   Hamiltonian * HamDMRG = new Hamiltonian(nOrbDMRG, SymmInfo.getGroupNumber(), iHandler->getIrrepOfEachDMRGorbital());
   int N = Nelectrons;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){ N -= 2*iHandler->getNOCC(irrep); }
   Problem * Prob = new Problem(HamDMRG, TwoS, N, Irrep);
   Prob->SetupReorderD2h(); //Doesn't matter if the group isn't d2h, Prob checks it.
   
   int maxlinsize = 0;
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ if (iHandler->getNORB(cnt) > maxlinsize){ maxlinsize = iHandler->getNORB(cnt); } }
   const bool doBlockWise = (maxlinsize <= CheMPS2::CASSCF_maxlinsizeCutoff) ? false : true; //Only if bigger, do we want to work blockwise
   int maxBlockSize = maxlinsize;
   if (doBlockWise){
      int factor   = (int) (ceil( (1.0 * maxlinsize) / CheMPS2::CASSCF_maxlinsizeCutoff ) + 0.01);
      maxBlockSize = (int) (ceil( (1.0 * maxlinsize) / factor ) + 0.01);
      cout << "CASSCF info: the max. # orb per irrep = " << maxlinsize << " and was truncated to " << maxBlockSize << " for the 2-body rotation." << endl;
   }
   
   //One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
   int maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!!!
   int sizeWorkmem1 = max( maxBSpower4 , maxlinsize*maxlinsize*3 ); //Second argument for updateUnitary
       sizeWorkmem1 = max( sizeWorkmem1 , nOrbDMRG*nOrbDMRG ); //Second argument for eigenvectors 1DM (calcNOON)
   int sizeWorkmem2 = max( maxBSpower4 , maxlinsize*maxlinsize*2 );
   if (CheMPS2::CASSCF_rotate2DMtoNO){
      sizeWorkmem2 = max( sizeWorkmem2 , nOrbDMRG*nOrbDMRG*nOrbDMRG*nOrbDMRG ); //Second argument to rotate 2DM
   }
   double * mem1 = new double[sizeWorkmem1];
   double * mem2 = new double[sizeWorkmem2];
   double * mem3 = NULL;
   if (doBlockWise){ mem3 = new double[maxBSpower4]; }
   
   if (CheMPS2::CASSCF_storeUnitary){
   
      struct stat stFileInfo;
      int intStat = stat(CheMPS2::CASSCF_unitaryStorageName.c_str(),&stFileInfo);
      if (intStat==0){ unitary->loadU(); }
   
   }

   while (gradNorm > CheMPS2::CASSCF_gradientNormThreshold){
   
      //Update the unitary transformations based on the previous unitary transformation and the xmatrix
      unitary->updateUnitary(mem1, mem2);
      
      if ((CheMPS2::CASSCF_storeUnitary) && (gradNorm!=1.0)){ unitary->saveU(); }
   
      //Setup rotated Hamiltonian matrix elements based on unitary transformations
      if (doBlockWise){ fillRotatedHamInMemoryBlockWise(mem1, mem2, mem3, maxBlockSize); }
      else{             fillRotatedHamAllInMemory(mem1, mem2); }
   
      //Fill HamDMRG based on the HamRotated --> requires QmatrixOcc to be filled
      buildQmatrixOCC();
      fillHamDMRG(HamDMRG);
      
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
      
      calcNOON(mem1);
      if (CheMPS2::CASSCF_rotate2DMtoNO){
         rotate2DMand1DM(N,mem1, mem2);
         unitary->rotateUnitaryNOeigenvecs(mem1, mem2);
         //Unitary update implies matrix element update. --> Can actually be done for DMRG orbitals alone = faster --> TODO
         if (doBlockWise){ fillRotatedHamInMemoryBlockWise(mem1, mem2, mem3, maxBlockSize); }
         else{             fillRotatedHamAllInMemory(mem1, mem2); }
      }
      
      //In order to construct the gradient and Hessian faster, the following three matrices need to be updated (in this order!)
      if (CheMPS2::CASSCF_rotate2DMtoNO){ buildQmatrixOCC(); }
      buildQmatrixACT();
      buildFmat();

      gradNorm = updateXmatrixAugmentedHessianNR(); //updateXmatrixNewtonRaphson();
      
      //PrintCoeff_C2(theDMRG); //Print coeff for C2 to discern (for 1Ag) ^1Sigma_g^+ <--> ^1Delta_g   &   (for 3B1u) ^3Sigma_u^+ <--> ^3Delta_u
      
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
      cout << "CASSCF :: Eigenvalues of the Hessian = [ ";
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
   cout << "CASSCF :: Norm of the gradient = " << gradNorm << endl;
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

   int irrep1 = 0;
   while (index1 >= iHandler->getOrigNOCCstart(irrep1+1)){ irrep1++; }
   int irrep2 = 0;
   while (index2 >= iHandler->getOrigNOCCstart(irrep2+1)){ irrep2++; }
   
   if (irrep1 != irrep2){ return 0.0; } //From now on: irrep1 == irrep2
   
   int irrep3 = 0;
   while (index3 >= iHandler->getOrigNOCCstart(irrep3+1)){ irrep3++; }
   int irrep4 = 0;
   while (index4 >= iHandler->getOrigNOCCstart(irrep4+1)){ irrep4++; }
   
   if (irrep3 != irrep4){ return 0.0; } //From now on: irrep3 == irrep4
   
   double value = 0.0;
   
   if ( (irrep1 == irrep3) && (index2 == index3) ) { //irrep1 == irrep2 == irrep3 == irrep4
      value = Fmat(index1,index4) + Fmat(index4,index1);
   }
   
   if (index1 >= iHandler->getOrigNVIRTstart(irrep1)){ return value; }
   if (index3 >= iHandler->getOrigNVIRTstart(irrep3)){ return value; } //index1 and index3 are now certainly not virtual!
   
   if (index1 < iHandler->getOrigNDMRGstart(irrep1)){
      if (index3 < iHandler->getOrigNDMRGstart(irrep3)){
      
         // (index1,index3) (occupied,occupied) --> (alpha,beta) can be (occupied,occupied) or (active,active)
         if (index1==index3){
         
            // Part of one-body matrix elements if 1-RDM occ indices
            value += 4 * ( HamRotated->getTmat(index2,index4) + QmatOCC(index2,index4) + QmatACT(index2, index4) );
            value += 12 * HamRotated->getVmat(index2,index4,index1,index1) - 4 * HamRotated->getVmat(index2,index1,index4,index1);
            return value;
            
         } else { //index1 != index3
         
            value += 16 * HamRotated->getVmat(index2,index4,index1,index3)
                   - 4  * HamRotated->getVmat(index2,index1,index4,index3)
                   - 4  * HamRotated->getVmat(index2,index4,index3,index1);
            return value;
            
         }
      
      } else {
      
         const int DMRGindex3 = iHandler->getDMRGcumulative(irrep3) + index3 - iHandler->getOrigNDMRGstart(irrep3);
      
         // (index1,index3) (occupied,active) --> (alpha,beta) can be (active,occupied) or (occupied,active)
         for (int alpha_index=iHandler->getOrigNDMRGstart(irrep3); alpha_index<iHandler->getOrigNVIRTstart(irrep3); alpha_index++){
            const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep3) + alpha_index - iHandler->getOrigNDMRGstart(irrep3);
            value += 2 * DMRG1DM[ DMRGindex3 + nOrbDMRG * DMRGindexALPHA ] * ( 4 * HamRotated->getVmat(index2,index4,index1,alpha_index)
                                                                                 - HamRotated->getVmat(index2,index4,alpha_index,index1)
                                                                                 - HamRotated->getVmat(index2,index1,index4,alpha_index) );
         }
         return value;
         
      }
   } else {
      
      const int DMRGindex1 = iHandler->getDMRGcumulative(irrep1) + index1 - iHandler->getOrigNDMRGstart(irrep1);
   
      if (index3 < iHandler->getOrigNDMRGstart(irrep3)){
      
         // (index1,index3) (active,occupied) --> (alpha,beta) can be (active,occupied) or (occupied,active)
         for (int alpha_index=iHandler->getOrigNDMRGstart(irrep1); alpha_index<iHandler->getOrigNVIRTstart(irrep1); alpha_index++){
            const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep1) + alpha_index - iHandler->getOrigNDMRGstart(irrep1);
            value += 2 * DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindexALPHA ] * ( 4 * HamRotated->getVmat(index2,index4,alpha_index,index3)
                                                                                 - HamRotated->getVmat(index2,alpha_index,index4,index3)
                                                                                 - HamRotated->getVmat(index2,index4,index3,alpha_index) );
         }
         return value;
         
      } else {
      
         const int DMRGindex3 = iHandler->getDMRGcumulative(irrep3) + index3 - iHandler->getOrigNDMRGstart(irrep3);
      
         // (index1,index3) (active,active) --> (alpha,beta) can be (occupied,occupied) or (active,active)
         // Case1: (alpha,beta)==(occ,occ) --> alpha == beta
         // Part of one-body matrix elements if 1-RDM active indices
         value += 2 * DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindex3 ] * ( HamRotated->getTmat(index2,index4) + QmatOCC(index2,index4) );
         
         // Case2: (alpha,beta)==(act,act)
         const int productIrrep = SymmInfo.directProd(irrep1,irrep3);
         for (int irrep_alpha=0; irrep_alpha<numberOfIrreps; irrep_alpha++){
            int irrep_beta = SymmInfo.directProd(productIrrep,irrep_alpha);
            for (int alpha_index=iHandler->getOrigNDMRGstart(irrep_alpha); alpha_index<iHandler->getOrigNVIRTstart(irrep_alpha); alpha_index++){
               const int DMRGindexALPHA = iHandler->getDMRGcumulative(irrep_alpha) + alpha_index - iHandler->getOrigNDMRGstart(irrep_alpha);
               for (int beta_index=iHandler->getOrigNDMRGstart(irrep_beta); beta_index<iHandler->getOrigNVIRTstart(irrep_beta); beta_index++){
                  const int DMRGindexBETA = iHandler->getDMRGcumulative(irrep_beta) + beta_index - iHandler->getOrigNDMRGstart(irrep_beta);
                  const double TwoDMval1 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindexALPHA + nOrbDMRG * ( DMRGindex1 + nOrbDMRG * DMRGindexBETA ) ) ];
                  const double TwoDMval2 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindexALPHA + nOrbDMRG * ( DMRGindexBETA + nOrbDMRG * DMRGindex1 ) ) ];
                  const double TwoDMval3 = DMRG2DM[ DMRGindex3 + nOrbDMRG * ( DMRGindex1 + nOrbDMRG * ( DMRGindexBETA + nOrbDMRG * DMRGindexALPHA ) ) ];
                  value += 2 * (   TwoDMval1 * HamRotated->getVmat(index2,alpha_index,index4,beta_index)
                               + ( TwoDMval2 + TwoDMval3 ) * HamRotated->getVmat(index2,index4,alpha_index,beta_index) );
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

   int irrep1 = 0;
   while (index1 >= iHandler->getOrigNOCCstart(irrep1+1)){ irrep1++; }
   int irrep2 = 0;
   while (index2 >= iHandler->getOrigNOCCstart(irrep2+1)){ irrep2++; }
   
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   return Fmatrix[irrep1][ index1 - iHandler->getOrigNOCCstart(irrep1) + iHandler->getNORB(irrep1) * ( index2 - iHandler->getOrigNOCCstart(irrep1) ) ];

}

double CheMPS2::CASSCF::FmatHelper(const int index1, const int index2) const{
   
   int irrep1 = 0;
   while (index1 >= iHandler->getOrigNOCCstart(irrep1+1)){ irrep1++; }
   int irrep2 = 0;
   while (index2 >= iHandler->getOrigNOCCstart(irrep2+1)){ irrep2++; }
   
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   if (index1 >= iHandler->getOrigNVIRTstart(irrep1)){ return 0.0; } //index1 virtual: 1/2DM returns 0.0
   
   double value = 0.0;
   
   if (index1 < iHandler->getOrigNDMRGstart(irrep1)){ //index1 occupied
   
      value += 2 * ( HamRotated->getTmat(index2,index1) + QmatOCC(index2,index1) + QmatACT(index2,index1) );
   
   } else { //index1 active
   
      const int DMRGindex1 = iHandler->getDMRGcumulative(irrep1) + index1 - iHandler->getOrigNDMRGstart(irrep1);
   
      //1DM will only return non-zero if r_index also active, and corresponds to the same irrep
      for (int r_index=iHandler->getOrigNDMRGstart(irrep1); r_index<iHandler->getOrigNVIRTstart(irrep1); r_index++){
         const int DMRGindexR = iHandler->getDMRGcumulative(irrep1) + r_index - iHandler->getOrigNDMRGstart(irrep1);
         value += DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindexR ] * ( HamRotated->getTmat(index2,r_index) + QmatOCC(index2,r_index) );
      }
      
      //All the summation indices are active
      for (int irrep_r=0; irrep_r<numberOfIrreps; irrep_r++){
         const int productIrrep = SymmInfo.directProd(irrep1,irrep_r);
         for (int irrep_s=0; irrep_s<numberOfIrreps; irrep_s++){
            int irrep_t = SymmInfo.directProd(productIrrep,irrep_s);
            for (int r_index=iHandler->getOrigNDMRGstart(irrep_r); r_index<iHandler->getOrigNVIRTstart(irrep_r); r_index++){
               const int DMRGindexR = iHandler->getDMRGcumulative(irrep_r) + r_index - iHandler->getOrigNDMRGstart(irrep_r);
               for (int s_index=iHandler->getOrigNDMRGstart(irrep_s); s_index<iHandler->getOrigNVIRTstart(irrep_s); s_index++){
                  const int DMRGindexS = iHandler->getDMRGcumulative(irrep_s) + s_index - iHandler->getOrigNDMRGstart(irrep_s);
                  for (int t_index=iHandler->getOrigNDMRGstart(irrep_t); t_index<iHandler->getOrigNVIRTstart(irrep_t); t_index++){
                     const int DMRGindexT = iHandler->getDMRGcumulative(irrep_t) + t_index - iHandler->getOrigNDMRGstart(irrep_t);
                     const double TwoDMvalue = DMRG2DM[ DMRGindex1 + nOrbDMRG * ( DMRGindexR + nOrbDMRG * ( DMRGindexS + nOrbDMRG * DMRGindexT ) ) ];
                     value += TwoDMvalue * HamRotated->getVmat(index2, r_index, s_index, t_index);
                  }
               }
            }
         }
      }
   
   }
   
   return value;

}


