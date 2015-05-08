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

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <assert.h>

#include "DMRG.h"
#include "Lapack.h"
#include "TensorK.h"
#include "TensorM.h"
#include "TensorGYZ.h"
#include "Gsl.h"
#include "Heff.h"

using std::cout;
using std::endl;

void CheMPS2::DMRG::calc2DMandCorrelations(){

   //First get the whole MPS into left-canonical form
   int index = L-2;
   Sobject * denS = new Sobject(index,denBK->gIrrep(index),denBK->gIrrep(index+1),denBK);
   denS->Join(MPS[index],MPS[index+1]);
   Heff Solver(denBK, Prob);
   double Energy = 0.0;
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
      for (int cnt=0; cnt<nStates-1; cnt++){ delete [] VeffTilde[cnt]; }
      delete [] VeffTilde;
   }
   Energy += Prob->gEconst();
   if (Energy<TotalMinEnergy){ TotalMinEnergy = Energy; }
   denS->Split(MPS[index],MPS[index+1],OptScheme->getD(OptScheme->getNInstructions()-1),true,true);
   delete denS;
   
   cout << "**************************************" << endl;
   cout << "** 2DM and Correlations calculation **" << endl;
   cout << "**************************************" << endl;
   updateMovingRightSafe(index);
   
   TensorDiag * Norm = new TensorDiag(L, denBK);
   MPS[L-1]->QR(Norm);
   delete Norm;
   
   //Allocate space for the 2DM
   if (the2DMallocated){
      delete the2DM;
      the2DMallocated = false;
   }
   the2DM = new TwoDM(denBK, Prob);
   the2DMallocated = true;
   
   //Then calculate step by step the 2DM
   for (int siteindex=L-1; siteindex>=0; siteindex--){
      the2DM->FillSite(MPS[siteindex], Ltensors, F0tensors, F1tensors, S0tensors, S1tensors);
      if (siteindex>0){
         TensorDiag * Left = new TensorDiag(siteindex, denBK);
         MPS[siteindex]->LQ(Left);
         MPS[siteindex-1]->RightMultiply(Left);
         delete Left;
         updateMovingLeftSafe2DM(siteindex-1);
      }
   }
   
   //Then perform two checks: double trace & energy
   double NtimesNminus1 = the2DM->doubletrace2DMA();
   cout << "   N(N-1) = " << denBK->gN() * (denBK->gN() - 1) << " and calculated by double trace of the 2DM-A = " << NtimesNminus1 << endl;
   
   double Energy2DMA = the2DM->calcEnergy();
   cout << "   Energy obtained by Heffective at edge = " << Energy << " and as Econst + 0.5*trace(2DM-A*Ham) = " << Energy2DMA << endl;

   //Now the MPS has the gauge form CRRRRRRRRR
   //Allocate space for the Correlations
   if (theCorrAllocated){
      delete theCorr;
      theCorrAllocated = false;
   }
   theCorr = new Correlations(denBK, Prob, the2DM);
   theCorrAllocated = true;
   
   //Then calculate step by step the mutual information.
   //Define the following tensor arrays thereto. The native DMRG ones are at the edge and are hence small.
   TensorGYZ ** Gtensors = new TensorGYZ * [L-1];
   TensorGYZ ** Ytensors = new TensorGYZ * [L-1];
   TensorGYZ ** Ztensors = new TensorGYZ * [L-1];
   TensorK   ** Ktensors = new TensorK   * [L-1];
   TensorM   ** Mtensors = new TensorM   * [L-1];
   
   //Do the actual work
   for (int siteindex=1; siteindex<L; siteindex++){
   
      //Switch MPS gauge
      TensorDiag * Right = new TensorDiag(siteindex, denBK);
      MPS[siteindex-1]->QR(Right);
      MPS[siteindex]->LeftMultiply(Right);
      delete Right;
      
      //Update the tensors
      const int dimL = denBK->gMaxDimAtBound(siteindex-1);
      const int dimR = denBK->gMaxDimAtBound(siteindex);
      double * workmemLR = new double[dimL*dimR];
      for (int previousindex=0; previousindex<siteindex-1; previousindex++){
         TensorGYZ * newG = new TensorGYZ(siteindex, 'G', denBK);
         TensorGYZ * newY = new TensorGYZ(siteindex, 'Y', denBK);
         TensorGYZ * newZ = new TensorGYZ(siteindex, 'Z', denBK);
         TensorK   * newK = new TensorK(siteindex, denBK->gIrrep(previousindex), denBK);
         TensorM   * newM = new TensorM(siteindex, denBK->gIrrep(previousindex), denBK);

         newG->update(MPS[siteindex-1], Gtensors[previousindex], workmemLR);
         newY->update(MPS[siteindex-1], Ytensors[previousindex], workmemLR);
         newZ->update(MPS[siteindex-1], Ztensors[previousindex], workmemLR);
         newK->update(Ktensors[previousindex], MPS[siteindex-1], workmemLR, false);
         newM->update(Mtensors[previousindex], MPS[siteindex-1], workmemLR, false);
         
         delete Gtensors[previousindex];
         delete Ytensors[previousindex];
         delete Ztensors[previousindex];
         delete Ktensors[previousindex];
         delete Mtensors[previousindex];
         
         Gtensors[previousindex] = newG;
         Ytensors[previousindex] = newY;
         Ztensors[previousindex] = newZ;
         Ktensors[previousindex] = newK;
         Mtensors[previousindex] = newM;
      }
      delete [] workmemLR;
      
      //Construct the new tensors
      Gtensors[siteindex-1] = new TensorGYZ(siteindex, 'G', denBK);
      Ytensors[siteindex-1] = new TensorGYZ(siteindex, 'Y', denBK);
      Ztensors[siteindex-1] = new TensorGYZ(siteindex, 'Z', denBK);
      Ktensors[siteindex-1] = new TensorK(siteindex, denBK->gIrrep(siteindex-1), denBK);
      Mtensors[siteindex-1] = new TensorM(siteindex, denBK->gIrrep(siteindex-1), denBK);
      
      Gtensors[siteindex-1]->construct(MPS[siteindex-1]);
      Ytensors[siteindex-1]->construct(MPS[siteindex-1]);
      Ztensors[siteindex-1]->construct(MPS[siteindex-1]);
      Ktensors[siteindex-1]->construct(MPS[siteindex-1]);
      Mtensors[siteindex-1]->construct(MPS[siteindex-1]);
      
      //Use the tensors to fill in the Correlations
      theCorr->FillSite(MPS[siteindex], Gtensors, Ytensors, Ztensors, Ktensors, Mtensors);
      
   }
   
   //Clean-up
   for (int previousindex=0; previousindex<L-1; previousindex++){
      delete Gtensors[previousindex];
      delete Ytensors[previousindex];
      delete Ztensors[previousindex];
      delete Ktensors[previousindex];
      delete Mtensors[previousindex];
   }
   
   delete [] Gtensors;
   delete [] Ytensors;
   delete [] Ztensors;
   delete [] Ktensors;
   delete [] Mtensors;
   
   cout << "   Single-orbital entropies (Hamiltonian index order is used!) = [ ";
   for (int index=0; index < L-1; index++){ cout << theCorr->SingleOrbitalEntropy_HAM(index) << " , "; }
   cout << theCorr->SingleOrbitalEntropy_HAM(L-1) << " ]." << endl;
   
   for (int power=0; power<=2; power++){
      cout << "   Idistance(" << power << ") = " << theCorr->MutualInformationDistance((double)power) << endl;
   }
   cout << "**************************************" << endl;

}

CheMPS2::TwoDM * CheMPS2::DMRG::get2DM(){ return the2DM; }

CheMPS2::Correlations * CheMPS2::DMRG::getCorrelations(){ return theCorr; }

double CheMPS2::DMRG::getSpecificCoefficient(int * coeff) const{
   
   int * alpha = new int[ L ];
   int * beta  = new int[ L ];
   for ( int orb=0; orb<L; orb++ ){
      if ( coeff[orb] == 0 ){ alpha[orb] = 0; beta[orb] = 0; }
      if ( coeff[orb] == 1 ){ alpha[orb] = 1; beta[orb] = 0; }
      if ( coeff[orb] == 2 ){ alpha[orb] = 1; beta[orb] = 1; }
   }
   const double FCIcoeff = getFCIcoefficient( alpha, beta );
   delete [] alpha;
   delete [] beta;
   return FCIcoeff;

}

double CheMPS2::DMRG::getFCIcoefficient(int * alpha, int * beta) const{ //DMRGcoeff = coeff[Hamindex = Prob->gf2(DMRGindex)]

   //Check if it's possible
   {
      int nTot = 0;
      int twoSz = 0;
      int iTot = 0;
      for (int DMRGindex=0; DMRGindex<L; DMRGindex++){
         int HamIndex = (Prob->gReorderD2h()) ? Prob->gf2(DMRGindex) : DMRGindex;
         assert( ( alpha[HamIndex] == 0 ) || ( alpha[HamIndex] == 1 ) );
         assert( (  beta[HamIndex] == 0 ) || (  beta[HamIndex] == 1 ) );
         nTot  += alpha[HamIndex] + beta[HamIndex];
         twoSz += alpha[HamIndex] - beta[HamIndex];
         if ((alpha[HamIndex]+beta[HamIndex])==1){ iTot = Irreps::directProd(iTot,denBK->gIrrep(DMRGindex)); }
      }
      if ( Prob->gN() != nTot ){
         cout << "Ndesired = " << Prob->gN() << " and Ntotal in alpha and beta strings = " << nTot << endl;
         return 0.0;
      }
      // 2Sz can be -Prob->2S() ; -Prob->2S()+2 ; -Prob->2S()+4 ; ... ; Prob->2S()
      if ( ( Prob->gTwoS() < twoSz ) || ( twoSz < -Prob->gTwoS() ) || ( ( Prob->gTwoS() - twoSz ) % 2 != 0 ) ){
         cout << "2Sdesired = " << Prob->gTwoS() << " and Sz in alpha and beta strings = " << twoSz << endl;
         return 0.0;
      }
      if ( Prob->gIrrep() != iTot ){
         cout << "Idesired = " << Prob->gIrrep() << " and Irrep of alpha and beta strings = " << iTot << endl;
         return 0.0;
      }
   }
   
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
      const int HamIndex = (Prob->gReorderD2h()) ? Prob->gf2(DMRGindex) : DMRGindex;
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
               }
            }
            if ( indexSL != -1 ){
               int dimL = jumpL[ indexSL+1 ] - jumpL[ indexSL ];
               double * Tblock = MPS[DMRGindex]->gStorage(NL,TwoSLvalue,IL,NR,TwoSRvalue,IR);
               double prefactor = sqrt( TwoSRvalue + 1 )
                                * gsl_sf_coupling_3j(TwoSLvalue, spread, TwoSRvalue, twoSLz, twoSzloc, -twoSRz)
                                * Heff::phase( -TwoSLvalue + spread - twoSRz );
               double add2array = 1.0;
               char notrans = 'N';
               dgemm_( &notrans, &notrans, &dimFirst, &dimR, &dimL, &prefactor, arrayL + jumpL[indexSL], &dimFirst, Tblock, &dimL,
                                                                    &add2array, arrayR + jumpR[cntSR], &dimFirst);
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
   
   const double theCoeff = arrayL[0];
   
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
   
   return theCoeff;

}

void CheMPS2::DMRG::calcVeffTilde(double * result, Sobject * currentS, int state_number){

   int dimTot = currentS->gKappa2index(currentS->gNKappa());
   for (int cnt=0; cnt<dimTot; cnt++){ result[cnt] = 0.0; }
   int index = currentS->gIndex();
   
   const int dimL = std::max(denBK->gMaxDimAtBound(index),   Exc_BKs[state_number]->gMaxDimAtBound(index)   );
   const int dimR = std::max(denBK->gMaxDimAtBound(index+2), Exc_BKs[state_number]->gMaxDimAtBound(index+2) );
   double * workmem = new double[dimL * dimR];
   
   //Construct Sup
   Sobject * Sup = new Sobject(index,Exc_BKs[state_number]->gIrrep(index),Exc_BKs[state_number]->gIrrep(index+1),Exc_BKs[state_number]);
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

void CheMPS2::DMRG::calcOverlapsWithLowerStates(){

   for (int state=0; state<nStates-1; state++){
 
      //Don't do updatemovingRightSafe here, as with storeRenormOp some left boundary operators are stored and removed from mem.
      const int cnt = L-2;
      if (isAllocated[cnt]==2){
         deleteTensors(cnt, false);
         isAllocated[cnt]=0;
      }
      if (isAllocated[cnt]==0){
         allocateTensors(cnt, true);
         isAllocated[cnt]=1;
      }
      updateMovingRight(cnt);
 
      TensorO * Otemp = new TensorO(L, true, Exc_BKs[state], denBK, Prob);
      Otemp->update(Exc_MPSs[state][L-1], MPS[L-1], Exc_Overlaps[state][L-2]);
      double overlap = Otemp->gStorage()[0];
      delete Otemp;
      
      cout << "The overlap between the current state and state " << state << " is : " << overlap << endl;

      if (isAllocated[cnt]==1){
         deleteTensors(cnt, true);
         isAllocated[cnt]=0;
      }
   
   }

}

void CheMPS2::DMRG::calcOverlapsWithLowerStatesDuringSweeps_debug(double ** VeffTilde, Sobject * denS){

   denS->prog2symm();
   int dim = denS->gKappa2index(denS->gNKappa());

   for (int state=0; state<nStates-1; state++){
      
      int inc = 1;
      double overlap = ddot_(&dim,denS->gStorage(),&inc,VeffTilde[state],&inc) / sqrt(Exc_Eshifts[state]);
      cout << "   Debugger: The overlap (based on VeffTilde) between the current state and state " << state << " is : " << overlap << endl;
   
   }
   
   denS->symm2prog();

}

