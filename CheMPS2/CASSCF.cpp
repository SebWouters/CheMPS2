/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013 Sebastian Wouters

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
#include <sstream>

#include "CASSCF.h"
#include "Lapack.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

CheMPS2::CASSCF::CASSCF(const string filename){

   HamOrig = new Hamiltonian(filename);
   
   int * orb2irrepHamOrig = new int[HamOrig->getL()];
   for (int cnt=0; cnt<HamOrig->getL(); cnt++){ orb2irrepHamOrig[cnt] = HamOrig->getOrbitalIrrep(cnt); }
   HamRotated = new Hamiltonian(HamOrig->getL(), HamOrig->getNGroup(), orb2irrepHamOrig);
   delete [] orb2irrepHamOrig;
   
   L = HamOrig->getL();
   
   SymmInfo.setGroup(HamOrig->getNGroup());
   
   numberOfIrreps = SymmInfo.getNumberOfIrreps();
   
   OrbPerIrrep = new int[numberOfIrreps];
   
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ OrbPerIrrep[cnt] = 0; }
   for (int cnt=0; cnt<L; cnt++){ OrbPerIrrep[HamOrig->getOrbitalIrrep(cnt)]++; }

   allocateAndFillOCC(filename);
   
   checkHF();
   
   setupStartCalled = false;
   
}

CheMPS2::CASSCF::CASSCF(Hamiltonian * HamIn, int * DOCCin, int * SOCCin){

   L = HamIn->getL();
   int SyGroup = HamIn->getNGroup();
   int * OrbIrreps = new int[L];
   for (int cnt=0; cnt<L; cnt++){ OrbIrreps[cnt] = HamIn->getOrbitalIrrep(cnt); }

   HamOrig    = new Hamiltonian(L, SyGroup, OrbIrreps);
   HamRotated = new Hamiltonian(L, SyGroup, OrbIrreps);
   
   delete [] OrbIrreps;
   
   SymmInfo.setGroup(SyGroup);
   
   HamOrig->setEconst(HamIn->getEconst());
   for (int cnt=0; cnt<L; cnt++){
      int I1 = HamIn->getOrbitalIrrep(cnt);
      for (int cnt2=cnt; cnt2<L; cnt2++){
         int I2 = HamIn->getOrbitalIrrep(cnt2);
         if (I1==I2){
            HamOrig->setTmat(cnt,cnt2,HamIn->getTmat(cnt,cnt2));
         }
         for (int cnt3=cnt; cnt3<L; cnt3++){
            int I3 = HamIn->getOrbitalIrrep(cnt3);
            for (int cnt4=cnt2; cnt4<L; cnt4++){
               int I4 = HamIn->getOrbitalIrrep(cnt4);
               if (SymmInfo.directProd(I1,I2) == SymmInfo.directProd(I3,I4)){
                  HamOrig->setVmat(cnt,cnt2,cnt3,cnt4,HamIn->getVmat(cnt,cnt2,cnt3,cnt4));
               }
            }
         }
      }
   }
   
   numberOfIrreps = SymmInfo.getNumberOfIrreps();
   
   OrbPerIrrep = new int[numberOfIrreps];
   
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ OrbPerIrrep[cnt] = 0; }
   for (int cnt=0; cnt<L; cnt++){ OrbPerIrrep[HamOrig->getOrbitalIrrep(cnt)]++; }

   allocateAndFillOCC(DOCCin, SOCCin);
   
   checkHF();
   
   setupStartCalled = false;
   
}

CheMPS2::CASSCF::~CASSCF(){

   delete HamOrig;
   delete HamRotated;
   
   delete [] DOCC;
   delete [] SOCC;
   delete [] OrbPerIrrep;
   
   if (setupStartCalled){
   
      delete [] Nocc;
      delete [] NDMRG;
      delete [] Nvirt;
      
      delete [] irrepOfEachDMRGOrbital;
      delete [] listDMRG;
      delete [] listCondensed;
      
      for (int irrep=0; irrep<numberOfIrreps; irrep++){
         for (int case_=0; case_<3; case_++){ delete [] xmatrix[irrep][case_]; }
         delete [] xmatrix[irrep];
      }
      delete [] xmatrix;
      
      for (int irrep=0; irrep<numberOfIrreps; irrep++){ delete [] unitary[irrep]; }
      delete [] unitary;
      
      delete [] jumpsHamOrig;
      delete [] DMRG1DM;
      delete [] DMRG2DM;
      delete [] x_firstindex;
      delete [] x_secondindex;
      
      for (int cnt=0; cnt<numberOfIrreps; cnt++){ delete [] Fmatrix[cnt]; }
      delete [] Fmatrix;
      
   }

}

int CheMPS2::CASSCF::getNumberOfIrreps(){ return numberOfIrreps; }

void CheMPS2::CASSCF::copyXsolutionBack(double * vector){

   for (int irrep=0; irrep<numberOfIrreps; irrep++){
         
      for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
         for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){
            int index1 = jumpsHamOrig[irrep] + Nocc[irrep] + cntDMRG;
            int index2 = jumpsHamOrig[irrep] + cntOcc;
            int xsolindex = x_tolin(index1,index2);
            if (xsolindex==-1){ cout << "Error: xsolindex==-1" << endl; }
            xmatrix[irrep][0][ cntDMRG + NDMRG[irrep] * cntOcc ] = vector[xsolindex];
         }
      }
      for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){
         for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){
            int index1 = jumpsHamOrig[irrep] + Nocc[irrep] + NDMRG[irrep] + cntVirt;
            int index2 = jumpsHamOrig[irrep] + Nocc[irrep] + cntDMRG;
            int xsolindex = x_tolin(index1,index2);
            if (xsolindex==-1){ cout << "Error: xsolindex==-1" << endl; }
            xmatrix[irrep][1][ cntVirt + Nvirt[irrep] * cntDMRG ] = vector[xsolindex];
         }
      }
      for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
         for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){
            int index1 = jumpsHamOrig[irrep] + Nocc[irrep] + NDMRG[irrep] + cntVirt;
            int index2 = jumpsHamOrig[irrep] + cntOcc;
            int xsolindex = x_tolin(index1,index2);
            if (xsolindex==-1){ cout << "Error: xsolindex==-1" << endl; }
            xmatrix[irrep][2][ cntVirt + Nvirt[irrep] * cntOcc ] = vector[xsolindex];
         }
      }
   }

}

int CheMPS2::CASSCF::x_tolin(const int p_index, const int q_index){

   for (int cnt=0; cnt<x_linearlength; cnt++){
      if ((p_index == x_firstindex[cnt]) && (q_index == x_secondindex[cnt])){ return cnt; }
   }
   
   return -1;

}

double CheMPS2::CASSCF::get2DMrotated(const int index1, const int index2, const int index3, const int index4) const{

   if ((index1<0) || (index1>=L)){ return 0.0; }
   if ((index2<0) || (index2>=L)){ return 0.0; }
   if ((index3<0) || (index3>=L)){ return 0.0; }
   if ((index4<0) || (index4>=L)){ return 0.0; } //Within bounds now

   const int irrep1 = HamOrig->getOrbitalIrrep(index1);
   const int irrep2 = HamOrig->getOrbitalIrrep(index2);
   const int irrep3 = HamOrig->getOrbitalIrrep(index3);
   const int irrep4 = HamOrig->getOrbitalIrrep(index4);
   
   if (SymmInfo.directProd(irrep1,irrep2) != SymmInfo.directProd(irrep3,irrep4)){ return 0.0; } //Matrix element is symmetry allowed
   
   int cnt1 = index1 - jumpsHamOrig[irrep1];
   int cnt2 = index2 - jumpsHamOrig[irrep2];
   int cnt3 = index3 - jumpsHamOrig[irrep3];
   int cnt4 = index4 - jumpsHamOrig[irrep4];
   
   if (cnt1 >= Nocc[irrep1] + NDMRG[irrep1]){ return 0.0; }
   if (cnt2 >= Nocc[irrep2] + NDMRG[irrep2]){ return 0.0; }
   if (cnt3 >= Nocc[irrep3] + NDMRG[irrep3]){ return 0.0; }
   if (cnt4 >= Nocc[irrep4] + NDMRG[irrep4]){ return 0.0; } //No virtual orbitals now
   
   const int numberOcc = ((cnt1 < Nocc[irrep1]) ? 1 : 0) + ((cnt2 < Nocc[irrep2]) ? 1 : 0) + ((cnt3 < Nocc[irrep3]) ? 1 : 0) + ((cnt4 < Nocc[irrep4]) ? 1 : 0);
   
   if ((numberOcc%2)!=0){ return 0.0; } //Now it is allowed --> need to find out its value.
   
   if (numberOcc==4){
      if (index3==index4){
         if ((index1==index2) && (index1==index3)){ return 2.0; } //all indices equal
         else { return 0.0; }
      }
      return (((index1 == index3) && (index2 == index4)) ? 4.0 : 0.0 ) - (((index1 == index4) && (index2 == index3)) ? 2.0 : 0.0 );
   }
   
   if (numberOcc==0){
   
      int DMRGindex1 = 0;
      int DMRGindex2 = 0;
      int DMRGindex3 = 0;
      int DMRGindex4 = 0;
      
      for (int bla=0; bla<irrep1; bla++){ DMRGindex1 += NDMRG[bla]; }
      for (int bla=0; bla<irrep2; bla++){ DMRGindex2 += NDMRG[bla]; }
      for (int bla=0; bla<irrep3; bla++){ DMRGindex3 += NDMRG[bla]; }
      for (int bla=0; bla<irrep4; bla++){ DMRGindex4 += NDMRG[bla]; }
      
      DMRGindex1 += cnt1 - Nocc[irrep1];
      DMRGindex2 += cnt2 - Nocc[irrep2];
      DMRGindex3 += cnt3 - Nocc[irrep3];
      DMRGindex4 += cnt4 - Nocc[irrep4];
      
      return DMRG2DM[DMRGindex1 + nOrbDMRG * ( DMRGindex2 + nOrbDMRG * (DMRGindex3 + nOrbDMRG * DMRGindex4 ) ) ];
      
   }
   
   //Now two indices are occ, two are active
   if (cnt1<Nocc[irrep1]){ // index1 in occ
   
      if (cnt3<Nocc[irrep3]){
         
         if (cnt1 != cnt3){
            return 0.0;
         } else {
            int DMRGindex2 = 0;
            int DMRGindex4 = 0;
            for (int bla=0; bla<irrep2; bla++){ DMRGindex2 += NDMRG[bla]; }
            for (int bla=0; bla<irrep4; bla++){ DMRGindex4 += NDMRG[bla]; }
            DMRGindex2 += cnt2 - Nocc[irrep2];
            DMRGindex4 += cnt4 - Nocc[irrep4];
            
            return 2 * DMRG1DM[DMRGindex2 + nOrbDMRG*DMRGindex4]; //2 * get1DMrotated(index2,index4);
         }
         
      }
      if (cnt4<Nocc[irrep4]){
      
         if (cnt1 != cnt4){
            return 0.0;
         } else {
            int DMRGindex2 = 0;
            int DMRGindex3 = 0;
            for (int bla=0; bla<irrep2; bla++){ DMRGindex2 += NDMRG[bla]; }
            for (int bla=0; bla<irrep3; bla++){ DMRGindex3 += NDMRG[bla]; }
            DMRGindex2 += cnt2 - Nocc[irrep2];
            DMRGindex3 += cnt3 - Nocc[irrep3];
            
            return - DMRG1DM[DMRGindex2 + nOrbDMRG*DMRGindex3]; // - get1DMrotated(index2,index3);
         }
      
      }
   
   } else { // index1 in act
   
      if (cnt3>=Nocc[irrep3]){
      
         if (cnt2 != cnt4){
            return 0.0;
         } else {
            int DMRGindex1 = 0;
            int DMRGindex3 = 0;
            for (int bla=0; bla<irrep1; bla++){ DMRGindex1 += NDMRG[bla]; }
            for (int bla=0; bla<irrep3; bla++){ DMRGindex3 += NDMRG[bla]; }
            DMRGindex1 += cnt1 - Nocc[irrep1];
            DMRGindex3 += cnt3 - Nocc[irrep3];
            
            return 2 * DMRG1DM[DMRGindex1 + nOrbDMRG*DMRGindex3]; //2 * get1DMrotated(index1,index3);
         }
      
      }
      if (cnt4>=Nocc[irrep4]){
      
         if (cnt2 != cnt3){
            return 0.0;
         } else {
            int DMRGindex1 = 0;
            int DMRGindex4 = 0;
            for (int bla=0; bla<irrep1; bla++){ DMRGindex1 += NDMRG[bla]; }
            for (int bla=0; bla<irrep4; bla++){ DMRGindex4 += NDMRG[bla]; }
            DMRGindex1 += cnt1 - Nocc[irrep1];
            DMRGindex4 += cnt4 - Nocc[irrep4];
            
            return - DMRG1DM[DMRGindex1 + nOrbDMRG*DMRGindex4]; //- get1DMrotated(index1,index4);
         }
      
      }
   
   }
   
   return 0.0; //Two occ, two active; but index1 and index2 are in the same space (occ/act) --> zero

}

void CheMPS2::CASSCF::copy2DMover(TwoDM * theDMRG2DM){

   for (int i1=0; i1<nOrbDMRG; i1++){
      for (int i2=0; i2<nOrbDMRG; i2++){
         for (int i3=0; i3<nOrbDMRG; i3++){
            for (int i4=0; i4<nOrbDMRG; i4++){
               DMRG2DM[i1 + nOrbDMRG * ( i2 + nOrbDMRG * (i3 + nOrbDMRG * i4 ) ) ] = theDMRG2DM->getTwoDMA_HAM(i1, i2, i3, i4);
            }
         }
      }
   }

}

double CheMPS2::CASSCF::get1DMrotated(const int index1, const int index2) const{

   if ((index1<0) || (index1>=L)){ return 0.0; }
   if ((index2<0) || (index2>=L)){ return 0.0; } //Within bounds now.

   const int irrep = HamOrig->getOrbitalIrrep(index1);

   if (irrep != HamOrig->getOrbitalIrrep(index2)){ return 0.0; } //Of same irrep now.
   
   int cnt1 = index1 - jumpsHamOrig[irrep];
   int cnt2 = index2 - jumpsHamOrig[irrep];
   
   if (cnt1 < Nocc[irrep]){ //If in occupied HF orbitals
      if (cnt1==cnt2){ return 2.0; }
      else{ return 0.0; }
   }
   
   cnt1 -= Nocc[irrep];
   cnt2 -= Nocc[irrep];
   
   if ((cnt1 >= 0) && (cnt1 < NDMRG[irrep]) && (cnt2 >= 0) && (cnt2 < NDMRG[irrep])){ //If both in DMRG orbitals
      int jump = 0;
      for (int bla=0; bla<irrep; bla++){ jump += NDMRG[bla]; }
      int DMRGindex1 = jump + cnt1;
      int DMRGindex2 = jump + cnt2;
      return DMRG1DM[DMRGindex1 + nOrbDMRG*DMRGindex2];
   }
   
   return 0.0; //All other cases : 0.0

}

void CheMPS2::CASSCF::setDMRG1DM(const int N){

   const double prefactor = 1.0/(N-1.0);

   for (int cnt1=0; cnt1<nOrbDMRG; cnt1++){
      for (int cnt2=cnt1; cnt2<nOrbDMRG; cnt2++){
         DMRG1DM[cnt1 + nOrbDMRG*cnt2] = 0.0;
         for (int cnt3=0; cnt3<nOrbDMRG; cnt3++){ DMRG1DM[cnt1 + nOrbDMRG*cnt2] += DMRG2DM[cnt1 + nOrbDMRG * (cnt3 + nOrbDMRG * (cnt2 + nOrbDMRG * cnt3 ) ) ]; }
         DMRG1DM[cnt1 + nOrbDMRG*cnt2] *= prefactor;
         DMRG1DM[cnt2 + nOrbDMRG*cnt1] = DMRG1DM[cnt1 + nOrbDMRG*cnt2];
      }
   }

}

void CheMPS2::CASSCF::calcNOON(){

   int size = nOrbDMRG * nOrbDMRG;
   double * copy = new double[size];
   double * eigenval = new double[nOrbDMRG];
   double * work = new double[size];
      
   for (int cnt=0; cnt<size; cnt++){ copy[cnt] = DMRG1DM[cnt]; }
      
   char jobz = 'N';
   char uplo = 'U';
   int info;
   dsyev_(&jobz,&uplo,&nOrbDMRG,copy,&nOrbDMRG,eigenval,work,&size,&info);
      
   cout << "CASSCF :: DMRG 1DM eigenvalues [NOON] = [ ";
   for (int cnt=0; cnt<nOrbDMRG-1; cnt++){ cout << eigenval[cnt] << " , "; }
   cout << eigenval[nOrbDMRG-1] << " ]." << endl;
      
   delete [] work;
   delete [] eigenval;
   delete [] copy;

}

void CheMPS2::CASSCF::calcNOON(double * eigenvecs){

   int size = nOrbDMRG * nOrbDMRG;
   double * eigenval = new double[nOrbDMRG];
   double * work = new double[size];

   for (int cnt=0; cnt<size; cnt++){ eigenvecs[cnt] = DMRG1DM[cnt]; }

   char jobz = 'V';
   char uplo = 'U';
   int info;
   int passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){

      if (NDMRG[irrep] > 0){

         //Calculate the eigenvectors and values per block
         dsyev_(&jobz, &uplo, NDMRG+irrep, eigenvecs + passed*(1+nOrbDMRG) ,&nOrbDMRG, eigenval + passed, work, &size, &info);

         //Print the NOON
         cout << "CASSCF :: DMRG 1DM eigenvalues [NOON] of irrep " << irrep << " = [ ";
         for (int cnt=0; cnt<NDMRG[irrep]-1; cnt++){ cout << eigenval[passed + NDMRG[irrep]-1-cnt] << " , "; }
         cout << eigenval[passed + 0] << " ]." << endl;

         //Sort the eigenvecs
         for (int col=0; col<NDMRG[irrep]/2; col++){
            for (int row=0; row<NDMRG[irrep]; row++){
               double temp = eigenvecs[passed + row + nOrbDMRG * (passed + NDMRG[irrep] - 1 - col)];
               eigenvecs[passed + row + nOrbDMRG * (passed + NDMRG[irrep] - 1 - col)] = eigenvecs[passed + row + nOrbDMRG * (passed + col)];
               eigenvecs[passed + row + nOrbDMRG * (passed + col)] = temp;
            }
         }

      }

      //Update the number of passed DMRG orbitals
      passed += NDMRG[irrep];

   }

   delete [] work;
   delete [] eigenval;

}

void CheMPS2::CASSCF::rotateUnitaryAnd2DMand1DM(const int N, double * eigenvecs, double * work){

   char notr = 'N';
   char tran = 'T';
   double alpha = 1.0;
   double beta = 0.0;

   int power1 = nOrbDMRG;
   int power2 = nOrbDMRG*nOrbDMRG;
   int power3 = nOrbDMRG*nOrbDMRG*nOrbDMRG;

   //2DM: Gamma_{ijkl} --> Gamma_{ajkl}
   dgemm_(&tran,&notr,&power1,&power3,&power1,&alpha,eigenvecs,&power1,DMRG2DM,&power1,&beta,work,&power1);
   //2DM: Gamma_{ajkl} --> Gamma_{ajkd}
   dgemm_(&notr,&notr,&power3,&power1,&power1,&alpha,work,&power3,eigenvecs,&power1,&beta,DMRG2DM,&power3);
   //2DM: Gamma_{ajkd} --> Gamma_{ajcd}
   for (int cnt=0; cnt<nOrbDMRG; cnt++){
      dgemm_(&notr,&notr,&power2,&power1,&power1,&alpha,DMRG2DM + cnt*power3,&power2,eigenvecs,&power1,&beta,work + cnt*power3,&power2);
   }
   //2DM: Gamma_{ajcd} --> Gamma_{abcd}
   for (int cnt=0; cnt<power2; cnt++){
      dgemm_(&notr,&notr,&power1,&power1,&power1,&alpha,work + cnt*power2,&power1,eigenvecs,&power1,&beta,DMRG2DM + cnt*power2,&power1);
   }

   //Update 1DM
   setDMRG1DM(N);

   //Rotate unitary...
   int passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){

      if (NDMRG[irrep] > 1){

         int rotationlinsize = NDMRG[irrep];
         int blocklinsize = OrbPerIrrep[irrep];

         double * temp1 = work;
         double * temp2 = work + rotationlinsize*blocklinsize;
         double * BlockEigen = eigenvecs + passed * (nOrbDMRG + 1);
      
         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               temp1[row + rotationlinsize*col] = unitary[irrep][ Nocc[irrep] + row + blocklinsize * col ];
            }
         }

         dgemm_(&tran,&notr,&rotationlinsize,&blocklinsize,&rotationlinsize,&alpha,BlockEigen,&nOrbDMRG,temp1,&rotationlinsize,&beta,temp2,&rotationlinsize);

         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               unitary[irrep][ Nocc[irrep] + row + blocklinsize * col ] = temp2[row + rotationlinsize*col];
            }
         }

      }

      if (CheMPS2::CASSCF_debugPrint){

         int linsize = OrbPerIrrep[irrep];
         dgemm_(&tran,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep],&linsize,unitary[irrep],&linsize,&beta,work,&linsize);
         double value = 0.0;
         for (int cnt=0; cnt<linsize; cnt++){
            value += (work[cnt*(1+linsize)]-1.0) * (work[cnt*(1+linsize)]-1.0);
            for (int cnt2=cnt+1; cnt2<linsize; cnt2++){
               value += work[cnt + cnt2*linsize] * work[cnt + cnt2*linsize] + work[cnt2 + cnt*linsize] * work[cnt2 + cnt*linsize];
            }
         }
         value = sqrt(value);
         cout << "Two-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] - I = " << value << endl;

      }

      passed += NDMRG[irrep];

   }

}

void CheMPS2::CASSCF::updateUnitary(double * temp1, double * temp2){

   //Per irrep
   for (int irrep=0; irrep<numberOfIrreps; irrep++){

      int linsize = OrbPerIrrep[irrep];
      int size = linsize * linsize;
      
      if (linsize>1){ //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok
         
         //Construct the anti-symmetric x-matrix
         double * xblock = temp1;
         for (int cnt=0; cnt<size; cnt++){ xblock[cnt] = 0.0; }
         for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
            for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){
               xblock[ Nocc[irrep] + cntDMRG + linsize*cntOcc ] = xmatrix[irrep][0][ cntDMRG + NDMRG[irrep] * cntOcc];
               xblock[ cntOcc + linsize*(Nocc[irrep] + cntDMRG) ] = - xblock[ Nocc[irrep] + cntDMRG + linsize*cntOcc ];
            }
         }
         for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){
            for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){
               xblock[ Nocc[irrep] + NDMRG[irrep] + cntVirt + linsize*(Nocc[irrep] + cntDMRG )] = xmatrix[irrep][1][ cntVirt + Nvirt[irrep] * cntDMRG ];
               xblock[ Nocc[irrep] + cntDMRG + linsize*(Nocc[irrep] + NDMRG[irrep] + cntVirt )]
                   = - xblock[ Nocc[irrep] + NDMRG[irrep] + cntVirt + linsize*( Nocc[irrep] + cntDMRG ) ];
            }
         }
         for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
            for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){
               xblock[ Nocc[irrep] + NDMRG[irrep] + cntVirt + linsize*cntOcc ] = xmatrix[irrep][2][ cntVirt + Nvirt[irrep] * cntOcc ];
               xblock[ cntOcc + linsize*( Nocc[irrep] + NDMRG[irrep] + cntVirt ) ]
                   = - xblock[ Nocc[irrep] + NDMRG[irrep] + cntVirt + linsize*cntOcc ];
            }
         }
         
         //Calculate its eigenvalues and eigenvectors
         double * Bmat = temp2;
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0; //SET !!!
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,xblock,&linsize,&beta,Bmat,&linsize); //Bmat = xblock * xblock
         
         char uplo = 'U';
         char jobz = 'V';
         double * eigenval = temp1 + 2*size;
         double * work = temp1 + size;
         int info;
         dsyev_(&jobz, &uplo, &linsize, Bmat, &linsize, eigenval, work, &size, &info); // xblock * xblock = Bmat * eigenval * Bmat^T
         
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,Bmat,&linsize,&beta,work,&linsize);
         char trans = 'T';
         double * work2 = temp2 + size;
         dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work,&linsize,&beta,work2,&linsize); //work2 = Bmat^T * xblock * Bmat
         
         if (CheMPS2::CASSCF_debugPrint){
            cout << "lambdas of irrep block " << irrep << " : " << endl;
            for (int cnt=0; cnt<linsize/2; cnt++){
               cout << "   block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
               cout << "           [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
            }
         }
         
         for (int cnt=0; cnt<linsize/2; cnt++){
            eigenval[cnt] = 0.5*( work2[2*cnt + linsize*(2*cnt+1)] - work2[2*cnt+1 + linsize*(2*cnt)] );
            work2[2*cnt + linsize*(2*cnt+1)] -= eigenval[cnt];
            work2[2*cnt+1 + linsize*(2*cnt)] += eigenval[cnt];
         }
         
         if (CheMPS2::CASSCF_debugPrint){
            double TwoNormResidual = 0.0;
            for (int cnt=0; cnt<size; cnt++){ TwoNormResidual += work2[cnt] * work2[cnt]; }
            TwoNormResidual = sqrt(TwoNormResidual);
            cout << "TwoNormResidual of irrep block " << irrep << " = " << TwoNormResidual << endl;
         }
         
         //Calculate exp(x)
         for (int cnt=0; cnt<size; cnt++){ work2[cnt] = 0.0; }
         for (int cnt=0; cnt<linsize/2; cnt++){
            double cosine = cos(eigenval[cnt]);
            double sine = sin(eigenval[cnt]);
            work2[2*cnt   + linsize*(2*cnt  )] = cosine;
            work2[2*cnt+1 + linsize*(2*cnt+1)] = cosine;
            work2[2*cnt   + linsize*(2*cnt+1)] = sine;
            work2[2*cnt+1 + linsize*(2*cnt  )] = - sine;
         }
         for (int cnt=2*(linsize/2); cnt<linsize; cnt++){
            work2[cnt*(linsize + 1)] = 1.0;
         }
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work2,&linsize,&beta,work,&linsize);
         dgemm_(&notr,&trans,&linsize,&linsize,&linsize,&alpha,work,&linsize,Bmat,&linsize,&beta,work2,&linsize); //work2 = exp(xblock)
         
         //U <-- exp(x) * U
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,work2,&linsize,unitary[irrep],&linsize,&beta,work,&linsize);
         int inc = 1;
         dcopy_(&size, work, &inc, unitary[irrep], &inc);
         
         //How unitary is the unitary
         if (CheMPS2::CASSCF_debugPrint){
         
            dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep],&linsize,unitary[irrep],&linsize,&beta,work2,&linsize);
            double value = 0.0;
            for (int cnt=0; cnt<linsize; cnt++){
               value += (work2[cnt*(1+linsize)]-1.0) * (work2[cnt*(1+linsize)]-1.0);
               for (int cnt2=cnt+1; cnt2<linsize; cnt2++){
                  value += work2[cnt + cnt2*linsize] * work2[cnt + cnt2*linsize] + work2[cnt2 + cnt*linsize] * work2[cnt2 + cnt*linsize];
               }
            }
            value = sqrt(value);
            cout << "Two-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] - I = " << value << endl;
         
         }
      
      }
   
   }

}

void CheMPS2::CASSCF::fillHamDMRG(Hamiltonian * HamDMRG){

   //Calculate the constant part of the energy.
   double Econst = HamRotated->getEconst();
   for (int cnt=0; cnt<nCondensed; cnt++){
      Econst += 2*HamRotated->getTmat(listCondensed[cnt],listCondensed[cnt]);
      for (int cnt2=0; cnt2<nCondensed; cnt2++){
         Econst += 2*HamRotated->getVmat(listCondensed[cnt],listCondensed[cnt2],listCondensed[cnt],listCondensed[cnt2])
                 -   HamRotated->getVmat(listCondensed[cnt],listCondensed[cnt2],listCondensed[cnt2],listCondensed[cnt]);
      }
   }
   HamDMRG->setEconst(Econst);
   
   //Calculate the one-body and two-body matrix elements.
   for (int cnt=0; cnt<nOrbDMRG; cnt++){
      for (int cnt2=0; cnt2<nOrbDMRG; cnt2++){
         const int productSymm = SymmInfo.directProd(irrepOfEachDMRGOrbital[cnt],irrepOfEachDMRGOrbital[cnt2]);
         if (productSymm == SymmInfo.getTrivialIrrep()){
            double value = HamRotated->getTmat(listDMRG[cnt], listDMRG[cnt2]);
            for (int cnt3=0; cnt3<nCondensed; cnt3++){
               value += 2 * HamRotated->getVmat(listDMRG[cnt],listCondensed[cnt3],listDMRG[cnt2],listCondensed[cnt3])
                      -     HamRotated->getVmat(listDMRG[cnt],listCondensed[cnt3],listCondensed[cnt3],listDMRG[cnt2]);
            }
            HamDMRG->setTmat(cnt,cnt2,value);
         }
         for (int cnt3=0; cnt3<nOrbDMRG; cnt3++){
            for (int cnt4=0; cnt4<nOrbDMRG; cnt4++){
               if (productSymm == SymmInfo.directProd(irrepOfEachDMRGOrbital[cnt3],irrepOfEachDMRGOrbital[cnt4])){
                  HamDMRG->setVmat(cnt,cnt2,cnt3,cnt4 , HamRotated->getVmat(listDMRG[cnt],listDMRG[cnt2],listDMRG[cnt3],listDMRG[cnt4]));
               }
            }
         }
      }
   }

}

void CheMPS2::CASSCF::allocateAndFillOCC(int * DOCCin, int * SOCCin){
   
   DOCC = new int[numberOfIrreps];
   SOCC = new int[numberOfIrreps];
   
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
      DOCC[cnt] = DOCCin[cnt];
      SOCC[cnt] = SOCCin[cnt];
   }
   
   cout << "Norb = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << OrbPerIrrep[cnt] << " , "; }
   cout << OrbPerIrrep[numberOfIrreps-1] << " ]" << endl;
   cout << "DOCC = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << DOCC[cnt] << " , "; }
   cout << DOCC[numberOfIrreps-1] << " ]" << endl;
   cout << "SOCC = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << SOCC[cnt] << " , "; }
   cout << SOCC[numberOfIrreps-1] << " ]" << endl;
   
}

void CheMPS2::CASSCF::allocateAndFillOCC(const string filename){

   string line, part;
   int pos, pos2;
   
   ifstream inputfile(filename.c_str());
   
   //First go to the start of the integral dump.
   bool stop = false;
   string start = "****  Molecular Integrals For CheMPS Start Here";
   do{
      getline(inputfile,line);
      pos = line.find(start);
      if (pos==0) stop = true;
   } while (!stop);
   
   //Get the group name and convert it to the group number
   getline(inputfile,line);
   getline(inputfile,line);
   getline(inputfile,line);
   getline(inputfile,line);
   getline(inputfile,line);
   getline(inputfile,line);
   DOCC = new int[numberOfIrreps];
   SOCC = new int[numberOfIrreps];
   
   //Doubly occupied orbitals
   getline(inputfile,line);
   
   pos = line.find("=") + 3;
   pos2 = line.find(" ",pos);
   part = line.substr(pos, pos2-pos);
   
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){
      DOCC[cnt] = atoi(part.c_str());
      pos = pos2 + 2;
      pos2 = line.find(" ",pos);
      part = line.substr(pos,pos2-pos);
   }
   DOCC[numberOfIrreps-1] = atoi(part.c_str());
   
   //Singly occupied orbitals
   getline(inputfile,line);
   
   pos = line.find("=") + 3;
   pos2 = line.find(" ",pos);
   part = line.substr(pos, pos2-pos);
   
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){
      SOCC[cnt] = atoi(part.c_str());
      pos = pos2 + 2;
      pos2 = line.find(" ",pos);
      part = line.substr(pos,pos2-pos);
   }
   SOCC[numberOfIrreps-1] = atoi(part.c_str());
   
   cout << "Norb = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << OrbPerIrrep[cnt] << " , "; }
   cout << OrbPerIrrep[numberOfIrreps-1] << " ]" << endl;
   cout << "DOCC = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << DOCC[cnt] << " , "; }
   cout << DOCC[numberOfIrreps-1] << " ]" << endl;
   cout << "SOCC = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << SOCC[cnt] << " , "; }
   cout << SOCC[numberOfIrreps-1] << " ]" << endl;
   
   inputfile.close();
   
}

void CheMPS2::CASSCF::setupStart(int * NoccIn, int * NDMRGIn, int * NvirtIn){

   setupStartCalled = true;

   //Copy the partitioning
   Nocc  = new int[numberOfIrreps];
   NDMRG = new int[numberOfIrreps];
   Nvirt = new int[numberOfIrreps];
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
   
      if ((NoccIn[cnt]<0) || (NDMRGIn[cnt]<0) || (NvirtIn[cnt]<0)){ cout << "One of the entries for the setup was smaller than zero." << endl; }
      if (NoccIn[cnt] + NDMRGIn[cnt] + NvirtIn[cnt] != OrbPerIrrep[cnt]){ cout << "The sum of the partition for the setup doesn't match the number of orbitals." << endl; }
      
      Nocc[cnt] = NoccIn[cnt];
      NDMRG[cnt] = NDMRGIn[cnt];
      Nvirt[cnt] = NvirtIn[cnt];
   
   }
   
   //Number of orbitals for DMRG, their irrep, and their number in terms of the original Hamiltonian
   nOrbDMRG = 0;
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ nOrbDMRG += NDMRG[cnt]; }
   irrepOfEachDMRGOrbital = new int[nOrbDMRG];
   listDMRG = new int[nOrbDMRG];
   int passed = 0;
   int passed2 = 0;
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
      for (int cnt2=0; cnt2<NDMRG[cnt]; cnt2++){
         irrepOfEachDMRGOrbital[passed+cnt2] = cnt;
         listDMRG[passed+cnt2] = passed2 + Nocc[cnt] + cnt2;
      }
      passed += NDMRG[cnt];
      passed2 += OrbPerIrrep[cnt];
   }
   
   //How many and which orbitals are condensed now.
   nCondensed = 0;
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ nCondensed += Nocc[cnt]; }
   listCondensed = new int[nCondensed];
   passed = 0;
   int irrep = 0;
   int counter = 0;
   for (int cnt=0; cnt<L; cnt++){
      if (cnt-passed < Nocc[irrep]){
         listCondensed[counter] = cnt;
         counter++;
      }
      if ((cnt-passed == OrbPerIrrep[irrep]-1) && (irrep<numberOfIrreps-1)){
         passed += OrbPerIrrep[irrep];
         irrep++;
      }
   }
   
   //Allocate the xmatrix
   x_linearlength = 0;
   xmatrix = new double**[numberOfIrreps];
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      xmatrix[irrep] = new double*[3];
      for (int case_=0; case_<3; case_++){
         int size = 0;
         if (case_==0){ size = Nocc[irrep] * NDMRG[irrep]; }
         if (case_==1){ size = NDMRG[irrep] * Nvirt[irrep]; }
         if (case_==2){ size = Nocc[irrep] * Nvirt[irrep]; }
         xmatrix[irrep][case_] = new double[size];
         for (int cnt=0; cnt<size; cnt++){ xmatrix[irrep][case_][cnt] = 0.0; }
         x_linearlength += size;
      }
   }
   
   //Allocate the unitary
   unitary = new double*[numberOfIrreps];
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      int size = OrbPerIrrep[irrep] * OrbPerIrrep[irrep];
      unitary[irrep] = new double[size];
      for (int cnt=0; cnt<size; cnt++){ unitary[irrep][cnt] = 0.0; }
      for (int cnt=0; cnt<OrbPerIrrep[irrep]; cnt++){ unitary[irrep][cnt*(1+OrbPerIrrep[irrep])] = 1.0; }
   }
   
   //Allocate jumpsHamOrig
   jumpsHamOrig = new int[numberOfIrreps+1];
   jumpsHamOrig[0] = 0;
   for (int cnt=0; cnt<numberOfIrreps; cnt++){ jumpsHamOrig[cnt+1] = jumpsHamOrig[cnt] + OrbPerIrrep[cnt]; }
   
   //Allocate space for the DMRG 1DM and 2DM
   DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
   DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];
   
   //Find the corresponding indices
   x_firstindex = new int[x_linearlength];
   x_secondindex = new int[x_linearlength];
   x_linearlength = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      for (int case_=0; case_<3; case_++){
         if (case_==0){
            for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
               for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){ //DMRG is the row index, hence fast moving one
                  x_firstindex[x_linearlength] = jumpsHamOrig[irrep] + Nocc[irrep] + cntDMRG;
                  x_secondindex[x_linearlength] = jumpsHamOrig[irrep] + cntOcc;
                  x_linearlength++;
               }
            }
         }
         if (case_==1){
            for (int cntDMRG=0; cntDMRG<NDMRG[irrep]; cntDMRG++){
               for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength] = jumpsHamOrig[irrep] + Nocc[irrep] + NDMRG[irrep] + cntVirt;
                  x_secondindex[x_linearlength] = jumpsHamOrig[irrep] + Nocc[irrep] + cntDMRG;
                  x_linearlength++;
               }
            }
         }
         if (case_==2){
            for (int cntOcc=0; cntOcc<Nocc[irrep]; cntOcc++){
               for (int cntVirt=0; cntVirt<Nvirt[irrep]; cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength] = jumpsHamOrig[irrep] + Nocc[irrep] + NDMRG[irrep] + cntVirt;
                  x_secondindex[x_linearlength] = jumpsHamOrig[irrep] + cntOcc;
                  x_linearlength++;
               }
            }
         }
      }
   }
   
   //To calculate the F-matrix elements only once, and to store them for future access
   Fmatrix = new double*[numberOfIrreps];
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
      Fmatrix[cnt] = new double[OrbPerIrrep[cnt] * OrbPerIrrep[cnt]];
      for (int cnt2=0; cnt2<OrbPerIrrep[cnt] * OrbPerIrrep[cnt]; cnt2++){ Fmatrix[cnt][cnt2] = 0.0; }
   }
   
   //Print what we have just set up.
   cout << "Nocc  = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << Nocc[cnt] << " , "; }
   cout << Nocc[numberOfIrreps-1] << " ]" << endl;
   
   cout << "Ndmrg = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << NDMRG[cnt] << " , "; }
   cout << NDMRG[numberOfIrreps-1] << " ]" << endl;
   
   cout << "Nvirt = [ ";
   for (int cnt=0; cnt<numberOfIrreps-1; cnt++){ cout << Nvirt[cnt] << " , "; }
   cout << Nvirt[numberOfIrreps-1] << " ]" << endl;
   
   cout << "Condensed orbitals = [ ";
   for (int cnt=0; cnt<nCondensed-1; cnt++){ cout << listCondensed[cnt] << " , "; }
   cout << listCondensed[nCondensed-1] << " ]" << endl;
      
   cout << "DMRG orbital irreps = [ ";
   for (int cnt=0; cnt<nOrbDMRG-1; cnt++){ cout << irrepOfEachDMRGOrbital[cnt] << " , "; }
   cout << irrepOfEachDMRGOrbital[nOrbDMRG-1] << " ]" << endl;
      
   cout << "DMRG orbital numbers = [ ";
   for (int cnt=0; cnt<nOrbDMRG-1; cnt++){ cout << listDMRG[cnt] << " , "; }
   cout << listDMRG[nOrbDMRG-1] << " ]" << endl;
   
   cout << "Number of variables in the x-matrix = " << x_linearlength << endl;

}

