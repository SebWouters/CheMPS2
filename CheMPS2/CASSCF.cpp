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

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "CASSCF.h"
#include "Lapack.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

CheMPS2::CASSCF::CASSCF(const string filename){

   HamOrig = new Hamiltonian(filename);
   shouldHamOrigBeDeleted = true;
   
   L = HamOrig->getL();
   SymmInfo.setGroup(HamOrig->getNGroup());
   numberOfIrreps = SymmInfo.getNumberOfIrreps();
   
   allocateAndFillOCC(filename);
   setupStartCalled = false;

}

CheMPS2::CASSCF::CASSCF(Hamiltonian * HamIn, int * DOCCin, int * SOCCin){

   HamOrig = HamIn;
   shouldHamOrigBeDeleted = false;
   
   L = HamOrig->getL();
   SymmInfo.setGroup(HamOrig->getNGroup());
   numberOfIrreps = SymmInfo.getNumberOfIrreps();
   
   allocateAndFillOCC(DOCCin, SOCCin);
   setupStartCalled = false;

}

CheMPS2::CASSCF::~CASSCF(){

   if (shouldHamOrigBeDeleted){ delete HamOrig; }
   
   delete [] DOCC;
   delete [] SOCC;
   
   if (setupStartCalled){
   
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

}

int CheMPS2::CASSCF::getNumberOfIrreps(){ return numberOfIrreps; }

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

void CheMPS2::CASSCF::setDMRG1DM(const int nDMRGelectrons, const int totOrbDMRG, double * localDMRG1DM, double * localDMRG2DM){

   const double prefactor = 1.0/(nDMRGelectrons-1.0);

   for (int cnt1=0; cnt1<totOrbDMRG; cnt1++){
      for (int cnt2=cnt1; cnt2<totOrbDMRG; cnt2++){
         localDMRG1DM[cnt1 + totOrbDMRG*cnt2] = 0.0;
         for (int cnt3=0; cnt3<totOrbDMRG; cnt3++){ localDMRG1DM[cnt1 + totOrbDMRG*cnt2] += localDMRG2DM[cnt1 + totOrbDMRG * (cnt3 + totOrbDMRG * (cnt2 + totOrbDMRG * cnt3 ) ) ]; }
         localDMRG1DM[cnt1 + totOrbDMRG*cnt2] *= prefactor;
         localDMRG1DM[cnt2 + totOrbDMRG*cnt1] = localDMRG1DM[cnt1 + totOrbDMRG*cnt2];
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

void CheMPS2::CASSCF::calcNOON(DMRGSCFindices * localIdx, double * eigenvecs, double * workmem, double * localDMRG1DM){

   const int numIrreps = localIdx->getNirreps();
   int totOrbDMRG = localIdx->getDMRGcumulative( numIrreps );
   int size = totOrbDMRG * totOrbDMRG;
   double * eigenval = workmem + size;

   for (int cnt=0; cnt<size; cnt++){ eigenvecs[cnt] = localDMRG1DM[cnt]; }

   char jobz = 'V';
   char uplo = 'U';
   int info;
   int passed = 0;
   for (int irrep=0; irrep<numIrreps; irrep++){

      int NDMRG = localIdx->getNDMRG(irrep);
      if (NDMRG > 0){

         //Calculate the eigenvectors and values per block
         dsyev_(&jobz, &uplo, &NDMRG, eigenvecs + passed*(1+totOrbDMRG), &totOrbDMRG, eigenval + passed, workmem, &size, &info);

         //Print the NOON
         if (irrep==0){ cout << "DMRGSCF::calcNOON : DMRG 1DM eigenvalues [NOON] of irrep " << irrep << " = [ "; }
         else {         cout << "                    DMRG 1DM eigenvalues [NOON] of irrep " << irrep << " = [ "; }
         for (int cnt=0; cnt<NDMRG-1; cnt++){ cout << eigenval[passed + NDMRG-1-cnt] << " , "; }
         cout << eigenval[passed + 0] << " ]." << endl;

         //Sort the eigenvecs
         for (int col=0; col<NDMRG/2; col++){
            for (int row=0; row<NDMRG; row++){
               double temp = eigenvecs[passed + row + totOrbDMRG * (passed + NDMRG - 1 - col)];
               eigenvecs[passed + row + totOrbDMRG * (passed + NDMRG - 1 - col)] = eigenvecs[passed + row + totOrbDMRG * (passed + col)];
               eigenvecs[passed + row + totOrbDMRG * (passed + col)] = temp;
            }
         }

      }

      //Update the number of passed DMRG orbitals
      passed += NDMRG;

   }

}

void CheMPS2::CASSCF::rotate2DMand1DM(const int nDMRGelectrons, int totOrbDMRG, double * eigenvecs, double * work, double * localDMRG1DM, double * localDMRG2DM){

   char notr = 'N';
   char tran = 'T';
   double alpha = 1.0;
   double beta = 0.0;

   int power1 = totOrbDMRG;
   int power2 = totOrbDMRG*totOrbDMRG;
   int power3 = totOrbDMRG*totOrbDMRG*totOrbDMRG;

   //2DM: Gamma_{ijkl} --> Gamma_{ajkl}
   dgemm_(&tran,&notr,&power1,&power3,&power1,&alpha,eigenvecs,&power1,localDMRG2DM,&power1,&beta,work,&power1);
   //2DM: Gamma_{ajkl} --> Gamma_{ajkd}
   dgemm_(&notr,&notr,&power3,&power1,&power1,&alpha,work,&power3,eigenvecs,&power1,&beta,localDMRG2DM,&power3);
   //2DM: Gamma_{ajkd} --> Gamma_{ajcd}
   for (int cnt=0; cnt<totOrbDMRG; cnt++){
      dgemm_(&notr,&notr,&power2,&power1,&power1,&alpha,localDMRG2DM + cnt*power3,&power2,eigenvecs,&power1,&beta,work + cnt*power3,&power2);
   }
   //2DM: Gamma_{ajcd} --> Gamma_{abcd}
   for (int cnt=0; cnt<power2; cnt++){
      dgemm_(&notr,&notr,&power1,&power1,&power1,&alpha,work + cnt*power2,&power1,eigenvecs,&power1,&beta,localDMRG2DM + cnt*power2,&power1);
   }

   //Update 1DM
   setDMRG1DM(nDMRGelectrons, totOrbDMRG, localDMRG1DM, localDMRG2DM);

}

void CheMPS2::CASSCF::rotateOldToNew(DMRGSCFmatrix * myMatrix){

   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
   
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

  for (int irrepQ = 0; irrepQ < numberOfIrreps; irrepQ++){
   
      const int linearsizeQ = iHandler->getNORB(irrepQ);
      const int numberOfUniqueIndices = (linearsizeQ * (linearsizeQ + 1))/2;
      
      #pragma omp parallel for schedule(static)
      for (int combinedindex = 0; combinedindex < numberOfUniqueIndices; combinedindex++){
      
         int colQ = 1;
         while ( (colQ*(colQ+1))/2 <= combinedindex ){ colQ++; }
         colQ -= 1;
         int rowQ = combinedindex - (colQ*(colQ+1))/2;
         
         const int HamIndexI = iHandler->getOrigNOCCstart(irrepQ) + rowQ;
         const int HamIndexJ = iHandler->getOrigNOCCstart(irrepQ) + colQ;
         
         double theValue = 0.0;
         
         for (int irrepN = 0; irrepN < numberOfIrreps; irrepN++){
            const int linearsizeN = iHandler->getNORB( irrepN );
            for (int rowN = 0; rowN < linearsizeN; rowN++){
            
               const int HamIndexS = iHandler->getOrigNOCCstart( irrepN ) + rowN;
               theValue += densityMatrix->get(irrepN, rowN, rowN) * ( HamOrig->getVmat(HamIndexI,HamIndexS,HamIndexJ,HamIndexS)
                                                              - 0.5 * HamOrig->getVmat(HamIndexI,HamIndexJ,HamIndexS,HamIndexS) );
               
               for (int colN = rowN+1; colN < linearsizeN; colN++){
               
                  const int HamIndexT = iHandler->getOrigNOCCstart( irrepN ) + colN;
                  theValue += densityMatrix->get(irrepN, rowN, colN) * ( 2 * HamOrig->getVmat(HamIndexI,HamIndexS,HamIndexJ,HamIndexT)
                                                                     - 0.5 * HamOrig->getVmat(HamIndexI,HamIndexJ,HamIndexS,HamIndexT) 
                                                                     - 0.5 * HamOrig->getVmat(HamIndexI,HamIndexJ,HamIndexT,HamIndexS) );
               
               }
            }
         }
         
         resultMatrix->set( irrepQ, rowQ, colQ, theValue );
         resultMatrix->set( irrepQ, colQ, rowQ, theValue );
      
      }
   }

}

void CheMPS2::CASSCF::buildQmatOCC(){

   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
   
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

   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
   
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

   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
      const int NumORB = iHandler->getNORB(irrep);
      for (int row = 0; row < NumORB; row++){
         const int HamIndexRow = iHandler->getOrigNOCCstart(irrep) + row;
         for (int col = 0; col < NumORB; col++){
            const int HamIndexCol = iHandler->getOrigNOCCstart(irrep) + col;
            theTmatrix->set( irrep, row, col, HamOrig->getTmat(HamIndexRow, HamIndexCol) );
         }
      }
   }
   
   rotateOldToNew( theTmatrix );

}

void CheMPS2::CASSCF::fillConstAndTmatDMRG(Hamiltonian * HamDMRG) const{

   //Constant part of the energy
   double value = HamOrig->getEconst();
   for (int irrep = 0; irrep < numberOfIrreps; irrep++){
      for (int orb = 0; orb < iHandler->getNOCC(irrep); orb++){
         value += 2 * theTmatrix->get(irrep, orb, orb) + theQmatOCC->get(irrep, orb, orb);
      }
   }
   HamDMRG->setEconst(value);
   
   //One-body terms: diagonal in the irreps
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
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

void CheMPS2::CASSCF::allocateAndFillOCC(int * DOCCin, int * SOCCin){
   
   DOCC = new int[numberOfIrreps];
   SOCC = new int[numberOfIrreps];
   
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
      DOCC[cnt] = DOCCin[cnt];
      SOCC[cnt] = SOCCin[cnt];
   }
   
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
   
   iHandler = new DMRGSCFindices(L, SymmInfo.getGroupNumber(), NoccIn, NDMRGIn, NvirtIn);
   unitary  = new DMRGSCFunitary(iHandler);
   theDIIS = NULL;
   theRotatedTEI = new DMRGSCFintegrals( iHandler );
   
   //Allocate space for the DMRG 1DM and 2DM
   nOrbDMRG = iHandler->getDMRGcumulative(numberOfIrreps);
   DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
   DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];
   
   //To calculate the F-matrix and Q-matrix(occ,act) elements only once, and to store them for future access
   theFmatrix = new DMRGSCFmatrix( iHandler ); theFmatrix->clear();
   theQmatOCC = new DMRGSCFmatrix( iHandler ); theQmatOCC->clear();
   theQmatACT = new DMRGSCFmatrix( iHandler ); theQmatACT->clear();
   theQmatWORK= new DMRGSCFmatrix( iHandler );theQmatWORK->clear();
   theTmatrix = new DMRGSCFmatrix( iHandler ); theTmatrix->clear();
   
   //To calculate the w_tilde elements only once, and store them for future access
   wmattilde = new DMRGSCFwtilde( iHandler );
   
   //Print the MO info. This requires the indexHandler to be created...
   checkHF();
   
   //Print what we have just set up.
   iHandler->Print();
   
   cout << "DMRGSCF::setupStart : Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;

}

