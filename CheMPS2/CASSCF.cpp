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

   allocateAndFillOCC(filename);
   
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

   allocateAndFillOCC(DOCCin, SOCCin);
   
   setupStartCalled = false;
   
}

CheMPS2::CASSCF::~CASSCF(){

   delete HamOrig;
   delete HamRotated;
   
   delete [] DOCC;
   delete [] SOCC;
   
   if (setupStartCalled){
   
      delete unitary; //First delete the unitary as it requires the iHandler in its destructor.
      delete iHandler;

      delete [] DMRG1DM;
      delete [] DMRG2DM;
      
      for (int irrep=0; irrep<numberOfIrreps; irrep++){
         delete [] Fmatrix[irrep];
         delete [] QmatrixOCC[irrep];
         delete [] QmatrixACT[irrep];
      }
      delete [] Fmatrix;
      delete [] QmatrixOCC;
      delete [] QmatrixACT;
      
   }

}

int CheMPS2::CASSCF::getNumberOfIrreps(){ return numberOfIrreps; }

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

      int NDMRG = iHandler->getNDMRG(irrep);
      if (NDMRG > 0){

         //Calculate the eigenvectors and values per block
         dsyev_(&jobz, &uplo, &NDMRG, eigenvecs + passed*(1+nOrbDMRG) ,&nOrbDMRG, eigenval + passed, work, &size, &info);

         //Print the NOON
         cout << "CASSCF :: DMRG 1DM eigenvalues [NOON] of irrep " << irrep << " = [ ";
         for (int cnt=0; cnt<NDMRG-1; cnt++){ cout << eigenval[passed + NDMRG-1-cnt] << " , "; }
         cout << eigenval[passed + 0] << " ]." << endl;

         //Sort the eigenvecs
         for (int col=0; col<NDMRG/2; col++){
            for (int row=0; row<NDMRG; row++){
               double temp = eigenvecs[passed + row + nOrbDMRG * (passed + NDMRG - 1 - col)];
               eigenvecs[passed + row + nOrbDMRG * (passed + NDMRG - 1 - col)] = eigenvecs[passed + row + nOrbDMRG * (passed + col)];
               eigenvecs[passed + row + nOrbDMRG * (passed + col)] = temp;
            }
         }

      }

      //Update the number of passed DMRG orbitals
      passed += NDMRG;

   }

   delete [] work;
   delete [] eigenval;

}

void CheMPS2::CASSCF::rotate2DMand1DM(const int N, double * eigenvecs, double * work){

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

}

void CheMPS2::CASSCF::buildQmatrixOCC(){

   for (int irrep_ij=0; irrep_ij<numberOfIrreps; irrep_ij++){
      for (int indexi=0; indexi<iHandler->getNORB(irrep_ij); indexi++){
         const int origIndexI = indexi + iHandler->getOrigNOCCstart(irrep_ij);
         for (int indexj=indexi; indexj<iHandler->getNORB(irrep_ij); indexj++){
            const int origIndexJ = indexj + iHandler->getOrigNOCCstart(irrep_ij);
            double value = 0.0;
            for (int irrep_occ=0; irrep_occ<numberOfIrreps; irrep_occ++){
               for (int origIndexOcc=iHandler->getOrigNOCCstart(irrep_occ); origIndexOcc<iHandler->getOrigNDMRGstart(irrep_occ); origIndexOcc++){
                  value += 2 * HamRotated->getVmat(origIndexI, origIndexOcc, origIndexJ, origIndexOcc)
                             - HamRotated->getVmat(origIndexI, origIndexJ, origIndexOcc, origIndexOcc);
               }
            }
            QmatrixOCC[irrep_ij][indexi + iHandler->getNORB(irrep_ij) * indexj] = value;
            QmatrixOCC[irrep_ij][indexj + iHandler->getNORB(irrep_ij) * indexi] = value;
         }
      }
   }

}

double CheMPS2::CASSCF::QmatOCC(const int index1, const int index2) const{

   int irrep1 = 0;
   while (index1 >= iHandler->getOrigNOCCstart(irrep1+1)){ irrep1++; }
   int irrep2 = 0;
   while (index2 >= iHandler->getOrigNOCCstart(irrep2+1)){ irrep2++; }
   
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   return QmatrixOCC[irrep1][ index1 - iHandler->getOrigNOCCstart(irrep1) + iHandler->getNORB(irrep1) * ( index2 - iHandler->getOrigNOCCstart(irrep1) ) ];

}

void CheMPS2::CASSCF::buildQmatrixACT(){

   for (int irrep_ij=0; irrep_ij<numberOfIrreps; irrep_ij++){
      for (int indexi=0; indexi<iHandler->getNORB(irrep_ij); indexi++){
         const int origIndexI = indexi + iHandler->getOrigNOCCstart(irrep_ij);
         for (int indexj=indexi; indexj<iHandler->getNORB(irrep_ij); indexj++){
            const int origIndexJ = indexj + iHandler->getOrigNOCCstart(irrep_ij);
            double value = 0.0;
            for (int irrep_act=0; irrep_act<numberOfIrreps; irrep_act++){
               for (int origIndexAct1=iHandler->getOrigNDMRGstart(irrep_act); origIndexAct1<iHandler->getOrigNVIRTstart(irrep_act); origIndexAct1++){
                  const int DMRGindex1 = iHandler->getDMRGcumulative(irrep_act) + origIndexAct1 - iHandler->getOrigNDMRGstart(irrep_act);
                  for (int origIndexAct2=iHandler->getOrigNDMRGstart(irrep_act); origIndexAct2<iHandler->getOrigNVIRTstart(irrep_act); origIndexAct2++){
                     const int DMRGindex2 = iHandler->getDMRGcumulative(irrep_act) + origIndexAct2 - iHandler->getOrigNDMRGstart(irrep_act);
                     value += DMRG1DM[ DMRGindex1 + nOrbDMRG * DMRGindex2 ] * ( HamRotated->getVmat(origIndexI, origIndexAct1, origIndexJ, origIndexAct2)
                                                                        - 0.5 * HamRotated->getVmat(origIndexI, origIndexJ, origIndexAct1, origIndexAct2) );
                  }
               }
            }
            QmatrixACT[irrep_ij][indexi + iHandler->getNORB(irrep_ij) * indexj] = value;
            QmatrixACT[irrep_ij][indexj + iHandler->getNORB(irrep_ij) * indexi] = value;
         }
      }
   }

}

double CheMPS2::CASSCF::QmatACT(const int index1, const int index2) const{

   int irrep1 = 0;
   while (index1 >= iHandler->getOrigNOCCstart(irrep1+1)){ irrep1++; }
   int irrep2 = 0;
   while (index2 >= iHandler->getOrigNOCCstart(irrep2+1)){ irrep2++; }
   
   if (irrep1 != irrep2){ return 0.0; } //From now on: both irreps are the same.
   
   return QmatrixACT[irrep1][ index1 - iHandler->getOrigNOCCstart(irrep1) + iHandler->getNORB(irrep1) * ( index2 - iHandler->getOrigNOCCstart(irrep1) ) ];

}

void CheMPS2::CASSCF::fillHamDMRG(Hamiltonian * HamDMRG){

   //Calculate the constant part of the energy.
   double Econst = HamRotated->getEconst();
   for (int irrep1=0; irrep1<numberOfIrreps; irrep1++){
      for (int orb1=iHandler->getOrigNOCCstart(irrep1); orb1<iHandler->getOrigNDMRGstart(irrep1); orb1++){
         Econst += 2 * HamRotated->getTmat(orb1,orb1) + QmatOCC(orb1,orb1);
      }
   }
   HamDMRG->setEconst(Econst);
   
   for (int irrep1=0; irrep1<numberOfIrreps; irrep1++){
      //Calculate the one-body matrix elements.
      for (int orb1=iHandler->getOrigNDMRGstart(irrep1); orb1<iHandler->getOrigNVIRTstart(irrep1); orb1++){
         const int DMRGorb1 = iHandler->getDMRGcumulative(irrep1) + orb1 - iHandler->getOrigNDMRGstart(irrep1);
         for (int orb2=orb1; orb2<iHandler->getOrigNVIRTstart(irrep1); orb2++){
            double value = HamRotated->getTmat(orb1, orb2) + QmatOCC(orb1, orb2);
            const int DMRGorb2 = iHandler->getDMRGcumulative(irrep1) + orb2 - iHandler->getOrigNDMRGstart(irrep1);
            HamDMRG->setTmat(DMRGorb1, DMRGorb2, value);
         }
      }
      
      //Calculate the two-body matrix elements.
      for (int irrep2=irrep1; irrep2<numberOfIrreps; irrep2++){
         const int irrep_product = SymmInfo.directProd(irrep1, irrep2);
         for (int irrep3=irrep1; irrep3<numberOfIrreps; irrep3++){
            const int irrep4 = SymmInfo.directProd(irrep_product, irrep3);
            if (irrep4 >= irrep2){
               for (int orb1=iHandler->getOrigNDMRGstart(irrep1); orb1<iHandler->getOrigNVIRTstart(irrep1); orb1++){
                  const int DMRGorb1 = iHandler->getDMRGcumulative(irrep1) + orb1 - iHandler->getOrigNDMRGstart(irrep1);
                  for (int orb2=iHandler->getOrigNDMRGstart(irrep2); orb2<iHandler->getOrigNVIRTstart(irrep2); orb2++){
                     const int DMRGorb2 = iHandler->getDMRGcumulative(irrep2) + orb2 - iHandler->getOrigNDMRGstart(irrep2);
                     for (int orb3=iHandler->getOrigNDMRGstart(irrep3); orb3<iHandler->getOrigNVIRTstart(irrep3); orb3++){
                        const int DMRGorb3 = iHandler->getDMRGcumulative(irrep3) + orb3 - iHandler->getOrigNDMRGstart(irrep3);
                        for (int orb4=iHandler->getOrigNDMRGstart(irrep4); orb4<iHandler->getOrigNVIRTstart(irrep4); orb4++){
                           const int DMRGorb4 = iHandler->getDMRGcumulative(irrep4) + orb4 - iHandler->getOrigNDMRGstart(irrep4);
                           HamDMRG->setVmat(DMRGorb1, DMRGorb2, DMRGorb3, DMRGorb4, HamRotated->getVmat(orb1, orb2, orb3, orb4));
                        }
                     }
                  }
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
   unitary = new DMRGSCFunitary(iHandler);
   
   //Allocate space for the DMRG 1DM and 2DM
   nOrbDMRG = iHandler->getDMRGcumulative(numberOfIrreps);
   DMRG1DM = new double[nOrbDMRG * nOrbDMRG];
   DMRG2DM = new double[nOrbDMRG * nOrbDMRG * nOrbDMRG * nOrbDMRG];
   
   //To calculate the F-matrix and Q-matrix(occ,act) elements only once, and to store them for future access
   Fmatrix    = new double*[numberOfIrreps];
   QmatrixOCC = new double*[numberOfIrreps];
   QmatrixACT = new double*[numberOfIrreps];
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      const int size = iHandler->getNORB(irrep) * iHandler->getNORB(irrep);
      Fmatrix[irrep]    = new double[size];
      QmatrixOCC[irrep] = new double[size];
      QmatrixACT[irrep] = new double[size];
      for (int cnt2=0; cnt2<size; cnt2++){
         Fmatrix[irrep][cnt2]    = 0.0;
         QmatrixOCC[irrep][cnt2] = 0.0;
         QmatrixACT[irrep][cnt2] = 0.0;
      }
   }
   
   //Print the MO info. This requires the indexHandler to be created...
   checkHF();
   
   //Print what we have just set up.
   iHandler->Print();
   
   cout << "Number of variables in the x-matrix = " << unitary->getNumVariablesX() << endl;

}

