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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "Heff.h"
#include "Davidson.h"
#include "Lapack.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::Heff::Heff(const SyBookkeeper * denBKIn, const Problem * ProbIn){

   denBK = denBKIn;
   Prob = ProbIn;

}

CheMPS2::Heff::~Heff(){

}

void CheMPS2::Heff::makeHeff(double * memS, double * memHeff, const Sobject * denS, TensorL *** Ltensors, TensorA **** Atensors, TensorB **** Btensors, TensorC **** Ctensors, TensorD **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   const int DIM = max(denBK->gMaxDimAtBound(indexS), denBK->gMaxDimAtBound(indexS+2));
   
   //PARALLEL
   #pragma omp parallel
   {
   
      double * temp  = new double[DIM*DIM];
      double * temp2 = new double[DIM*DIM];
   
      #pragma omp for schedule(dynamic)
      for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
      
         const int ikappa = denS->gReorder(ikappaBIS);

         for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeff[cnt] = 0.0; }
         
         addDiagram1C(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
         addDiagram1D(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
         addDiagram2dall(ikappa, memS, memHeff, denS);
         addDiagram3Eand3H(ikappa, memS, memHeff, denS);
         addDiagramExcitations(ikappa, memS, memHeff, denS, nLower, VeffTilde);
         
         if (!atLeft){

            addDiagram1A(ikappa, memS, memHeff, denS, Xtensors[indexS-1]);
            
            addDiagram2b1and2b2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][0]);
            addDiagram2c1and2c2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][1]);
            addDiagram2b3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][0]);
            addDiagram2c3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][1]);
            addDiagram2b3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][0]);
            addDiagram2c3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][1]);
            
            addDiagram3Aand3D(ikappa, memS, memHeff, denS, Qtensors[indexS-1][0], Ltensors[indexS-1], temp);
            addDiagram3Band3I(ikappa, memS, memHeff, denS, Qtensors[indexS-1][1], Ltensors[indexS-1], temp);
            
            addDiagram4A1and4A2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1][1][0]);
            addDiagram4A1and4A2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1][1][0]);
            addDiagram4A3and4A4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][1][0]);
            addDiagram4A3and4A4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][1][0]);
            addDiagram4D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp);
            addDiagram4I(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp);
         
         }
         
         if (!atRight){
         
            addDiagram1B(ikappa, memS, memHeff, denS, Xtensors[indexS+1]);
            
            addDiagram2e1and2e2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][1]);
            addDiagram2f1and2f2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][0]);
            addDiagram2e3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][1]);
            addDiagram2f3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][0]);
            addDiagram2e3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][1]);
            addDiagram2f3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][0]);
            
            addDiagram3Kand3F(ikappa, memS, memHeff, denS, Qtensors[indexS+1][1], Ltensors[indexS+1], temp);
            addDiagram3Land3G(ikappa, memS, memHeff, denS, Qtensors[indexS+1][0], Ltensors[indexS+1], temp);
            
            addDiagram4F(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp);
            addDiagram4G(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp);
            addDiagram4J1and4J2spin0(ikappa, memS, memHeff, denS, Atensors[indexS+1][1][0]);
            addDiagram4J1and4J2spin1(ikappa, memS, memHeff, denS, Btensors[indexS+1][1][0]);
            addDiagram4J3and4J4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][1][0]);
            addDiagram4J3and4J4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][1][0]);
            
         }
         
         if ((!atLeft) && (!atRight)){
         
            addDiagram2a1spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp);
            addDiagram2a2spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp);
            addDiagram2a1spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp);
            addDiagram2a2spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp);
            addDiagram2a3spin0(ikappa, memS, memHeff, denS, Ctensors, F0tensors, temp);
            addDiagram2a3spin1(ikappa, memS, memHeff, denS, Dtensors, F1tensors, temp);
            
            addDiagram3C(ikappa, memS, memHeff, denS, Qtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram3J(ikappa, memS, memHeff, denS, Qtensors[indexS+1], Ltensors[indexS-1], temp);
            
            addDiagram4B1and4B2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B1and4B2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B3and4B4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B3and4B4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C1and4C2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C1and4C2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C3and4C4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C3and4C4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram4H(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram4K1and4K2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp);
            addDiagram4L1and4L2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp);
            addDiagram4K1and4K2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp);
            addDiagram4L1and4L2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp);
            addDiagram4K3and4K4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp);
            addDiagram4L3and4L4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp);
            addDiagram4K3and4K4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp);
            addDiagram4L3and4L4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp);
            
            addDiagram5A(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5B(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5C(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5F(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
                  
         }
         
      }
      
      delete [] temp;
      delete [] temp2;
   
   }

}

void CheMPS2::Heff::fillHeffDiag(double * memHeffDiag, const Sobject * denS, TensorC **** Ctensors, TensorD **** Dtensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   
   /*for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
      const int ikappa = denS->gReorder(ikappaBIS);
      cout << "Size of block " << ikappa << " is " << denS->gKappa2index(ikappa+1) - denS->gKappa2index(ikappa) << endl;
   }*/
   
   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
   
      const int ikappa = denS->gReorder(ikappaBIS);

      for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeffDiag[cnt] = 0.0; }
      
      addDiagonal1C(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
      addDiagonal1D(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
      addDiagonal2d3all(ikappa, memHeffDiag, denS);
      if (nLower>0){ addDiagonalExcitations(ikappa, memHeffDiag, denS, nLower, VeffTilde); }
      
      if (!atLeft){
         addDiagonal1A(ikappa, memHeffDiag, denS, Xtensors[indexS-1]);
         addDiagonal2b3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][0]);
         addDiagonal2c3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][1]);
         addDiagonal2b3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][0]);
         addDiagonal2c3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][1]);
      }
      
      if (!atRight){
         addDiagonal1B(ikappa, memHeffDiag, denS, Xtensors[indexS+1]);
         addDiagonal2e3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][1]);
         addDiagonal2f3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][0]);
         addDiagonal2e3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][1]);
         addDiagonal2f3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][0]);
      }
      
      if ((!atLeft) && (!atRight)){
         addDiagonal2a3spin0(ikappa, memHeffDiag, denS, Ctensors, F0tensors);
         addDiagonal2a3spin1(ikappa, memHeffDiag, denS, Dtensors, F1tensors);
      }
      
   }
   
}

double CheMPS2::Heff::SolveDAVIDSON(Sobject * denS, TensorL *** Ltensors, TensorA **** Atensors, TensorB **** Btensors, TensorC **** Ctensors, TensorD **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   int inc1 = 1;

   int veclength     = denS->gKappa2index( denS->gNKappa() );
   const double RTOL = CheMPS2::HEFF_DAVIDSON_RTOL_BASE * sqrt( 1.0 * veclength );

   Davidson deBoskabouter( veclength, CheMPS2::HEFF_DAVIDSON_NUM_VEC, CheMPS2::HEFF_DAVIDSON_NUM_VEC_KEEP, RTOL, CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF, CheMPS2::HEFF_debugPrint );
   double ** whichpointers = new double*[2];

   char instruction = deBoskabouter.FetchInstruction( whichpointers );
   assert( instruction == 'A' );
   denS->prog2symm(); // Convert mem of Sobject to symmetric conventions
   dcopy_(&veclength, denS->gStorage(), &inc1, whichpointers[0], &inc1); // Starting vector for Davidson is the current state of the Sobject in symmetric conventions
   fillHeffDiag(whichpointers[1], denS, Ctensors, Dtensors, F0tensors, F1tensors, Xtensors, nLower, VeffTilde);

   instruction = deBoskabouter.FetchInstruction( whichpointers );
   while ( instruction == 'B' ){
      makeHeff(whichpointers[0], whichpointers[1], denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
      instruction = deBoskabouter.FetchInstruction( whichpointers );
   }

   assert( instruction == 'C' );
   dcopy_( &veclength, whichpointers[0], &inc1, denS->gStorage(), &inc1 ); // Copy the solution in symmetric conventions back
   denS->symm2prog(); // Convert mem of Sobject to program conventions
   const double eigenvalue = whichpointers[1][0];
   if (CheMPS2::HEFF_debugPrint){ cout << "   Stats: nIt(DAVIDSON) = " << deBoskabouter.GetNumMultiplications() << endl; }
   delete [] whichpointers;
   return eigenvalue;

}

int CheMPS2::Heff::phase(const int TwoTimesPower){

   return (((TwoTimesPower/2)%2)!=0)?-1:1;

}



