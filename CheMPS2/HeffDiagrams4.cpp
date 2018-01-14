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

#include <math.h>
#include <stdlib.h>

#include "Heff.h"
#include "Lapack.h"
#include "MPIchemps2.h"
#include "Wigner.h"

void CheMPS2::Heff::addDiagram4A1and4A2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atens) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   int ILdown = Irreps::directProd(IL,Atens->get_irrep());
   
   char trans = 'T';
   char notrans = 'N';
   double beta = 1.0; //add
   
   //4A1A.spin0
   if ((N1==0) && (N2==0)){
   
      int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         double * Ablock = Atens->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
         int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,ILdown);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A1B.spin0
   if ((N1==1) && (N2==0)){
   
      int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,2,1,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         double * Ablock = Atens->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
         int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,ILdown);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A1C.spin0
   if ((N1==0) && (N2==1)){
   
      int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,1,2,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         double * Ablock = Atens->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
         int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,ILdown);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A1D.spin0 and 4A2A.spin0
   if ((N1==1) && (N2==1) && (TwoJ==0)){
   
      int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,2,2,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = -1.0;
         double * Ablock = Atens->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
         int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,ILdown);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
      memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,0,0,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         double * Ablock = Atens->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
         int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A2B.spin0
   if ((N1==2) && (N2==1)){
   
      int memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,1,0,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         double * Ablock = Atens->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
         int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A2C.spin0
   if ((N1==1) && (N2==2)){
   
      int memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,0,1,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         double * Ablock = Atens->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
         int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,ILdown);

         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A2D.spin0
   if ((N1==2) && (N2==2)){
   
      int memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = -1.0;
         double * Ablock = Atens->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
         int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,ILdown);

         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Ablock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
}

void CheMPS2::Heff::addDiagram4A1and4A2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Btens) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   int ILdown = Irreps::directProd(IL,Btens->get_irrep());
   
   char trans = 'T';
   char notrans = 'N';
   double beta = 1.0; //add
   
   //4A1A.spin1
   if ((N1==0) && (N2==0)){ //which means TwoSL==TwoSR --> no checker for TwoSLdown
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if (TwoSLdown>=0){
   
            int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,1,1,2,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               double factor = 1.0;
               double * Bblock = Btens->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
               int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSLdown,ILdown);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
            }
         }
      }
   }
   
   //4A1B.spin1
   if ((N1==1) && (N2==0)){
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,2,1,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSLdown + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(1,1,2,TwoSL,TwoSLdown,TwoSR);
               double * Bblock = Btens->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
               int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSLdown,ILdown);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A1C.spin1
   if ((N1==0) && (N2==1)){
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,1,2,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSLdown + TwoSR + 3);
               double factor = fase * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(1,1,2,TwoSL,TwoSLdown,TwoSR);
               double * Bblock = Btens->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
               int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSLdown,ILdown);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A1D.spin1 and 4A2A.spin1
   if ((N1==1) && (N2==1) && (TwoJ==2)){
   
      int TwoSLdown = TwoSR;
   
      int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,2,2,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         int fase = phase(TwoSLdown - TwoSL);
         double factor = fase * sqrt((TwoSL+1.0)/(TwoSLdown+1.0));
         double * Bblock = Btens->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
         int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSLdown,ILdown);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
      memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,0,0,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
   
         double factor = 1.0;
         double * Bblock = Btens->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
         int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSLdown,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A2B.spin1
   if ((N1==2) && (N2==1)){
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,1,0,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSL + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSLdown+1)) * Wigner::wigner6j(1,1,2,TwoSL,TwoSLdown,TwoSR);
               double * Bblock = Btens->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
               int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSLdown,ILdown);
         
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A2C.spin1
   if ((N1==1) && (N2==2)){
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,0,1,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSL + TwoSR + 3);
               double factor = fase * sqrt(3.0*(TwoSLdown+1)) * Wigner::wigner6j(1,1,2,TwoSL,TwoSLdown,TwoSR);
               double * Bblock = Btens->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
               int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSLdown,ILdown);
         
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A2D.spin1
   if ((N1==2) && (N2==2)){ //TwoSL==TwoSR --> no extra check for TwoSLdown needed
   
      for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
      
         if (TwoSLdown>=0){
   
            int memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,1,1,2,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSLdown - TwoSL);
               double factor = fase * sqrt((TwoSLdown+1.0)/(TwoSL+1.0));
               double * Bblock = Btens->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
               int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSLdown,ILdown);
         
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,Bblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4A3and4A4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctens) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   int ILdown = Irreps::directProd(IL,Ctens->get_irrep());
   
   char trans = 'T';
   char notrans = 'N';
   double beta = 1.0; //add
   
   //4A3A.spin0
   if ((N1==1) && (N2==0)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,0,1,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = sqrt(0.5);
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);

         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }  
   }
   
   //4A3B.spin0
   if ((N1==1) && (N2==1) && (TwoJ==0)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,0,2,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A3C.spin0
   if ((N1==2) && (N2==0)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A3D.spin0
   if ((N1==2) && (N2==1)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,1,2,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A4A.spin0
   if ((N1==0) && (N2==1)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,1,0,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = sqrt(0.5);
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A4B.spin0
   if ((N1==0) && (N2==2)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A4C.spin0
   if ((N1==1) && (N2==1) && (TwoJ==0)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,2,0,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = 1.0;
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A4D.spin0
   if ((N1==1) && (N2==2)){
   
      int memSkappa = denS->gKappa(NL,TwoSL,ILdown,2,1,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         double factor = - sqrt(0.5);
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSL,ILdown);
         double * ptr = Ctens->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
         
         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }

}

void CheMPS2::Heff::addDiagram4A3and4A4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtens) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   int ILdown = Irreps::directProd(IL,Dtens->get_irrep());
   
   char trans = 'T';
   char notrans = 'N';
   double beta = 1.0; //add
   
   //4A3A.spin1
   if ((N1==1) && (N2==0)){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,0,1,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSLdown + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(1, 1, 2, TwoSL, TwoSLdown, TwoSR);
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A3B.spin1
   if ((N1==1) && (N2==1) && (TwoJ==2)){
   
      int TwoSLdown = TwoSR;
   
      int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,0,2,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
      
         int fase = phase(TwoSLdown - TwoSL);
         double factor = fase * sqrt((TwoSL+1.0)/(TwoSLdown+1.0));
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
         double * ptr = Dtens->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
         
         dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
      }
   }
   
   //4A3C.spin1
   if ((N1==2) && (N2==0)){ //TwoSL==TwoSR --> no extra bounds needed on TwoSLdown
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if (TwoSLdown>=0){
   
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,1,1,2,NR,TwoSR,IR);
            if (memSkappa!=-1){
         
               double factor = -1.0;
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A3D.spin1
   if ((N1==2) && (N2==1)){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,1,2,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSLdown + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(1, 1, 2, TwoSL, TwoSLdown, TwoSR);
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
         
               dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A4A.spin1
   if ((N1==0) && (N2==1)){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
         
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,1,0,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSL + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSLdown+1)) * Wigner::wigner6j(1, 1, 2, TwoSL, TwoSLdown, TwoSR);
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
         
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A4B.spin1
   if ((N1==0) && (N2==2)){ //TwoSL==TwoSR --> no extra bounds needed on TwoSLdown
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if (TwoSLdown>=0){
   
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,1,1,2,NR,TwoSR,IR);
            if (memSkappa!=-1){
         
               int fase = phase(TwoSL - TwoSLdown);
               double factor = fase * sqrt((TwoSLdown+1.0)/(TwoSL+1.0));
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
            
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   //4A4C.spin1
   if ((N1==1) && (N2==1) && (TwoJ==2)){

      int TwoSLdown = TwoSR;

      int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,2,0,0,NR,TwoSR,IR);
      if (memSkappa!=-1){

         double factor = -1.0;
         int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
         double * ptr = Dtens->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);

         dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);

      }
   }
   
   //4A4D.spin1
   if ((N1==1) && (N2==2)){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=1) && (TwoSLdown>=0)){
   
            int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,2,1,1,NR,TwoSR,IR);
            if (memSkappa!=-1){
      
               int fase = phase(TwoSL + TwoSR + 1);
               double factor = fase * sqrt(3.0*(TwoSLdown+1)) * Wigner::wigner6j(1, 1, 2, TwoSL, TwoSLdown, TwoSR);
               int dimLdown = denBK->gCurrentDim(theindex,NL,TwoSLdown,ILdown);
               double * ptr = Dtens->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
         
               dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,ptr,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4B1and4B2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Aleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4B1A.spin0
   if (N1==0){
      
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSR + TwoSL + 1 + TwoJdown + 2*TwoS2);
               const double factor = fase * sqrt(0.5*(TwoSR+1)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL-2, TwoSL, ILdown, 1, N2, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown, &beta,temp,&dimLdown);
                     
                        double * Ablock = Aleft[l_index-theindex][0]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                   
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B1B.spin0
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
         
            int fase = phase(TwoSR + TwoSL + 2 + TwoJ + 2*TwoS2);
            const double factor = fase * sqrt(0.5*(TwoSR+1)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, TwoS2, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex][0]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL-2, TwoSL, ILdown, 2, N2, TwoS2, NR-1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                 
                     int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                     double * Ablock = Aleft[l_index-theindex][0]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4B2A.spin0
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
         
            int fase = phase(TwoSRdown + TwoSL + 1 + TwoJ + 2*TwoS2);
            const double factor = fase * sqrt(0.5*(TwoSRdown+1)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, TwoS2, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex][0]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL+2, TwoSL, ILdown, 0, N2, TwoS2, NR+1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                 
                     int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                     double * Ablock = Aleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4B2B.spin0
   if (N1==2){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSRdown + TwoSL + 2 + TwoJdown + 2*TwoS2);
               const double factor = fase * sqrt(0.5*(TwoSRdown+1)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL+2, TwoSL, ILdown, 1, N2, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * Ablock = Aleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4B1and4B2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Bleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4B1A.spin1
   if (N1==0){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = (TwoS2==0)?1:-1;
                  const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoSL+1)*(TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS2);
   
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex][0]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL-2, TwoSLdown, ILdown, 1, N2, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown, Lblock,&dimRdown, &beta,temp, &dimLdown);
                     
                           double * Bblock = Bleft[l_index-theindex][0]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B1B.spin1
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSR - TwoSRdown + TwoSL + 3 - TwoSLdown + 2*TwoS2);
               const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoSL+1)*(TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS2);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL-2, TwoSLdown, ILdown, 2, N2, TwoS2, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp, &dimLdown);
                     
                        double * Bblock = Bleft[l_index-theindex][0]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B2A.spin1
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = (TwoS2==0)?1:-1;
               const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoSLdown+1)*(TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS2);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL+2, TwoSLdown, ILdown, 0, N2, TwoS2, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * Bblock = Bleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B2B.spin1
   if (N1==2){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSLdown + 3 - TwoSL + TwoSRdown - TwoSR + 2*TwoS2);
                  const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoJdown+1)*(TwoSLdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS2);
   
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( theindex, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex][0]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL+2, TwoSLdown, ILdown, 1, N2, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                 
                           int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp, &dimLdown);
                     
                           double * Bblock = Bleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4B3and4B4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Cleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4B3A.spin0
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
         
            int fase = phase(TwoSR + TwoSL + TwoJ + 2*TwoS2);
            const double factor = fase * sqrt(0.5*(TwoSR+1)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, TwoS2, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex][0]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL, TwoSL, ILdown, 0, N2, TwoS2, NR-1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                  
                     int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                     double * ptr = Cleft[l_index-theindex][0]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
                  
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4B3B.spin0
   if (N1==2){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSR + TwoSL + 1 + TwoJdown + 2*TwoS2);
               const double factor = fase * sqrt(0.5*(TwoSR+1)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSL, ILdown, 1, N2, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp, &dimLdown);
                     
                        double * ptr = Cleft[l_index-theindex][0]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
                  
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B4A.spin0
   if (N1==0){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSRdown + TwoSL + TwoJdown + 2*TwoS2);
               const double factor = fase * sqrt(0.5*(TwoSRdown+1)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSL, ILdown, 1, N2, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * ptr = Cleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
                  
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B4B.spin0
   if (N1==1){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
         
            int fase = phase(TwoSRdown + TwoSL + 1 + TwoJ + 2*TwoS2);
            const double factor = fase * sqrt(0.5*(TwoSRdown+1)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, TwoS2, 1, TwoSRdown, TwoSR, TwoSL);
    
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex][0]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL, TwoSL, ILdown, 2, N2, TwoS2, NR+1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                  
                     int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                     double * ptr = Cleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
                  
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4B3and4B4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Dleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4B3A.spin1
   if (N1==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSL-TwoSLdown + TwoSR - TwoSRdown + 3 + 2*TwoS2);
               const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoJ+1)*(TwoSL+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS2);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, 0, N2, TwoS2, NR-1, TwoSRdown, IRdown);
                
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp, &dimLdown);
                     
                        double * ptr = Dleft[l_index-theindex][0]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
                        
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B3B.spin1
   if (N1==2){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = (TwoS2==0)?-1:1;
                  const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoJdown+1)*(TwoSL+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS2);
    
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex][0]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, 1, N2, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp, &dimLdown);
                     
                           double * ptr = Dleft[l_index-theindex][0]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
                  
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B4A.spin1
   if (N1==0){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSRdown - TwoSR + TwoSLdown - TwoSL + 3 + 2*TwoS2);
                  const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoJdown+1)*(TwoSLdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS2);
    
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex][0]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, 1, N2, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp, &dimLdown);
                     
                           double * ptr = Dleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
                  
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4B4B.spin1
   if (N1==1){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = (TwoS2==0)?-1:1;
               const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoJ+1)*(TwoSLdown+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS2);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex][0]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, 2, N2, TwoS2, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * ptr = Dleft[l_index-theindex][0]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
                  
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4C1and4C2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Aleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4C1A.spin0
   if (N2==0){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               const double factor = phase(TwoSR + TwoSL + 2 + TwoS1) * sqrt(0.5*(TwoSR+1)*(TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown, TwoS1, 1, TwoSR, TwoSRdown, TwoSL);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL-2, TwoSL, ILdown, N1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                        double * Ablock = Aleft[l_index-theindex-1][1]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C1B.spin0
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
         
            const double factor = phase(TwoSR + TwoSL + 3 + TwoS1) * sqrt(0.5*(TwoSR+1)*(TwoJ+1))
                                * Wigner::wigner6j(TwoJ, TwoS1, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex-1][1]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL-2, TwoSL, ILdown, N1, 2, TwoS1, NR-1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                 
                     int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                     double * Ablock = Aleft[l_index-theindex-1][1]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4C2A.spin0
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
         
            const double factor = phase(TwoSRdown + TwoSL + 2 + TwoS1) * sqrt(0.5*(TwoSRdown+1)*(TwoJ+1))
                                * Wigner::wigner6j(TwoJ, TwoS1, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex-1][1]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL+2, TwoSL, ILdown, N1, 0, TwoS1, NR+1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                 
                     int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                     double * Ablock = Aleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4C2B.spin0
   if (N2==2){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               const double factor = phase(TwoSRdown + TwoSL + 3 + TwoS1) * sqrt(0.5*(TwoSRdown+1)*(TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown, TwoS1, 1, TwoSR, TwoSRdown, TwoSL);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Aleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL+2, TwoSL, ILdown, N1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * Ablock = Aleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Ablock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4C1and4C2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Bleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4C1A.spin1
   if (N2==0){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  const double factor = phase(1 + TwoS1 - TwoJdown)
                                      * sqrt(3.0*(TwoSR+1)*(TwoSL+1)*(TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS1);
   
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex-1][1]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL-2, TwoSLdown, ILdown, N1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                           double * Bblock = Bleft[l_index-theindex-1][1]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C1B.spin1
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               const double factor = phase(TwoSR - TwoSRdown + TwoSL - TwoSLdown + TwoS1 - TwoJ)
                                   * sqrt(3.0*(TwoSR+1)*(TwoSL+1)*(TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS1);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL-2, TwoSLdown, ILdown, N1, 2, TwoS1, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL-2, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                        double * Bblock = Bleft[l_index-theindex-1][1]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C2A.spin1
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               const double factor = phase(1 + TwoS1 - TwoJ)
                                   * sqrt(3.0*(TwoSRdown+1)*(TwoSLdown+1)*(TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS1);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL+2, TwoSLdown, ILdown, N1, 0, TwoS1, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                 
                        int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * Bblock = Bleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C2B.spin1
   if (N2==2){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  const double factor = phase(TwoSLdown - TwoSL + TwoSRdown - TwoSR + TwoS1 - TwoJdown)
                                      * sqrt(3.0*(TwoSRdown+1)*(TwoJdown+1)*(TwoSLdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS1);
   
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( theindex+1, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Bleft[l_index-theindex-1][1]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL+2, TwoSLdown, ILdown, N1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                 
                           int dimLdown = denBK->gCurrentDim(theindex  , NL+2, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                           double * Bblock = Bleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,Bblock,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4C3and4C4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Cleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4C3A.spin0
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
         
            const double factor = phase(TwoSR + TwoSL + 1 + TwoS1) * sqrt(0.5*(TwoSR+1)*(TwoJ+1))
                                * Wigner::wigner6j(TwoJ, TwoS1, 1, TwoSRdown, TwoSR, TwoSL);
   
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex-1][1]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL, TwoSL, ILdown, N1, 0, TwoS1, NR-1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                  
                     int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                     double * ptr = Cleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
                  
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4C3B.spin0
   if (N2==2){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               const double factor = phase(TwoSR + TwoSL + 2 + TwoS1) * sqrt(0.5*(TwoSR+1)*(TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown, TwoS1, 1, TwoSR, TwoSRdown, TwoSL);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSL, ILdown, N1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                        double * ptr = Cleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
                  
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C4A.spin0
   if (N2==0){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSL-TwoSRdown)<=TwoJdown) && (TwoSRdown>=0)){
         
               const double factor = phase(TwoSRdown + TwoSL + 1 + TwoS1) * sqrt(0.5*(TwoSRdown+1)*(TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown, TwoS1, 1, TwoSR, TwoSRdown, TwoSL);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSL, ILdown, N1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * ptr = Cleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
                  
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C4B.spin0
   if (N2==1){

      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
         
            const double factor = phase(TwoSRdown+TwoSL+2+TwoS1) * sqrt(0.5*(TwoSRdown+1)*(TwoJ+1))
                                * Wigner::wigner6j(TwoJ, TwoS1, 1, TwoSRdown, TwoSR, TwoSL);
    
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, Cleft[l_index-theindex-1][1]->get_irrep() );
                  int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                  int memSkappa = denS->gKappa(NL, TwoSL, ILdown, N1, 2, TwoS1, NR+1, TwoSRdown, IRdown);
                  
                  if (memSkappa!=-1){
                  
                     int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSL,     ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                     double * ptr = Cleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
                  
                     alpha = factor;
                     beta = 1.0; //add
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4C3and4C4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator *** Dleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4C3A.spin1
   if (N2==1){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSL-TwoSLdown + TwoSR - TwoSRdown + TwoS1 - TwoJ);
               const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoJ+1)*(TwoSL+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS1);
   
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, N1, 0, TwoS1, NR-1, TwoSRdown, IRdown);
                
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                        double * ptr = Dleft[l_index-theindex-1][1]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
                        
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C3B.spin1
   if (N2==2){
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(3 + TwoS1 - TwoJdown);
                  const double factor = fase * sqrt(3.0*(TwoSR+1)*(TwoJdown+1)*(TwoSL+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS1);
    
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex-1][1]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, N1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRdown,&beta,temp,&dimLdown);
                     
                           double * ptr = Dleft[l_index-theindex-1][1]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
                           
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C4A.spin1
   if (N2==0){
   
      for (int TwoSRdown = TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSRdown - TwoSR + TwoSLdown - TwoSL + TwoS1 - TwoJdown);
                  const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoJdown+1)*(TwoSLdown+1))
                                      * Wigner::wigner9j(2, TwoSLdown, TwoSL, 1, TwoSRdown, TwoSR, 1, TwoJdown, TwoS1);
    
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex-1][1]->get_irrep() );
                        int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                        int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, N1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                  
                        if (memSkappa!=-1){
                  
                           int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                           double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                           double * ptr = Dleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
                           
                           alpha = factor;
                           beta = 1.0; //add
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4C4B.spin1
   if (N2==1){

      for (int TwoSRdown = TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown = TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(3 + TwoS1 - TwoJ);
               const double factor = fase * sqrt(3.0*(TwoSRdown+1)*(TwoJ+1)*(TwoSLdown+1))
                                   * Wigner::wigner9j(2, TwoSL, TwoSLdown, 1, TwoSR, TwoSRdown, 1, TwoJ, TwoS1);
    
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), theindex+1, l_index ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, Dleft[l_index-theindex-1][1]->get_irrep() );
                     int IRdown = Irreps::directProd(IR, denBK->gIrrep(l_index) );
                     int memSkappa = denS->gKappa(NL, TwoSLdown, ILdown, N1, 2, TwoS1, NR+1, TwoSRdown, IRdown);
                  
                     if (memSkappa!=-1){
                  
                        int dimLdown = denBK->gCurrentDim(theindex  , NL,   TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                  
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        double alpha = 1.0;
                        double beta = 0.0; //set
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,Lblock,&dimRup,&beta,temp,&dimLdown);
                     
                        double * ptr = Dleft[l_index-theindex-1][1]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
                        
                        alpha = factor;
                        beta = 1.0; //add
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,ptr,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int inc = 1;
   double beta = 1.0; //add
   int ILdown = Irreps::directProd(IL, denBK->gIrrep(theindex+1));
   
   //4D1A and 4D1B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4D1AB ) == MPIRANK ) && (N1==0) && (N2>0)){
   #else
   if ((N1==0) && (N2>0)){
   #endif
   
      int TwoS2down = (N2==1)?0:1;
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if (abs(TwoSLdown-TwoSR)<=TwoS2down){
      
            int dimLdown = denBK->gCurrentDim(theindex,NL-1,TwoSLdown,ILdown);
            if (dimLdown>0){
            
               int size = dimLup * dimLdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
         
               int number = 0;
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(l_index,theindex+1,theindex,theindex);
                     double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
         
               if (number>0){
            
                  double factor = -1.0;
                  if (N2==1){
                     int fase = phase(TwoSR+1-TwoSL);
                     factor = fase * sqrt((TwoSL+1.0)/(TwoSR+1.0));
                  }
                  int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2-1,TwoS2down,NR,TwoSR,IR);
                  dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
               }
            }
         }
      }
   }
   
   //4D2A and 4D2B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4D2AB ) == MPIRANK ) && (N1==2) && (N2<2)){
   #else
   if ((N1==2) && (N2<2)){
   #endif
      
      int TwoS2down = (N2==0)?1:0;
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if (abs(TwoSLdown-TwoSR)<=TwoS2down){
      
            int dimLdown = denBK->gCurrentDim(theindex,NL+1,TwoSLdown,ILdown);
            if (dimLdown>0){
            
               int size = dimLup * dimLdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
         
               int number = 0;
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(l_index,theindex+1,theindex,theindex);
                     double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
         
               if (number>0){
            
                  double factor = -1.0;
                  if (N2==0){
                     int fase = phase(TwoSR+1-TwoSLdown);
                     factor = fase * sqrt((TwoSLdown+1.0)/(TwoSR+1.0));
                  }
                  int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,0,N2+1,TwoS2down,NR,TwoSR,IR);
                  dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
               }
            }
         }
      }
   }
   
   //4D3A and 4D3B and 4D3C and 4D3D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4D3ABCD ) == MPIRANK ) && (N1>0) && (N2<2)){
   #else
   if ((N1>0) && (N2<2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex,NL-1,TwoSLdown,ILdown);
         if (dimLdown>0){
         
            int TwoSdownSum = ((N1==1)?1:0) + ((N2==0)?1:0);
            int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart ; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int size = dimLup * dimLdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
                  double alpha_fact = 0.0;
                  if ((N1==1) && (N2==0)){ //4D3A
                     int fase = phase(TwoSLdown + TwoSR + 2);
                     alpha_fact = fase * sqrt((TwoSL+1.0)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown,1,1,TwoSL,TwoSLdown,TwoSR);
                  }
                  if ((N1==1) && (N2==1)){ //4D3B
                     int fase = phase(TwoSLdown + TwoSR + 3 + TwoJ);
                     alpha_fact = fase * sqrt((TwoSL+1.0)*(TwoJ+1)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,1,1,TwoSLdown);
                  }
                  if ((N1==2) && (N2==0)){ //4D3C
                     alpha_fact = -1.0;
                  }
                  if ((N1==2) && (N2==1)){ //4D3D
                     int fase = phase(TwoSL+1-TwoSR);
                     alpha_fact = fase * sqrt((TwoSL+1.0)/(TwoSR+1.0));
                  }
         
                  int number = 0;
                  for (int l_index=0; l_index<theindex; l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                  
                        double alpha = 0.0;
                        if ((N1==1) && (N2==0)){ //4D3A
                           alpha = alpha_fact * ( Prob->gMxElement(l_index, theindex, theindex, theindex+1) + ((TwoJdown==0)?1:-1) * Prob->gMxElement(l_index, theindex, theindex+1, theindex) );
                        }
                        if ((N1==1) && (N2==1)){ //4D3B
                           alpha = alpha_fact * Prob->gMxElement(l_index,theindex,theindex+1,theindex);
                           if (TwoJ==0){
                              alpha += sqrt(2.0) * Prob->gMxElement(l_index,theindex,theindex,theindex+1);
                           }
                        }
                        if (N1==2){ //4D3C and 4D3D
                           alpha = alpha_fact * ( Prob->gMxElement(l_index, theindex, theindex, theindex+1) - 2 * Prob->gMxElement(l_index, theindex, theindex+1, theindex) );
                        }
                     
                        double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
         
                  if (number>0){
            
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,N2+1,TwoJdown,NR,TwoSR,IR);
                     dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   //4D4A and 4D4B and 4D4C and 4D4D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4D4ABCD ) == MPIRANK ) && (N1>0) && (N2>0)){
   #else
   if ((N1>0) && (N2>0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex,NL+1,TwoSLdown,ILdown);
         if (dimLdown>0){
         
            //int N1down = N1;
            //int N2down = N2-1;
            int TwoSdownSum = ((N1==1)?1:0) + ((N2==2)?1:0);
            int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart ; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int size = dimLup * dimLdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
                  double alpha_fact = 0.0;
                  if ((N1==1) && (N2==1)){ //4D4A
                     int fase = phase(TwoSL + TwoSR + 2);
                     alpha_fact = fase * sqrt((TwoSLdown+1.0)*(TwoJ+1)) * Wigner::wigner6j(TwoJ,1,1,TwoSLdown,TwoSL,TwoSR);
                  }
                  if ((N1==1) && (N2==2)){ //4D4B
                     int fase = phase(TwoSL + TwoSR + 3 + TwoJdown);
                     alpha_fact = fase * sqrt((TwoSLdown+1.0)*(TwoJdown+1)) * Wigner::wigner6j(TwoSLdown,TwoSR,TwoJdown,1,1,TwoSL);
                  }
                  if ((N1==2) && (N2==1)){ //4D4C
                     alpha_fact = -1.0;
                  }
                  if ((N1==2) && (N2==2)){ //4D4D
                     int fase = phase(TwoSLdown+1-TwoSR);
                     alpha_fact = fase * sqrt((TwoSLdown+1.0)/(TwoSR+1.0));
                  }
         
                  int number = 0;
                  for (int l_index=0; l_index<theindex; l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                  
                        double alpha = 0.0;
                        if ((N1==1) && (N2==1)){ //4D4A
                           alpha = alpha_fact * ( Prob->gMxElement(l_index, theindex, theindex, theindex+1) + ((TwoJ==0)?1:-1) * Prob->gMxElement(l_index, theindex, theindex+1, theindex) );
                        }
                        if ((N1==1) && (N2==2)){ //4D4B
                           alpha = alpha_fact * Prob->gMxElement(l_index,theindex,theindex+1,theindex);
                           if (TwoJdown==0){
                              alpha += sqrt(2.0) * Prob->gMxElement(l_index,theindex,theindex,theindex+1);
                           }
                        }
                        if (N1==2){ //4D4C and 4D4D
                           alpha = alpha_fact * ( Prob->gMxElement(l_index, theindex, theindex, theindex+1) - 2 * Prob->gMxElement(l_index, theindex, theindex+1, theindex) );
                        }
                     
                        double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
         
                  if (number>0){
            
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,N2-1,TwoJdown,NR,TwoSR,IR);
                     dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4E(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int inc = 1;
   
   //4E1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E1 ) == MPIRANK ) && (N1==0)){
   #else
   if (N1==0){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSL+TwoSR-TwoS2);
               const double factor = fase * sqrt((TwoSL+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoS2,TwoSRdown,TwoSLdown,1);
         
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
            
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                     
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){
                     
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (Irrep == denBK->gIrrep(l_beta)){
                                    double * LblockRight = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    double prefact = Prob->gMxElement(l_alpha,l_beta,theindex,theindex);
                                    daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                 }
                              }
                           
                              int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2,TwoS2,NR+1,TwoSRdown,IRdown);
                              double alpha = factor;
                              double beta = 0.0; //set
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                           
                              alpha = 1.0;
                              beta = 1.0; //add
                              double * LblockLeft = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                              dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4E2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E2 ) == MPIRANK ) && (N1==2)){
   #else
   if (N1==2){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSLdown+TwoSRdown-TwoS2);
               const double factor = fase * sqrt((TwoSLdown+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoS2,TwoSR,TwoSL,1);
         
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
            
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                        if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                        if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                     
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){
                     
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                 if (Irrep == denBK->gIrrep(l_delta)){
                                    double * LblockRight = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    double prefact = Prob->gMxElement(l_gamma,l_delta,theindex,theindex);
                                    daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                 }
                              }
                           
                              int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,0,N2,TwoS2,NR-1,TwoSRdown,IRdown);
                              double alpha = factor;
                              double beta = 0.0; //set
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                           
                              alpha = 1.0;
                              beta = 1.0; //add
                              double * LblockLeft = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4E3A
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E3A ) == MPIRANK ) && (N1==1)){
   #else
   if (N1==1){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSL + TwoSR + TwoJ + TwoSLdown + TwoSRdown + 1 - TwoS2);
                  const double factor1 = fase * sqrt((TwoJ+1)*(TwoJdown+1)*(TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSL, TwoSRdown, TwoS2, TwoJdown, 1, TwoSLdown) * Wigner::wigner6j(TwoJ, 1, TwoS2, TwoSRdown, TwoSL, TwoSR);
               
                  double factor2 = 0.0;
                  if (TwoJ == TwoJdown){
                     fase = phase(TwoSL+TwoSRdown+TwoJ+3+2*TwoS2);
                     factor2 = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSRdown, TwoJ, TwoSR, TwoSL, 1);
                  }
            
                  for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
            
                     int ILdown = Irreps::directProd(IL,Irrep);
                     int IRdown = Irreps::directProd(IR,Irrep);
                     int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
               
                     if ((dimLdown>0) && (dimRdown>0)){
                        bool isPossibleLeft = false;
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                        }
                        bool isPossibleRight = false;
                        for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                           if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                        }
                        if ( (isPossibleLeft) && (isPossibleRight) ){
                     
                           for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                              if (Irrep == denBK->gIrrep(l_alpha)){
                     
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                                 for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                    if (Irrep == denBK->gIrrep(l_delta)){
                                       double * LblockRight = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                       double prefact = factor1 * Prob->gMxElement(l_alpha,theindex,theindex,l_delta);
                                       if (TwoJ == TwoJdown){ prefact += factor2 * Prob->gMxElement(l_alpha,theindex,l_delta,theindex); }
                                       daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                    }
                                 }
                           
                                 int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR-1,TwoSRdown,IRdown);
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockLeft = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                                 dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4E3B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E3B ) == MPIRANK ) && (N1==2)){
   #else
   if (N1==2){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSL + TwoSRdown - TwoS2 + 3);
               const double factor = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSRdown, TwoS2, TwoSR, TwoSL, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                        if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                     
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){
                     
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                 if (Irrep == denBK->gIrrep(l_delta)){
                                    double * LblockRight = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    double prefact = Prob->gMxElement(l_alpha,theindex,theindex,l_delta) - 2 * Prob->gMxElement(l_alpha,theindex,l_delta,theindex);
                                    daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                 }
                              }
                             
                              int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2,TwoS2,NR-1,TwoSRdown,IRdown);
                              double alpha = factor;
                              double beta = 0.0; //set
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              alpha = 1.0;
                              beta = 1.0; //add
                              double * LblockLeft = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                              dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4E4A
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E4A ) == MPIRANK ) && (N1==1)){
   #else
   if (N1==1){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSL + TwoSR + TwoJdown + TwoSLdown + TwoSRdown + 1 - TwoS2);
                  const double factor1 = fase * sqrt((TwoJ+1)*(TwoJdown+1)*(TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSR, TwoS2, TwoJ, 1, TwoSL) * Wigner::wigner6j(TwoJdown, 1, TwoS2, TwoSR, TwoSLdown, TwoSRdown);
                  
                  double factor2 = 0.0;
                  if (TwoJ == TwoJdown){
                     fase = phase(TwoSLdown+TwoSR+TwoJ+3+2*TwoS2);
                     factor2 = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL, TwoSR, TwoJ, TwoSRdown, TwoSLdown, 1);
                  }
                  
                  for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                     int ILdown = Irreps::directProd(IL,Irrep);
                     int IRdown = Irreps::directProd(IR,Irrep);
                     int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                        bool isPossibleLeft = false;
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                        }
                        bool isPossibleRight = false;
                        for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                           if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                        }
                        if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                           for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                              if (Irrep == denBK->gIrrep(l_gamma)){
                        
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (Irrep == denBK->gIrrep(l_beta)){
                                       double * LblockRight = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                       double prefact = factor1 * Prob->gMxElement(l_gamma,theindex,theindex,l_beta);
                                       if (TwoJ == TwoJdown){ prefact += factor2 * Prob->gMxElement(l_gamma,theindex,l_beta,theindex); }
                                       daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                    }
                                 }
                              
                                 int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR+1,TwoSRdown,IRdown);
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                                 
                                 beta = 1.0; //add
                                 double * LblockLeft = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                                 dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4E4B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4E4B ) == MPIRANK ) && (N1==2)){
   #else
   if (N1==2){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSLdown + TwoSR - TwoS2 + 3);
               const double factor = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL, TwoSR, TwoS2, TwoSRdown, TwoSLdown, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                        if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                     
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){
                     
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (Irrep == denBK->gIrrep(l_beta)){
                                    double * LblockRight = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    double prefact = Prob->gMxElement(l_gamma,theindex,theindex,l_beta) - 2 * Prob->gMxElement(l_gamma,theindex,l_beta,theindex);
                                    daxpy_(&size,&prefact,LblockRight,&inc,temp,&inc);
                                 }
                              }
                             
                              int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,2,N2,TwoS2,NR+1,TwoSRdown,IRdown);
                              double alpha = factor;
                              double beta = 0.0; //set
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              alpha = 1.0;
                              beta = 1.0; //add
                              double * LblockLeft = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockLeft,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int inc = 1;
   double beta = 1.0; //add
   int IRdown = Irreps::directProd(IR, denBK->gIrrep(theindex+1)); //I_{L} must be equal to I_{i+1}
   
   //4F1A and 4F1B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4F1AB ) == MPIRANK ) && (N1==2) && (N2<2)){
   #else
   if ((N1==2) && (N2<2)){
   #endif
   
      int TwoS2down = (N2==1)?0:1;
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=TwoS2down){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
            if (dimRdown>0){
               
               int size = dimRup * dimRdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
            
               int number = 0;
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(theindex,theindex,theindex+1,l_index);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
            
               if (number>0){
               
                  double factor = 0.0;
                  if (N2==1){
                     factor = - sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
                  } else {
                     factor = phase(TwoSR+1-TwoSRdown);
                  }
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,0,N2+1,TwoS2down,NR-1,TwoSRdown,IRdown);
                  dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
               }
            }
         }
      }
   }
   
   //4F2A and 4F2B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4F2AB ) == MPIRANK ) && (N1==0) && (N2>0)){
   #else
   if ((N1==0) && (N2>0)){
   #endif
   
      int TwoS2down = (N2==1)?0:1;
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=TwoS2down){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
            if (dimRdown>0){
               
               int size = dimRup * dimRdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
            
               int number = 0;
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(theindex,theindex,theindex+1,l_index);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
            
               if (number>0){
               
                  double factor = 0.0;
                  if (N2==2){
                     factor = - sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
                  } else {
                     factor = phase(TwoSRdown+1-TwoSR);
                  }
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,2,N2-1,TwoS2down,NR+1,TwoSRdown,IRdown);
                  dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
               }
            }
         }
      }
   }
   
   //4F3A and 4F3B and 4F3C and 4F3D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4F3ABCD ) == MPIRANK ) && (N1>0) && (N2>0)){
   #else
   if ((N1>0) && (N2>0)){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
         if (dimRdown>0){
         
            //int N1down = N1;
            //int N2down = N2-1;
            int TwoSdownSum = ((N1==1)?1:0) + ((N2==2)?1:0);
            int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown=TwoJstart ; TwoJdown<=TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
               
                  int size = dimRup * dimRdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                  
                  double factor = 0.0;
                  double factor2 = 0.0;
                  if ((N1==1) && (N2==1)){ // 4F3A
                     int fase = phase(TwoSL+TwoSR+2);
                     factor = fase * sqrt((TwoSR+1.0)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, 1, 1, TwoSRdown, TwoSR, TwoSL);
                  }
                  if ((N1==1) && (N2==2)){ // 4F3B
                     int fase = phase(TwoSL+TwoSR+3);
                     factor = fase * sqrt((TwoSR+1.0)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, 1, 1, TwoSR, TwoSRdown, TwoSL);
                     factor2 = (TwoJdown==0) ? sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)) : 0.0;
                  }
                  if ((N1==2) && (N2==1)){ // 4F3C
                     factor = sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
                  }
                  if ((N1==2) && (N2==2)){ // 4F3D
                     factor = phase(TwoSR + 1 - TwoSRdown);
                  }
            
                  int number = 0;
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                        double prefact = 0.0;
                        if ((N1==1) && (N2==1)){ // 4F3A
                           prefact = factor * ( Prob->gMxElement(theindex,theindex+1,theindex,l_index)
                           + ((TwoJ==0)?1:-1) * Prob->gMxElement(theindex,theindex+1,l_index,theindex) );
                        }
                        if ((N1==1) && (N2==2)){ // 4F3B
                                             prefact  = factor  * Prob->gMxElement(theindex,theindex+1,theindex,l_index);
                           if (TwoJdown==0){ prefact += factor2 * Prob->gMxElement(theindex,theindex+1,l_index,theindex); }
                        }
                        if (N1==2){ // 4F3C and 4F3D
                           prefact = factor * ( 2 * Prob->gMxElement(theindex,theindex+1,theindex,l_index)
                                                  - Prob->gMxElement(theindex,theindex+1,l_index,theindex) );
                        }
                        
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        daxpy_(&size,&prefact,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
            
                  if (number>0){
               
                     double alpha = 1.0;
                     int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,N2-1,TwoJdown,NR-1,TwoSRdown,IRdown);
                     dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            
                  }
               }
            }
         }
      }
   }
   
   //4F4A and 4F4B and 4F4C and 4F4D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4F4ABCD ) == MPIRANK ) && (N1>0) && (N2<2)){
   #else
   if ((N1>0) && (N2<2)){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
         if (dimRdown>0){
         
            //int N1down = N1;
            //int N2down = N2+1;
            int TwoSdownSum = ((N1==1)?1:0) + ((N2==0)?1:0);
            int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
               
                  int size = dimRup * dimRdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                  
                  double factor = 0.0;
                  double factor2 = 0.0;
                  if ((N1==1) && (N2==0)){ // 4F3A
                     int fase = phase(TwoSL+TwoSRdown+2);
                     factor = fase * sqrt((TwoSRdown+1.0)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown, 1, 1, TwoSR, TwoSRdown, TwoSL);
                  }
                  if ((N1==1) && (N2==1)){ // 4F3B
                     int fase = phase(TwoSL+TwoSRdown+3);
                     factor = fase * sqrt((TwoSRdown+1.0)*(TwoJ+1)) * Wigner::wigner6j(TwoJ, 1, 1, TwoSRdown, TwoSR, TwoSL);
                     factor2 = (TwoJ==0) ? sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)) : 0.0;
                  }
                  if ((N1==2) && (N2==0)){ // 4F3C
                     factor = sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
                  }
                  if ((N1==2) && (N2==1)){ // 4F3D
                     factor = phase(TwoSRdown + 1 - TwoSR);
                  }
            
                  int number = 0;
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                        double prefact = 0.0;
                        if ((N1==1) && (N2==0)){ // 4F3A
                              prefact = factor * ( Prob->gMxElement(theindex,theindex+1,theindex,l_index)
                          + ((TwoJdown==0)?1:-1) * Prob->gMxElement(theindex,theindex+1,l_index,theindex) );
                        }
                        if ((N1==1) && (N2==1)){ // 4F3B
                                         prefact  = factor  * Prob->gMxElement(theindex,theindex+1,theindex,l_index);
                           if (TwoJ==0){ prefact += factor2 * Prob->gMxElement(theindex,theindex+1,l_index,theindex); }
                        }
                        if (N1==2){ // 4F3C and 4F3D
                           prefact = factor * ( 2 * Prob->gMxElement(theindex,theindex+1,theindex,l_index)
                                                  - Prob->gMxElement(theindex,theindex+1,l_index,theindex) );
                        }
                        
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        daxpy_(&size,&prefact,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
            
                  if (number>0){
               
                     double alpha = 1.0;
                     int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,N2+1,TwoJdown,NR+1,TwoSRdown,IRdown);
                     dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4G(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int inc = 1;
   double beta = 1.0; //add
   int IRdown = Irreps::directProd(IR, denBK->gIrrep(theindex)); //I_{L} must be equal to I_{i+1}
   
   //4G1A and 4G1B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4G1AB ) == MPIRANK ) && (N1<2) && (N2==2)){
   #else
   if ((N1<2) && (N2==2)){
   #endif
   
      int TwoS1down = (N1==1)?0:1;
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=TwoS1down){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
            if (dimRdown>0){
               
               int size = dimRup * dimRdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
            
               int number = 0;
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     double alpha = Prob->gMxElement(theindex,l_index,theindex+1,theindex+1);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
            
               if (number>0){
               
                  double factor = 0.0;
                  if (N1==1){ //4G1B
                     factor = - sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
                  } else {
                     factor = phase(TwoSR+1-TwoSRdown);
                  }
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,N1+1,0,TwoS1down,NR-1,TwoSRdown,IRdown);
                  dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
               }
            }
         }
      }
   }
   
   //4G2A and 4G2B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4G2AB ) == MPIRANK ) && (N1>0) && (N2==0)){
   #else
   if ((N1>0) && (N2==0)){
   #endif
   
      int TwoS1down = (N1==1)?0:1;
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=TwoS1down){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
            if (dimRdown>0){
               
               int size = dimRup * dimRdown;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
            
               int number = 0;
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     double alpha = Prob->gMxElement(theindex,l_index,theindex+1,theindex+1);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     number++;
                  }
               }
            
               if (number>0){

                  double factor = 0.0;
                  if (N1==2){ //4G2B --> bug fixed
                     factor = - sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
                  } else { //4G2A --> bug fixed
                     factor = phase(TwoSRdown+1-TwoSR);
                  }
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,N1-1,2,TwoS1down,NR+1,TwoSRdown,IRdown);
                  dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   }
   
   //4G3A and 4G3B and 4G3C and 4G3D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4G3ABCD ) == MPIRANK ) && (N1<2) && (N2>0)){
   #else
   if ((N1<2) && (N2>0)){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
         if (dimRdown>0){
            
            //int N1down = N1+1;
            //int N2down = N2;
            int TwoSdownSum = ((N1==1)?0:1) + ((N2==1)?1:0);
            int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoSdownSum<2))? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
               
                  int size = dimRup * dimRdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
                  double alpha_prefact = 0.0;
                  double alpha_prefact2 = 0.0;
                  if ((N1==0) && (N2==1)){
                     int fase = phase(TwoSL + TwoSRdown + 2);
                     alpha_prefact = fase * sqrt((TwoJdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoJdown, 1, 1, TwoSR, TwoSRdown, TwoSL);
                  }
                  if ((N1==1) && (N2==1)){
                     int fase = phase(TwoSL + TwoSRdown + TwoJ + 3);
                     alpha_prefact = fase * sqrt((TwoJ+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoJ, 1, 1, TwoSRdown, TwoSR, TwoSL);
                     alpha_prefact2 = (TwoJ==0)? sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)) : 0.0;
                  }
                  if ((N1==0) && (N2==2)){
                     alpha_prefact = - sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
                  }
                  if ((N1==1) && (N2==2)){
                     alpha_prefact = phase(TwoSR + 1 - TwoSRdown);
                  }
            
                  int number = 0;
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                        
                        double alpha = 0.0;
                        if ((N1==0) && (N2==1)){
                           alpha = alpha_prefact * ( Prob->gMxElement(theindex,theindex+1,theindex+1,l_index) + ((TwoJdown==0)?1:-1) * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1) );
                        }
                        if ((N1==1) && (N2==1)){
                           alpha = alpha_prefact * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1);
                           if (TwoJ==0){ alpha += alpha_prefact2 * Prob->gMxElement(theindex,theindex+1,theindex+1,l_index); }
                        }
                        if (N2==2){
                           alpha = alpha_prefact * (Prob->gMxElement(theindex,theindex+1,theindex+1,l_index) - 2 * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1) );
                        }
                        
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
            
                  if (number>0){
               
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL,TwoSL,IL,N1+1,N2,TwoJdown,NR+1,TwoSRdown,IRdown);
                     dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  }
               }
            }
         }
      }
   }
   
   //4G4A and 4G4B and 4G4C and 4G4D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4G4ABCD ) == MPIRANK ) && (N1>0) && (N2>0)){
   #else
   if ((N1>0) && (N2>0)){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
         if (dimRdown>0){
            
            //int N1down = N1-1;
            //int N2down = N2;
            int TwoSdownSum = ((N1==1)?0:1) + ((N2==1)?1:0);
            int TwoJstart = ((TwoSL!=TwoSRdown) || (TwoSdownSum < 2)) ? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart ; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
            
                  int size = dimRup * dimRdown;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
                  double alpha_prefact = 0.0;
                  double alpha_prefact2 = 0.0;
                  if ((N1==1) && (N2==1)){
                     int fase = phase(TwoSL + TwoSR + 2);
                     alpha_prefact = fase * sqrt((TwoJ+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoJ, 1, 1, TwoSRdown, TwoSR, TwoSL);
                  }
                  if ((N1==2) && (N2==1)){
                     int fase = phase(TwoSL + TwoSR + TwoJdown + 3);
                     alpha_prefact = fase * sqrt((TwoJdown+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoJdown, 1, 1, TwoSR, TwoSRdown, TwoSL);
                     alpha_prefact2 = (TwoJdown==0) ? sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)) : 0.0;
                  }
                  if ((N1==1) && (N2==2)){
                     alpha_prefact = - sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
                  }
                  if ((N1==2) && (N2==2)){
                     alpha_prefact = phase(TwoSRdown + 1 - TwoSR);
                  }
            
                  int number = 0;
                  for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                        
                        double alpha = 0.0;
                        if ((N1==1) && (N2==1)){
                           alpha = alpha_prefact * ( Prob->gMxElement(theindex,theindex+1,theindex+1,l_index) + ((TwoJ==0)?1:-1) * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1) );
                        }
                        if ((N1==2) && (N2==1)){
                           alpha = alpha_prefact * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1);
                           if (TwoJdown==0){ alpha += alpha_prefact2 * Prob->gMxElement(theindex,theindex+1,theindex+1,l_index); }
                        }
                        if (N2==2){
                           alpha = alpha_prefact * (Prob->gMxElement(theindex,theindex+1,theindex+1,l_index) - 2 * Prob->gMxElement(theindex,theindex+1,l_index,theindex+1) );
                        }
                        
                        double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                        number++;
                     }
                  }
            
                  if (number>0){
               
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL,TwoSL,IL,N1-1,N2,TwoJdown,NR-1,TwoSRdown,IRdown);
                     dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4H(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4H1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H1 ) == MPIRANK ) && (N2==2)){
   #else
   if (N2==2){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               int fase = phase(TwoSLdown + TwoSRdown - TwoS1);
               const double factor = fase * sqrt((TwoSLdown+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSL, TwoSR, TwoS1, TwoSRdown, TwoSLdown, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                        if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                        if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                 if (Irrep == denBK->gIrrep(l_delta)){
                                    double fact = factor * Prob->gMxElement(l_gamma,l_delta,theindex+1,theindex+1);
                                    double * LblockR = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    int inc = 1;
                                    daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              double * LblockL = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                              
                              int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,0,TwoS1,NR-1,TwoSRdown,IRdown);
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta, temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4H2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H2 ) == MPIRANK ) && (N2==0)){
   #else
   if (N2==0){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSL + TwoSR - TwoS1);
               const double factor = fase * sqrt((TwoSL+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSRdown, TwoS1, TwoSR, TwoSL, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (Irrep == denBK->gIrrep(l_beta)){
                                    double fact = factor * Prob->gMxElement(l_alpha,l_beta,theindex+1,theindex+1);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    int inc = 1;
                                    daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                              
                              int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,2,TwoS1,NR+1,TwoSRdown,IRdown);
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta, temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4H3A
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H3A ) == MPIRANK ) && (N2==1)){
   #else
   if (N2==1){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSL+TwoSR+TwoSLdown+TwoSRdown+TwoJdown+1-TwoS1);
                  const double factor1 = fase * sqrt((TwoJ+1)*(TwoJdown+1)*(TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSL, TwoSRdown, TwoS1, TwoJdown, 1, TwoSLdown) * Wigner::wigner6j(TwoJ, 1, TwoS1, TwoSRdown, TwoSL, TwoSR);
                  
                  double factor2 = 0.0;
                  if (TwoJ == TwoJdown){
                     fase = phase(TwoSL+TwoSRdown+TwoJ+3+2*TwoS1);
                     factor2 = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSRdown, TwoJ, TwoSR, TwoSL, 1);
                  }
                  
                  for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                     int ILdown = Irreps::directProd(IL,Irrep);
                     int IRdown = Irreps::directProd(IR,Irrep);
                     int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                        bool isPossibleLeft = false;
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                        }
                        bool isPossibleRight = false;
                        for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                           if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                        }
                        if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                           for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                              if (Irrep == denBK->gIrrep(l_alpha)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                    if (Irrep == denBK->gIrrep(l_delta)){
                                       double fact = factor1 * Prob->gMxElement(l_alpha,theindex+1,theindex+1,l_delta);
                                       if (TwoJ==TwoJdown){ fact += factor2 * Prob->gMxElement(l_alpha,theindex+1,l_delta,theindex+1); }
                                       double * LblockR = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                       int inc = 1;
                                       daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                              
                                 int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR-1,TwoSRdown,IRdown);
                                 dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta, temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4H3B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H3B ) == MPIRANK ) && (N2==2)){
   #else
   if (N2==2){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSL+TwoSRdown+3-TwoS1);
               const double factor = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSRdown, TwoS1, TwoSR, TwoSL, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL-1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (Irrep == denBK->gIrrep(l_alpha)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                        if (Irrep == denBK->gIrrep(l_delta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (Irrep == denBK->gIrrep(l_alpha)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                                 if (Irrep == denBK->gIrrep(l_delta)){
                                    double fact = factor * ( Prob->gMxElement(l_alpha,theindex+1,theindex+1,l_delta) - 2 * Prob->gMxElement(l_alpha,theindex+1,l_delta,theindex+1) );
                                    double * LblockR = Lright[l_delta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    int inc = 1;
                                    daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                              
                              int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,2,TwoS1,NR-1,TwoSRdown,IRdown);
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta, temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4H4A
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H4A ) == MPIRANK ) && (N2==1)){
   #else
   if (N2==1){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){

            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSL+TwoSR+TwoSLdown+TwoSRdown+TwoJ+1-TwoS1);
                  const double factor1 = fase * sqrt((TwoJ+1)*(TwoJdown+1)*(TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSLdown, TwoSR, TwoS1, TwoJ, 1, TwoSL) * Wigner::wigner6j(TwoJdown, 1, TwoS1, TwoSR, TwoSLdown, TwoSRdown);
                  
                  double factor2 = 0.0;
                  if (TwoJ == TwoJdown){
                     fase = phase(TwoSLdown+TwoSR+TwoJ+3+2*TwoS1);
                     factor2 = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL, TwoSR, TwoJ, TwoSRdown, TwoSLdown, 1);
                  }
                  
                  for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                     int ILdown = Irreps::directProd(IL,Irrep);
                     int IRdown = Irreps::directProd(IR,Irrep);
                     int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                        bool isPossibleLeft = false;
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                        }
                        bool isPossibleRight = false;
                        for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                           if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                        }
                        if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                           for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                              if (Irrep == denBK->gIrrep(l_gamma)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (Irrep == denBK->gIrrep(l_beta)){
                                       double fact = factor1 * Prob->gMxElement(l_gamma,theindex+1,theindex+1,l_beta);
                                       if (TwoJ==TwoJdown){ fact += factor2 * Prob->gMxElement(l_gamma,theindex+1,l_beta,theindex+1); }
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                       int inc = 1;
                                       daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 double * LblockL = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                              
                                 int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR+1,TwoSRdown,IRdown);
                                 dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta, temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4H4B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4H4B ) == MPIRANK ) && (N2==2)){
   #else
   if (N2==2){
   #endif
   
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSLdown+TwoSR+3-TwoS1);
               const double factor = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL, TwoSR, TwoS1, TwoSRdown, TwoSLdown, 1);
                  
               for (int Irrep=0; Irrep < (denBK->getNumberOfIrreps()); Irrep++){
               
                  int ILdown = Irreps::directProd(IL,Irrep);
                  int IRdown = Irreps::directProd(IR,Irrep);
                  int dimLdown = denBK->gCurrentDim(theindex,  NL+1,TwoSLdown,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+1,TwoSRdown,IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                     bool isPossibleLeft = false;
                     for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                        if (Irrep == denBK->gIrrep(l_gamma)){ isPossibleLeft = true; }
                     }
                     bool isPossibleRight = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (Irrep == denBK->gIrrep(l_beta)){ isPossibleRight = true; }
                     }
                     if ( (isPossibleLeft) && (isPossibleRight) ){
                        
                        for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                           if (Irrep == denBK->gIrrep(l_gamma)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (Irrep == denBK->gIrrep(l_beta)){
                                    double fact = factor * ( Prob->gMxElement(l_gamma,theindex+1,theindex+1,l_beta) - 2 * Prob->gMxElement(l_gamma,theindex+1,l_beta,theindex+1) );
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    int inc = 1;
                                    daxpy_(&size,&fact,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              double * LblockL = Lleft[theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                              
                              int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,2,TwoS1,NR+1,TwoSRdown,IRdown);
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta, temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                              
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4I(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int ILdown = Irreps::directProd(IL,denBK->gIrrep(theindex));
   int inc = 1;
   double beta = 1.0; //add
   
   //4I1A and 4I1B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4I1AB ) == MPIRANK ) && (N1>0) && (N2==0)){
   #else
   if ((N1>0) && (N2==0)){
   #endif
   
      int TwoJdown = ((N1==2)?1:0);
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if (abs(TwoSLdown-TwoSR)<=TwoJdown){
      
            int dimLdown = denBK->gCurrentDim(theindex,  NL-1, TwoSLdown, ILdown);
            if (dimLdown>0){
               
               int size = dimLdown * dimLup;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
               int number = 0;
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     number++;
                     double alpha = Prob->gMxElement(l_index,theindex,theindex+1,theindex+1);
                     double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                  }
               }
               
               if (number>0){
               
                  double factor = -1.0; //4I1B
                  if (N1==1){
                     int fase = phase(TwoSR + 1 - TwoSL);
                     factor = fase * sqrt((TwoSL+1.0)/(TwoSR+1.0));
                  }
                  int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, N1-1, 2, TwoJdown, NR, TwoSR, IR);
                  dgemm_(&trans, &notrans, &dimLup, &dimR, &dimLdown, &factor, temp, &dimLdown, memS+denS->gKappa2index(memSkappa), &dimLdown, &beta, memHeff+denS->gKappa2index(ikappa), &dimLup);
               }
            }
         }
      }
   }
   
   //4I2A and 4I2B
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4I2AB ) == MPIRANK ) && (N1<2) && (N2==2)){
   #else
   if ((N1<2) && (N2==2)){
   #endif
   
      int TwoJdown = ((N1==0)?1:0);
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if (abs(TwoSLdown-TwoSR)<=TwoJdown){
         
            int dimLdown = denBK->gCurrentDim(theindex,  NL+1, TwoSLdown, ILdown);
            if (dimLdown>0){
               
               int size = dimLdown * dimLup;
               for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
               
               int number = 0;
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     number++;
                     double alpha = Prob->gMxElement(l_index,theindex,theindex+1,theindex+1);
                     double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                  }
               }
               
               if (number>0){
               
                  double factor = -1.0; //4I2B
                  if (N1==0){ //4I2A
                     int fase = phase(TwoSR + 1 - TwoSLdown);
                     factor = fase * sqrt((TwoSLdown+1.0)/(TwoSR+1.0));
                  }
                  int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, N1+1, 0, TwoJdown, NR, TwoSR, IR);
                  dgemm_(&notrans, &notrans, &dimLup, &dimR, &dimLdown, &factor, temp, &dimLup, memS+denS->gKappa2index(memSkappa), &dimLdown, &beta, memHeff+denS->gKappa2index(ikappa), &dimLup);
               }
            }
         }
      }
   }
   
   //4I3A and 4I3B and 4I3C and 4I3D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4I3ABCD ) == MPIRANK ) && (N1<2) && (N2>0)){
   #else
   if ((N1<2) && (N2>0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){

         int dimLdown = denBK->gCurrentDim(theindex,  NL-1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            //int N1down = N1+1;
            //int N2down = N2;
            int TwoSdownSum = ((N1==0)?1:0) + ((N2==1)?1:0);
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart ; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int size = dimLdown * dimLup;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                  
                  double prefact = 0.0;
                  if ((N1==0)&&(N2==1)){
                     int fase = phase(TwoSLdown+TwoSR+2);
                     prefact = fase * sqrt((TwoSL+1.0)*(TwoJdown+1)) * Wigner::wigner6j(TwoJdown,1,1,TwoSL,TwoSLdown,TwoSR);
                  }
                  if ((N1==1)&&(N2==1)){
                     int fase = phase(TwoSLdown+TwoSR+3);
                     prefact = fase * sqrt((TwoJ+1)*(TwoSL+1.0)) * Wigner::wigner6j(TwoJ, 1, 1, TwoSLdown, TwoSL, TwoSR);
                  }
                  if ((N1==0)&&(N2==2)){
                     prefact = 1.0;
                  }
                  if ((N1==1)&&(N2==2)){
                     int fase = phase(TwoSR+1-TwoSL);
                     prefact = fase * sqrt((TwoSL+1.0)/(TwoSR+1.0));
                  }
               
                  int number = 0;
                  for (int l_index=0; l_index<theindex; l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                        number++;
                        double alpha = 0.0;
                        if ((N1==0)&&(N2==1)){
                                alpha = prefact * ( Prob->gMxElement(l_index,theindex+1,theindex,theindex+1)
                           + ((TwoJdown==0)?1:-1) * Prob->gMxElement(l_index,theindex+1,theindex+1,theindex) );
                        }
                        if ((N1==1)&&(N2==1)){
                                         alpha  =   prefact * Prob->gMxElement(l_index,theindex+1,theindex,theindex+1);
                           if (TwoJ==0){ alpha += sqrt(2.0) * Prob->gMxElement(l_index,theindex+1,theindex+1,theindex); }
                        }
                        if (N2==2){
                           alpha = prefact * ( 2 * Prob->gMxElement(l_index,theindex+1,theindex,theindex+1)
                                                 - Prob->gMxElement(l_index,theindex+1,theindex+1,theindex) );
                        }
                        
                        double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     }
                  }
               
                  if (number>0){
      
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, N1+1, N2, TwoJdown, NR, TwoSR, IR);
                     dgemm_(&trans, &notrans, &dimLup, &dimR, &dimLdown, &factor, temp, &dimLdown, memS+denS->gKappa2index(memSkappa), &dimLdown, &beta, memHeff+denS->gKappa2index(ikappa), &dimLup);
                  }
               }   
            }
         }
      }
   }
   
   //4I4A and 4I4B and 4I4C and 4I4D
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_4I4ABCD ) == MPIRANK ) && (N1>0) && (N2>0)){
   #else
   if ((N1>0) && (N2>0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){

         int dimLdown = denBK->gCurrentDim(theindex,  NL+1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            //int N1down = N1-1;
            //int N2down = N2;
            int TwoSdownSum = ((N1==2)?1:0) + ((N2==1)?1:0);
            int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoSdownSum<2)) ? TwoSdownSum : 0;
            for (int TwoJdown = TwoJstart ; TwoJdown <= TwoSdownSum ; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
               
                  int size = dimLdown * dimLup;
                  for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                  
                  double prefact = 0.0;
                  if ((N1==1)&&(N2==1)){
                     int fase = phase(TwoSL+TwoSR+2);
                     prefact = fase * sqrt((TwoSLdown+1.0)*(TwoJ+1)) * Wigner::wigner6j(TwoJ,1,1,TwoSLdown,TwoSL,TwoSR);
                  }
                  if ((N1==2)&&(N2==1)){
                     int fase = phase(TwoSL+TwoSR+3);
                     prefact = fase * sqrt((TwoJdown+1)*(TwoSLdown+1.0)) * Wigner::wigner6j(TwoJdown, 1, 1, TwoSL, TwoSLdown, TwoSR);
                  }
                  if ((N1==1)&&(N2==2)){
                     prefact = 1.0;
                  }
                  if ((N1==2)&&(N2==2)){
                     int fase = phase(TwoSR+1-TwoSLdown);
                     prefact = fase * sqrt((TwoSLdown+1.0)/(TwoSR+1.0));
                  }
               
                  int number = 0;
                  for (int l_index=0; l_index<theindex; l_index++){
                     if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                        number++;
                        double alpha = 0.0;
                        if ((N1==1)&&(N2==1)){
                            alpha = prefact * ( Prob->gMxElement(l_index,theindex+1,theindex,theindex+1)
                           + ((TwoJ==0)?1:-1) * Prob->gMxElement(l_index,theindex+1,theindex+1,theindex) );
                        }
                        if ((N1==2)&&(N2==1)){
                                             alpha  =   prefact * Prob->gMxElement(l_index,theindex+1,theindex,theindex+1);
                           if (TwoJdown==0){ alpha += sqrt(2.0) * Prob->gMxElement(l_index,theindex+1,theindex+1,theindex); }
                        }
                        if (N2==2){
                           alpha = prefact * ( 2 * Prob->gMxElement(l_index,theindex+1,theindex,theindex+1)
                                                 - Prob->gMxElement(l_index,theindex+1,theindex+1,theindex) );
                        }
                        
                        double * Lblock = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                        daxpy_(&size,&alpha,Lblock,&inc,temp,&inc);
                     }
                  }
               
                  if (number>0){
      
                     double factor = 1.0;
                     int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, N1-1, N2, TwoJdown, NR, TwoSR, IR);
                     dgemm_(&notrans, &notrans, &dimLup, &dimR, &dimLdown, &factor, temp, &dimLup, memS+denS->gKappa2index(memSkappa), &dimLdown, &beta, memHeff+denS->gKappa2index(ikappa), &dimLup);
                  }
               }   
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram4J1and4J2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Aright) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int IRdown = Irreps::directProd(IR,Aright->get_irrep());

   //4J1A.spin0
   if ((N1==0) && (N2==0)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR+2,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }
   
   //4J1B.spin0
   if ((N1==1) && (N2==0)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,2,1,1,NR+2,TwoSR,IRdown);
         double alpha = - sqrt(0.5);
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }
   
   //4J1C.spin0
   if ((N1==0) && (N2==1)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,2,1,NR+2,TwoSR,IRdown);
         double alpha = - sqrt(0.5);
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }
   
   //4J1D.spin0 and 4J2A.spin0
   if ((N1==1) && (N2==1) && (TwoJ==0)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,2,2,0,NR+2,TwoSR,IRdown);
         double alpha = -1.0;
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
      
      dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,0,0,0,NR-2,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }

   //4J2B.spin0
   if ((N1==2) && (N2==1)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,0,1,NR-2,TwoSR,IRdown);
         double alpha = - sqrt(0.5);
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }
   
   //4J2C.spin0
   if ((N1==1) && (N2==2)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,0,1,1,NR-2,TwoSR,IRdown);
         double alpha = - sqrt(0.5);
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }

   //4J2D.spin0
   if ((N1==2) && (N2==2)){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
      if (dimRdown>0){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR-2,TwoSR,IRdown);
         double alpha = -1.0;
         double beta = 1.0;
         double * Ablock = Aright->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Ablock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   }

}

void CheMPS2::Heff::addDiagram4J1and4J2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Bright) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int IRdown = Irreps::directProd(IR,Bright->get_irrep());

   //4J1A.spin1
   if ((N1==0) && (N2==0)){ //TwoSL = TwoSR --> TwoSR can be what it wants to be.
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
   
         int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
         if (dimRdown>0){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,2,NR+2,TwoSRdown,IRdown);
            double alpha = sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
            double beta = 1.0;
            double * Bblock = Bright->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
            dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
         
         }
      }
   }
   
   //4J1B.spin1
   if ((N1==1) && (N2==0)){
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=1){
   
            int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
            if (dimRdown>0){
         
               int memSkappa = denS->gKappa(NL,TwoSL,IL,2,1,1,NR+2,TwoSRdown,IRdown);
               int fase = phase(TwoSL+TwoSRdown+3);
               double alpha = fase * sqrt(3.0 * (TwoSRdown+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0;
               double * Bblock = Bright->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
               dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
            }
         }
      }
   }
   
   //4J1C.spin1
   if ((N1==0) && (N2==1)){
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=1){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
            if (dimRdown>0){
         
               int memSkappa = denS->gKappa(NL,TwoSL,IL,1,2,1,NR+2,TwoSRdown,IRdown);
               int fase = phase(TwoSL+TwoSRdown+1);
               double alpha = fase * sqrt(3.0 * (TwoSRdown+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0;
               double * Bblock = Bright->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
               dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
            }
         }
      }
   }
   
   //4J1D.spin1 and 4J2A.spin1
   if ((N1==1) && (N2==1) && (TwoJ==2)){
   
      int TwoSRdown = TwoSL;
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
      if (dimRdown>0){ //4J1D.spin1
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,2,2,0,NR+2,TwoSRdown,IRdown);
         double alpha = phase(TwoSR-TwoSRdown);
         double beta = 1.0;
         double * Bblock = Bright->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }

      
      dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
      if (dimRdown>0){ //4J2A.spin1
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,0,0,0,NR-2,TwoSRdown,IRdown);
         double alpha = sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
         double beta = 1.0;
         double * Bblock = Bright->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
         
      }

   }

   //4J2B.spin1
   if ((N1==2) && (N2==1)){
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=1){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
            if (dimRdown>0){
         
               int memSkappa = denS->gKappa(NL,TwoSL,IL,1,0,1,NR-2,TwoSRdown,IRdown);
               int fase = phase(TwoSL+TwoSR+3);
               double alpha = fase * sqrt(3.0 * (TwoSR+1)) * Wigner::wigner6j(1,1,2,TwoSRdown,TwoSR,TwoSL);
               double beta = 1.0;
               double * Bblock = Bright->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
               dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
            }
         }
      }
   }
   
   //4J2C.spin1
   if ((N1==1) && (N2==2)){
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         if (abs(TwoSL-TwoSRdown)<=1){
         
            int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
            if (dimRdown>0){
         
               int memSkappa = denS->gKappa(NL,TwoSL,IL,0,1,1,NR-2,TwoSRdown,IRdown);
               int fase = phase(TwoSL+TwoSR+1);
               double alpha = fase * sqrt(3.0 * (TwoSR+1)) * Wigner::wigner6j(1,1,2,TwoSRdown,TwoSR,TwoSL);
               double beta = 1.0;
               double * Bblock = Bright->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
               dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
            
            }
         }
      }
   }

   //4J2D.spin1
   if ((N1==2) && (N2==2)){ //TwoSL == TwoSR --> TwoSRdown can be what it wants to be.
   
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
   
         int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
         if (dimRdown>0){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,2,NR-2,TwoSRdown,IRdown);
            double alpha = phase(TwoSR-TwoSRdown);
            double beta = 1.0;
            double * Bblock = Bright->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
            dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Bblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4J3and4J4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Cright) const{

   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);

   int IRdown = Irreps::directProd(IR,Cright->get_irrep());
   int theindex = denS->gIndex();
   
   int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSR,IRdown);
   if (dimRdown>0){

      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
   
      int N1 = denS->gN1(ikappa);
      int N2 = denS->gN2(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
   
      int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
      int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
      char trans = 'T';
      char notrans = 'N';

      //4J3A.spin0
      if ((N1==1) && (N2==0)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,0,1,1,NR,TwoSR,IRdown);
         double alpha = sqrt(0.5);
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
         
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
         
      }
   
      //4J3B.spin0
      if ((N1==2) && (N2==0)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
         
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   
      //4J3C.spin0
      if ((N1==1) && (N2==1) && (TwoJ==0)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,0,2,0,NR,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
         
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   
      //4J3D.spin0
      if ((N1==2) && (N2==1)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,2,1,NR,TwoSR,IRdown);
         double alpha = -sqrt(0.5);
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
         
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }

      //4J4A.spin0
      if ((N1==0) && (N2==1)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,0,1,NR,TwoSR,IRdown);
         double alpha = sqrt(0.5);
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
         
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   
      //4J4B.spin0
      if ((N1==1) && (N2==1) && (TwoJ==0)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,2,0,0,NR,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
         
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }

      //4J4C.spin0
      if ((N1==0) && (N2==2)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR,TwoSR,IRdown);
         double alpha = 1.0;
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
         
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   
      //4J4D.spin0
      if ((N1==1) && (N2==2)){
      
         int memSkappa = denS->gKappa(NL,TwoSL,IL,2,1,1,NR,TwoSR,IRdown);
         double alpha = -sqrt(0.5);
         double beta = 1.0;
         double * ptr = Cright->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
         
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
      }
   
   }
   
}

void CheMPS2::Heff::addDiagram4J3and4J4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dright) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   
   int theindex = denS->gIndex();
   int dimL   = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';

   int IRdown = Irreps::directProd(IR,Dright->get_irrep());
   
   for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
   
      int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSRdown,IRdown);
      if (dimRdown>0){

         //4J3A.spin1
         if ((N1==1) && (N2==0) && (abs(TwoSL-TwoSRdown)<=1)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,0,1,1,NR,TwoSRdown,IRdown);
            int fase = phase(TwoSL+TwoSRdown+3);
            double alpha = fase * sqrt(3.0*(TwoSRdown+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
         
            dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
         
         }
   
         //4J3C.spin1
         if ((N1==2) && (N2==0)){ //TwoSL==TwoSR and hence TwoSRdown can be what it wants
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,2,NR,TwoSRdown,IRdown);
            double alpha = - sqrt((TwoSRdown+1.0)/(TwoSR+1.0));
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
         
            dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
   
         //4J3B.spin1
         if ((N1==1) && (N2==1) && (TwoJ==2) && (TwoSL==TwoSRdown)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,0,2,0,NR,TwoSRdown,IRdown);
            double alpha = phase(TwoSR-TwoSRdown);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
         
            dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
   
         //4J3D.spin1
         if ((N1==2) && (N2==1) && (abs(TwoSL-TwoSRdown)<=1)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,2,1,NR,TwoSRdown,IRdown);
            int fase = phase(TwoSL+TwoSRdown+3);
            double alpha = fase * sqrt(3.0*(TwoSRdown+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
         
            dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }

         //4J4A.spin1
         if ((N1==0) && (N2==1) && (abs(TwoSL-TwoSRdown)<=1)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,0,1,NR,TwoSRdown,IRdown);
            int fase = phase(TwoSL+TwoSR+3);
            double alpha = fase * sqrt(3.0*(TwoSR+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
         
            dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
   
         //4J4B.spin1
         if ((N1==1) && (N2==1) && (TwoJ==2) && (TwoSL==TwoSRdown)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,2,0,0,NR,TwoSRdown,IRdown);
            double alpha = - sqrt((TwoSR+1.0)/(TwoSRdown+1.0));
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
         
            dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }

         //4J4C.spin1
         if ((N1==0) && (N2==2)){ // TwoSL == TwoSR --> TwoSRdown can be what it wants to be.
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,2,NR,TwoSRdown,IRdown);
            double alpha = phase(TwoSR-TwoSRdown);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
         
            dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
   
         //4J4D.spin1
         if ((N1==1) && (N2==2) && (abs(TwoSL-TwoSRdown)<=1)){
      
            int memSkappa = denS->gKappa(NL,TwoSL,IL,2,1,1,NR,TwoSRdown,IRdown);
            int fase = phase(TwoSL+TwoSR+3);
            double alpha = fase * sqrt(3.0*(TwoSR+1)) * Wigner::wigner6j(1,1,2,TwoSR,TwoSRdown,TwoSL);
            double beta = 1.0;
            double * ptr = Dright->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
         
            dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,ptr,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa), &dimL);
      
         }
      }
   }
   
}


void CheMPS2::Heff::addDiagram4K1and4K2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Aright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4K1A.spin0
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
   
            int fase = phase(TwoSLdown + TwoSR + TwoJ + 1 + 2*TwoS1);
            const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJ+1)) * Wigner::wigner6j(TwoS1,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Aright[theindex+1-l_index][0]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,0,TwoS1,NR-2,TwoSR,IRdown);
                     double * blockA = Aright[theindex+1-l_index][0]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRdown,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4K1B.spin0
   if (N2==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSLdown + TwoSR + TwoJdown + 2 + 2*TwoS1);
               const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Aright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR-2,TwoSR,IRdown);
                        double * blockA = Aright[theindex+1-l_index][0]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K2A.spin0
   if (N2==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSL + TwoSR + TwoJdown + 1 + 2*TwoS1);
               const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Aright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR+2,TwoSR,IRdown);
                        double * blockA = Aright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K2B.spin0
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
   
            int fase = phase(TwoSL + TwoSR + TwoJ + 2 + 2*TwoS1);
            const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJ+1)) * Wigner::wigner6j(TwoJ,TwoS1,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Aright[theindex+1-l_index][0]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,2,TwoS1,NR+2,TwoSR,IRdown);
                     double * blockA = Aright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRup,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4L1and4L2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Aright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4L1A.spin0
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
   
            int fase = phase(TwoSLdown + TwoSR + 2 + TwoS2);
            const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJ+1))
                                * Wigner::wigner6j(TwoS2,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Aright[theindex-l_index][1]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,0,N2,TwoS2,NR-2,TwoSR,IRdown);
                     double * blockA = Aright[theindex-l_index][1]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRdown,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
   //4L1B.spin0
   if (N1==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSLdown + TwoSR + 3 + TwoS2);
               const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Aright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR-2,TwoSR,IRdown);
                        double * blockA = Aright[theindex-l_index][1]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L2A.spin0
   if (N1==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSL + TwoSR + 2 + TwoS2);
               const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJdown+1))
                                   * Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Aright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR+2,TwoSR,IRdown);
                        double * blockA = Aright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L2B.spin0
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
         
            int fase = phase(TwoSL + TwoSR + 3 + TwoS2);
            const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJ+1))
                                * Wigner::wigner6j(TwoJ,TwoS2,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Aright[theindex-l_index][1]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,2,N2,TwoS2,NR+2,TwoSR,IRdown);
                     double * blockA = Aright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockA,&dimRup,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4K1and4K2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Bright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4K1A.spin1
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = (TwoS1==1)?-1:1;
               const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSL+1) * (TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS1);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Bright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,0,TwoS1,NR-2,TwoSRdown,IRdown);
                        double * blockB = Bright[theindex+1-l_index][0]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K1B.spin1
   if (N2==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSR-TwoSRdown+TwoSLdown-TwoSL+3+2*TwoS1);
                  const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSL+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS1);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Bright[theindex+1-l_index][0]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR-2,TwoSRdown,IRdown);
                           double * blockB = Bright[theindex+1-l_index][0]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRdown,&beta, temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K2A.spin1
   if (N2==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = (TwoS1==1)?-1:1;
                  const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSLdown+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS1);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Bright[theindex+1-l_index][0]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR+2,TwoSRdown,IRdown);
                           double * blockB = Bright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRup,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K2B.spin1
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSRdown-TwoSR+TwoSL-TwoSLdown+3+2*TwoS1);
               const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSLdown+1) * (TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS1);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Bright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,2,TwoS1,NR+2,TwoSRdown,IRdown);
                        double * blockB = Bright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4L1and4L2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Bright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4L1A.spin1
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(1+TwoS2-TwoJ);
               const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSL+1) * (TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS2);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Bright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,0,N2,TwoS2,NR-2,TwoSRdown,IRdown);
                        double * blockB = Bright[theindex-l_index][1]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L1B.spin1
   if (N1==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){

                  int fase = phase(TwoSR-TwoSRdown+TwoSLdown-TwoSL+TwoS2-TwoJdown); //bug fixed
                  const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSL+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS2);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Bright[theindex-l_index][1]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-2, TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR-2,TwoSRdown,IRdown);
                           double * blockB = Bright[theindex-l_index][1]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRdown,&beta, temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L2A.spin1
   if (N1==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(1+TwoS2-TwoJdown);
                  const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSLdown+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS2);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Bright[theindex-l_index][1]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR+2,TwoSRdown,IRdown);
                           double * blockB = Bright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRup,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L2B.spin1
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSRdown-TwoSR+TwoSL-TwoSLdown+TwoS2-TwoJ);
               const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSLdown+1) * (TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS2);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_absigma( l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Bright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+2, TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,2,N2,TwoS2,NR+2,TwoSRdown,IRdown);
                        double * blockB = Bright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,blockB,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::Heff::addDiagram4K3and4K4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Cright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4K3A.spin0
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
         
            int fase = phase(TwoSL + TwoSR + TwoJ + 2*TwoS1);
            const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJ+1)) * Wigner::wigner6j(TwoJ,TwoS1,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Cright[theindex+1-l_index][0]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,0,TwoS1,NR,TwoSR,IRdown);
                     double * ptr = Cright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
                     
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   //4K3B.spin0
   if (N2==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSL + TwoSR + TwoJdown + 1 + 2*TwoS1);
               const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Cright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSR,IRdown);
                        double * ptr = Cright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
                     
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K4A.spin0
   if (N2==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS1==0)) ? TwoS1+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSLdown + TwoSR + TwoJdown + 2*TwoS1);
               const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Cright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSR,IRdown);
                        double * ptr = Cright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
                     
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K4B.spin0
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
         
            int fase = phase(TwoSLdown + TwoSR + 1 + TwoJ + 2*TwoS1);
            const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJ+1)) * Wigner::wigner6j(TwoJ,TwoS1,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Cright[theindex+1-l_index][0]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,2,TwoS1,NR,TwoSR,IRdown);
                     double * ptr = Cright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
                     
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                   
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
  
}

void CheMPS2::Heff::addDiagram4L3and4L4spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Cright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4L3A.spin0
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
         
            int fase = phase(TwoSL + TwoSR + 1 + TwoS2);
            const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJ+1)) * Wigner::wigner6j(TwoJ,TwoS2,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Cright[theindex-l_index][1]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,0,N2,TwoS2,NR,TwoSR,IRdown);
                     double * ptr = Cright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
                     
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   //4L3B.spin0
   if (N1==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSL + TwoSR + 2 + TwoS2);
               const double factor = fase * sqrt(0.5 * (TwoSLdown+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Cright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSR,IRdown);
                        double * ptr = Cright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
                     
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L4A.spin0
   if (N1==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         
         int TwoJstart = ((TwoSLdown!=TwoSR) || (TwoS2==0)) ? TwoS2+1 : 0;
         for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
            if ((abs(TwoSLdown-TwoSR)<=TwoJdown) && (TwoSLdown>=0)){
            
               int fase = phase(TwoSLdown + TwoSR + 1 + TwoS2);
               const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJdown+1)) * Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Cright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSR,IRdown);
                        double * ptr = Cright[theindex-l_index][1]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
                     
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L4B.spin0
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
         
            int fase = phase(TwoSLdown + TwoSR + 2 + TwoS2);
            const double factor = fase * sqrt(0.5 * (TwoSL+1) * (TwoJ+1)) * Wigner::wigner6j(TwoJ,TwoS2,1,TwoSLdown,TwoSL,TwoSR);
         
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR, Cright[theindex-l_index][1]->get_irrep() );
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSR,     IRdown);
                  
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2,TwoS2,NR,TwoSR,IRdown);
                     double * ptr = Cright[theindex-l_index][1]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
                     
                     double beta = 0.0; //set
                     double alpha = factor;
                     
                     dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                     beta = 1.0; //add
                     alpha = 1.0;
                     double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                   
                     dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
  
}

void CheMPS2::Heff::addDiagram4K3and4K4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Dright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS1 = (N1==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4K3A.spin1
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSL - TwoSLdown + 1 + 2*TwoS1);
               const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSLdown+1) * (TwoJ+1))
                                   * Wigner::wigner9j(2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS1);
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Dright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,0,TwoS1,NR,TwoSRdown,IRdown);
                        double * ptr = Dright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K3B.spin1
   if (N2==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSR - TwoSRdown + 2*TwoS1);
                  const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSLdown+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS1);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Dright[theindex+1-l_index][0]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSRdown,IRdown);
                           double * ptr = Dright[theindex+1-l_index][0]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K4A.spin1
   if (N2==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS1==0)) ? TwoS1+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS1+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSLdown + 1 - TwoSL + 2*TwoS1);
                  const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSL+1) * (TwoJdown+1))
                                      * Wigner::wigner9j(2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS1);
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Dright[theindex+1-l_index][0]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSRdown,IRdown);
                           double * ptr = Dright[theindex+1-l_index][0]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4K4B.spin1
   if (N2==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSR - TwoSRdown + 2*TwoS1);
               const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSL+1) * (TwoJ+1))
                                   * Wigner::wigner9j( 2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS1 );
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex+1 ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Dright[theindex+1-l_index][0]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                   
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,2,TwoS1,NR,TwoSRdown,IRdown);
                        double * ptr = Dright[theindex+1-l_index][0]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                   
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
  
}

void CheMPS2::Heff::addDiagram4L3and4L4spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorOperator *** Dright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
  
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int TwoS2 = (N2==1)?1:0;
   
   int theindex = denS->gIndex();
   int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   
   //4L3A.spin1
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSL - TwoSLdown + 2 + TwoS2 - TwoJ);
               const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSLdown+1) * (TwoJ+1))
                                   * Wigner::wigner9j( 2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS2 );
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Dright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,0, N2,TwoS2,NR,TwoSRdown,IRdown);
                        double * ptr = Dright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L3B.spin1
   if (N1==2){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSR - TwoSRdown + 1 + TwoS2 - TwoJdown);
                  const double factor = fase * sqrt(3.0 * (TwoSR+1) * (TwoSLdown+1) * (TwoJdown+1))
                                      * Wigner::wigner9j( 2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS2 );
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Dright[theindex-l_index][1]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSRdown,IRdown);
                           double * ptr = Dright[theindex-l_index][1]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRup,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLup,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L4A.spin1
   if (N1==0){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            
            int TwoJstart = ((TwoSLdown!=TwoSRdown) || (TwoS2==0)) ? TwoS2+1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=TwoS2+1; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  int fase = phase(TwoSLdown + 2 - TwoSL + TwoS2 - TwoJdown);
                  const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSL+1) * (TwoJdown+1))
                                      * Wigner::wigner9j( 2, TwoSRdown, TwoSR, 1, TwoSLdown, TwoSL, 1, TwoJdown, TwoS2 );
         
                  for (int l_index=0; l_index<theindex; l_index++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                        int IRdown = Irreps::directProd(IR, Dright[theindex-l_index][1]->get_irrep() );
                  
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                  
                        if ((dimLdown>0) && (dimRdown>0)){
                  
                           int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSRdown,IRdown);
                           double * ptr = Dright[theindex-l_index][1]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double beta = 0.0; //set
                           double alpha = factor;
                     
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                           beta = 1.0; //add
                           alpha = 1.0;
                           double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //4L4B.spin1
   if (N1==1){
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=TwoS2) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSR - TwoSRdown + 1 + TwoS2 - TwoJ);
               const double factor = fase * sqrt(3.0 * (TwoSRdown+1) * (TwoSL+1) * (TwoJ+1))
                                   * Wigner::wigner9j( 2, TwoSR, TwoSRdown, 1, TwoSL, TwoSLdown, 1, TwoJ, TwoS2 );
         
               for (int l_index=0; l_index<theindex; l_index++){
               
                  #ifdef CHEMPS2_MPI_COMPILATION
                  if ( MPIchemps2::owner_cdf( Prob->gL(), l_index, theindex ) == MPIRANK )
                  #endif
                  {
                     int ILdown = Irreps::directProd(IL, denBK->gIrrep(l_index));
                     int IRdown = Irreps::directProd(IR, Dright[theindex-l_index][1]->get_irrep() );
                  
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR,   TwoSRdown, IRdown);
                   
                     if ((dimLdown>0) && (dimRdown>0)){
                  
                        int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2,TwoS2,NR,TwoSRdown,IRdown);
                        double * ptr = Dright[theindex-l_index][1]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
                        double beta = 0.0; //set
                        double alpha = factor;
                     
                        dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,ptr,&dimRdown,&beta,temp,&dimLdown);
                     
                        beta = 1.0; //add
                        alpha = 1.0;
                        double * blockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                   
                        dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,blockL,&dimLdown,temp,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     }
                  }
               }
            }
         }
      }
   }
  
}


