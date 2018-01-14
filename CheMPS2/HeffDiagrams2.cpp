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

void CheMPS2::Heff::addDiagram2a1spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Atensors, TensorS0 **** S0tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int l_alpha=0; l_alpha<theindex; l_alpha++){
         for (int l_beta=l_alpha; l_beta<theindex; l_beta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( l_alpha, l_beta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,S0tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->get_irrep());
               int IRdown = Irreps::directProd(IR,Atensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->get_irrep());
               int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,N1,N2,TwoJ,NR-2,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL-2,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
               
                  double * BlockS0 = S0tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                  double * BlockA = Atensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockS0,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockA,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   
   } else {
      
      for (int l_gamma=theindex+2; l_gamma<Prob->gL(); l_gamma++){
         for (int l_delta=l_gamma; l_delta<Prob->gL(); l_delta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( l_gamma, l_delta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,Atensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->get_irrep());
               int IRdown = Irreps::directProd(IR,S0tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->get_irrep());
               int memSkappa = denS->gKappa(NL-2,TwoSL,ILdown,N1,N2,TwoJ,NR-2,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL-2,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IRdown);
               
                  double * BlockA = Atensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->gStorage(NL-2,TwoSL,ILdown,NL,TwoSL,IL);
                  double * BlockS0 = S0tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->gStorage(NR-2,TwoSR,IRdown,NR,TwoSR,IR);

                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockA,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockS0,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2a2spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Atensors, TensorS0 **** S0tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int l_alpha=0; l_alpha<theindex; l_alpha++){
         for (int l_beta=l_alpha; l_beta<theindex; l_beta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( l_alpha, l_beta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,S0tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->get_irrep());
               int IRdown = Irreps::directProd(IR,Atensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->get_irrep());
               int memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,N1,N2,TwoJ,NR+2,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL+2,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
               
                  double * BlockS0 = S0tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                  double * BlockA = Atensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockS0,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockA,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   
   } else {
      
      for (int l_gamma=theindex+2; l_gamma<Prob->gL(); l_gamma++){
         for (int l_delta=l_gamma; l_delta<Prob->gL(); l_delta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( l_gamma, l_delta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,Atensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->get_irrep());
               int IRdown = Irreps::directProd(IR,S0tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->get_irrep());
               int memSkappa = denS->gKappa(NL+2,TwoSL,ILdown,N1,N2,TwoJ,NR+2,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL+2,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IRdown);
               
                  double * BlockA = Atensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->gStorage(NL,TwoSL,IL,NL+2,TwoSL,ILdown);
                  double * BlockS0 = S0tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IRdown);

                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockA,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockS0,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2a1spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Btensors, TensorS1 **** S1tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSRdown+TwoSL+TwoJ+2);
               const double thefactor = fase * sqrt((TwoSR + 1)*(TwoSL + 1.0)) * Wigner::wigner6j( TwoSLdown, TwoSRdown, TwoJ, TwoSR, TwoSL, 2 );
      
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  for (int l_beta=l_alpha+1; l_beta<theindex; l_beta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_alpha, l_beta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,S1tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->get_irrep());
                        int IRdown = Irreps::directProd(IR,Btensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->get_irrep());
                        int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,N1,N2,TwoJ,NR-2,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL-2,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
               
                           double * BlockS1 = S1tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                           double * BlockB = Btensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);
            
                           char trans = 'T';
                           char notrans = 'N';
                           double alpha = thefactor;
                           double beta = 0.0;
                  
                           dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockS1,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           alpha = beta = 1.0;
                  
                           dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockB,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   
   } else {
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSRdown+TwoSL+TwoJ+2);
               const double thefactor = fase * sqrt((TwoSR + 1)*(TwoSL + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_gamma=theindex+2; l_gamma<Prob->gL(); l_gamma++){
                  for (int l_delta=l_gamma+1; l_delta<Prob->gL(); l_delta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_gamma, l_delta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,Btensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->get_irrep());
                        int IRdown = Irreps::directProd(IR,S1tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->get_irrep());
                        int memSkappa = denS->gKappa(NL-2,TwoSLdown,ILdown,N1,N2,TwoJ,NR-2,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL-2,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSRdown,IRdown);
               
                           double * BlockB = Btensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->gStorage(NL-2,TwoSLdown,ILdown,NL,TwoSL,IL);
                           double * BlockS1 = S1tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->gStorage(NR-2,TwoSRdown,IRdown,NR,TwoSR,IR);

                           char trans = 'T';
                           char notrans = 'N';
                           double alpha = thefactor;
                           double beta = 0.0;
                        
                           dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockB,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                        
                           alpha = beta = 1.0;
                        
                           dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockS1,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2a2spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Btensors, TensorS1 **** S1tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSLdown+TwoSR+TwoJ+2);
               const double thefactor = fase * sqrt((TwoSRdown + 1)*(TwoSLdown + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
         
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  for (int l_beta=l_alpha+1; l_beta<theindex; l_beta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_alpha, l_beta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,S1tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->get_irrep());
                        int IRdown = Irreps::directProd(IR,Btensors[theindex+1][l_beta-l_alpha][theindex+1-l_beta]->get_irrep());
                        int memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,N1,N2,TwoJ,NR+2,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
                
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL+2,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
               
                           double * BlockS1 = S1tensors[theindex-1][l_beta-l_alpha][theindex-1-l_beta]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                           double * BlockB  = Btensors[ theindex+1][l_beta-l_alpha][theindex+1-l_beta]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);
            
                           char trans = 'T';
                           char notr = 'N';
                           double alpha = thefactor;
                           double beta = 0.0;
                  
                           dgemm_(&notr,&notr,&dimL,&dimRdown,&dimLdown,&alpha,BlockS1,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           alpha = beta = 1.0;
                  
                           dgemm_(&notr,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockB,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   
   } else {
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSLdown+TwoSR+TwoJ+2);
               const double thefactor = fase * sqrt((TwoSRdown + 1)*(TwoSLdown + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_gamma=theindex+2; l_gamma<Prob->gL(); l_gamma++){
                  for (int l_delta=l_gamma+1; l_delta<Prob->gL(); l_delta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_absigma( l_gamma, l_delta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,Btensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->get_irrep());
                        int IRdown = Irreps::directProd(IR,S1tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->get_irrep());
                        int memSkappa = denS->gKappa(NL+2,TwoSLdown,ILdown,N1,N2,TwoJ,NR+2,TwoSRdown,IRdown);
                
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL+2,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSRdown,IRdown);
               
                           double * BlockB = Btensors[theindex-1][l_delta-l_gamma][l_gamma-theindex]->gStorage(NL,TwoSL,IL,NL+2,TwoSLdown,ILdown);
                           double * BlockS1 = S1tensors[theindex+1][l_delta-l_gamma][l_gamma-theindex-2]->gStorage(NR,TwoSR,IR,NR+2,TwoSRdown,IRdown);

                           char trans = 'T';
                           char notrans = 'N';
                           double alpha = thefactor;
                           double beta = 0.0;
                  
                           dgemm_(&notrans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockB,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           alpha = beta = 1.0;
                  
                           dgemm_(&notrans,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockS1,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                    
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2a3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Ctensors, TensorF0 **** F0tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int l_gamma=0; l_gamma<theindex; l_gamma++){
         for (int l_alpha=l_gamma+1; l_alpha<theindex; l_alpha++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_gamma, l_alpha ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,F0tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->get_irrep());
               int IRdown = Irreps::directProd(IR,Ctensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->get_irrep());
               int memSkappa = denS->gKappa(NL,TwoSL,ILdown,N1,N2,TwoJ,NR,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSR,IRdown);
               
                  //no transpose
                  double * ptr = Ctensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
                  double * BlockF0 = F0tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockF0,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,ptr,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
      
      for (int l_alpha=0; l_alpha<theindex; l_alpha++){
         for (int l_gamma=l_alpha; l_gamma<theindex; l_gamma++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_alpha, l_gamma ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,F0tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->get_irrep());
               int IRdown = Irreps::directProd(IR,Ctensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->get_irrep());
               int memSkappa = denS->gKappa(NL,TwoSL,ILdown,N1,N2,TwoJ,NR,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSR,IRdown);
               
                  //transpose
                  double * ptr = Ctensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
                  double * BlockF0 = F0tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,BlockF0,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   
   } else {
      
      for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
         for (int l_beta=l_delta+1; l_beta<Prob->gL(); l_beta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_delta, l_beta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,Ctensors[theindex-1][l_beta-l_delta][l_delta-theindex]->get_irrep());
               int IRdown = Irreps::directProd(IR,F0tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->get_irrep());
               int memSkappa = denS->gKappa(NL,TwoSL,ILdown,N1,N2,TwoJ,NR,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSR,IRdown);
               
                  //no transpose
                  double * ptr = Ctensors[theindex-1][l_beta-l_delta][l_delta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSL,ILdown);
                  double * BlockF0 = F0tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSR,IRdown);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,ptr,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&trans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockF0,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
      
      for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
         for (int l_delta=l_beta; l_delta<Prob->gL(); l_delta++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_beta, l_delta ) == MPIRANK )
            #endif
            {
               int ILdown = Irreps::directProd(IL,Ctensors[theindex-1][l_delta-l_beta][l_beta-theindex]->get_irrep());
               int IRdown = Irreps::directProd(IR,F0tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->get_irrep());
               int memSkappa = denS->gKappa(NL,TwoSL,ILdown,N1,N2,TwoJ,NR,TwoSR,IRdown);
               
               if (memSkappa!=-1){
               
                  int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSL,ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSR,IRdown);
               
                  //transpose
                  double * ptr = Ctensors[theindex-1][l_delta-l_beta][l_beta-theindex]->gStorage(NL,TwoSL,ILdown,NL,TwoSL,IL);
                  double * BlockF0 = F0tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->gStorage(NR,TwoSR,IRdown,NR,TwoSR,IR);
            
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0;
                  
                  dgemm_(&trans,&notrans,&dimL,&dimRdown,&dimLdown,&alpha,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                  beta = 1.0;
                  
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimRdown,&alpha,workspace,&dimL,BlockF0,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2a3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator **** Dtensors, TensorF1 **** F1tensors, double * workspace) const{

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
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   const bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSLdown+TwoSRdown+TwoJ+2);
               double prefactor = fase * sqrt((TwoSR + 1)*(TwoSLdown + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_gamma=0; l_gamma<theindex; l_gamma++){
                  for (int l_alpha=l_gamma+1; l_alpha<theindex; l_alpha++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_gamma, l_alpha ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,F1tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->get_irrep());
                        int IRdown = Irreps::directProd(IR,Dtensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->get_irrep());
                        int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,N1,N2,TwoJ,NR,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSRdown,IRdown);
               
                           //no transpose
                           double * ptr = Dtensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
                           double * BlockF1 = F1tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
            
                           char trans = 'T';
                           char notr = 'N';
                           double beta = 0.0;
                  
                           dgemm_(&notr,&notr,&dimL,&dimRdown,&dimLdown,&prefactor,BlockF1,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           beta = 1.0;
                  
                           dgemm_(&notr,&trans,&dimL,&dimR,&dimRdown,&beta,workspace,&dimL,ptr,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
               
               fase = phase(TwoSL+TwoSR+TwoJ+2);
               prefactor = fase * sqrt((TwoSRdown + 1)*(TwoSL + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  for (int l_gamma=l_alpha; l_gamma<theindex; l_gamma++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_alpha, l_gamma ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,F1tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->get_irrep());
                        int IRdown = Irreps::directProd(IR,Dtensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->get_irrep());
                        int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,N1,N2,TwoJ,NR,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSRdown,IRdown);
               
                           //transpose
                           double * ptr = Dtensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
                           double * BlockF1 = F1tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
            
                           char trans = 'T';
                           char notr = 'N';
                           double beta = 0.0;
                           
                           dgemm_(&trans,&notr,&dimL,&dimRdown,&dimLdown,&prefactor,BlockF1,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           beta = 1.0;
                  
                           dgemm_(&notr,&notr,&dimL,&dimR,&dimRdown,&beta,workspace,&dimL,ptr,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   
   } else {
   
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
         
            if ((TwoSLdown>=0) && (TwoSRdown>=0) && (abs(TwoSLdown-TwoSRdown)<=TwoJ)){
            
               int fase = phase(TwoSLdown+TwoSRdown+TwoJ+2);
               double prefactor = fase * sqrt((TwoSR + 1)*(TwoSLdown + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_delta=theindex+2; l_delta<Prob->gL(); l_delta++){
                  for (int l_beta=l_delta+1; l_beta<Prob->gL(); l_beta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_delta, l_beta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,Dtensors[theindex-1][l_beta-l_delta][l_delta-theindex]->get_irrep());
                        int IRdown = Irreps::directProd(IR,F1tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->get_irrep());
                        int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,N1,N2,TwoJ,NR,TwoSRdown,IRdown);
               
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSRdown,IRdown);
               
                           //no transpose
                           double * ptr = Dtensors[theindex-1][l_beta-l_delta][l_delta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSLdown,ILdown);
                           double * BlockF1 = F1tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSRdown,IRdown);
            
                           char trans = 'T';
                           char notr = 'N';
                           double beta = 0.0;
                  
                           dgemm_(&notr,&notr,&dimL,&dimRdown,&dimLdown,&prefactor,ptr,&dimL,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           beta = 1.0;
                  
                           dgemm_(&notr,&trans,&dimL,&dimR,&dimRdown,&beta,workspace,&dimL,BlockF1,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
               
               fase = phase(TwoSL+TwoSR+TwoJ+2);
               prefactor = fase * sqrt((TwoSRdown + 1)*(TwoSL + 1.0)) * Wigner::wigner6j(TwoSLdown,TwoSRdown,TwoJ,TwoSR,TwoSL,2);
      
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  for (int l_delta=l_beta; l_delta<Prob->gL(); l_delta++){
                  
                     #ifdef CHEMPS2_MPI_COMPILATION
                     if ( MPIchemps2::owner_cdf( Prob->gL(), l_beta, l_delta ) == MPIRANK )
                     #endif
                     {
                        int ILdown = Irreps::directProd(IL,Dtensors[theindex-1][l_delta-l_beta][l_beta-theindex]->get_irrep());
                        int IRdown = Irreps::directProd(IR,F1tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->get_irrep());
                        int memSkappa = denS->gKappa(NL,TwoSLdown,ILdown,N1,N2,TwoJ,NR,TwoSRdown,IRdown);
                        
                        if (memSkappa!=-1){
               
                           int dimLdown = denBK->gCurrentDim(theindex  ,NL,TwoSLdown,ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2,NR,TwoSRdown,IRdown);
               
                           //transpose
                           double * ptr = Dtensors[theindex-1][l_delta-l_beta][l_beta-theindex]->gStorage(NL,TwoSLdown,ILdown,NL,TwoSL,IL);
                           double * BlockF1 = F1tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->gStorage(NR,TwoSRdown,IRdown,NR,TwoSR,IR);
            
                           char trans = 'T';
                           char notr = 'N';
                           double beta = 0.0;
                  
                           dgemm_(&trans,&notr,&dimL,&dimRdown,&dimLdown,&prefactor,ptr,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,workspace,&dimL);
                  
                           beta = 1.0;
                  
                           dgemm_(&notr,&notr,&dimL,&dimR,&dimRdown,&beta,workspace,&dimL,BlockF1,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2b1and2b2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==0){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,IL);
      
      if (dimLdown>0){
   
         int NR = denS->gNR(ikappa);
         int TwoSR = denS->gTwoSR(ikappa);
         int IR = denS->gIR(ikappa);
   
         int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
         int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
         
         int memSkappa = denS->gKappa(NL-2,TwoSL,IL,2, denS->gN2(ikappa), denS->gTwoJ(ikappa), NR,TwoSR,IR);
         
         if (memSkappa!=-1){
         
            double * BlockA = Atensor->gStorage(NL-2,TwoSL,IL,NL,TwoSL,IL);
            char trans = 'T';
            char notrans = 'N';
            double alpha = sqrt(2.0);
            double beta = 1.0;
            
            dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,BlockA,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
         }
      }
   }
   
   if (N1==2){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,IL);
      
      if (dimLdown>0){
   
         int NR = denS->gNR(ikappa);
         int TwoSR = denS->gTwoSR(ikappa);
         int IR = denS->gIR(ikappa);
   
         int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
         int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
         
         int memSkappa = denS->gKappa(NL+2,TwoSL,IL,0, denS->gN2(ikappa), denS->gTwoJ(ikappa), NR,TwoSR,IR);
         
         if (memSkappa!=-1){
         
            double * BlockA = Atensor->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
            char notrans = 'N';
            double alpha = sqrt(2.0);
            double beta = 1.0;
            
            dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,BlockA,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2c1and2c2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==0){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimLdown = denBK->gCurrentDim(theindex,NL-2,TwoSL,IL);
      
      if (dimLdown>0){
   
         int NR = denS->gNR(ikappa);
         int TwoSR = denS->gTwoSR(ikappa);
         int IR = denS->gIR(ikappa);
   
         int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
         int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
         
         int memSkappa = denS->gKappa(NL-2,TwoSL,IL, denS->gN1(ikappa), 2, denS->gTwoJ(ikappa), NR,TwoSR,IR);
         
         if (memSkappa!=-1){
         
            double * BlockA = Atensor->gStorage(NL-2,TwoSL,IL,NL,TwoSL,IL);
            char trans = 'T';
            char notrans = 'N';
            double alpha = sqrt(2.0);
            double beta = 1.0;
            
            dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,BlockA,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
         }
      }
   }
   
   if (N2==2){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimLdown = denBK->gCurrentDim(theindex,NL+2,TwoSL,IL);
      
      if (dimLdown>0){
   
         int NR = denS->gNR(ikappa);
         int TwoSR = denS->gTwoSR(ikappa);
         int IR = denS->gIR(ikappa);
   
         int dimLup = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
         int dimR   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
         
         int memSkappa = denS->gKappa(NL+2,TwoSL,IL, denS->gN1(ikappa), 0, denS->gTwoJ(ikappa), NR,TwoSR,IR);
         
         if (memSkappa!=-1){
         
            double * BlockA = Atensor->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
            char notrans = 'N';
            double alpha = sqrt(2.0);
            double beta = 1.0;
            
            dgemm_(&notrans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,BlockA,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
         
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2dall(const int ikappa, double * memS, double * memHeff, const Sobject * denS) const{

   const int N1 = denS->gN1(ikappa);
   const int N2 = denS->gN2(ikappa);
   const int theindex = denS->gIndex();
   int size = denBK->gCurrentDim(theindex,  denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa))
            * denBK->gCurrentDim(theindex+2,denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
   int inc = 1;
   
   if ((N1==2) && (N2==0)){ //2d1
   
      int memSkappa = denS->gKappa(denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa), 0, 2, 0, denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      if (memSkappa!=-1){
         double factor = Prob->gMxElement(theindex, theindex, theindex+1, theindex+1);
         daxpy_(&size,&factor,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   }
   
   if ((N1==0) && (N2==2)){ //2d2
   
      int memSkappa = denS->gKappa(denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa), 2, 0, 0, denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      if (memSkappa!=-1){
         double factor = Prob->gMxElement(theindex, theindex, theindex+1, theindex+1);
         daxpy_(&size,&factor,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   }
   
   if ((N1==2) && (N2==2)){ //2d3a
   
      double factor = 4 * Prob->gMxElement(theindex, theindex+1, theindex, theindex+1)
                    - 2 * Prob->gMxElement(theindex, theindex+1, theindex+1, theindex);
      daxpy_(&size,&factor,memS+denS->gKappa2index(ikappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
   
   }
   
   if ((N1==1) && (N2==1)){ //2d3b
   
      int fase = (denS->gTwoJ(ikappa) == 0)? 1: -1;
      double factor = Prob->gMxElement(theindex, theindex+1, theindex, theindex+1)
             + fase * Prob->gMxElement(theindex, theindex+1, theindex+1, theindex);
      daxpy_(&size,&factor,memS+denS->gKappa2index(ikappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
   
   }
   
   if ((N1==2) && (N2==1)){ //2d3c
   
      double factor = 2 * Prob->gMxElement(theindex, theindex+1, theindex, theindex+1)
                        - Prob->gMxElement(theindex, theindex+1, theindex+1, theindex);
      daxpy_(&size,&factor,memS+denS->gKappa2index(ikappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
   
   }
   
   if ((N1==1) && (N2==2)){ //2d3d
   
      double factor = 2 * Prob->gMxElement(theindex, theindex+1, theindex, theindex+1)
                        - Prob->gMxElement(theindex, theindex+1, theindex+1, theindex);
      daxpy_(&size,&factor,memS+denS->gKappa2index(ikappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
   
   }
   
}

void CheMPS2::Heff::addDiagram2e1and2e2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==2){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
      int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IR);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,  TwoSR,IR);
      
      int memSkappa = denS->gKappa(NL,TwoSL,IL, 0, denS->gN2(ikappa), denS->gTwoJ(ikappa), NR-2,TwoSR,IR);
      
      if (memSkappa!=-1){
         
         double * BlockA = Atensor->gStorage(NR-2,TwoSR,IR,NR,TwoSR,IR);
         char notrans = 'N';
         double alpha = sqrt(2.0);
         double beta = 1.0;
            
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,BlockA,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

      }
   }
   
   if (N1==0){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IR);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,  TwoSR,IR);
      
      int memSkappa = denS->gKappa(NL,TwoSL,IL, 2, denS->gN2(ikappa), denS->gTwoJ(ikappa), NR+2,TwoSR,IR);
      
      if (memSkappa!=-1){
         
         double * BlockA = Atensor->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IR);
         char notrans = 'N';
         char trans = 'T';
         double alpha = sqrt(2.0);
         double beta = 1.0;
            
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,BlockA,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

      }
   }
   
}

void CheMPS2::Heff::addDiagram2f1and2f2(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Atensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==2){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
      int dimRdown = denBK->gCurrentDim(theindex+2,NR-2,TwoSR,IR);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,  TwoSR,IR);
      
      int memSkappa = denS->gKappa(NL,TwoSL,IL, denS->gN1(ikappa), 0, denS->gTwoJ(ikappa), NR-2,TwoSR,IR);
      
      if (memSkappa!=-1){
         
         double * BlockA = Atensor->gStorage(NR-2,TwoSR,IR,NR,TwoSR,IR);
         char notrans = 'N';
         double alpha = sqrt(2.0);
         double beta = 1.0;
            
         dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,BlockA,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

      }
   }
   
   if (N2==0){

      int theindex = denS->gIndex();
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
      int dimRdown = denBK->gCurrentDim(theindex+2,NR+2,TwoSR,IR);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,  TwoSR,IR);
      
      int memSkappa = denS->gKappa(NL,TwoSL,IL, denS->gN1(ikappa), 2, denS->gTwoJ(ikappa), NR+2,TwoSR,IR);
      
      if (memSkappa!=-1){
         
         double * BlockA = Atensor->gStorage(NR,TwoSR,IR,NR+2,TwoSR,IR);
         char notrans = 'N';
         char trans = 'T';
         double alpha = sqrt(2.0);
         double beta = 1.0;
            
         dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,BlockA,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

      }
   }
   
}

void CheMPS2::Heff::addDiagram2b3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1!=0){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,               TwoSL,               IL);
      int dimR     = denBK->gCurrentDim(theindex+2,denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      double * Cblock = Ctensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);

      char trans = 'T';
      char notrans = 'N';
      double alpha = ((N1==2)?1.0:0.5)*sqrt(2.0);
      double beta = 1.0;
            
      dgemm_(&trans,&notrans,&dimL,&dimR,&dimL,&alpha,Cblock,&dimL,memS+denS->gKappa2index(ikappa),&dimL,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

   }
   
}

void CheMPS2::Heff::addDiagram2c3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2!=0){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,               TwoSL,               IL);
      int dimR     = denBK->gCurrentDim(theindex+2,denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      double * Cblock = Ctensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);

      char trans = 'T';
      char notrans = 'N';
      double alpha = ((N2==2)?1.0:0.5)*sqrt(2.0);
      double beta = 1.0;
            
      dgemm_(&trans,&notrans,&dimL,&dimR,&dimL,&alpha,Cblock,&dimL,memS+denS->gKappa2index(ikappa),&dimL,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

   }
   
}

void CheMPS2::Heff::addDiagram2e3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1!=0){

      int theindex = denS->gIndex();
      
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimR     = denBK->gCurrentDim(theindex+2,NR,               TwoSR,               IR);
      int dimL     = denBK->gCurrentDim(theindex  ,denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa));

      double * Cblock = Ctensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);

      char notrans = 'N';
      double alpha = ((N1==2)?1.0:0.5)*sqrt(2.0);
      double beta = 1.0;
            
      dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,memS+denS->gKappa2index(ikappa),&dimL,Cblock,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

   }
   
}

void CheMPS2::Heff::addDiagram2f3spin0(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Ctensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2!=0){

      int theindex = denS->gIndex();
      
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimR     = denBK->gCurrentDim(theindex+2,NR,               TwoSR,               IR);
      int dimL     = denBK->gCurrentDim(theindex  ,denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa));
      
      double * Cblock = Ctensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
      
      char notrans = 'N';
      double alpha = ((N2==2)?1.0:0.5)*sqrt(2.0);
      double beta = 1.0;
            
      dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,memS+denS->gKappa2index(ikappa),&dimL,Cblock,&dimR,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);

   }
   
}

void CheMPS2::Heff::addDiagram2b3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==1){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N2 = denS->gN2(ikappa);
      
      int dimLup   = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL, TwoSLdown, IL);
         
         if (dimLdown>0){
         
            double * Dblock = Dtensor->gStorage(NL,TwoSLdown,IL,NL,TwoSL,IL);

            int TwoS2 = (N2==1)?1:0;
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
         
                  int memSkappa = denS->gKappa(NL, TwoSLdown, IL, N1, N2, TwoJdown, NR, TwoSR, IR);
      
                  if (memSkappa!=-1){
            
                     int fase = phase(TwoSLdown + TwoSR + TwoJ + TwoS2 + TwoJdown - 1);
                     double alpha = fase * sqrt(3.0*(TwoJ+1)*(TwoJdown+1)*(TwoSL+1)) * Wigner::wigner6j(TwoJdown,TwoJ,2,1,1,TwoS2) * Wigner::wigner6j(TwoJdown,TwoJ,2,TwoSL,TwoSLdown,TwoSR);
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 1.0;
               
                     dgemm_(&trans,&notra,&dimLup,&dimR,&dimLdown,&alpha,Dblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2c3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==1){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N1 = denS->gN1(ikappa);
      
      int dimLup   = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      for (int TwoSLdown=TwoSL-2; TwoSLdown<=TwoSL+2; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL, TwoSLdown, IL);
         
         if (dimLdown>0){
         
            double * Dblock = Dtensor->gStorage(NL,TwoSLdown,IL,NL,TwoSL,IL);
            
            int TwoS1 = (N1==1)?1:0;
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
         
                  int memSkappa = denS->gKappa(NL, TwoSLdown, IL, N1, N2, TwoJdown, NR, TwoSR, IR);
      
                  if (memSkappa!=-1){
            
                     int fase = phase(TwoSLdown + TwoSR + 2*TwoJ + TwoS1 - 1);
                     double alpha = fase * sqrt(3.0*(TwoJ+1)*(TwoJdown+1)*(TwoSL+1)) * Wigner::wigner6j(TwoJdown,TwoJ,2,1,1,TwoS1) * Wigner::wigner6j(TwoJdown,TwoJ,2,TwoSL,TwoSLdown,TwoSR);
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 1.0;
               
                     dgemm_(&trans,&notra,&dimLup,&dimR,&dimLdown,&alpha,Dblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                     
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2e3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==1){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N2 = denS->gN2(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR, TwoSRdown, IR);
         
         if (dimRdown>0){

            double * Dblock = Dtensor->gStorage(NR,TwoSRdown,IR,NR,TwoSR,IR);
            
            int TwoS2 = (N2==1)?1:0;
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
         
                  int memSkappa = denS->gKappa(NL, TwoSL, IL, N1, N2, TwoJdown, NR, TwoSRdown, IR);
      
                  if (memSkappa!=-1){
            
                     int fase = phase(TwoSRdown + TwoSL + 2*TwoJ + TwoS2 + 1);
                     double alpha = fase * sqrt(3.0*(TwoJ+1)*(TwoJdown+1)*(TwoSRdown+1)) * Wigner::wigner6j(TwoJdown,TwoJ,2,1,1,TwoS2) * Wigner::wigner6j(TwoJdown,TwoJ,2,TwoSR,TwoSRdown,TwoSL);
                     char notr = 'N';
                     double beta = 1.0;
               
                     dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Dblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                     
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram2f3spin1(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorOperator * Dtensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==1){

      int theindex = denS->gIndex();
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N1 = denS->gN1(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimRup   = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      for (int TwoSRdown=TwoSR-2; TwoSRdown<=TwoSR+2; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR, TwoSRdown, IR);
         
         if (dimRdown>0){
         
            double * Dblock = Dtensor->gStorage(NR,TwoSRdown,IR,NR,TwoSR,IR);
            
            int TwoS1 = (N1==1)?1:0;
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
               
                  int memSkappa = denS->gKappa(NL, TwoSL, IL, N1, N2, TwoJdown, NR, TwoSRdown, IR);
      
                  if (memSkappa!=-1){
            
                     int fase = phase(TwoSRdown + TwoSL + TwoJ + TwoS1 + TwoJdown + 1);
                     double alpha = fase * sqrt(3.0*(TwoJ+1)*(TwoJdown+1)*(TwoSRdown+1)) * Wigner::wigner6j(TwoJdown,TwoJ,2,1,1,TwoS1) * Wigner::wigner6j(TwoJdown,TwoJ,2,TwoSR,TwoSRdown,TwoSL);
                     char notr = 'N';
                     double beta = 1.0;
               
                     dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimL,Dblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                     
                  }
               }
            }
         }
      }
   }
   
}



