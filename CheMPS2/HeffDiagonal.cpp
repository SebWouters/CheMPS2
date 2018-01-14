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

#include "Heff.h"
#include "Lapack.h"
#include "MPIchemps2.h"
#include "Wigner.h"

void CheMPS2::Heff::addDiagonal1A(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorX * Xleft) const{
   int dimL = denBK->gCurrentDim(denS->gIndex(), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa));
   int dimR = denBK->gCurrentDim(denS->gIndex()+2, denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa));
   double * BlockX = Xleft->gStorage( denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa) );
   int ptr = denS->gKappa2index(ikappa);
   
   for (int cnt=0; cnt<dimL; cnt++){
      for (int cnt2=0; cnt2<dimR; cnt2++){
         memHeffDiag[ptr + cnt + dimL*cnt2] += BlockX[cnt*(dimL+1)];
      }
   }

}

void CheMPS2::Heff::addDiagonal1B(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorX * Xright) const{
   int dimL = denBK->gCurrentDim(denS->gIndex(), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa));
   int dimR = denBK->gCurrentDim(denS->gIndex()+2, denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa));
   double * BlockX = Xright->gStorage( denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa), denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa) );
   int ptr = denS->gKappa2index(ikappa);
   
   for (int cnt=0; cnt<dimL; cnt++){
      for (int cnt2=0; cnt2<dimR; cnt2++){
         memHeffDiag[ptr + cnt + dimL*cnt2] += BlockX[cnt2*(dimR+1)];
      }
   }
   
}

void CheMPS2::Heff::addDiagonal1C(const int ikappa, double * memHeffDiag, const Sobject * denS, const double Helem_links) const{
   if (denS->gN1(ikappa)==2){
      int ptr = denS->gKappa2index(ikappa);
      int dim = denS->gKappa2index(ikappa+1) - ptr;
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += Helem_links; }
   }
}

void CheMPS2::Heff::addDiagonal1D(const int ikappa, double * memHeffDiag, const Sobject * denS, const double Helem_rechts) const{
   if (denS->gN2(ikappa)==2){
      int ptr = denS->gKappa2index(ikappa);
      int dim = denS->gKappa2index(ikappa+1) - ptr;
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += Helem_rechts; }
   }
}

void CheMPS2::Heff::addDiagonal2d3all(const int ikappa, double * memHeffDiag, const Sobject * denS) const{

   if ((denS->gN1(ikappa)==2)&&(denS->gN2(ikappa)==2)){ //2d3a
   
      const int theindex = denS->gIndex();
      const int ptr = denS->gKappa2index(ikappa);
      const int dim = denS->gKappa2index(ikappa+1) - ptr;
      const double factor = 4 * Prob->gMxElement(theindex,theindex+1,theindex,theindex+1)
                          - 2 * Prob->gMxElement(theindex,theindex+1,theindex+1,theindex);
      
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += factor; }
      
   }
   
   if ((denS->gN1(ikappa)==1)&&(denS->gN2(ikappa)==1)){ //2d3b
   
      const int theindex = denS->gIndex();
      const int ptr = denS->gKappa2index(ikappa);
      const int dim = denS->gKappa2index(ikappa+1) - ptr;
      const int fase = (denS->gTwoJ(ikappa) == 2)? -1: 1;
      const double factor = Prob->gMxElement(theindex,theindex+1,theindex,theindex+1)
                   + fase * Prob->gMxElement(theindex,theindex+1,theindex+1,theindex);
      
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += factor; }
      
   }
   
   if ((denS->gN1(ikappa)==2)&&(denS->gN2(ikappa)==1)){ //2d3c
   
      const int theindex = denS->gIndex();
      const int ptr = denS->gKappa2index(ikappa);
      const int dim = denS->gKappa2index(ikappa+1) - ptr;
      const double factor = 2 * Prob->gMxElement(theindex,theindex+1,theindex,theindex+1)
                              - Prob->gMxElement(theindex,theindex+1,theindex+1,theindex);
      
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += factor; }
      
   }
   
   if ((denS->gN1(ikappa)==1)&&(denS->gN2(ikappa)==2)){ //2d3d
   
      const int theindex = denS->gIndex();
      const int ptr = denS->gKappa2index(ikappa);
      const int dim = denS->gKappa2index(ikappa+1) - ptr;
      const double factor = 2 * Prob->gMxElement(theindex,theindex+1,theindex,theindex+1)
                              - Prob->gMxElement(theindex,theindex+1,theindex+1,theindex);
      
      for (int cnt=0; cnt<dim; cnt++){ memHeffDiag[ptr + cnt] += factor; }
      
   }
   
}

void CheMPS2::Heff::addDiagonal2b3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1!=0){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      double sqrt0p5 = sqrt(0.5);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,               TwoSL,               IL);
      int dimR     = denBK->gCurrentDim(theindex+2,denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      double * Cblock = Ctensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += N1*sqrt0p5*Cblock[(dimL+1)*cntL];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2c3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2!=0){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      double sqrt0p5 = sqrt(0.5);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,               TwoSL,               IL);
      int dimR     = denBK->gCurrentDim(theindex+2,denS->gNR(ikappa),denS->gTwoSR(ikappa),denS->gIR(ikappa));
      
      double * Cblock = Ctensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += N2*sqrt0p5*Cblock[(dimL+1)*cntL];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2e3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1!=0){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      double sqrt0p5 = sqrt(0.5);
      
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimR     = denBK->gCurrentDim(theindex+2,NR,               TwoSR,               IR);
      int dimL     = denBK->gCurrentDim(theindex  ,denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa));
      
      double * Cblock = Ctensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += N1*sqrt0p5*Cblock[(dimR+1)*cntR];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2f3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Ctensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2!=0){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      double sqrt0p5 = sqrt(0.5);
      
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      
      int dimR     = denBK->gCurrentDim(theindex+2,NR,               TwoSR,               IR);
      int dimL     = denBK->gCurrentDim(theindex  ,denS->gNL(ikappa),denS->gTwoSL(ikappa),denS->gIL(ikappa));
      
      double * Cblock = Ctensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += N2*sqrt0p5*Cblock[(dimR+1)*cntR];
         }
      }
      
   }
   
}

void CheMPS2::Heff::addDiagonal2b3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==1){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N2 = denS->gN2(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      int fase = phase(TwoSL + TwoSR + 2*TwoJ + ((N2==1)?1:0) - 1);
      const double alpha = fase * (TwoJ+1) * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(TwoJ,TwoJ,2,1,1,((N2==1)?1:0)) * Wigner::wigner6j(TwoJ,TwoJ,2,TwoSL,TwoSL,TwoSR);
      
      double * Dblock = Dtensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimL+1)*cntL];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2c3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==1){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N1 = denS->gN1(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      int fase = phase(TwoSL + TwoSR + 2*TwoJ + ((N1==1)?1:0) - 1);
      const double alpha = fase * (TwoJ+1) * sqrt(3.0*(TwoSL+1)) * Wigner::wigner6j(TwoJ,TwoJ,2,1,1,((N1==1)?1:0)) * Wigner::wigner6j(TwoJ,TwoJ,2,TwoSL,TwoSL,TwoSR);
      
      double * Dblock = Dtensor->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimL+1)*cntL];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2e3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const{

   int N1 = denS->gN1(ikappa);
   
   if (N1==1){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N2 = denS->gN2(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      int fase = phase(TwoSR + TwoSL + 2*TwoJ + ((N2==1)?1:0) + 1);
      const double alpha = fase * (TwoJ+1) * sqrt(3.0*(TwoSR+1)) * Wigner::wigner6j(TwoJ,TwoJ,2,1,1,((N2==1)?1:0)) * Wigner::wigner6j(TwoJ,TwoJ,2,TwoSR,TwoSR,TwoSL);
      
      double * Dblock = Dtensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimR+1)*cntR];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2f3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator * Dtensor) const{

   int N2 = denS->gN2(ikappa);
   
   if (N2==1){

      int theindex = denS->gIndex();
      int ptr = denS->gKappa2index(ikappa);
      
      int NL = denS->gNL(ikappa);
      int TwoSL = denS->gTwoSL(ikappa);
      int IL = denS->gIL(ikappa);
      int NR = denS->gNR(ikappa);
      int TwoSR = denS->gTwoSR(ikappa);
      int IR = denS->gIR(ikappa);
      int TwoJ = denS->gTwoJ(ikappa);
      int N1 = denS->gN1(ikappa);
      
      int dimL     = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
      int dimR     = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
      
      int fase = phase(TwoSR + TwoSL + 2*TwoJ + ((N1==1)?1:0) + 1);
      const double alpha = fase * (TwoJ+1) * sqrt(3.0*(TwoSR+1)) * Wigner::wigner6j(TwoJ,TwoJ,2,1,1,((N1==1)?1:0)) * Wigner::wigner6j(TwoJ,TwoJ,2,TwoSR,TwoSR,TwoSL);
      
      double * Dblock = Dtensor->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
      for (int cntR=0; cntR<dimR; cntR++){
         for (int cntL=0; cntL<dimL; cntL++){
            memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimR+1)*cntR];
         }
      }

   }
   
}

void CheMPS2::Heff::addDiagonal2a3spin0(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator **** Ctensors, TensorF0 **** F0tensors) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();
   int ptr = denS->gKappa2index(ikappa);
   
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int l_gamma=0; l_gamma<theindex; l_gamma++){
         for (int l_alpha=l_gamma+1; l_alpha<theindex; l_alpha++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_gamma, l_alpha ) == MPIRANK )
            #endif
            {
               if (denBK->gIrrep(l_alpha) == denBK->gIrrep(l_gamma)){
            
                  double * Cblock = Ctensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
                  double * BlockF0 = F0tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += BlockF0[(dimL+1)*cntL] * Cblock[(dimR+1)*cntR];
                     }
                  }
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
               if (denBK->gIrrep(l_alpha) == denBK->gIrrep(l_gamma)){
            
                  double * Cblock = Ctensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
                  double * BlockF0 = F0tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += BlockF0[(dimL+1)*cntL] * Cblock[(dimR+1)*cntR];
                     }
                  }
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
               if (denBK->gIrrep(l_delta) == denBK->gIrrep(l_beta)){
            
                  double * Cblock = Ctensors[theindex-1][l_beta-l_delta][l_delta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
                  double * BlockF0 = F0tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += Cblock[(dimL+1)*cntL] * BlockF0[(dimR+1)*cntR];
                     }
                  }
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
               if (denBK->gIrrep(l_delta) == denBK->gIrrep(l_beta)){
            
                  double * Cblock = Ctensors[theindex-1][l_delta-l_beta][l_beta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
                  double * BlockF0 = F0tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += Cblock[(dimL+1)*cntL] * BlockF0[(dimR+1)*cntR];
                     }
                  }
               }
            }
         }
      }
   
   }
   
}

void CheMPS2::Heff::addDiagonal2a3spin1(const int ikappa, double * memHeffDiag, const Sobject * denS, TensorOperator **** Dtensors, TensorF1 **** F1tensors) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int TwoJ = denS->gTwoJ(ikappa);
   const int fase = phase(TwoSL+TwoSR+TwoJ+2);
   const double alpha = fase * sqrt((TwoSR + 1)*(TwoSL + 1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,TwoSR,TwoSL,2);
   
   int theindex = denS->gIndex();
   int ptr = denS->gKappa2index(ikappa);
   
   int dimL = denBK->gCurrentDim(theindex  ,NL,TwoSL,IL);
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   bool leftSum = ( theindex < Prob->gL()*0.5 )?true:false;
   
   if (leftSum){
      
      for (int l_gamma=0; l_gamma<theindex; l_gamma++){
         for (int l_alpha=l_gamma+1; l_alpha<theindex; l_alpha++){
         
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), l_gamma, l_alpha ) == MPIRANK )
            #endif
            {
               if (denBK->gIrrep(l_alpha) == denBK->gIrrep(l_gamma)){
            
                  double * Dblock = Dtensors[theindex+1][l_alpha-l_gamma][theindex+1-l_alpha]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
                  double * BlockF1 = F1tensors[theindex-1][l_alpha-l_gamma][theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += alpha * BlockF1[(dimL+1)*cntL] * Dblock[(dimR+1)*cntR];
                     }
                  }
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
         
              if (denBK->gIrrep(l_alpha) == denBK->gIrrep(l_gamma)){
            
                 double * Dblock = Dtensors[theindex+1][l_gamma-l_alpha][theindex+1-l_gamma]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
                 double * BlockF1 = F1tensors[theindex-1][l_gamma-l_alpha][theindex-1-l_gamma]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += alpha * BlockF1[(dimL+1)*cntL] * Dblock[(dimR+1)*cntR];
                     }
                  }
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
               if (denBK->gIrrep(l_delta) == denBK->gIrrep(l_beta)){
            
                  double * Dblock = Dtensors[theindex-1][l_beta-l_delta][l_delta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
                  double * BlockF1 = F1tensors[theindex+1][l_beta-l_delta][l_delta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimL+1)*cntL] * BlockF1[(dimR+1)*cntR];
                     }
                  }
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
               if (denBK->gIrrep(l_delta) == denBK->gIrrep(l_beta)){
            
                  double * Dblock = Dtensors[theindex-1][l_delta-l_beta][l_beta-theindex]->gStorage(NL,TwoSL,IL,NL,TwoSL,IL);
                  double * BlockF1 = F1tensors[theindex+1][l_delta-l_beta][l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               
                  for (int cntL=0; cntL<dimL; cntL++){
                     for (int cntR=0; cntR<dimR; cntR++){
                        memHeffDiag[ptr + cntL + dimL*cntR] += alpha * Dblock[(dimL+1)*cntL] * BlockF1[(dimR+1)*cntR];
                     }
                  }
               }
            }
         }
      }
   
   }
   
}

void CheMPS2::Heff::addDiagonalExcitations(const int ikappa, double * memHeffDiag, const Sobject * denS, int nLower, double ** VeffTilde) const{

   const int loc = denS->gKappa2index(ikappa);
   const int dim = denS->gKappa2index(ikappa+1) - loc;
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   for (int state=0; state<nLower; state++){
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_specific_excitation( Prob->gL(), state ) == MPIRANK )
      #endif
      {
         for (int cnt=0; cnt<dim; cnt++){
            memHeffDiag[loc + cnt] += VeffTilde[state][loc+cnt] * VeffTilde[state][loc+cnt];
         }
      }
   }

}


