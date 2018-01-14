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

void CheMPS2::Heff::addDiagram1A(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorX * Xleft) const{
   int dimL = denBK->gCurrentDim(denS->gIndex(), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa));
   int dimR = denBK->gCurrentDim(denS->gIndex()+2, denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa));
   double * BlockX = Xleft->gStorage( denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa) );
   
   double one = 1.0;
   char notr = 'N';
   dgemm_(&notr,&notr,&dimL,&dimR,&dimL,&one,BlockX,&dimL,memS+denS->gKappa2index(ikappa),&dimL,&one,memHeff+denS->gKappa2index(ikappa),&dimL);
}

void CheMPS2::Heff::addDiagram1B(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorX * Xright) const{
   int dimL = denBK->gCurrentDim(denS->gIndex(), denS->gNL(ikappa), denS->gTwoSL(ikappa), denS->gIL(ikappa));
   int dimR = denBK->gCurrentDim(denS->gIndex()+2, denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa));
   double * BlockX = Xright->gStorage( denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa), denS->gNR(ikappa), denS->gTwoSR(ikappa), denS->gIR(ikappa) );
   
   double one = 1.0;
   char notr = 'N';
   char trans = 'T';
   dgemm_(&notr,&trans,&dimL,&dimR,&dimR,&one,memS+denS->gKappa2index(ikappa),&dimL,BlockX,&dimR,&one,memHeff+denS->gKappa2index(ikappa),&dimL);
}

void CheMPS2::Heff::addDiagram1C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, double Helem_links) const{
   if (denS->gN1(ikappa)==2){
      int inc = 1;
      int ptr = denS->gKappa2index(ikappa);
      int dim = denS->gKappa2index(ikappa+1) - ptr;
      daxpy_(&dim,&Helem_links,memS+ptr,&inc,memHeff+ptr,&inc);
   }
}

void CheMPS2::Heff::addDiagram1D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, double Helem_rechts) const{
   if (denS->gN2(ikappa)==2){
      int inc = 1;
      int ptr = denS->gKappa2index(ikappa);
      int dim = denS->gKappa2index(ikappa+1) - ptr;
      daxpy_(&dim,&Helem_rechts,memS+ptr,&inc,memHeff+ptr,&inc);
   }
}

void CheMPS2::Heff::addDiagramExcitations(const int ikappa, double * memS, double * memHeff, const Sobject * denS, int nLower, double ** VeffTilde) const{

   int dimTotal = denS->gKappa2index(denS->gNKappa());
   int ptr = denS->gKappa2index(ikappa);
   int dimBlock = denS->gKappa2index(ikappa+1) - ptr;
   int inc = 1;
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   for (int state=0; state<nLower; state++){
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_specific_excitation( Prob->gL(), state ) == MPIRANK )
      #endif
      {
         double alpha = ddot_(&dimTotal, memS, &inc, VeffTilde[state], &inc);
         daxpy_(&dimBlock,&alpha,VeffTilde[state]+ptr,&inc,memHeff+ptr,&inc);
      }
   }

}

