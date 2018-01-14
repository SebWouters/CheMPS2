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

void CheMPS2::Heff::addDiagram3Aand3D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qleft, TensorL ** Lleft, double * temp) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();
   int ILdown = Irreps::directProd(IL,denBK->gIrrep(theindex));
   int TwoS2 = (N2==1)?1:0;
   
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   int dimLup = denBK->gCurrentDim(theindex,NL,TwoSL,IL);

   if (N1==2){ //3A1A and 3D1
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL+1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSR,IR);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSR+2+TwoS2);
                     double factor = sqrt((TwoJdown+1)*(TwoSLdown+1.0))*fase*Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
                     double beta = 1.0; //add
                     char notr = 'N';
                  
                     double * BlockQ = Qleft->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     int inc = 1;
                     int size = dimLup * dimLdown;
                     dcopy_(&size, BlockQ, &inc, temp, &inc);
                  
                     for (int l_index=0; l_index<theindex; l_index++){
                        if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                           double alpha = Prob->gMxElement(l_index,theindex,theindex,theindex);
                           double * BlockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                           daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                        }
                     }
                  
                     dgemm_(&notr,&notr,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   if (N1==1){ //3A1B
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
            int dimLdown = denBK->gCurrentDim(theindex, NL+1, TwoSLdown, ILdown);
            int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,0,N2,TwoS2,NR,TwoSR,IR);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSR+1+TwoS2);
               double factor = sqrt((TwoSLdown+1)*(TwoJ+1.0))*fase*Wigner::wigner6j(TwoS2,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
               double beta = 1.0;
               char notr = 'N';
               double * BlockQ = Qleft->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
               dgemm_(&notr,&notr,&dimLup,&dimR,&dimLdown,&factor,BlockQ,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   if (N1==0){ //3A2A
      for (int TwoSLdown=TwoSL-1;TwoSLdown<=TwoSL+1;TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL-1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,1,N2,TwoJdown,NR,TwoSR,IR);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSLdown+TwoSR+1+TwoS2);
                     double factor = fase*sqrt((TwoSL+1)*(TwoJdown+1.0))*Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSL,TwoSLdown,TwoSR);
                     double beta = 1.0;
                     char notr = 'N';
                     char trans = 'T';
                     double * BlockQ = Qleft->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     dgemm_(&trans,&notr,&dimLup,&dimR,&dimLdown,&factor,BlockQ,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   if (N1==1){ //3A2B ans 3D2
      for (int TwoSLdown=TwoSL-1;TwoSLdown<=TwoSL+1;TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS2) && (TwoSLdown>=0)){
            int dimLdown = denBK->gCurrentDim(theindex, NL-1, TwoSLdown, ILdown);
            int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,2,N2,TwoS2,NR,TwoSR,IR);
            if (memSkappa!=-1){
               int fase = phase(TwoSLdown+TwoSR+2+TwoS2);
               double factor = fase*sqrt((TwoSL+1)*(TwoJ+1.0))*Wigner::wigner6j(TwoS2,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
               double beta = 1.0;
               char notr = 'N';
               char trans = 'T';
            
               double * BlockQ = Qleft->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
               int inc = 1;
               int size = dimLup * dimLdown;
               dcopy_(&size, BlockQ, &inc, temp, &inc);
               
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     double alpha = Prob->gMxElement(l_index,theindex,theindex,theindex);
                     double * BlockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                  }
               }
               
               dgemm_(&trans,&notr,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram3Band3I(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qleft, TensorL ** Lleft, double * temp) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();
   int ILdown = Irreps::directProd(IL,denBK->gIrrep(theindex+1));
   int TwoS1 = (N1==1)?1:0;
   
   int dimR = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   int dimLup = denBK->gCurrentDim(theindex,NL,TwoSL,IL);

   if (N2==2){ //3B1A and 3I2
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL+1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSR,IR);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSR+3-TwoJdown);
                     double factor = sqrt((TwoJdown+1)*(TwoSLdown+1.0))*fase*Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
                     double beta = 1.0; //add
                     char notr = 'N';
                  
                     double * BlockQ = Qleft->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     int inc = 1;
                     int size = dimLup * dimLdown;
                     dcopy_(&size, BlockQ, &inc, temp, &inc);
                  
                     for (int l_index=0; l_index<theindex; l_index++){
                        if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                           double alpha = Prob->gMxElement(l_index,theindex+1,theindex+1,theindex+1);
                           double * BlockL = Lleft[theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                           daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                        }
                     }
                  
                     dgemm_(&notr,&notr,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   if (N2==1){ //3B1B
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
            int dimLdown = denBK->gCurrentDim(theindex, NL+1, TwoSLdown, ILdown);
            int memSkappa = denS->gKappa(NL+1,TwoSLdown,ILdown,N1,0,TwoS1,NR,TwoSR,IR);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSR+2-TwoJ);
               double factor = sqrt((TwoSLdown+1)*(TwoJ+1.0))*fase*Wigner::wigner6j(TwoS1,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
               double beta = 1.0;
               char notr = 'N';
               double * BlockQ = Qleft->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
               dgemm_(&notr,&notr,&dimLup,&dimR,&dimLdown,&factor,BlockQ,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
   if (N2==0){ //3B2A
      for (int TwoSLdown=TwoSL-1;TwoSLdown<=TwoSL+1;TwoSLdown+=2){
      
         int dimLdown = denBK->gCurrentDim(theindex, NL-1, TwoSLdown, ILdown);
         if (dimLdown>0){
         
            int TwoJstart = ((TwoSR!=TwoSLdown) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSLdown-TwoSR)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,1,TwoJdown,NR,TwoSR,IR);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSLdown+TwoSR+2-TwoJdown);
                     double factor = fase*sqrt((TwoSL+1)*(TwoJdown+1.0))*Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSL,TwoSLdown,TwoSR);
                     double beta = 1.0;
                     char notr = 'N';
                     char trans = 'T';
                     double * BlockQ = Qleft->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     dgemm_(&trans,&notr,&dimLup,&dimR,&dimLdown,&factor,BlockQ,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  }
               }
            }
         }
      }
   }
   
   if (N2==1){ //3B2B and 3I1
      for (int TwoSLdown=TwoSL-1;TwoSLdown<=TwoSL+1;TwoSLdown+=2){
         if ((abs(TwoSLdown-TwoSR)<=TwoS1) && (TwoSLdown>=0)){
            int dimLdown = denBK->gCurrentDim(theindex, NL-1, TwoSLdown, ILdown);
            int memSkappa = denS->gKappa(NL-1,TwoSLdown,ILdown,N1,2,TwoS1,NR,TwoSR,IR);
            if (memSkappa!=-1){
               int fase = phase(TwoSLdown+TwoSR+3-TwoJ);
               double factor = fase*sqrt((TwoSL+1)*(TwoJ+1.0))*Wigner::wigner6j(TwoS1,TwoJ,1,TwoSL,TwoSLdown,TwoSR);
               double beta = 1.0;
               char notr = 'N';
               char trans = 'T';
            
               double * BlockQ = Qleft->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
               int inc = 1;
               int size = dimLup * dimLdown;
               dcopy_(&size, BlockQ, &inc, temp, &inc);
               
               for (int l_index=0; l_index<theindex; l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(l_index,theindex+1,theindex+1,theindex+1);
                     double * BlockL = Lleft[theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                  }
               }
            
               dgemm_(&trans,&notr,&dimLup,&dimR,&dimLdown,&factor,temp,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram3C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ ** Qleft, TensorL ** Lright, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();

   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   int dimLup = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
   
   //First do 3C1
   for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
      
            int fase = phase(TwoSLdown+TwoSR+TwoJ+1 + ((N1==1)?2:0) + ((N2==1)?2:0) );
            const double factor = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,TwoSRdown,TwoSLdown,1);
      
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_q( Prob->gL(), l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL,denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR,denBK->gIrrep(l_index));
                  int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, N1, N2, TwoJ, NR+1, TwoSRdown, IRdown);
                  if (memSkappa!=-1){
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  
                     double * Qblock = Qleft[ l_index-theindex  ]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                  
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 0.0; //set
                     double alpha = factor;
                     dgemm_(&notra,&notra,&dimLup,&dimRdown,&dimLdown,&alpha,Qblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,temp,&dimLup);
                  
                     beta = 1.0; //add
                     alpha = 1.0;
                     dgemm_(&notra,&trans,&dimLup,&dimRup,&dimRdown,&alpha,temp,&dimLup,Lblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }

   //Then do 3C2
   for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
      
            int fase = phase(TwoSL+TwoSRdown+TwoJ+1 + ((N1==1)?2:0) + ((N2==1)?2:0) );
            const double factor = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,TwoSRdown,TwoSLdown,1);
      
            for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_q( Prob->gL(), l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL,denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR,denBK->gIrrep(l_index));
                  int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, N1, N2, TwoJ, NR-1, TwoSRdown, IRdown);
                  if (memSkappa!=-1){
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  
                     double * Qblock = Qleft[ l_index-theindex  ]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     double * Lblock = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                  
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 0.0; //set
                     double alpha = factor;
                     dgemm_(&trans,&notra,&dimLup,&dimRdown,&dimLdown,&alpha,Qblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,temp,&dimLup);
                  
                     beta = 1.0; //add
                     alpha = 1.0;
                     dgemm_(&notra,&notra,&dimLup,&dimRup,&dimRdown,&alpha,temp,&dimLup,Lblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram3Eand3H(const int ikappa, double * memS, double * memHeff, const Sobject * denS) const{ //TwoJ = TwoJdown

   int theindex = denS->gIndex();

   if (denBK->gIrrep(theindex) != denBK->gIrrep(theindex+1)){ return; }

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int size = ( denBK->gCurrentDim(theindex,NL,TwoSL,IL) ) * ( denBK->gCurrentDim(theindex+2,NR,TwoSR,IR) );
   int inc = 1;
   
   if ((N1==2) && (N2==0)){ //3E1A
   
      int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = sqrt(2.0) * Prob->gMxElement(theindex,theindex,theindex,theindex+1);
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   
   }
   
   if ((N1==2) && (N2==1)){ //3E1B and 3H1B
   
      int memSkappa = denS->gKappa(NL,TwoSL,IL,1,2,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = - ( Prob->gMxElement(theindex,theindex,theindex,theindex+1) + Prob->gMxElement(theindex,theindex+1,theindex+1,theindex+1) );
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   
   }
   
   if ((N1==1) && (N2==1) && (TwoJ==0)){ //3E2A and 3H1A
   
      int memSkappa = denS->gKappa(NL,TwoSL,IL,2,0,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = sqrt(2.0) * Prob->gMxElement(theindex,theindex,theindex,theindex+1);
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
      
      memSkappa = denS->gKappa(NL,TwoSL,IL,0,2,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = sqrt(2.0) * Prob->gMxElement(theindex,theindex+1,theindex+1,theindex+1);
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   
   }
   
   if ((N1==1) && (N2==2)){ //3E2B and 3H2B
   
      int memSkappa = denS->gKappa(NL,TwoSL,IL,2,1,1,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = - ( Prob->gMxElement(theindex,theindex,theindex,theindex+1) + Prob->gMxElement(theindex,theindex+1,theindex+1,theindex+1) );
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   
   }
   
   if ((N1==0) && (N2==2)){ //3H2A
   
      int memSkappa = denS->gKappa(NL,TwoSL,IL,1,1,0,NR,TwoSR,IR);
      if (memSkappa!=-1){
         double alpha = sqrt(2.0) * Prob->gMxElement(theindex,theindex+1,theindex+1,theindex+1);
         daxpy_(&size,&alpha,memS+denS->gKappa2index(memSkappa),&inc,memHeff+denS->gKappa2index(ikappa),&inc);
      }
   
   }
   
}

void CheMPS2::Heff::addDiagram3Kand3F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qright, TensorL ** Lright, double * temp) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();
   int IRdown = Irreps::directProd(IR,denBK->gIrrep(theindex));
   int TwoS2 = (N2==1)?1:0;
   
   int dimL   = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);

   if (N1==1){ //3K1A
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
            int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
            int memSkappa = denS->gKappa(NL,TwoSL,IL,0,N2,TwoS2,NR-1,TwoSRdown,IRdown);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSR+TwoJ+2*TwoS2);
               double factor = sqrt((TwoJ+1)*(TwoSR+1.0)) * fase * Wigner::wigner6j(TwoS2,TwoJ,1,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0; //add
               char notr = 'N';
               double * BlockQ = Qright->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
               dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,BlockQ,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            }
         }
      }
   }
   
   if (N1==2){ //3K1B and 3F1
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
         if (dimRdown>0){
         
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,1,N2,TwoJdown,NR-1,TwoSRdown,IRdown);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSR+TwoJdown+1+2*TwoS2);
                     double factor = sqrt((TwoJdown+1)*(TwoSR+1.0)) * fase * Wigner::wigner6j( TwoJdown, TwoS2, 1, TwoSR, TwoSRdown, TwoSL );
                     double beta = 1.0; //add
                     char notr = 'N';
                  
                     double * BlockQ = Qright->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     int inc = 1;
                     int size = dimRup * dimRdown;
                     dcopy_(&size,BlockQ,&inc,temp,&inc);
                  
                     for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                        if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                           double alpha = Prob->gMxElement(theindex,theindex,theindex,l_index);
                           double * BlockL = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                        }
                     }
                  
                     dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  }
               }
            }
         }
      }
   }
   
   if (N1==0){ //3K2A
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
         if (dimRdown>0){
         
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS2==0)) ? 1 + TwoS2 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS2; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,1,N2,TwoJdown,NR+1,TwoSRdown,IRdown);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSRdown+TwoJdown+2*TwoS2);
                     double factor = sqrt((TwoJdown+1)*(TwoSRdown+1.0)) * fase * Wigner::wigner6j(TwoJdown,TwoS2,1,TwoSR,TwoSRdown,TwoSL);
                     double beta = 1.0; //add
                     char notr = 'N';
                     char tran = 'T';
                     double * BlockQ = Qright->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     dgemm_(&notr,&tran,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,BlockQ,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  }
               }
            }
         }
      }
   }
   
   if (N1==1){ //3K2B and 3F2
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS2) && (TwoSRdown>=0)){
            int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
            int memSkappa = denS->gKappa(NL,TwoSL,IL,2,N2,TwoS2,NR+1,TwoSRdown,IRdown);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSRdown+TwoJ+1+2*TwoS2);
               double factor = sqrt((TwoJ+1)*(TwoSRdown+1.0)) * fase * Wigner::wigner6j(TwoS2,TwoJ,1,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0; //add
               char notr = 'N';
               char tran = 'T';
            
               double * BlockQ = Qright->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
               int inc = 1;
               int size = dimRup * dimRdown;
               dcopy_(&size,BlockQ,&inc,temp,&inc);
            
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex)){
                     double alpha = Prob->gMxElement(theindex,theindex,theindex,l_index);
                     double * BlockL = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                  }
               }
            
               dgemm_(&notr,&tran,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram3Land3G(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ * Qright, TensorL ** Lright, double * temp) const{

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();
   int IRdown = Irreps::directProd(IR,denBK->gIrrep(theindex+1));
   int TwoS1 = (N1==1)?1:0;
   
   int dimL   = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);

   if (N2==1){ //3L1A
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
            int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
            int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,0,TwoS1,NR-1,TwoSRdown,IRdown);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSR+TwoS1+1);
               double factor = sqrt((TwoJ+1)*(TwoSR+1.0)) * fase * Wigner::wigner6j(TwoS1,TwoJ,1,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0; //add
               char notr = 'N';
               double * BlockQ = Qright->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
               dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,BlockQ,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            }
         }
      }
   }
   
   if (N2==2){ //3L1B and 3G1
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
         if (dimRdown>0){
         
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,1,TwoJdown,NR-1,TwoSRdown,IRdown);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSR+TwoS1+2);
                     double factor = sqrt((TwoJdown+1)*(TwoSR+1.0)) * fase * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSR,TwoSRdown,TwoSL);
                     double beta = 1.0; //add
                     char notr = 'N';
                  
                     double * BlockQ = Qright->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                     int inc = 1;
                     int size = dimRup * dimRdown;
                     dcopy_(&size,BlockQ,&inc,temp,&inc);
                  
                     for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                        if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                           double alpha = Prob->gMxElement(theindex+1,theindex+1,theindex+1,l_index);
                           double * BlockL = Lright[l_index-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                           daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                        }
                     }
                  
                     dgemm_(&notr,&notr,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  }
               }
            }
         }
      }
   }
   
   if (N2==0){ //3L2A
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
      
         int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
         if (dimRdown>0){
         
            int TwoJstart = ((TwoSRdown!=TwoSL) || (TwoS1==0)) ? 1 + TwoS1 : 0;
            for (int TwoJdown=TwoJstart; TwoJdown<=1+TwoS1; TwoJdown+=2){
               if (abs(TwoSL-TwoSRdown)<=TwoJdown){
            
                  int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,1,TwoJdown,NR+1,TwoSRdown,IRdown);
                  if (memSkappa!=-1){
                     int fase = phase(TwoSL+TwoSRdown+TwoS1+1);
                     double factor = sqrt((TwoJdown+1)*(TwoSRdown+1.0)) * fase * Wigner::wigner6j(TwoJdown,TwoS1,1,TwoSR,TwoSRdown,TwoSL);
                     double beta = 1.0; //add
                     char notr = 'N';
                     char tran = 'T';
                     double * BlockQ = Qright->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     dgemm_(&notr,&tran,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,BlockQ,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
                  }
               }
            }
         }
      }
   }
   
   if (N2==1){ //3L2B and 3G2
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSL-TwoSRdown)<=TwoS1) && (TwoSRdown>=0)){
            int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
            int memSkappa = denS->gKappa(NL,TwoSL,IL,N1,2,TwoS1,NR+1,TwoSRdown,IRdown);
            if (memSkappa!=-1){
               int fase = phase(TwoSL+TwoSRdown+TwoS1+2);
               double factor = sqrt((TwoJ+1)*(TwoSRdown+1.0)) * fase * Wigner::wigner6j(TwoS1,TwoJ,1,TwoSR,TwoSRdown,TwoSL);
               double beta = 1.0; //add
               char notr = 'N';
               char tran = 'T';
            
               double * BlockQ = Qright->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
               int inc = 1;
               int size = dimRup * dimRdown;
               dcopy_(&size,BlockQ,&inc,temp,&inc);
            
               for (int l_index=theindex+2; l_index<Prob->gL(); l_index++){
                  if (denBK->gIrrep(l_index) == denBK->gIrrep(theindex+1)){
                     double alpha = Prob->gMxElement(theindex+1,theindex+1,theindex+1,l_index);
                     double * BlockL = Lright[l_index-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                     daxpy_(&size, &alpha, BlockL, &inc, temp, &inc);
                  }
               }
            
               dgemm_(&notr,&tran,&dimL,&dimRup,&dimRdown,&factor,memS+denS->gKappa2index(memSkappa),&dimL,temp,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimL);
            }
         }
      }
   }
   
}

void CheMPS2::Heff::addDiagram3J(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorQ ** Qright, TensorL ** Lleft, double * temp) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif

   int NL = denS->gNL(ikappa);
   int TwoSL = denS->gTwoSL(ikappa);
   int IL = denS->gIL(ikappa);
   int N1 = denS->gN1(ikappa);
   int N2 = denS->gN2(ikappa);
   int TwoJ = denS->gTwoJ(ikappa);
   int NR = denS->gNR(ikappa);
   int TwoSR = denS->gTwoSR(ikappa);
   int IR = denS->gIR(ikappa);
   
   int theindex = denS->gIndex();

   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   int dimLup = denBK->gCurrentDim(theindex,  NL,TwoSL,IL);
   
   //First do 3J2
   for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
      
            int fase = phase(TwoSLdown+TwoSR+TwoJ+1 + ((N1==1)?2:0) + ((N2==1)?2:0) );
            const double factor = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,TwoSRdown,TwoSLdown,1);
      
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_q( Prob->gL(), l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL,denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR,denBK->gIrrep(l_index));
                  int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, N1, N2, TwoJ, NR+1, TwoSRdown, IRdown);
                  if (memSkappa!=-1){
               
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                  
                     double * Lblock = Lleft[ theindex-1-l_index]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                     double * Qblock = Qright[theindex+1-l_index]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                  
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 0.0; //set
                     double alpha = factor;
                     dgemm_(&notra,&notra,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLup,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,temp,&dimLup);
                  
                     beta = 1.0; //add
                     alpha = 1.0;
                     dgemm_(&notra,&trans,&dimLup,&dimRup,&dimRdown,&alpha,temp,&dimLup,Qblock,&dimRup,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }

   //Then do 3J1
   for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
      for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         if ((abs(TwoSLdown-TwoSRdown)<=TwoJ) && (TwoSLdown>=0) && (TwoSRdown>=0)){
      
            int fase = phase(TwoSL+TwoSRdown+TwoJ+1 + ((N1==1)?2:0) + ((N2==1)?2:0) );
            const double factor = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j(TwoSL,TwoSR,TwoJ,TwoSRdown,TwoSLdown,1);
      
            for (int l_index=0; l_index<theindex; l_index++){
            
               #ifdef CHEMPS2_MPI_COMPILATION
               if ( MPIchemps2::owner_q( Prob->gL(), l_index ) == MPIRANK )
               #endif
               {
                  int ILdown = Irreps::directProd(IL,denBK->gIrrep(l_index));
                  int IRdown = Irreps::directProd(IR,denBK->gIrrep(l_index));
                  int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, N1, N2, TwoJ, NR-1, TwoSRdown, IRdown);
                  if (memSkappa!=-1){
                     int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  
                     double * Lblock = Lleft[ theindex-1-l_index]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                     double * Qblock = Qright[theindex+1-l_index]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                  
                     char trans = 'T';
                     char notra = 'N';
                     double beta = 0.0; //set
                     double alpha = factor;
                     dgemm_(&trans,&notra,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,memS+denS->gKappa2index(memSkappa),&dimLdown,&beta,temp,&dimLup);
                  
                     beta = 1.0; //add
                     alpha = 1.0;
                     dgemm_(&notra,&notra,&dimLup,&dimRup,&dimRdown,&alpha,temp,&dimLup,Qblock,&dimRdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                  
                  }
               }
            }
         }
      }
   }
   
}


