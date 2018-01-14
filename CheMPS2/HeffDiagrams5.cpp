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

void CheMPS2::Heff::addDiagram5A(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char trans = 'T';
   char notrans = 'N';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5A1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5A1 ) == MPIRANK ) && (N1==0) && (N2==0)){
   #else
   if ((N1==0) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
                  int fase = phase(TwoSLdown+TwoSRdown+2);
                  const double factor = fase * sqrt((TwoJdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac =             factor * ( Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                     + ((TwoJdown==0)?1:-1) * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex) );
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                              
                                 dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

   //5A2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5A2 ) == MPIRANK ) && (N1==1) && (N2==0)){
   #else
   if ((N1==1) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
            if (((TwoSL == TwoSRdown) || (TwoSLdown == TwoSR)) && (abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               double factor1 = 0.0;
               double factor2 = 0.0;
            
               if (TwoSL == TwoSRdown){ factor2 = phase(TwoSRdown + 1 - TwoSR); }
               if (TwoSLdown == TwoSR){
                  int fase = phase(TwoSLdown + 1 - TwoSL);
                  factor1 = fase * sqrt(((TwoSRdown+1)*(TwoSL+1.0))/((TwoSR+1)*(TwoSLdown+1.0)));
               }
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                  + factor2 * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 2, 1, 1, NR+1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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
   
   //5A3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5A3 ) == MPIRANK ) && (N1==0) && (N2==1)){
   #else
   if ((N1==0) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
            if (((TwoSL == TwoSRdown) || (TwoSLdown == TwoSR)) && (abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               double factor1 = 0.0;
               double factor2 = 0.0;
            
               if (TwoSL == TwoSRdown){ factor1 = phase(TwoSRdown + 1 - TwoSR); }
               if (TwoSLdown == TwoSR){
                  int fase = phase(TwoSLdown + 1 - TwoSL);
                  factor2 = fase * sqrt(((TwoSRdown+1)*(TwoSL+1.0))/((TwoSR+1)*(TwoSLdown+1.0)));
               }
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                  + factor2 * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 2, 1, NR+1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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
   
   //5A4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5A4 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
            
            int fase = (((TwoSL+1)%2)!=0)?-1:1;
            const double factor = fase * sqrt((TwoJ+1)*(TwoSL+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
            
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
              
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                     
               if ((leftOK) && (rightOK)){
                     
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                              
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                       
                        if ((dimLdown>0) && (dimRdown>0)){
                         
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                             
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac =         factor * ( Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                               + ((TwoJ==0)?1:-1) * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex) );
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                                 
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 2, 2, 0, NR+1, TwoSRdown, IRdown);
                                 
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                                 
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

void CheMPS2::Heff::addDiagram5B(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char notrans = 'N';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5B1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5B1 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
               
            int fase = phase(TwoSL+TwoSR+2);
            const double factor = fase * sqrt((TwoJ+1)*(TwoSR+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
           
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                  
               if ((leftOK) && (rightOK)){
                  
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                              
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                              
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                             
                        if ((dimLdown>0) && (dimRdown>0)){
                             
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac =         factor * ( Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                               + ((TwoJ==0)?1:-1) * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex) );
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                                 
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 0, 0, 0, NR-1, TwoSRdown, IRdown);
                                
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                                 
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   //5B2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5B2 ) == MPIRANK ) && (N1==2) && (N2==1)){
   #else
   if ((N1==2) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){

            if (((TwoSL == TwoSRdown) || (TwoSLdown == TwoSR)) && (abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               double factor1 = 0.0;
               double factor2 = 0.0;
            
               if (TwoSLdown == TwoSR){ factor2 = phase(TwoSR + 1 - TwoSRdown); }
               if (TwoSL == TwoSRdown){
                  int fase = phase(TwoSL + 1 - TwoSLdown);
                  factor1 = fase * sqrt(((TwoSR+1)*(TwoSLdown+1.0))/((TwoSRdown+1)*(TwoSL+1.0)));
               }
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                  + factor2 * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 0, 1, NR-1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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
   
   //5B3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5B3 ) == MPIRANK ) && (N1==1) && (N2==2)){
   #else
   if ((N1==1) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
         
            if (((TwoSL == TwoSRdown) || (TwoSLdown == TwoSR)) && (abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
               double factor1 = 0.0;
               double factor2 = 0.0;
            
               if (TwoSLdown == TwoSR){ factor1 = phase(TwoSR + 1 - TwoSRdown); }
               if (TwoSL == TwoSRdown){
                  factor2 = phase(TwoSL + 1 - TwoSLdown) * sqrt(((TwoSR+1)*(TwoSLdown+1.0))/((TwoSRdown+1)*(TwoSL+1.0)));
               }
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                  + factor2 * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 0, 1, 1, NR-1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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
   
   //5B4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5B4 ) == MPIRANK ) && (N1==2) && (N2==2)){
   #else
   if ((N1==2) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
                  int fase = (((TwoSLdown+1)%2)!=0)?-1:1;
                  const double factor = fase * sqrt((TwoJdown+1)*(TwoSLdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
   
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac =             factor * ( Prob->gMxElement(l_alpha, l_beta, theindex, theindex+1) 
                                                     + ((TwoJdown==0)?1:-1) * Prob->gMxElement(l_alpha, l_beta, theindex+1, theindex) );
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                              
                                 dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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
   
}

void CheMPS2::Heff::addDiagram5C(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char notrans = 'N';
   char trans = 'T';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5C1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5C1 ) == MPIRANK ) && (N1==1) && (N2==0)){
   #else
   if ((N1==1) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSL+TwoSRdown);
               const double factor2 = fase * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j( TwoSLdown, TwoSRdown, 1, TwoSR, TwoSL, 1 );
               const double factor1 = (TwoSL==TwoSRdown)?sqrt((TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
           
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
   
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha,theindex,theindex+1,l_beta)
                                                  + factor2 * Prob->gMxElement(l_alpha,theindex,l_beta,theindex+1);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 0, 1, 1, NR-1, TwoSRdown, IRdown);
                             
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

   //5C2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5C2 ) == MPIRANK ) && (N1==2) && (N2==0)){
   #else
   if ((N1==2) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSR + TwoSLdown + 3 + TwoJdown);
                  const double factor1 = fase * sqrt((TwoSR+1)*(TwoJdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
                  const double factor2 = (TwoJdown==0)?sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                                     + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                              
                                 dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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
   
   //5C3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5C3 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
            
            int fase = phase(TwoSR + TwoSRdown + 3 + TwoJ);
            double factor1 = fase * sqrt((TwoSL+1.0)*(TwoJ+1.0)*(TwoSR+1.0)/(TwoSRdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
            double factor2 = (TwoJ==0)?sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
            
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                   
               if ((leftOK) && (rightOK)){
                  
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                           
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                        if ((dimLdown>0) && (dimRdown>0)){
                             
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                               + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                                 
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 0, 2, 0, NR-1, TwoSRdown, IRdown);
                               
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                                 
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //5C4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5C4 ) == MPIRANK ) && (N1==2) && (N2==1)){
   #else
   if ((N1==2) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor1 = (TwoSLdown==TwoSR) ? phase(TwoSL-TwoSRdown) * sqrt((TwoSL+1.0)/(TwoSLdown+1.0)) : 0.0;
               const double factor2 = phase(TwoSL+TwoSRdown+2) * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j( TwoSLdown, TwoSRdown, 1, TwoSR, TwoSL, 1 );
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                                  + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 2, 1, NR-1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

void CheMPS2::Heff::addDiagram5D(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char notrans = 'N';
   char trans = 'T';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5D1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5D1 ) == MPIRANK ) && (N1==0) && (N2==1)){
   #else
   if ((N1==0) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               int fase = phase(TwoSLdown+TwoSR);
               const double factor2 = fase * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRdown, TwoSLdown, 1 );
               const double factor1 = (TwoSLdown==TwoSR)?sqrt((TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
           
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha,theindex,theindex+1,l_beta)
                                                  + factor2 * Prob->gMxElement(l_alpha,theindex,l_beta,theindex+1);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 0, 1, NR+1, TwoSRdown, IRdown);
                             
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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

   //5D2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5D2 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
            
            int fase = phase(TwoSRdown + TwoSL + 3 + TwoJ);
            const double factor1 = fase * sqrt((TwoSRdown+1)*(TwoJ+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
            const double factor2 = (TwoJ==0)?sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
               
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
                 
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
                  
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                     
               if ((leftOK) && (rightOK)){
                     
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                              
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                              
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                              
                        if ((dimLdown>0) && (dimRdown>0)){
                                 
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                                 
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                               + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                                
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 2, 0, 0, NR+1, TwoSRdown, IRdown);
                                 
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                                 
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //5D3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5D3 ) == MPIRANK ) && (N1==0) && (N2==2)){
   #else
   if ((N1==0) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSR + TwoSRdown + 3 + TwoJdown);
                  double factor1 = fase * sqrt((TwoSLdown+1.0)*(TwoJdown+1.0)*(TwoSRdown+1.0)/(TwoSR+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
                  double factor2 = (TwoJdown==0)?sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                                     + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                            
                                 dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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
   
   //5D4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5D4 ) == MPIRANK ) && (N1==1) && (N2==2)){
   #else
   if ((N1==1) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor1 = (TwoSL==TwoSRdown) ? phase(TwoSLdown-TwoSR) * sqrt((TwoSLdown+1.0)/(TwoSL+1.0)) : 0.0;
               const double factor2 = phase(TwoSLdown+TwoSR+2) * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRdown, TwoSLdown, 1 );
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, theindex, theindex+1, l_beta) 
                                                  + factor2 * Prob->gMxElement(l_alpha, theindex, l_beta, theindex+1);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 2, 1, 1, NR+1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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

void CheMPS2::Heff::addDiagram5E(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char notrans = 'N';
   char trans = 'T';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5E1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5E1 ) == MPIRANK ) && (N1==0) && (N2==1)){
   #else
   if ((N1==0) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor2 = phase(TwoSL+TwoSRdown) * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j( TwoSLdown, TwoSRdown, 1, TwoSR, TwoSL, 1 );
               const double factor1 = (TwoSL==TwoSRdown)?sqrt((TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
           
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
                
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                   
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha,theindex+1,theindex,l_beta)
                                                  + factor2 * Prob->gMxElement(l_alpha,theindex+1,l_beta,theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 0, 1, NR-1, TwoSRdown, IRdown);
                             
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

   //5E2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5E2 ) == MPIRANK ) && (N1==0) && (N2==2)){
   #else
   if ((N1==0) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
         
                  int fase = phase(TwoSR + TwoSLdown + 3);
                  const double factor1 = fase * sqrt((TwoSR+1)*(TwoJdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
                  const double factor2 = (TwoJdown==0)?sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                                     + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR-1, TwoSRdown, IRdown);
                              
                                 dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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
   
   //5E3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5E3 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
            
            int fase = phase(TwoSR + TwoSRdown + 3);
            double factor1 = fase * sqrt((TwoSL+1.0)*(TwoJ+1.0)*(TwoSR+1.0)/(TwoSRdown+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
            double factor2 = (TwoJ==0)?sqrt(2.0*(TwoSR+1.0)/(TwoSRdown+1.0)):0.0;
               
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
                  
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
                  
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                     
               if ((leftOK) && (rightOK)){
                     
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                              
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                              
                        int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                              
                        if ((dimLdown>0) && (dimRdown>0)){
                                 
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                                 
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                               + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                                 
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 2, 0, 0, NR-1, TwoSRdown, IRdown);
                               
                           dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                                 
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
                           dgemm_(&trans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLdown,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //5E4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5E4 ) == MPIRANK ) && (N1==1) && (N2==2)){
   #else
   if ((N1==1) && (N2==2)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor1 = (TwoSLdown==TwoSR) ? phase(TwoSL-TwoSRdown) * sqrt((TwoSL+1.0)/(TwoSLdown+1.0)) : 0.0;
               const double factor2 = phase(TwoSL+TwoSRdown+2) * sqrt((TwoSL+1)*(TwoSR+1.0)) * Wigner::wigner6j( TwoSLdown, TwoSRdown, 1, TwoSR, TwoSL, 1 );
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);

                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR-1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                                  + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR-1,TwoSRdown,IRdown,NR,TwoSR,IR);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                               
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL-1, TwoSLdown, ILdown, 2, 1, 1, NR-1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&notrans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRdown,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL-1,TwoSLdown,ILdown,NL,TwoSL,IL);
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

void CheMPS2::Heff::addDiagram5F(const int ikappa, double * memS, double * memHeff, const Sobject * denS, TensorL ** Lleft, TensorL ** Lright, double * temp, double * temp2) const{

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
   int dimRup = denBK->gCurrentDim(theindex+2,NR,TwoSR,IR);
   
   char notrans = 'N';
   char trans = 'T';
   int inc = 1;
   
   int IprodMID = Irreps::directProd(denBK->gIrrep(theindex),denBK->gIrrep(theindex+1));
   
   //5F1
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5F1 ) == MPIRANK ) && (N1==1) && (N2==0)){
   #else
   if ((N1==1) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor2 = phase(TwoSLdown+TwoSR) * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRdown, TwoSLdown, 1 );
               const double factor1 = (TwoSLdown==TwoSR)?sqrt((TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
           
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha,theindex+1,theindex,l_beta)
                                                  + factor2 * Prob->gMxElement(l_alpha,theindex+1,l_beta,theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 0, 1, 1, NR+1, TwoSRdown, IRdown);
                             
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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

   //5F2
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5F2 ) == MPIRANK ) && (N1==1) && (N2==1)){
   #else
   if ((N1==1) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         if ((TwoSLdown>=0) && (abs(TwoSR-TwoSLdown)<=1)){
            int TwoSRdown = TwoSLdown;
            
            const double factor1 = phase(TwoSRdown + TwoSL + 3) * sqrt((TwoSRdown+1)*(TwoJ+1.0)) * Wigner::wigner6j( 1, 1, TwoJ, TwoSL, TwoSR, TwoSRdown );
            const double factor2 = (TwoJ==0)?sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
               
            for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
                 
               int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
               bool leftOK = false;
               for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                  if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
               }
               bool rightOK = false;
               for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                  if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
               }
                     
               if ((leftOK) && (rightOK)){
                  
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){
                              
                        int ILdown = Irreps::directProd(IL, Irrep);
                        int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                          
                        int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                        int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                        
                        if ((dimLdown>0) && (dimRdown>0)){
                           
                           int size = dimRup * dimRdown;
                           for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                            
                           for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                              if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                 double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                               + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                 double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                 daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                              }
                           }
                             
                           double alpha = 1.0;
                           double beta = 0.0; //set
                           int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 0, 2, 0, NR+1, TwoSRdown, IRdown);
                              
                           dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                           beta = 1.0; //add
                           double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
                           dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimLdown,&alpha,LblockL,&dimLup,temp2,&dimLdown,&beta,memHeff+denS->gKappa2index(ikappa),&dimLup);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   //5F3
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5F3 ) == MPIRANK ) && (N1==2) && (N2==0)){
   #else
   if ((N1==2) && (N2==0)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            for (int TwoJdown=0; TwoJdown<=2; TwoJdown+=2){
               if ((abs(TwoSLdown-TwoSRdown)<=TwoJdown) && (TwoSLdown>=0) && (TwoSRdown>=0)){
               
                  double factor1 = phase(TwoSR + TwoSRdown + 3) * sqrt((TwoSLdown+1.0)*(TwoJdown+1.0)*(TwoSRdown+1.0)/(TwoSR+1.0)) * Wigner::wigner6j( 1, 1, TwoJdown, TwoSLdown, TwoSRdown, TwoSR );
                  double factor2 = (TwoJdown==0)?sqrt(2.0*(TwoSRdown+1.0)/(TwoSR+1.0)):0.0;
            
                  for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                     int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                     bool leftOK = false;
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                     }
                     bool rightOK = false;
                     for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                        if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                     }
                  
                     if ((leftOK) && (rightOK)){
                  
                        for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                           if (denBK->gIrrep(l_alpha) == Irrep){
                           
                              int ILdown = Irreps::directProd(IL, Irrep);
                              int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                              int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                              int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                              if ((dimLdown>0) && (dimRdown>0)){
                              
                                 int size = dimRup * dimRdown;
                                 for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                                 for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                    if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                       double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                                     + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                       double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                       daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                    }
                                 }
                              
                                 double alpha = 1.0;
                                 double beta = 0.0; //set
                                 int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 1, TwoJdown, NR+1, TwoSRdown, IRdown);
                            
                                 dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                                 beta = 1.0; //add
                                 double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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
   
   //5F4
   #ifdef CHEMPS2_MPI_COMPILATION
   if (( MPIchemps2::owner_specific_diagram( Prob->gL(), MPI_CHEMPS2_5F4 ) == MPIRANK ) && (N1==2) && (N2==1)){
   #else
   if ((N1==2) && (N2==1)){
   #endif
   
      for (int TwoSLdown=TwoSL-1; TwoSLdown<=TwoSL+1; TwoSLdown+=2){
         for (int TwoSRdown=TwoSR-1; TwoSRdown<=TwoSR+1; TwoSRdown+=2){
            if ((abs(TwoSLdown-TwoSRdown)<=1) && (TwoSLdown>=0) && (TwoSRdown>=0)){
            
               const double factor1 = (TwoSL==TwoSRdown) ? phase(TwoSLdown-TwoSR) * sqrt((TwoSLdown+1.0)/(TwoSL+1.0)) : 0.0;
               const double factor2 = phase(TwoSLdown+TwoSR+2) * sqrt((TwoSLdown+1)*(TwoSRdown+1.0)) * Wigner::wigner6j( TwoSL, TwoSR, 1, TwoSRdown, TwoSLdown, 1 );
            
               for (int Irrep=0; Irrep<(denBK->getNumberOfIrreps()); Irrep++){
               
                  int IrrepTimesMid = Irreps::directProd(Irrep,IprodMID);
               
                  bool leftOK = false;
                  for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                     if (denBK->gIrrep(l_alpha) == Irrep){ leftOK = true; }
                  }
                  bool rightOK = false;
                  for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                     if (denBK->gIrrep(l_beta) == IrrepTimesMid){ rightOK = true; }
                  }
                  
                  if ((leftOK) && (rightOK)){
                  
                     for (int l_alpha=0; l_alpha<theindex; l_alpha++){
                        if (denBK->gIrrep(l_alpha) == Irrep){
                           
                           int ILdown = Irreps::directProd(IL, Irrep);
                           int IRdown = Irreps::directProd(IR, IrrepTimesMid);
                           
                           int dimLdown = denBK->gCurrentDim(theindex,   NL+1, TwoSLdown, ILdown);
                           int dimRdown = denBK->gCurrentDim(theindex+2, NR+1, TwoSRdown, IRdown);
                           
                           if ((dimLdown>0) && (dimRdown>0)){
                              
                              int size = dimRup * dimRdown;
                              for (int cnt=0; cnt<size; cnt++){ temp[cnt] = 0.0; }
                              
                              for (int l_beta=theindex+2; l_beta<Prob->gL(); l_beta++){
                                 if (denBK->gIrrep(l_beta) == IrrepTimesMid){
                                    double prefac = factor1 * Prob->gMxElement(l_alpha, theindex+1, theindex, l_beta) 
                                                  + factor2 * Prob->gMxElement(l_alpha, theindex+1, l_beta, theindex);
                                    double * LblockR = Lright[l_beta-theindex-2]->gStorage(NR,TwoSR,IR,NR+1,TwoSRdown,IRdown);
                                    daxpy_(&size,&prefac,LblockR,&inc,temp,&inc);
                                 }
                              }
                              
                              double alpha = 1.0;
                              double beta = 0.0; //set
                              int memSkappa = denS->gKappa(NL+1, TwoSLdown, ILdown, 1, 2, 1, NR+1, TwoSRdown, IRdown);
                              
                              dgemm_(&notrans,&trans,&dimLdown,&dimRup,&dimRdown,&alpha,memS+denS->gKappa2index(memSkappa),&dimLdown,temp,&dimRup,&beta,temp2,&dimLdown);
                              
                              beta = 1.0; //add
                              double * LblockL = Lleft[theindex-1-l_alpha]->gStorage(NL,TwoSL,IL,NL+1,TwoSLdown,ILdown);
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



