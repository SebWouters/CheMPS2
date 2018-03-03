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

#include <stdlib.h> /*rand*/
#include <algorithm>
#include <math.h>

#include "TensorT.h"
#include "Lapack.h"

using std::min;

CheMPS2::TensorT::TensorT( const int site_index, const SyBookkeeper * denBK ) : Tensor(){

   this->index = site_index; //left boundary = index ; right boundary = index+1
   this->denBK = denBK;

   AllocateAllArrays();

}

void CheMPS2::TensorT::AllocateAllArrays(){

   nKappa = 0;
   for ( int NL = denBK->gNmin( index ); NL <= denBK->gNmax( index ); NL++ ){
      for ( int TwoSL = denBK->gTwoSmin( index, NL ); TwoSL <= denBK->gTwoSmax( index, NL ); TwoSL += 2 ){
         for ( int IL = 0; IL < denBK->getNumberOfIrreps(); IL++ ){
            const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
            if ( dimL > 0 ){
               for ( int NR = NL; NR <= NL+2; NR++ ){
                  const int TwoJ = (( NR == NL + 1 ) ? 1 : 0 );
                  for ( int TwoSR = TwoSL - TwoJ; TwoSR <= TwoSL + TwoJ; TwoSR += 2 ){
                     if ( TwoSR >= 0 ){
                        int IR = (( NR == NL + 1 ) ? Irreps::directProd( IL, denBK->gIrrep( index ) ) : IL );
                        const int dimR = denBK->gCurrentDim( index + 1, NR, TwoSR, IR );
                        if ( dimR > 0 ){
                           nKappa++;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   sectorNL    = new int[ nKappa ];
   sectorNR    = new int[ nKappa ];
   sectorIL    = new int[ nKappa ];
   sectorIR    = new int[ nKappa ];
   sectorTwoSL = new int[ nKappa ];
   sectorTwoSR = new int[ nKappa ];
   kappa2index = new int[ nKappa + 1 ];
   kappa2index[ 0 ] = 0;

   nKappa = 0;
   for ( int NL = denBK->gNmin( index ); NL <= denBK->gNmax( index ); NL++ ){
      for ( int TwoSL = denBK->gTwoSmin( index, NL ); TwoSL <= denBK->gTwoSmax( index, NL ); TwoSL += 2 ){
         for ( int IL = 0; IL < denBK->getNumberOfIrreps(); IL++ ){
            const int dimL = denBK->gCurrentDim( index, NL, TwoSL, IL );
            if ( dimL > 0 ){
               for ( int NR = NL; NR <= NL+2; NR++ ){
                  const int TwoJ = (( NR == NL + 1 ) ? 1 : 0 );
                  for ( int TwoSR = TwoSL - TwoJ; TwoSR <= TwoSL + TwoJ; TwoSR += 2 ){
                     if ( TwoSR >= 0 ){
                        int IR = (( NR == NL + 1 ) ? Irreps::directProd( IL, denBK->gIrrep( index ) ) : IL );
                        const int dimR = denBK->gCurrentDim( index + 1, NR, TwoSR, IR );
                        if ( dimR > 0 ){
                           sectorNL[ nKappa ] = NL;
                           sectorNR[ nKappa ] = NR;
                           sectorIL[ nKappa ] = IL;
                           sectorIR[ nKappa ] = IR;
                           sectorTwoSL[ nKappa ] = TwoSL;
                           sectorTwoSR[ nKappa ] = TwoSR;
                           kappa2index[ nKappa + 1 ] = kappa2index[ nKappa ] + dimL * dimR;
                           nKappa++;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   storage = new double[ kappa2index[ nKappa ] ];

}

CheMPS2::TensorT::~TensorT(){

   DeleteAllArrays();

}

void CheMPS2::TensorT::DeleteAllArrays(){

   delete [] sectorNL;
   delete [] sectorNR;
   delete [] sectorIL;
   delete [] sectorIR;
   delete [] sectorTwoSL;
   delete [] sectorTwoSR;
   delete [] kappa2index;
   delete [] storage;

}

void CheMPS2::TensorT::Reset(){

   DeleteAllArrays();
   AllocateAllArrays();

}

int CheMPS2::TensorT::gNKappa() const { return nKappa; }

double * CheMPS2::TensorT::gStorage() { return storage; }

int CheMPS2::TensorT::gKappa( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 ) const{

   for ( int cnt = 0; cnt < nKappa; cnt++ ){
      if (( sectorNL[ cnt ] == N1 ) &&
          ( sectorNR[ cnt ] == N2 ) &&
          ( sectorIL[ cnt ] == I1 ) &&
          ( sectorIR[ cnt ] == I2 ) &&
          ( sectorTwoSL[ cnt ] == TwoS1 ) &&
          ( sectorTwoSR[ cnt ] == TwoS2 )){ return cnt; }
   }

   return -1;

}

int CheMPS2::TensorT::gKappa2index( const int kappa ) const{ return kappa2index[ kappa ]; }

double * CheMPS2::TensorT::gStorage( const int N1, const int TwoS1, const int I1, const int N2, const int TwoS2, const int I2 ){

   int kappa = gKappa( N1, TwoS1, I1, N2, TwoS2, I2 );
   if ( kappa == -1 ){ return NULL; }
   return storage + kappa2index[ kappa ];

}

int CheMPS2::TensorT::gIndex() const { return index; }

const CheMPS2::SyBookkeeper * CheMPS2::TensorT::gBK() const{ return denBK; }

void CheMPS2::TensorT::sBK( const SyBookkeeper * newBK ){ denBK = newBK; }

void CheMPS2::TensorT::random(){

   for ( int cnt = 0; cnt < kappa2index[ nKappa ]; cnt++ ){
      storage[ cnt ] = ( 2 * ( (double) rand() ) / RAND_MAX ) - 1.0; // Value in [-1,1[
   }

}

void CheMPS2::TensorT::number_operator( const double alpha, const double beta ){

   #pragma omp parallel for schedule(dynamic)
   for ( int ikappa = 0; ikappa < nKappa; ikappa++ ){
      int size = kappa2index[ ikappa + 1 ] - kappa2index[ ikappa ];
      double * array = storage + kappa2index[ ikappa ];
      double factor = beta + alpha * ( sectorNR[ ikappa ] - sectorNL[ ikappa ] );
      int inc1 = 1;
      dscal_( &size, &factor, array, &inc1 );
   }

}

void CheMPS2::TensorT::QR(Tensor * Rstorage){

   //Left normalization occurs in T-convention: no pre or after multiplication
   //Work per right symmetry sector
   
   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int NR=denBK->gNmin(index+1); NR<=denBK->gNmax(index+1); NR++){
      for (int TwoSR=denBK->gTwoSmin(index+1,NR); TwoSR<=denBK->gTwoSmax(index+1,NR); TwoSR+=2){
         for (int IR=0; IR<denBK->getNumberOfIrreps(); IR++){
            int dimR = denBK->gCurrentDim(index+1,NR,TwoSR,IR);
            if (dimR>0){
            
               //Find out the total left dimension
               int dimLtotal = 0;
               for (int ikappa=0; ikappa<nKappa; ikappa++){
                  if ((NR==sectorNR[ikappa])&&(TwoSR==sectorTwoSR[ikappa])&&(IR==sectorIR[ikappa])){
                     dimLtotal += denBK->gCurrentDim(index,sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
                  }
               }
               
               if (dimLtotal>0){ //Due to the initial truncation, it is possible that dimLtotal is temporarily smaller than dimR ...
               
                  double * mem = new double[dimLtotal*dimR];
                  //Copy the relevant parts from storage to mem
                  int dimLtotal2 = 0;
                  for (int ikappa=0; ikappa<nKappa; ikappa++){
                     if ((NR==sectorNR[ikappa])&&(TwoSR==sectorTwoSR[ikappa])&&(IR==sectorIR[ikappa])){
                        int dimL = denBK->gCurrentDim(index,sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
                        if (dimL>0){
                           for (int l=0; l<dimL; l++){
                              for (int r=0; r<dimR; r++){
                                 mem[dimLtotal2 + l + dimLtotal * r] = storage[kappa2index[ikappa] + l + dimL * r];
                              }
                           }
                           dimLtotal2 += dimL;
                        }
                     }
                  }
               
                  //QR mem --> m = dimLtotal ; n = dimR
                  int info;
                  int minofdims = min(dimR,dimLtotal);
                  double * tau = new double[minofdims];
                  double * work = new double[dimR];
                  dgeqrf_(&dimLtotal,&dimR,mem,&dimLtotal,tau,work,&dimR,&info);
               
                  //Copy R to Rstorage
                  double * wheretoput = Rstorage->gStorage(NR,TwoSR,IR,NR,TwoSR,IR); //dimR x dimR
               
                  for (int irow = 0; irow<minofdims; irow++){
                     for (int icol = 0; icol<irow; icol++){
                        wheretoput[irow + dimR * icol] = 0.0;
                     }
                     for (int icol = irow; icol<dimR; icol++){
                        wheretoput[irow + dimR * icol] = mem[irow + dimLtotal * icol];
                     }
                  }
                  for (int irow = minofdims; irow<dimR; irow++){
                     for (int icol = 0; icol<dimR; icol++){
                        wheretoput[irow + dimR * icol] = 0.0;
                     }
                  }
               
                  //Construct Q
                  dorgqr_(&dimLtotal,&minofdims,&minofdims,mem,&dimLtotal,tau,work,&dimR,&info);
                  if (dimLtotal < dimR){ //if number of cols larger than number of rows, rest of cols zero.
                     for (int irow=0; irow<dimLtotal; irow++){
                        for (int icol=dimLtotal; icol<dimR; icol++){
                           mem[irow + dimLtotal * icol] = 0.0;
                        }
                     }
                  }
              
                  //Copy from mem to storage
                  dimLtotal2 = 0;
                  for (int ikappa=0; ikappa<nKappa; ikappa++){
                     if ((NR==sectorNR[ikappa])&&(TwoSR==sectorTwoSR[ikappa])&&(IR==sectorIR[ikappa])){
                        int dimL = denBK->gCurrentDim(index,sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
                        if (dimL>0){
                           for (int l=0; l<dimL; l++){
                              for (int r=0; r<dimR; r++){
                                 storage[kappa2index[ikappa] + l + dimL * r] = mem[dimLtotal2 + l + dimLtotal * r];
                              }
                           }
                           dimLtotal2 += dimL;
                        }
                     }
                  }
                
                  //Clear the memory
                  delete [] work;
                  delete [] tau;
                  delete [] mem;
               
               }
            }
         }
      }
   }

}

void CheMPS2::TensorT::LQ(Tensor * Lstorage){

   //Right normalization occurs in U-convention: pre-multiplication with sqrt{2jR+1/2jL+1} and after multiplication with sqrt{2jL+1/2jR+1}
   //Work per left symmetry sector
   
   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int NL=denBK->gNmin(index); NL<=denBK->gNmax(index); NL++){
      for (int TwoSL=denBK->gTwoSmin(index,NL); TwoSL<=denBK->gTwoSmax(index,NL); TwoSL+=2){
         for (int IL=0; IL<denBK->getNumberOfIrreps(); IL++){
            int dimL = denBK->gCurrentDim(index,NL,TwoSL,IL);
            if (dimL>0){
            
               //Find out the total right dimension
               int dimRtotal = 0;
               for (int ikappa=0; ikappa<nKappa; ikappa++){
                  if ((NL==sectorNL[ikappa])&&(TwoSL==sectorTwoSL[ikappa])&&(IL==sectorIL[ikappa])){
                     dimRtotal += denBK->gCurrentDim(index+1,sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
                  }
               }
               
               if (dimRtotal>0){ //Due to the initial truncation, it is possible that dimRtotal is temporarily smaller than dimL ...
               
                  double * mem = new double[dimRtotal*dimL];
                  //Copy the relevant parts from storage to mem & multiply with factor !!
                  int dimRtotal2 = 0;
                  for (int ikappa=0; ikappa<nKappa; ikappa++){
                     if ((NL==sectorNL[ikappa])&&(TwoSL==sectorTwoSL[ikappa])&&(IL==sectorIL[ikappa])){
                        int dimR = denBK->gCurrentDim(index+1,sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
                        if (dimR>0){
                           double factor = sqrt((sectorTwoSR[ikappa]+1.0)/(TwoSL+1.0));
                           for (int l=0; l<dimL; l++){
                              for (int r=0; r<dimR; r++){
                                 mem[l + dimL * (dimRtotal2 + r)] = factor * storage[kappa2index[ikappa] + l + dimL * r];
                              }
                           }
                           dimRtotal2 += dimR;
                        }
                     }
                  }
               
                  //LQ mem --> m = dimL ; n = dimRtotal
                  int info;
                  int minofdims = min(dimL,dimRtotal);
                  double * tau = new double[minofdims];
                  double * work = new double[dimL];
                  dgelqf_(&dimL,&dimRtotal,mem,&dimL,tau,work,&dimL,&info);

                  //Copy L to Lstorage
                  double * wheretoput = Lstorage->gStorage(NL,TwoSL,IL,NL,TwoSL,IL); //dimL x dimL
               
                  for (int irow = 0; irow<dimL; irow++){
                     for (int icol = 0; icol<min(irow+1,dimRtotal); icol++){ //icol can be max. irow and max. dimRtotal-1
                        wheretoput[irow + dimL * icol] = mem[irow + dimL * icol];
                     }
                     for (int icol = min(irow+1,dimRtotal); icol<dimL; icol++){
                        wheretoput[irow + dimL * icol] = 0.0;
                     }
                  }
                  
                  //Construct Q
                  dorglq_(&minofdims,&dimRtotal,&minofdims,mem,&dimL,tau,work,&dimL,&info);
                  if (dimRtotal < dimL){ //if number of rows larger than number of cols, rest of rows zero.
                     for (int irow=dimRtotal; irow<dimL; irow++){
                        for (int icol=0; icol<dimRtotal; icol++){
                           mem[irow + dimL * icol] = 0.0;
                        }
                     }
                  }
              
                  //Copy from mem to storage & multiply with factor !!
                  dimRtotal2 = 0;
                  for (int ikappa=0; ikappa<nKappa; ikappa++){
                     if ((NL==sectorNL[ikappa])&&(TwoSL==sectorTwoSL[ikappa])&&(IL==sectorIL[ikappa])){
                        int dimR = denBK->gCurrentDim(index+1,sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
                        if (dimR>0){
                           double factor = sqrt((TwoSL+1.0)/(sectorTwoSR[ikappa]+1.0));
                           for (int l=0; l<dimL; l++){
                              for (int r=0; r<dimR; r++){
                                 storage[kappa2index[ikappa] + l + dimL * r] = factor * mem[l + dimL * (r + dimRtotal2)];
                              }
                           }
                           dimRtotal2 += dimR;
                        }
                     }
                  }
                
                  //Clear the memory
                  delete [] work;
                  delete [] tau;
                  delete [] mem;
               
               }
            }
         }
      }
   }

}

void CheMPS2::TensorT::LeftMultiply(Tensor * Mx){

   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimL = denBK->gCurrentDim(index,sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
      int dimR = denBK->gCurrentDim(index+1,sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
      double * MxBlock = Mx->gStorage(sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa],sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
      char notrans = 'N';
      double one = 1.0;
      double zero = 0.0;
      int dim = dimL*dimR;
      double * mem = new double[dim];
      dgemm_(&notrans,&notrans,&dimL,&dimR,&dimL,&one,MxBlock,&dimL,storage+kappa2index[ikappa],&dimL,&zero,mem,&dimL);
      int inc = 1;
      dcopy_(&dim,mem,&inc,storage+kappa2index[ikappa],&inc);
      delete [] mem;
   }
   
}

void CheMPS2::TensorT::RightMultiply(Tensor * Mx){

   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int ikappa=0; ikappa<nKappa; ikappa++){
      int dimL = denBK->gCurrentDim(index,sectorNL[ikappa],sectorTwoSL[ikappa],sectorIL[ikappa]);
      int dimR = denBK->gCurrentDim(index+1,sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
      double * MxBlock = Mx->gStorage(sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa],sectorNR[ikappa],sectorTwoSR[ikappa],sectorIR[ikappa]);
      char notrans = 'N';
      double one = 1.0;
      double zero = 0.0;
      int dim = dimL*dimR;
      double * mem = new double[dim];
      dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&one,storage+kappa2index[ikappa],&dimL,MxBlock,&dimR,&zero,mem,&dimL);
      int inc = 1;
      dcopy_(&dim,mem,&inc,storage+kappa2index[ikappa],&inc);
      delete [] mem;
   }
   
}

bool CheMPS2::TensorT::CheckLeftNormal() const{

   bool isLeftNormal = true;

   for (int NR=denBK->gNmin(index+1); NR<=denBK->gNmax(index+1); NR++){
      for (int TwoSR=denBK->gTwoSmin(index+1,NR); TwoSR<=denBK->gTwoSmax(index+1,NR); TwoSR+=2){
         for (int IR=0; IR<denBK->getNumberOfIrreps(); IR++){
            int dimR = denBK->gCurrentDim(index+1,NR,TwoSR,IR);
            if (dimR>0){
               double * result = new double[dimR*dimR];
               bool firsttime = true;
               for (int NL=NR-2; NL<=NR; NL++){
                  for (int TwoSL=TwoSR-((NR==NL+1)?1:0); TwoSL<TwoSR+2; TwoSL+=2){
                     int IL = (NR==NL+1)?(Irreps::directProd(denBK->gIrrep(index),IR)):IR;
                     int dimL = denBK->gCurrentDim(index,NL,TwoSL,IL);
                     if (dimL>0){
                        double * Block = storage + kappa2index[gKappa(NL,TwoSL,IL,NR,TwoSR,IR)];
                        char trans = 'T';
                        char notrans = 'N';
                        double one = 1.0;
                        double beta = (firsttime)?0.0:1.0;
                        dgemm_(&trans,&notrans,&dimR,&dimR,&dimL,&one,Block,&dimL,Block,&dimL,&beta,result,&dimR);
                        firsttime = false;
                     }
                  }
               }
               for (int cnt=0; cnt<dimR; cnt++) result[(dimR+1)*cnt] -= 1.0;
               char norm = 'F'; //Frobenius norm
               char uplo = 'U'; //Doesn't matter as result is fully filled
               double TwoNorm = dlansy_(&norm,&uplo,&dimR,result,&dimR,result); //last result is work space, not referenced as norm='F';
               if (TwoNorm>CheMPS2::TENSORT_orthoComparison) isLeftNormal = false;
               delete [] result;
            }
         }
      }
   }
   
   return isLeftNormal;
            
}

bool CheMPS2::TensorT::CheckRightNormal() const{

   bool isRightNormal = true;

   for (int NL=denBK->gNmin(index); NL<=denBK->gNmax(index); NL++){
      for (int TwoSL=denBK->gTwoSmin(index,NL); TwoSL<=denBK->gTwoSmax(index,NL); TwoSL+=2){
         for (int IL=0; IL<denBK->getNumberOfIrreps(); IL++){
            int dimL = denBK->gCurrentDim(index,NL,TwoSL,IL);
            if (dimL>0){
               double * result = new double[dimL*dimL];
               bool firsttime = true;
               for (int NR=NL; NR<=NL+2; NR++){
                  for (int TwoSR=TwoSL-((NR==NL+1)?1:0); TwoSR<TwoSL+2; TwoSR+=2){
                     int IR = (NR==NL+1)?(Irreps::directProd(denBK->gIrrep(index),IL)):IL;
                     int dimR = denBK->gCurrentDim(index+1,NR,TwoSR,IR);
                     if (dimR>0){
                        double * Block = storage + kappa2index[gKappa(NL,TwoSL,IL,NR,TwoSR,IR)];
                        char trans = 'T';
                        char notrans = 'N';
                        double alpha = (TwoSR+1.0)/(TwoSL+1.0);
                        double beta = (firsttime)?0.0:1.0;
                        
                        dgemm_(&notrans,&trans,&dimL,&dimL,&dimR,&alpha,Block,&dimL,Block,&dimL,&beta,result,&dimL);
                        firsttime = false;
                     }
                  }
               }
               for (int cnt=0; cnt<dimL; cnt++) result[(dimL+1)*cnt] -= 1.0;
               char norm = 'F'; //Frobenius norm
               char uplo = 'U'; //Doesn't matter as result is fully filled
               double TwoNorm = dlansy_(&norm,&uplo,&dimL,result,&dimL,result); //last result is work space, not referenced as norm='F';
               if (TwoNorm>CheMPS2::TENSORT_orthoComparison) isRightNormal = false;
               delete [] result;
            }
         }
      }
   }
   
   return isRightNormal;
            
}


