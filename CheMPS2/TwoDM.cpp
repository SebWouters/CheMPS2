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
#include <algorithm>

#include "TwoDM.h"
#include "Lapack.h"
#include "Gsl.h"
#include "Options.h"

using std::max;

CheMPS2::TwoDM::TwoDM(const SyBookkeeper * denBKIn, const Problem * ProbIn){

   denBK = denBKIn;
   Prob = ProbIn;
   L = denBK->gL();
   
   orb2IndexSy = new int[L];
   irrep2num_orb = new int[denBK->getNumberOfIrreps()];
   for (int cnt=0; cnt<denBK->getNumberOfIrreps(); cnt++){ irrep2num_orb[cnt] = 0; }
   for (int orb=0; orb<L; orb++){
      const int irrep = Prob->gIrrep( orb ); //Prob assumes you use DMRG orbs...
      orb2IndexSy[ orb ] = irrep2num_orb[ irrep ];
      irrep2num_orb[ irrep ] += 1;
   }
   TwoDMA = new TwoDMstorage( Prob->gSy() , irrep2num_orb );
   TwoDMB = new TwoDMstorage( Prob->gSy() , irrep2num_orb );

}

CheMPS2::TwoDM::~TwoDM(){

   delete TwoDMA;
   delete TwoDMB;
   delete [] orb2IndexSy;
   delete [] irrep2num_orb;

}

void CheMPS2::TwoDM::setTwoDMA_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value){

   //Prob assumes you use DMRG orbs...
   //Irrep sanity checks are performed in TwoDM::FillSite
   TwoDMA->set( Prob->gIrrep(cnt1), Prob->gIrrep(cnt2), Prob->gIrrep(cnt3), Prob->gIrrep(cnt4),
                orb2IndexSy[cnt1],  orb2IndexSy[cnt2],  orb2IndexSy[cnt3],  orb2IndexSy[cnt4],  value );
   TwoDMA->set( Prob->gIrrep(cnt2), Prob->gIrrep(cnt1), Prob->gIrrep(cnt4), Prob->gIrrep(cnt3),
                orb2IndexSy[cnt2],  orb2IndexSy[cnt1],  orb2IndexSy[cnt4],  orb2IndexSy[cnt3],  value );
   TwoDMA->set( Prob->gIrrep(cnt3), Prob->gIrrep(cnt4), Prob->gIrrep(cnt1), Prob->gIrrep(cnt2),
                orb2IndexSy[cnt3],  orb2IndexSy[cnt4],  orb2IndexSy[cnt1],  orb2IndexSy[cnt2],  value );
   TwoDMA->set( Prob->gIrrep(cnt4), Prob->gIrrep(cnt3), Prob->gIrrep(cnt2), Prob->gIrrep(cnt1),
                orb2IndexSy[cnt4],  orb2IndexSy[cnt3],  orb2IndexSy[cnt2],  orb2IndexSy[cnt1],  value );


}

void CheMPS2::TwoDM::setTwoDMB_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4, const double value){

   //Prob assumes you use DMRG orbs...
   //Irrep sanity checks are performed in TwoDM::FillSite
   TwoDMB->set( Prob->gIrrep(cnt1), Prob->gIrrep(cnt2), Prob->gIrrep(cnt3), Prob->gIrrep(cnt4),
                orb2IndexSy[cnt1],  orb2IndexSy[cnt2],  orb2IndexSy[cnt3],  orb2IndexSy[cnt4],  value );
   TwoDMB->set( Prob->gIrrep(cnt2), Prob->gIrrep(cnt1), Prob->gIrrep(cnt4), Prob->gIrrep(cnt3),
                orb2IndexSy[cnt2],  orb2IndexSy[cnt1],  orb2IndexSy[cnt4],  orb2IndexSy[cnt3],  value );
   TwoDMB->set( Prob->gIrrep(cnt3), Prob->gIrrep(cnt4), Prob->gIrrep(cnt1), Prob->gIrrep(cnt2),
                orb2IndexSy[cnt3],  orb2IndexSy[cnt4],  orb2IndexSy[cnt1],  orb2IndexSy[cnt2],  value );
   TwoDMB->set( Prob->gIrrep(cnt4), Prob->gIrrep(cnt3), Prob->gIrrep(cnt2), Prob->gIrrep(cnt1),
                orb2IndexSy[cnt4],  orb2IndexSy[cnt3],  orb2IndexSy[cnt2],  orb2IndexSy[cnt1],  value );

}

double CheMPS2::TwoDM::getTwoDMA_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   const int irrep3 = Prob->gIrrep(cnt3);
   const int irrep4 = Prob->gIrrep(cnt4);
   if ( Irreps::directProd(irrep1, irrep2) == Irreps::directProd(irrep3, irrep4) ){
      return TwoDMA->get( irrep1, irrep2, irrep3, irrep4, orb2IndexSy[cnt1], orb2IndexSy[cnt2], orb2IndexSy[cnt3], orb2IndexSy[cnt4] );
   }
   
   return 0.0;

}

double CheMPS2::TwoDM::getTwoDMB_DMRG(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs...
   const int irrep1 = Prob->gIrrep(cnt1);
   const int irrep2 = Prob->gIrrep(cnt2);
   const int irrep3 = Prob->gIrrep(cnt3);
   const int irrep4 = Prob->gIrrep(cnt4);
   if ( Irreps::directProd(irrep1, irrep2) == Irreps::directProd(irrep3, irrep4) ){
      return TwoDMB->get( irrep1, irrep2, irrep3, irrep4, orb2IndexSy[cnt1], orb2IndexSy[cnt2], orb2IndexSy[cnt3], orb2IndexSy[cnt4] );
   }
   
   return 0.0;

}

double CheMPS2::TwoDM::getTwoDMA_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorderD2h() ){
      return getTwoDMA_DMRG( Prob->gf1(cnt1), Prob->gf1(cnt2), Prob->gf1(cnt3), Prob->gf1(cnt4) );
   }
   return getTwoDMA_DMRG( cnt1, cnt2, cnt3, cnt4 );

}

double CheMPS2::TwoDM::getTwoDMB_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorderD2h() ){
      return getTwoDMB_DMRG( Prob->gf1(cnt1), Prob->gf1(cnt2), Prob->gf1(cnt3), Prob->gf1(cnt4) );
   }
   return getTwoDMB_DMRG( cnt1, cnt2, cnt3, cnt4 );

}

double CheMPS2::TwoDM::doubletrace2DMA(){

   double val = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         val += getTwoDMA_DMRG(cnt1,cnt2,cnt1,cnt2);
      }
   }
   return val;

}

double CheMPS2::TwoDM::calcEnergy(){

   double val = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         for (int cnt3=0; cnt3<L; cnt3++){
            for (int cnt4=0; cnt4<L; cnt4++){
               val += getTwoDMA_DMRG(cnt1,cnt2,cnt3,cnt4) * Prob->gMxElement(cnt1,cnt2,cnt3,cnt4);
            }
         }
      }
   }
   val *= 0.5;
   return val + Prob->gEconst();

}

void CheMPS2::TwoDM::save(){

   TwoDMA->save(CheMPS2::TWODM_2DM_A_storagename);
   TwoDMB->save(CheMPS2::TWODM_2DM_B_storagename);

}

void CheMPS2::TwoDM::read(){

   TwoDMA->read(CheMPS2::TWODM_2DM_A_storagename);
   TwoDMB->read(CheMPS2::TWODM_2DM_B_storagename);

}

int CheMPS2::TwoDM::trianglefunction(const int k, const int glob){

   int cnt2tilde = 1;
   while(cnt2tilde*(cnt2tilde+1)/2 <= glob){ cnt2tilde++; }
   return k - cnt2tilde;
   
}


void CheMPS2::TwoDM::FillSite(TensorT * denT, TensorL *** Ltens, TensorF0 **** F0tens, TensorF1 **** F1tens, TensorS0 **** S0tens, TensorS1 **** S1tens){

   const int theindex = denT->gIndex();
   const int DIM = max(denBK->gMaxDimAtBound(theindex), denBK->gMaxDimAtBound(theindex+1));
   const double prefactorSpin = 1.0/(Prob->gTwoS() + 1.0);
   
   //Diagram 1
   const double d1 = doD1(denT) * prefactorSpin;
   setTwoDMA_DMRG(theindex,theindex,theindex,theindex, 2*d1);
   setTwoDMB_DMRG(theindex,theindex,theindex,theindex,-2*d1);
   
   #pragma omp parallel
   {
   
      double * workmem = new double[DIM*DIM];
      double * workmem2 = new double[DIM*DIM];

      //Diagram 2
      #pragma omp for schedule(static) nowait
      for (int j_index=theindex+1; j_index<L; j_index++){
         if (denBK->gIrrep(j_index) == denBK->gIrrep(theindex)){
            const double d2 = doD2(denT, Ltens[theindex][j_index-theindex-1], workmem) * prefactorSpin;
            setTwoDMA_DMRG(theindex,j_index,theindex,theindex, 2*d2);
            setTwoDMB_DMRG(theindex,j_index,theindex,theindex,-2*d2);
         }
      }

      /*for (int j_index=theindex+1; j_index<L; j_index++){
         for (int k_index=j_index; k_index<L; k_index++){*/
      const int dimTriangle = L - theindex - 1;
      const int upperboundTriangle = dimTriangle*(dimTriangle+1)/2;
      #pragma omp for schedule(static) nowait
      for (int global=0; global<upperboundTriangle; global++){
         const int row = trianglefunction(dimTriangle, global);
         const int col = global - (dimTriangle-row)*(dimTriangle-1-row)/2;
         const int j_index = theindex + 1 + row;
         const int k_index = j_index + col;
         if (denBK->gIrrep(j_index) == denBK->gIrrep(k_index)){
         
            //Diagram 3
            const double d3 = doD3(denT, S0tens[theindex][k_index-j_index][j_index-theindex-1], workmem) * prefactorSpin;
            setTwoDMA_DMRG(theindex,theindex,j_index,k_index, 2*d3);
            setTwoDMB_DMRG(theindex,theindex,j_index,k_index,-2*d3);
            
            //Diagrams 4,5 and 6
            const double d4 = doD4(denT, F0tens[theindex][k_index-j_index][j_index-theindex-1], workmem) * prefactorSpin;
            const double d5 = doD5(denT, F0tens[theindex][k_index-j_index][j_index-theindex-1], workmem) * prefactorSpin;
            const double d6 = doD6(denT, F1tens[theindex][k_index-j_index][j_index-theindex-1], workmem) * prefactorSpin;
            setTwoDMA_DMRG(theindex,j_index,k_index,theindex, -2*d4 - 2*d5 - 3*d6);
            setTwoDMB_DMRG(theindex,j_index,k_index,theindex, -2*d4 - 2*d5 +   d6);
            setTwoDMA_DMRG(theindex,j_index,theindex,k_index,  4*d4 + 4*d5);
            setTwoDMB_DMRG(theindex,j_index,theindex,k_index,  2*d6);
            
         }
      }

      //Diagram 7
      #pragma omp for schedule(static) nowait
      for (int g_index=0; g_index<theindex; g_index++){
         if (denBK->gIrrep(g_index) == denBK->gIrrep(theindex)){
            const double d7 = doD7(denT, Ltens[theindex-1][theindex-g_index-1], workmem) * prefactorSpin;
            setTwoDMA_DMRG(g_index,theindex,theindex,theindex, 2*d7);
            setTwoDMB_DMRG(g_index,theindex,theindex,theindex,-2*d7);
         }
      }

      /*for (int g_index=0; g_index<theindex; g_index++){
         for (int j_index=theindex+1; j_index<L; j_index++){*/
      const int globalsize8to12 = theindex * ( L - 1 - theindex );
      #pragma omp for schedule(static) nowait
      for (int gj_index=0; gj_index<globalsize8to12; gj_index++){
         const int g_index = gj_index % theindex;
         const int j_index = ( gj_index / theindex ) + theindex + 1;
         const int I_g = denBK->gIrrep(g_index);
         if (denBK->gIrrep(g_index) == denBK->gIrrep(j_index)){
            //Diagrams 8,9,10 and 11
            const double d8 = doD8(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, I_g) * prefactorSpin;
            double d9, d10, d11;
            doD9andD10andD11(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, &d9, &d10, &d11, I_g);
            d9 *= prefactorSpin;
            d10 *= prefactorSpin;
            d11 *= prefactorSpin;
            setTwoDMA_DMRG(g_index,theindex,j_index,theindex, -4*d8-d9);
            setTwoDMA_DMRG(g_index,theindex,theindex,j_index, 2*d8 + d11);
            setTwoDMB_DMRG(g_index,theindex,j_index,theindex, d9 - 2*d10);
            setTwoDMB_DMRG(g_index,theindex,theindex,j_index, 2*d8 + 2*d10 - d11);
            
            //Diagram 12
            const double d12 = doD12(denT, Ltens[theindex-1][theindex-g_index-1], Ltens[theindex][j_index-theindex-1], workmem, workmem2, I_g) * prefactorSpin;
            setTwoDMA_DMRG(g_index,j_index,theindex,theindex, 2*d12);
            setTwoDMB_DMRG(g_index,j_index,theindex,theindex,-2*d12);
         }
      }

      /*for (int g_index=0; g_index<theindex; g_index++){
         for (int j_index=theindex+1; j_index<L; j_index++){
            for (int k_index=j_index; k_index<L; k_index++){*/
      const int globalsize = theindex * upperboundTriangle;
      #pragma omp for schedule(static) nowait
      for (int gjk_index=0; gjk_index<globalsize; gjk_index++){
         const int g_index = gjk_index % theindex;
         const int global  = gjk_index / theindex;
         const int row = trianglefunction(dimTriangle, global);
         const int col = global - (dimTriangle-row)*(dimTriangle-1-row)/2;
         const int j_index = theindex + 1 + row;
         const int k_index = j_index + col;
         const int I_g = denBK->gIrrep(g_index);

         if (Irreps::directProd(I_g, denBK->gIrrep(theindex)) == Irreps::directProd(denBK->gIrrep(j_index), denBK->gIrrep(k_index))){
            //Diagrams 13,14,15 and 16
            const double d13 = doD13(denT, Ltens[theindex-1][theindex-g_index-1], S0tens[theindex][k_index-j_index][j_index-theindex-1],
                                     workmem, workmem2, I_g) * prefactorSpin;
            const double d14 = doD14(denT, Ltens[theindex-1][theindex-g_index-1], S0tens[theindex][k_index-j_index][j_index-theindex-1],
                                     workmem, workmem2, I_g) * prefactorSpin;
            double d15 = 0.0;
            double d16 = 0.0;
            if (k_index>j_index){
               d15 = doD15(denT, Ltens[theindex-1][theindex-g_index-1], S1tens[theindex][k_index-j_index][j_index-theindex-1], workmem, workmem2, I_g) * prefactorSpin;
               d16 = doD16(denT, Ltens[theindex-1][theindex-g_index-1], S1tens[theindex][k_index-j_index][j_index-theindex-1], workmem, workmem2, I_g) * prefactorSpin;
            }
            setTwoDMA_DMRG(g_index,theindex,j_index,k_index, 2*d13 + 2*d14 + 3*d15 + 3*d16);
            setTwoDMA_DMRG(g_index,theindex,k_index,j_index, 2*d13 + 2*d14 - 3*d15 - 3*d16);
            setTwoDMB_DMRG(g_index,theindex,j_index,k_index,-2*d13 - 2*d14 +   d15 +   d16);
            setTwoDMB_DMRG(g_index,theindex,k_index,j_index,-2*d13 - 2*d14 -   d15 -   d16);
            
            //Diagrams 17,18,19 and 20
            const double d17 = doD17orD21(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, true) * prefactorSpin;
            const double d18 = doD18orD22(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, true) * prefactorSpin;
            const double d19 = doD19orD23(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, true) * prefactorSpin;
            const double d20 = doD20orD24(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, true) * prefactorSpin;
            setTwoDMA_DMRG(g_index,j_index,k_index,theindex, -2*d17 - 2*d18 - 3*d19 - 3*d20);
            setTwoDMA_DMRG(g_index,j_index,theindex,k_index,  4*d17 + 4*d18                );
            setTwoDMB_DMRG(g_index,j_index,k_index,theindex, -2*d17 - 2*d18 +   d19 +   d20);
            setTwoDMB_DMRG(g_index,j_index,theindex,k_index,                  2*d19 + 2*d20);
            
            //Diagrams 21,22,23 and 24
            const double d21 = doD17orD21(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, false) * prefactorSpin;
            const double d22 = doD18orD22(denT, Ltens[theindex-1][theindex-g_index-1], F0tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, false) * prefactorSpin;
            const double d23 = doD19orD23(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, false) * prefactorSpin;
            const double d24 = doD20orD24(denT, Ltens[theindex-1][theindex-g_index-1], F1tens[theindex][k_index-j_index][j_index-theindex-1],
                                          workmem, workmem2, I_g, false) * prefactorSpin;
            setTwoDMA_DMRG(g_index,k_index,j_index,theindex, -2*d21 - 2*d22 - 3*d23 - 3*d24);
            setTwoDMA_DMRG(g_index,k_index,theindex,j_index,  4*d21 + 4*d22                );
            setTwoDMB_DMRG(g_index,k_index,j_index,theindex, -2*d21 - 2*d22 +   d23 +   d24);
            setTwoDMB_DMRG(g_index,k_index,theindex,j_index,                  2*d23 + 2*d24);
         }
      }

      delete [] workmem;
      delete [] workmem2;
   
   }

}

double CheMPS2::TwoDM::doD1(TensorT * denT){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            int dimL = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimR = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            if ((dimL>0) && (dimR>0)){
            
               double * Tblock = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
               
               int length = dimL*dimR;
               int inc = 1;
               total += (TwoSL+1) * ddot_(&length, Tblock, &inc, Tblock, &inc);
               
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::TwoDM::doD2(TensorT * denT, TensorL * Lright, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSR = TwoSL-1; TwoSR<=TwoSL+1; TwoSR+=2){
            
               int IRup = Irreps::directProd(IL, denBK->gIrrep(theindex)); 
               
               int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
               int dimRdown = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
               int dimRup   = denBK->gCurrentDim(theindex+1,NL+1,TwoSR,IRup);
            
               if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
                  double * Tdown = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL  );
                  double * Tup   = denT->gStorage(NL,TwoSL,IL,NL+1,TwoSR,IRup);
                  double * Lblock = Lright->gStorage(NL+1,TwoSR,IRup,NL+2,TwoSL,IL);
               
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&notrans,&trans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,Lblock,&dimRup,&beta,workmem,&dimL);
                  
                  int fase = ((((TwoSL+1-TwoSR)/2)%2)!=0)?-1:1;
                  double factor = fase * 0.5 * sqrt((TwoSL+1)*(TwoSR+1.0));
                  
                  int length = dimL * dimRup;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

               }
            }
         }
      }
   }
   
   return total;

}

double CheMPS2::TwoDM::doD3(TensorT * denT, TensorS0 * S0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
               
            int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimRdown = denBK->gCurrentDim(theindex+1,NL,  TwoSL,IL);
            int dimRup   = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            
            if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
               double * Tdown = denT->gStorage(NL,TwoSL,IL,NL,  TwoSL,IL);
               double * Tup   = denT->gStorage(NL,TwoSL,IL,NL+2,TwoSL,IL);
               double * S0block = S0right->gStorage(NL, TwoSL, IL, NL+2, TwoSL, IL);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //set
               dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,S0block,&dimRdown,&beta,workmem,&dimL);
               
               double factor = sqrt(0.5) * (TwoSL+1);
               
               int length = dimL * dimRup;
               int inc = 1;
               total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD4(TensorT * denT, TensorF0 * F0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
               
            int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
            int dimR     = denBK->gCurrentDim(theindex+1,NL+2,TwoSL,IL);
            
            if ((dimL>0) && (dimR>0)){
               
               double * Tblock  =    denT->gStorage(NL,   TwoSL, IL, NL+2, TwoSL, IL);
               double * F0block = F0right->gStorage(NL+2, TwoSL, IL, NL+2, TwoSL, IL);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //set
               dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,Tblock,&dimL,F0block,&dimR,&beta,workmem,&dimL);
               
               double factor = sqrt(0.5) * (TwoSL+1);
               
               int length = dimL * dimR;
               int inc = 1;
               total += factor * ddot_(&length, workmem, &inc, Tblock, &inc);

            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD5(TensorT * denT, TensorF0 * F0right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSR=TwoSL-1; TwoSR<=TwoSL+1; TwoSR+=2){
               
               int IR = Irreps::directProd(IL,denBK->gIrrep(theindex));
               int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,IL);
               int dimR     = denBK->gCurrentDim(theindex+1,NL+1,TwoSR,IR);
            
               if ((dimL>0) && (dimR>0)){
               
                  double * Tblock  =    denT->gStorage(NL,   TwoSL, IL, NL+1, TwoSR, IR);
                  double * F0block = F0right->gStorage(NL+1, TwoSR, IR, NL+1, TwoSR, IR);
               
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&notrans,&notrans,&dimL,&dimR,&dimR,&alpha,Tblock,&dimL,F0block,&dimR,&beta,workmem,&dimL);
               
                  double factor = 0.5 * sqrt(0.5) * (TwoSR+1);
                  
                  int length = dimL * dimR;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tblock, &inc);
                  
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD6(TensorT * denT, TensorF1 * F1right, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSL = denBK->gTwoSmin(theindex,NL); TwoSL<= denBK->gTwoSmax(theindex,NL); TwoSL+=2){
         for (int IL = 0; IL<denBK->getNumberOfIrreps(); IL++){
            for (int TwoSRup=TwoSL-1; TwoSRup<=TwoSL+1; TwoSRup+=2){
               for (int TwoSRdown=TwoSL-1; TwoSRdown<=TwoSL+1; TwoSRdown+=2){
               
                  int IR = Irreps::directProd(IL,denBK->gIrrep(theindex));
                  int dimL     = denBK->gCurrentDim(theindex,  NL,  TwoSL,    IL);
                  int dimRup   = denBK->gCurrentDim(theindex+1,NL+1,TwoSRup,  IR);
                  int dimRdown = denBK->gCurrentDim(theindex+1,NL+1,TwoSRdown,IR);
            
                  if ((dimL>0) && (dimRup>0) && (dimRdown>0)){
               
                     double * Tup   =      denT->gStorage(NL,   TwoSL,     IL, NL+1, TwoSRup,   IR);
                     double * Tdown =      denT->gStorage(NL,   TwoSL,     IL, NL+1, TwoSRdown, IR);
                     double * F1block = F1right->gStorage(NL+1, TwoSRdown, IR, NL+1, TwoSRup,   IR);
               
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&notrans,&notrans,&dimL,&dimRup,&dimRdown,&alpha,Tdown,&dimL,F1block,&dimRdown,&beta,workmem,&dimL);
               
                     int fase = ((((TwoSL + TwoSRdown - 1)/2)%2)!=0)?-1:1;
                     double factor = sqrt((TwoSRup+1)/3.0) * (TwoSRdown+1) * fase * gsl_sf_coupling_6j(1,1,2,TwoSRup,TwoSRdown,TwoSL);
                     
                     int length = dimL * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem, &inc, Tup, &inc);
                     
                  }
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD7(TensorT * denT, TensorL * Lleft, double * workmem){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NLup = denBK->gNmin(theindex); NLup<=denBK->gNmax(theindex); NLup++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NLup); TwoSLup<= denBK->gTwoSmax(theindex,NLup); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NLup, TwoSLup, ILup);
         
            for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
               
               int IR = Irreps::directProd(ILup,denBK->gIrrep(theindex));
               int dimLdown = denBK->gCurrentDim(theindex,   NLup-1, TwoSLdown, IR);
               int dimR     = denBK->gCurrentDim(theindex+1, NLup+1, TwoSLdown, IR);
            
               if ((dimLup>0) && (dimLdown>0) && (dimR>0)){
               
                  double * Tup   =   denT->gStorage(NLup,   TwoSLup,   ILup, NLup+1, TwoSLdown, IR);
                  double * Tdown =   denT->gStorage(NLup-1, TwoSLdown, IR,   NLup+1, TwoSLdown, IR);
                  double * Lblock = Lleft->gStorage(NLup-1, TwoSLdown, IR,   NLup,   TwoSLup,   ILup);
               
                  char trans = 'T';
                  char notrans = 'N';
                  double alpha = 1.0;
                  double beta = 0.0; //set
                  dgemm_(&trans,&notrans,&dimLup,&dimR,&dimLdown,&alpha,Lblock,&dimLdown,Tdown,&dimLdown,&beta,workmem,&dimLup);
               
                  int fase = ((((TwoSLup - TwoSLdown + 3)/2)%2)!=0)?-1:1;
                  double factor = 0.5 * sqrt((TwoSLdown+1)*(TwoSLup+1.0)) * fase;
                  
                  int length = dimLup * dimR;
                  int inc = 1;
                  total += factor * ddot_(&length, workmem, &inc, Tup, &inc);

               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::TwoDM::doD8(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g){

   int theindex = denT->gIndex();
   
   double total = 0.0;

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int Idown = Irreps::directProd(ILup,Irrep_g);
                  
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                  
                     double * Tup       =   denT->gStorage(NL,   TwoSLup,   ILup,  NL+2, TwoSLup,   ILup);
                     double * Tdown     =   denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup, ILup);
                     double * LrightBlk = Lright->gStorage(NL+1, TwoSLdown, Idown, NL+2, TwoSLup, ILup);
                  
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //set
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,Tdown,&dimLdown,&beta,workmem,&dimLup);
                  
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRdown,&beta,workmem2,&dimLup);
                  
                     double factor = -0.5 * (TwoSLup+1);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem2, &inc, Tup, &inc);
                  
                  }
               }
            }
         }
      }
   }

   return total;

}

void CheMPS2::TwoDM::doD9andD10andD11(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, double * d9, double * d10, double * d11, int Irrep_g){

   d9[0]  = 0.0;
   d10[0] = 0.0;
   d11[0] = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NL, TwoSLup, ILup);
            if (dimLup>0){
            
               int IRup   = Irreps::directProd(ILup,   denBK->gIrrep(theindex));
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                     for (int TwoSRdown=TwoSRup-1; TwoSRdown<=TwoSRup+1; TwoSRdown+=2){
                        if ((TwoSLdown>=0) && (TwoSRup>=0) && (TwoSRdown>=0) && (abs(TwoSLdown - TwoSRdown)<=1)){
                        
                           int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                           int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                           int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
               
                           if ((dimLdown>0) && (dimRup>0) && (dimRdown>0)){
                           
                              double * T_up      =   denT->gStorage(NL,   TwoSLup,   ILup,   NL+1, TwoSRup,   IRup);
                              double * T_down    =   denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown);
                              double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup,   ILup);
                              double * LrightBlk = Lright->gStorage(NL,   TwoSRdown, IRdown, NL+1, TwoSRup,   IRup);
                              
                              char trans = 'T';
                              char notrans = 'N';
                              double alpha = 1.0;
                              double beta = 0.0; //SET
                              
                              dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                              dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRdown,&beta,workmem2,&dimLup);
                              
                              int length = dimLup * dimRup;
                              int inc = 1;
                              double value = ddot_(&length, workmem2, &inc, T_up, &inc);
                              
                              int fase = ((((TwoSLup + TwoSRdown + 2)/2)%2)!=0) ? -1 : 1;
                              double fact1 = fase * (TwoSRup+1) * sqrt((TwoSRdown+1)*(TwoSLup+1.0)) * gsl_sf_coupling_6j(TwoSRup,1,TwoSLup,TwoSLdown,1,TwoSRdown);
                              double fact2 = 2 * (TwoSRup+1) * sqrt((TwoSRdown+1)*(TwoSLup+1.0)) * gsl_sf_coupling_6j(TwoSRup, TwoSLdown, 2, 1, 1, TwoSLup) * gsl_sf_coupling_6j(TwoSRup, TwoSLdown, 2, 1, 1, TwoSRdown);
                              double fact3 = (TwoSRdown == TwoSLup) ? TwoSRup+1.0 : 0.0 ;
                              
                              d9[0] += fact1  * value;
                              d10[0] += fact2 * value;
                              d11[0] += fact3 * value;
                           
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

double CheMPS2::TwoDM::doD12(TensorT * denT, TensorL * Lleft, TensorL * Lright, double * workmem, double * workmem2, int Irrep_g){

   double d12 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL, TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL, TwoSLup, ILup);
            if ((dimLup>0) && (dimRup>0)){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                        
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                           
                     double * T_up      =   denT->gStorage(NL,   TwoSLup,   ILup,  NL,   TwoSLup,   ILup);
                     double * T_down    =   denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * LleftBlk  =  Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                     double * LrightBlk = Lright->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, Idown);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,LleftBlk,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&trans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,LrightBlk,&dimRup,&beta,workmem2,&dimLup);
               
                     int fase = ((((TwoSLdown+1-TwoSLup)/2)%2)!=0) ? -1 : 1;
                     double factor = fase * 0.5 * sqrt((TwoSLup+1)*(TwoSLdown+1.0));
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     d12 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d12;
   
}

double CheMPS2::TwoDM::doD13(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g){

   double d13 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                        
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL+2, TwoSLup, ILup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup, IRdown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup, ILup);
                     double * S0block = S0right->gStorage(NL,   TwoSLup,   IRdown, NL+2, TwoSLup, ILup);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S0block,&dimRdown,&beta,workmem2,&dimLup);
               
                     double factor = -0.5 * sqrt(0.5) * (TwoSLup+1);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     d13 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d13;
   
}

double CheMPS2::TwoDM::doD14(TensorT * denT, TensorL * Lleft, TensorS0 * S0right, double * workmem, double * workmem2, int Irrep_g){

   double d14 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, IRup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL-1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, IRup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL-1, TwoSLdown, Idown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                     double * S0block = S0right->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, IRup);
                              
                     char trans = 'T';
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S0block,&dimRdown,&beta,workmem2,&dimLup);
               
                     int fase = ((((TwoSLdown + 1 - TwoSLup)/2)%2)!=0) ? -1 : 1;
                     double factor = fase * 0.5 * sqrt(0.5 * (TwoSLup+1) * (TwoSLdown+1));
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     d14 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return d14;
   
}

double CheMPS2::TwoDM::doD15(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g){

   double d15 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            int dimRup = denBK->gCurrentDim(theindex+1, NL+2, TwoSLup, ILup);
            
            if ((dimLup>0) && (dimRup>0)){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRdown=TwoSLdown-1; TwoSRdown<=TwoSLdown+1; TwoSRdown+=2){
                        
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
               
                     if ((dimLdown>0) && (dimRdown>0)){
                           
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL+2, TwoSLup,   ILup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSRdown, IRdown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL,   TwoSLup,   ILup);
                        double * S1block = S1right->gStorage(NL,   TwoSRdown, IRdown, NL+2, TwoSLup,   ILup);
                              
                        char trans = 'T';
                        char notrans = 'N';
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                              
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S1block,&dimRdown,&beta,workmem2,&dimLup);
               
                        int fase = ((((TwoSLdown + TwoSLup + 1)/2)%2)!=0) ? -1 : 1;
                        double factor = fase * (TwoSLup+1) * sqrt((TwoSRdown+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSLup,TwoSRdown,TwoSLdown);
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        d15 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return d15;
   
}

double CheMPS2::TwoDM::doD16(TensorT * denT, TensorL * Lleft, TensorS1 * S1right, double * workmem, double * workmem2, int Irrep_g){

   double d16 = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL-1, TwoSLdown, Idown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSRup,   IRup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL-1, TwoSLdown, Idown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                        double * S1block = S1right->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSRup,   IRup);
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&notrans,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,S1block,&dimRdown,&beta,workmem2,&dimLup);
                  
                        int fase = ((((TwoSRup+TwoSLdown+2)/2)%2)!=0) ? -1 : 1;
                        double factor = fase * (TwoSRup+1) * sqrt((TwoSLup+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSRup,TwoSLdown,TwoSLup);
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        d16 += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return d16;
   
}

double CheMPS2::TwoDM::doD17orD21(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD17){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   ILup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   IRdown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup, ILup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, IRdown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup, ILup);
                     double * F0block = (shouldIdoD17) ? F0right->gStorage(NL,   TwoSLup,   IRdown, NL, TwoSLup, ILup)
                                                       : F0right->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup, IRdown) ;
                              
                     char trans = 'T';
                     char notrans = 'N';
                     char var = (shouldIdoD17) ? notrans : trans;
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     int dimvar = (shouldIdoD17) ? dimRdown : dimRup ;
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F0block,&dimvar,&beta,workmem2,&dimLup);
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += sqrt(0.5) * 0.5 * (TwoSLup+1) * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD18orD22(TensorT * denT, TensorL * Lleft, TensorF0 * F0right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD18){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  
                  int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, IRup);
                  int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                  int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
               
                  if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                           
                     double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSLdown, IRup);
                     double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                     double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL, TwoSLup, ILup);
                     double * F0block = (shouldIdoD18) ? F0right->gStorage(NL+1, TwoSLdown, Idown, NL+1, TwoSLdown, IRup)
                                                       : F0right->gStorage(NL+1, TwoSLdown, IRup,  NL+1, TwoSLdown, Idown) ;
                              
                     char trans = 'T';
                     char notrans = 'N';
                     char var = (shouldIdoD18) ? notrans : trans;
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     int dimvar = (shouldIdoD18) ? dimRdown : dimRup ;
                              
                     dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
               
                     dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F0block,&dimvar,&beta,workmem2,&dimLup);
               
                     int fase = ((((TwoSLdown + 1 - TwoSLup)/2)%2)!=0) ? -1 : 1;
                     double factor = fase * 0.5 * sqrt(0.5*(TwoSLup+1)*(TwoSLdown+1));
                     
                     int length = dimLup * dimRup;
                     int inc = 1;
                     total += factor * ddot_(&length, workmem2, &inc, T_up, &inc);

                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD19orD23(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD19){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex, NL, TwoSLup, ILup);
            
            if (dimLup>0){
            
               int ILdown = Irreps::directProd(ILup,   Irrep_g);
               int IRdown = Irreps::directProd(ILdown, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRdown=TwoSLdown-1; TwoSRdown<=TwoSLdown+1; TwoSRdown+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL,   TwoSLup,   ILup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, ILdown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL,   TwoSRdown, IRdown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,   NL, TwoSLup,   ILup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSRdown, IRdown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, ILdown, NL, TwoSLup,   ILup);
                        double * F1block = (shouldIdoD19) ? F1right->gStorage(NL, TwoSRdown, IRdown, NL, TwoSLup,   ILup)
                                                          : F1right->gStorage(NL, TwoSLup,   ILup,   NL, TwoSRdown, IRdown) ;
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        char var = (shouldIdoD19) ? notrans : trans;
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                        int dimvar = (shouldIdoD19) ? dimRdown : dimRup ;
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F1block,&dimvar,&beta,workmem2,&dimLup);
                  
                        double factor = 0.0;
                        if (shouldIdoD19){
                           int fase = ((((TwoSLdown + TwoSRdown - 1)/2)%2)!=0) ? -1 : 1;
                           factor = fase * (TwoSRdown+1) * sqrt((TwoSLup+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSLup,TwoSRdown,TwoSLdown);
                        } else {
                           int fase = ((((TwoSLdown + TwoSLup - 1)/2)%2)!=0) ? -1 : 1;
                           factor = fase * (TwoSLup+1) * sqrt((TwoSRdown+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSLup,TwoSRdown,TwoSLdown);
                        }
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        total += factor * ddot_(&length,workmem2,&inc,T_up,&inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;
   
}

double CheMPS2::TwoDM::doD20orD24(TensorT * denT, TensorL * Lleft, TensorF1 * F1right, double * workmem, double * workmem2, int Irrep_g, bool shouldIdoD20){

   double total = 0.0;

   int theindex = denT->gIndex();

   for (int NL = denBK->gNmin(theindex); NL<=denBK->gNmax(theindex); NL++){
      for (int TwoSLup = denBK->gTwoSmin(theindex,NL); TwoSLup<= denBK->gTwoSmax(theindex,NL); TwoSLup+=2){
         for (int ILup = 0; ILup<denBK->getNumberOfIrreps(); ILup++){
         
            int dimLup = denBK->gCurrentDim(theindex,   NL,   TwoSLup, ILup);
            
            if (dimLup>0){
            
               int Idown = Irreps::directProd(ILup, Irrep_g);
               int IRup  = Irreps::directProd(ILup, denBK->gIrrep(theindex));
         
               for (int TwoSLdown=TwoSLup-1; TwoSLdown<=TwoSLup+1; TwoSLdown+=2){
                  for (int TwoSRup=TwoSLup-1; TwoSRup<=TwoSLup+1; TwoSRup+=2){
                  
                     int dimRup   = denBK->gCurrentDim(theindex+1, NL+1, TwoSRup,   IRup);
                     int dimLdown = denBK->gCurrentDim(theindex,   NL-1, TwoSLdown, Idown);
                     int dimRdown = denBK->gCurrentDim(theindex+1, NL+1, TwoSLdown, Idown);
                  
                     if ((dimLdown>0) && (dimRdown>0) && (dimRup>0)){
                              
                        double * T_up    =    denT->gStorage(NL,   TwoSLup,   ILup,  NL+1, TwoSRup,   IRup);
                        double * T_down  =    denT->gStorage(NL-1, TwoSLdown, Idown, NL+1, TwoSLdown, Idown);
                        double * Lblock  =   Lleft->gStorage(NL-1, TwoSLdown, Idown, NL,   TwoSLup,   ILup);
                        double * F1block = (shouldIdoD20) ? F1right->gStorage(NL+1, TwoSLdown, Idown, NL+1, TwoSRup,   IRup)
                                                          : F1right->gStorage(NL+1, TwoSRup,   IRup,  NL+1, TwoSLdown, Idown) ;
                                 
                        char trans = 'T';
                        char notrans = 'N';
                        char var = (shouldIdoD20) ? notrans : trans;
                        double alpha = 1.0;
                        double beta = 0.0; //SET
                        int dimvar = (shouldIdoD20) ? dimRdown : dimRup ;
                                 
                        dgemm_(&trans,&notrans,&dimLup,&dimRdown,&dimLdown,&alpha,Lblock,&dimLdown,T_down,&dimLdown,&beta,workmem,&dimLup);
                  
                        dgemm_(&notrans,&var,&dimLup,&dimRup,&dimRdown,&alpha,workmem,&dimLup,F1block,&dimvar,&beta,workmem2,&dimLup);
                  
                        double factor = 0.0;
                        if (shouldIdoD20){
                           int fase = (((TwoSLup)%2)!=0) ? -1 : 1;
                           factor = fase * sqrt((TwoSLup+1)*(TwoSRup+1)*(TwoSLdown+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSRup,TwoSLdown,TwoSLup);
                        } else {
                           int fase = ((((2*TwoSLup + TwoSRup - TwoSLdown)/2)%2)!=0) ? -1 : 1;
                           factor = fase * (TwoSRup+1) * sqrt((TwoSLup+1)/3.0) * gsl_sf_coupling_6j(1,1,2,TwoSRup,TwoSLdown,TwoSLup);
                        }
                        
                        int length = dimLup * dimRup;
                        int inc = 1;
                        total += factor * ddot_(&length, workmem2, &inc, T_up, &inc);
                     
                     }
                  }
               }
            }
         }
      }
   }
   
   return total;
   
}




