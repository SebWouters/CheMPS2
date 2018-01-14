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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "Correlations.h"
#include "Lapack.h"
#include "MPIchemps2.h"

using std::cout;
using std::endl;
using std::max;
using std::min;

CheMPS2::Correlations::Correlations(const SyBookkeeper * denBKIn, const Problem * ProbIn, TwoDM * the2DMin){

   denBK = denBKIn;
   Prob = ProbIn;
   the2DM = the2DMin;
   
   L = denBK->gL();
   
   Cspin     = new double[L*L];
   Cdens     = new double[L*L];
   Cspinflip = new double[L*L];
   Cdirad    = new double[L*L];
   MutInfo   = new double[L*L];
   
   for (int cnt=0; cnt<L*L; cnt++){ Cspin[cnt]     = 0.0; }
   for (int cnt=0; cnt<L*L; cnt++){ Cdens[cnt]     = 0.0; }
   for (int cnt=0; cnt<L*L; cnt++){ Cspinflip[cnt] = 0.0; }
   for (int cnt=0; cnt<L*L; cnt++){ Cdirad[cnt]    = 0.0; }
   for (int cnt=0; cnt<L*L; cnt++){ MutInfo[cnt]   = 0.0; }
   
   FillSpinDensSpinflip();

}

CheMPS2::Correlations::~Correlations(){

   delete [] Cspin;
   delete [] Cdens;
   delete [] Cspinflip;
   delete [] Cdirad;
   delete [] MutInfo;

}

void CheMPS2::Correlations::FillSpinDensSpinflip(){
   
   //Spin
   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         Cspin[row + L*col] = the2DM->getTwoDMB_DMRG(row,col,row,col);
      }
      Cspin[row + L*row] += the2DM->get1RDM_DMRG(row, row);
   }
   
   //Density
   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         Cdens[row + L*col] = the2DM->getTwoDMA_DMRG(row,col,row,col) - the2DM->get1RDM_DMRG(row, row) * the2DM->get1RDM_DMRG(col, col);
      }
      Cdens[row + L*row] += the2DM->get1RDM_DMRG(row, row);
   }
   
   //Spin-flip
   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         Cspinflip[row + L*col] = 0.5 * ( the2DM->getTwoDMB_DMRG(row,col,col,row) - the2DM->getTwoDMA_DMRG(row,col,col,row) );
      }
      Cspinflip[row + L*row] += the2DM->get1RDM_DMRG(row, row);
   }
   
   //Singlet diradical: partial fill
   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         Cdirad[row + L*col] = - 0.5 * ( the2DM->get1RDM_DMRG(row, row) - the2DM->getTwoDMA_DMRG(row,row,row,row) )
                                     * ( the2DM->get1RDM_DMRG(col, col) - the2DM->getTwoDMA_DMRG(col,col,col,col) );
      }
   }

}

double CheMPS2::Correlations::getCspin_DMRG(const int row, const int col) const{ return Cspin[row + L*col]; }

double CheMPS2::Correlations::getCspin_HAM(const int row, const int col) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getCspin_DMRG( Prob->gf1(row), Prob->gf1(col) );
   }
   return getCspin_DMRG( row, col );

}

double CheMPS2::Correlations::getCdens_DMRG(const int row, const int col) const{ return Cdens[row + L*col]; }

double CheMPS2::Correlations::getCdens_HAM(const int row, const int col) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getCdens_DMRG( Prob->gf1(row), Prob->gf1(col) );
   }
   return getCdens_DMRG( row, col );

}

double CheMPS2::Correlations::getCspinflip_DMRG(const int row, const int col) const{ return Cspinflip[row + L*col]; }

double CheMPS2::Correlations::getCspinflip_HAM(const int row, const int col) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getCspinflip_DMRG( Prob->gf1(row), Prob->gf1(col) );
   }
   return getCspinflip_DMRG( row, col );

}

double CheMPS2::Correlations::getCdirad_DMRG(const int row, const int col) const{ return Cdirad[row + L*col]; }

double CheMPS2::Correlations::getCdirad_HAM(const int row, const int col) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getCdirad_DMRG( Prob->gf1(row), Prob->gf1(col) );
   }
   return getCdirad_DMRG( row, col );

}

double CheMPS2::Correlations::getMutualInformation_DMRG(const int row, const int col) const{ return MutInfo[row + L*col]; }

double CheMPS2::Correlations::getMutualInformation_HAM(const int row, const int col) const{

   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   if ( Prob->gReorder() ){
      return getMutualInformation_DMRG( Prob->gf1(row), Prob->gf1(col) );
   }
   return getMutualInformation_DMRG( row, col );

}

double CheMPS2::Correlations::SingleOrbitalEntropy_DMRG(const int index) const{

   const double val4  = 0.5 * the2DM->getTwoDMA_DMRG(index,index,index,index);
   const double val23 = 0.5 * ( the2DM->get1RDM_DMRG(index, index) - the2DM->getTwoDMA_DMRG(index,index,index,index) );
   const double val1  = 1.0 - val4 - 2*val23;
   
   double entropy = 0.0;
   if (val1  > CheMPS2::CORRELATIONS_discardEig){ entropy -=     val1  * log(val1 ); }
   if (val23 > CheMPS2::CORRELATIONS_discardEig){ entropy -= 2 * val23 * log(val23); }
   if (val4  > CheMPS2::CORRELATIONS_discardEig){ entropy -=     val4  * log(val4 ); }
   return entropy;

}

double CheMPS2::Correlations::SingleOrbitalEntropy_HAM(const int index) const{

   if ( Prob->gReorder() ){
      return SingleOrbitalEntropy_DMRG( Prob->gf1(index) );
   }
   return SingleOrbitalEntropy_DMRG( index );

}

double CheMPS2::Correlations::MutualInformationDistance(const double power) const{

   double Idist = 0.0;
   for (int row=0; row<L; row++){
      for (int col=0; col<L; col++){
         if (row != col){
            Idist += MutInfo[row + L*col] * pow(abs(row - col), power);
         }
      }
   }
   return Idist;

}

#ifdef CHEMPS2_MPI_COMPILATION
void CheMPS2::Correlations::mpi_broadcast(){

   int size = L*L;
   MPIchemps2::broadcast_array_double( MutInfo, size, MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_array_double( Cdirad,  size, MPI_CHEMPS2_MASTER );

}
#endif

void CheMPS2::Correlations::FillSite(TensorT * denT, TensorGYZ ** Gtensors, TensorGYZ ** Ytensors, TensorGYZ ** Ztensors, TensorKM ** Ktensors, TensorKM ** Mtensors){

   const int theindex = denT->gIndex();
   const int MAXDIM = max(denBK->gMaxDimAtBound(theindex), denBK->gMaxDimAtBound(theindex+1));
   double * workmem  = new double[MAXDIM*MAXDIM];
   int lindimRDM = 16;
   double * RDM = new double[lindimRDM*lindimRDM];
   int lwork = 3*lindimRDM - 1;
   double * work = new double[lwork];
   double * eigs = new double[lindimRDM];
   
   const double prefactorSpin = 1.0/(Prob->gTwoS() + 1.0);
   const double sqrt_one_half = sqrt(0.5);
   
   for (int previousindex=0; previousindex<theindex; previousindex++){
   
      const bool equalIrreps = (denBK->gIrrep(previousindex) == denBK->gIrrep(theindex)) ? true : false;
      const double diag1  = diagram3(denT, Gtensors[previousindex], workmem) * prefactorSpin * 0.5 * sqrt_one_half;
      const double diag2  = 0.125 * (   the2DM->getTwoDMB_DMRG(previousindex,theindex,theindex,previousindex)
                                      - the2DM->getTwoDMA_DMRG(previousindex,theindex,theindex,previousindex) );
      
      const double val1 = diagram1(denT, Ytensors[previousindex], workmem) * prefactorSpin;                                  //1x1 block N=0, Sz=0
      const double val2 = diagram2(denT, Ztensors[previousindex], workmem) * prefactorSpin;                                  //1x1 block N=4, Sz=0
      const double val3 = diag1 + diag2;                                                                                     //1x1 block N=2, Sz=2*sigma
      const double val4 = diagram1(denT, Gtensors[previousindex], workmem) * prefactorSpin * sqrt_one_half;                  //2x2 block N=1, alpha_LL
      const double val5 = diagram3(denT, Ytensors[previousindex], workmem) * prefactorSpin * 0.5;                            //2x2 block N=1, alpha_RR
      const double val6 = (equalIrreps) ? ( diagram4(denT, Ktensors[previousindex], workmem) * prefactorSpin * 0.5 ) : 0.0 ; //2x2 block N=1, alpha_LR
      const double val7 = diagram2(denT, Gtensors[previousindex], workmem) * prefactorSpin * sqrt_one_half;                  //2x2 block N=3, alpha_LL
      const double val8 = diagram3(denT, Ztensors[previousindex], workmem) * prefactorSpin * 0.5;                            //2x2 block N=3, alpha_RR
      const double val9 = (equalIrreps) ? ( diagram5(denT, Mtensors[previousindex], workmem) * prefactorSpin * 0.5 ) : 0.0 ; //2x2 block N=3, alpha_LR
      
      //4x4 block N=2, Sz=0
      const double alpha   = diagram2(denT, Ytensors[previousindex], workmem) * prefactorSpin;
      const double gamma   = diagram1(denT, Ztensors[previousindex], workmem) * prefactorSpin;
      const double beta    = diag1 - diag2;
      const double lambda  = 2*diag2;
      const double delta   = (equalIrreps) ? ( - diagram5(denT, Ktensors[previousindex], workmem) * prefactorSpin * 0.5 ) : 0.0;
      const double epsilon = (equalIrreps) ? (   diagram4(denT, Mtensors[previousindex], workmem) * prefactorSpin * 0.5 ) : 0.0;
      const double kappa   = 0.5 * the2DM->getTwoDMA_DMRG(previousindex,previousindex,theindex,theindex);
      
      /*
      
         [ val1                                                                                                     ]
         [       val4  val6                                                                                         ]
         [       val6  val5                                                                                         ]
         [                   val4  val6                                                                             ]
         [                   val6  val5                                                                             ]
         [                               val3                                                                       ]
         [                                     alpha  delta   -delta    kappa                                       ]
         [                                     delta  beta     lambda   epsilon                                     ]
         [                                    -delta  lambda   beta    -epsilon                                     ]
         [                                     kappa  epsilon -epsilon  gamma                                       ]
         [                                                                       val3                               ]
         [                                                                             val7  val9                   ]
         [                                                                             val9  val8                   ]
         [                                                                                         val7  val9       ]
         [                                                                                         val9  val8       ]
         [                                                                                                     val2 ]
      
      */
      
      for (int cnt=0; cnt<lindimRDM*lindimRDM; cnt++){ RDM[cnt] = 0.0; }
      RDM[0  + lindimRDM * 0 ] = val1;
      RDM[15 + lindimRDM * 15] = val2;
      RDM[5  + lindimRDM * 5 ] = RDM[10 + lindimRDM * 10] = val3;
      RDM[1  + lindimRDM * 1 ] = RDM[3  + lindimRDM * 3 ] = val4;
      RDM[2  + lindimRDM * 2 ] = RDM[4  + lindimRDM * 4 ] = val5;
      RDM[1  + lindimRDM * 2 ] = RDM[2  + lindimRDM * 1 ] = RDM[3  + lindimRDM * 4 ] = RDM[4  + lindimRDM * 3 ] = val6;
      RDM[11 + lindimRDM * 11] = RDM[13 + lindimRDM * 13] = val7;
      RDM[12 + lindimRDM * 12] = RDM[14 + lindimRDM * 14] = val8;
      RDM[11 + lindimRDM * 12] = RDM[12 + lindimRDM * 11] = RDM[13 + lindimRDM*14] = RDM[14 + lindimRDM*13] = val9;
      RDM[6  + lindimRDM * 6 ] = alpha;
      RDM[7  + lindimRDM * 7 ] = RDM[8  + lindimRDM * 8 ] = beta;
      RDM[9  + lindimRDM * 9 ] = gamma;
      RDM[6  + lindimRDM * 7 ] = RDM[7  + lindimRDM * 6 ] = delta;
      RDM[6  + lindimRDM * 8 ] = RDM[8  + lindimRDM * 6 ] = - delta;
      RDM[7  + lindimRDM * 9 ] = RDM[9  + lindimRDM * 7 ] = epsilon;
      RDM[8  + lindimRDM * 9 ] = RDM[9  + lindimRDM * 8 ] = - epsilon;
      RDM[6  + lindimRDM * 9 ] = RDM[9  + lindimRDM * 6 ] = kappa;
      RDM[7  + lindimRDM * 8 ] = RDM[8  + lindimRDM * 7 ] = lambda;
      
      if (CheMPS2::CORRELATIONS_debugPrint){
         const double RDM_1orb_prev_4  = 0.5 * the2DM->getTwoDMA_DMRG(previousindex,previousindex,previousindex,previousindex);
         const double RDM_1orb_prev_23 = 0.5 * (   the2DM->get1RDM_DMRG(previousindex, previousindex)
                                                 - the2DM->getTwoDMA_DMRG(previousindex,previousindex,previousindex,previousindex) );
         const double RDM_1orb_prev_1  = 1.0 - RDM_1orb_prev_4 - 2*RDM_1orb_prev_23;
         
         const double RDM_1orb_curr_4  = 0.5 * the2DM->getTwoDMA_DMRG(theindex,theindex,theindex,theindex);
         const double RDM_1orb_curr_23 = 0.5 * ( the2DM->get1RDM_DMRG(theindex, theindex) - the2DM->getTwoDMA_DMRG(theindex,theindex,theindex,theindex) );
         const double RDM_1orb_curr_1  = 1.0 - RDM_1orb_curr_4 - 2*RDM_1orb_curr_23;
         
         cout << "   Correlations::FillSite : Looking at DMRG sites (" << previousindex << "," << theindex << ")." << endl;
         
         //Check 1 : full trace
         double fulltrace = 0.0;
         for (int cnt=0; cnt<lindimRDM; cnt++){ fulltrace += RDM[ cnt * ( 1 + lindimRDM ) ]; }
         cout << "                            1 - trace(2-orb RDM) = " << 1.0 - fulltrace << endl;
         
         //Check 2 : trace over orb previousindex
         double getal1 = RDM[0  + lindimRDM * 0 ] + RDM[1  + lindimRDM * 1 ] + RDM[3  + lindimRDM * 3 ] + RDM[9  + lindimRDM * 9 ] - RDM_1orb_curr_1;
         double getal2 = RDM[2  + lindimRDM * 2 ] + RDM[5  + lindimRDM * 5 ] + RDM[8  + lindimRDM * 8 ] + RDM[12 + lindimRDM * 12] - RDM_1orb_curr_23;
         double getal3 = RDM[4  + lindimRDM * 4 ] + RDM[7  + lindimRDM * 7 ] + RDM[10 + lindimRDM * 10] + RDM[14 + lindimRDM * 14] - RDM_1orb_curr_23;
         double getal4 = RDM[6  + lindimRDM * 6 ] + RDM[11 + lindimRDM * 11] + RDM[13 + lindimRDM * 13] + RDM[15 + lindimRDM * 15] - RDM_1orb_curr_4;
         double RMS = sqrt(getal1*getal1 + getal2*getal2 + getal3*getal3 + getal4*getal4);
         cout << "                            2-norm difference of one-orb RDM of the CURRENT site via 2DM and via trace(2-orb RDM) = " << RMS << endl;
         
         //Check 3 : trace over orb currentindex
         getal1 = RDM[0  + lindimRDM * 0 ] + RDM[2  + lindimRDM * 2 ] + RDM[4  + lindimRDM * 4 ] + RDM[6  + lindimRDM * 6 ] - RDM_1orb_prev_1;
         getal2 = RDM[1  + lindimRDM * 1 ] + RDM[5  + lindimRDM * 5 ] + RDM[7  + lindimRDM * 7 ] + RDM[11 + lindimRDM * 11] - RDM_1orb_prev_23;
         getal3 = RDM[3  + lindimRDM * 3 ] + RDM[8  + lindimRDM * 8 ] + RDM[10 + lindimRDM * 10] + RDM[13 + lindimRDM * 13] - RDM_1orb_prev_23;
         getal4 = RDM[9  + lindimRDM * 9 ] + RDM[12 + lindimRDM * 12] + RDM[14 + lindimRDM * 14] + RDM[15 + lindimRDM * 15] - RDM_1orb_prev_4;
         RMS = sqrt(getal1*getal1 + getal2*getal2 + getal3*getal3 + getal4*getal4);
         cout << "                            2-norm difference of one-orb RDM of the PREVIOUS site via 2DM and via trace(2-orb RDM) = " << RMS << endl;
      }
      
      char jobz = 'N'; //eigenvalues only
      char uplo = 'U';
      int info;
      dsyev_(&jobz, &uplo, &lindimRDM, RDM, &lindimRDM, eigs, work, &lwork, &info);
      
      double entropy = 0.0;
      for (int cnt=0; cnt<lindimRDM; cnt++){
         if (eigs[cnt] > CheMPS2::CORRELATIONS_discardEig){ //With discardEig = 1e-100, the discarded contributions are smaller than 2.4e-98
            entropy -= eigs[cnt] * log(eigs[cnt]);
         }
      }
      
      const double thisMutInfo = 0.5 * ( SingleOrbitalEntropy_DMRG(previousindex) + SingleOrbitalEntropy_DMRG(theindex) - entropy );
      MutInfo[previousindex + L * theindex] = MutInfo[theindex + L * previousindex] = thisMutInfo;
      Cdirad[previousindex + L * theindex] += 2 * beta;
      Cdirad[theindex + L * previousindex] += 2 * beta;
      
   }
   
   delete [] eigs;
   delete [] work;
   delete [] RDM;
   delete [] workmem;

}

double CheMPS2::Correlations::diagram1(TensorT * denT, TensorGYZ * denY, double * workmem) const{ //checked

   const int theindex = denT->gIndex();
   
   double total = 0.0;
   
   for (int NR = denBK->gNmin(theindex+1); NR <= denBK->gNmax(theindex+1); NR++){
      for (int TwoSR = denBK->gTwoSmin(theindex+1,NR); TwoSR <= denBK->gTwoSmax(theindex+1,NR); TwoSR += 2){
         for (int IR = 0; IR < denBK->getNumberOfIrreps(); IR++){
         
            int dimL = denBK->gCurrentDim(theindex,   NR, TwoSR, IR);
            int dimR = denBK->gCurrentDim(theindex+1, NR, TwoSR, IR);
            
            if ((dimL>0)&&(dimR>0)){
            
               double * blockT = denT->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               double * blockY = denY->gStorage(NR,TwoSR,IR,NR,TwoSR,IR);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //SET
               dgemm_(&notrans, &notrans, &dimL, &dimR, &dimL, &alpha, blockY, &dimL, blockT, &dimL, &beta, workmem, &dimL);
               
               int length = dimL*dimR;
               int inc = 1;
               total += (TwoSR + 1.0) * ddot_(&length, workmem, &inc, blockT, &inc);
            
            }
         }
      }
   }

   return total;

}

double CheMPS2::Correlations::diagram2(TensorT * denT, TensorGYZ * denZ, double * workmem) const{ //checked

   const int theindex = denT->gIndex();
   
   double total = 0.0;
   
   for (int NR = denBK->gNmin(theindex+1); NR <= denBK->gNmax(theindex+1); NR++){
      for (int TwoSR = denBK->gTwoSmin(theindex+1,NR); TwoSR <= denBK->gTwoSmax(theindex+1,NR); TwoSR += 2){
         for (int IR = 0; IR < denBK->getNumberOfIrreps(); IR++){
         
            int dimL = denBK->gCurrentDim(theindex,   NR-2, TwoSR, IR);
            int dimR = denBK->gCurrentDim(theindex+1, NR,   TwoSR, IR);
            
            if ((dimL>0)&&(dimR>0)){
            
               double * blockT = denT->gStorage(NR-2,TwoSR,IR,NR,  TwoSR,IR);
               double * blockZ = denZ->gStorage(NR-2,TwoSR,IR,NR-2,TwoSR,IR);
               
               char notrans = 'N';
               double alpha = 1.0;
               double beta = 0.0; //SET
               dgemm_(&notrans, &notrans, &dimL, &dimR, &dimL, &alpha, blockZ, &dimL, blockT, &dimL, &beta, workmem, &dimL);
               
               int length = dimL*dimR;
               int inc = 1;
               total += (TwoSR + 1.0) * ddot_(&length, workmem, &inc, blockT, &inc);
            
            }
         }
      }
   }

   return total;

}

double CheMPS2::Correlations::diagram3(TensorT * denT, TensorGYZ * denG, double * workmem) const{ //checked

   const int theindex = denT->gIndex();
   
   double total = 0.0;
   
   for (int NR = denBK->gNmin(theindex+1); NR <= denBK->gNmax(theindex+1); NR++){
      for (int TwoSR = denBK->gTwoSmin(theindex+1,NR); TwoSR <= denBK->gTwoSmax(theindex+1,NR); TwoSR += 2){
         for (int IR = 0; IR < denBK->getNumberOfIrreps(); IR++){
         
            int dimR = denBK->gCurrentDim(theindex+1, NR, TwoSR, IR);
            const int IL = Irreps::directProd( IR, denBK->gIrrep(theindex) );
            
            if (dimR>0){
            
               for (int TwoSL = TwoSR-1; TwoSL <= TwoSR+1; TwoSL+=2){
               
                  int dimL = denBK->gCurrentDim(theindex, NR-1, TwoSL, IL);
            
                  if (dimL>0){
                  
                     double * blockT = denT->gStorage(NR-1,TwoSL,IL,NR,  TwoSR,IR);
                     double * blockG = denG->gStorage(NR-1,TwoSL,IL,NR-1,TwoSL,IL);
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     dgemm_(&notrans, &notrans, &dimL, &dimR, &dimL, &alpha, blockG, &dimL, blockT, &dimL, &beta, workmem, &dimL);
                     
                     int length = dimL*dimR;
                     int inc = 1;
                     total += (TwoSR + 1.0) * ddot_(&length, workmem, &inc, blockT, &inc);
                  
                  }
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::Correlations::diagram4(TensorT * denT, TensorKM * denK, double * workmem) const{ //checked

   const int theindex = denT->gIndex();
   
   double total = 0.0;
   
   for (int NR = denBK->gNmin(theindex+1); NR <= denBK->gNmax(theindex+1); NR++){
      for (int TwoSR = denBK->gTwoSmin(theindex+1,NR); TwoSR <= denBK->gTwoSmax(theindex+1,NR); TwoSR += 2){
         for (int IR = 0; IR < denBK->getNumberOfIrreps(); IR++){
         
            int dimR   = denBK->gCurrentDim(theindex+1, NR, TwoSR, IR);
            int dimLup = denBK->gCurrentDim(theindex,   NR, TwoSR, IR);
            const int ILdown = Irreps::directProd( IR, denBK->gIrrep(theindex) );
            
            if ((dimR>0) && (dimLup>0)){
            
               for (int TwoSLdown = TwoSR-1; TwoSLdown <= TwoSR+1; TwoSLdown += 2){
               
                  int dimLdown = denBK->gCurrentDim(theindex, NR-1, TwoSLdown, ILdown);
            
                  if (dimLdown>0){
                  
                     double * blockTup   = denT->gStorage(NR,   TwoSR,     IR,     NR, TwoSR, IR);
                     double * blockTdown = denT->gStorage(NR-1, TwoSLdown, ILdown, NR, TwoSR, IR);
                     double * blockK     = denK->gStorage(NR-1, TwoSLdown, ILdown, NR, TwoSR, IR);
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     dgemm_(&notrans, &notrans, &dimLdown, &dimR, &dimLup, &alpha, blockK, &dimLdown, blockTup, &dimLup, &beta, workmem, &dimLdown);
                     
                     int length = dimLdown*dimR;
                     int inc = 1;
                     total += (TwoSR + 1.0) * ddot_(&length, workmem, &inc, blockTdown, &inc);
                  
                  }
               }
            }
         }
      }
   }

   return total;

}

double CheMPS2::Correlations::diagram5(TensorT * denT, TensorKM * denM, double * workmem) const{ //checked

   const int theindex = denT->gIndex();
   
   double total = 0.0;
   
   for (int NR = denBK->gNmin(theindex+1); NR <= denBK->gNmax(theindex+1); NR++){
      for (int TwoSR = denBK->gTwoSmin(theindex+1,NR); TwoSR <= denBK->gTwoSmax(theindex+1,NR); TwoSR += 2){
         for (int IR = 0; IR < denBK->getNumberOfIrreps(); IR++){
         
            int dimR     = denBK->gCurrentDim(theindex+1, NR,   TwoSR, IR);
            int dimLdown = denBK->gCurrentDim(theindex,   NR-2, TwoSR, IR);
            const int ILup = Irreps::directProd( IR, denBK->gIrrep(theindex) );
            
            if ((dimR>0) && (dimLdown>0)){
            
               for (int TwoSLup = TwoSR-1; TwoSLup <= TwoSR+1; TwoSLup += 2){
               
                  int dimLup = denBK->gCurrentDim(theindex, NR-1, TwoSLup, ILup);
            
                  if (dimLup>0){
                  
                     double * blockTdown = denT->gStorage(NR-2, TwoSR,     IR,   NR,   TwoSR,   IR);
                     double * blockTup   = denT->gStorage(NR-1, TwoSLup,   ILup, NR,   TwoSR,   IR);
                     double * blockM     = denM->gStorage(NR-2, TwoSR,     IR,   NR-1, TwoSLup, ILup);
                     
                     char notrans = 'N';
                     double alpha = 1.0;
                     double beta = 0.0; //SET
                     dgemm_(&notrans, &notrans, &dimLdown, &dimR, &dimLup, &alpha, blockM, &dimLdown, blockTup, &dimLup, &beta, workmem, &dimLdown);
                     
                     int length = dimLdown*dimR;
                     int inc = 1;
                     const int fase = ((((TwoSLup + 1 - TwoSR)/2)%2)!=0)?-1:1;
                     total += fase * sqrt((TwoSLup+1.0)*(TwoSR+1)) * ddot_(&length, workmem, &inc, blockTdown, &inc);
                  
                  }
               }
            }
         }
      }
   }

   return total;

}

void CheMPS2::Correlations::PrintTableNice(const double * table, const int sPrecision, const int columnsPerLine) const{

   std::stringstream thestream;
   thestream.precision(sPrecision);
   thestream << std::fixed;
   
   int numGroups = L / columnsPerLine;
   if (numGroups * columnsPerLine < L){ numGroups++; }
   
   std::string prefix = "   ";
   
   for (int groups=0; groups<numGroups; groups++){
      const int startCol = groups * columnsPerLine + 1;
      const int stopCol  = min( (groups + 1) * columnsPerLine , L );
      thestream << prefix << "Columns " << startCol << " to " << stopCol << "\n \n";
      for (int row=0; row<L; row++){
         for (int col=startCol-1; col<stopCol; col++){
            if ((row==col) && (table==MutInfo)){
               thestream << prefix << " 0 ";
               for (int cnt=0; cnt<sPrecision; cnt++){ thestream << " "; }
            } else {
               int row2 = (Prob->gReorder()) ? Prob->gf1(row) : row;
               int col2 = (Prob->gReorder()) ? Prob->gf1(col) : col;
               if (table[row2 + L*col2] < 0.0){
                  thestream << prefix        << table[row2 + L*col2];
               } else {
                  thestream << prefix << " " << table[row2 + L*col2];
               }
            }
         }
         thestream << "\n";
      }
      thestream << "\n";
   
   }
   
   cout << thestream.str();

}

void CheMPS2::Correlations::Print(const int precision, const int columnsPerLine) const{
   
   cout << "--------------------------------------------------------" << endl;
   cout << "Spin correlation function = 4 * ( < S_i^z S_j^z > - < S_i^z > * < S_j^z > ) \nHamiltonian index order is used!\n" << endl;
   PrintTableNice( Cspin , precision, columnsPerLine );
   cout << "--------------------------------------------------------" << endl;
   cout << "Spin-flip correlation function = < S_i^+ S_j^- > + < S_i^- S_j^+ > \nHamiltonian index order is used!\n" << endl;
   PrintTableNice( Cspinflip , precision, columnsPerLine);
   cout << "--------------------------------------------------------" << endl;
   cout << "Density correlation function = < n_i n_j > - < n_i > * < n_j > \nHamiltonian index order is used!\n" << endl;
   PrintTableNice( Cdens , precision, columnsPerLine );
   cout << "--------------------------------------------------------" << endl;
   cout << "Singlet diradical correlation function = < d_i,up d_j,down > + < d_i,down d_j,up > - < d_i,up > * < d_j,down > - < d_i,down > * < d_j,up > \nHamiltonian index order is used!\n" << endl;
   PrintTableNice( Cdirad , precision, columnsPerLine );
   cout << "--------------------------------------------------------" << endl;
   cout << "Two-orbital mutual information = 0.5 * ( s1(i) + s1(j) - s2(i,j) ) * ( 1 - delta(i,j) ) \nHamiltonian index order is used!\n" << endl;
   PrintTableNice( MutInfo , precision, columnsPerLine );
   cout << "--------------------------------------------------------" << endl;

}


