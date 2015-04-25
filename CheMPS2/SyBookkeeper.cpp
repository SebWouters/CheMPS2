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

#include <algorithm>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <math.h>

#include "SyBookkeeper.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

CheMPS2::SyBookkeeper::SyBookkeeper(const Problem * Probin, const int Din) : Irreps(Probin->gSy()){

   Prob = Probin;
   
   //Set the min and max particle number and spin
   Nmin = new int[gL()+1];
   Nmax = new int[gL()+1];
   TwoSmin = new int*[gL()+1];
   TwoSmax = new int*[gL()+1];
   for (int bound=0; bound<=gL(); bound++){
      Nmin[bound] = max( max(0, gN()+2*(bound - gL())) , bound - gL() + (gN() + gTwoS())/2);
      Nmax[bound] = min( min(2*bound, gN() ), bound + (gN() - gTwoS())/2);
      TwoSmin[bound] = new int[Nmax[bound]-Nmin[bound]+1];
      TwoSmax[bound] = new int[Nmax[bound]-Nmin[bound]+1];
      for (int N=Nmin[bound]; N<=Nmax[bound]; N++){
         TwoSmin[bound][N-Nmin[bound]] = max(N%2, gTwoS() - (gL() - bound - abs(gN() - N - gL() + bound)));
         TwoSmax[bound][N-Nmin[bound]] = min( bound - abs(bound - N), gTwoS() + (gL() - bound - abs(gN() - N - gL() + bound)));
      }
   }
   
   //FCIdim & CurrentDim memory allocation
   FCIdim = new int***[gL()+1];
   CurrentDim = new int***[gL()+1];
   for (int bound=0; bound<=gL(); bound++){
      FCIdim[bound] = new int**[gNmax(bound)-gNmin(bound)+1];
      CurrentDim[bound] = new int**[gNmax(bound)-gNmin(bound)+1];
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         FCIdim[bound][N-gNmin(bound)] = new int*[(gTwoSmax(bound,N)-gTwoSmin(bound,N))/2+1];
         CurrentDim[bound][N-gNmin(bound)] = new int*[(gTwoSmax(bound,N)-gTwoSmin(bound,N))/2+1];
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            FCIdim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2] = new int[getNumberOfIrreps()];
            CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2] = new int[getNumberOfIrreps()];
         }
      }
   }
   
   //Fill the FCI dims & copy it to the Current dims
   fillFCIdim();
   
   //Scale the Current dims
   ScaleCurrentDim(Din);
   
   if (CheMPS2::SYBK_debugPrint) print();
   assert( IsPossible() );

}

CheMPS2::SyBookkeeper::~SyBookkeeper(){

   for (int bound=0; bound<=gL(); bound++){
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            delete [] FCIdim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2];
            delete [] CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2];
         }
         delete [] FCIdim[bound][N-gNmin(bound)];
         delete [] CurrentDim[bound][N-gNmin(bound)];
      }
      delete [] FCIdim[bound];
      delete [] CurrentDim[bound];
   }
   delete [] FCIdim;
   delete [] CurrentDim;

   for (int bound=0; bound<=gL(); bound++){
      delete [] TwoSmin[bound];
      delete [] TwoSmax[bound];
   }
   delete [] TwoSmin;
   delete [] TwoSmax;
   delete [] Nmin;
   delete [] Nmax;

}

int CheMPS2::SyBookkeeper::gL() const{ return Prob->gL(); }

int CheMPS2::SyBookkeeper::gIrrep(const int nOrb) const{ return Prob->gIrrep(nOrb); }
      
int CheMPS2::SyBookkeeper::gTwoS() const{ return Prob->gTwoS(); }

int CheMPS2::SyBookkeeper::gN() const{ return Prob->gN(); }

int CheMPS2::SyBookkeeper::gIrrep() const{ return Prob->gIrrep(); }

int CheMPS2::SyBookkeeper::gNmin(const int bound) const{ return Nmin[bound]; }

int CheMPS2::SyBookkeeper::gNmax(const int bound) const{ return Nmax[bound]; }

int CheMPS2::SyBookkeeper::gTwoSmin(const int bound, const int N) const{ return TwoSmin[bound][N-Nmin[bound]]; }

int CheMPS2::SyBookkeeper::gTwoSmax(const int bound, const int N) const{ return TwoSmax[bound][N-Nmin[bound]]; }

int CheMPS2::SyBookkeeper::gFCIdim(const int bound, const int N, const int TwoS, const int Icnt) const{ return gDimPrivate(FCIdim, bound, N, TwoS, Icnt); }

int CheMPS2::SyBookkeeper::gCurrentDim(const int bound, const int N, const int TwoS, const int Icnt) const{ return gDimPrivate(CurrentDim, bound, N, TwoS, Icnt); }

void CheMPS2::SyBookkeeper::SetDim(const int bound, const int N, const int TwoS, const int Icnt, const int val){

   if (gFCIdim(bound,N,TwoS,Icnt)!=0){
      CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = val;
   }
   
}

void CheMPS2::SyBookkeeper::fillFCIdim(){

   //First fill FCIdim from left
   for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++) FCIdim[0][0][0][Icnt] = 0;
   FCIdim[0][0][0][0] = 1;
   
   for (int bound=1; bound<=gL(); bound++){
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
               FCIdim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = min( CheMPS2::SYBK_dimensionCutoff, gFCIdim(bound-1,N,TwoS,Icnt) + gFCIdim(bound-1,N-2,TwoS,Icnt) + gFCIdim(bound-1,N-1,TwoS+1,directProd(Icnt,gIrrep(bound-1))) + gFCIdim(bound-1,N-1,TwoS-1,directProd(Icnt,gIrrep(bound-1))) );
            }
         }
      }
   }

   //Then allocate FCIdimRight
   int **** FCIdimRight = new int***[gL()+1];
   for (int bound=0; bound<=gL(); bound++){
      FCIdimRight[bound] = new int**[gNmax(bound)-gNmin(bound)+1];
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         FCIdimRight[bound][N-gNmin(bound)] = new int*[(gTwoSmax(bound,N)-gTwoSmin(bound,N))/2+1];
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            FCIdimRight[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2] = new int[getNumberOfIrreps()];
         }
      }
   }
   
   //Then calculate FCI dim from right -->FCIdimRight
   for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++) FCIdimRight[gL()][0][0][Icnt] = 0;
   FCIdimRight[gL()][0][0][gIrrep()] = 1;
   
   for (int bound=gL()-1; bound>=0; bound--){
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
               FCIdimRight[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = min( CheMPS2::SYBK_dimensionCutoff, gDimPrivate(FCIdimRight,bound+1,N,TwoS,Icnt) + gDimPrivate(FCIdimRight,bound+1,N+2,TwoS,Icnt) + gDimPrivate(FCIdimRight,bound+1,N+1,TwoS+1,directProd(Icnt,gIrrep(bound))) + gDimPrivate(FCIdimRight,bound+1,N+1,TwoS-1,directProd(Icnt,gIrrep(bound))) );
            }
         }
      }
   }
   
   //Then take min from FCIdim and FCIdimRight and store in FCIdim; and copy it to CurrentDim
   for (int bound=0; bound<=gL(); bound++){
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
               FCIdim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = min( gFCIdim(bound,N,TwoS,Icnt), gDimPrivate(FCIdimRight,bound,N,TwoS,Icnt) );
               CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = gFCIdim(bound,N,TwoS,Icnt);
            }
         }
      }
   }
   
   //Deallocate FCIdimRight
   for (int bound=0; bound<=gL(); bound++){
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            delete [] FCIdimRight[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2];
         }
         delete [] FCIdimRight[bound][N-gNmin(bound)];
      }
      delete [] FCIdimRight[bound];
   }
   delete [] FCIdimRight;

}

void CheMPS2::SyBookkeeper::ScaleCurrentDim(const int virtualD){

   for (int bound=1; bound<gL(); bound++){
      int totaldim = 0;
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
               totaldim += gCurrentDim(bound,N,TwoS,Icnt);
            }
         }
      }
      
      if (totaldim > virtualD){
         double factor = (1.0 * virtualD) / totaldim;
         for (int N=gNmin(bound); N<=gNmax(bound); N++){
            for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
               for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
                  CurrentDim[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt] = (ceil( factor * gCurrentDim(bound,N,TwoS,Icnt)) + 0.1);
               }
            }
         }
      }
      
      if (CheMPS2::SYBK_debugPrint){
         cout << "Bound = " << bound << endl;
         cout << "   Totaldim (FCI)        = " << totaldim << endl;
         totaldim = 0;
         for (int N=gNmin(bound); N<=gNmax(bound); N++){
            for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
               for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
                  totaldim += gCurrentDim(bound,N,TwoS,Icnt);
               }
            }
         }
         cout << "   Totaldim (rescaled)   = " << totaldim << endl;
      }
   }
   
   if (CheMPS2::SYBK_debugPrint){
      for (int bound=0; bound<=gL(); bound++){
         for (int N=gNmin(bound); N<=gNmax(bound); N++){
            for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
               for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
                  if ((gFCIdim(bound,N,TwoS,Icnt)!=0) && (gCurrentDim(bound,N,TwoS,Icnt)==0)){
                     cout << "Error at print: " << endl;
                     cout << "       gFCIdim(" << bound << "," << N << "," << TwoS << "," << Icnt << ") = " << gFCIdim(bound,N,TwoS,Icnt) << endl;
                     cout << "   gCurrentDim(" << bound << "," << N << "," << TwoS << "," << Icnt << ") = " << gCurrentDim(bound,N,TwoS,Icnt) << endl;
                  }
               }
            }
         }
      }
   }

}

int CheMPS2::SyBookkeeper::gDimPrivate(int **** storage, const int bound, const int N, const int TwoS, const int Icnt) const{

   if ((bound<0) || (bound>gL())) return 0;
   if ((N>gNmax(bound)) || (N<gNmin(bound))) return 0;
   if ((TwoS%2) != (gTwoSmin(bound,N)%2)) return 0;
   if ((TwoS < gTwoSmin(bound,N)) || (TwoS > gTwoSmax(bound,N))) return 0;
   if ((Icnt<0) || (Icnt>getNumberOfIrreps())) return 0;
   return storage[bound][N-gNmin(bound)][(TwoS-gTwoSmin(bound,N))/2][Icnt];

}

int CheMPS2::SyBookkeeper::gMaxDimAtBound(const int iBound) const{

   int maxDim = 0;
   for (int N=gNmin(iBound); N<=gNmax(iBound); N++){
      for (int TwoS=gTwoSmin(iBound,N); TwoS<=gTwoSmax(iBound,N); TwoS+=2){
         for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
            int dim = gCurrentDim(iBound,N,TwoS,Icnt);
            if (dim>maxDim) maxDim = dim;
         }
      }
   }
   return maxDim;

}

void CheMPS2::SyBookkeeper::print() const{

   for (int bound=0; bound<=gL(); bound++){
      int totaldim = 0;
      int totaldim2 = 0;
      cout << "   Nmin[" << bound << "] = " << gNmin(bound) << " and Nmax[" << bound << "] = " << gNmax(bound) << endl;
      for (int N=gNmin(bound); N<=gNmax(bound); N++){
         cout << "      2Smin[" << bound << "][" << N << "] = " << gTwoSmin(bound,N) << " and 2Smax[" << bound << "][" << N << "] = " << gTwoSmax(bound,N) << endl;
         for (int TwoS=gTwoSmin(bound,N); TwoS<=gTwoSmax(bound,N); TwoS+=2){
            cout << "         bound = " << bound << " and N = " << N << " and TwoS = " << TwoS << endl;
            for (int Icnt=0; Icnt<getNumberOfIrreps(); Icnt++){
               totaldim += gCurrentDim(bound,N,TwoS,Icnt);
               totaldim2 += (TwoS + 1) * gCurrentDim(bound,N,TwoS,Icnt);
               if ( gFCIdim(bound,N,TwoS,Icnt) != 0 ){
                  cout << "            gFCIdim(" << bound << "," << N << "," << TwoS << "," << Icnt << ") = " << gFCIdim(bound,N,TwoS,Icnt) << endl;
                  cout << "            gCurrentDim(" << bound << "," << N << "," << TwoS << "," << Icnt << ") = " << gCurrentDim(bound,N,TwoS,Icnt) << endl;
               }
            }
         }
      }
      cout << "CurrentDim (reduced  multiplets) at bound = " << bound << " is " << totaldim << endl;
      cout << "CurrentDim (complete multiplets) at bound = " << bound << " is " << totaldim2 << endl;
   }

}

bool CheMPS2::SyBookkeeper::IsPossible() const{

   return (gCurrentDim(gL(), gN(), gTwoS(), gIrrep())==1);

}



