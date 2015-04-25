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

#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "Problem.h"

using std::cout;
using std::endl;

CheMPS2::Problem::Problem(const Hamiltonian * Hamin, const int TwoSin, const int Nin, const int Irrepin){

   Ham = Hamin;
   TwoS = TwoSin;
   N = Nin;
   OneOverNMinusOne = 1.0/(N-1);
   Irrep = Irrepin;
   bReorderD2h = false;
   
   checkConsistency();

}

CheMPS2::Problem::~Problem(){

   if (bReorderD2h){
      delete [] f1;
      delete [] f2;
   }

}

void CheMPS2::Problem::SetupReorderD2h(){

   if (gSy()==7){ //Only if D2h of course
   
      bReorderD2h = true;
      f1 = new int[Ham->getL()];
      f2 = new int[Ham->getL()];
      
      int DMRGirrepOrder[8];
      DMRGirrepOrder[0] = 0; //Ag  sigma
      DMRGirrepOrder[1] = 5; //B1u sigma^*
      DMRGirrepOrder[2] = 7; //B3u pi_x
      DMRGirrepOrder[3] = 2; //B2g pi_x^*
      DMRGirrepOrder[4] = 6; //B2u pi_y
      DMRGirrepOrder[5] = 3; //B3g pi_y^*
      DMRGirrepOrder[6] = 1; //B1g
      DMRGirrepOrder[7] = 4; //Au
      
      int DMRGOrb = 0;
      for (int irrep=0; irrep<8; irrep++){
         for (int HamOrb=0; HamOrb<Ham->getL(); HamOrb++){
            if (Ham->getOrbitalIrrep(HamOrb)==DMRGirrepOrder[irrep]){
               f1[HamOrb] = DMRGOrb;
               f2[DMRGOrb] = HamOrb;
               DMRGOrb++;
            }
         }
      }
      assert( DMRGOrb==Ham->getL() );
      
   }

}

int CheMPS2::Problem::gL() const{ return Ham->getL(); }
int CheMPS2::Problem::gSy() const{ return Ham->getNGroup(); }

int CheMPS2::Problem::gIrrep(const int nOrb) const{
   
   if (!bReorderD2h){
      return Ham->getOrbitalIrrep(nOrb);
   }
   
   return Ham->getOrbitalIrrep(f2[nOrb]);

}

int CheMPS2::Problem::gTwoS() const{ return TwoS; }
int CheMPS2::Problem::gN() const{ return N; }
int CheMPS2::Problem::gIrrep() const{ return Irrep; }
double CheMPS2::Problem::gEconst() const{ return Ham->getEconst(); }

bool CheMPS2::Problem::gReorderD2h() const{ return bReorderD2h; }
int CheMPS2::Problem::gf1(const int HamOrb) const{ return (bReorderD2h)?f1[HamOrb]:-1; }
int CheMPS2::Problem::gf2(const int DMRGOrb) const{ return (bReorderD2h)?f2[DMRGOrb]:-1; }

double CheMPS2::Problem::gMxElement(const int alpha, const int beta, const int gamma, const int delta) const{

   if (!bReorderD2h){
      return Ham->getVmat(alpha, beta, gamma, delta) + OneOverNMinusOne*(((alpha==gamma)?Ham->getTmat(beta,delta):0) + ((beta==delta)?Ham->getTmat(alpha,gamma):0));
   }
   
   return Ham->getVmat(f2[alpha],f2[beta],f2[gamma],f2[delta]) + OneOverNMinusOne*(((alpha==gamma)?Ham->getTmat(f2[beta],f2[delta]):0) + ((beta==delta)?Ham->getTmat(f2[alpha],f2[gamma]):0));

}

bool CheMPS2::Problem::checkConsistency() const{

   Irreps SymmInfo(gSy());
   if ((gIrrep()<0) || (gIrrep()>=SymmInfo.getNumberOfIrreps())){
      cout << "Problem::Problem() : Irrep out of bound : Irrep = " << gIrrep() << endl;
      return false;
   }
   if (gTwoS()<0){
      cout << "Problem::checkConsistency() : TwoS = " << gTwoS() << endl;
      return false;
   }
   if (gN()<0){
      cout << "Problem::checkConsistency() : N = " << gN() << endl;
      return false;
   }
   if (gL()<0){
      cout << "Problem::checkConsistency() : L = " << gL() << endl;
      return false;
   }
   if (gN()>2*gL()){
      cout << "Problem::checkConsistency() : N > 2*L ; N = " << gN() << " and L = " << gL() << endl;
      return false;
   }
   if ( (gN()%2) != (gTwoS()%2) ){
      cout << "Problem::checkConsistency() : N%2 != TwoS%2 ; N = " << gN() << " and TwoS = " << gTwoS() << endl;
      return false;
   }
   if ( gTwoS() > gL() - abs(gN() - gL()) ){
      cout << "Problem::checkConsistency() : TwoS > L - |N-L| ; N = " << gN() << " and TwoS = " << gTwoS() << " and L = " << gL() << endl;
      return false;
   }
   
   return true;

}


