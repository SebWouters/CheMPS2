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
#include <string>

#include "Irreps.h"

using std::string;
using std::cout;
using std::endl;

CheMPS2::Irreps::Irreps(){

   isActivated = false;

}

CheMPS2::Irreps::Irreps(const int nGroup){

   if ((nGroup >= 0) && (nGroup <= 7)){
      isActivated = true;
      groupNumber = nGroup;
      nIrreps = getNumberOfIrreps(groupNumber);
   } else {
      isActivated = false;
   }

}

CheMPS2::Irreps::~Irreps(){ }

bool CheMPS2::Irreps::setGroup(const int nGroup){

   if ((nGroup >= 0) && (nGroup <= 7)){
      isActivated = true;
      groupNumber = nGroup;
      nIrreps = getNumberOfIrreps(groupNumber);
   } else {
      isActivated = false;
   }
   
   return isActivated;

}

bool CheMPS2::Irreps::getIsActivated() const{ return isActivated; }

int CheMPS2::Irreps::getGroupNumber() const{ return isActivated ? groupNumber : -1 ; }

string CheMPS2::Irreps::getGroupName() const{

   return isActivated ? getGroupNamePrivate(groupNumber) : "error" ;

}

string CheMPS2::Irreps::getGroupName(const int nGroup){

   return ((nGroup>=0)&&(nGroup<=7)) ? getGroupNamePrivate(nGroup) : "error" ;

}

string CheMPS2::Irreps::getGroupNamePrivate(const int nGroup){

   if (nGroup==0) return "c1";
   if (nGroup==1) return "ci";
   if (nGroup==2) return "c2";
   if (nGroup==3) return "cs";
   if (nGroup==4) return "d2";
   if (nGroup==5) return "c2v";
   if (nGroup==6) return "c2h";
   if (nGroup==7) return "d2h";
   return "error";

}

int CheMPS2::Irreps::getNumberOfIrreps() const{

   return isActivated ? nIrreps : -1 ;

}

int CheMPS2::Irreps::getNumberOfIrreps(const int nGroup){
   
   if ((nGroup < 0) || ( nGroup > 7)){ return -1; }
   if (nGroup == 0) return 1;
   if (nGroup <= 3) return 2;
   if (nGroup <= 6) return 4;
   return 8;

}

string CheMPS2::Irreps::getIrrepName(const int irrepNumber) const{

   if (!isActivated) return "error1";
   
   if ( (irrepNumber<0) || (irrepNumber >= nIrreps) ) return "error2";
   
   return getIrrepNamePrivate(groupNumber,irrepNumber);
   
}
   
string CheMPS2::Irreps::getIrrepNamePrivate(const int nGroup, const int nIrrep){
   
   if (nGroup == 0){
      if (nIrrep == 0) return "A";
   }
   
   if (nGroup == 1){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "Au";
   }
   
   if (nGroup == 2){
      if (nIrrep == 0) return "A";
      if (nIrrep == 1) return "B";
   }
   
   if (nGroup == 3){
      if (nIrrep == 0) return "Ap";
      if (nIrrep == 1) return "App";
   }
   
   if (nGroup == 4){
      if (nIrrep == 0) return "A";
      if (nIrrep == 1) return "B1";
      if (nIrrep == 2) return "B2";
      if (nIrrep == 3) return "B3";
   }
   
   if (nGroup == 5){
      if (nIrrep == 0) return "A1";
      if (nIrrep == 1) return "A2";
      if (nIrrep == 2) return "B1";
      if (nIrrep == 3) return "B2";
   }
   
   if (nGroup == 6){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "Bg";
      if (nIrrep == 2) return "Au";
      if (nIrrep == 3) return "Bu";
   }
   
   if (nGroup == 7){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "B1g";
      if (nIrrep == 2) return "B2g";
      if (nIrrep == 3) return "B3g";
      if (nIrrep == 4) return "Au";
      if (nIrrep == 5) return "B1u";
      if (nIrrep == 6) return "B2u";
      if (nIrrep == 7) return "B3u";
   }
   
   return "error2";

}

void CheMPS2::Irreps::symm_psi2molpro( int * psi2molpro ) const{

   if (!isActivated) return;
   symm_psi2molpro( psi2molpro, getGroupName() );

}

void CheMPS2::Irreps::symm_psi2molpro( int * psi2molpro, const string SymmLabel ){

   if ( SymmLabel.compare("c1")==0 ){
      psi2molpro[0] = 1;
   }
   if ( ( SymmLabel.compare("ci")==0 ) || ( SymmLabel.compare("c2")==0 ) || ( SymmLabel.compare("cs")==0 ) ){
      psi2molpro[0] = 1;
      psi2molpro[1] = 2;
   }
   if ( ( SymmLabel.compare("d2")==0 ) ){
      psi2molpro[0] = 1;
      psi2molpro[1] = 4;
      psi2molpro[2] = 3;
      psi2molpro[3] = 2;
   }
   if ( ( SymmLabel.compare("c2v")==0 ) || ( SymmLabel.compare("c2h")==0 ) ){
      psi2molpro[0] = 1;
      psi2molpro[1] = 4;
      psi2molpro[2] = 2;
      psi2molpro[3] = 3;
   }
   if ( ( SymmLabel.compare("d2h")==0 ) ){
      psi2molpro[0] = 1;
      psi2molpro[1] = 4;
      psi2molpro[2] = 6;
      psi2molpro[3] = 7;
      psi2molpro[4] = 8;
      psi2molpro[5] = 5;
      psi2molpro[6] = 3;
      psi2molpro[7] = 2;
   }

}

void CheMPS2::Irreps::printAll(){

   for (int thegroup=0; thegroup<8; thegroup++){
      cout << "######################################################" << endl;
      cout << "Name = " << getGroupNamePrivate(thegroup) << endl;
      cout << "nIrreps = " << getNumberOfIrreps(thegroup) << endl;
      cout << "Multiplication table :" << endl;
      for (int irrep1=-1; irrep1<getNumberOfIrreps(thegroup); irrep1++){
         for (int irrep2=-1; irrep2<getNumberOfIrreps(thegroup); irrep2++){
            if ((irrep1 == -1) && (irrep2 == -1)) cout << "\t";
            if ((irrep1 == -1) && (irrep2 >= 0 )) cout << getIrrepNamePrivate(thegroup, irrep2) << "\t";
            if ((irrep2 == -1) && (irrep1 >= 0 )) cout << getIrrepNamePrivate(thegroup, irrep1) << "\t";
            if ((irrep2 >=  0) && (irrep1 >= 0 )) cout << getIrrepNamePrivate(thegroup, directProd(irrep1, irrep2)) << "\t";
         }
         cout << endl;
      }
   }
   cout << "######################################################" << endl;
   
}


