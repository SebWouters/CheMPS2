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
#include <math.h>

#include "CASSCF.h"

using std::cout;
using std::endl;

void CheMPS2::CASSCF::checkHF(){
   
   double EnergyHF = HamOrig->getEconst();
      
   int passed = 0;
   int irrep = 0;
   cout << "Single particle energy levels : " << endl;
   for (int cnt=0; cnt<L; cnt++){
   
      double SPenergy = HamOrig->getTmat(cnt,cnt);
   
      if (cnt-passed < DOCC[irrep]){ EnergyHF += 2*HamOrig->getTmat(cnt,cnt); }
      else{
         if (cnt-passed < SOCC[irrep] + DOCC[irrep]) EnergyHF +=   HamOrig->getTmat(cnt,cnt);
      }
      
      int passed2 = 0;
      int irrep2 = 0;
      for (int cnt2=0; cnt2<L; cnt2++){
      
         if (cnt2-passed2 < DOCC[irrep2]){ SPenergy += 2*HamOrig->getVmat(cnt,cnt2,cnt,cnt2) - HamOrig->getVmat(cnt,cnt,cnt2,cnt2); }
         else {
            if (cnt2-passed2 < SOCC[irrep2] + DOCC[irrep2]) SPenergy += HamOrig->getVmat(cnt,cnt2,cnt,cnt2) - HamOrig->getVmat(cnt,cnt,cnt2,cnt2);
         }
         
         if ((cnt-passed < DOCC[irrep]) && (cnt2-passed2 < DOCC[irrep2])){
            EnergyHF +=   2*HamOrig->getVmat(cnt,cnt2,cnt,cnt2) -     HamOrig->getVmat(cnt,cnt,cnt2,cnt2);
         }
         if ((cnt-passed >= DOCC[irrep]) && (cnt-passed < SOCC[irrep] + DOCC[irrep]) && (cnt2-passed2 < DOCC[irrep2])){
            EnergyHF +=     HamOrig->getVmat(cnt,cnt2,cnt,cnt2) - 0.5*HamOrig->getVmat(cnt,cnt,cnt2,cnt2);
         }
         if ((cnt-passed < DOCC[irrep]) && (cnt2-passed2 >= DOCC[irrep2]) && (cnt2-passed2 < SOCC[irrep2] + DOCC[irrep2])){
            EnergyHF +=     HamOrig->getVmat(cnt,cnt2,cnt,cnt2) - 0.5*HamOrig->getVmat(cnt,cnt,cnt2,cnt2);
         }
         if ((cnt-passed >= DOCC[irrep]) && (cnt-passed < SOCC[irrep] + DOCC[irrep]) && (cnt2-passed2 >= DOCC[irrep2]) && (cnt2-passed2 < SOCC[irrep2] + DOCC[irrep2])){
            EnergyHF += 0.5*HamOrig->getVmat(cnt,cnt2,cnt,cnt2) - 0.5*HamOrig->getVmat(cnt,cnt,cnt2,cnt2);
         }
         
         if (cnt2-passed2 == iHandler->getNORB(irrep2)-1){
            passed2 += iHandler->getNORB(irrep2);
            irrep2++;
         }
            
      }
      
      cout << "   Orb " << cnt << " : "<< cnt-passed+1 << SymmInfo.getIrrepName(irrep) << " = " << SPenergy << endl;
         
      if ((cnt-passed == iHandler->getNORB(irrep)-1) && (irrep<numberOfIrreps-1)){
         passed += iHandler->getNORB(irrep);
         irrep++;
      }
         
   }
 
   cout << "HF energy = " << EnergyHF << endl;

}

void CheMPS2::CASSCF::PrintCoeff_C2(DMRG * theDMRG){

   //cc-pVDZ  full active space : Ag 7 / B1g 1 / B2g 3 / B3g 3 / Au 1 / B1u 7 / B2u 3 / B3u 3 ==> total 28 orbitals
   //cc-pCVDZ full active space : Ag 9 / B1g 1 / B2g 4 / B3g 4 / Au 1 / B1u 9 / B2u 4 / B3u 4 ==> total 36 orbitals

   int * coeff1 = new int[nOrbDMRG]; //1pi_x^2
   int * coeff2 = new int[nOrbDMRG]; //1pi_y^2
   //Fill coeff1 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B3u^2 3Ag^2 > or | 1pi_x^2 >
   //Fill coeff2 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B2u^2 3Ag^2 > or | 1pi_y^2 >
   
   for (int cnt=0; cnt<nOrbDMRG; cnt++){ coeff1[cnt] = coeff2[cnt] = 0; } //No particles at all
   int passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      if (irrep==0){ //Ag
         for (int orb=iHandler->getNOCC(irrep); orb<3; orb++){ //1-3 Ag : double
            coeff1[passed+orb] = coeff2[passed+orb] = 2;
         }
      }
      if (irrep==5){ //B1u
         for (int orb=iHandler->getNOCC(irrep); orb<2; orb++){ //1-2 B1u : double
            coeff1[passed+orb] = coeff2[passed+orb] = 2;
         }
      }
      if (irrep==7){ //B3u
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B3u orbitals not allowed.
         coeff1[passed+0] = 2; //1 B3u : double for | 1pi_x^2 >
      }
      if (irrep==6){ //B2u
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B2u orbitals not allowed.
         coeff2[passed+0] = 2; //1 B2u : double for | 1pi_y^2 >
      }
      passed += iHandler->getNDMRG(irrep);
   }
   
   int * coeff3 = new int[nOrbDMRG]; //1pi_x^1 1pi_x^{*1}
   int * coeff4 = new int[nOrbDMRG]; //1pi_y^1 1pi_y^{*1}
   //Fill coeff3 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B3u^1 3Ag^2 1B2g^1 > or | 1pi_x^1 1pi_x^{*1} >
   //Fill coeff4 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B2u^1 3Ag^2 1B3g^1 > or | 1pi_y^1 1pi_y^{*1} >
   
   for (int cnt=0; cnt<nOrbDMRG; cnt++){ coeff3[cnt] = coeff4[cnt] = 0; } //No particles at all
   passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      if (irrep==0){ //Ag
         for (int orb=iHandler->getNOCC(irrep); orb<3; orb++){ //1-3 Ag : double
            coeff3[passed+orb] = coeff4[passed+orb] = 2;
         }
      }
      if (irrep==5){ //B1u
         for (int orb=iHandler->getNOCC(irrep); orb<2; orb++){ //1-2 B1u : double
            coeff3[passed+orb] = coeff4[passed+orb] = 2;
         }
      }
      if (irrep==7){ //B3u
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B3u orbitals not allowed.
         coeff3[passed+0] = 1; //1 B3u : single for | 1pi_x^1 1pi_x^{*1} >
      }
      if (irrep==6){ //B2u
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B2u orbitals not allowed.
         coeff4[passed+0] = 1; //1 B2u : single for | 1pi_y^1 1pi_y^{*1} >
      }
      if (irrep==2){ //B2g
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B2g orbitals not allowed.
         coeff3[passed+0] = 1; //1 B2g : single for | 1pi_x^1 1pi_x^{*1} >
      }
      if (irrep==3){ //B3g
         assert( iHandler->getNOCC(irrep)==0 ); //Condensed B3g orbitals not allowed.
         coeff4[passed+0] = 1; //1 B3g : single for | 1pi_y^1 1pi_y^{*1} >
      }
      passed += iHandler->getNDMRG(irrep);
   }
      
   double thecoeff1 = theDMRG->getSpecificCoefficient(coeff1);
   cout << "Coeff of 1pi_x^2 = " << thecoeff1 << endl;
   double thecoeff2 = theDMRG->getSpecificCoefficient(coeff2);
   cout << "Coeff of 1pi_y^2 = " << thecoeff2 << endl;
   double thecoeff3 = theDMRG->getSpecificCoefficient(coeff3);
   cout << "Coeff of 1pi_x^1 1pi_x^{*1} = " << thecoeff3 << endl;
   double thecoeff4 = theDMRG->getSpecificCoefficient(coeff4);
   cout << "Coeff of 1pi_y^1 1pi_y^{*1} = " << thecoeff4 << endl;
         
   delete [] coeff1;
   delete [] coeff2;
   delete [] coeff3;
   delete [] coeff4;

}


