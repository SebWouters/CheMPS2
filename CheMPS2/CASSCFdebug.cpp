/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013 Sebastian Wouters

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
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <hdf5.h>

#include "CASSCF.h"
#include "Lapack.h"

using std::string;
using std::ifstream;
using std::cout;
using std::cerr;
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
         
         if (cnt2-passed2 == OrbPerIrrep[irrep2]-1){
            passed2 += OrbPerIrrep[irrep2];
            irrep2++;
         }
            
      }
      
      cout << "   Orb " << cnt << " : "<< cnt-passed+1 << SymmInfo.getIrrepName(irrep) << " = " << SPenergy << endl;
         
      if ((cnt-passed == OrbPerIrrep[irrep]-1) && (irrep<numberOfIrreps-1)){
         passed += OrbPerIrrep[irrep];
         irrep++;
      }
         
   }
 
   cout << "HF energy = " << EnergyHF << endl;

}

void CheMPS2::CASSCF::check1DMand2DMrotated(const int N){

   //Check particle number by tracing the 1DMrotated --> IS OK
   double SingleTrace1DM = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      SingleTrace1DM += get1DMrotated(cnt1,cnt1);
   }
   cout << "CASSCF :: Single trace 1DMrotated = " << SingleTrace1DM << " and N = " << N << endl;
   
   //Compare the 1DM eigenvalues with the NOON from DMRG --> IS OK
   double * OneDM = new double[L*L];
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         OneDM[cnt1+L*cnt2] = get1DMrotated(cnt1,cnt2);
      }
   }
   double * eigenval = new double[L];
   int lwork = 3*L;
   double * work = new double[lwork];
   char jobz = 'N';
   char uplo = 'U';
   int info;
   dsyev_(&jobz,&uplo,&L,OneDM,&L,eigenval,work,&lwork,&info);
   cout << "CASSCF :: CASSCF 1DM eigenvalues= [ ";
   for (int cnt=0; cnt<L-1; cnt++){ cout << eigenval[cnt] << " , "; }
   cout << eigenval[L-1] << " ]." << endl;
   delete [] OneDM;
   delete [] eigenval;
   delete [] work;

   //Compare the double trace of 2DMrotated with N*(N-1) --> IS OK
   double DoubleTrace = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         DoubleTrace += get2DMrotated(cnt1,cnt2,cnt1,cnt2);
      }
   }
   cout << "CASSCF :: Double trace 2DMrotated = " << DoubleTrace << " and N*(N-1) = " << N*(N-1) << endl;
   
   //Check the 2-norm of 1DMrotated and 1/(N-1)*trace(2DMrotated) --> IS OK
   double prefactor = 1.0/(N-1.0);
   double * SingleTrace = new double[L*L];
   double OneDMerror = 0.0;
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         SingleTrace[cnt1 + L*cnt2] = 0.0;
         for (int cnt3=0; cnt3<L; cnt3++){ SingleTrace[cnt1 + L*cnt2] += get2DMrotated(cnt1,cnt3,cnt2,cnt3); }
         SingleTrace[cnt1 + L*cnt2] *= prefactor;
         OneDMerror += ( SingleTrace[cnt1+L*cnt2] - get1DMrotated(cnt1,cnt2) ) * ( SingleTrace[cnt1+L*cnt2] - get1DMrotated(cnt1,cnt2) );
      }
   }
   OneDMerror = sqrt(OneDMerror);
   cout << "CASSCF :: 2-norm of (1DMrotated - (Single trace 2DMrotated)/(N-1)) = " << OneDMerror << endl;
   delete [] SingleTrace;
   
   //Check the energy based on the 1DMrotated and 2DMrotated --> IS OK
   double Energy = HamRotated->getEconst();
   for (int cnt1=0; cnt1<L; cnt1++){
      for (int cnt2=0; cnt2<L; cnt2++){
         Energy += HamRotated->getTmat(cnt1,cnt2) * get1DMrotated(cnt1,cnt2);
         for (int cnt3=0; cnt3<L; cnt3++){
            for (int cnt4=0; cnt4<L; cnt4++){
               Energy += 0.5 * HamRotated->getVmat(cnt1,cnt2,cnt3,cnt4) * get2DMrotated(cnt1,cnt2,cnt3,cnt4);
            }
         }
      }
   }
   cout << "CASSCF :: Energy (with 1DM and 2DM) = " << Energy << endl;

}

void CheMPS2::CASSCF::saveU(){

   hid_t file_id = H5Fcreate(CheMPS2::CASSCF_unitaryStorageName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
   
      std::stringstream irrepname;
      irrepname << "irrep_" << cnt;
      
      hsize_t dimarray      = OrbPerIrrep[cnt] * OrbPerIrrep[cnt];
      hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
      hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[cnt]);

      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::CASSCF::loadU(){

   hid_t file_id = H5Fopen(CheMPS2::CASSCF_unitaryStorageName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
   for (int cnt=0; cnt<numberOfIrreps; cnt++){
   
      std::stringstream irrepname;
      irrepname << "irrep_" << cnt;

      hid_t dataset_id = H5Dopen(group_id, irrepname.str().c_str(), H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[cnt]);
         
      H5Dclose(dataset_id);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::CASSCF::deleteStoredUnitary(){

   std::stringstream temp;
   temp << "rm " << CheMPS2::CASSCF_unitaryStorageName;
   int info = system(temp.str().c_str());
   cout << "Info on CASSCF::Unitary rm call to system: " << info << endl;

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
         for (int orb=Nocc[irrep]; orb<3; orb++){ //1-3 Ag : double
            coeff1[passed+orb] = coeff2[passed+orb] = 2;
         }
      }
      if (irrep==5){ //B1u
         for (int orb=Nocc[irrep]; orb<2; orb++){ //1-2 B1u : double
            coeff1[passed+orb] = coeff2[passed+orb] = 2;
         }
      }
      if (irrep==7){ //B3u
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B3u orbitals not allowed." << endl; }
         coeff1[passed+0] = 2; //1 B3u : double for | 1pi_x^2 >
      }
      if (irrep==6){ //B2u
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B2u orbitals not allowed." << endl; }
         coeff2[passed+0] = 2; //1 B2u : double for | 1pi_y^2 >
      }
      passed += NDMRG[irrep];
   }
   
   int * coeff3 = new int[nOrbDMRG]; //1pi_x^1 1pi_x^{*1}
   int * coeff4 = new int[nOrbDMRG]; //1pi_y^1 1pi_y^{*1}
   //Fill coeff3 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B3u^1 3Ag^2 1B2g^1 > or | 1pi_x^1 1pi_x^{*1} >
   //Fill coeff4 with | 1Ag^2 1B1u^2 2Ag^2 2B1u^2 1B2u^1 3Ag^2 1B3g^1 > or | 1pi_y^1 1pi_y^{*1} >
   
   for (int cnt=0; cnt<nOrbDMRG; cnt++){ coeff3[cnt] = coeff4[cnt] = 0; } //No particles at all
   passed = 0;
   for (int irrep=0; irrep<numberOfIrreps; irrep++){
      if (irrep==0){ //Ag
         for (int orb=Nocc[irrep]; orb<3; orb++){ //1-3 Ag : double
            coeff3[passed+orb] = coeff4[passed+orb] = 2;
         }
      }
      if (irrep==5){ //B1u
         for (int orb=Nocc[irrep]; orb<2; orb++){ //1-2 B1u : double
            coeff3[passed+orb] = coeff4[passed+orb] = 2;
         }
      }
      if (irrep==7){ //B3u
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B3u orbitals not allowed." << endl; }
         coeff3[passed+0] = 1; //1 B3u : single for | 1pi_x^1 1pi_x^{*1} >
      }
      if (irrep==6){ //B2u
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B2u orbitals not allowed." << endl; }
         coeff4[passed+0] = 1; //1 B2u : single for | 1pi_y^1 1pi_y^{*1} >
      }
      if (irrep==2){ //B2g
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B2g orbitals not allowed." << endl; }
         coeff3[passed+0] = 1; //1 B2g : single for | 1pi_x^1 1pi_x^{*1} >
      }
      if (irrep==3){ //B3g
         if (Nocc[irrep] > 0){ cerr << "Coeff printer: condensed B3g orbitals not allowed." << endl; }
         coeff4[passed+0] = 1; //1 B3g : single for | 1pi_y^1 1pi_y^{*1} >
      }
      passed += NDMRG[irrep];
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


