/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013, 2014 Sebastian Wouters

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
#include <string>
#include <fstream>

#include "Irreps.h"
#include "TwoIndex.h"
#include "FourIndex.h"
#include "Hamiltonian.h"
#include "MyHDF5.h"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;

CheMPS2::Hamiltonian::Hamiltonian(const int Norbitals, const int nGroup, const int * OrbIrreps){

   L = Norbitals;
   assert( nGroup>=0 );
   assert( nGroup<=7 );
   SymmInfo.setGroup(nGroup);
   
   orb2irrep = new int[L];
   orb2indexSy = new int[L];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[nIrreps];
   for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
   for (int cnt=0; cnt<L; cnt++){
      assert( OrbIrreps[cnt]>=0 );
      assert( OrbIrreps[cnt]<nIrreps );
      orb2irrep[cnt] = OrbIrreps[cnt];
      orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
      irrep2num_orb[orb2irrep[cnt]]++;
   }
   
   Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);

}

CheMPS2::Hamiltonian::Hamiltonian(const string file_psi4text){

   CreateAndFillFromPsi4dump( file_psi4text );
   
}

CheMPS2::Hamiltonian::Hamiltonian(const bool fileh5, const string main_file, const string file_tmat, const string file_vmat){

   if (fileh5){
      CreateAndFillFromH5( main_file, file_tmat, file_vmat );
   } else {
      CreateAndFillFromPsi4dump( main_file );
   }
   
}

CheMPS2::Hamiltonian::~Hamiltonian(){
   
   delete [] orb2irrep;
   delete [] orb2indexSy;
   delete [] irrep2num_orb;
   delete Tmat;
   delete Vmat;
   
}

int CheMPS2::Hamiltonian::getL() const{ return L; }

int CheMPS2::Hamiltonian::getNGroup() const{ return SymmInfo.getGroupNumber(); }

int CheMPS2::Hamiltonian::getOrbitalIrrep(const int nOrb) const{ return orb2irrep[nOrb]; }

void CheMPS2::Hamiltonian::setEconst(const double val){ Econst = val; }

double CheMPS2::Hamiltonian::getEconst() const{ return Econst; }

void CheMPS2::Hamiltonian::setTmat(const int index1, const int index2, const double val){

   assert( orb2irrep[index1]==orb2irrep[index2] );
   Tmat->set(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2], val);

}

double CheMPS2::Hamiltonian::getTmat(const int index1, const int index2) const{

   if (orb2irrep[index1]==orb2irrep[index2]){
      return Tmat->get(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2]);
   }

   return 0.0;
   
}

void CheMPS2::Hamiltonian::setVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   assert( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) );
   Vmat->set(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);

}

void CheMPS2::Hamiltonian::addToVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   assert( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) );
   Vmat->add(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);

}

double CheMPS2::Hamiltonian::getVmat(const int index1, const int index2, const int index3, const int index4) const{

   if ( Irreps::directProd(orb2irrep[index1],orb2irrep[index2]) == Irreps::directProd(orb2irrep[index3],orb2irrep[index4]) ){
      return Vmat->get(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4]);
   }

   return 0.0;
   
}

void CheMPS2::Hamiltonian::save(const string file_parent, const string file_tmat, const string file_vmat) const{

   Tmat->save(file_tmat);
   Vmat->save(file_vmat);

   //The hdf5 file
   hid_t file_id = H5Fcreate(file_parent.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
         //The chain length
         hsize_t dimarray       = 1;
         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "L", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);
         
         //The group number
         hsize_t dimarray2      = 1;
         hid_t dataspace_id2    = H5Screate_simple(1, &dimarray2, NULL);
         hid_t dataset_id2      = H5Dcreate(group_id, "nGroup", H5T_STD_I32LE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         int nGroup = SymmInfo.getGroupNumber();
         H5Dwrite(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
         
         //orb2irrep
         hsize_t dimarray3      = L;
         hid_t dataspace_id3    = H5Screate_simple(1, &dimarray3, NULL);
         hid_t dataset_id3      = H5Dcreate(group_id, "orb2irrep", H5T_STD_I32LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);
         
         //Econst
         hsize_t dimarray4      = 1;
         hid_t dataspace_id4    = H5Screate_simple(1, &dimarray4, NULL);
         hid_t dataset_id4      = H5Dcreate(group_id, "Econst", H5T_IEEE_F64LE, dataspace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);
    
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);
         H5Dclose(dataset_id2);
         H5Sclose(dataspace_id2);
         H5Dclose(dataset_id3);
         H5Sclose(dataspace_id3);
         H5Dclose(dataset_id4);
         H5Sclose(dataspace_id4);

      H5Gclose(group_id);
      
   H5Fclose(file_id);

}

void CheMPS2::Hamiltonian::read(const string file_parent, const string file_tmat, const string file_vmat){

   Tmat->read(file_tmat);
   Vmat->read(file_vmat);

   //The hdf5 file
   hid_t file_id = H5Fopen(file_parent.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
         //The chain length
         hid_t dataset_id1 = H5Dopen(group_id, "L", H5P_DEFAULT);
         int Lagain;
         H5Dread(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lagain);
         assert( L==Lagain );
         
         //The group number
         hid_t dataset_id2 = H5Dopen(group_id, "nGroup", H5P_DEFAULT);
         int nGroup;
         H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
         assert( nGroup==SymmInfo.getGroupNumber() );
         
         //orb2irrep
         hid_t dataset_id3 = H5Dopen(group_id, "orb2irrep", H5P_DEFAULT);
         int * orb2irrepAgain = new int[L];
         H5Dread(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrepAgain);
         for (int cnt=0; cnt<L; cnt++){
            assert( orb2irrep[cnt]==orb2irrepAgain[cnt] );
         }
         delete [] orb2irrepAgain;
         
         //Econst
         hid_t dataset_id4 = H5Dopen(group_id, "Econst", H5P_DEFAULT);
         H5Dread(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);
    
         H5Dclose(dataset_id1);
         H5Dclose(dataset_id2);
         H5Dclose(dataset_id3);
         H5Dclose(dataset_id4);

      H5Gclose(group_id);
      
   H5Fclose(file_id);
   
   if (CheMPS2::HAMILTONIAN_debugPrint) debugcheck();

}

void CheMPS2::Hamiltonian::CreateAndFillFromH5(const string file_parent, const string file_tmat, const string file_vmat){

   //The hdf5 file
   hid_t file_id = H5Fopen(file_parent.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
         //The number of orbitals: Set L directly!
         hid_t dataset_id1 = H5Dopen(group_id, "L", H5P_DEFAULT);
         H5Dread(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);
         
         //The group number: SymmInfo.setGroup()
         hid_t dataset_id2 = H5Dopen(group_id, "nGroup", H5P_DEFAULT);
         int nGroup_LOADH5;
         H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup_LOADH5);
         SymmInfo.setGroup( nGroup_LOADH5 );
         
         //The irrep of each orbital: Create and fill orb2irrep directly!
         hid_t dataset_id3 = H5Dopen(group_id, "orb2irrep", H5P_DEFAULT);
         orb2irrep = new int[ L ];
         H5Dread(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);
    
         H5Dclose(dataset_id1);
         H5Dclose(dataset_id2);
         H5Dclose(dataset_id3);

      H5Gclose(group_id);
      
   H5Fclose(file_id);
   
   orb2indexSy = new int[ L ];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[ nIrreps ];
   for (int irrep = 0; irrep < nIrreps; irrep++){ irrep2num_orb[ irrep ] = 0; }
   for (int orb = 0; orb < L; orb++){
      orb2indexSy[ orb ] = irrep2num_orb[ orb2irrep[ orb ] ];
      irrep2num_orb[ orb2irrep[ orb ] ]++;
   }
   
   Tmat = new TwoIndex(  SymmInfo.getGroupNumber(), irrep2num_orb );
   Vmat = new FourIndex( SymmInfo.getGroupNumber(), irrep2num_orb );

   read(file_parent, file_tmat, file_vmat);

}

//Works for the file mointegrals/mointegrals.cc_PRINT which can be used as a plugin in psi4 beta5
void CheMPS2::Hamiltonian::CreateAndFillFromPsi4dump(const string filename){

   string line, part;
   int pos;
   
   ifstream inputfile(filename.c_str());
   
   //First go to the start of the integral dump.
   bool stop = false;
   string start = "****  Molecular Integrals For CheMPS Start Here";
   do{
      getline(inputfile,line);
      pos = line.find(start);
      if (pos==0) stop = true;
   } while (!stop);
   
   //Get the group name and convert it to the group number
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   int nGroup = 0;
   stop = false;
   do {
      if (part.compare(SymmInfo.getGroupName(nGroup))==0) stop = true;
      else nGroup += 1;
   } while (!stop);
   SymmInfo.setGroup(nGroup);
   //cout << "The group was found to be " << SymmInfo.getGroupName() << " ." << endl;
   
   //This line says how many irreps there are: skip.
   getline(inputfile,line);
   
   //This line contains the nuclear energy part.
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   Econst = atof(part.c_str());

   //This line contains the number of MO's.
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   L = atoi(part.c_str());
   
   //This line contains only text
   getline(inputfile,line);
   
   //This line contains the irrep numbers --> allocate, read in & set
   getline(inputfile,line);
   
   orb2irrep = new int[L];
   orb2indexSy = new int[L];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[nIrreps];
   
   pos = 0;
   do {
      orb2irrep[pos] = atoi(line.substr(2*pos,1).c_str());
      pos++;
   } while (2*pos < (int)line.size()-1);
   
   for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
   for (int cnt=0; cnt<L; cnt++){
      orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
      irrep2num_orb[orb2irrep[cnt]]++;
   }
   Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   
   //Skip three lines --> number of double occupations, single occupations and test line
   getline(inputfile,line);
   getline(inputfile,line);
   getline(inputfile,line);
   
   //Read in one-electron integrals
   getline(inputfile,line);
   int pos2, index1, index2;
   double value;
   while( (line.substr(0,1)).compare("*")!=0 ){
   
      pos = 0;
      pos2 = line.find(" ",pos);
      index1 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index2 = atoi(line.substr(pos,pos2-pos).c_str());
      
      value = atof(line.substr(pos2+1,line.size()-pos2-2).c_str());
      
      setTmat(index1,index2,value);
      
      getline(inputfile,line);
   
   }
   
   //Read in two-electron integrals --> in file: chemical notation; in Vmat: physics notation
   getline(inputfile,line);
   int index3, index4;
   while( (line.substr(0,1)).compare("*")!=0 ){
   
      pos = 0;
      pos2 = line.find(" ",pos);
      index1 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index2 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index3 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index4 = atoi(line.substr(pos,pos2-pos).c_str());
      
      value = atof(line.substr(pos2+1,line.size()-pos2-2).c_str());
      
      setVmat(index1, index3, index2, index4, value);
      
      getline(inputfile,line);
   
   }
   
   if (CheMPS2::HAMILTONIAN_debugPrint) debugcheck();
   
   inputfile.close();

}

void CheMPS2::Hamiltonian::debugcheck() const{

   cout << "Econst = " << Econst << endl;
   
   double test = 0.0;
   double test2 = 0.0;
   double test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getTmat(i,i);
      for (int j=0; j<L; j++){
         test += getTmat(i,j);
         if (i<=j) test2 += getTmat(i,j);
      }
   }
   cout << "1-electron integrals: Trace                  : " << test3 << endl;
   cout << "1-electron integrals: Sum over all elements  : " << test << endl;
   cout << "1-electron integrals: Sum over Tij with i<=j : " << test2 << endl;
      
   test = 0.0;
   test2 = 0.0;
   test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getVmat(i,i,i,i);
      for (int j=0; j<L; j++){
         for (int k=0; k<L; k++){
            for (int l=0; l<L; l++){
               test += getVmat(i,j,k,l);
               if ((i<=j) && (j<=k) && (k<=l)) test2 += getVmat(i,j,k,l);
            }
         }
      }
   }
   cout << "2-electron integrals: Trace                          : " << test3 << endl;
   cout << "2-electron integrals: Sum over all elements          : " << test << endl;
   cout << "2-electron integrals: Sum over Vijkl with i<=j<=k<=l : " << test2 << endl;

}


