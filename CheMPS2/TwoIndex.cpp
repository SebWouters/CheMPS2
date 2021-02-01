/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2021 Sebastian Wouters

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
#include <sstream>
#include <string>

#include <hdf5.h>

#include "TwoIndex.h"

using namespace std;

CheMPS2::TwoIndex::TwoIndex(const int nGroup, const int * IrrepSizes){

   SymmInfo.setGroup(nGroup);
   
   Isizes = new int[SymmInfo.getNumberOfIrreps()];
   storage = new double*[SymmInfo.getNumberOfIrreps()];
   
   for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
      Isizes[cnt] = IrrepSizes[cnt];
      if (Isizes[cnt]>0) storage[cnt] = new double[Isizes[cnt]*(Isizes[cnt]+1)/2];
   }
   
   Clear();

}

CheMPS2::TwoIndex::~TwoIndex(){
   
   for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++) if (Isizes[cnt]>0) delete [] storage[cnt];
   delete [] storage;
   delete [] Isizes;
   
}

void CheMPS2::TwoIndex::Clear(){

   for (int irrep = 0; irrep < SymmInfo.getNumberOfIrreps(); irrep++){
      const int loopsize = (Isizes[irrep]*(Isizes[irrep]+1))/2;
      for (int count = 0; count < loopsize; count++){ storage[irrep][count] = 0.0; }
   }

}

void CheMPS2::TwoIndex::set(const int irrep, const int i, const int j, const double val){

   if (i>j) storage[irrep][j + i*(i+1)/2] = val;
   else storage[irrep][i + j*(j+1)/2] = val;

}

double CheMPS2::TwoIndex::get(const int irrep, const int i, const int j) const{

   if (i>j) return storage[irrep][j + i*(i+1)/2];
   return storage[irrep][i + j*(j+1)/2];

}

void CheMPS2::TwoIndex::save(const std::string name) const{
 
   //The hdf5 file
   hid_t file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      //The metadata
      hid_t group_id = H5Gcreate(file_id, "/MetaData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
         //The IrrepSizes
         hsize_t dimarray       = SymmInfo.getNumberOfIrreps();
         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "IrrepSizes", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Isizes);
    
            //Attributes
            hid_t attribute_space_id1  = H5Screate(H5S_SCALAR);
            hid_t attribute_id1        = H5Acreate(dataset_id, "nGroup", H5T_STD_I32LE, attribute_space_id1, H5P_DEFAULT, H5P_DEFAULT);
            int nGroup                 = SymmInfo.getGroupNumber();
            H5Awrite(attribute_id1, H5T_NATIVE_INT, &nGroup); 

            hid_t attribute_space_id2  = H5Screate(H5S_SCALAR);
            hid_t attribute_id2        = H5Acreate(dataset_id, "nIrreps", H5T_STD_I32LE, attribute_space_id2, H5P_DEFAULT, H5P_DEFAULT);
            int nIrreps                = SymmInfo.getNumberOfIrreps();
            H5Awrite(attribute_id2, H5T_NATIVE_INT, &nIrreps); 

            H5Aclose(attribute_id1);
            H5Aclose(attribute_id2);
            H5Sclose(attribute_space_id1);
            H5Sclose(attribute_space_id2);
    
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
      
      //The object itself.
      for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
      
         if (Isizes[cnt]>0) {
      
            std::stringstream sstream;
            sstream << "/TwoIndex" << cnt ;
            hid_t group_id3 = H5Gcreate(file_id, sstream.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
               hsize_t dimarray3       = Isizes[cnt]*(Isizes[cnt]+1)/2;
               hid_t dataspace_id3     = H5Screate_simple(1, &dimarray3, NULL);
               hid_t dataset_id3       = H5Dcreate(group_id3, "Matrix elements", H5T_IEEE_F64LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
               H5Dwrite(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[cnt]);
            
               H5Dclose(dataset_id3);
               H5Sclose(dataspace_id3);

            H5Gclose(group_id3);
         
         }
         
      }
      
   H5Fclose(file_id);

}

void CheMPS2::TwoIndex::read(const std::string name){
 
   //The hdf5 file
   hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The metadata
      hid_t group_id = H5Gopen(file_id, "/MetaData",H5P_DEFAULT);
       
         //The IrrepSizes
         hid_t dataset_id = H5Dopen(group_id, "IrrepSizes", H5P_DEFAULT);
    
            //Attributes
            hid_t attribute_id1 = H5Aopen_by_name(group_id,"IrrepSizes", "nGroup", H5P_DEFAULT, H5P_DEFAULT);
            int nGroup;
            H5Aread(attribute_id1, H5T_NATIVE_INT, &nGroup);
            assert( nGroup==SymmInfo.getGroupNumber() );
            
            hid_t attribute_id2 = H5Aopen_by_name(group_id,"IrrepSizes", "nIrreps", H5P_DEFAULT, H5P_DEFAULT);
            int nIrreps;
            H5Aread(attribute_id2, H5T_NATIVE_INT, &nIrreps);
            assert( nIrreps==SymmInfo.getNumberOfIrreps() );

            H5Aclose(attribute_id1);
            H5Aclose(attribute_id2);
    
         int * IsizesAgain = new int[SymmInfo.getNumberOfIrreps()];
         H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
         for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
            assert( IsizesAgain[cnt]==Isizes[cnt] );
         }
         delete [] IsizesAgain;
         H5Dclose(dataset_id);

      H5Gclose(group_id);
      
      //The object itself.
      for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
      
         if (Isizes[cnt]>0) {
      
            std::stringstream sstream;
            sstream << "/TwoIndex" << cnt ;
            hid_t group_id3 = H5Gopen(file_id, sstream.str().c_str(), H5P_DEFAULT);
            
               hid_t dataset_id3 = H5Dopen(group_id3, "Matrix elements", H5P_DEFAULT);
               H5Dread(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[cnt]);
               H5Dclose(dataset_id3);

            H5Gclose(group_id3);
         
         }
         
      }
      
   H5Fclose(file_id);

}



