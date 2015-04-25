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
#include <string>

#include "TwoDMstorage.h"
#include "MyHDF5.h"

using std::cout;
using std::endl;
using std::string;

CheMPS2::TwoDMstorage::TwoDMstorage(const int nGroup, const int * IrrepSizes){

   SymmInfo.setGroup(nGroup);
   
   Isizes = new int[SymmInfo.getNumberOfIrreps()];
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){ Isizes[Icenter] = IrrepSizes[Icenter]; }
   
   storage = new long long**[SymmInfo.getNumberOfIrreps()];
   
   arrayLength = calcNumberOfElements(true); //true means allocate the storage!
   theElements = new double[arrayLength];
   
}

long long CheMPS2::TwoDMstorage::calcNumberOfElements(const bool allocate){

   //The object size: see text above storage in TwoDMstorage.h
   long long theTotalSize = 0;
   const int NumIrreps = SymmInfo.getNumberOfIrreps();
   
   for (int Icenter=0; Icenter<NumIrreps; Icenter++){
   
      if (allocate){ storage[Icenter] = new long long * [NumIrreps]; } //ALLOCATION
      
      for (int I_i=0; I_i<NumIrreps; I_i++){
         const int I_j = Irreps::directProd(Icenter,I_i);
         if ((I_i <= I_j) && (Isizes[I_i]>0) && (Isizes[I_j]>0)){ //I_i <= I_j and both I_i and I_j have orbitals
         
            if (allocate){ storage[Icenter][I_i] = new long long[NumIrreps - I_i]; } //ALLOCATION
         
            for (int I_k=I_i; I_k<NumIrreps; I_k++){ //I_i <= I_k
               const int I_l = Irreps::directProd(Icenter,I_k);
               if ((I_i <= I_l) && (Isizes[I_k]>0) && (Isizes[I_l]>0)){ //I_i <= I_l and both I_k and I_l have orbitals
               
                  storage[Icenter][I_i][I_k - I_i] = theTotalSize;
                  theTotalSize += Isizes[I_i] * Isizes[I_j] * Isizes[I_k] * Isizes[I_l];
               
               }
            }
            if (!allocate){ delete [] storage[Icenter][I_i]; } //DEALLOCATION
         }
      }
      if (!allocate){ delete [] storage[Icenter]; } //DEALLOCATION
   }
   
   return theTotalSize;
   
}

CheMPS2::TwoDMstorage::~TwoDMstorage(){
   
   arrayLength = calcNumberOfElements(false); //false means delete the storage!
   delete [] theElements;
   delete [] Isizes;
   delete [] storage;
   
}

void CheMPS2::TwoDMstorage::set(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val){

   theElements[getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)] = val;

}

double CheMPS2::TwoDMstorage::get(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const{

   return theElements[getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)];

}

long long CheMPS2::TwoDMstorage::getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {

   const int Icenter = Irreps::directProd(irrep_i,irrep_j);

   if ( Icenter == Irreps::directProd(irrep_k,irrep_l) ){

      if ((irrep_i <= irrep_j) && (irrep_i <= irrep_k) && (irrep_i <= irrep_l)){ // Irrep I_i is the smallest : (ijkl)
         return storage[Icenter][irrep_i][irrep_k - irrep_i] + i + Isizes[irrep_i] * ( j + Isizes[irrep_j] * ( k + Isizes[irrep_k] * l ) );
      }
   
      if ((irrep_j <= irrep_i) && (irrep_j <= irrep_k) && (irrep_j <= irrep_l)){ // Irrep I_j is the smallest : (jilk)
         return storage[Icenter][irrep_j][irrep_l - irrep_j] + j + Isizes[irrep_j] * ( i + Isizes[irrep_i] * ( l + Isizes[irrep_l] * k ) );
      }
      
      if ((irrep_k <= irrep_i) && (irrep_k <= irrep_j) && (irrep_k <= irrep_l)){ // Irrep I_k is the smallest : (klij)
         return storage[Icenter][irrep_k][irrep_i - irrep_k] + k + Isizes[irrep_k] * ( l + Isizes[irrep_l] * ( i + Isizes[irrep_i] * j ) );
      }
      
      if ((irrep_l <= irrep_i) && (irrep_l <= irrep_j) && (irrep_l <= irrep_k)){ // Irrep I_l is the smallest : (lkji)
         return storage[Icenter][irrep_l][irrep_j - irrep_l] + l + Isizes[irrep_l] * ( k + Isizes[irrep_k] * ( j + Isizes[irrep_j] * i ) );
      }
        
   }
   
   return -1;

}

void CheMPS2::TwoDMstorage::save(const std::string name) const{
 
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
            
            hid_t attribute_space_id3  = H5Screate(H5S_SCALAR);
            hid_t attribute_id3        = H5Acreate(dataset_id, "theTotalSize", H5T_STD_I64LE, attribute_space_id3, H5P_DEFAULT, H5P_DEFAULT);
            H5Awrite(attribute_id3, H5T_NATIVE_LLONG, &arrayLength); 

            H5Aclose(attribute_id1);
            H5Aclose(attribute_id2);
            H5Aclose(attribute_id3);
            H5Sclose(attribute_space_id1);
            H5Sclose(attribute_space_id2);
            H5Sclose(attribute_space_id3);
    
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
      
      //The object itself
      hid_t group_id7 = H5Gcreate(file_id, "/TwoDMstorageObject", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
         hsize_t dimarray7       = arrayLength; //hsize_t is defined by default as unsigned long long, so no problem
         hid_t dataspace_id7     = H5Screate_simple(1, &dimarray7, NULL);
         hid_t dataset_id7       = H5Dcreate(group_id7, "TwoDM elements", H5T_IEEE_F64LE, dataspace_id7, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theElements);
             
         H5Dclose(dataset_id7);
         H5Sclose(dataspace_id7);

      H5Gclose(group_id7);
      
   H5Fclose(file_id);

}

void CheMPS2::TwoDMstorage::read(const std::string name){
 
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
            
            hid_t attribute_id3 = H5Aopen_by_name(group_id,"IrrepSizes", "theTotalSize", H5P_DEFAULT, H5P_DEFAULT);
            long long theTotalSize;
            H5Aread(attribute_id3, H5T_NATIVE_LLONG, &theTotalSize);
            assert( theTotalSize==arrayLength );

            H5Aclose(attribute_id1);
            H5Aclose(attribute_id2);
            H5Aclose(attribute_id3);
    
         int * IsizesAgain = new int[SymmInfo.getNumberOfIrreps()];
         H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
         for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
            assert( IsizesAgain[cnt]==Isizes[cnt] );
         }
         delete [] IsizesAgain;
         H5Dclose(dataset_id);

      H5Gclose(group_id);
      
      std::cout << "TwoDMstorage::read : loading " << arrayLength << " doubles." << std::endl;
      
      //The object itself.
      hid_t group_id7 = H5Gopen(file_id, "/TwoDMstorageObject", H5P_DEFAULT);

      hid_t dataset_id7 = H5Dopen(group_id7, "TwoDM elements", H5P_DEFAULT);
      H5Dread(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, theElements);
      H5Dclose(dataset_id7);

      H5Gclose(group_id7);
      
   H5Fclose(file_id);
   
   std::cout << "TwoDMstorage::read : everything loaded!" << std::endl;

}

