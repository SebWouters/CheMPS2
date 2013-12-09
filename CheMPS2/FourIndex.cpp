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

#include <hdf5.h>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "FourIndex.h"
#include "Lapack.h"

using namespace std;

CheMPS2::FourIndex::FourIndex(const int nGroup, const int * IrrepSizes){

   SymmInfo.setGroup(nGroup);
   
   Isizes = new int[SymmInfo.getNumberOfIrreps()];
   storage = new double*****[SymmInfo.getNumberOfIrreps()];
   
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      Isizes[Icenter] = IrrepSizes[Icenter];
   }
   
   //For the following: see text above storage in Fourindex.h
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      storage[Icenter] = new double****[SymmInfo.getNumberOfIrreps()];
      for (int I_i=0; I_i<SymmInfo.getNumberOfIrreps(); I_i++){
         int I_j = SymmInfo.directProd(Icenter,I_i);
         if ((Isizes[I_i]>0)&&(Isizes[I_j]>0)){
            storage[Icenter][I_i] = new double***[SymmInfo.getNumberOfIrreps()];
            for (int I_k=I_i; I_k<SymmInfo.getNumberOfIrreps(); I_k++){
               int I_l = SymmInfo.directProd(Icenter,I_k);
               if ((Isizes[I_k]>0)&&(Isizes[I_l]>0)){
                  if ((I_i <= I_j) && (I_j <= I_l)){
                     if (Icenter == 0){ // I_i = I_j and I_k = I_l
                        if (I_i == I_k){
                           storage[Icenter][I_i][I_k] = new double**[Isizes[I_i]*(Isizes[I_i]+1)/2];
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 storage[Icenter][I_i][I_k][i + k*(k+1)/2] = new double*[Isizes[I_j]-i];
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    if (i==j) storage[Icenter][I_i][I_k][i + k*(k+1)/2][j-i] = new double[Isizes[I_l]-k]; 
                                    else      storage[Icenter][I_i][I_k][i + k*(k+1)/2][j-i] = new double[Isizes[I_l]-j]; 
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           storage[Icenter][I_i][I_k] = new double**[Isizes[I_i]*Isizes[I_k]];
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 storage[Icenter][I_i][I_k][i + k*Isizes[I_i]] = new double*[Isizes[I_j]-i];
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    if (i==j) storage[Icenter][I_i][I_k][i+k*Isizes[I_i]][j-i] = new double[Isizes[I_l]-k];
                                    else      storage[Icenter][I_i][I_k][i+k*Isizes[I_i]][j-i] = new double[Isizes[I_l]];
                                 }
                              }
                           }
                        }
                     } else { //Icenter !=0 ; I_i < I_j and I_k != I_l
                        if (I_i == I_k){
                           storage[Icenter][I_i][I_k] = new double**[Isizes[I_i]*(Isizes[I_i]+1)/2];
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 storage[Icenter][I_i][I_k][i + k*(k+1)/2] = new double*[Isizes[I_j]];
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    storage[Icenter][I_i][I_k][i + k*(k+1)/2][j] = new double[Isizes[I_l]-j];
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           storage[Icenter][I_i][I_k] = new double**[Isizes[I_i]*Isizes[I_k]];
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 storage[Icenter][I_i][I_k][i + k*Isizes[I_i]] = new double*[Isizes[I_j]];
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j] = new double[Isizes[I_l]];
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
   }
   
}

CheMPS2::FourIndex::~FourIndex(){
   
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      for (int I_i=0; I_i<SymmInfo.getNumberOfIrreps(); I_i++){
         int I_j = SymmInfo.directProd(Icenter,I_i);
         if ((Isizes[I_i]>0)&&(Isizes[I_j]>0)){
            for (int I_k=I_i; I_k<SymmInfo.getNumberOfIrreps(); I_k++){
               int I_l = SymmInfo.directProd(Icenter,I_k);
               if ((Isizes[I_k]>0)&&(Isizes[I_l]>0)){
                  if ((I_i <= I_j) && (I_j <= I_l)){
                     if (Icenter == 0){
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++) delete [] storage[Icenter][I_i][I_k][i + k*(k+1)/2][j-i];
                                 delete [] storage[Icenter][I_i][I_k][i + k*(k+1)/2];
                              }
                           }
                        } else {
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++) delete [] storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j-i];
                                 delete [] storage[Icenter][I_i][I_k][i + k*Isizes[I_i]];
                              }
                           }
                        }
                     } else {
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++) delete [] storage[Icenter][I_i][I_k][i + k*(k+1)/2][j];
                                 delete [] storage[Icenter][I_i][I_k][i + k*(k+1)/2];
                              }
                           }
                        } else {
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++) delete [] storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j];
                                 delete [] storage[Icenter][I_i][I_k][i + k*Isizes[I_i]];
                              }
                           }
                        }
                     }
                     delete [] storage[Icenter][I_i][I_k];
                  }
               }
            }
            delete [] storage[Icenter][I_i];
         }
      }
      delete [] storage[Icenter];
   }
   
   delete [] Isizes;
   delete [] storage;
   
}

void CheMPS2::FourIndex::set(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val){

   getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)[0] = val;

}

void CheMPS2::FourIndex::add(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l, const double val){

   getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)[0] += val;

}

double CheMPS2::FourIndex::get(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const{

   return getPointer(irrep_i, irrep_j, irrep_k, irrep_l, i, j, k, l)[0];

}

double * CheMPS2::FourIndex::getPointer(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {

   if (SymmInfo.directProd(irrep_i,irrep_j)==SymmInfo.directProd(irrep_k,irrep_l)){

      if ((irrep_i <= irrep_j) && (irrep_i <= irrep_k) && (irrep_j <= irrep_l)){ // (ijkl irrep ordering)
         return getPtrIrrepOrderOK(irrep_i,irrep_j,irrep_k,irrep_l,i,j,k,l);
      }
   
      if ((irrep_j <= irrep_i) && (irrep_j <= irrep_l) && (irrep_i <= irrep_k)){ // (jilk irrep ordering)
         return getPtrIrrepOrderOK(irrep_j,irrep_i,irrep_l,irrep_k,j,i,l,k);
      }
      
      if ((irrep_k <= irrep_j) && (irrep_k <= irrep_i) && (irrep_j <= irrep_l)){ // (kjil irrep ordering)
         return getPtrIrrepOrderOK(irrep_k,irrep_j,irrep_i,irrep_l,k,j,i,l);
      }   
      
      if ((irrep_j <= irrep_k) && (irrep_j <= irrep_l) && (irrep_k <= irrep_i)){ // (jkli irrep ordering)
         return getPtrIrrepOrderOK(irrep_j,irrep_k,irrep_l,irrep_i,j,k,l,i);
      }   
      
      if ((irrep_i <= irrep_l) && (irrep_i <= irrep_k) && (irrep_l <= irrep_j)){ // (ilkj irrep ordering)
         return getPtrIrrepOrderOK(irrep_i,irrep_l,irrep_k,irrep_j,i,l,k,j);
      }   
      
      if ((irrep_l <= irrep_i) && (irrep_l <= irrep_j) && (irrep_i <= irrep_k)){ // (lijk irrep ordering)
         return getPtrIrrepOrderOK(irrep_l,irrep_i,irrep_j,irrep_k,l,i,j,k);
      }   
      
      if ((irrep_k <= irrep_l) && (irrep_k <= irrep_i) && (irrep_l <= irrep_j)){ // (klij irrep ordering)
         return getPtrIrrepOrderOK(irrep_k,irrep_l,irrep_i,irrep_j,k,l,i,j);
      }   
      
      if ((irrep_l <= irrep_k) && (irrep_l <= irrep_j) && (irrep_k <= irrep_i)){ // (lkji irrep ordering)
         return getPtrIrrepOrderOK(irrep_l,irrep_k,irrep_j,irrep_i,l,k,j,i);
      }
        
   }
   
   return NULL;

}

double * CheMPS2::FourIndex::getPtrIrrepOrderOK(const int irrep_i, const int irrep_j, const int irrep_k, const int irrep_l, const int i, const int j, const int k, const int l) const {

   //I_i <= I_j <= I_l and I_i <= I_k
   int Icenter = SymmInfo.directProd(irrep_i,irrep_j);
   
   if (Icenter>0){ // I_i < I_j and I_k != I_l
   
      if (irrep_i == irrep_k){ //and hence I_j == I_l
      
         if (k>=i){
            if (l>=j) return getPtrAllOK(5, Icenter, irrep_i, irrep_k, i, j, k, l);
            else      return getPtrAllOK(5, Icenter, irrep_i, irrep_k, i, l, k, j);
         } else {
            if (l>=j) return getPtrAllOK(5, Icenter, irrep_k, irrep_i, k, j, i, l);
            else      return getPtrAllOK(5, Icenter, irrep_k, irrep_i, k, l, i, j);
         }
      
      } else return getPtrAllOK(6, Icenter, irrep_i, irrep_k, i, j, k, l);
   
   } else {
   
      if (irrep_i == irrep_k){ // all irreps the same
      
         //i en j
         if ((i <  j) && (i <= k) && (j <= l)) return getPtrAllOK(2, Icenter, irrep_i, irrep_k, i, j, k, l); // (ijkl ordering) 
         if ((i == j) && (i <= k) && (j <= l)){
            if (l>=k) return getPtrAllOK(1, Icenter, irrep_i, irrep_k, i, j, k, l); // (ijkl ordering)
            else      return getPtrAllOK(1, Icenter, irrep_j, irrep_l, j, i, l, k); // (jilk ordering)
         }
         if ((j <  i) && (j <= l) && (i <= k)) return getPtrAllOK(2, Icenter, irrep_j, irrep_l, j, i, l, k); // (jilk ordering)
         
         //k en l
         if ((k <  l) && (k <= i) && (l <= j)) return getPtrAllOK(2, Icenter, irrep_k, irrep_i, k, l, i, j); // (klij ordering)
         if ((k == l) && (k <= i) && (l <= j)){
            if (j>=i) return getPtrAllOK(1, Icenter, irrep_k, irrep_i, k, l, i, j); // (klij ordering)
            else      return getPtrAllOK(1, Icenter, irrep_l, irrep_j, l, k, j, i); // (lkji ordering)
         }
         if ((l <  k) && (l <= j) && (k <= i)) return getPtrAllOK(2, Icenter, irrep_l, irrep_j, l, k, j, i); // (lkji ordering)

         // k en j
         if ((k <  j) && (k <= i) && (j <= l)) return getPtrAllOK(2, Icenter, irrep_k, irrep_i, k, j, i, l); // (kjil ordering)
         if ((k == j) && (k <= i) && (j <= l)){
            if (l>=i) return getPtrAllOK(1, Icenter, irrep_k, irrep_i, k, j, i, l); // (kjil ordering)
            else      return getPtrAllOK(1, Icenter, irrep_j, irrep_l, j, k, l, i); // (jkli ordering)
         }
         if ((j <  k) && (j <= l) && (k <= i)) return getPtrAllOK(2, Icenter, irrep_j, irrep_l, j, k, l, i); // (jkli ordering)
         
         // i en l
         if ((i <  l) && (i <= k) && (l <= j)) return getPtrAllOK(2, Icenter, irrep_i, irrep_k, i, l, k, j); // (ilkj ordering)
         if ((i == l) && (i <= k) && (l <= j)){
            if (j>=k) return getPtrAllOK(1, Icenter, irrep_i, irrep_k, i, l, k, j); // (ilkj ordering)
            else      return getPtrAllOK(1, Icenter, irrep_l, irrep_j, l, i, j, k); // (lijk ordering)
         }
         if ((l <  i) && (l <= j) && (i <= k)) return getPtrAllOK(2, Icenter, irrep_l, irrep_j, l, i, j, k); // (lijk ordering) 

      } else {
      
         if (j==i){
            if (l>=k) return getPtrAllOK(3, Icenter, irrep_i, irrep_k, i, j, k, l);
            else      return getPtrAllOK(3, Icenter, irrep_j, irrep_l, j, i, l, k);
         } else {
            if (j>i)  return getPtrAllOK(4, Icenter, irrep_i, irrep_k, i, j, k, l);
            else      return getPtrAllOK(4, Icenter, irrep_j, irrep_l, j, i, l, k);
         }
         
      }
      
   }
   
   return NULL;
   
}

double * CheMPS2::FourIndex::getPtrAllOK(const int number, const int Icent, const int irrep_i, const int irrep_k, const int i, const int j, const int k, const int l) const {

   switch (number){
      case 1:
         return storage[Icent][irrep_i][irrep_k][i + k*(k+1)/2][j-i] + l-k;
      case 2:
         return storage[Icent][irrep_i][irrep_k][i + k*(k+1)/2][j-i] + l-j;
      case 3:
         return storage[Icent][irrep_i][irrep_k][i + Isizes[irrep_i]*k][j-i] + l-k;
      case 4:
         return storage[Icent][irrep_i][irrep_k][i + Isizes[irrep_i]*k][j-i] + l;
      case 5:
         return storage[Icent][irrep_i][irrep_k][i + k*(k+1)/2][j] + l-j;
      case 6:
         return storage[Icent][irrep_i][irrep_k][i + Isizes[irrep_i]*k][j] + l;
   }
   
   return NULL;

}

void CheMPS2::FourIndex::save(const std::string name) const{

   //The object size.
   long long theTotalSize = 0;
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      for (int I_i=0; I_i<SymmInfo.getNumberOfIrreps(); I_i++){
         int I_j = SymmInfo.directProd(Icenter,I_i);
         if ((Isizes[I_i]>0)&&(Isizes[I_j]>0)){
            for (int I_k=I_i; I_k<SymmInfo.getNumberOfIrreps(); I_k++){
               int I_l = SymmInfo.directProd(Icenter,I_k);
               if ((Isizes[I_k]>0)&&(Isizes[I_l]>0)){
                  if ((I_i <= I_j) && (I_j <= I_l)){
                     if (Icenter == 0){ // I_i = I_j and I_k = I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    theTotalSize += Isizes[I_l] - ((i==j)?k:j);
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    theTotalSize += Isizes[I_l] - ((i==j)?k:0);
                                 }
                              }
                           }
                        }
                     } else { //Icenter !=0 ; I_i < I_j and I_k != I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    theTotalSize += Isizes[I_l] - j;
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    theTotalSize += Isizes[I_l];
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
   }
   
   //Space to copy it: needed for fast HDF5 write.
   double * writeDump = new double[theTotalSize];
   
   //Copy the object.
   long long theTotalSize2 = 0;
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      for (int I_i=0; I_i<SymmInfo.getNumberOfIrreps(); I_i++){
         int I_j = SymmInfo.directProd(Icenter,I_i);
         if ((Isizes[I_i]>0)&&(Isizes[I_j]>0)){
            for (int I_k=I_i; I_k<SymmInfo.getNumberOfIrreps(); I_k++){
               int I_l = SymmInfo.directProd(Icenter,I_k);
               if ((Isizes[I_k]>0)&&(Isizes[I_l]>0)){
                  if ((I_i <= I_j) && (I_j <= I_l)){
                     if (Icenter == 0){ // I_i = I_j and I_k = I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    int mySize = Isizes[I_l] - ((i==j)?k:j);
                                    int inc = 1;
                                    dcopy_(&mySize, storage[Icenter][I_i][I_k][i + k*(k+1)/2][j-i], &inc, writeDump + theTotalSize2, &inc);
                                    theTotalSize2 += mySize;
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                    int mySize = Isizes[I_l] - ((i==j)?k:0);
                                    int inc = 1;
                                    dcopy_(&mySize, storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j-i], &inc, writeDump + theTotalSize2, &inc);
                                    theTotalSize2 += mySize;
                                 }
                              }
                           }
                        }
                     } else { //Icenter !=0 ; I_i < I_j and I_k != I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    int mySize = Isizes[I_l] - j;
                                    int inc = 1;
                                    dcopy_(&mySize, storage[Icenter][I_i][I_k][i + k*(k+1)/2][j], &inc, writeDump + theTotalSize2, &inc);
                                    theTotalSize2 += mySize;
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                    int mySize = Isizes[I_l];
                                    int inc = 1;
                                    dcopy_(&mySize, storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j], &inc, writeDump + theTotalSize2, &inc);
                                    theTotalSize2 += mySize;
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
   }
   
   if (theTotalSize != theTotalSize2){ std::cerr << "FourIndex::save : mismatch of theTotalSize and theTotalSize2" << std::endl;}
 
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
            H5Awrite(attribute_id3, H5T_NATIVE_LLONG, &theTotalSize); 

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
      hid_t group_id7 = H5Gcreate(file_id, "/FourIndexObject", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            
         hsize_t dimarray7       = theTotalSize; //hsize_t is defined by default as unsigned long long, so no problem
         hid_t dataspace_id7     = H5Screate_simple(1, &dimarray7, NULL);
         hid_t dataset_id7       = H5Dcreate(group_id7, "Matrix elements", H5T_IEEE_F64LE, dataspace_id7, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, writeDump);
             
         H5Dclose(dataset_id7);
         H5Sclose(dataspace_id7);

      H5Gclose(group_id7);
      
   H5Fclose(file_id);
   
   delete [] writeDump;

}

void CheMPS2::FourIndex::read(const std::string name){
 
   //The hdf5 file
   hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The metadata
      hid_t group_id = H5Gopen(file_id, "/MetaData",H5P_DEFAULT);
       
         //The IrrepSizes
         hid_t dataset_id = H5Dopen(group_id, "IrrepSizes", H5P_DEFAULT);
    
            //Attributes
            hid_t attribute_id1 = H5Aopen_name(dataset_id, "nGroup");
            int nGroup;
            H5Aread(attribute_id1, H5T_NATIVE_INT, &nGroup);
            if (nGroup != SymmInfo.getGroupNumber()) std::cout << "Error at FourIndex::read : nGroup doesn't match." << std::endl;
            
            hid_t attribute_id2 = H5Aopen_name(dataset_id, "nIrreps");
            int nIrreps;
            H5Aread(attribute_id2, H5T_NATIVE_INT, &nIrreps);
            if (nGroup != SymmInfo.getGroupNumber()) std::cout << "Error at FourIndex::read : nIrreps doesn't match." << std::endl;
            
            hid_t attribute_id3 = H5Aopen_name(dataset_id, "theTotalSize");
            long long theTotalSize;
            H5Aread(attribute_id3, H5T_NATIVE_LLONG, &theTotalSize);

            H5Aclose(attribute_id1);
            H5Aclose(attribute_id2);
            H5Aclose(attribute_id3);
    
         int * IsizesAgain = new int[SymmInfo.getNumberOfIrreps()];
         H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
         for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
            if (IsizesAgain[cnt]!=Isizes[cnt]) std::cout << "Error at FourIndex::read : One of the Isizes doesn't match." << std::endl;
         }
         delete [] IsizesAgain;
         H5Dclose(dataset_id);

      H5Gclose(group_id);
      
      std::cout << "FourIndex::read : loading " << theTotalSize << " doubles." << std::endl;
      double * readDump = new double[theTotalSize];
      
      //The object itself.
      hid_t group_id7 = H5Gopen(file_id, "/FourIndexObject", H5P_DEFAULT);

      hid_t dataset_id7 = H5Dopen(group_id7, "Matrix elements", H5P_DEFAULT);
      H5Dread(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, readDump);
      H5Dclose(dataset_id7);

      H5Gclose(group_id7);
      
   H5Fclose(file_id);
   
   //Copy the matrix elements back
   long long theTotalSize2 = 0;
   for (int Icenter=0; Icenter<SymmInfo.getNumberOfIrreps(); Icenter++){
      for (int I_i=0; I_i<SymmInfo.getNumberOfIrreps(); I_i++){
         int I_j = SymmInfo.directProd(Icenter,I_i);
         if ((Isizes[I_i]>0)&&(Isizes[I_j]>0)){
            for (int I_k=I_i; I_k<SymmInfo.getNumberOfIrreps(); I_k++){
               int I_l = SymmInfo.directProd(Icenter,I_k);
               if ((Isizes[I_k]>0)&&(Isizes[I_l]>0)){
                  if ((I_i <= I_j) && (I_j <= I_l)){
                     if (Icenter == 0){ // I_i = I_j and I_k = I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                 
                                    double * location = storage[Icenter][I_i][I_k][i + k*(k+1)/2][j-i];
                                    int mySize = Isizes[I_l] - ((i==j)?k:j);
                                    int inc = 1;
                                    dcopy_(&mySize, readDump + theTotalSize2, &inc, location, &inc);
                                    theTotalSize2 += mySize;

                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=i; j<Isizes[I_j]; j++){
                                 
                                    double * location = storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j-i];
                                    int mySize = Isizes[I_l] - ((i==j)?k:0);
                                    int inc = 1;
                                    dcopy_(&mySize, readDump + theTotalSize2, &inc, location, &inc);
                                    theTotalSize2 += mySize;

                                 }
                              }
                           }
                        }
                     } else { //Icenter !=0 ; I_i < I_j and I_k != I_l
                        if (I_i == I_k){
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=i; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                 
                                    double * location = storage[Icenter][I_i][I_k][i + k*(k+1)/2][j];
                                    int mySize = Isizes[I_l] - j;
                                    int inc = 1;
                                    dcopy_(&mySize, readDump + theTotalSize2, &inc, location, &inc);
                                    theTotalSize2 += mySize;
                           
                                 }
                              }
                           }
                        } else { // I_i < I_k
                           for (int i=0; i<Isizes[I_i]; i++){
                              for (int k=0; k<Isizes[I_k]; k++){
                                 for (int j=0; j<Isizes[I_j]; j++){
                                 
                                    double * location = storage[Icenter][I_i][I_k][i + k*Isizes[I_i]][j];
                                    int mySize = Isizes[I_l];
                                    int inc = 1;
                                    dcopy_(&mySize, readDump + theTotalSize2, &inc, location, &inc);
                                    theTotalSize2 += mySize;
                           
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
   }
   
   if (theTotalSize != theTotalSize2){ std::cerr << "FourIndex::read : mismatch of theTotalSize and theTotalSize2" << std::endl; }
   
   delete [] readDump;

}

