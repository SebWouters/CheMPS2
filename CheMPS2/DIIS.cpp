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
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include <hdf5.h>

#include "Lapack.h"
#include "DIIS.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

CheMPS2::DIIS::DIIS(const int numVarsParamIn, const int numVarsErrorIn, const int numVecsIn){

   numVarsParam = numVarsParamIn;
   numVarsError = numVarsErrorIn;
   numVecs = numVecsIn;
   
   errorVectors = new double*[numVecs];
   paramVectors = new double*[numVecs];
   currentNumVecs = 0;
   
   lastLinco = new double[numVarsParam];

}

CheMPS2::DIIS::~DIIS(){

   for (int cnt=0; cnt<currentNumVecs; cnt++){
      delete [] errorVectors[cnt];
      delete [] paramVectors[cnt];
   }

   delete [] errorVectors;
   delete [] paramVectors;
   delete [] lastLinco;
      
}

int CheMPS2::DIIS::getNumVarsParam() const{ return numVarsParam; }

int CheMPS2::DIIS::getNumVarsError() const{ return numVarsError; }
         
int CheMPS2::DIIS::getNumVecs() const{ return numVecs; }
         
int CheMPS2::DIIS::getCurrentNumVecs() const{ return currentNumVecs; }

double * CheMPS2::DIIS::getLastLinco(){ return lastLinco; }

void CheMPS2::DIIS::appendNew(double * newError, double * newParam){

   if (currentNumVecs==numVecs){
   
      double * ptrE = errorVectors[0];
      double * ptrP = paramVectors[0];
      
      for (int cnt=1; cnt<numVecs; cnt++){
         errorVectors[cnt-1] = errorVectors[cnt];
         paramVectors[cnt-1] = paramVectors[cnt];
      }
      
      errorVectors[numVecs-1] = ptrE;
      paramVectors[numVecs-1] = ptrP;
   
   } else {
   
      errorVectors[currentNumVecs] = new double[numVarsError];
      paramVectors[currentNumVecs] = new double[numVarsParam];
      currentNumVecs += 1;
   
   }

   int inc = 1;
   dcopy_(&numVarsError, newError, &inc, errorVectors[currentNumVecs-1], &inc);
   dcopy_(&numVarsParam, newParam, &inc, paramVectors[currentNumVecs-1], &inc);

}

void CheMPS2::DIIS::calculateParam(double * newParam){

   int lindim = currentNumVecs + 1;
   int inc = 1;
   
   //The symmetric Pulay matrix
   double * matrix = new double[ lindim * lindim ];
   matrix[currentNumVecs*(1+lindim)] = 0.0;
   for (int cnt=0; cnt<currentNumVecs; cnt++){
      matrix[currentNumVecs + cnt*lindim] = 1.0 ;
      matrix[cnt + currentNumVecs*lindim] = 1.0 ;
      for (int cnt2=cnt; cnt2<currentNumVecs; cnt2++){
         matrix[cnt + cnt2*lindim] = ddot_(&numVarsError, errorVectors[cnt], &inc, errorVectors[cnt2], &inc);
         matrix[cnt2 + cnt*lindim] = matrix[cnt + cnt2*lindim];
      }
   }
   
   //Eigenvalue decomposition of the Pulay matrix
   char jobz = 'V';
   char uplo = 'U';
   double * array = new double[ lindim ]; //eigs
   int lwork = 3 * lindim;
   double * work = new double[ lwork ];
   int info;
   dsyev_(&jobz, &uplo, &lindim, matrix, &lindim, array, work, &lwork, &info);
   
   //The desired coefficients: [vec(c), lambda] = V (1/eigs) V^T [vec(0), 1]
   //Step 1: [vec(0), 1] --> work
   for (int cnt=0; cnt<currentNumVecs; cnt++){ work[cnt] = 0.0; }
   work[currentNumVecs] = 1.0;
   
   //Step 2: V^T [vec(0), 1] --> work+lindim
   char trans = 'T';
   char notra = 'N';
   int one = 1;
   double alpha = 1.0;
   double beta = 0.0;
   dgemm_(&trans, &notra, &lindim, &one, &lindim, &alpha, matrix, &lindim, work, &lindim, &beta, work+lindim, &lindim);
   
   //Step 3: (1/eigs) V^T [vec(0), 1] --> work+lindim
   for (int cnt=0; cnt<lindim; cnt++){ work[lindim + cnt] /= array[cnt]; }
   
   //Step 4: V (1/eigs) V^T [vec(0), 1] --> work
   dgemm_(&notra, &notra, &lindim, &one, &lindim, &alpha, matrix, &lindim, work+lindim, &lindim, &beta, work, &lindim);
   
   //Fill newParam and make a copy in lastLinco
   for (int cnt=0; cnt<numVarsParam; cnt++){ newParam[cnt] = 0.0; }
   for (int cnt=0; cnt<currentNumVecs; cnt++){ daxpy_(&numVarsParam, work+cnt, paramVectors[cnt], &inc, newParam, &inc); }
   dcopy_(&numVarsParam, newParam, &inc, lastLinco, &inc);
   
   //Print out the DIIS coefficients
   cout << "   DIIS::calculateParam : coefficients (newer vectors --> older vectors) : ";
   for (int cnt=0; cnt<currentNumVecs; cnt++){ cout << work[currentNumVecs-1-cnt] << "\t"; }
   cout << endl;
   
   delete [] matrix;
   delete [] array;
   delete [] work;

}


void CheMPS2::DIIS::saveDIIS(const string filename) const{

   hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   
   //The number of vectors saved
   hsize_t dimarray1       = 1;
   hid_t dataspace_id1     = H5Screate_simple(1, &dimarray1, NULL);
   hid_t dataset_id1       = H5Dcreate(group_id, "currentNumVecs", H5T_STD_I32LE, dataspace_id1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &currentNumVecs);
   H5Dclose(dataset_id1);
   H5Sclose(dataspace_id1);
   
   //The parameter vector sizes
   hsize_t dimarray2       = 1;
   hid_t dataspace_id2     = H5Screate_simple(1, &dimarray2, NULL);
   hid_t dataset_id2       = H5Dcreate(group_id, "numVarsParam", H5T_STD_I32LE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numVarsParam);
   H5Dclose(dataset_id2);
   H5Sclose(dataspace_id2);
   
   //The error vector sizes
   hsize_t dimarray5       = 1;
   hid_t dataspace_id5     = H5Screate_simple(1, &dimarray5, NULL);
   hid_t dataset_id5       = H5Dcreate(group_id, "numVarsError", H5T_STD_I32LE, dataspace_id5, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(dataset_id5, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numVarsError);
   H5Dclose(dataset_id5);
   H5Sclose(dataspace_id5);
   
   for (int cnt=0; cnt<currentNumVecs; cnt++){
   
      //The error vectors
      std::stringstream nameE;
      nameE << "error_" << cnt;
      hsize_t dimarray3      = numVarsError;
      hid_t dataspace_id3    = H5Screate_simple(1, &dimarray3, NULL);
      hid_t dataset_id3      = H5Dcreate(group_id, nameE.str().c_str(), H5T_IEEE_F64LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, errorVectors[cnt]);
      H5Dclose(dataset_id3);
      H5Sclose(dataspace_id3);
      
      //The parameter vectors
      std::stringstream nameP;
      nameP << "param_" << cnt;
      hsize_t dimarray4      = numVarsParam;
      hid_t dataspace_id4    = H5Screate_simple(1, &dimarray4, NULL);
      hid_t dataset_id4      = H5Dcreate(group_id, nameP.str().c_str(), H5T_IEEE_F64LE, dataspace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, paramVectors[cnt]);
      H5Dclose(dataset_id4);
      H5Sclose(dataspace_id4);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::DIIS::loadDIIS(const string filename){

   hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
   
   //The parameter vector sizes
   int numVarsBIS;
   hid_t dataset_id2 = H5Dopen(group_id, "numVarsParam", H5P_DEFAULT);
   H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numVarsBIS);
   H5Dclose(dataset_id2);
   assert( numVarsParam==numVarsBIS );
   
   //The error vector sizes
   hid_t dataset_id5 = H5Dopen(group_id, "numVarsError", H5P_DEFAULT);
   H5Dread(dataset_id5, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numVarsBIS);
   H5Dclose(dataset_id5);
   assert( numVarsError==numVarsBIS );
   
   //The number of vectors saved
   int currentNumVecsBIS;
   hid_t dataset_id1 = H5Dopen(group_id, "currentNumVecs", H5P_DEFAULT);
   H5Dread(dataset_id1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &currentNumVecsBIS);
   H5Dclose(dataset_id1);
   assert( currentNumVecsBIS<=numVecs );
   assert( currentNumVecsBIS>=0 );
   
   //Create just enough storage for the vectors to be loaded
   if (currentNumVecs < currentNumVecsBIS){
      for (int cnt=currentNumVecs; cnt<currentNumVecsBIS; cnt++){
         errorVectors[cnt] = new double[numVarsError];
         paramVectors[cnt] = new double[numVarsParam];
      }
      currentNumVecs = currentNumVecsBIS;
   }
   
   if (currentNumVecs > currentNumVecsBIS){
      for (int cnt=currentNumVecs; cnt>currentNumVecsBIS; cnt--){
         delete [] errorVectors[cnt-1];
         delete [] paramVectors[cnt-1];
      }
      currentNumVecs = currentNumVecsBIS;
   }
       
   for (int cnt=0; cnt<currentNumVecs; cnt++){
   
      //The error vectors
      std::stringstream nameE;
      nameE << "error_" << cnt;
      hid_t dataset_id3 = H5Dopen(group_id, nameE.str().c_str(), H5P_DEFAULT);
      H5Dread(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, errorVectors[cnt]);
      H5Dclose(dataset_id3);
      
      //The parameter vectors
      std::stringstream nameP;
      nameP << "param_" << cnt;
      hid_t dataset_id4 = H5Dopen(group_id, nameP.str().c_str(), H5P_DEFAULT);
      H5Dread(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, paramVectors[cnt]);
      H5Dclose(dataset_id4);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

