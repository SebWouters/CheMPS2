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

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>

#include "DMRG.h"

void CheMPS2::DMRG::saveMPS(const std::string name, TensorT ** MPSlocation, SyBookkeeper * BKlocation, bool isConverged) const{
 
   //The hdf5 file
   hid_t file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      //Whether the MPS was converged or not
      hid_t group_id = H5Gcreate(file_id, "/Convergence", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
         hsize_t dimarray     = 1; //One integer
         hid_t dataspace_id   = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id     = H5Dcreate(group_id, "Converged_yn", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         int toWrite = (isConverged)?1:0;
         H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &toWrite);
    
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);

      H5Gclose(group_id);
      
      //The current virtual dimensions
      for (int bound=0; bound<=BKlocation->gL(); bound++){
         for (int N=BKlocation->gNmin(bound); N<=BKlocation->gNmax(bound); N++){
            for (int TwoS=BKlocation->gTwoSmin(bound,N); TwoS<=BKlocation->gTwoSmax(bound,N); TwoS+=2){
               for (int Irrep=0; Irrep<BKlocation->getNumberOfIrreps(); Irrep++){
     
                  std::stringstream sstream;
                  sstream << "/VirtDim_" << bound << "_" << N << "_" << TwoS << "_" << Irrep;
                  hid_t group_id2 = H5Gcreate(file_id, sstream.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                  
                     hsize_t dimarray2     = 1; //One integer
                     hid_t dataspace_id2   = H5Screate_simple(1, &dimarray2, NULL);
                     hid_t dataset_id2     = H5Dcreate(group_id2, "Value", H5T_STD_I32LE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                     int toWrite2 = BKlocation->gCurrentDim(bound,N,TwoS,Irrep);
                     H5Dwrite(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &toWrite2);
                
                     H5Dclose(dataset_id2);
                     H5Sclose(dataspace_id2);

                  H5Gclose(group_id2);
                  
               }
            }
         }
      }
      
      //The MPS
      for (int site=0; site<BKlocation->gL(); site++){
         
         std::stringstream sstream;
         sstream << "/MPS_" << site;
         hid_t group_id3 = H5Gcreate(file_id, sstream.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         
            hsize_t dimarray3     = MPSlocation[site]->gKappa2index(MPSlocation[site]->gNKappa()); //An array of doubles
            hid_t dataspace_id3   = H5Screate_simple(1, &dimarray3, NULL);
            hid_t dataset_id3     = H5Dcreate(group_id3, "Values", H5T_IEEE_F64LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MPSlocation[site]->gStorage());
       
            H5Dclose(dataset_id3);
            H5Sclose(dataspace_id3);

         H5Gclose(group_id3);
         
      }
      
      
   H5Fclose(file_id);

}

void CheMPS2::DMRG::loadDIM(const std::string name, SyBookkeeper * BKlocation){

   //The hdf5 file
   hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The current virtual dimensions
      for (int bound=0; bound<=BKlocation->gL(); bound++){
         for (int N=BKlocation->gNmin(bound); N<=BKlocation->gNmax(bound); N++){
            for (int TwoS=BKlocation->gTwoSmin(bound,N); TwoS<=BKlocation->gTwoSmax(bound,N); TwoS+=2){
               for (int Irrep=0; Irrep<BKlocation->getNumberOfIrreps(); Irrep++){
     
                  std::stringstream sstream;
                  sstream << "/VirtDim_" << bound << "_" << N << "_" << TwoS << "_" << Irrep;
                  hid_t group_id2 = H5Gopen(file_id, sstream.str().c_str(), H5P_DEFAULT);
                  
                     hid_t dataset_id2 = H5Dopen(group_id2, "Value", H5P_DEFAULT);
                     int toRead;
                     H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &toRead);
                     BKlocation->SetDim(bound, N, TwoS, Irrep, toRead);
                     H5Dclose(dataset_id2);
                  
                  H5Gclose(group_id2);
                  
               }
            }
         }
      }
      
   H5Fclose(file_id);

}

void CheMPS2::DMRG::loadMPS(const std::string name, TensorT ** MPSlocation, bool * isConverged){

   //The hdf5 file
   hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //Whether the MPS was converged or not
      hid_t group_id = H5Gopen(file_id, "/Convergence", H5P_DEFAULT);
      
         hid_t dataset_id = H5Dopen(group_id, "Converged_yn", H5P_DEFAULT);
         int toRead;
         H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &toRead);
         isConverged[0] = (toRead==0)?false:true;
         H5Dclose(dataset_id);

      H5Gclose(group_id);
      
      //The MPS
      for (int site=0; site<L; site++){
         
         std::stringstream sstream;
         sstream << "/MPS_" << site;
         hid_t group_id3 = H5Gopen(file_id, sstream.str().c_str(), H5P_DEFAULT);
         
            hid_t dataset_id3     = H5Dopen(group_id3, "Values", H5P_DEFAULT);
            H5Dread(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MPSlocation[site]->gStorage());
            H5Dclose(dataset_id3);

         H5Gclose(group_id3);
         
      }
      
      
   H5Fclose(file_id);

}

void CheMPS2::DMRG::deleteStoredMPS(){

   std::stringstream thestream;
   thestream << "rm " << CheMPS2::DMRG_MPS_storage_prefix << "*.h5";
   int info = system(thestream.str().c_str());
   std::cout << "Info on DMRG::MPS rm call to system: " << info << std::endl;

}


