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
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>

#include "MyHDF5.h"
#include "Lapack.h"
#include "DMRGSCFunitary.h"

using std::string;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;

CheMPS2::DMRGSCFunitary::DMRGSCFunitary(DMRGSCFindices * iHandlerIn){

   this->iHandler = iHandlerIn;
   
   //Allocate the xmatrix
   x_linearlength = 0;
   xmatrix = new double**[ iHandler->getNirreps() ];
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      xmatrix[irrep] = new double*[3];
      for (int geval=0; geval<3; geval++){
         int size = 0;
         if (geval==0){ size = iHandler->getNOCC( irrep) * iHandler->getNDMRG(irrep); }
         if (geval==1){ size = iHandler->getNDMRG(irrep) * iHandler->getNVIRT(irrep); }
         if (geval==2){ size = iHandler->getNOCC( irrep) * iHandler->getNVIRT(irrep); }
         xmatrix[irrep][geval] = new double[size];
         for (int cnt=0; cnt<size; cnt++){ xmatrix[irrep][geval][cnt] = 0.0; }
         x_linearlength += size;
      }
   }
   
   //Allocate the unitary
   unitary = new double*[ iHandler->getNirreps() ];
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      const int linsize = iHandler->getNORB(irrep);
      const int size = linsize * linsize;
      unitary[irrep] = new double[size];
      for (int cnt=0; cnt<size; cnt++){ unitary[irrep][cnt] = 0.0; }
      for (int cnt=0; cnt<linsize; cnt++){ unitary[irrep][cnt*(1+linsize)] = 1.0; }
   }
   
   //Find the corresponding indices
   x_firstindex  = new int[x_linearlength];
   x_secondindex = new int[x_linearlength];
   int x_linearlength2 = 0;
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      for (int geval=0; geval<3; geval++){
         if (geval==0){
            for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
               for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){ //DMRG is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNDMRGstart(irrep)+ cntDMRG;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNOCCstart(irrep) + cntOcc;
                  x_linearlength2++;
               }
            }
         }
         if (geval==1){
            for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
               for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNDMRGstart(irrep) + cntDMRG;
                  x_linearlength2++;
               }
            }
         }
         if (geval==2){
            for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
               for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNOCCstart(irrep) + cntOcc;
                  x_linearlength2++;
               }
            }
         }
      }
   }
   if (x_linearlength != x_linearlength2){ cerr << "DMRGSCFunitary::DMRGSCFunitary : The number of variables is different!" << endl; }

}

CheMPS2::DMRGSCFunitary::~DMRGSCFunitary(){

   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      for (int geval=0; geval<3; geval++){ delete [] xmatrix[irrep][geval]; }
      delete [] xmatrix[irrep];
   }
   delete [] xmatrix;
      
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){ delete [] unitary[irrep]; }
   delete [] unitary;
      
   delete [] x_firstindex;
   delete [] x_secondindex;
      
}

int CheMPS2::DMRGSCFunitary::getNumVariablesX() const{ return x_linearlength; }

int CheMPS2::DMRGSCFunitary::getLinearIndex(const int p_index, const int q_index) const{

   for (int cnt=0; cnt<x_linearlength; cnt++){
      if ((p_index == x_firstindex[cnt]) && (q_index == x_secondindex[cnt])){ return cnt; }
   }
   
   return -1;

}

int CheMPS2::DMRGSCFunitary::getFirstIndex(const int linearindex) const{ return x_firstindex[linearindex]; }

int CheMPS2::DMRGSCFunitary::getSecondIndex(const int linearindex) const{ return x_secondindex[linearindex]; }

double * CheMPS2::DMRGSCFunitary::getBlock(const int irrep){ return unitary[irrep]; }

void CheMPS2::DMRGSCFunitary::copyXsolutionBack(double * vector){

   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
         
      for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
         for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
            int index1 = iHandler->getOrigNDMRGstart(irrep) + cntDMRG;
            int index2 = iHandler->getOrigNOCCstart(irrep) + cntOcc;
            int xsolindex = getLinearIndex(index1,index2);
            if (xsolindex==-1){ cerr << "DMRGSCFunitary::copyXsolutionBack : xsolindex==-1" << endl; }
            xmatrix[irrep][0][ cntDMRG + iHandler->getNDMRG(irrep) * cntOcc ] = vector[xsolindex];
         }
      }
      for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
         for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){
            int index1 = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
            int index2 = iHandler->getOrigNDMRGstart(irrep) + cntDMRG;
            int xsolindex = getLinearIndex(index1,index2);
            if (xsolindex==-1){ cerr << "DMRGSCFunitary::copyXsolutionBack : xsolindex==-1" << endl; }
            xmatrix[irrep][1][ cntVirt + iHandler->getNVIRT(irrep) * cntDMRG ] = vector[xsolindex];
         }
      }
      for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
         for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){
            int index1 = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
            int index2 = iHandler->getOrigNOCCstart(irrep) + cntOcc;
            int xsolindex = getLinearIndex(index1,index2);
            if (xsolindex==-1){ cerr << "DMRGSCFunitary::copyXsolutionBack : xsolindex==-1" << endl; }
            xmatrix[irrep][2][ cntVirt + iHandler->getNVIRT(irrep) * cntOcc ] = vector[xsolindex];
         }
      }
      
   }

}

void CheMPS2::DMRGSCFunitary::updateUnitary(double * temp1, double * temp2){

   //Per irrep
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){

      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;
      
      const int NOCC  = iHandler->getNOCC( irrep);
      const int NDMRG = iHandler->getNDMRG(irrep);
      const int NVIRT = iHandler->getNVIRT(irrep);
      
      if (linsize>1){ //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok
         
         //Construct the anti-symmetric x-matrix
         double * xblock = temp1;
         for (int cnt=0; cnt<size; cnt++){ xblock[cnt] = 0.0; }
         for (int cntOcc=0; cntOcc<NOCC; cntOcc++){
            for (int cntDMRG=0; cntDMRG<NDMRG; cntDMRG++){
               xblock[ NOCC + cntDMRG + linsize*cntOcc ]   =   xmatrix[irrep][0][ cntDMRG + NDMRG * cntOcc];
               xblock[ cntOcc + linsize*(NOCC + cntDMRG) ] = - xblock[ NOCC + cntDMRG + linsize*cntOcc ];
            }
         }
         for (int cntDMRG=0; cntDMRG<NDMRG; cntDMRG++){
            for (int cntVirt=0; cntVirt<NVIRT; cntVirt++){
               xblock[ NOCC + NDMRG + cntVirt + linsize*(NOCC + cntDMRG )] =   xmatrix[irrep][1][ cntVirt + NVIRT * cntDMRG ];
               xblock[ NOCC + cntDMRG + linsize*(NOCC + NDMRG + cntVirt )] = - xblock[ NOCC + NDMRG + cntVirt + linsize*( NOCC + cntDMRG ) ];
            }
         }
         for (int cntOcc=0; cntOcc<NOCC; cntOcc++){
            for (int cntVirt=0; cntVirt<NVIRT; cntVirt++){
               xblock[ NOCC + NDMRG + cntVirt + linsize*cntOcc ]     =   xmatrix[irrep][2][ cntVirt + NVIRT * cntOcc ];
               xblock[ cntOcc + linsize*( NOCC + NDMRG + cntVirt ) ] = - xblock[ NOCC + NDMRG + cntVirt + linsize*cntOcc ];
            }
         }
         
         //Calculate its eigenvalues and eigenvectors
         double * Bmat = temp2;
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0; //SET !!!
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,xblock,&linsize,&beta,Bmat,&linsize); //Bmat = xblock * xblock
         
         char uplo = 'U';
         char jobz = 'V';
         double * eigenval = temp1 + size;
         double * work = eigenval + linsize;
         int lwork = 2*size - linsize; //For linsize=2, lwork is 6 and should be 3*linsize-1=5, i.e. larger than linsize^2.
         int info;
         dsyev_(&jobz, &uplo, &linsize, Bmat, &linsize, eigenval, work, &lwork, &info); // xblock * xblock = Bmat * eigenval * Bmat^T
         
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,Bmat,&linsize,&beta,work,&linsize);
         char trans = 'T';
         double * work2 = temp2 + size;
         dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work,&linsize,&beta,work2,&linsize); //work2 = Bmat^T * xblock * Bmat
         
         if (CheMPS2::CASSCF_debugPrint){
            cout << "lambdas of irrep block " << irrep << " : " << endl;
            for (int cnt=0; cnt<linsize/2; cnt++){
               cout << "   block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
               cout << "           [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
            }
         }
         
         for (int cnt=0; cnt<linsize/2; cnt++){
            eigenval[cnt] = 0.5*( work2[2*cnt + linsize*(2*cnt+1)] - work2[2*cnt+1 + linsize*(2*cnt)] );
            work2[2*cnt + linsize*(2*cnt+1)] -= eigenval[cnt];
            work2[2*cnt+1 + linsize*(2*cnt)] += eigenval[cnt];
         }
         
         if (CheMPS2::CASSCF_debugPrint){
            double TwoNormResidual = 0.0;
            for (int cnt=0; cnt<size; cnt++){ TwoNormResidual += work2[cnt] * work2[cnt]; }
            TwoNormResidual = sqrt(TwoNormResidual);
            cout << "TwoNormResidual of irrep block " << irrep << " = " << TwoNormResidual << endl;
         }
         
         //Calculate exp(x)
         for (int cnt=0; cnt<size; cnt++){ work2[cnt] = 0.0; }
         for (int cnt=0; cnt<linsize/2; cnt++){
            double cosine = cos(eigenval[cnt]);
            double sine = sin(eigenval[cnt]);
            work2[2*cnt   + linsize*(2*cnt  )] = cosine;
            work2[2*cnt+1 + linsize*(2*cnt+1)] = cosine;
            work2[2*cnt   + linsize*(2*cnt+1)] = sine;
            work2[2*cnt+1 + linsize*(2*cnt  )] = - sine;
         }
         for (int cnt=2*(linsize/2); cnt<linsize; cnt++){
            work2[cnt*(linsize + 1)] = 1.0;
         }
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work2,&linsize,&beta,work,&linsize);
         dgemm_(&notr,&trans,&linsize,&linsize,&linsize,&alpha,work,&linsize,Bmat,&linsize,&beta,work2,&linsize); //work2 = exp(xblock)
         
         //U <-- exp(x) * U
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,work2,&linsize,unitary[irrep],&linsize,&beta,work,&linsize);
         int inc = 1;
         dcopy_(&size, work, &inc, unitary[irrep], &inc);

      }
   }
   
   if (CheMPS2::CASSCF_debugPrint){ CheckDeviationFromUnitary(temp2); }

}

void CheMPS2::DMRGSCFunitary::rotateUnitaryNOeigenvecs(double * eigenvecs, double * work){

   int passed = 0;
   int nOrbDMRG = iHandler->getDMRGcumulative(iHandler->getNirreps());
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){

      const int NDMRG = iHandler->getNDMRG(irrep);
      if (NDMRG > 1){

         int rotationlinsize = NDMRG;
         int blocklinsize = iHandler->getNORB(irrep);

         double * temp1 = work;
         double * temp2 = work + rotationlinsize*blocklinsize;
         double * BlockEigen = eigenvecs + passed * (nOrbDMRG + 1);
      
         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               temp1[row + rotationlinsize*col] = unitary[irrep][ iHandler->getNOCC(irrep) + row + blocklinsize * col ];
            }
         }

         char tran = 'T';
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0;
         dgemm_(&tran,&notr,&rotationlinsize,&blocklinsize,&rotationlinsize,&alpha,BlockEigen,&nOrbDMRG,temp1,&rotationlinsize,&beta,temp2,&rotationlinsize);

         for (int row = 0; row<rotationlinsize; row++){
            for (int col = 0; col<blocklinsize; col++){
               unitary[irrep][ iHandler->getNOCC(irrep) + row + blocklinsize * col ] = temp2[row + rotationlinsize*col];
            }
         }

      }

      passed += NDMRG;

   }
   
   if (CheMPS2::CASSCF_debugPrint){ CheckDeviationFromUnitary(work); }

}

void CheMPS2::DMRGSCFunitary::CheckDeviationFromUnitary(double * work) const{

   char tran = 'T';
   char notr = 'N';
   double alpha = 1.0;
   double beta = 0.0;
   
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      dgemm_(&tran,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep],&linsize,unitary[irrep],&linsize,&beta,work,&linsize);
      double value = 0.0;
      for (int cnt=0; cnt<linsize; cnt++){
         value += (work[cnt*(1+linsize)]-1.0) * (work[cnt*(1+linsize)]-1.0);
         for (int cnt2=cnt+1; cnt2<linsize; cnt2++){
            value += work[cnt + cnt2*linsize] * work[cnt + cnt2*linsize] + work[cnt2 + cnt*linsize] * work[cnt2 + cnt*linsize];
         }
      }
      value = sqrt(value);
      cout << "Two-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] - I = " << value << endl;
      
   }

}

void CheMPS2::DMRGSCFunitary::saveU() const{

   hid_t file_id = H5Fcreate(CheMPS2::CASSCF_unitaryStorageName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      std::stringstream irrepname;
      irrepname << "irrep_" << irrep;
      
      hsize_t dimarray      = iHandler->getNORB(irrep) * iHandler->getNORB(irrep);
      hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
      hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep]);

      H5Dclose(dataset_id);
      H5Sclose(dataspace_id);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::DMRGSCFunitary::loadU(){

   hid_t file_id = H5Fopen(CheMPS2::CASSCF_unitaryStorageName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      std::stringstream irrepname;
      irrepname << "irrep_" << irrep;

      hid_t dataset_id = H5Dopen(group_id, irrepname.str().c_str(), H5P_DEFAULT);
      H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep]);
         
      H5Dclose(dataset_id);
      
   }

   H5Gclose(group_id);
   H5Fclose(file_id);

}

void CheMPS2::DMRGSCFunitary::deleteStoredUnitary() const{

   std::stringstream temp;
   temp << "rm " << CheMPS2::CASSCF_unitaryStorageName;
   int info = system(temp.str().c_str());
   cout << "Info on CASSCF::Unitary rm call to system: " << info << endl;

}



