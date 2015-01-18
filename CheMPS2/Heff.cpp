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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "Heff.h"
#include "Lapack.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::Heff::Heff(const SyBookkeeper * denBKIn, const Problem * ProbIn){

   denBK = denBKIn;
   Prob = ProbIn;

}

CheMPS2::Heff::~Heff(){

}

void CheMPS2::Heff::makeHeff(double * memS, double * memHeff, const Sobject * denS, TensorL *** Ltensors, TensorA **** Atensors, TensorB **** Btensors, TensorC **** Ctensors, TensorD **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   const int DIM = max(denBK->gMaxDimAtBound(indexS), denBK->gMaxDimAtBound(indexS+2));
   
   //PARALLEL
   #pragma omp parallel
   {
   
      double * temp  = new double[DIM*DIM];
      double * temp2 = new double[DIM*DIM];
   
      #pragma omp for schedule(dynamic)
      for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
      
         const int ikappa = denS->gReorder(ikappaBIS);

         for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeff[cnt] = 0.0; }
         
         addDiagram1C(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
         addDiagram1D(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
         addDiagram2dall(ikappa, memS, memHeff, denS);
         addDiagram3Eand3H(ikappa, memS, memHeff, denS);
         addDiagramExcitations(ikappa, memS, memHeff, denS, nLower, VeffTilde);
         
         if (!atLeft){

            addDiagram1A(ikappa, memS, memHeff, denS, Xtensors[indexS-1]);
            
            addDiagram2b1and2b2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][0]);
            addDiagram2c1and2c2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][1]);
            addDiagram2b3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][0]);
            addDiagram2c3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][1]);
            addDiagram2b3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][0]);
            addDiagram2c3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][1]);
            
            addDiagram3Aand3D(ikappa, memS, memHeff, denS, Qtensors[indexS-1][0], Ltensors[indexS-1], temp);
            addDiagram3Band3I(ikappa, memS, memHeff, denS, Qtensors[indexS-1][1], Ltensors[indexS-1], temp);
            
            addDiagram4A1and4A2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1][1][0]);
            addDiagram4A1and4A2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1][1][0]);
            addDiagram4A3and4A4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][1][0]);
            addDiagram4A3and4A4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][1][0]);
            addDiagram4D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp);
            addDiagram4I(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp);
         
         }
         
         if (!atRight){
         
            addDiagram1B(ikappa, memS, memHeff, denS, Xtensors[indexS+1]);
            
            addDiagram2e1and2e2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][1]);
            addDiagram2f1and2f2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][0]);
            addDiagram2e3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][1]);
            addDiagram2f3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][0]);
            addDiagram2e3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][1]);
            addDiagram2f3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][0]);
            
            addDiagram3Kand3F(ikappa, memS, memHeff, denS, Qtensors[indexS+1][1], Ltensors[indexS+1], temp);
            addDiagram3Land3G(ikappa, memS, memHeff, denS, Qtensors[indexS+1][0], Ltensors[indexS+1], temp);
            
            addDiagram4F(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp);
            addDiagram4G(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp);
            addDiagram4J1and4J2spin0(ikappa, memS, memHeff, denS, Atensors[indexS+1][1][0]);
            addDiagram4J1and4J2spin1(ikappa, memS, memHeff, denS, Btensors[indexS+1][1][0]);
            addDiagram4J3and4J4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][1][0]);
            addDiagram4J3and4J4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][1][0]);
            
         }
         
         if ((!atLeft) && (!atRight)){
         
            addDiagram2a1spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp);
            addDiagram2a2spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp);
            addDiagram2a1spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp);
            addDiagram2a2spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp);
            addDiagram2a3spin0(ikappa, memS, memHeff, denS, Ctensors, F0tensors, temp);
            addDiagram2a3spin1(ikappa, memS, memHeff, denS, Dtensors, F1tensors, temp);
            
            addDiagram3C(ikappa, memS, memHeff, denS, Qtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram3J(ikappa, memS, memHeff, denS, Qtensors[indexS+1], Ltensors[indexS-1], temp);
            
            addDiagram4B1and4B2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B1and4B2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B3and4B4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4B3and4B4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C1and4C2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C1and4C2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C3and4C4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4C3and4C4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp);
            addDiagram4E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram4H(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram4K1and4K2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp);
            addDiagram4L1and4L2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp);
            addDiagram4K1and4K2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp);
            addDiagram4L1and4L2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp);
            addDiagram4K3and4K4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp);
            addDiagram4L3and4L4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp);
            addDiagram4K3and4K4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp);
            addDiagram4L3and4L4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp);
            
            addDiagram5A(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5B(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5C(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
            addDiagram5F(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2);
                  
         }
         
      }
      
      delete [] temp;
      delete [] temp2;
   
   }

}

void CheMPS2::Heff::fillHeffDiag(double * memHeffDiag, const Sobject * denS, TensorC **** Ctensors, TensorD **** Dtensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   
   /*for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
      const int ikappa = denS->gReorder(ikappaBIS);
      cout << "Size of block " << ikappa << " is " << denS->gKappa2index(ikappa+1) - denS->gKappa2index(ikappa) << endl;
   }*/
   
   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
   
      const int ikappa = denS->gReorder(ikappaBIS);

      for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeffDiag[cnt] = 0.0; }
      
      addDiagonal1C(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
      addDiagonal1D(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
      addDiagonal2d3all(ikappa, memHeffDiag, denS);
      if (nLower>0){ addDiagonalExcitations(ikappa, memHeffDiag, denS, nLower, VeffTilde); }
      
      if (!atLeft){
         addDiagonal1A(ikappa, memHeffDiag, denS, Xtensors[indexS-1]);
         addDiagonal2b3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][0]);
         addDiagonal2c3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][1]);
         addDiagonal2b3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][0]);
         addDiagonal2c3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][1]);
      }
      
      if (!atRight){
         addDiagonal1B(ikappa, memHeffDiag, denS, Xtensors[indexS+1]);
         addDiagonal2e3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][1]);
         addDiagonal2f3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][0]);
         addDiagonal2e3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][1]);
         addDiagonal2f3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][0]);
      }
      
      if ((!atLeft) && (!atRight)){
         addDiagonal2a3spin0(ikappa, memHeffDiag, denS, Ctensors, F0tensors);
         addDiagonal2a3spin1(ikappa, memHeffDiag, denS, Dtensors, F1tensors);
      }
      
   }
   
}

double CheMPS2::Heff::SolveDAVIDSON(Sobject * denS, TensorL *** Ltensors, TensorA **** Atensors, TensorB **** Btensors, TensorC **** Ctensors, TensorD **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   //From people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter11.pdf : algorithm 11.1, with instead of line (16), equation (11.3).
   const int DAVIDSON_NUM_VEC      = CheMPS2::HEFF_DAVIDSON_NUM_VEC;
   const int DAVIDSON_NUM_VEC_KEEP = CheMPS2::HEFF_DAVIDSON_NUM_VEC_KEEP;

   //Convert mem of Sobject to symmetric conventions
   denS->prog2symm();

   int length_vec = denS->gKappa2index(denS->gNKappa());
   int num_vec = 0;
   double ** vecs  = new double*[DAVIDSON_NUM_VEC];
   double ** Hvecs = new double*[DAVIDSON_NUM_VEC];
   int num_allocated = 0;
   
   double * mxM = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   double * mxM_eigs = new double[DAVIDSON_NUM_VEC];
   double * mxM_vecs = new double[DAVIDSON_NUM_VEC * DAVIDSON_NUM_VEC];
   int mxM_lwork = 3*DAVIDSON_NUM_VEC-1;
   double * mxM_work = new double[mxM_lwork];
   
   double rtol = CheMPS2::HEFF_DAVIDSON_RTOL_BASE * sqrt(length_vec);
   double rnorm = 10*rtol;
   
   double * t_vec = new double[length_vec];
   double * u_vec = new double[length_vec];
   double * work_vec = new double[length_vec];
   int inc1 = 1;
   dcopy_(&length_vec,denS->gStorage(),&inc1,t_vec,&inc1); //starting vector for Davidson is the current state of the Sobject in symmetric conventions.
   
   //Checking whether the S-object contains anything
   double Sobjectnorm = 0.0;
   for (int cnt=0; cnt<length_vec; cnt++){ Sobjectnorm += t_vec[cnt]*t_vec[cnt]; }
   Sobjectnorm = sqrt(Sobjectnorm);
   if (Sobjectnorm==0.0){
      for (int cnt=0; cnt<length_vec; cnt++){ t_vec[cnt] = ((double) rand())/RAND_MAX; }
      if (CheMPS2::HEFF_debugPrint){
         cout << "WARNING AT HEFF : S-object with zero norm was replaced with random vector. This should be a rare (but possible) warning." << endl;
      }
   }
   //End checking whether the S-object contains anything
   
   double * HeffDiag = new double[length_vec];
   fillHeffDiag(HeffDiag, denS, Ctensors, Dtensors, F0tensors, F1tensors, Xtensors, nLower, VeffTilde);
   
   double * Reortho_Lowdin = NULL;
   double * Reortho_Overlap_eigs = NULL;
   double * Reortho_Overlap = NULL;
   double * Reortho_Eigenvecs = NULL;
   bool Reortho_Allocated = false;
   
   int nIterations = 0;
   
   while (rnorm > rtol){

      //1. Orthogonalize the new t_vec w.r.t. the old basis
      for (int cnt=0; cnt<num_vec; cnt++){
         double min_overlap = - ddot_(&length_vec,t_vec,&inc1,vecs[cnt],&inc1);
         daxpy_(&length_vec,&min_overlap,vecs[cnt],&inc1,t_vec,&inc1);
      }
   
      //2. Normalize the t_vec
      char norm = 'F';
      double alpha = 1.0/dlange_(&norm,&length_vec,&inc1,t_vec,&length_vec,t_vec); //work not referenced as Frobenius norm
      dscal_(&length_vec,&alpha,t_vec,&inc1);
      
      //3. T_vec becomes part of vecs
      if (num_vec<num_allocated){
         double * temp = vecs[num_vec];
         vecs[num_vec] = t_vec;
         t_vec = temp;
      } else {
         vecs[num_allocated] = t_vec;
         Hvecs[num_allocated] = new double[length_vec];
         t_vec = new double[length_vec];
         num_allocated++;
      }
      makeHeff(vecs[num_vec], Hvecs[num_vec], denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
      nIterations++;
      
      //4. mxM contains the Hamiltonian in the basis "vecs"
      for (int cnt=0; cnt<num_vec; cnt++){
         mxM[cnt + DAVIDSON_NUM_VEC * num_vec] = ddot_(&length_vec,vecs[num_vec],&inc1,Hvecs[cnt],&inc1);
         mxM[num_vec + DAVIDSON_NUM_VEC * cnt] = mxM[cnt + DAVIDSON_NUM_VEC * num_vec];
      }
      mxM[num_vec + DAVIDSON_NUM_VEC * num_vec] = ddot_(&length_vec,vecs[num_vec],&inc1,Hvecs[num_vec],&inc1);
      
      //5. When t-vec was added to vecs, the number of vecs was actually increased by one. For convenience (doing 4.), only now the number is incremented.
      num_vec++;
      
      //6. Calculate the eigenvalues and vectors of mxM
      char jobz = 'V';
      char uplo = 'U';
      int info;
      for (int cnt1=0; cnt1<num_vec; cnt1++){
         for (int cnt2=0; cnt2<num_vec; cnt2++){
            mxM_vecs[cnt1 + DAVIDSON_NUM_VEC * cnt2] = mxM[cnt1 + DAVIDSON_NUM_VEC * cnt2];
         }
      }
      int lda = DAVIDSON_NUM_VEC;
      dsyev_(&jobz,&uplo,&num_vec,mxM_vecs,&lda,mxM_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
      
      //7. Calculate u and r. r is stored in t_vec, u in u_vec.
      for (int cnt=0; cnt<length_vec; cnt++){
         t_vec[cnt] = 0.0;
         u_vec[cnt] = 0.0;
      }
      for (int cnt=0; cnt<num_vec; cnt++){
         double alpha = mxM_vecs[cnt]; //eigenvector with lowest eigenvalue, hence mxM_vecs[cnt + DAVIDSON_NUM_VEC * 0]
         daxpy_(&length_vec,&alpha,Hvecs[cnt],&inc1,t_vec,&inc1);
         daxpy_(&length_vec,&alpha, vecs[cnt],&inc1,u_vec,&inc1);
      }
      alpha = -mxM_eigs[0];
      daxpy_(&length_vec,&alpha,u_vec,&inc1,t_vec,&inc1);
      
      //8. Calculate the norm of r
      rnorm = dlange_(&norm,&length_vec,&inc1,t_vec,&length_vec,t_vec);
      
      //9. In case convergence is not yet reached: prepare for the following iteration
      if (rnorm > rtol){
      
         //9a. Calculate the new t_vec based on the residual of the lowest eigenvalue, to add to the vecs.
         for (int cnt=0; cnt<length_vec; cnt++){
            if (fabs(HeffDiag[cnt] - mxM_eigs[0])> CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF ){
               work_vec[cnt] = u_vec[cnt]/(HeffDiag[cnt] - mxM_eigs[0]); // work_vec = K^(-1) u_vec
            } else {
               work_vec[cnt] = u_vec[cnt]/ CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF ;
               if (CheMPS2::HEFF_debugPrint) cout << "|(HeffDiag[" << cnt << "] - mxM_eigs[0])| = " << fabs(HeffDiag[cnt] - mxM_eigs[0]) << endl;
            }
         }
         alpha = - ddot_(&length_vec,work_vec,&inc1,t_vec,&inc1)/ddot_(&length_vec,work_vec,&inc1,u_vec,&inc1); // alpha = - (u^T K^(-1) r) / (u^T K^(-1) u)
         daxpy_(&length_vec,&alpha,u_vec,&inc1,t_vec,&inc1); // t_vec = r - (u^T K^(-1) r) / (u^T K^(-1) u) u
         for (int cnt=0; cnt<length_vec; cnt++){
            if (fabs(HeffDiag[cnt] - mxM_eigs[0])> CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF ){
               t_vec[cnt] = - t_vec[cnt]/(HeffDiag[cnt] - mxM_eigs[0]); //t_vec = - K^(-1) (r - (u^T K^(-1) r) / (u^T K^(-1) u) u)
            } else {
               t_vec[cnt] = - t_vec[cnt]/ CheMPS2::HEFF_DAVIDSON_PRECOND_CUTOFF ;
            }
         }
         
         // 9b. When the maximum number of vectors is reached: construct the one with lowest eigenvalue & restart
         if (num_vec == DAVIDSON_NUM_VEC){
         
            if (DAVIDSON_NUM_VEC_KEEP<=1){
            
               alpha = 1.0/dlange_(&norm,&length_vec,&inc1,u_vec,&length_vec,u_vec); //work not referenced as Frobenius norm
               dscal_(&length_vec,&alpha,u_vec,&inc1);
               dcopy_(&length_vec,u_vec,&inc1,vecs[0],&inc1);
               makeHeff(vecs[0], Hvecs[0], denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
               nIterations++;
               mxM[0] = ddot_(&length_vec,vecs[0],&inc1,Hvecs[0],&inc1);
            
               num_vec = 1;
            
            } else {
            
               //Construct the lowest DAVIDSON_NUM_VEC_KEEP eigenvectors
               if (!Reortho_Allocated) Reortho_Eigenvecs = new double[length_vec * DAVIDSON_NUM_VEC_KEEP];
               dcopy_(&length_vec,u_vec,&inc1,Reortho_Eigenvecs,&inc1);
               for (int cnt=1; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  for (int irow=0; irow<length_vec; irow++){
                     Reortho_Eigenvecs[irow + length_vec * cnt] = 0.0;
                     for (int ivec=0; ivec<DAVIDSON_NUM_VEC; ivec++){
                        Reortho_Eigenvecs[irow + length_vec * cnt] += vecs[ivec][irow] * mxM_vecs[ivec + DAVIDSON_NUM_VEC * cnt];
                     }
                  }
               }
               
               //Reorthonormalize them
               //Reortho: Calculate the overlap matrix
               if (!Reortho_Allocated) Reortho_Overlap = new double[DAVIDSON_NUM_VEC_KEEP * DAVIDSON_NUM_VEC_KEEP];
               char trans = 'T';
               char notr = 'N';
               int DVDS_KEEP = DAVIDSON_NUM_VEC_KEEP;
               double one = 1.0;
               double zero = 0.0; //set
               dgemm_(&trans,&notr,&DVDS_KEEP,&DVDS_KEEP,&length_vec,&one,Reortho_Eigenvecs,&length_vec, Reortho_Eigenvecs,&length_vec,&zero,Reortho_Overlap,&DVDS_KEEP);
               
               //Reortho: Calculate the Lowdin tfo
               if (!Reortho_Allocated) Reortho_Overlap_eigs = new double[DVDS_KEEP];
               dsyev_(&jobz,&uplo,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap_eigs,mxM_work,&mxM_lwork,&info); //ascending order of eigs
               for (int icnt=0; icnt<DVDS_KEEP; icnt++){
                  Reortho_Overlap_eigs[icnt] = pow(Reortho_Overlap_eigs[icnt],-0.25);
                  dscal_(&DVDS_KEEP, Reortho_Overlap_eigs+icnt, Reortho_Overlap+DVDS_KEEP*icnt, &inc1);
               }
               if (!Reortho_Allocated) Reortho_Lowdin = new double[DVDS_KEEP*DVDS_KEEP];
               dgemm_(&notr,&trans,&DVDS_KEEP,&DVDS_KEEP,&DVDS_KEEP,&one,Reortho_Overlap,&DVDS_KEEP,Reortho_Overlap,&DVDS_KEEP,&zero,Reortho_Lowdin,&DVDS_KEEP);
               
               //Reortho: Put the Lowdin tfo eigenvecs in vecs
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  dscal_(&length_vec,&zero,vecs[ivec],&inc1);
                  for (int ivec2=0; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     daxpy_(&length_vec,Reortho_Lowdin + ivec2 + DAVIDSON_NUM_VEC_KEEP * ivec,Reortho_Eigenvecs + length_vec * ivec2, &inc1, vecs[ivec], &inc1);
                  }
               }
               
               if (!Reortho_Allocated) Reortho_Allocated = true;
               
               //Construct the H*vecs
               for (int cnt=0; cnt<DAVIDSON_NUM_VEC_KEEP; cnt++){
                  makeHeff(vecs[cnt], Hvecs[cnt], denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
                  nIterations++;
               }
               
               //Build MxM
               for (int ivec=0; ivec<DAVIDSON_NUM_VEC_KEEP; ivec++){
                  for (int ivec2=ivec; ivec2<DAVIDSON_NUM_VEC_KEEP; ivec2++){
                     mxM[ivec + DAVIDSON_NUM_VEC * ivec2] = ddot_(&length_vec, vecs[ivec], &inc1, Hvecs[ivec2], &inc1);
                     mxM[ivec2 + DAVIDSON_NUM_VEC * ivec] = mxM[ivec + DAVIDSON_NUM_VEC * ivec2];
                  }
               }
               
               //Set num_vec
               num_vec = DAVIDSON_NUM_VEC_KEEP;
            
            }

         }
      }
      
   }
   
   if (CheMPS2::HEFF_debugPrint) cout << "   Stats: nIt(DAVIDSON) = " << nIterations << endl;
   
   double eigenvalue = mxM_eigs[0]; //mxM_eigs[0] and u_vec are eigenvalue and vector
   dcopy_(&length_vec,u_vec,&inc1,denS->gStorage(),&inc1);
   
   for (int cnt=0; cnt<num_allocated; cnt++){
      delete [] vecs[cnt];
      delete [] Hvecs[cnt];
   }
   delete [] vecs;
   delete [] Hvecs;
   delete [] t_vec;
   delete [] u_vec;
   delete [] work_vec;
   delete [] mxM;
   delete [] mxM_eigs;
   delete [] mxM_vecs;
   delete [] mxM_work;
   delete [] HeffDiag;
   
   if (Reortho_Allocated){
      delete [] Reortho_Eigenvecs;
      delete [] Reortho_Overlap;
      delete [] Reortho_Overlap_eigs;
      delete [] Reortho_Lowdin;
   }
   
   //convert denS->gStorage() to program conventions
   denS->symm2prog();
   
   return eigenvalue;

}

int CheMPS2::Heff::phase(const int TwoTimesPower){

   return (((TwoTimesPower/2)%2)!=0)?-1:1;

}



