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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "Heff.h"
#include "Davidson.h"
#include "Lapack.h"
#include "MPIchemps2.h"

CheMPS2::Heff::Heff(const SyBookkeeper * denBKIn, const Problem * ProbIn, const double dvdson_rtol_in){

   denBK = denBKIn;
   Prob = ProbIn;
   dvdson_rtol = dvdson_rtol_in;

}

CheMPS2::Heff::~Heff(){

}

void CheMPS2::Heff::makeHeff(double * memS, double * memHeff, const Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   const int DIM = std::max(denBK->gMaxDimAtBound(indexS), denBK->gMaxDimAtBound(indexS+2));
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   //PARALLEL
   #pragma omp parallel
   {
   
      double * temp  = new double[DIM*DIM];
      double * temp2 = new double[DIM*DIM];
   
      #pragma omp for schedule(dynamic)
      for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
      
         const int ikappa = denS->gReorder(ikappaBIS);
         for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeff[cnt] = 0.0; }
         
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_1cd2d3eh() == MPIRANK )
         #endif
         {
            addDiagram1C(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
            addDiagram1D(ikappa, memS,memHeff,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
            addDiagram2dall(ikappa, memS, memHeff, denS);
            addDiagram3Eand3H(ikappa, memS, memHeff, denS);
         }
         addDiagramExcitations(ikappa, memS, memHeff, denS, nLower, VeffTilde); //The MPI check occurs in this function
         
         if (!atLeft){

            /*********************
            *  Diagrams group 1  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_x() == MPIRANK )
            #endif
            {  addDiagram1A(ikappa, memS, memHeff, denS, Xtensors[indexS-1]); }

            /*********************
            *  Diagrams group 2  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS, indexS ) == MPIRANK )
            #endif
            {  addDiagram2b1and2b2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][0]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS+1, indexS+1 ) == MPIRANK )
            #endif
            { addDiagram2c1and2c2(ikappa, memS, memHeff, denS, Atensors[indexS-1][0][1]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS ) == MPIRANK )
            #endif
            {  addDiagram2b3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][0]);
               addDiagram2b3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][0]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS+1, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram2c3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][0][1]);
               addDiagram2c3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][0][1]); }

            /*********************
            *  Diagrams group 3  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_q( Prob->gL(), indexS ) == MPIRANK )
            #endif
            {  addDiagram3Aand3D(ikappa, memS, memHeff, denS, Qtensors[indexS-1][0], Ltensors[indexS-1], temp); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_q( Prob->gL(), indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram3Band3I(ikappa, memS, memHeff, denS, Qtensors[indexS-1][1], Ltensors[indexS-1], temp); }

            /*********************
            *  Diagrams group 4  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram4A1and4A2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1][1][0]);
               addDiagram4A1and4A2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1][1][0]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram4A3and4A4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1][1][0]);
               addDiagram4A3and4A4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1][1][0]); }
            addDiagram4D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp); //The MPI check occurs in this function
            addDiagram4I(ikappa, memS, memHeff, denS, Ltensors[indexS-1], temp); //The MPI check occurs in this function

         }
         
         if (!atRight){

            /*********************
            *  Diagrams group 1  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_x() == MPIRANK )
            #endif
            {  addDiagram1B(ikappa, memS, memHeff, denS, Xtensors[indexS+1]); }

            /*********************
            *  Diagrams group 2  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS, indexS ) == MPIRANK )
            #endif
            {  addDiagram2e1and2e2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][1]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS+1, indexS+1 ) == MPIRANK )
            #endif
            { addDiagram2f1and2f2(ikappa, memS, memHeff, denS, Atensors[indexS+1][0][0]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS ) == MPIRANK )
            #endif
            {  addDiagram2e3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][1]);
               addDiagram2e3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][1]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS+1, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram2f3spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][0][0]);
               addDiagram2f3spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][0][0]); }

            /*********************
            *  Diagrams group 3  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_q( Prob->gL(), indexS ) == MPIRANK )
            #endif
            {  addDiagram3Kand3F(ikappa, memS, memHeff, denS, Qtensors[indexS+1][1], Ltensors[indexS+1], temp); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_q( Prob->gL(), indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram3Land3G(ikappa, memS, memHeff, denS, Qtensors[indexS+1][0], Ltensors[indexS+1], temp); }

            /*********************
            *  Diagrams group 4  *
            *********************/
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_absigma( indexS, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram4J1and4J2spin0(ikappa, memS, memHeff, denS, Atensors[indexS+1][1][0]);
               addDiagram4J1and4J2spin1(ikappa, memS, memHeff, denS, Btensors[indexS+1][1][0]); }
            #ifdef CHEMPS2_MPI_COMPILATION
            if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS+1 ) == MPIRANK )
            #endif
            {  addDiagram4J3and4J4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS+1][1][0]);
               addDiagram4J3and4J4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS+1][1][0]); }
            addDiagram4F(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4G(ikappa, memS, memHeff, denS, Ltensors[indexS+1], temp); //The MPI check occurs in this function

         }
         
         if ((!atLeft) && (!atRight)){
         
            addDiagram2a1spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp); //The MPI check occurs in this function
            addDiagram2a2spin0(ikappa, memS, memHeff, denS, Atensors, S0tensors, temp); //The MPI check occurs in this function
            addDiagram2a1spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp); //The MPI check occurs in this function
            addDiagram2a2spin1(ikappa, memS, memHeff, denS, Btensors, S1tensors, temp); //The MPI check occurs in this function
            addDiagram2a3spin0(ikappa, memS, memHeff, denS, Ctensors, F0tensors, temp); //The MPI check occurs in this function
            addDiagram2a3spin1(ikappa, memS, memHeff, denS, Dtensors, F1tensors, temp); //The MPI check occurs in this function
            
            addDiagram3C(ikappa, memS, memHeff, denS, Qtensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram3J(ikappa, memS, memHeff, denS, Qtensors[indexS+1], Ltensors[indexS-1], temp); //The MPI check occurs in this function
            
            addDiagram4B1and4B2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4B1and4B2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4B3and4B4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4B3and4B4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4C1and4C2spin0(ikappa, memS, memHeff, denS, Atensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4C1and4C2spin1(ikappa, memS, memHeff, denS, Btensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4C3and4C4spin0(ikappa, memS, memHeff, denS, Ctensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4C3and4C4spin1(ikappa, memS, memHeff, denS, Dtensors[indexS-1], Ltensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram4H(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram4K1and4K2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4L1and4L2spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Atensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4K1and4K2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4L1and4L2spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Btensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4K3and4K4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4L3and4L4spin0(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ctensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4K3and4K4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp); //The MPI check occurs in this function
            addDiagram4L3and4L4spin1(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Dtensors[indexS+1], temp); //The MPI check occurs in this function
            
            addDiagram5A(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram5B(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram5C(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram5D(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram5E(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
            addDiagram5F(ikappa, memS, memHeff, denS, Ltensors[indexS-1], Ltensors[indexS+1], temp, temp2); //The MPI check occurs in this function
                  
         }
         
      }
      
      delete [] temp;
      delete [] temp2;
   
   }

}

void CheMPS2::Heff::fillHeffDiag(double * memHeffDiag, const Sobject * denS, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   const int indexS = denS->gIndex();
   const bool atLeft  = (indexS==0)?true:false;
   const bool atRight = (indexS==Prob->gL()-2)?true:false;
   #ifdef CHEMPS2_MPI_COMPILATION
   const int MPIRANK = MPIchemps2::mpi_rank();
   #endif
   
   //PARALLEL
   #pragma omp parallel for schedule(dynamic)
   for (int ikappaBIS=0; ikappaBIS<denS->gNKappa(); ikappaBIS++){
   
      const int ikappa = denS->gReorder(ikappaBIS);
      for (int cnt=denS->gKappa2index(ikappa); cnt<denS->gKappa2index(ikappa+1); cnt++){ memHeffDiag[cnt] = 0.0; }
      
      #ifdef CHEMPS2_MPI_COMPILATION
      if ( MPIchemps2::owner_1cd2d3eh() == MPIRANK )
      #endif
      {  addDiagonal1C(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS,indexS,indexS,indexS));
         addDiagonal1D(ikappa, memHeffDiag,denS,Prob->gMxElement(indexS+1,indexS+1,indexS+1,indexS+1));
         addDiagonal2d3all(ikappa, memHeffDiag, denS); }
      if (nLower>0){ addDiagonalExcitations(ikappa, memHeffDiag, denS, nLower, VeffTilde); } //The MPI check occurs in this function
      
      if (!atLeft){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_x() == MPIRANK )
         #endif
         {  addDiagonal1A(ikappa, memHeffDiag, denS, Xtensors[indexS-1]); }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS ) == MPIRANK )
         #endif
         {  addDiagonal2b3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][0]);
            addDiagonal2b3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][0]); }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf( Prob->gL(), indexS+1, indexS+1 ) == MPIRANK )
         #endif
         {  addDiagonal2c3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS-1][0][1]);
            addDiagonal2c3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS-1][0][1]); }
      }
      
      if (!atRight){
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_x() == MPIRANK )
         #endif
         {  addDiagonal1B(ikappa, memHeffDiag, denS, Xtensors[indexS+1]); }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf( Prob->gL(), indexS, indexS ) == MPIRANK )
         #endif
         {  addDiagonal2e3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][1]);
            addDiagonal2e3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][1]); }
         #ifdef CHEMPS2_MPI_COMPILATION
         if ( MPIchemps2::owner_cdf( Prob->gL(), indexS+1, indexS+1 ) == MPIRANK )
         #endif
         {  addDiagonal2f3spin0(ikappa, memHeffDiag, denS, Ctensors[indexS+1][0][0]);
            addDiagonal2f3spin1(ikappa, memHeffDiag, denS, Dtensors[indexS+1][0][0]); }
      }
      
      if ((!atLeft) && (!atRight)){
         addDiagonal2a3spin0(ikappa, memHeffDiag, denS, Ctensors, F0tensors); //The MPI check occurs in this function
         addDiagonal2a3spin1(ikappa, memHeffDiag, denS, Dtensors, F1tensors); //The MPI check occurs in this function
      }
      
   }
   
}

double CheMPS2::Heff::SolveDAVIDSON(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   #ifdef CHEMPS2_MPI_COMPILATION
   if ( MPIchemps2::mpi_rank() == MPI_CHEMPS2_MASTER ){
      return SolveDAVIDSON_main(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
   } else {
      return SolveDAVIDSON_help(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
   }
   #else
      return SolveDAVIDSON_main(denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
   #endif

}

double CheMPS2::Heff::SolveDAVIDSON_main(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   int inc1 = 1;
   int veclength = denS->gKappa2index( denS->gNKappa() );

   Davidson deBoskabouter( veclength, CheMPS2::DAVIDSON_NUM_VEC,
                                      CheMPS2::DAVIDSON_NUM_VEC_KEEP,
                                      // CheMPS2::DAVIDSON_DMRG_RTOL,
                                      dvdson_rtol,
                                      CheMPS2::DAVIDSON_PRECOND_CUTOFF, CheMPS2::HEFF_debugPrint );
   double ** whichpointers = new double*[2];

   char instruction = deBoskabouter.FetchInstruction( whichpointers );
   assert( instruction == 'A' );
   denS->prog2symm(); // Convert mem of Sobject to symmetric conventions
   dcopy_(&veclength, denS->gStorage(), &inc1, whichpointers[0], &inc1); // Starting vector for Davidson is the current state of the Sobject in symmetric conventions
   #ifdef CHEMPS2_MPI_COMPILATION
      double * workspace = new double[ veclength ];
      fillHeffDiag(workspace, denS, Ctensors, Dtensors, F0tensors, F1tensors, Xtensors, nLower, VeffTilde);
      MPIchemps2::reduce_array_double( workspace, whichpointers[1], veclength, MPI_CHEMPS2_MASTER );
   #else
      fillHeffDiag(whichpointers[1], denS, Ctensors, Dtensors, F0tensors, F1tensors, Xtensors, nLower, VeffTilde);
   #endif

   instruction = deBoskabouter.FetchInstruction( whichpointers );
   while ( instruction == 'B' ){
   
      #ifdef CHEMPS2_MPI_COMPILATION
      {
         int mpi_instruction = 2;
         MPIchemps2::broadcast_array_int( &mpi_instruction, 1, MPI_CHEMPS2_MASTER );
         MPIchemps2::broadcast_array_double( whichpointers[0], veclength, MPI_CHEMPS2_MASTER );
         makeHeff(whichpointers[0], workspace, denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
         MPIchemps2::reduce_array_double( workspace, whichpointers[1], veclength, MPI_CHEMPS2_MASTER );
      }
      #else
         makeHeff(whichpointers[0], whichpointers[1], denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
      #endif
      instruction = deBoskabouter.FetchInstruction( whichpointers );
   }

   assert( instruction == 'C' );
   dcopy_( &veclength, whichpointers[0], &inc1, denS->gStorage(), &inc1 ); // Copy the solution in symmetric conventions back
   denS->symm2prog(); // Convert mem of Sobject to program conventions
   double eigenvalue = whichpointers[1][0];
   if (CheMPS2::HEFF_debugPrint){ std::cout << "   Stats: nIt(DAVIDSON) = " << deBoskabouter.GetNumMultiplications() << std::endl; }
   delete [] whichpointers;
   #ifdef CHEMPS2_MPI_COMPILATION
      delete [] workspace;
      int mpi_instruction = 3;
      MPIchemps2::broadcast_array_int( &mpi_instruction, 1, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_double( &eigenvalue, 1, MPI_CHEMPS2_MASTER );
   #endif
   return eigenvalue;

}

#ifdef CHEMPS2_MPI_COMPILATION
double CheMPS2::Heff::SolveDAVIDSON_help(Sobject * denS, TensorL *** Ltensors, TensorOperator **** Atensors, TensorOperator **** Btensors, TensorOperator **** Ctensors, TensorOperator **** Dtensors, TensorS0 **** S0tensors, TensorS1 **** S1tensors, TensorF0 **** F0tensors, TensorF1 **** F1tensors, TensorQ *** Qtensors, TensorX ** Xtensors, int nLower, double ** VeffTilde) const{

   int veclength = denS->gKappa2index( denS->gNKappa() );
   double * vecin  = new double[ veclength ];
   double * vecout = new double[ veclength ];
   int mpi_instruction = -1;
   
   fillHeffDiag( vecout, denS, Ctensors, Dtensors, F0tensors, F1tensors, Xtensors, nLower, VeffTilde );
   MPIchemps2::reduce_array_double( vecout, vecin, veclength, MPI_CHEMPS2_MASTER );
   MPIchemps2::broadcast_array_int( &mpi_instruction, 1, MPI_CHEMPS2_MASTER );
   
   while ( mpi_instruction == 2 ){ // Mat Vec
   
      MPIchemps2::broadcast_array_double( vecin, veclength, MPI_CHEMPS2_MASTER );
      makeHeff(vecin, vecout, denS, Ltensors, Atensors, Btensors, Ctensors, Dtensors, S0tensors, S1tensors, F0tensors, F1tensors, Qtensors, Xtensors, nLower, VeffTilde);
      MPIchemps2::reduce_array_double( vecout, vecin, veclength, MPI_CHEMPS2_MASTER );
      MPIchemps2::broadcast_array_int( &mpi_instruction, 1, MPI_CHEMPS2_MASTER );
   
   }
   
   assert( mpi_instruction == 3 ); // Receive energy
   double eigenvalue = 0.0;
   MPIchemps2::broadcast_array_double( &eigenvalue, 1, MPI_CHEMPS2_MASTER );
   delete [] vecin;
   delete [] vecout;
   
   return eigenvalue; // The eigenvalue is correct on each process, denS not

}
#endif


