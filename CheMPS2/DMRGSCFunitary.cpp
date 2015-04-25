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
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>

#include "MyHDF5.h"
#include "Lapack.h"
#include "DMRGSCFunitary.h"

using std::string;
using std::ifstream;
using std::cout;
using std::endl;

CheMPS2::DMRGSCFunitary::DMRGSCFunitary(DMRGSCFindices * iHandlerIn){

   this->iHandler = iHandlerIn;
   
   //Allocate the unitary and set to I
   unitary = new double*[ iHandler->getNirreps() ];
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      const int linsize = iHandler->getNORB(irrep);
      const int size = linsize * linsize;
      unitary[irrep] = new double[size];
      for (int cnt=0; cnt<size; cnt++){ unitary[irrep][cnt] = 0.0; }
      for (int cnt=0; cnt<linsize; cnt++){ unitary[irrep][cnt*(1+linsize)] = 1.0; }
   }
   
   //Find the unique indices for OCC-ACT, OCC-VIRT, and ACT-VIRT rotations
   x_linearlength = 0;
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      x_linearlength += iHandler->getNOCC(irrep)*iHandler->getNDMRG(irrep) + iHandler->getNDMRG(irrep)*iHandler->getNVIRT(irrep) + iHandler->getNOCC(irrep)*iHandler->getNVIRT(irrep);
   }
   if (x_linearlength==0){ return; }
   x_firstindex  = new int[x_linearlength];
   x_secondindex = new int[x_linearlength];
   jumper = new int*[ iHandler->getNirreps() ];
   int jumper_previous = 0;
   int x_linearlength2 = 0;
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
      jumper[ irrep ] = new int[ 3 ];
      for (int geval=0; geval<3; geval++){
         if (geval==0){
            for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
               for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){ //DMRG is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNDMRGstart(irrep)+ cntDMRG;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNOCCstart(irrep) + cntOcc;
                  x_linearlength2++;
               }
            }
            jumper[irrep][geval] = jumper_previous;
            jumper_previous += iHandler->getNOCC(irrep)*iHandler->getNDMRG(irrep);
         }
         if (geval==1){
            for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
               for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNDMRGstart(irrep) + cntDMRG;
                  x_linearlength2++;
               }
            }
            jumper[irrep][geval] = jumper_previous;
            jumper_previous += iHandler->getNDMRG(irrep)*iHandler->getNVIRT(irrep);
         }
         if (geval==2){
            for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
               for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){ //Virt is the row index, hence fast moving one
                  x_firstindex[x_linearlength2]  = iHandler->getOrigNVIRTstart(irrep) + cntVirt;
                  x_secondindex[x_linearlength2] = iHandler->getOrigNOCCstart(irrep) + cntOcc;
                  x_linearlength2++;
               }
            }
            jumper[irrep][geval] = jumper_previous;
            jumper_previous += iHandler->getNOCC(irrep)*iHandler->getNVIRT(irrep);
         }
      }
   }
   assert( x_linearlength==x_linearlength2 );
   assert( x_linearlength==jumper_previous );

}

CheMPS2::DMRGSCFunitary::~DMRGSCFunitary(){
      
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){ delete [] unitary[irrep]; }
   delete [] unitary;
   
   if (x_linearlength!=0){
      delete [] x_firstindex;
      delete [] x_secondindex;
      for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){ delete [] jumper[irrep]; }
      delete [] jumper;
   }
      
}

int CheMPS2::DMRGSCFunitary::getNumVariablesX() const{ return x_linearlength; }

int CheMPS2::DMRGSCFunitary::getJumper( const int irrep, const int geval ) const{ return jumper[ irrep ][ geval ]; }

int CheMPS2::DMRGSCFunitary::getLinearIndex(const int p_index, const int q_index) const{

   int irrep_p = iHandler->getNirreps() - 1;
   int irrep_q = iHandler->getNirreps() - 1;
   
   while ( p_index < iHandler->getOrigNOCCstart( irrep_p ) ){ irrep_p--; }
   while ( q_index < iHandler->getOrigNOCCstart( irrep_q ) ){ irrep_q--; }
   
   assert( irrep_p == irrep_q );
   const int irrep = irrep_p;
   
   const bool p_virt = ( p_index >= iHandler->getOrigNVIRTstart( irrep )) ? true : false;
   const bool p_dmrg = ((p_index >= iHandler->getOrigNDMRGstart( irrep )) && (!( p_virt ))) ? true : false;
   
   const bool q_occ  = ( q_index < iHandler->getOrigNDMRGstart( irrep )) ? true : false;
   const bool q_dmrg = ((q_index < iHandler->getOrigNVIRTstart( irrep )) && (!( q_occ ))) ? true : false;
   
   // 0 : p DMRG q OCC
   // 1 : p VIRT q DMRG
   // 2 : p VIRT q OCC
   
   if ((p_dmrg) && (q_occ)){ // geval 0
      
      const int p_rel = p_index - iHandler->getOrigNDMRGstart( irrep );
      const int q_rel = q_index - iHandler->getOrigNOCCstart(  irrep );
      const int theLinIndex = getJumper( irrep, 0 ) + p_rel + iHandler->getNDMRG( irrep ) * q_rel;
      assert( p_index == x_firstindex[  theLinIndex ] );
      assert( q_index == x_secondindex[ theLinIndex ] );
      return theLinIndex;
      
   }
   
   if ((p_virt) && (q_dmrg)){ // geval 1
   
      const int p_rel = p_index - iHandler->getOrigNVIRTstart( irrep );
      const int q_rel = q_index - iHandler->getOrigNDMRGstart( irrep );
      const int theLinIndex = getJumper( irrep, 1 ) + p_rel + iHandler->getNVIRT( irrep ) * q_rel;
      assert( p_index == x_firstindex[  theLinIndex ] );
      assert( q_index == x_secondindex[ theLinIndex ] );
      return theLinIndex;
   
   }
   
   if ((p_virt) && (q_occ)){ // geval 2
   
      const int p_rel = p_index - iHandler->getOrigNVIRTstart( irrep );
      const int q_rel = q_index - iHandler->getOrigNOCCstart(  irrep );
      const int theLinIndex = getJumper( irrep, 2 ) + p_rel + iHandler->getNVIRT( irrep ) * q_rel;
      assert( p_index == x_firstindex[  theLinIndex ] );
      assert( q_index == x_secondindex[ theLinIndex ] );
      return theLinIndex;
   
   }
   
   assert( 0 == 1 );
   return -1;

}

int CheMPS2::DMRGSCFunitary::getFirstIndex(const int linearindex) const{ return x_firstindex[linearindex]; }

int CheMPS2::DMRGSCFunitary::getSecondIndex(const int linearindex) const{ return x_secondindex[linearindex]; }

double * CheMPS2::DMRGSCFunitary::getBlock(const int irrep){ return unitary[irrep]; }

void CheMPS2::DMRGSCFunitary::buildSkewSymmX(const int irrep, double * result, double * Xelem, const bool compact) const{

   const int linsize = iHandler->getNORB(irrep);
   for (int cnt=0; cnt<linsize*linsize; cnt++){ result[cnt] = 0.0; }

   if (compact){

      for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
         for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
            const int xsolindex = getJumper( irrep , 0 ) + cntDMRG + iHandler->getNDMRG( irrep ) * cntOcc;
            const int index1 = iHandler->getNOCC(irrep) + cntDMRG; //Index within irrep block
            result[ index1 + linsize*cntOcc ] =   Xelem[xsolindex];
            result[ cntOcc + linsize*index1 ] = - Xelem[xsolindex];
         }
      }
      for (int cntDMRG=0; cntDMRG<iHandler->getNDMRG(irrep); cntDMRG++){
         for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){
            const int xsolindex = getJumper( irrep , 1 ) + cntVirt + iHandler->getNVIRT( irrep ) * cntDMRG;
            const int index1 = iHandler->getNOCC(irrep) + iHandler->getNDMRG(irrep) + cntVirt; //Index within irrep block
            const int index2 = iHandler->getNOCC(irrep) + cntDMRG; //Index within irrep block
            result[ index1 + linsize*index2 ] =   Xelem[xsolindex];
            result[ index2 + linsize*index1 ] = - Xelem[xsolindex];
         }
      }
      for (int cntOcc=0; cntOcc<iHandler->getNOCC(irrep); cntOcc++){
         for (int cntVirt=0; cntVirt<iHandler->getNVIRT(irrep); cntVirt++){
            const int xsolindex = getJumper( irrep , 2 ) + cntVirt + iHandler->getNVIRT( irrep ) * cntOcc;
            const int index1 = iHandler->getNOCC(irrep) + iHandler->getNDMRG(irrep) + cntVirt; //Index within irrep block
            result[ index1 + linsize*cntOcc ] =   Xelem[xsolindex];
            result[ cntOcc + linsize*index1 ] = - Xelem[xsolindex];
         }
      }
         
   } else { //NOT compact
   
      int jump = 0;
      for (int cnt=0; cnt<irrep; cnt++){
         int linsizeCNT = iHandler->getNORB(cnt);
         jump += linsizeCNT * (linsizeCNT-1) / 2;
      }
      
      for (int row=0; row<linsize; row++){
         for (int col=row+1; col<linsize; col++){
            result[ row + linsize * col ] =   Xelem[ jump + row + col*(col-1)/2 ];
            result[ col + linsize * row ] = - Xelem[ jump + row + col*(col-1)/2 ];
         }
      }
   
   }

}

void CheMPS2::DMRGSCFunitary::updateUnitary(double * temp1, double * temp2, double * vector, const bool multiply, const bool compact){

   //Per irrep
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){

      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;
      
      if (linsize>1){ //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok
      
         double * xblock    = temp1;            //linsize*linsize
         double * Bmat      = temp1 +   size;   //linsize*linsize
         double * work1     = temp1 + 2*size;   //linsize*linsize
         double * work2     = temp1 + 3*size;   //linsize*linsize
         
         double * workLARGE = temp2;  //4*size
         int     lworkLARGE = 4*size; //4*size = 4*linsize*linsize > 3*linsize-1
         
         //Construct the antisymmetric x-matrix
         buildSkewSymmX(irrep, xblock, vector, compact);
         
         //Bmat <= xblock * xblock
         char notr = 'N';
         double alpha = 1.0;
         double beta = 0.0; //SET !!!
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,xblock,&linsize,&beta,Bmat,&linsize);
         
         //Bmat * work1 * Bmat^T <= xblock * xblock
         char uplo = 'U';
         char jobz = 'V';
         int info;
         dsyev_(&jobz, &uplo, &linsize, Bmat, &linsize, work1, workLARGE, &lworkLARGE, &info);
         
         //work2 <= Bmat^T * xblock * Bmat
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,Bmat,&linsize,&beta,work1,&linsize);
         char trans = 'T';
         dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work1,&linsize,&beta,work2,&linsize);
         
         if (CheMPS2::DMRGSCF_debugPrint){
            cout << "   DMRGSCFunitary::updateUnitary : Lambdas of irrep block " << irrep << " : " << endl;
            for (int cnt=0; cnt<linsize/2; cnt++){
               cout << "      block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
               cout << "              [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
            }
         }
         
         //work1 <= values of the antisymmetric 2x2 blocks
         for (int cnt=0; cnt<linsize/2; cnt++){ work1[cnt] = 0.5*( work2[2*cnt + linsize*(2*cnt+1)] - work2[2*cnt+1 + linsize*(2*cnt)] ); }
         
         if (CheMPS2::DMRGSCF_debugPrint){
            for (int cnt=0; cnt<linsize/2; cnt++){
               work2[2*cnt + linsize*(2*cnt+1)] -= work1[cnt];
               work2[2*cnt+1 + linsize*(2*cnt)] += work1[cnt];
            }
            double RMSdeviation = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSdeviation += work2[cnt] * work2[cnt]; }
            RMSdeviation = sqrt(RMSdeviation);
            cout << "   DMRGSCFunitary::updateUnitary : RMSdeviation of irrep block " << irrep << " (should be 0.0) = " << RMSdeviation << endl;
         }
         
         //work2 <= exp(Bmat^T * xblock * Bmat)
         for (int cnt=0; cnt<size; cnt++){ work2[cnt] = 0.0; }
         for (int cnt=0; cnt<linsize/2; cnt++){
            double cosine = cos(work1[cnt]);
            double sine   = sin(work1[cnt]);
            work2[2*cnt   + linsize*(2*cnt  )] = cosine;
            work2[2*cnt+1 + linsize*(2*cnt+1)] = cosine;
            work2[2*cnt   + linsize*(2*cnt+1)] = sine;
            work2[2*cnt+1 + linsize*(2*cnt  )] = -sine;
         }
         for (int cnt=2*(linsize/2); cnt<linsize; cnt++){ work2[cnt*(linsize + 1)] = 1.0; }
         
         //work2 <= Bmat * exp(Bmat^T * xblock * Bmat) * Bmat^T = exp(xblock)
         dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work2,&linsize,&beta,work1,&linsize);
         dgemm_(&notr,&trans,&linsize,&linsize,&linsize,&alpha,work1,&linsize,Bmat,&linsize,&beta,work2,&linsize);
         
         int inc = 1;
         if (multiply){ //U <-- exp(x) * U
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,work2,&linsize,unitary[irrep],&linsize,&beta,work1,&linsize);
            dcopy_(&size, work1, &inc, unitary[irrep], &inc);
         } else { //U <-- exp(x)
            dcopy_(&size, work2, &inc, unitary[irrep], &inc);
         }

      }
   }
   
   if (CheMPS2::DMRGSCF_debugPrint){ CheckDeviationFromUnitary(temp2); }

}

void CheMPS2::DMRGSCFunitary::rotateActiveSpaceVectors(double * eigenvecs, double * work){

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
   
   if (CheMPS2::DMRGSCF_debugPrint){ CheckDeviationFromUnitary(work); }

}

void CheMPS2::DMRGSCFunitary::CheckDeviationFromUnitary(double * work) const{

   char tran = 'T';
   char notr = 'N';
   double alpha = 1.0;
   double beta  = 0.0;
   
   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      if (linsize>0){
         dgemm_(&tran,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep],&linsize,unitary[irrep],&linsize,&beta,work,&linsize);
         double value = 0.0;
         for (int cnt=0; cnt<linsize; cnt++){
            value += (work[cnt*(1+linsize)]-1.0) * (work[cnt*(1+linsize)]-1.0);
            for (int cnt2=cnt+1; cnt2<linsize; cnt2++){
               value += work[cnt + cnt2*linsize] * work[cnt + cnt2*linsize] + work[cnt2 + cnt*linsize] * work[cnt2 + cnt*linsize];
            }
         }
         value = sqrt(value);
         cout << "   DMRGSCFunitary::CheckDeviationFromUnitary : 2-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] - I = " << value << endl;
      }
      
   }

}

void CheMPS2::DMRGSCFunitary::getLog(double * vector, double * temp1, double * temp2) const{

   int jump = 0;

   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;
      
      if (linsize>1){
         /* linsize >= 2; hence temp1 is at least of size 4*linsize*linsize
            if linsize <= 1; there corresponds no block in xmatrix to it */
         
         double * work1     = temp1;          //linsize * linsize
         double * work2     = temp1 +   size; //linsize * linsize
         double * work3     = temp1 + 2*size; //linsize * linsize
         double * work4     = temp1 + 3*size; //linsize * linsize
         
         double * workLARGE = temp2;
         int lworkLARGE     = 4*size;
         
         //S = U + U^T --> work1
         for (int row=0; row<linsize; row++){
            for (int col=0; col<linsize; col++){
               work1[row + linsize * col] = unitary[irrep][row + linsize * col] + unitary[irrep][col + linsize * row];
            }
         }
         
         //Get the eigenvectors of S = U + U^T: eigvals in work4, eigvecs in work1
         char jobz = 'V'; //compute eigenvectors
         char uplo = 'U';
         int info;
         dsyev_(&jobz, &uplo, &linsize, work1, &linsize, work4, workLARGE, &lworkLARGE, &info);
         
         //Transform the original orthogonal matrix with the real orthonormal eigenbasis of S --> the result V^T U V is stored in work3
         char trans = 'T';
         char notra = 'N';
         double alpha = 1.0;
         double beta = 0.0; //set
         dgemm_(&trans, &notra, &linsize, &linsize, &linsize, &alpha, work1, &linsize, unitary[irrep], &linsize, &beta, workLARGE, &linsize);
         dgemm_(&notra, &notra, &linsize, &linsize, &linsize, &alpha, workLARGE, &linsize, work1, &linsize, &beta, work3, &linsize);
         
         //Work3 should contain 1x1 blocks containing [+/-1] and 2x2 blocks containing [[cos(u) sin(u)][-sin(u) cos(u)]] and is hence TRIDIAGONAL.
         //If det(work3) = 1, a real-valued logarithm can be computed. If det(work3) = -1, the logarithm is complex valued.
         double f_low = 1.0;
         double f_high = work3[0];
         for (int cnt=1; cnt<linsize; cnt++){
            double temp = work3[cnt * (1 + linsize)] * f_high - work3[cnt + linsize*(cnt-1)] * work3[cnt-1 + linsize*cnt] * f_low;
            f_low = f_high;
            f_high = temp;
         }
         assert( f_high>0.0 );
         if (CheMPS2::DMRGSCF_debugPrint){
            cout << "   DMRGSCFunitary::getLog : Determinant of V^T*U*V for irrep " << irrep << " (should be +1) = " << f_high << endl;
         }
         
         //Fill work2 with ln(V^T U V) = ln(work3).
         for (int cnt=0; cnt<size; cnt++){ work2[cnt] = 0.0; }
         for (int cnt=0; cnt<linsize/2; cnt++){ //2x2 blocks
            const double cosinus = 0.5 * ( work3[2*cnt + linsize*2*cnt    ] + work3[2*cnt+1 + linsize*(2*cnt+1)] );
            const double sinus   = 0.5 * ( work3[2*cnt + linsize*(2*cnt+1)] - work3[2*cnt+1 + linsize*2*cnt    ] );
            const double theta   = atan2( sinus, cosinus );
            if (CheMPS2::DMRGSCF_debugPrint){
               cout << "   DMRGSCFunitary::getLog : block = [ " << work3[2*cnt   + linsize*2*cnt] << " , " << work3[2*cnt   + linsize*(2*cnt+1)] << " ]" << endl;
               cout << "                                    [ " << work3[2*cnt+1 + linsize*2*cnt] << " , " << work3[2*cnt+1 + linsize*(2*cnt+1)]
                                                                << " ] and corresponds to theta = " << theta << endl;
            }
            work3[2*cnt   + linsize*2*cnt    ] -= cosinus;
            work3[2*cnt+1 + linsize*(2*cnt+1)] -= cosinus;
            work3[2*cnt   + linsize*(2*cnt+1)] -= sinus;
            work3[2*cnt+1 + linsize*2*cnt    ] += sinus;
            work2[2*cnt   + linsize*(2*cnt+1)] = theta;
            work2[2*cnt+1 + linsize*2*cnt    ] = -theta;
         } //The rest are 1x1 blocks, corresponding to eigenvalue 1 --> ln(1) = 0 --> work2 does not need to be updated.
         for (int cnt=2*(linsize/2); cnt<linsize; cnt++){ work3[cnt*(linsize+1)] -= 1; }
         
         //Calculate the 2-norm of updated work3 (should be 0.0)
         if (CheMPS2::DMRGSCF_debugPrint){
            double RMSdeviation = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSdeviation += work3[cnt]*work3[cnt]; }
            RMSdeviation = sqrt(RMSdeviation);
            cout << "   DMRGSCFunitary::getLog : 2-norm of [ V^T*U*V - assumed blocks ] for irrep " << irrep << " (should be 0.0) = " << RMSdeviation << endl;
         }
         
         //Calculate V * ln(blocks) * V^T --> work4
         dgemm_(&notra, &notra, &linsize, &linsize, &linsize, &alpha, work1, &linsize, work2, &linsize, &beta, workLARGE, &linsize); //V*ln(blocks)
         dgemm_(&notra, &trans, &linsize, &linsize, &linsize, &alpha, workLARGE, &linsize, work1, &linsize, &beta, work4, &linsize); //V*ln(blocks)*V^T
         
         //Fill the vector with the just calculated ln(U) = work4
         for (int row=0; row<linsize; row++){
            for (int col=row+1; col<linsize; col++){
               vector[ jump + row + col*(col-1)/2 ] = 0.5 * ( work4[ row + linsize * col ] - work4[ col + linsize * row ] );
            }
         }
         
         jump += linsize*(linsize-1)/2;
         
      }
   
   }
   
   if (true){
      DMRGSCFunitary temporaryU = DMRGSCFunitary(iHandler);
      temporaryU.updateUnitary(temp1, temp2, vector, false, false); //multiply = compact = false
      double TwoNormDifference = 0.0;
      for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
         const int linsize = iHandler->getNORB(irrep);
         for (int cnt=0; cnt<linsize*linsize; cnt++){
            TwoNormDifference += ( temporaryU.getBlock(irrep)[cnt] - unitary[irrep][cnt] ) * ( temporaryU.getBlock(irrep)[cnt] - unitary[irrep][cnt] );
         }
      }
      TwoNormDifference = sqrt(TwoNormDifference);
      cout << "   DMRGSCFunitary::getLog : 2-norm of [ U - exp(ln(U)) ] (should be 0.0) = " << TwoNormDifference << endl;
   }

}

void CheMPS2::DMRGSCFunitary::BCH(double * Xprev, double * step, double * Xnew, double * temp1, double * temp2) const{

   int jump = 0;

   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;
      
      if (linsize>1){
         /* linsize >= 2; hence temp1 and temp2 are at least of size 4*linsize*linsize
            if linsize <= 1; the corresponding blocks in x are zero */
         
         double * theY      = temp1;          //linsize * linsize --> will contain Xprev in unfolded form
         double * xblock    = temp1 +   size; //linsize * linsize --> will contain step in unfolded form
         double * work1     = temp1 + 2*size; //linsize * linsize
         double * work2     = temp1 + 3*size; //linsize * linsize
         double * work3     = temp2;          //linsize * linsize
         double * work4     = temp2 +   size; //linsize * linsize
         //double * work5     = temp2 + 2*size; //linsize * linsize
         double * result    = temp2 + 3*size; //linsize * linsize
         
         //Copy Xprev to theY (right hand side of BCH)
         buildSkewSymmX(irrep, theY, Xprev, false); //NOT compact
         buildSkewSymmX(irrep, xblock, step, true); //compact
         
         //Z1 = X+Y
         int inc = 1;
         dcopy_(&size, xblock, &inc, result, &inc);
         double alpha = 1.0;
         daxpy_(&size, &alpha, theY, &inc, result, &inc);
         
         //Z2 = 0.5[X,Y]
         char notrans = 'N';
         double beta = 0.0; //SET
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, xblock, &linsize, theY, &linsize, &beta, work1, &linsize); //X*Y --> work1
         alpha = -1.0;
         beta = 1.0; //ADD
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, theY, &linsize, xblock, &linsize, &beta, work1, &linsize); //work1 = XY - YX = [X,Y]
         alpha = 0.5;
         daxpy_(&size, &alpha, work1, &inc, result, &inc);
         if (CheMPS2::DMRGSCF_debugPrint){
            double RMSmat = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSmat += work1[cnt] * work1[cnt]; }         
            RMSmat = sqrt(0.25 * RMSmat);
            cout << "   DMRGSCFunitary::BCH : for irrep " << irrep << endl;
            cout << "                         the 2-norm of the 2nd order term is " << RMSmat << endl;
         }
         
         //Z3 = 1/12 * ( [X,[X,Y]] - [Y,[X,Y]] )
         beta = 0.0; //SET
         alpha = 1.0;
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, xblock, &linsize, work1, &linsize, &beta, work2, &linsize); //X * [X,Y] --> work2
         alpha = -1.0;
         beta = 1.0; //ADD
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, work1, &linsize, xblock, &linsize, &beta, work2, &linsize); //work2 = [X,[X,Y]]
         beta = 0.0; //SET
         alpha = 1.0;
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, theY, &linsize, work1, &linsize, &beta, work3, &linsize); //Y * [X,Y] --> work3
         alpha = -1.0;
         beta = 1.0; //ADD
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, work1, &linsize, theY, &linsize, &beta, work3, &linsize); //work3 = [Y,[X,Y]]
         alpha = 1 / 12.0;
         daxpy_(&size, &alpha, work2, &inc, result, &inc);
         alpha = - 1 / 12.0;
         daxpy_(&size, &alpha, work3, &inc, result, &inc);
         if (CheMPS2::DMRGSCF_debugPrint){
            double RMSmat = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSmat += ( work2[cnt] - work3[cnt] ) * ( work2[cnt] - work3[cnt] ); }
            RMSmat = sqrt(RMSmat / 144.0);
            cout << "                         the 2-norm of the 3rd order term is " << RMSmat << endl;
         }
         
         //Z4 = -1/24 [Y,[X,[X,Y]]]
         beta = 0.0; //SET
         alpha = 1.0;
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, theY, &linsize, work2, &linsize, &beta, work4, &linsize); //Y * [X,[X,Y]] --> work4
         alpha = -1.0;
         beta = 1.0; //ADD
         dgemm_(&notrans, &notrans, &linsize, &linsize, &linsize, &alpha, work2, &linsize, theY, &linsize, &beta, work4, &linsize); //work4 = [Y,[X,[X,Y]]]
         alpha = - 1 / 24.0;
         daxpy_(&size, &alpha, work4, &inc, result, &inc);
         if (CheMPS2::DMRGSCF_debugPrint){
            double RMSmat = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSmat += work4[cnt] * work4[cnt]; }
            RMSmat = sqrt(RMSmat / 576.0);
            cout << "                         the 2-norm of the 4th order term is " << RMSmat << endl;
         }
         
         //Copy back the result
         for (int row=0; row<linsize; row++){
            for (int col=row+1; col<linsize; col++){
               Xnew[ jump + row + col*(col-1)/2 ] = 0.5 * ( result[ row + linsize * col ] - result[ col + linsize * row ] );
            }
         }
         
         if (CheMPS2::DMRGSCF_debugPrint){
            for (int row=0; row<linsize; row++){
               for (int col=row+1; col<linsize; col++){
                  result[ row + linsize * col ] -= Xnew[ jump + row + col*(col-1)/2 ];
                  result[ col + linsize * row ] += Xnew[ jump + row + col*(col-1)/2 ];
               }
            }
            double RMSmat = 0.0;
            for (int cnt=0; cnt<size; cnt++){ RMSmat += result[cnt] * result[cnt]; }
            RMSmat = sqrt(RMSmat);
            cout << "                         the 2-norm of the adapted result (should be 0.0) = " << RMSmat << endl;
         }
         
         jump += linsize*(linsize-1)/2;
         
      }
   }
   
   if (CheMPS2::DMRGSCF_debugPrint){
      DMRGSCFunitary temporaryU = DMRGSCFunitary(iHandler);
      temporaryU.updateUnitary(temp1, temp2, Xnew, false, false); //multiply = compact = false
      double TwoNormDifference = 0.0;
      for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
         const int linsize = iHandler->getNORB(irrep);
         for (int cnt=0; cnt<linsize*linsize; cnt++){
            TwoNormDifference += ( temporaryU.getBlock(irrep)[cnt] - unitary[irrep][cnt] ) * ( temporaryU.getBlock(irrep)[cnt] - unitary[irrep][cnt] );
         }
      }
      TwoNormDifference = sqrt(TwoNormDifference);
      cout << "   DMRGSCFunitary::getLog : 2-norm of [ U - exp( Xnew ) ] = " << TwoNormDifference << endl;
   }

}

void CheMPS2::DMRGSCFunitary::makeSureAllBlocksDetOne(double * temp1, double * temp2){

   for (int irrep=0; irrep<iHandler->getNirreps(); irrep++){
   
      int linsize = iHandler->getNORB(irrep);
      int size = linsize * linsize;
      
      if (linsize>1){
         /* linsize >= 2; hence temp1 is at least of size 4*linsize*linsize
            if linsize <= 1; there corresponds no block in xmatrix to it */
         
         double * work1     = temp1;          //linsize * linsize
         //double * work2     = temp1 +   size; //linsize * linsize
         double * work3     = temp1 + 2*size; //linsize * linsize
         double * work4     = temp1 + 3*size; //linsize * linsize
         
         double * workLARGE = temp2;
         int lworkLARGE     = 4*size;
         
         //S = U + U^T --> work1
         for (int row=0; row<linsize; row++){
            for (int col=0; col<linsize; col++){
               work1[row + linsize * col] = unitary[irrep][row + linsize * col] + unitary[irrep][col + linsize * row];
            }
         }
         
         //Get the eigenvectors of S = U + U^T: eigvals in work4, eigvecs in work1.
         char jobz = 'V'; //compute eigenvectors
         char uplo = 'U';
         int info;
         dsyev_(&jobz, &uplo, &linsize, work1, &linsize, work4, workLARGE, &lworkLARGE, &info);
         
         //Transform the original orthogonal matrix with the real orthonormal eigenbasis of S --> the result V^T U V is stored in work3.
         char trans = 'T';
         char notra = 'N';
         double alpha = 1.0;
         double beta = 0.0; //set
         dgemm_(&trans, &notra, &linsize, &linsize, &linsize, &alpha, work1, &linsize, unitary[irrep], &linsize, &beta, workLARGE, &linsize);
         dgemm_(&notra, &notra, &linsize, &linsize, &linsize, &alpha, workLARGE, &linsize, work1, &linsize, &beta, work3, &linsize);
         
         /*  Work3 should contain
                * 1x1 blocks containing [-1]
                * 2x2 blocks containing [[cos(u) sin(u)][-sin(u) cos(u)]]
                * 1x1 blocks containing [+1]
             and is hence TRIDIAGONAL. Its determinant is easy to compute: */
         double f_low = 1.0;
         double f_high = work3[0];
         for (int cnt=1; cnt<linsize; cnt++){
            double temp = work3[cnt * (1 + linsize)] * f_high - work3[cnt + linsize*(cnt-1)] * work3[cnt-1 + linsize*cnt] * f_low;
            f_low = f_high;
            f_high = temp;
         }
         if (f_high < 0.0){ //f_high = det = -1
            // U <-- diag(-1, 1, 1, 1, ...) * U : First row of U changes sign!
            for (int cnt=0; cnt<linsize; cnt++){ unitary[irrep][0 + linsize * cnt] *= -1; }
         }
         
      }
   
   }

}

void CheMPS2::DMRGSCFunitary::saveU(const string filename) const{

   hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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

void CheMPS2::DMRGSCFunitary::loadU(const string filename){

   hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
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

void CheMPS2::DMRGSCFunitary::deleteStoredUnitary(const string filename) const{

   std::stringstream temp;
   temp << "rm " << filename;
   int info = system(temp.str().c_str());
   cout << "Info on DMRGSCF::Unitary rm call to system: " << info << endl;

}



