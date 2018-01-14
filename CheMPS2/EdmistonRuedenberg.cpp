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

#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <algorithm>

#include "EdmistonRuedenberg.h"
#include "DMRGSCFrotations.h"
#include "Lapack.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::EdmistonRuedenberg::EdmistonRuedenberg( const FourIndex * Vmat, const int group, const int printLevelIn ){

   VMAT_ORIG = Vmat;
   printLevel = printLevelIn;
   SymmInfo.setGroup( group );
   
   int * Isizes = new int[ SymmInfo.getNumberOfIrreps() ];
   int * Zeroes = new int[ SymmInfo.getNumberOfIrreps() ];
   int L = 0;
   for ( int irrep = 0; irrep < SymmInfo.getNumberOfIrreps(); irrep++ ){
      Isizes[ irrep ] = VMAT_ORIG->get_irrep_size( irrep );
      Zeroes[ irrep ] = 0;
      L += Isizes[ irrep ];
   }
   
   iHandler = new DMRGSCFindices( L, group, Zeroes, Isizes, Zeroes ); //Supposes all orbitals are active
   unitary  = new DMRGSCFunitary( iHandler );
   VmatRotated = new FourIndex( group, Isizes );
   
   delete [] Zeroes;
   delete [] Isizes;

}

CheMPS2::EdmistonRuedenberg::~EdmistonRuedenberg(){

   delete unitary; //unitary needs iHandler in its destructor!
   delete VmatRotated;
   delete iHandler;
      
}

double CheMPS2::EdmistonRuedenberg::Optimize(double * temp1, double * temp2, const bool startFromRandomUnitary, const double gradThreshold, const int maxIter){

   //Clear the unitary
   unitary->identity();

   //Setting up the variables for the gradient
   double gradNorm = 1.0;
   const int numVariables = iHandler->getROTparamsize();
   double * gradient = new double[numVariables];
   for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] = 0.0; }

   //Randomize the unitary if asked
   if (startFromRandomUnitary){
      for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] = ((double) rand())/RAND_MAX - 0.5; }
      unitary->updateUnitary(temp1, temp2, gradient, false, false); //multiply = compact = false
      for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] = 0.0; }
   }

   const int mem_size = iHandler->getL() * iHandler->getL() * iHandler->getL() * iHandler->getL();
   DMRGSCFrotations::rotate( VMAT_ORIG, VmatRotated, NULL, 'F', 'F', 'F', 'F', iHandler, unitary, temp1, temp2, mem_size, "edmistonruedenberg" );

   //Setting up the variables for the cost function
   double Icost = costFunction();
   if (printLevel>0){ cout << "   EdmistonRuedenberg::Optimize : Cost function at start = " << Icost << endl; }
   double Icost_previous = 0.0;
   
   int nIterations = 0;
   
   while ((gradNorm > gradThreshold) && (nIterations < maxIter)){
   
      nIterations++;
   
      //Update the unitary
      unitary->updateUnitary(temp1, temp2, gradient, true, false); //multiply = true; compact = false
   
      //Rotate the Vmat
      Icost_previous = Icost;
      DMRGSCFrotations::rotate( VMAT_ORIG, VmatRotated, NULL, 'F', 'F', 'F', 'F', iHandler, unitary, temp1, temp2, mem_size, "edmistonruedenberg" );
      Icost = costFunction();
      
      /* What if the cost function has dimished? Then make the rotation step a bit smaller!
         Rotate back ever closer to unitary_effective = I; untill Icost_previous < Icost
         1 - 0.5 - 0.5^2 - 0.5^3 - 0.5^4 - ... = 2 - ( 1 + 0.5 + 0.5^2 + ... ) = 2 - 1/(1-0.5) = 0  */
      if (Icost_previous > Icost){
         if (printLevel>1){ cout << "                                     WARNING : Icost = " << Icost << " < Icost_previous = " << Icost_previous << endl; }
         for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] = -gradient[cnt]; } //Switch the sign of the update : rotate back!
         int nIterationsBACK = 0;
         while ((Icost_previous > Icost) && (nIterationsBACK < CheMPS2::EDMISTONRUED_maxIterBackTfo)){
            nIterationsBACK++;
            for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] *= 0.5; }
            unitary->updateUnitary(temp1, temp2, gradient, true, false); //multiply = true; compact = false
            DMRGSCFrotations::rotate( VMAT_ORIG, VmatRotated, NULL, 'F', 'F', 'F', 'F', iHandler, unitary, temp1, temp2, mem_size, "edmistonruedenberg" );
            Icost = costFunction();
         }
         if (printLevel>1){ cout << "                                     WARNING : Rotated back a bit. Now Icost = " << Icost << endl; }
      }
      
      //Calculate the gradient and Hessian and update
      gradNorm = augmentedHessianNewtonRaphson(gradient, temp1, temp2);
      
      if (printLevel>1){
         cout << "                                  After iteration " << nIterations << " : Icost = " << Icost << " has gradNorm = " << gradNorm << endl; 
      }
   
   }
   
   delete [] gradient;
   
   if (printLevel>0){
      cout << "                                  Cost function at stop  = " << Icost << endl;
      cout << "                                  Gradient norm = " << gradNorm << " after " << nIterations << " iterations." << endl;
   }
   
   return Icost;

}

CheMPS2::DMRGSCFunitary * CheMPS2::EdmistonRuedenberg::getUnitary(){ return unitary; }

double CheMPS2::EdmistonRuedenberg::costFunction() const{

   double Cost = 0.0;
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      for (int orb=0; orb<iHandler->getNORB(irrep); orb++){
         Cost += VmatRotated->get(irrep,irrep,irrep,irrep,orb,orb,orb,orb);
      }
   }
   return Cost;

}

double CheMPS2::EdmistonRuedenberg::calcGradientValue(const int irrep, const int p, const int q) const{

   return 4 * ( VmatRotated->get(irrep, irrep, irrep, irrep, p, p, p, q) - VmatRotated->get(irrep, irrep, irrep, irrep, q, q, q, p) );

}

double CheMPS2::EdmistonRuedenberg::calcHessianValue(const int irrep, const int p, const int q, const int r, const int s) const{

   double value = 0.0;
   if (p==r){
      value += ( 8*VmatRotated->get(irrep, irrep, irrep, irrep, p, p, q, s)
               + 4*VmatRotated->get(irrep, irrep, irrep, irrep, p, q, p, s)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, q, q, q, s)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, s, s, s, q) );
   }
   if (q==s){
      value += ( 8*VmatRotated->get(irrep, irrep, irrep, irrep, q, q, p, r)
               + 4*VmatRotated->get(irrep, irrep, irrep, irrep, q, p, q, r)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, p, p, p, r)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, r, r, r, p) );
   }
   if (p==s){
      value -= ( 8*VmatRotated->get(irrep, irrep, irrep, irrep, p, p, q, r)
               + 4*VmatRotated->get(irrep, irrep, irrep, irrep, p, q, p, r)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, q, q, q, r)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, r, r, r, q) );
   }
   if (q==r){
      value -= ( 8*VmatRotated->get(irrep, irrep, irrep, irrep, q, q, p, s)
               + 4*VmatRotated->get(irrep, irrep, irrep, irrep, q, p, q, s)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, p, p, p, s)
               - 2*VmatRotated->get(irrep, irrep, irrep, irrep, s, s, s, p) );
   }
   return value;

}

double CheMPS2::EdmistonRuedenberg::augmentedHessianNewtonRaphson(double * gradient, double * temp1, double * temp2) const{

   int jump = 0;
   double gradNorm = 0.0;
   
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      
      int linsize = iHandler->getNORB(irrep);
      
      if (linsize>1){
      
         int hessianlinearsize = linsize*(linsize-1)/2 + 1;
         
         /*   linsize    linsize^4    hesslinsize    hesslinsize^2   3*hesslinsize-1   4*hesslinsize-1   (hesslinsize+1)*hesslinsize
              2          16           2              4               5                 7                 6
              3          81           4              16              11                15                20
              4          256          7              49              20                27                56
              5          625          11             121             32                43                132
              6          1296         16             256             47                63                272
              
              Conclusion : - temp1 is large enough for hessian (hessianlinearsize^2)
                           - temp2 is large enough for workspace (max(hesslinsize^2, 3*hesslinsize-1)) AND eigenvecs (hesslinsize)
         */
         
         double * hessian = temp1;
         double * eigen   = temp2;
         double * work    = temp2 + hessianlinearsize;
         int lwork        = linsize*linsize*linsize*linsize - hessianlinearsize;
      
         /*   Maximizing the cost function  O(x) =  O(0) + g^T * x + x^T * H * x / 2
         ==   Minimizing the cost function -O(x) = -O(0) - g^T * x - x^T * H * x / 2 */
         for (int row=0; row<linsize; row++){
            for (int col=row+1; col<linsize; col++){
               const double waarde = calcGradientValue(irrep, row, col);
               gradient[jump + row + col*(col-1)/2] = waarde;
               gradNorm += waarde * waarde;
               for (int row2=0; row2<linsize; row2++){
                  for (int col2=row2+1; col2<linsize; col2++){
                     const double value = calcHessianValue(irrep, row, col, row2, col2);
                     hessian[row + col*(col-1)/2 + hessianlinearsize * (row2 + col2*(col2-1)/2)] = - value; //MINUS SIGN FOR OTHER DIRECTION
                  }
               }
            }
         }
         
         for (int cnt=0; cnt<hessianlinearsize-1; cnt++){
            hessian[hessianlinearsize-1 + hessianlinearsize * cnt]   = - gradient[jump + cnt]; //MINUS SIGN FOR OTHER DIRECTION
            hessian[cnt + hessianlinearsize * (hessianlinearsize-1)] = - gradient[jump + cnt]; //MINUS SIGN FOR OTHER DIRECTION
         }
         hessian[hessianlinearsize*hessianlinearsize-1] = 0.0;
   
         //Find its lowest eigenvalue and vector
         char jobz = 'V';
         char uplo = 'U';
         int info;
         dsyev_(&jobz,&uplo,&hessianlinearsize,hessian,&hessianlinearsize,eigen,work,&lwork,&info);
   
         double scalar = 1.0/hessian[hessianlinearsize-1];
         int inc = 1;
         dscal_(&hessianlinearsize,&scalar,hessian,&inc);
         
         for (int cnt=0; cnt<hessianlinearsize-1; cnt++){ gradient[jump + cnt] = hessian[cnt]; }
         
         jump += hessianlinearsize-1;
      
      }
   
   }
   
   gradNorm = sqrt(gradNorm);
   
   return gradNorm;

}

void CheMPS2::EdmistonRuedenberg::Fiedler(const int irrep, int * reorder, double * laplacian, double * temp2){

   //For information on the Fiedler vector: see http://web.eecs.utk.edu/~mberry/order/node9.html

   int linsize = iHandler->getNORB(irrep);
   assert( linsize>=2 );
   
   //Preamble: linsize>=2 at this point
   double * work = temp2;             //temp2 at least of size 4*linsize*linsize
   int lwork     = 3*linsize*linsize; //3*linsize*linsize > 3*linsize-1 for linsize>=2
   double * eigs = temp2 + lwork;     //has size linsize*linsize > linsize for linsize>=2
   
   //Calculate the eigenspectrum of the Laplacian
   char jobz = 'V';
   char uplo = 'U';
   int info;
   dsyev_(&jobz, &uplo, &linsize, laplacian, &linsize, eigs, work, &lwork, &info);
   if (printLevel>1){
      cout << "   EdmistonRuedenberg::Fiedler : Smallest eigs(Laplacian[" << irrep << "]) = [ " << eigs[0] << "  ,  " << eigs[1] << " ]." << endl;
   }
   
   //The second eigenvector is the Fiedler vector
   double * FiedlerVector     = laplacian + linsize;
   double * FiedlerVectorCopy = laplacian;
   for (int cnt=0; cnt<linsize; cnt++){ FiedlerVectorCopy[cnt] = FiedlerVector[cnt]; }
   
   //Find the Fiedler ordering
   const double upperBound = 2.0; //Sum_i FiedlerVector[i]^2 = 1 --> fabs(FiedlerVector[i]) <= 1.0
   for (int value=0; value<linsize; value++){
      int index=0;
      for (int cnt=1; cnt<linsize; cnt++){
         if (FiedlerVectorCopy[cnt] < FiedlerVectorCopy[index]){ index = cnt; }
      } //Index now corresponds to the smallest value in FiedlerVectorCopy
      reorder[value] = index;
      FiedlerVectorCopy[index] = upperBound; //Remove from the interesting entries in FiedlerVectorCopy
   }
   
   if (printLevel>1){
      bool isOK = true;
      for (int cnt=0; cnt<linsize-1; cnt++){
         if (FiedlerVector[reorder[cnt]] > FiedlerVector[reorder[cnt+1]]){ isOK = false; }
      }
      assert( isOK );
      cout << "                                 Reorder[" << irrep << "] = [ ";
      for (int cnt=0; cnt<linsize-1; cnt++){ cout << reorder[cnt] << "  ,  "; }
      cout << reorder[linsize-1] << " ]." << endl;
   }
   
   //Reorder the vectors in the unitary
   double * blockU = unitary->getBlock(irrep);
   for (int row=0; row<linsize; row++){
      for (int col=0; col<linsize; col++){
         work[row + linsize * col] = blockU[reorder[row] + linsize * col];
      }
   }
   
   int size = linsize*linsize;
   int inc = 1;
   dcopy_(&size, work, &inc, blockU, &inc);
   
   if (printLevel>1){
      char trans = 'T';
      char notrans = 'N';
      double alpha = 1.0;
      double beta = 0.0;
      dgemm_(&trans, &notrans, &linsize, &linsize, &linsize, &alpha, blockU, &linsize, blockU, &linsize, &beta, work, &linsize);
      double sum = 0.0;
      for (int row=0; row<linsize; row++){
         sum += (work[row+linsize*row]-1.0)*(work[row+linsize*row]-1.0);
         for (int col=row+1; col<linsize; col++){
            sum += work[row+linsize*col]*work[row+linsize*col] + work[col+linsize*row]*work[col+linsize*row];
         }
      }
      sum = sqrt(sum);
      cout << "                                 2-norm of Unitary[" << irrep << "]^T * Unitary[" << irrep << "] - I = " << sum << endl;
   }

}

double CheMPS2::EdmistonRuedenberg::FiedlerExchangeCost() const{

   double Cost = 0.0;
   
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      const int linsize = iHandler->getNORB(irrep);
      if (linsize>1){
         for (int row=0; row<linsize; row++){
            for (int col=row+1; col<linsize; col++){
               Cost += 2 * VmatRotated->get(irrep,irrep,irrep,irrep,row,col,col,row) * (col-row) * (col-row);
            }
         }
      }
   }
   
   return Cost;

}

void CheMPS2::EdmistonRuedenberg::FiedlerExchange(const int maxlinsize, double * temp1, double * temp2){

   //For information on the Fiedler vector: see http://web.eecs.utk.edu/~mberry/order/node9.html

   const int mem_size = iHandler->getL() * iHandler->getL() * iHandler->getL() * iHandler->getL();
   DMRGSCFrotations::rotate( VMAT_ORIG, VmatRotated, NULL, 'F', 'F', 'F', 'F', iHandler, unitary, temp1, temp2, mem_size, "edmistonruedenberg" );

   if (printLevel>0){ cout << "   EdmistonRuedenberg::FiedlerExchange : Cost function at start = " << FiedlerExchangeCost() << endl; }
   
   int * reorder = new int[maxlinsize];

   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      
      const int linsize = iHandler->getNORB(irrep);
      if (linsize>1){
      
         //temp1 and temp2 both at least of size 4*linsize*linsize
         double * laplacian = temp1;
         
         //Fill the weighted graph Laplacian
         for (int row=0; row<linsize; row++){
            laplacian[row*(1+linsize)] = 0.0;
            for (int col=row+1; col<linsize; col++){
               laplacian[row + linsize*col]  = - VmatRotated->get(irrep,irrep,irrep,irrep,row,col,col,row); //minus the exchange matrix
               laplacian[col + linsize*row]  = laplacian[row + linsize*col]; //Symmetric matrix
               laplacian[row + linsize*row] -= laplacian[row + linsize*col]; //On diagonal : Minus the sum of the other terms on that row
            }
            for (int col=0; col<row; col++){
               laplacian[row + linsize*row] -= laplacian[row + linsize*col]; //On diagonal : Minus the sum of the other terms on that row
            }
         }
         
         Fiedler(irrep, reorder, laplacian, temp2);
         
      }
   }
   
   delete [] reorder;

   DMRGSCFrotations::rotate( VMAT_ORIG, VmatRotated, NULL, 'F', 'F', 'F', 'F', iHandler, unitary, temp1, temp2, mem_size, "edmistonruedenberg" );

   if (printLevel>0){ cout << "   EdmistonRuedenberg::FiedlerExchange : Cost function at end   = " << FiedlerExchangeCost() << endl; }

}

double CheMPS2::EdmistonRuedenberg::FiedlerGlobalCost( const DMRGSCFindices * idx, const FourIndex * VMAT_LOCAL, int * dmrg2ham ){

   double cost = 0.0;

   for ( int dmrg_row = 0; dmrg_row < idx->getL(); dmrg_row++ ){
      for ( int dmrg_col = 0; dmrg_col < idx->getL(); dmrg_col++ ){
         const int ham_row = dmrg2ham[ dmrg_row ];
         const int ham_col = dmrg2ham[ dmrg_col ];
         const int irrep_row = idx->getOrbitalIrrep( ham_row );
         const int irrep_col = idx->getOrbitalIrrep( ham_col );
         const int rel_row = ham_row - idx->getOrigNOCCstart( irrep_row );
         const int rel_col = ham_col - idx->getOrigNOCCstart( irrep_col );
         cost += VMAT_LOCAL->get( irrep_row, irrep_col, irrep_col, irrep_row, rel_row, rel_col, rel_col, rel_row ) * ( dmrg_row - dmrg_col ) * ( dmrg_row - dmrg_col );
      }
   }

   return cost;

}

void CheMPS2::EdmistonRuedenberg::FiedlerGlobal( int * dmrg2ham ) const{

   // For information on the Fiedler vector: see http://web.eecs.utk.edu/~mberry/order/node9.html

   for ( int orb = 0; orb < iHandler->getL(); orb++ ){ dmrg2ham[ orb ] = orb; }
   if ( printLevel > 0 ){ cout << "   EdmistonRuedenberg::FiedlerGlobal : Cost function at start = " << FiedlerGlobalCost( iHandler, VMAT_ORIG, dmrg2ham ) << endl; }

   // Build the Laplacian
   double * laplacian = new double[ iHandler->getL() * iHandler->getL() ];
   for ( int ham_row = 0; ham_row < iHandler->getL(); ham_row++ ){
      double sum_over_column = 0.0;
      for ( int ham_col = 0; ham_col < iHandler->getL(); ham_col++ ){
         if ( ham_row != ham_col ){
            const int irrep_row = iHandler->getOrbitalIrrep( ham_row );
            const int irrep_col = iHandler->getOrbitalIrrep( ham_col );
            const int rel_row = ham_row - iHandler->getOrigNOCCstart( irrep_row );
            const int rel_col = ham_col - iHandler->getOrigNOCCstart( irrep_col );
            const double value = fabs( VMAT_ORIG->get( irrep_row, irrep_col, irrep_col, irrep_row, rel_row, rel_col, rel_col, rel_row ) );
            laplacian[ ham_row + iHandler->getL() * ham_col ] = - value;
            sum_over_column += value;
         } else {
            laplacian[ ham_row + iHandler->getL() * ham_col ] = 0.0;
         }
      }
      laplacian[ ham_row + iHandler->getL() * ham_row ] = sum_over_column;
   }

   // Calculate the eigenspectrum of the Laplacian
   int lwork     = 3 * iHandler->getL() * iHandler->getL();
   double * work = new double[ lwork ];
   double * eigs = new double[ iHandler->getL() ];
   char jobz     = 'V';
   char uplo     = 'U';
   int linsize   = iHandler->getL();
   int info;
   dsyev_( &jobz, &uplo, &linsize, laplacian, &linsize, eigs, work, &lwork, &info );
   delete [] work;
   delete [] eigs;

   // Fill dmrg2ham
   double * fiedler_vec = laplacian + linsize;
   for ( int dummy = 0; dummy < linsize; dummy++ ){
      int index = 0;
      for ( int orb = 1; orb < linsize; orb++ ){
         if ( fiedler_vec[ orb ] < fiedler_vec[ index ] ){ index = orb; }
      }
      dmrg2ham[ dummy ] = index;
      fiedler_vec[ index ] = 2.0; // Eigenvectors are normalized to 1.0, so certainly OK
   }

   delete [] laplacian;

   if ( printLevel > 0 ){
      cout << "   EdmistonRuedenberg::FiedlerGlobal : Cost function at end   = " << FiedlerGlobalCost( iHandler, VMAT_ORIG, dmrg2ham ) << endl;
      cout << "   EdmistonRuedenberg::FiedlerGlobal : Reordering = [ ";
      for ( int orb = 0; orb < linsize - 1; orb++ ){ cout << dmrg2ham[ orb ] << ", "; }
      cout << dmrg2ham[ linsize - 1 ] << " ]." << endl;
   }

}


