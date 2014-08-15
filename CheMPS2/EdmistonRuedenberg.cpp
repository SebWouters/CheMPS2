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
#include <math.h>
#include <algorithm>

#include "EdmistonRuedenberg.h"
#include "DMRGSCFVmatRotations.h"
#include "Lapack.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::EdmistonRuedenberg::EdmistonRuedenberg(Hamiltonian * HamIn, const int printLevelIn){

   Ham = HamIn;
   printLevel = printLevelIn;
   SymmInfo.setGroup(Ham->getNGroup());
   
   int * Isizes = new int[SymmInfo.getNumberOfIrreps()];
   int * Zeroes = new int[SymmInfo.getNumberOfIrreps()];
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){ Isizes[irrep] = Zeroes[irrep] = 0; }
   for (int orb=0; orb<Ham->getL(); orb++){ Isizes[Ham->getOrbitalIrrep(orb)]++; }
   
   iHandler = new DMRGSCFindices(Ham->getL(), Ham->getNGroup(), Zeroes, Isizes, Zeroes); //Supposes all orbitals are active
   unitary  = new DMRGSCFunitary(iHandler);
   VmatRotated = new FourIndex(Ham->getNGroup(), Isizes);
   
   delete [] Zeroes;
   delete [] Isizes;

}

CheMPS2::EdmistonRuedenberg::~EdmistonRuedenberg(){

   delete unitary; //unitary needs iHandler in its destructor!
   delete VmatRotated;
   delete iHandler;
      
}

double CheMPS2::EdmistonRuedenberg::Optimize(double * temp1, double * temp2, const double gradThreshold, const int maxIter){

   DMRGSCFVmatRotations theRotator(Ham, iHandler);

   //Setting up the variables for the gradient
   double gradNorm = 1.0;
   int numVariables = 0;
   for (int irrep=0; irrep<SymmInfo.getNumberOfIrreps(); irrep++){
      if (iHandler->getNORB(irrep) > 1){
         numVariables += iHandler->getNORB(irrep) * (iHandler->getNORB(irrep) - 1) / 2;
      }
   }
   double * gradient = new double[numVariables];
   for (int cnt=0; cnt<numVariables; cnt++){ gradient[cnt] = 0.0; }
   
   //Setting up the variables for the cost function
   double Icost = 0.0;
   for (int orb=0; orb<Ham->getL(); orb++){ Icost += Ham->getVmat(orb,orb,orb,orb); }
   if (printLevel>0){ cout << "   EdmistonRuedenberg::Optimize : Cost function at start = " << Icost << endl; }
   double Icost_previous = 0.0;
   
   int nIterations = 0;
   
   while ((gradNorm > gradThreshold) && (nIterations < maxIter)){
   
      nIterations++;
   
      //Update the unitary
      unitary->updateUnitary(temp1, temp2, gradient, true, false); //multiply = true; compact = false
   
      //Rotate the Vmat
      Icost_previous = Icost;
      theRotator.fillVmatRotated(VmatRotated, unitary, temp1, temp2);
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
            theRotator.fillVmatRotated(VmatRotated, unitary, temp1, temp2);
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


