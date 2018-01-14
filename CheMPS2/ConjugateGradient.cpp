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

#include "ConjugateGradient.h"

using std::cout;
using std::endl;

CheMPS2::ConjugateGradient::ConjugateGradient( const int veclength_in, const double RTOL_in, const double DIAG_CUTOFF_in, const bool print_in ){

   veclength = veclength_in;
   RTOL = RTOL_in;
   DIAG_CUTOFF = DIAG_CUTOFF_in;
   print = print_in;

   state = 'I';
   num_matvec = 0;

   XVEC   = new double[ veclength ];
   PRECON = new double[ veclength ];
   RHS    = new double[ veclength ];
   WORK   = new double[ veclength ];
   RESID  = new double[ veclength ];
   PVEC   = new double[ veclength ];
   OPVEC  = new double[ veclength ];

}

CheMPS2::ConjugateGradient::~ConjugateGradient(){

   delete [] XVEC;
   delete [] PRECON;
   delete [] RHS;
   delete [] WORK;
   delete [] RESID;
   delete [] PVEC;
   delete [] OPVEC;

}

int CheMPS2::ConjugateGradient::get_num_matvec() const{ return num_matvec; }

char CheMPS2::ConjugateGradient::step( double ** pointers ){

   /*
      Possible states:
       - I : just created the class
       - G : the guess has been set in XVEC, the diagonal in PRECON, and the right-hand side in RHS
       - H : the PRECON, RESID, and XVEC have been set to start the iterations of ( PRECON * operator * PRECON ) * XVEC = RESID = PRECON * RHS
       - J : at start-up OPVEC contains operator * PRECON * XVEC
       - K : OPVEC, RESID, PVEC have just been set, as well as rnorm and rkT_rk
       - L : OPVEC contains operator * PRECON * x_k
       - Y : XVEC contains x = operator^{-1} * rhs, and OPVEC contains operator * XVEC
       - Z : the converged signal has been given to the user, nothing remains to be done
      
      Possible instructions:
       - A : copy the guess to pointers[0], the diagonal of the operator to pointers[1], and the right-hand side of the problem to pointers[2]
       - B : perform pointers[1] = operator * pointers[0]
       - C : pointers[0] contains the solution; pointers[1][0] the residual norm
       - D : there was an error
   */

   if ( state == 'I' ){
      pointers[0] = XVEC;
      pointers[1] = PRECON;
      pointers[2] = RHS;
      state = 'G';
      return 'A';
   }

   if ( state == 'G' ){
      stepG2H();
      state = 'H';
   }

   if ( state == 'H' ){
      apply_precon( XVEC, WORK );
      pointers[0] = WORK;
      pointers[1] = OPVEC;
      state = 'J';
      num_matvec++;
      return 'B';
   }

   if ( state == 'J' ){
      stepJ2K();
      state = 'K';
   }

   if ( state == 'L' ){
      stepL2K();
      state = 'K';
   }

   if ( state == 'K' ){
      if ( rnorm >= RTOL ){
         apply_precon( PVEC, WORK ); // WORK = PRECON * PVEC
         pointers[0] = WORK;
         pointers[1] = OPVEC;
         state = 'L';
      } else {
         apply_precon( XVEC );
         pointers[0] = XVEC;
         pointers[1] = OPVEC;
         state = 'Y';
      }
      num_matvec++;
      return 'B';
   }

   if ( state == 'Y' ){
      stepY2Z();
      pointers[0] = XVEC;
      pointers[1] = WORK;
      pointers[1][0] = rnorm;
      state = 'Z';
      return 'C';
   }

   return 'D';

}

void CheMPS2::ConjugateGradient::stepL2K(){

   apply_precon( OPVEC );                                    // OPVEC_old = ( PRECON * operator * PRECON ) * PVEC_old
   const double alpha = rdotr / inprod( PVEC, OPVEC );       // alpha = RESID_old^T * RESID_old / ( PVEC_old^T * ( PRECON * operator * PRECON ) * PVEC_old )
   for ( int elem = 0; elem < veclength; elem++ ){
      XVEC[ elem ] = XVEC[ elem ] + alpha * PVEC[ elem ];    // XVEC_new <-- XVEC_old + alpha * PVEC_old
   }
   for ( int elem = 0; elem < veclength; elem++ ){
      RESID[ elem ] = RESID[ elem ] - alpha * OPVEC[ elem ]; // RESID_new <-- RESID_old - alpha * ( PRECON * operator * PRECON ) * PVEC_old
   }
   const double new_rdotr = inprod( RESID );
   const double beta = new_rdotr / rdotr;                    // beta = RESID_new^T * RESID_new / ( RESID_old^T * RESID_old )
   for ( int elem = 0; elem < veclength; elem++ ){
      PVEC[ elem ] = RESID[ elem ] + beta * PVEC[ elem ];    // PVEC_new = RESID_new + beta * PVEC_old
   }
   rdotr = new_rdotr;
   rnorm = sqrt( rdotr );
   if ( print ){ cout << "ConjugateGradient : After " << num_matvec << " matrix-vector products, the residual of p*O*p * x = p*RHS is " << rnorm << endl; }

}

void CheMPS2::ConjugateGradient::stepY2Z(){

   rnorm = 0.0;
   for ( int elem = 0; elem < veclength; elem++ ){
      const double diff = OPVEC[ elem ] - RHS[ elem ];
      rnorm += diff * diff;
   }
   rnorm = sqrt( rnorm );
   if ( print ){ cout << "ConjugateGradient : At convergence the residual of O * x = RHS is " << rnorm << endl; }

}

void CheMPS2::ConjugateGradient::stepJ2K(){

   apply_precon( OPVEC );                            // OPVEC = ( PRECON * operator * PRECON ) * XVEC
   for ( int elem = 0; elem < veclength; elem++ ){
      RESID[ elem ] = RESID[ elem ] - OPVEC[ elem ]; // RESID = ( precon * RHS ) - ( precon * operator * precon ) * XVEC
   }
   for ( int elem = 0; elem < veclength; elem++ ){
      PVEC[ elem ] = RESID[ elem ];                  // PVEC = RESID
   }
   rdotr = inprod( RESID );
   rnorm = sqrt( rdotr );

}

void CheMPS2::ConjugateGradient::stepG2H(){

   // PRECON = 1 / sqrt( diag ( operator ) )
   for ( int elem = 0; elem < veclength; elem++ ){
      if ( PRECON[ elem ] < DIAG_CUTOFF ){ PRECON[ elem ] = DIAG_CUTOFF; }
      PRECON[ elem ] = 1.0 / sqrt( PRECON[ elem ] );
   }

   // RESID = PRECON * RHS
   apply_precon( RHS, RESID );

   // XVEC = guess / PRECON
   for ( int elem = 0; elem < veclength; elem++ ){
      XVEC[ elem ] = XVEC[ elem ] / PRECON[ elem ];
   }

}

double CheMPS2::ConjugateGradient::inprod( double * vector ){

   double inproduct = 0.0;
   for ( int elem = 0; elem < veclength; elem++ ){
      inproduct += vector[ elem ] * vector[ elem ];
   }
   return inproduct;

}

double CheMPS2::ConjugateGradient::inprod( double * vector, double * othervector ){

   double inproduct = 0.0;
   for ( int elem = 0; elem < veclength; elem++ ){
      inproduct += vector[ elem ] * othervector[ elem ];
   }
   return inproduct;

}

void CheMPS2::ConjugateGradient::apply_precon( double * vector ){

   for ( int elem = 0; elem < veclength; elem++ ){
      vector[ elem ] = PRECON[ elem ] * vector[ elem ];
   }

}

void CheMPS2::ConjugateGradient::apply_precon( double * vector, double * result ){

   for ( int elem = 0; elem < veclength; elem++ ){
      result[ elem ] = PRECON[ elem ] * vector[ elem ];
   }

}

