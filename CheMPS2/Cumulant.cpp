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
#include <sys/time.h>
#include <unistd.h>
#include <iostream>

#include "Cumulant.h"

CheMPS2::Cumulant::Cumulant(){ }

CheMPS2::Cumulant::~Cumulant(){ }

void CheMPS2::Cumulant::gamma4_fock_contract_ham_slow(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, double * fock, double * result){

   struct timeval start, end;
   gettimeofday(&start, NULL);
   const int L = prob->gL();
   
   /*
   for ( int i = 0; i < L; i++ ){
      for ( int j = 0; j < L; j++ ){
         for ( int k = 0; k < L; k++ ){
            for ( int p = 0; p < L; p++ ){
               for ( int q = 0; q < L; q++ ){
                  for ( int r = 0; r < L; r++ ){
                     double value = 0.0;
                     for ( int l = 0; l < L; l++ ){
                        for ( int s = 0; s < L; s++ ){
                           value += fock[ l + L * s ] * gamma4_ham(prob, the3DM, the2DM, i, j, k, l, p, q, r, s);
                        }
                     }
                     result[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r ) ) ) ) ] = value;
                  }
               }
            }
         }
      }
   }
   */
   
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         for ( int k = i; k < L; k++ ){
            for ( int p = i; p < L; p++ ){
               for ( int q = i; q < L; q++ ){
                  for ( int r = q; r < L; r++ ){
                     double value = 0.0;
                     #pragma omp parallel for schedule(static) reduction(+:value)
                     for ( int ls = 0; ls < L*L; ls++ ){
                        const int l_val = ls % L;
                        const int s_val = ls / L;
                           value += fock[ ls ] * gamma4_ham(prob, the3DM, the2DM, i, j, k, l_val, p, q, r, s_val);
                     }
                     result[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r ) ) ) ) ] = value;
                     result[ i + L * ( k + L * ( j + L * ( p + L * ( r + L * q ) ) ) ) ] = value;
                     result[ j + L * ( i + L * ( k + L * ( q + L * ( p + L * r ) ) ) ) ] = value;
                     result[ k + L * ( i + L * ( j + L * ( r + L * ( p + L * q ) ) ) ) ] = value;
                     result[ j + L * ( k + L * ( i + L * ( q + L * ( r + L * p ) ) ) ) ] = value;
                     result[ k + L * ( j + L * ( i + L * ( r + L * ( q + L * p ) ) ) ) ] = value;
                     result[ p + L * ( q + L * ( r + L * ( i + L * ( j + L * k ) ) ) ) ] = value;
                     result[ p + L * ( r + L * ( q + L * ( i + L * ( k + L * j ) ) ) ) ] = value;
                     result[ q + L * ( p + L * ( r + L * ( j + L * ( i + L * k ) ) ) ) ] = value;
                     result[ r + L * ( p + L * ( q + L * ( k + L * ( i + L * j ) ) ) ) ] = value;
                     result[ q + L * ( r + L * ( p + L * ( j + L * ( k + L * i ) ) ) ) ] = value;
                     result[ r + L * ( q + L * ( p + L * ( k + L * ( j + L * i ) ) ) ) ] = value;
                  }
               }
            }
         }
      }
   }
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   std::cout << "Cumulant :: Contraction of cu(4)-4RDM with CASPT2 Fock operator took " << elapsed << " seconds." << std::endl;

}

double CheMPS2::Cumulant::lambda2_ham(const TwoDM * the2DM, const int i, const int j, const int p, const int q){

   const double value = the2DM->getTwoDMA_HAM( i, j, p, q )
                      - the2DM->get1RDM_HAM( i, p ) * the2DM->get1RDM_HAM( j, q )
                      + the2DM->get1RDM_HAM( i, q ) * the2DM->get1RDM_HAM( j, p ) * 0.5;
   return value;

}

double CheMPS2::Cumulant::gamma4_ham(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, const int i, const int j, const int k, const int l,
                                                                                                         const int p, const int q, const int r, const int s){
   
   //Prob assumes you use DMRG orbs... f1 converts HAM orbs to DMRG orbs
   const int irrep_i = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( i ) : i );
   const int irrep_j = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( j ) : j );
   const int irrep_k = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( k ) : k );
   const int irrep_l = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( l ) : l );
   const int irrep_p = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( p ) : p );
   const int irrep_q = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( q ) : q );
   const int irrep_r = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( r ) : r );
   const int irrep_s = prob->gIrrep(( prob->gReorderD2h() ) ? prob->gf1( s ) : s );
   
   const int irrep_ij = Irreps::directProd( irrep_i, irrep_j );
   const int irrep_kl = Irreps::directProd( irrep_k, irrep_l );
   const int irrep_pq = Irreps::directProd( irrep_p, irrep_q );
   const int irrep_rs = Irreps::directProd( irrep_r, irrep_s );
   
   if ( Irreps::directProd( irrep_ij, irrep_kl ) != Irreps::directProd( irrep_pq, irrep_rs ) ){ return 0.0; }
   
   const double part1 = (       the3DM->get_ham_index(i,j,k,p,q,r) * the2DM->get1RDM_HAM(l,s)
                        - 0.5 * the3DM->get_ham_index(i,j,k,s,q,r) * the2DM->get1RDM_HAM(l,p)
                        - 0.5 * the3DM->get_ham_index(i,j,k,p,s,r) * the2DM->get1RDM_HAM(l,q)
                        - 0.5 * the3DM->get_ham_index(i,j,k,p,q,s) * the2DM->get1RDM_HAM(l,r)
                        +       the3DM->get_ham_index(i,j,l,p,q,s) * the2DM->get1RDM_HAM(k,r)
                        - 0.5 * the3DM->get_ham_index(i,j,l,r,q,s) * the2DM->get1RDM_HAM(k,p)
                        - 0.5 * the3DM->get_ham_index(i,j,l,p,r,s) * the2DM->get1RDM_HAM(k,q)
                        - 0.5 * the3DM->get_ham_index(i,j,l,p,q,r) * the2DM->get1RDM_HAM(k,s)
                        +       the3DM->get_ham_index(i,k,l,p,r,s) * the2DM->get1RDM_HAM(j,q)
                        - 0.5 * the3DM->get_ham_index(i,k,l,q,r,s) * the2DM->get1RDM_HAM(j,p)
                        - 0.5 * the3DM->get_ham_index(i,k,l,p,q,s) * the2DM->get1RDM_HAM(j,r)
                        - 0.5 * the3DM->get_ham_index(i,k,l,p,r,q) * the2DM->get1RDM_HAM(j,s)
                        +       the3DM->get_ham_index(j,k,l,q,r,s) * the2DM->get1RDM_HAM(i,p)
                        - 0.5 * the3DM->get_ham_index(j,k,l,p,r,s) * the2DM->get1RDM_HAM(i,q)
                        - 0.5 * the3DM->get_ham_index(j,k,l,q,p,s) * the2DM->get1RDM_HAM(i,r)
                        - 0.5 * the3DM->get_ham_index(j,k,l,q,r,p) * the2DM->get1RDM_HAM(i,s) );

   const double part2 = ( the2DM->getTwoDMA_HAM(i,j,p,q) * the2DM->getTwoDMA_HAM(k,l,r,s)
                        - the2DM->getTwoDMA_HAM(i,j,p,r) * the2DM->getTwoDMA_HAM(k,l,q,s) * 0.5
                        - the2DM->getTwoDMA_HAM(i,j,p,s) * the2DM->getTwoDMA_HAM(k,l,r,q) * 0.5
                        - the2DM->getTwoDMA_HAM(i,j,r,q) * the2DM->getTwoDMA_HAM(k,l,p,s) * 0.5
                        - the2DM->getTwoDMA_HAM(i,j,s,q) * the2DM->getTwoDMA_HAM(k,l,r,p) * 0.5
                        + the2DM->getTwoDMA_HAM(i,j,r,s) * the2DM->getTwoDMA_HAM(k,l,p,q) / 3.0
                        + the2DM->getTwoDMA_HAM(i,j,r,s) * the2DM->getTwoDMA_HAM(k,l,q,p) / 6.0
                        + the2DM->getTwoDMA_HAM(i,j,s,r) * the2DM->getTwoDMA_HAM(k,l,p,q) / 6.0
                        + the2DM->getTwoDMA_HAM(i,j,s,r) * the2DM->getTwoDMA_HAM(k,l,q,p) / 3.0
                        + the2DM->getTwoDMA_HAM(i,k,p,r) * the2DM->getTwoDMA_HAM(j,l,q,s)
                        - the2DM->getTwoDMA_HAM(i,k,p,q) * the2DM->getTwoDMA_HAM(j,l,r,s) * 0.5
                        - the2DM->getTwoDMA_HAM(i,k,p,s) * the2DM->getTwoDMA_HAM(j,l,q,r) * 0.5
                        - the2DM->getTwoDMA_HAM(i,k,q,r) * the2DM->getTwoDMA_HAM(j,l,p,s) * 0.5
                        - the2DM->getTwoDMA_HAM(i,k,s,r) * the2DM->getTwoDMA_HAM(j,l,q,p) * 0.5
                        + the2DM->getTwoDMA_HAM(i,k,q,s) * the2DM->getTwoDMA_HAM(j,l,p,r) / 3.0
                        + the2DM->getTwoDMA_HAM(i,k,s,q) * the2DM->getTwoDMA_HAM(j,l,p,r) / 6.0
                        + the2DM->getTwoDMA_HAM(i,k,q,s) * the2DM->getTwoDMA_HAM(j,l,r,p) / 6.0
                        + the2DM->getTwoDMA_HAM(i,k,s,q) * the2DM->getTwoDMA_HAM(j,l,r,p) / 3.0
                        + the2DM->getTwoDMA_HAM(i,l,p,s) * the2DM->getTwoDMA_HAM(k,j,r,q)
                        - the2DM->getTwoDMA_HAM(i,l,p,r) * the2DM->getTwoDMA_HAM(k,j,s,q) * 0.5
                        - the2DM->getTwoDMA_HAM(i,l,p,q) * the2DM->getTwoDMA_HAM(k,j,r,s) * 0.5
                        - the2DM->getTwoDMA_HAM(i,l,r,s) * the2DM->getTwoDMA_HAM(k,j,p,q) * 0.5
                        - the2DM->getTwoDMA_HAM(i,l,q,s) * the2DM->getTwoDMA_HAM(k,j,r,p) * 0.5
                        + the2DM->getTwoDMA_HAM(i,l,r,q) * the2DM->getTwoDMA_HAM(k,j,p,s) / 3.0
                        + the2DM->getTwoDMA_HAM(i,l,q,r) * the2DM->getTwoDMA_HAM(k,j,p,s) / 6.0
                        + the2DM->getTwoDMA_HAM(i,l,r,q) * the2DM->getTwoDMA_HAM(k,j,s,p) / 6.0
                        + the2DM->getTwoDMA_HAM(i,l,q,r) * the2DM->getTwoDMA_HAM(k,j,s,p) / 3.0 );
                                    
   const double part3 = ( lambda2_ham(the2DM,i,j,p,q) * lambda2_ham(the2DM,k,l,r,s)
                        - lambda2_ham(the2DM,i,j,p,r) * lambda2_ham(the2DM,k,l,q,s) * 0.5
                        - lambda2_ham(the2DM,i,j,p,s) * lambda2_ham(the2DM,k,l,r,q) * 0.5
                        - lambda2_ham(the2DM,i,j,r,q) * lambda2_ham(the2DM,k,l,p,s) * 0.5
                        - lambda2_ham(the2DM,i,j,s,q) * lambda2_ham(the2DM,k,l,r,p) * 0.5
                        + lambda2_ham(the2DM,i,j,r,s) * lambda2_ham(the2DM,k,l,p,q) / 3.0
                        + lambda2_ham(the2DM,i,j,r,s) * lambda2_ham(the2DM,k,l,q,p) / 6.0
                        + lambda2_ham(the2DM,i,j,s,r) * lambda2_ham(the2DM,k,l,p,q) / 6.0
                        + lambda2_ham(the2DM,i,j,s,r) * lambda2_ham(the2DM,k,l,q,p) / 3.0
                        + lambda2_ham(the2DM,i,k,p,r) * lambda2_ham(the2DM,j,l,q,s)
                        - lambda2_ham(the2DM,i,k,p,q) * lambda2_ham(the2DM,j,l,r,s) * 0.5
                        - lambda2_ham(the2DM,i,k,p,s) * lambda2_ham(the2DM,j,l,q,r) * 0.5
                        - lambda2_ham(the2DM,i,k,q,r) * lambda2_ham(the2DM,j,l,p,s) * 0.5
                        - lambda2_ham(the2DM,i,k,s,r) * lambda2_ham(the2DM,j,l,q,p) * 0.5
                        + lambda2_ham(the2DM,i,k,q,s) * lambda2_ham(the2DM,j,l,p,r) / 3.0
                        + lambda2_ham(the2DM,i,k,s,q) * lambda2_ham(the2DM,j,l,p,r) / 6.0
                        + lambda2_ham(the2DM,i,k,q,s) * lambda2_ham(the2DM,j,l,r,p) / 6.0
                        + lambda2_ham(the2DM,i,k,s,q) * lambda2_ham(the2DM,j,l,r,p) / 3.0
                        + lambda2_ham(the2DM,i,l,p,s) * lambda2_ham(the2DM,k,j,r,q)
                        - lambda2_ham(the2DM,i,l,p,r) * lambda2_ham(the2DM,k,j,s,q) * 0.5
                        - lambda2_ham(the2DM,i,l,p,q) * lambda2_ham(the2DM,k,j,r,s) * 0.5
                        - lambda2_ham(the2DM,i,l,r,s) * lambda2_ham(the2DM,k,j,p,q) * 0.5
                        - lambda2_ham(the2DM,i,l,q,s) * lambda2_ham(the2DM,k,j,r,p) * 0.5
                        + lambda2_ham(the2DM,i,l,r,q) * lambda2_ham(the2DM,k,j,p,s) / 3.0
                        + lambda2_ham(the2DM,i,l,q,r) * lambda2_ham(the2DM,k,j,p,s) / 6.0
                        + lambda2_ham(the2DM,i,l,r,q) * lambda2_ham(the2DM,k,j,s,p) / 6.0
                        + lambda2_ham(the2DM,i,l,q,r) * lambda2_ham(the2DM,k,j,s,p) / 3.0 );

   return ( part1 - part2 + 2 * part3 );

}




