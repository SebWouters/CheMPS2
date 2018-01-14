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
#include <sys/time.h>
#include <unistd.h>
#include <iostream>

#include "Cumulant.h"

/*void CheMPS2::Cumulant::gamma4_fock_contract_ham_slow(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, double * fock, double * result){

   struct timeval start, end;
   gettimeofday(&start, NULL);
   const int L = prob->gL();

   for ( int cnt = 0; cnt < L*L*L*L*L*L; cnt++ ){ result[ cnt ] = 0.0; }

   int * irreps = new int[ L ];
   for ( int orb = 0; orb < L; orb++ ){ irreps[ orb ] = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( orb ) : orb ); }

   #pragma omp parallel for schedule(dynamic)
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         for ( int k = i; k < L; k++ ){
            const int irrep_ijk = Irreps::directProd( Irreps::directProd( irreps[ i ], irreps[ j ] ), irreps[ k ] );
            for ( int p = i; p < L; p++ ){
               for ( int q = i; q < L; q++ ){
                  for ( int r = q; r < L; r++ ){
                     const int irrep_pqr = Irreps::directProd( Irreps::directProd( irreps[ p ], irreps[ q ] ), irreps[ r ] );
                     if ( irrep_ijk == irrep_pqr ){
                        double value = 0.0;
                        for ( int l = 0; l < L; l++ ){
                           for ( int s = 0; s < L; s++ ){
                              if ( irreps[ l ] == irreps[ s ] ){
                                 value += fock[ l + L * s ] * gamma4_ham(prob, the3DM, the2DM, i, j, k, l, p, q, r, s);
                              }
                           }
                        }
                        result[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] = value;
                        result[ i + L * ( k + L * ( j + L * ( p + L * ( r + L * q )))) ] = value;
                        result[ j + L * ( i + L * ( k + L * ( q + L * ( p + L * r )))) ] = value;
                        result[ k + L * ( i + L * ( j + L * ( r + L * ( p + L * q )))) ] = value;
                        result[ j + L * ( k + L * ( i + L * ( q + L * ( r + L * p )))) ] = value;
                        result[ k + L * ( j + L * ( i + L * ( r + L * ( q + L * p )))) ] = value;
                        result[ p + L * ( q + L * ( r + L * ( i + L * ( j + L * k )))) ] = value;
                        result[ p + L * ( r + L * ( q + L * ( i + L * ( k + L * j )))) ] = value;
                        result[ q + L * ( p + L * ( r + L * ( j + L * ( i + L * k )))) ] = value;
                        result[ r + L * ( p + L * ( q + L * ( k + L * ( i + L * j )))) ] = value;
                        result[ q + L * ( r + L * ( p + L * ( j + L * ( k + L * i )))) ] = value;
                        result[ r + L * ( q + L * ( p + L * ( k + L * ( j + L * i )))) ] = value;
                     }
                  }
               }
            }
         }
      }
   }
   
   delete [] irreps;
   
   gettimeofday(&end, NULL);
   const double elapsed = (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
   std::cout << "Cumulant :: Contraction of cu(4)-4RDM with CASPT2 Fock operator took " << elapsed << " seconds." << std::endl;

}*/

void CheMPS2::Cumulant::gamma4_fock_contract_ham(const Problem * prob, const ThreeDM * the3DM, const TwoDM * the2DM, double * fock, double * result){

   struct timeval start, end;
   gettimeofday(&start, NULL);
   const int L = prob->gL();
   
   /* Clear result */
   for ( int cnt = 0; cnt < L*L*L*L*L*L; cnt++ ){ result[ cnt ] = 0.0; }
   
   /* Construct an array with the orbital irreps in Hamiltonian indices */
   int * irreps = new int[ L ];
   for ( int orb = 0; orb < L; orb++ ){ irreps[ orb ] = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( orb ) : orb ); }
   
   /* Helper arrays with partial (mult) or full (dot) contractions of objects with the CASPT2 Fock operator */
   double * G3dotF  = new double[ L*L*L*L ];
   double * lambda2 = new double[ L*L*L*L ];
   double * G2multF = new double[ L*L*L*L ];
   double * L2multF = new double[ L*L*L*L ];
   double * gamma1  = new double[ L*L ];
   double * G1multF = new double[ L*L ];
   double * G2dotF  = new double[ L*L ];
   double * L2dotF  = new double[ L*L ];
   for ( int cnt = 0; cnt < L*L*L*L; cnt++ ){  G3dotF[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L*L*L; cnt++ ){ lambda2[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L*L*L; cnt++ ){ G2multF[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L*L*L; cnt++ ){ L2multF[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L;     cnt++ ){  gamma1[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L;     cnt++ ){ G1multF[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L;     cnt++ ){  G2dotF[ cnt ] = 0.0; }
   for ( int cnt = 0; cnt < L*L;     cnt++ ){  L2dotF[ cnt ] = 0.0; }
   double G1dotF = 0.0;
   
   /* Fill gamma1 */
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         const double value = the2DM->get1RDM_HAM( i, j );
         gamma1[ i + L * j ] = value;
         gamma1[ j + L * i ] = value;
      }
   }
   
   /* Build  G3dotF[i,j,p,q] = sum_[l,s] Gamma3[i,j,l,p,q,s] F[l,s]
      Build lambda2[i,j,p,q] = Lambda2[i,j,p,q]                     */
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         const int irrep_ij = Irreps::directProd( irreps[ i ], irreps[ j ] );
         for ( int p = i; p < L; p++ ){
            for ( int q = i; q < L; q++ ){
               const int irrep_pq = Irreps::directProd( irreps[ p ], irreps[ q ] );
               if ( irrep_ij == irrep_pq ){
                  { // sum_[l,s] Gamma3[i,j,l,p,q,s] F[l,s]
                     double value = 0.0;
                     for ( int l = 0; l < L; l++ ){
                        for ( int s = 0; s < L; s++ ){
                           if ( irreps[ l ] == irreps[ s ] ){
                              value += the3DM->get_ham_index(i,j,l,p,q,s) * fock[ l + L * s ];
                           }
                        }
                     }
                     G3dotF[ i + L * ( j + L * ( p + L * q )) ] = value;
                     G3dotF[ j + L * ( i + L * ( q + L * p )) ] = value;
                     G3dotF[ p + L * ( q + L * ( i + L * j )) ] = value;
                     G3dotF[ q + L * ( p + L * ( j + L * i )) ] = value;
                  }
                  { // Lambda2[i,j,p,q]
                     const double val_lambda = the2DM->getTwoDMA_HAM( i, j, p, q )
                                             - gamma1[ i + L * p ] * gamma1[ j + L * q ]
                                             + gamma1[ i + L * q ] * gamma1[ j + L * p ] * 0.5;
                     lambda2[ i + L * ( j + L * ( p + L * q )) ] = val_lambda;
                     lambda2[ j + L * ( i + L * ( q + L * p )) ] = val_lambda;
                     lambda2[ p + L * ( q + L * ( i + L * j )) ] = val_lambda;
                     lambda2[ q + L * ( p + L * ( j + L * i )) ] = val_lambda;
                  }
               }
            }
         }
      }
   }
   
   /* Build G1multF[i,j] = sum_[p]   Gamma1[i,p] F[p,j]
      Build       G1dotF = sum_[i,j] Gamma1[i,j] F[j,i] */
   for ( int i = 0; i < L; i++ ){
      for ( int j = 0; j < L; j++ ){
         if ( irreps[ i ] == irreps[ j ] ){
            double value = 0.0;
            for ( int p = 0; p < L; p++ ){
               if ( irreps[ j ] == irreps[ p ] ){
                  value += gamma1[ i + L * p ] * fock[ p + L * j ];
               }
            }
            G1multF[ i + L * j ] = value;
         }
      }
      G1dotF += G1multF[ i * ( L + 1 ) ];
   }
   
   /* Build G2dotF[i,j] = sum_[p,q]  Gamma2[i,p,j,q] F[p,q]
      Build L2dotF[i,j] = sum_[p,q] Lambda2[i,p,j,q] F[p,q] */
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         if ( irreps[ i ] == irreps[ j ] ){
            double val_gamma  = 0.0;
            double val_lambda = 0.0;
            for ( int p = 0; p < L; p++ ){
               for ( int q = 0; q < L; q++ ){
                  if ( irreps[ p ] == irreps[ q ] ){
                     val_gamma  += fock[ p + L * q ] * the2DM->getTwoDMA_HAM( i, p, j, q );
                     val_lambda += fock[ p + L * q ] * lambda2[ i + L * ( p + L * ( j + L * q )) ];
                  }
               }
            }
            G2dotF[ i + L * j ] = val_gamma;
            G2dotF[ j + L * i ] = val_gamma;
            L2dotF[ i + L * j ] = val_lambda;
            L2dotF[ j + L * i ] = val_lambda;
         }
      }
   }
   
   /* Build G2multF[i,s,j,k] = sum_[l]  Gamma2[i,l,j,k] F[l,s]
      Build L2multF[i,s,j,k] = sum_[l] Lambda2[i,l,j,k] F[l,s] */
   for ( int i = 0; i < L; i++ ){
      for ( int j = 0; j < L; j++ ){
         const int irrep_ij = Irreps::directProd( irreps[ i ], irreps[ j ] );
         for ( int k = 0; k < L; k++ ){
            for ( int s = 0; s < L; s++ ){
               const int irrep_ks = Irreps::directProd( irreps[ k ], irreps[ s ] );
               if ( irrep_ij == irrep_ks ){
                  double val_gamma  = 0.0;
                  double val_lambda = 0.0;
                  for ( int l = 0; l < L; l++ ){
                     if ( irreps[ l ] == irreps[ s ] ){
                        val_gamma  += fock[ l + L * s ] * the2DM->getTwoDMA_HAM( i, l, j, k );
                        val_lambda += fock[ l + L * s ] * lambda2[ i + L * ( l + L * ( j + L * k ) ) ];
                     }
                  }
                  G2multF[ i + L * ( s + L * ( j + L * k ) ) ] = val_gamma;
                  L2multF[ i + L * ( s + L * ( j + L * k ) ) ] = val_lambda;
               }
            }
         }
      }
   }
   
   /* Fill result */
   #pragma omp parallel for schedule(dynamic)
   for ( int i = 0; i < L; i++ ){
      for ( int j = i; j < L; j++ ){
         for ( int k = i; k < L; k++ ){
            const int irrep_ijk = Irreps::directProd( Irreps::directProd( irreps[ i ], irreps[ j ] ), irreps[ k ] );
            for ( int p = i; p < L; p++ ){
               for ( int q = i; q < L; q++ ){
                  for ( int r = q; r < L; r++ ){
                     const int irrep_pqr = Irreps::directProd( Irreps::directProd( irreps[ p ], irreps[ q ] ), irreps[ r ] );
                     if ( irrep_ijk == irrep_pqr ){
                     
                        double dm3_contribution = 0.0;
                        double gamma2_part1 = 0.0;
                        double gamma2_part2 = 0.0;
                        double lambda2_part1 = 0.0;
                        double lambda2_part2 = 0.0;
                        
                        for ( int ls = 0; ls < L; ls++ ){
                        
                           dm3_contribution += ( the3DM->get_ham_index( ls, j, k, p, q, r ) * G1multF[ i + L * ls ]
                                               + the3DM->get_ham_index( i, ls, k, p, q, r ) * G1multF[ j + L * ls ]
                                               + the3DM->get_ham_index( i, j, ls, p, q, r ) * G1multF[ k + L * ls ]
                                               + the3DM->get_ham_index( i, j, k, ls, q, r ) * G1multF[ p + L * ls ]
                                               + the3DM->get_ham_index( i, j, k, p, ls, r ) * G1multF[ q + L * ls ]
                                               + the3DM->get_ham_index( i, j, k, p, q, ls ) * G1multF[ r + L * ls ] );
                           
                           gamma2_part1 += ( the2DM->getTwoDMA_HAM( i, j, p, ls ) * G2multF[ k + L * ( ls + L * ( r + L * q )) ]
                                           + the2DM->getTwoDMA_HAM( i, j, ls, q ) * G2multF[ k + L * ( ls + L * ( r + L * p )) ]
                                           + the2DM->getTwoDMA_HAM( i, k, p, ls ) * G2multF[ j + L * ( ls + L * ( q + L * r )) ]
                                           + the2DM->getTwoDMA_HAM( i, k, ls, r ) * G2multF[ j + L * ( ls + L * ( q + L * p )) ]
                                           + the2DM->getTwoDMA_HAM( k, j, ls, q ) * G2multF[ i + L * ( ls + L * ( p + L * r )) ]
                                           + the2DM->getTwoDMA_HAM( k, j, r, ls ) * G2multF[ i + L * ( ls + L * ( p + L * q )) ] );
                           
                           gamma2_part2 += ( the2DM->getTwoDMA_HAM( i, j, r, ls ) * ( G2multF[ k + L * ( ls + L * ( p + L * q )) ] 
                                                                              + 0.5 * G2multF[ k + L * ( ls + L * ( q + L * p )) ] )
                                           + the2DM->getTwoDMA_HAM( i, j, ls, r ) * ( G2multF[ k + L * ( ls + L * ( q + L * p )) ]
                                                                              + 0.5 * G2multF[ k + L * ( ls + L * ( p + L * q )) ] )
                                           + the2DM->getTwoDMA_HAM( i, k, q, ls ) * ( G2multF[ j + L * ( ls + L * ( p + L * r )) ]
                                                                              + 0.5 * G2multF[ j + L * ( ls + L * ( r + L * p )) ] )
                                           + the2DM->getTwoDMA_HAM( i, k, ls, q ) * ( G2multF[ j + L * ( ls + L * ( r + L * p )) ]
                                                                              + 0.5 * G2multF[ j + L * ( ls + L * ( p + L * r )) ] )
                                           + the2DM->getTwoDMA_HAM( k, j, p, ls ) * ( G2multF[ i + L * ( ls + L * ( r + L * q )) ]
                                                                              + 0.5 * G2multF[ i + L * ( ls + L * ( q + L * r )) ] )
                                           + the2DM->getTwoDMA_HAM( k, j, ls, p ) * ( G2multF[ i + L * ( ls + L * ( q + L * r )) ]
                                                                              + 0.5 * G2multF[ i + L * ( ls + L * ( r + L * q )) ] ) );
                           
                           lambda2_part1 += ( lambda2[ i + L * ( j + L * ( p + L * ls )) ] * L2multF[ k + L * ( ls + L * ( r + L * q )) ]
                                            + lambda2[ i + L * ( j + L * ( ls + L * q )) ] * L2multF[ k + L * ( ls + L * ( r + L * p )) ]
                                            + lambda2[ i + L * ( k + L * ( p + L * ls )) ] * L2multF[ j + L * ( ls + L * ( q + L * r )) ]
                                            + lambda2[ i + L * ( k + L * ( ls + L * r )) ] * L2multF[ j + L * ( ls + L * ( q + L * p )) ]
                                            + lambda2[ k + L * ( j + L * ( ls + L * q )) ] * L2multF[ i + L * ( ls + L * ( p + L * r )) ]
                                            + lambda2[ k + L * ( j + L * ( r + L * ls )) ] * L2multF[ i + L * ( ls + L * ( p + L * q )) ] );
                           
                           lambda2_part2 += ( lambda2[ i + L * ( j + L * ( r + L * ls )) ] * ( L2multF[ k + L * ( ls + L * ( p + L * q )) ] 
                                                                                       + 0.5 * L2multF[ k + L * ( ls + L * ( q + L * p )) ] )
                                            + lambda2[ i + L * ( j + L * ( ls + L * r )) ] * ( L2multF[ k + L * ( ls + L * ( q + L * p )) ]
                                                                                       + 0.5 * L2multF[ k + L * ( ls + L * ( p + L * q )) ] )
                                            + lambda2[ i + L * ( k + L * ( q + L * ls )) ] * ( L2multF[ j + L * ( ls + L * ( p + L * r )) ]
                                                                                       + 0.5 * L2multF[ j + L * ( ls + L * ( r + L * p )) ] )
                                            + lambda2[ i + L * ( k + L * ( ls + L * q )) ] * ( L2multF[ j + L * ( ls + L * ( r + L * p )) ]
                                                                                       + 0.5 * L2multF[ j + L * ( ls + L * ( p + L * r )) ] )
                                            + lambda2[ k + L * ( j + L * ( p + L * ls )) ] * ( L2multF[ i + L * ( ls + L * ( r + L * q )) ]
                                                                                       + 0.5 * L2multF[ i + L * ( ls + L * ( q + L * r )) ] )
                                            + lambda2[ k + L * ( j + L * ( ls + L * p )) ] * ( L2multF[ i + L * ( ls + L * ( q + L * r )) ]
                                                                                       + 0.5 * L2multF[ i + L * ( ls + L * ( r + L * q )) ] ) );
                        }
                     
                        const double contracted_value = ( the3DM->get_ham_index( i, j, k, p, q, r ) * G1dotF

                                                  +       G3dotF[ i + L * ( j + L * ( p + L * q )) ] * gamma1[ k + L * r ]
                                                  - 0.5 * G3dotF[ i + L * ( j + L * ( r + L * q )) ] * gamma1[ k + L * p ]
                                                  - 0.5 * G3dotF[ i + L * ( j + L * ( p + L * r )) ] * gamma1[ k + L * q ]
                                                  +       G3dotF[ i + L * ( k + L * ( p + L * r )) ] * gamma1[ j + L * q ]
                                                  - 0.5 * G3dotF[ i + L * ( k + L * ( q + L * r )) ] * gamma1[ j + L * p ]
                                                  - 0.5 * G3dotF[ i + L * ( k + L * ( p + L * q )) ] * gamma1[ j + L * r ]
                                                  +       G3dotF[ j + L * ( k + L * ( q + L * r )) ] * gamma1[ i + L * p ]
                                                  - 0.5 * G3dotF[ j + L * ( k + L * ( p + L * r )) ] * gamma1[ i + L * q ]
                                                  - 0.5 * G3dotF[ j + L * ( k + L * ( q + L * p )) ] * gamma1[ i + L * r ]

                                                  - 0.5 * dm3_contribution

                                                  -       the2DM->getTwoDMA_HAM( i, j, p, q ) * G2dotF[ k + L * r ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( i, j, p, r ) * G2dotF[ k + L * q ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( i, j, r, q ) * G2dotF[ k + L * p ]
                                                  -       the2DM->getTwoDMA_HAM( i, k, p, r ) * G2dotF[ j + L * q ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( i, k, p, q ) * G2dotF[ j + L * r ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( i, k, q, r ) * G2dotF[ j + L * p ]
                                                  -       the2DM->getTwoDMA_HAM( k, j, r, q ) * G2dotF[ i + L * p ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( k, j, p, q ) * G2dotF[ i + L * r ]
                                                  + 0.5 * the2DM->getTwoDMA_HAM( k, j, r, p ) * G2dotF[ i + L * q ]

                                                  + 0.5 * gamma2_part1
                                                  -       gamma2_part2 / 3.0

                                                  + 2 * lambda2[ i + L * ( j + L * ( p + L * q )) ] * L2dotF[ k + L * r ]
                                                  -     lambda2[ i + L * ( j + L * ( p + L * r )) ] * L2dotF[ k + L * q ]
                                                  -     lambda2[ i + L * ( j + L * ( r + L * q )) ] * L2dotF[ k + L * p ]
                                                  + 2 * lambda2[ i + L * ( k + L * ( p + L * r )) ] * L2dotF[ j + L * q ]
                                                  -     lambda2[ i + L * ( k + L * ( p + L * q )) ] * L2dotF[ j + L * r ]
                                                  -     lambda2[ i + L * ( k + L * ( q + L * r )) ] * L2dotF[ j + L * p ]
                                                  + 2 * lambda2[ k + L * ( j + L * ( r + L * q )) ] * L2dotF[ i + L * p ]
                                                  -     lambda2[ k + L * ( j + L * ( p + L * q )) ] * L2dotF[ i + L * r ]
                                                  -     lambda2[ k + L * ( j + L * ( r + L * p )) ] * L2dotF[ i + L * q ]

                                                  - lambda2_part1
                                                  + lambda2_part2 / 1.5 );
                              
                        result[ i + L * ( j + L * ( k + L * ( p + L * ( q + L * r )))) ] = contracted_value;
                        result[ i + L * ( k + L * ( j + L * ( p + L * ( r + L * q )))) ] = contracted_value;
                        result[ j + L * ( i + L * ( k + L * ( q + L * ( p + L * r )))) ] = contracted_value;
                        result[ k + L * ( i + L * ( j + L * ( r + L * ( p + L * q )))) ] = contracted_value;
                        result[ j + L * ( k + L * ( i + L * ( q + L * ( r + L * p )))) ] = contracted_value;
                        result[ k + L * ( j + L * ( i + L * ( r + L * ( q + L * p )))) ] = contracted_value;
                        result[ p + L * ( q + L * ( r + L * ( i + L * ( j + L * k )))) ] = contracted_value;
                        result[ p + L * ( r + L * ( q + L * ( i + L * ( k + L * j )))) ] = contracted_value;
                        result[ q + L * ( p + L * ( r + L * ( j + L * ( i + L * k )))) ] = contracted_value;
                        result[ r + L * ( p + L * ( q + L * ( k + L * ( i + L * j )))) ] = contracted_value;
                        result[ q + L * ( r + L * ( p + L * ( j + L * ( k + L * i )))) ] = contracted_value;
                        result[ r + L * ( q + L * ( p + L * ( k + L * ( j + L * i )))) ] = contracted_value;
                     }
                  }
               }
            }
         }
      }
   }
   
   delete [] G3dotF;
   delete [] lambda2;
   delete [] G2multF;
   delete [] L2multF;
   delete [] gamma1;
   delete [] G1multF;
   delete [] G2dotF;
   delete [] L2dotF;
   delete [] irreps;
   
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
   const int irrep_i = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( i ) : i );
   const int irrep_j = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( j ) : j );
   const int irrep_k = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( k ) : k );
   const int irrep_l = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( l ) : l );
   const int irrep_p = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( p ) : p );
   const int irrep_q = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( q ) : q );
   const int irrep_r = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( r ) : r );
   const int irrep_s = prob->gIrrep(( prob->gReorder() ) ? prob->gf1( s ) : s );
   
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




