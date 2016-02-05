/*
   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2016 Sebastian Wouters

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
#include <assert.h>
#include <algorithm>
#include <sys/time.h>

#include "CASPT2.h"
#include "Lapack.h"
#include "Options.h"
#include "ConjugateGradient.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

CheMPS2::CASPT2::CASPT2( DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei, DMRGSCFmatrix * fock_in, double * one_dm, double * two_dm, double * three_dm, double * contract_4dm ){

   indices    = idx;
   fock       = fock_in;
   one_rdm    = one_dm;
   two_rdm    = two_dm;
   three_rdm  = three_dm;
   f_dot_4dm  = contract_4dm;
   num_irreps = indices->getNirreps();
   create_f_dots(); // Sets E_FOCK

   vector_helper(); // Needs to be called BEFORE make_S**()!

   make_SAA_SCC();
   make_SDD();
   make_SEE_SGG();
   make_SBB_SFF_singlet();
   make_SBB_SFF_triplet();

   construct_rhs( oei, ints ); // Needs to be called AFTER make_S**()!

   // The following Fock operator constructors need to be called AFTER make_S**()!
   make_FAA_FCC();
   make_FDD();
   make_FEE_FGG();
   make_FBB_FFF_singlet();
   make_FBB_FFF_triplet();

   make_FAD_FCD();
   make_FEH_FGH();
   make_FAB_FCF_singlet();
   make_FAB_FCF_triplet();
   make_FBE_FFG_singlet();
   make_FBE_FFG_triplet();
   make_FDE_FDG();

   delete [] f_dot_3dm;
   delete [] f_dot_2dm;

   recreate(); // Remove the overlap matrices

   /*{

      cout << "MOLCAS test8 CASPT2-N= " << -0.1599978130 << endl;
      cout << "MOLCAS test8 CASPT2-D= " << -0.1596306078 << endl;
      const double energy_caspt2_n = solve();
      cout << "MOLCAS test8 CASPT2-N= " << -0.1599978130 << endl;
      cout << "MOLCAS test8 CASPT2-D= " << -0.1596306078 << endl;

      if ( false ){
         int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
         double * vector = new double[ total_size ];
         double * matrix = new double[ total_size * total_size ];
         double * diag_fock = new double[ total_size ];
         diagonal( diag_fock, -E_FOCK );

         for ( int col = 0; col < total_size; col++ ){
            for ( int row = 0; row < total_size; row++ ){ vector[ row ] = 0.0; }
            vector[ col ] = 1.0;
            matvec( vector, matrix + col * total_size, diag_fock );
         }

         double rms = 0.0;
         for ( int col = 0; col < total_size; col++ ){
            for ( int row = col+1; row < total_size; row++ ){
               const double diff = matrix[ row + total_size * col ] - matrix[ col + total_size * row ];
               rms += diff * diff;
               if ( fabs( diff ) > 1e-6 ){ cout << " matrix[" << row << "," << col << "] - matrix[" << col << "," << row << "] = " <<  diff << endl; }
            }
         }
         rms = sqrt( rms );
         cout << "RMS deviation from symmetric = " << rms << endl;

         delete [] vector;
         delete [] matrix;
         delete [] diag_fock;
      }
   }*/

}

CheMPS2::CASPT2::~CASPT2(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      delete [] FAA[ irrep ];
      delete [] FCC[ irrep ];
      delete [] FDD[ irrep ];
      delete [] FEE[ irrep ];
      delete [] FGG[ irrep ];
      delete [] FBB_singlet[ irrep ];
      delete [] FBB_triplet[ irrep ];
      delete [] FFF_singlet[ irrep ];
      delete [] FFF_triplet[ irrep ];
   }
   delete [] FAA;
   delete [] FCC;
   delete [] FDD;
   delete [] FEE;
   delete [] FGG;
   delete [] FBB_singlet;
   delete [] FBB_triplet;
   delete [] FFF_singlet;
   delete [] FFF_triplet;

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){
      for ( int irrep_right = 0; irrep_right < num_irreps; irrep_right++ ){
         const int irrep_w = Irreps::directProd( irrep_left, irrep_right );
         const int num_w   = indices->getNDMRG( irrep_w );
         for ( int w = 0; w < num_w; w++ ){
            delete [] FAD[ irrep_left ][ irrep_right ][ w ];
            delete [] FCD[ irrep_left ][ irrep_right ][ w ];
            delete [] FAB_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FAB_triplet[ irrep_left ][ irrep_right ][ w ];
            delete [] FCF_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FCF_triplet[ irrep_left ][ irrep_right ][ w ];
            delete [] FBE_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FBE_triplet[ irrep_left ][ irrep_right ][ w ];
            delete [] FFG_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FFG_triplet[ irrep_left ][ irrep_right ][ w ];
            delete [] FDE_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FDE_triplet[ irrep_left ][ irrep_right ][ w ];
            delete [] FDG_singlet[ irrep_left ][ irrep_right ][ w ];
            delete [] FDG_triplet[ irrep_left ][ irrep_right ][ w ];
         }
         delete [] FAD[ irrep_left ][ irrep_right ];
         delete [] FCD[ irrep_left ][ irrep_right ];
         delete [] FAB_singlet[ irrep_left ][ irrep_right ];
         delete [] FAB_triplet[ irrep_left ][ irrep_right ];
         delete [] FCF_singlet[ irrep_left ][ irrep_right ];
         delete [] FCF_triplet[ irrep_left ][ irrep_right ];
         delete [] FBE_singlet[ irrep_left ][ irrep_right ];
         delete [] FBE_triplet[ irrep_left ][ irrep_right ];
         delete [] FFG_singlet[ irrep_left ][ irrep_right ];
         delete [] FFG_triplet[ irrep_left ][ irrep_right ];
         delete [] FDE_singlet[ irrep_left ][ irrep_right ];
         delete [] FDE_triplet[ irrep_left ][ irrep_right ];
         delete [] FDG_singlet[ irrep_left ][ irrep_right ];
         delete [] FDG_triplet[ irrep_left ][ irrep_right ];
      }
      delete [] FAD[ irrep_left ];
      delete [] FCD[ irrep_left ];
      delete [] FAB_singlet[ irrep_left ];
      delete [] FAB_triplet[ irrep_left ];
      delete [] FCF_singlet[ irrep_left ];
      delete [] FCF_triplet[ irrep_left ];
      delete [] FBE_singlet[ irrep_left ];
      delete [] FBE_triplet[ irrep_left ];
      delete [] FFG_singlet[ irrep_left ];
      delete [] FFG_triplet[ irrep_left ];
      delete [] FDE_singlet[ irrep_left ];
      delete [] FDE_triplet[ irrep_left ];
      delete [] FDG_singlet[ irrep_left ];
      delete [] FDG_triplet[ irrep_left ];
   }
   delete [] FAD;
   delete [] FCD;
   delete [] FAB_singlet;
   delete [] FAB_triplet;
   delete [] FCF_singlet;
   delete [] FCF_triplet;
   delete [] FBE_singlet;
   delete [] FBE_triplet;
   delete [] FFG_singlet;
   delete [] FFG_triplet;
   delete [] FDE_singlet;
   delete [] FDE_triplet;
   delete [] FDG_singlet;
   delete [] FDG_triplet;

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int num_w = indices->getNDMRG( irrep );
      for ( int w = 0; w < num_w; w++ ){
         delete [] FEH[ irrep ][ w ];
         delete [] FGH[ irrep ][ w ];
      }
      delete [] FEH[ irrep ];
      delete [] FGH[ irrep ];
   }
   delete [] FEH;
   delete [] FGH;

   delete [] size_A;
   delete [] size_C;
   delete [] size_D;
   delete [] size_E;
   delete [] size_G;
   delete [] size_B_singlet;
   delete [] size_B_triplet;
   delete [] size_F_singlet;
   delete [] size_F_triplet;
   delete [] jump;
   delete [] vector_rhs;

}

double CheMPS2::CASPT2::solve( const bool diag_only ) const{

   int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];

   double * diag_fock = new double[ total_size ];
   diagonal( diag_fock, -E_FOCK );

   double energy_caspt2_d = 0.0;
   for ( int elem = 0; elem < total_size; elem++ ){
      energy_caspt2_d -= vector_rhs[ elem ] * vector_rhs[ elem ] / diag_fock[ elem ];
   }
   cout << "CASPT2 : E(CASPT2-D) = " << energy_caspt2_d << endl;
   
   double min_eig = diag_fock[ 0 ];
   for ( int elem = 1; elem < total_size; elem++ ){ min_eig = min( min_eig, diag_fock[ elem ] ); }
   cout << "CASPT2 : Minimum value on diagonal of 0th order operator = " << min_eig << endl;

   double best_energy = energy_caspt2_d;
   if ( diag_only == false ){

      ConjugateGradient CG( total_size, CheMPS2::CONJ_GRADIENT_RTOL, CheMPS2::CONJ_GRADIENT_PRECOND_CUTOFF, false );
      double ** pointers = new double*[ 3 ];
      char instruction = CG.step( pointers );
      assert( instruction == 'A' );
      for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 0 ][ elem ] = vector_rhs[ elem ] / diag_fock[ elem ]; } // Initial guess of F * x = V
      for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 1 ][ elem ] =  diag_fock[ elem ];                     } // Diagonal of the operator F
      for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 2 ][ elem ] = vector_rhs[ elem ];                     } // RHS of the linear problem F * x = V
      int inc1 = 1;
      instruction = CG.step( pointers );
      assert( instruction == 'B' );
      while ( instruction == 'B' ){
         matvec( pointers[ 0 ], pointers[ 1 ], diag_fock );
         instruction = CG.step( pointers );
      }
      assert( instruction == 'C' );
      best_energy = - ddot_( &total_size, pointers[ 0 ], &inc1, vector_rhs, &inc1 );
      const double rnorm = pointers[ 1 ][ 0 ];
      cout << "CASPT2 : Residual norm of F * c - V = " << rnorm << endl;
      const double weight = 1.0 / ( 1.0 + ddot_( &total_size, pointers[ 0 ], &inc1, pointers[ 0 ], & inc1 ) );
      cout << "CASPT2 : Fraction of CASSCF wfn in 1st order wfn = " << weight << endl;
      delete [] pointers;
      cout << "CASPT2 : E(CASPT2-N) = " << best_energy << endl;

   }
   delete [] diag_fock;
   return best_energy;

}

void CheMPS2::CASPT2::create_f_dots(){

   const int LAS = indices->getDMRGcumulative( num_irreps );
   f_dot_3dm = new double[ LAS * LAS * LAS * LAS ];
   f_dot_2dm = new double[ LAS * LAS ];
   f_dot_1dm = 0.0;

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NDMRG = indices->getNDMRG( irrep );
      const int NOCC  = indices->getNOCC( irrep );
      const int jumpx = indices->getDMRGcumulative( irrep );

      double value = 0.0;
      for ( int row = 0; row < NDMRG; row++ ){
         for ( int col = 0; col < NDMRG; col++ ){
            value += fock->get( irrep, NOCC + row, NOCC + col ) * one_rdm[ jumpx + row + LAS * ( jumpx + col ) ];
         }
      }
      f_dot_1dm += value;
   }

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      const int lower1 = indices->getDMRGcumulative( irrep1 );
      const int upper1 = lower1 + indices->getNDMRG( irrep1 );

      for ( int i1 = lower1; i1 < upper1; i1++ ){
         for ( int i2 = lower1; i2 < upper1; i2++ ){

            double value = 0.0;
            for ( int irrepx = 0; irrepx < num_irreps; irrepx++ ){
               const int NOCCx  = indices->getNOCC( irrepx );
               const int NDMRGx = indices->getNDMRG( irrepx );
               const int jumpx  = indices->getDMRGcumulative( irrepx );

               for ( int rowx = 0; rowx < NDMRGx; rowx++ ){
                  for ( int colx = 0; colx < NDMRGx; colx++ ){
                     value += ( fock->get( irrepx, NOCCx + rowx, NOCCx + colx )
                              * two_rdm[ i1 + LAS * ( jumpx + rowx + LAS * ( i2 + LAS * ( jumpx + colx ))) ] );
                  }
               }
            }
            f_dot_2dm[ i1 + LAS * i2 ] = value;
         }
      }
   }

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      const int lower1 = indices->getDMRGcumulative( irrep1 );
      const int upper1 = lower1 + indices->getNDMRG( irrep1 );
      for ( int irrep2 = 0; irrep2 < num_irreps; irrep2++ ){
         const int lower2 = indices->getDMRGcumulative( irrep2 );
         const int upper2 = lower2 + indices->getNDMRG( irrep2 );
         const int irr_12 = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = 0; irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( irr_12, irrep3 );
            const int lower3 = indices->getDMRGcumulative( irrep3 );
            const int lower4 = indices->getDMRGcumulative( irrep4 );
            const int upper3 = lower3 + indices->getNDMRG( irrep3 );
            const int upper4 = lower4 + indices->getNDMRG( irrep4 );

            for ( int i1 = lower1; i1 < upper1; i1++ ){
               for ( int i2 = lower2; i2 < upper2; i2++ ){
                  for ( int i3 = lower3; i3 < upper3; i3++ ){
                     for ( int i4 = lower4; i4 < upper4; i4++ ){

                        double value = 0.0;
                        for ( int irrepx = 0; irrepx < num_irreps; irrepx++ ){
                           const int jumpx  = indices->getDMRGcumulative( irrepx );
                           const int NOCCx  = indices->getNOCC( irrepx );
                           const int NDMRGx = indices->getNDMRG( irrepx );

                           for ( int rowx = 0; rowx < NDMRGx; rowx++ ){
                              for ( int colx = 0; colx < NDMRGx; colx++ ){
                                 value += ( fock->get( irrepx, NOCCx + rowx, NOCCx + colx )
                                          * three_rdm[ i1 + LAS*( i2 + LAS*( jumpx + rowx + LAS*( i3 + LAS*( i4 + LAS*( jumpx + colx ))))) ] );
                              }
                           }
                        }
                        f_dot_3dm[ i1 + LAS * ( i2 + LAS * ( i3 + LAS * i4 )) ] = value;
                     }
                  }
               }
            }
         }
      }
   }

   sum_f_kk = 0.0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      for ( int orb = 0; orb < NOCC; orb++ ){
         sum_f_kk += fock->get( irrep, orb, orb );
      }
   }

   E_FOCK = 2 * sum_f_kk + f_dot_1dm;
   cout << "CASPT2 : < 0 | F | 0 > = " << E_FOCK << endl;

}

int CheMPS2::CASPT2::vector_helper(){

   int * helper = new int[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];

   /*** Type A : c_tiuv E_ti E_uv | 0 >
                 c_tiuv = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + count_tuv + size_A[ irrep ] * count_i ]
        Type C : c_atuv E_at E_uv | 0 >
                 c_atuv = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + count_tuv + size_C[ irrep ] * count_a ]

        1/ (TYPE A) count_i = 0 .. NOCC[ irrep ]
           (TYPE C) count_a = 0 .. NVIRT[ irrep ]
        2/ jump_tuv = 0
           irrep_t = 0 .. num_irreps
              irrep_u = 0 .. num_irreps
                 irrep_v = irrep x irrep_t x irrep_u
                    ---> count_tuv = jump_tuv + t + NDMRG[ irrep_t ] * ( u + NDMRG[ irrep_u ] * v )
                    jump_tuv += NDMRG[ irrep_t ] * NDMRG[ irrep_u ] * NDMRG[ irrep_v ]
   */

   size_A = new int[ num_irreps ];
   size_C = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int linsize_AC = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         for ( int irrep_u = 0; irrep_u < num_irreps; irrep_u++ ){
            const int irrep_v = Irreps::directProd( Irreps::directProd( irrep, irrep_t ), irrep_u );
            linsize_AC += indices->getNDMRG( irrep_t ) * indices->getNDMRG( irrep_u ) * indices->getNDMRG( irrep_v );
         }
      }
      size_A[ irrep ] = linsize_AC;
      size_C[ irrep ] = linsize_AC;
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_A ] = size_A[ irrep ] * indices->getNOCC( irrep );
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_C ] = size_C[ irrep ] * indices->getNVIRT( irrep );
   }

   /*** Type D1 : c1_aitu E_ai E_tu | 0 >
        Type D2 : c2_tiau E_ti E_au | 0 >

                  c1_aitu = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + count_tu +          size_D[ irrep ] * count_ai ]
                  c2_tiau = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + count_tu + D2JUMP + size_D[ irrep ] * count_ai ]

                  D2JUMP = size_D[ irrep ] / 2

        1/ jump_ai = 0
           irrep_i = 0 .. num_irreps
              irrep_a = irrep_i x irrep
              ---> count_ai = jump_ai + i + NOCC[ irrep_i ] * a
              jump_ai += NOCC[ irrep_i ] * NVIRT[ irrep_a ]
        2/ jump_tu = 0
           irrep_t = 0 .. num_irreps
              irrep_u = irrep x irrep_t
              ---> count_tu = jump_tu + t + NDMRG[ irrep_t ] * u
              jump_tu += NDMRG[ irrep_t ] * NDMRG[ irrep_u ]
   */

   size_D = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int jump_tu = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         const int nact_t = indices->getNDMRG( irrep_t );
         const int nact_u = indices->getNDMRG( irrep_u );
         jump_tu += nact_t * nact_u;
      }
      size_D[ irrep ] = 2 * jump_tu;
      int jump_ai = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_a = Irreps::directProd( irrep_i, irrep );
         const int nocc_i  = indices->getNOCC( irrep_i );
         const int nvirt_a = indices->getNVIRT( irrep_a );
         jump_ai += nocc_i * nvirt_a;
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_D ] = jump_ai * size_D[ irrep ];
   }

   /*** Type B singlet : c_tiuj ( E_ti E_uj + E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 > with i <= j and t <= u
                         c_tiuj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + count_tu + size_B_singlet[ irrep ] * count_ij ]
        Type B triplet : c_tiuj ( E_ti E_uj - E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 > with i <  j and t <  u
                         c_tiuj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + count_tu + size_B_triplet[ irrep ] * count_ij ]
        Type F singlet : c_atbu ( E_at E_bu + E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 > with a <= b and t <= u
                         c_atbu = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + count_tu + size_F_singlet[ irrep ] * count_ab ]
        Type F triplet : c_atbu ( E_at E_bu - E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 > with a <  b and t <  u
                         c_atbu = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + count_tu + size_F_triplet[ irrep ] * count_ab ]

        1/ jump_ab = 0
           if ( irrep == 0 ):
              irrep_ab = 0 .. num_irreps
                 (SINGLET) ---> count_ab = jump_ab + a + b(b+1)/2
                 (SINGLET) jump_ab += NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] + 1 ) / 2
                 (TRIPLET) ---> count_ab = jump_ab + a + b(b-1)/2
                 (TRIPLET) jump_ab += NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] - 1 ) / 2
           else:
              irrep_a = 0 .. num_irreps
                 irrep_b = irrep x irrep_a
                 if ( irrep_a < irrep_b ):
                    ---> count_ab = jump_ab + a + NVIRT[ irrep_a ] * b
                    jump_ab += NVIRT[ irrep_i ] * NVIRT[ irrep_j ]

        2/ jump_ij = 0
           if irrep == 0:
              irrep_ij = 0 .. num_irreps
                 (SINGLET) ---> count_ij = jump_ij + i + j(j+1)/2
                 (SINGLET) jump_ij += NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] + 1 ) / 2
                 (TRIPLET) ---> count_ij = jump_ij + i + j(j-1)/2
                 (TRIPLET) jump_ij += NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] - 1 ) / 2
           else:
              irrep_i = 0 .. num_irreps
                 irrep_j = irrep x irrep_i
                 if ( irrep_i < irrep_j ):
                    ---> count_ij = jump_ij + i + NOCC[ irrep_i ] * j
                    jump_ij += NOCC[ irrep_i ] * NOCC[ irrep_j ]

        3/ jump_tu = 0
           if irrep == 0:
              irrep_tu = 0 .. num_irreps
                 (SINGLET) ---> count_tu = jump_tu + t + u(u+1)/2
                 (SINGLET) jump_tu += NDMRG[ irrep_tu ] * ( NDMRG[ irrep_tu ] + 1 ) / 2
                 (TRIPLET) ---> count_tu = jump_tu + t + u(u-1)/2
                 (TRIPLET) jump_tu += NDMRG[ irrep_tu ] * ( NDMRG[ irrep_tu ] - 1 ) / 2
           else:
              irrep_t = 0 .. num_irreps
                 irrep_u = irrep x irrep_t
                 if ( irrep_t < irrep_u ):
                    ---> count_tu = jump_tu + t + NDMRG[ irrep_t ] * u
                    jump_tu += NDMRG[ irrep_t ] * NDMRG[ irrep_u ]
   */

   size_B_singlet = new int[ num_irreps ];
   size_B_triplet = new int[ num_irreps ];
   size_F_singlet = new int[ num_irreps ];
   size_F_triplet = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int jump_tu_singlet = 0;
      int jump_tu_triplet = 0;
      if ( irrep == 0 ){
         for ( int irrep_tu = 0; irrep_tu < num_irreps; irrep_tu++ ){ // irrep_u == irrep_t
            const int nact_tu = indices->getNDMRG( irrep_tu );
            jump_tu_singlet += ( nact_tu * ( nact_tu + 1 )) / 2;
            jump_tu_triplet += ( nact_tu * ( nact_tu - 1 )) / 2;
         }
      } else {
         for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
            const int irrep_u = Irreps::directProd( irrep, irrep_t );
            if ( irrep_t < irrep_u ){
               const int nact_t = indices->getNDMRG( irrep_t );
               const int nact_u = indices->getNDMRG( irrep_u );
               jump_tu_singlet += nact_t * nact_u;
               jump_tu_triplet += nact_t * nact_u;
            }
         }
      }
      size_B_singlet[ irrep ] = jump_tu_singlet;
      size_B_triplet[ irrep ] = jump_tu_triplet;
      size_F_singlet[ irrep ] = jump_tu_singlet;
      size_F_triplet[ irrep ] = jump_tu_triplet;

      int linsize_B_singlet = 0;
      int linsize_B_triplet = 0;
      int linsize_F_singlet = 0;
      int linsize_F_triplet = 0;
      if ( irrep == 0 ){ // irrep_i == irrep_j    or    irrep_a == irrep_b
         for ( int irrep_ijab = 0; irrep_ijab < num_irreps; irrep_ijab++ ){
            const int nocc_ij = indices->getNOCC( irrep_ijab );
            linsize_B_singlet += ( nocc_ij * ( nocc_ij + 1 ))/2;
            linsize_B_triplet += ( nocc_ij * ( nocc_ij - 1 ))/2;
            const int nvirt_ab = indices->getNVIRT( irrep_ijab );
            linsize_F_singlet += ( nvirt_ab * ( nvirt_ab + 1 ))/2;
            linsize_F_triplet += ( nvirt_ab * ( nvirt_ab - 1 ))/2;
         }
      } else { // irrep_i < irrep_j = irrep_i x irrep    or    irrep_a < irrep_b = irrep_a x irrep
         for ( int irrep_ai = 0; irrep_ai < num_irreps; irrep_ai++ ){
            const int irrep_bj = Irreps::directProd( irrep, irrep_ai );
            if ( irrep_ai < irrep_bj ){
               const int nocc_i = indices->getNOCC( irrep_ai );
               const int nocc_j = indices->getNOCC( irrep_bj );
               linsize_B_singlet += nocc_i * nocc_j;
               linsize_B_triplet += nocc_i * nocc_j;
               const int nvirt_a = indices->getNVIRT( irrep_ai );
               const int nvirt_b = indices->getNVIRT( irrep_bj );
               linsize_F_singlet += nvirt_a * nvirt_b;
               linsize_F_triplet += nvirt_a * nvirt_b;
            }
         }
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] = linsize_B_singlet * size_B_singlet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] = linsize_B_triplet * size_B_triplet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] = linsize_F_singlet * size_F_singlet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] = linsize_F_triplet * size_F_triplet[ irrep ];
   }

   /*** Type E singlet : c_tiaj ( E_ti E_aj + E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 > with i <= j
                         c_tiaj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + count_t + size_E[ irrep ] * count_aij ]
        Type E triplet : c_tiaj ( E_ti E_aj - E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 > with i <  j
                         c_tiaj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + count_t + size_E[ irrep ] * count_aij ]

        1/ jump_aij = 0
           irrep_a = 0 .. num_irreps
              irrep_occ = irrep_a x irrep
              if ( irrep_occ == 0 ):
                 irrep_ij = 0 .. num_irreps
                    (SINGLET) ---> count_aij = jump_aij + a + NVIRT[ irrep_a ] * ( i + j(j+1)/2 )
                    (SINGLET) jump_aij += NVIRT[ irrep_a ] * NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] + 1 ) / 2
                    (TRIPLET) ---> count_aij = jump_aij + a + NVIRT[ irrep_a ] * ( i + j(j-1)/2 )
                    (TRIPLET) jump_aij += NVIRT[ irrep_a ] * NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] - 1 ) / 2
              else:
                 irrep_i = 0 .. num_irreps
                    irrep_j = irrep_occ x irrep_i
                       if ( irrep_i < irrep_j ):
                          ---> count_aij = jump_aij + a + NVIRT[ irrep_a ] * ( i + NOCC[ irrep_i ] * j )
                          jump_aij += NVIRT[ irrep_a ] * NOCC[ irrep_i ] * NOCC[ irrep_j ]
        2/ count_t = 0 .. NDMRG[ irrep ]
   */

   size_E = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      size_E[ irrep ] = indices->getNDMRG( irrep );
      int linsize_E_singlet = 0;
      int linsize_E_triplet = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int nvirt_a = indices->getNVIRT( irrep_a );
         const int irrep_occ = Irreps::directProd( irrep, irrep_a );
         if ( irrep_occ == 0 ){ // irrep_i == irrep_j
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               const int nocc_ij = indices->getNOCC( irrep_ij );
               linsize_E_singlet += ( nvirt_a * nocc_ij * ( nocc_ij + 1 )) / 2;
               linsize_E_triplet += ( nvirt_a * nocc_ij * ( nocc_ij - 1 )) / 2;
            }
         } else { // irrep_i < irrep_j = irrep_i x irrep_occ
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
               if ( irrep_i < irrep_j ){
                  const int nocc_i = indices->getNOCC( irrep_i );
                  const int nocc_j = indices->getNOCC( irrep_j );
                  linsize_E_singlet += nvirt_a * nocc_i * nocc_j;
                  linsize_E_triplet += nvirt_a * nocc_i * nocc_j;
               }
            }
         }
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] = linsize_E_singlet * size_E[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] = linsize_E_triplet * size_E[ irrep ];
   }

   /*** Type G singlet : c_aibt ( E_ai E_bt + E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 > with a <= b
                         c_aibt = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + count_t + size_G[ irrep ] * count_abi ]
        Type G triplet : c_aibt ( E_ai E_bt - E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 > with a <  b
                         c_aibt = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + count_t + size_G[ irrep ] * count_abi ]

        1/ jump_abi = 0
           irrep_i = 0 .. num_irreps
              irrep_virt = irrep_i x irrep
              if irrep_virt == 0:
                 irrep_ab = 0 .. num_irreps
                    (SINGLET) ---> count_abi = jump_abi + i + NOCC[ irrep_i ] * ( a + b(b+1)/2 )
                    (SINGLET) jump_abi += NOCC[ irrep_i ] * NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] + 1 ) / 2
                    (TRIPLET) ---> count_abi = jump_abi + i + NOCC[ irrep_i ] * ( a + b(b-1)/2 )
                    (TRIPLET) jump_abi += NOCC[ irrep_i ] * NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] - 1 ) / 2
              else:
                 irrep_a = 0 .. num_irreps
                    irrep_b = irrep_virt x irrep_a
                       if ( irrep_a < irrep_b ):
                          ---> count_abi = jump_abi + i + NOCC[ irrep_i ] * ( a + NVIRT[ irrep_a ] * b )
                          jump_abi += NOCC[ irrep_i ] * NVIRT[ irrep_a ] * NVIRT[ irrep_b ]
        2/ count_t = 0 .. NDMRG[ irrep ]
   */

   size_G = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      size_G[ irrep ] = indices->getNDMRG( irrep );
      int linsize_G_singlet = 0;
      int linsize_G_triplet = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int nocc_i = indices->getNOCC( irrep_i );
         const int irrep_virt = Irreps::directProd( irrep, irrep_i );
         if ( irrep_virt == 0 ){ // irrep_a == irrep_b
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int nvirt_ab = indices->getNVIRT( irrep_ab );
               linsize_G_singlet += ( nocc_i * nvirt_ab * ( nvirt_ab + 1 )) / 2;
               linsize_G_triplet += ( nocc_i * nvirt_ab * ( nvirt_ab - 1 )) / 2;
            }
         } else { // irrep_a < irrep_b = irrep_a x irrep_virt
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_virt, irrep_a );
               if ( irrep_a < irrep_b ){
                  const int nvirt_a = indices->getNVIRT( irrep_a );
                  const int nvirt_b = indices->getNVIRT( irrep_b );
                  linsize_G_singlet += nocc_i * nvirt_a * nvirt_b;
                  linsize_G_triplet += nocc_i * nvirt_a * nvirt_b;
               }
            }
         }
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] = linsize_G_singlet * size_G[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] = linsize_G_triplet * size_G[ irrep ];
   }

   /*** Type H singlet : c_aibj ( E_ai E_bj + E_bi E_aj ) / sqrt( ( 1 + delta_ij ) * ( 1 + delta_ab ) ) | 0 > with a <= b and i <= j
                         c_aibj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + count_aibj ]
        Type H triplet : c_aibj ( E_ai E_bj - E_bi E_aj ) / sqrt( ( 1 + delta_ij ) * ( 1 + delta_ab ) ) | 0 > with a <  b and i <  j
                         c_aibj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + count_aibj ]

        1/ irrep = 0 .. num_irreps
              if ( irrep == 0 ):
                 jump_aibj = 0
                 irrep_ij = 0 .. num_irreps
                    (SINGLET) linsize_ij = NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] + 1 ) / 2
                    (TRIPLET) linsize_ij = NOCC[ irrep_ij ] * ( NOCC[ irrep_ij ] - 1 ) / 2
                    irrep_ab = 0 .. num_irreps
                       (SINGLET) linsize_ab = NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] + 1 ) / 2
                       (TRIPLET) linsize_ab = NVIRT[ irrep_ab ] * ( NVIRT[ irrep_ab ] - 1 ) / 2
                       (SINGLET) ---> count_aibj = jump_aibj + i + j(j+1)/2 + linsize_ij * ( a + b(b+1)/2 )
                       (TRIPLET) ---> count_aibj = jump_aibj + i + j(j-1)/2 + linsize_ij * ( a + b(b-1)/2 )
                       jump_aibj += linsize_ij * linsize_ab
              else:
                 jump_aibj = 0
                 irrep_i = 0 .. num_irreps
                    irrep_j = irrep x irrep_i
                    if ( irrep_i < irrep_j ):
                       irrep_a = 0 .. num_irreps
                          irrep_b = irrep x irrep_a
                          if ( irrep_a < irrep_b ):
                             ---> count_aibj = jump_aibj + i + NOCC[ irrep_i ] * ( j + NOCC[ irrep_j ] * ( a + NVIRT[ irrep_a ] * b ) )
                             jump_aibj += NOCC[ irrep_i ] * NOCC[ irrep_j ] * NVIRT[ irrep_a ] * NVIRT[ irrep_b ]
   */

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int linsize_H_singlet = 0;
      int linsize_H_triplet = 0;
      if ( irrep == 0 ){ // irrep_i == irrep_j  and  irrep_a == irrep_b
         int linsize_ij_singlet = 0;
         int linsize_ij_triplet = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            linsize_ij_singlet += ( nocc_ij * ( nocc_ij + 1 )) / 2;
            linsize_ij_triplet += ( nocc_ij * ( nocc_ij - 1 )) / 2;
         }
         int linsize_ab_singlet = 0;
         int linsize_ab_triplet = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int nvirt_ab = indices->getNVIRT( irrep_ab );
            linsize_ab_singlet += ( nvirt_ab * ( nvirt_ab + 1 )) / 2;
            linsize_ab_triplet += ( nvirt_ab * ( nvirt_ab - 1 )) / 2;
         }
         linsize_H_singlet = linsize_ij_singlet * linsize_ab_singlet;
         linsize_H_triplet = linsize_ij_triplet * linsize_ab_triplet;
      } else { // irrep_i < irrep_j = irrep_i x irrep   and   irrep_a < irrep_b = irrep_a x irrep
         int linsize_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){ linsize_ij += indices->getNOCC( irrep_i ) * indices->getNOCC( irrep_j ); }
         }
         int linsize_ab = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int irrep_b = Irreps::directProd( irrep, irrep_a );
            if ( irrep_a < irrep_b ){ linsize_ab += indices->getNVIRT( irrep_a ) * indices->getNVIRT( irrep_b ); }
         }
         linsize_H_singlet = linsize_ij * linsize_ab;
         linsize_H_triplet = linsize_ij * linsize_ab;
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] = linsize_H_singlet;
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] = linsize_H_triplet;
   }

   jump = new int[ CHEMPS2_CASPT2_NUM_CASES * num_irreps + 1 ];
   jump[ 0 ] = 0;
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES * num_irreps; cnt++ ){ jump[ cnt+1 ] = jump[ cnt ] + helper[ cnt ]; }
   delete [] helper;
   const int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
   assert( total_size == debug_total_length() );
   cout << "CASPT2 :     Original size of the V_SD space = " << total_size << endl;
   return total_size;

}

long long CheMPS2::CASPT2::debug_total_length() const{

   long long length = 0;
   for ( int i1 = 0; i1 < num_irreps; i1++ ){
      const long long nocc1 = indices->getNOCC( i1 );
      const long long nact1 = indices->getNDMRG( i1 );
      const long long nvir1 = indices->getNVIRT( i1 );
      for ( int i2 = 0; i2 < num_irreps; i2++ ){
         const long long nocc2 = indices->getNOCC( i2 );
         const long long nact2 = indices->getNDMRG( i2 );
         for ( int i3 = 0; i3 < num_irreps; i3++ ){
            const int i4 = Irreps::directProd( Irreps::directProd( i1, i2 ), i3 );
            const long long nact3 = indices->getNDMRG( i3 );
            const long long nvir3 = indices->getNVIRT( i3 );
            const long long nocc4 = indices->getNOCC( i4 );
            const long long nact4 = indices->getNDMRG( i4 );
            const long long nvir4 = indices->getNVIRT( i4 );

            length +=     nocc1 * nact2 * nact3 * nact4; // A:  E_ti E_uv | 0 >
            length +=     nact1 * nact2 * nact3 * nvir4; // C:  E_at E_uv | 0 >
            length += 2 * nocc1 * nact2 * nact3 * nvir4; // D:  E_ai E_tu | 0 >  and  E_ti E_au | 0 >
            length +=     nocc1 * nocc2 * nact3 * nvir4; // E:  E_ti E_aj | 0 >
            length +=     nocc1 * nact2 * nvir3 * nvir4; // G:  E_ai E_bt | 0 >
            if ( i2 < i4 ){ // 2 < 4 and irrep_2 < irrep_4
               length += nact1 * nact3 * nocc2 * nocc4;  // B:  E_ti E_uj | 0 >
               length += nvir1 * nvir3 * nact2 * nact4;  // F:  E_at E_bu | 0 >
               length += nvir1 * nvir3 * nocc2 * nocc4;  // H:  E_ai E_bj | 0 >
            }
            if ( i2 == i4 ){ // i2 == i4 implies i1 == i3
               // 2 < 4 and irrep_2 == irrep_4
               length += ( nact1 * nact3 * nocc2 * ( nocc2 - 1 ) ) / 2; // B:  E_ti E_uj | 0 >
               length += ( nvir1 * nvir3 * nact2 * ( nact2 - 1 ) ) / 2; // F:  E_at E_bu | 0 >
               length += ( nvir1 * nvir3 * nocc2 * ( nocc2 - 1 ) ) / 2; // H:  E_ai E_bj | 0 >
               // 2 == 4 and 1 <= 3
               length += ( nact1 * ( nact3 + 1 ) * nocc2 ) / 2; // B:  E_ti E_uj | 0 >
               length += ( nvir1 * ( nvir3 + 1 ) * nact2 ) / 2; // F:  E_at E_bu | 0 >
               length += ( nvir1 * ( nvir3 + 1 ) * nocc2 ) / 2; // H:  E_ai E_bj | 0 >
            }
         }
      }
   }

   return length;

}

int CheMPS2::CASPT2::recreatehelper1( double * FOCK, double * OVLP, int SIZE, double * work, double * eigs, int lwork ){

   if ( SIZE == 0 ){ return SIZE; }

   // S = U_S eigs_S U_S^T
   char jobz = 'V';
   char uplo = 'U';
   int info;
   dsyev_( &jobz, &uplo, &SIZE, OVLP, &SIZE, eigs, work, &lwork, &info ); // eigs in ascending order

   // Discard smallest eigenvalues
   int skip = 0;
   bool ctu = true;
   while (( skip < SIZE ) && ( ctu )){
      if ( eigs[ skip ] < CheMPS2::CASPT2_OVLP_CUTOFF ){
         skip++;
      } else {
         ctu = false;
      }
   }
   int NEWSIZE = SIZE - skip;

   // OVLP  <---  U_S eigs_S^{-0.5}
   for ( int col = skip; col < SIZE; col++ ){
      const double prefactor = 1.0 / sqrt( eigs[ col ] );
      for ( int row = 0; row < SIZE; row++ ){
         OVLP[ row + SIZE * col ] *= prefactor;
      }
   }

   // FOCK  <---  FOCK_tilde = eigs_S^{-0.5} U_S^T FOCK U_S eigs_S^{-0.5}
   char notrans = 'N';
   char trans   = 'T';
   double one   = 1.0;
   double set   = 0.0;
   dgemm_( &notrans, &notrans, &SIZE,    &NEWSIZE, &SIZE, &one, FOCK, &SIZE, OVLP + skip * SIZE, &SIZE, &set, work, &SIZE    ); // work = FOCK * V * eigs^{-0.5}
   dgemm_( &trans,   &notrans, &NEWSIZE, &NEWSIZE, &SIZE, &one, OVLP + skip * SIZE, &SIZE, work, &SIZE, &set, FOCK, &NEWSIZE ); // FOCK = ( V * eigs^{-0.5} )^T * work

   // FOCK_tilde = U_F_tilde eigs_F_tilde U_F_tilde^T
   dsyev_( &jobz, &uplo, &NEWSIZE, FOCK, &NEWSIZE, eigs, work, &lwork, &info ); // eigs in ascending order

   // OVLP  <---  U_S eigs_S^{-0.5} U_F_tilde
   dgemm_( &notrans, &notrans, &SIZE, &NEWSIZE, &NEWSIZE, &one, OVLP + skip * SIZE, &SIZE, FOCK, &NEWSIZE, &set, work, &SIZE );
   int size_copy = SIZE * NEWSIZE;
   int inc1 = 1;
   dcopy_( &size_copy, work, &inc1, OVLP, &inc1 );

   // FOCK  <---  eigs_F_tilde
   dcopy_( &NEWSIZE, eigs, &inc1, FOCK, &inc1 );

   return NEWSIZE;

}

void CheMPS2::CASPT2::recreatehelper2( double * LEFT, double * RIGHT, double ** matrix, double * work, int OLD_LEFT, int NEW_LEFT, int OLD_RIGHT, int NEW_RIGHT, const int number ){

   double set = 0.0;
   double one = 1.0;
   char trans   = 'T';
   char notrans = 'N';
   for ( int count = 0; count < number; count++ ){
      // work  <---  LEFT[ old_left, new_left ]^T * matrix[ count ][ old_left, old_right ]
      dgemm_( &trans, &notrans, &NEW_LEFT, &OLD_RIGHT, &OLD_LEFT, &one, LEFT, &OLD_LEFT, matrix[ count ], &OLD_LEFT, &set, work, &NEW_LEFT );
      // matrix[ count ]  <---  work[ new_left, old_right ] * RIGHT[ old_right, new_right ]
      dgemm_( &notrans, &notrans, &NEW_LEFT, &NEW_RIGHT, &OLD_RIGHT, &one, work, &NEW_LEFT, RIGHT, &OLD_RIGHT, &set, matrix[ count ], &NEW_LEFT );
   }

}

void CheMPS2::CASPT2::recreatehelper3( double * OVLP, int OLDSIZE, int NEWSIZE, double * rhs_old, double * rhs_new, const int num_rhs ){

   // rhs <-- eigs^{-0.5} V^T rhs
   int inc1 = 1;
   double set = 0.0;
   double one = 1.0;
   char trans = 'T';
   for ( int sector = 0; sector < num_rhs; sector++ ){
      dgemv_( &trans, &OLDSIZE, &NEWSIZE, &one, OVLP, &OLDSIZE, rhs_old + OLDSIZE * sector, &inc1, &set, rhs_new + NEWSIZE * sector, &inc1 );
   }

}

void CheMPS2::CASPT2::recreate(){

   const int maxsize = get_maxsize(); // Minimum 3
   const int lwork = maxsize * maxsize;
   double * work = new double[ lwork ];
   double * eigs = new double[ maxsize ];

   int * newsize_A = new int[ num_irreps ];
   int * newsize_C = new int[ num_irreps ];
   int * newsize_D = new int[ num_irreps ];
   int * newsize_E = new int[ num_irreps ];
   int * newsize_G = new int[ num_irreps ];
   int * newsize_B_singlet = new int[ num_irreps ];
   int * newsize_B_triplet = new int[ num_irreps ];
   int * newsize_F_singlet = new int[ num_irreps ];
   int * newsize_F_triplet = new int[ num_irreps ];

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      newsize_A[ irrep ] = recreatehelper1( FAA[ irrep ], SAA[ irrep ], size_A[ irrep ], work, eigs, lwork );
      newsize_C[ irrep ] = recreatehelper1( FCC[ irrep ], SCC[ irrep ], size_C[ irrep ], work, eigs, lwork );
      newsize_D[ irrep ] = recreatehelper1( FDD[ irrep ], SDD[ irrep ], size_D[ irrep ], work, eigs, lwork );
      newsize_E[ irrep ] = recreatehelper1( FEE[ irrep ], SEE[ irrep ], size_E[ irrep ], work, eigs, lwork );
      newsize_G[ irrep ] = recreatehelper1( FGG[ irrep ], SGG[ irrep ], size_G[ irrep ], work, eigs, lwork );
      newsize_B_singlet[ irrep ] = recreatehelper1( FBB_singlet[ irrep ], SBB_singlet[ irrep ], size_B_singlet[ irrep ], work, eigs, lwork );
      newsize_B_triplet[ irrep ] = recreatehelper1( FBB_triplet[ irrep ], SBB_triplet[ irrep ], size_B_triplet[ irrep ], work, eigs, lwork );
      newsize_F_singlet[ irrep ] = recreatehelper1( FFF_singlet[ irrep ], SFF_singlet[ irrep ], size_F_singlet[ irrep ], work, eigs, lwork );
      newsize_F_triplet[ irrep ] = recreatehelper1( FFF_triplet[ irrep ], SFF_triplet[ irrep ], size_F_triplet[ irrep ], work, eigs, lwork );
   }

   for ( int IL = 0; IL < num_irreps; IL++ ){
      for ( int IR = 0; IR < num_irreps; IR++ ){

         const int irrep_w = Irreps::directProd( IL, IR );
         const int num_w   = indices->getNDMRG( irrep_w );

         if ( newsize_A[ IL ] * newsize_D[ IR ] * num_w > 0 ){
            recreatehelper2( SAA[ IL ], SDD[ IR ], FAD[ IL ][ IR ], work, size_A[ IL ], newsize_A[ IL ], size_D[ IR ], newsize_D[ IR ], num_w );
         }

         if ( newsize_C[ IL ] * newsize_D[ IR ] * num_w > 0 ){
            recreatehelper2( SCC[ IL ], SDD[ IR ], FCD[ IL ][ IR ], work, size_C[ IL ], newsize_C[ IL ], size_D[ IR ], newsize_D[ IR ], num_w );
         }

         if ( newsize_A[ IL ] * newsize_B_singlet[ IR ] * num_w > 0 ){
            recreatehelper2( SAA[ IL ], SBB_singlet[ IR ], FAB_singlet[ IL ][ IR ], work, size_A[ IL ], newsize_A[ IL ], size_B_singlet[ IR ], newsize_B_singlet[ IR ], num_w );
         }

         if ( newsize_A[ IL ] * newsize_B_triplet[ IR ] * num_w > 0 ){
            recreatehelper2( SAA[ IL ], SBB_triplet[ IR ], FAB_triplet[ IL ][ IR ], work, size_A[ IL ], newsize_A[ IL ], size_B_triplet[ IR ], newsize_B_triplet[ IR ], num_w );
         }

         if ( newsize_C[ IL ] * newsize_F_singlet[ IR ] * num_w > 0 ){
            recreatehelper2( SCC[ IL ], SFF_singlet[ IR ], FCF_singlet[ IL ][ IR ], work, size_C[ IL ], newsize_C[ IL ], size_F_singlet[ IR ], newsize_F_singlet[ IR ], num_w );
         }

         if ( newsize_C[ IL ] * newsize_F_triplet[ IR ] * num_w > 0 ){
            recreatehelper2( SCC[ IL ], SFF_triplet[ IR ], FCF_triplet[ IL ][ IR ], work, size_C[ IL ], newsize_C[ IL ], size_F_triplet[ IR ], newsize_F_triplet[ IR ], num_w );
         }

         if ( newsize_B_singlet[ IL ] * newsize_E[ IR ] * num_w > 0 ){
            recreatehelper2( SBB_singlet[ IL ], SEE[ IR ], FBE_singlet[ IL ][ IR ], work, size_B_singlet[ IL ], newsize_B_singlet[ IL ], size_E[ IR ], newsize_E[ IR ], num_w );
         }

         if ( newsize_B_triplet[ IL ] * newsize_E[ IR ] * num_w > 0 ){
            recreatehelper2( SBB_triplet[ IL ], SEE[ IR ], FBE_triplet[ IL ][ IR ], work, size_B_triplet[ IL ], newsize_B_triplet[ IL ], size_E[ IR ], newsize_E[ IR ], num_w );
         }

         if ( newsize_F_singlet[ IL ] * newsize_G[ IR ] * num_w > 0 ){
            recreatehelper2( SFF_singlet[ IL ], SGG[ IR ], FFG_singlet[ IL ][ IR ], work, size_F_singlet[ IL ], newsize_F_singlet[ IL ], size_G[ IR ], newsize_G[ IR ], num_w );
         }

         if ( newsize_F_triplet[ IL ] * newsize_G[ IR ] * num_w > 0 ){
            recreatehelper2( SFF_triplet[ IL ], SGG[ IR ], FFG_triplet[ IL ][ IR ], work, size_F_triplet[ IL ], newsize_F_triplet[ IL ], size_G[ IR ], newsize_G[ IR ], num_w );
         }

         if ( newsize_D[ IL ] * newsize_E[ IR ] * num_w > 0 ){
            recreatehelper2( SDD[ IL ], SEE[ IR ], FDE_singlet[ IL ][ IR ], work, size_D[ IL ], newsize_D[ IL ], size_E[ IR ], newsize_E[ IR ], num_w );
            recreatehelper2( SDD[ IL ], SEE[ IR ], FDE_triplet[ IL ][ IR ], work, size_D[ IL ], newsize_D[ IL ], size_E[ IR ], newsize_E[ IR ], num_w );
         }

         if ( newsize_D[ IL ] * newsize_G[ IR ] * num_w > 0 ){
            recreatehelper2( SDD[ IL ], SGG[ IR ], FDG_singlet[ IL ][ IR ], work, size_D[ IL ], newsize_D[ IL ], size_G[ IR ], newsize_G[ IR ], num_w );
            recreatehelper2( SDD[ IL ], SGG[ IR ], FDG_triplet[ IL ][ IR ], work, size_D[ IL ], newsize_D[ IL ], size_G[ IR ], newsize_G[ IR ], num_w );
         }

         if ( IR == 0 ){ // IL == irrep_w
            if ( newsize_E[ IL ] * num_w > 0 ){
               double one = 1.0;
               recreatehelper2( SEE[ IL ], &one, FEH[ IL ], work, size_E[ IL ], newsize_E[ IL ], 1, 1, num_w );
            }

            if ( newsize_G[ IL ] * num_w > 0 ){
               double one = 1.0;
               recreatehelper2( SGG[ IL ], &one, FGH[ IL ], work, size_G[ IL ], newsize_G[ IL ], 1, 1, num_w );
            }
         }
      }
   }

   delete [] work;
   delete [] eigs;

   double * tempvector_rhs = new double[ jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ] ];
   int * helper = new int[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ];
   for ( int ptr = 0; ptr < num_irreps * CHEMPS2_CASPT2_NUM_CASES; ptr++ ){ helper[ ptr ] = 0; }

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      if ( newsize_A[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_A;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_A[ irrep ];
         assert( num_rhs * size_A[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_A[ irrep ];
         recreatehelper3( SAA[ irrep ], size_A[ irrep ], newsize_A[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_B_singlet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_B_singlet[ irrep ];
         assert( num_rhs * size_B_singlet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_B_singlet[ irrep ];
         recreatehelper3( SBB_singlet[ irrep ], size_B_singlet[ irrep ], newsize_B_singlet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_B_triplet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_B_triplet[ irrep ];
         assert( num_rhs * size_B_triplet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_B_triplet[ irrep ];
         recreatehelper3( SBB_triplet[ irrep ], size_B_triplet[ irrep ], newsize_B_triplet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_C[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_C;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_C[ irrep ];
         assert( num_rhs * size_C[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_C[ irrep ];
         recreatehelper3( SCC[ irrep ], size_C[ irrep ], newsize_C[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_D[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_D;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_D[ irrep ];
         assert( num_rhs * size_D[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_D[ irrep ];
         recreatehelper3( SDD[ irrep ], size_D[ irrep ], newsize_D[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_E[ irrep ] > 0 ){
         const int ptr1 = irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET;
         const int num_rhs1 = ( jump[ ptr1 + 1 ] - jump[ ptr1 ] ) / size_E[ irrep ];
         assert( num_rhs1 * size_E[ irrep ] == jump[ ptr1 + 1 ] - jump[ ptr1 ] );
         helper[ ptr1 ] = num_rhs1 * newsize_E[ irrep ];
         recreatehelper3( SEE[ irrep ], size_E[ irrep ], newsize_E[ irrep ], vector_rhs + jump[ ptr1 ], tempvector_rhs + jump[ ptr1 ], num_rhs1 );
         const int ptr2 = irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET;
         const int num_rhs2 = ( jump[ ptr2 + 1 ] - jump[ ptr2 ] ) / size_E[ irrep ];
         assert( num_rhs2 * size_E[ irrep ] == jump[ ptr2 + 1 ] - jump[ ptr2 ] );
         helper[ ptr2 ] = num_rhs2 * newsize_E[ irrep ];
         recreatehelper3( SEE[ irrep ], size_E[ irrep ], newsize_E[ irrep ], vector_rhs + jump[ ptr2 ], tempvector_rhs + jump[ ptr2 ], num_rhs2 );
      }

      if ( newsize_F_singlet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_F_singlet[ irrep ];
         assert( num_rhs * size_F_singlet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_F_singlet[ irrep ];
         recreatehelper3( SFF_singlet[ irrep ], size_F_singlet[ irrep ], newsize_F_singlet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_F_triplet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_F_triplet[ irrep ];
         assert( num_rhs * size_F_triplet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_F_triplet[ irrep ];
         recreatehelper3( SFF_triplet[ irrep ], size_F_triplet[ irrep ], newsize_F_triplet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_G[ irrep ] > 0 ){
         const int ptr1 = irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET;
         const int num_rhs1 = ( jump[ ptr1 + 1 ] - jump[ ptr1 ] ) / size_G[ irrep ];
         assert( num_rhs1 * size_G[ irrep ] == jump[ ptr1 + 1 ] - jump[ ptr1 ] );
         helper[ ptr1 ] = num_rhs1 * newsize_G[ irrep ];
         recreatehelper3( SGG[ irrep ], size_G[ irrep ], newsize_G[ irrep ], vector_rhs + jump[ ptr1 ], tempvector_rhs + jump[ ptr1 ], num_rhs1 );
         const int ptr2 = irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET;
         const int num_rhs2 = ( jump[ ptr2 + 1 ] - jump[ ptr2 ] ) / size_G[ irrep ];
         assert( num_rhs2 * size_G[ irrep ] == jump[ ptr2 + 1 ] - jump[ ptr2 ] );
         helper[ ptr2 ] = num_rhs2 * newsize_G[ irrep ];
         recreatehelper3( SGG[ irrep ], size_G[ irrep ], newsize_G[ irrep ], vector_rhs + jump[ ptr2 ], tempvector_rhs + jump[ ptr2 ], num_rhs2 );
      }

      const int ptr1 = irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET;
      helper[ ptr1 ] = jump[ ptr1 + 1 ] - jump[ ptr1 ];
      int inc1 = 1;
      dcopy_( helper + ptr1, vector_rhs + jump[ ptr1 ], &inc1, tempvector_rhs + jump[ ptr1 ], &inc1 );
      const int ptr2 = irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET;
      helper[ ptr2 ] = jump[ ptr2 + 1 ] - jump[ ptr2 ];
      dcopy_( helper + ptr2, vector_rhs + jump[ ptr2 ], &inc1, tempvector_rhs + jump[ ptr2 ], &inc1 );

   }

   delete [] vector_rhs;
   delete [] size_A;         size_A = newsize_A;
   delete [] size_C;         size_C = newsize_C;
   delete [] size_D;         size_D = newsize_D;
   delete [] size_E;         size_E = newsize_E;
   delete [] size_G;         size_G = newsize_G;
   delete [] size_B_singlet; size_B_singlet = newsize_B_singlet;
   delete [] size_B_triplet; size_B_triplet = newsize_B_triplet;
   delete [] size_F_singlet; size_F_singlet = newsize_F_singlet;
   delete [] size_F_triplet; size_F_triplet = newsize_F_triplet;

   int * newjump = new int[ num_irreps * CHEMPS2_CASPT2_NUM_CASES + 1 ];
   newjump[ 0 ] = 0;
   for ( int cnt = 0; cnt < num_irreps * CHEMPS2_CASPT2_NUM_CASES; cnt++ ){
      newjump[ cnt + 1 ] = newjump[ cnt ] + helper[ cnt ];
   }
   vector_rhs = new double[ newjump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ] ];
   int inc1 = 1;
   for ( int cnt = 0; cnt < num_irreps * CHEMPS2_CASPT2_NUM_CASES; cnt++ ){
      dcopy_( helper + cnt, tempvector_rhs + jump[ cnt ], &inc1, vector_rhs + newjump[ cnt ], &inc1 );
   }
   delete [] tempvector_rhs;
   delete [] helper;
   delete [] jump;
   jump = newjump;

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      delete [] SAA[ irrep ];
      delete [] SCC[ irrep ];
      delete [] SDD[ irrep ];
      delete [] SEE[ irrep ];
      delete [] SGG[ irrep ];
      delete [] SBB_singlet[ irrep ];
      delete [] SBB_triplet[ irrep ];
      delete [] SFF_singlet[ irrep ];
      delete [] SFF_triplet[ irrep ];
   }
   delete [] SAA;
   delete [] SCC;
   delete [] SDD;
   delete [] SEE;
   delete [] SGG;
   delete [] SBB_singlet;
   delete [] SBB_triplet;
   delete [] SFF_singlet;
   delete [] SFF_triplet;

   const int total_size = jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ];
   cout << "CASPT2 : Nonredundant size of the V_SD space = " << total_size << endl;

}

int CheMPS2::CASPT2::get_maxsize() const{

   int maxsize = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      maxsize = max( max( max( max( max( max( max( max( max(
                                        size_A[irrep],
                                        size_C[irrep] ),
                                        size_D[irrep] ),
                                        size_E[irrep] ),
                                        size_G[irrep] ),
                                        size_B_singlet[irrep] ),
                                        size_B_triplet[irrep] ),
                                        size_F_singlet[irrep] ),
                                        size_F_triplet[irrep] ),
                                        maxsize );
   }
   if ( maxsize <= 2 ){ maxsize = 3; }
   return maxsize;

}

void CheMPS2::CASPT2::matvec( double * vector, double * result, double * diag_fock ) const{

   /*
         FOCK  | A  Bsinglet  Btriplet  C     D1     D2    Esinglet  Etriplet  Fsinglet  Ftriplet  Gsinglet  Gtriplet  Hsinglet  Htriplet
      ---------+-------------------------------------------------------------------------------------------------------------------------
      A        | OK    OK      OK       0     OK     OK    GRAD      GRAD      0         0         0         0         0         0
      Bsinglet | OK    OK      0        0     0      0     OK        0         0         0         0         0         0         0
      Btriplet | OK    0       OK       0     0      0     0         OK        0         0         0         0         0         0
      C        | 0     0       0        OK    OK     OK    0         0         OK        OK        GRAD      GRAD      0         0
      D1       | OK    0       0        OK    OK     OK    OK        OK        0         0         OK        OK        GRAD      GRAD
      D2       | OK    0       0        OK    OK     OK    OK        OK        0         0         OK        OK        GRAD      GRAD
      Esinglet | GRAD  OK      0        0     OK     OK    OK        0         0         0         0         0         OK        0
      Etriplet | GRAD  0       OK       0     OK     OK    0         OK        0         0         0         0         0         OK
      Fsinglet | 0     0       0        OK    0      0     0         0         OK        0         OK        0         0         0
      Ftriplet | 0     0       0        OK    0      0     0         0         0         OK        0         OK        0         0
      Gsinglet | 0     0       0        GRAD  OK     OK    0         0         OK        0         OK        0         OK        0
      Gtriplet | 0     0       0        GRAD  OK     OK    0         0         0         OK        0         OK        0         OK
      Hsinglet | 0     0       0        0     GRAD   GRAD  OK        0         0         0         OK        0         OK        0
      Htriplet | 0     0       0        0     GRAD   GRAD  0         OK        0         0         0         OK        0         OK
      
   */

   const int vectorlength = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
   for ( int elem = 0; elem < vectorlength; elem++ ){ result[ elem ] = diag_fock[ elem ] * vector[ elem ]; }

   const int maxlinsize = get_maxsize();
   double * workspace = new double[ maxlinsize * maxlinsize ];
   const double SQRT2 = sqrt( 2.0 );

   // FAD: < A(xjyz) E_wc D(aitu) > = delta_ac delta_ij FAD[ Ij ][ Ii x Ia ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ii == Ij == Ix x Iy x Iz
      int SIZE_L = size_A[ IL ];
      const int nocc_ij = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ia
         int SIZE_R = size_D[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Ia == Ic == Iw == IL x IR
         const int shift = shift_D_nonactive( indices, IL, Iw );
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         const int n_oa_w = nocc_w + nact_w;
         const int nvir_w = indices->getNVIRT( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nocc_ij * nact_w * nvir_w > 0 ){
            for ( int ac = 0; ac < nvir_w; ac++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_wc = fock->get( Iw, nocc_w + w, n_oa_w + ac );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_wc, FAD[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ij = 0; ij < nocc_ij; ij++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A ] + SIZE_L * ij;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + ij + nocc_ij * ac );
                  double one = 1.0;
                  char notrans = 'N';
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ij = 0; ij < nocc_ij; ij++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A ] + SIZE_L * ij;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + ij + nocc_ij * ac );
                  double one = 1.0;
                  char trans = 'T';
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }

   // FCD: < C(bxyz) E_kw D(aitu) > = delta_ik delta_ab FCD[ Ib ][ Ii x Ia ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ia == Ib
      int SIZE_L = size_C[ IL ];
      const int nvir_ab = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ia
         int SIZE_R = size_D[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Ii == Ik == Iw == IL x IR
         const int shift = shift_D_nonactive( indices, Iw, IL );
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nvir_ab * nact_w * nocc_w > 0 ){
            for ( int ik = 0; ik < nocc_w; ik++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_kw = fock->get( Iw, ik, nocc_w + w );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_kw, FCD[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ab = 0; ab < nvir_ab; ab++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C ] + SIZE_L * ab;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + ik + nocc_w * ab );
                  double one = 1.0;
                  char notrans = 'N';
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ab = 0; ab < nvir_ab; ab++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C ] + SIZE_L * ab;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + ik + nocc_w * ab );
                  double one = 1.0;
                  char trans = 'T';
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }

   // FAB singlet: < A(xlyz) E_kw SB_tiuj > = ( delta_ik delta_jl + delta_jk delta_il ) / sqrt( 1 + delta_ij ) * FAB_singlet[ Il ][ Ii x Ij ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Il == Ix x Iy x Iz
      int SIZE_L = size_A[ IL ];
      const int nocc_l = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ij
         int SIZE_R = size_B_singlet[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Ik
         const int shift = (( Iw < IL ) ? shift_B_nonactive( indices, Iw, IL, +1 ) : shift_B_nonactive( indices, IL, Iw, +1 ));
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nocc_l * nact_w * nocc_w > 0 ){
            for ( int k = 0; k < nocc_w; k++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_kw = fock->get( Iw, k, nocc_w + w );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_kw, FAB_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               if ( IR == 0 ){ // Ii == Ij  and  Ik == Il
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( k < l ) ? ( k + ( l * ( l + 1 ) ) / 2 ) : ( l + ( k * ( k + 1 ) ) / 2 ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     double factor = (( k == l ) ? SQRT2 : 1.0 );
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( k < l ) ? ( k + ( l * ( l + 1 ) ) / 2 ) : ( l + ( k * ( k + 1 ) ) / 2 ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     double factor = (( k == l ) ? SQRT2 : 1.0 );
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               } else {
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( Iw < IL ) ? ( k + nocc_w * l ) : ( l + nocc_l * k ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( Iw < IL ) ? ( k + nocc_w * l ) : ( l + nocc_l * k ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               }
            }
         }
      }
   }

   // FAB triplet: < A(xlyz) E_kw TB_tiuj > = ( delta_ik delta_jl - delta_jk delta_il ) * FAB_triplet[ Il ][ Ii x Ij ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Il == Ix x Iy x Iz
      int SIZE_L = size_A[ IL ];
      const int nocc_l = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ij
         int SIZE_R = size_B_triplet[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Ik
         const int shift = (( Iw < IL ) ? shift_B_nonactive( indices, Iw, IL, -1 ) : shift_B_nonactive( indices, IL, Iw, -1 ));
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nocc_l * nact_w * nocc_w > 0 ){
            for ( int k = 0; k < nocc_w; k++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_kw = fock->get( Iw, k, nocc_w + w );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_kw, FAB_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               if ( IR == 0 ){ // Ii == Ij  and  Ik == Il
                  for ( int l = 0; l < k; l++ ){
                     const int count = shift + l + ( k * ( k - 1 ) ) / 2;
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     double minus_one = -1.0; // ( k > l  --->  - delta_jk delta_il )
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &minus_one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int count = shift + k + ( l * ( l - 1 ) ) / 2;
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0; // ( k < l  --->  + delta_ik delta_jl )
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = 0; l < k; l++ ){
                     const int count = shift + l + ( k * ( k - 1 ) ) / 2;
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     double minus_one = -1.0; // ( k > l  --->  - delta_jk delta_il )
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &minus_one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int count = shift + k + ( l * ( l - 1 ) ) / 2;
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0; // ( k < l  --->  + delta_ik delta_jl )
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               } else {
                  double factor = (( Iw < IL ) ? 1.0 : -1.0 ); // ( k < l  --->  + delta_ik delta_jl ) and ( k > l  --->  - delta_jk delta_il )
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( Iw < IL ) ? ( k + nocc_w * l ) : ( l + nocc_l * k ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int count = shift + (( Iw < IL ) ? ( k + nocc_w * l ) : ( l + nocc_l * k ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               }
            }
         }
      }
   }

   // FCF singlet: < C(dxyz) E_wc SF_atbu > = ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) * FCF_singlet[ Id ][ Ia x Ib ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Id == Ix x Iy x Iz
      int SIZE_L = size_C[ IL ];
      const int nvir_d = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ia x Ib
         int SIZE_R = size_F_singlet[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Ic
         const int shift = (( Iw < IL ) ? shift_F_nonactive( indices, Iw, IL, +1 ) : shift_F_nonactive( indices, IL, Iw, +1 ));
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         const int n_oa_w = nocc_w + nact_w;
         const int nvir_w = indices->getNVIRT( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nvir_d * nact_w * nvir_w > 0 ){
            for ( int c = 0; c < nvir_w; c++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_wc = fock->get( Iw, nocc_w + w, n_oa_w + c );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_wc, FCF_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               if ( IR == 0 ){ // Ia == Ib  and  Ic == Id
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( c < d ) ? ( c + ( d * ( d + 1 ) ) / 2 ) : ( d + ( c * ( c + 1 ) ) / 2 ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     double factor = (( c == d ) ? SQRT2 : 1.0 );
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( c < d ) ? ( c + ( d * ( d + 1 ) ) / 2 ) : ( d + ( c * ( c + 1 ) ) / 2 ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     double factor = (( c == d ) ? SQRT2 : 1.0 );
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               } else {
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( Iw < IL ) ? ( c + nvir_w * d ) : ( d + nvir_d * c ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( Iw < IL ) ? ( c + nvir_w * d ) : ( d + nvir_d * c ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * count;
                     double one = 1.0;
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               }
            }
         }
      }
   }

   // FCF triplet: < C(dxyz) E_wc TF_atbu > = ( delta_ac delta_bd - delta_ad delta_bc ) * FCF_triplet[ Id ][ Ia x Ib ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Id == Ix x Iy x Iz
      int SIZE_L = size_C[ IL ];
      const int nvir_d = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ia x Ib
         int SIZE_R = size_F_triplet[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Ic
         const int shift = (( Iw < IL ) ? shift_F_nonactive( indices, Iw, IL, -1 ) : shift_F_nonactive( indices, IL, Iw, -1 ));
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         const int n_oa_w = nocc_w + nact_w;
         const int nvir_w = indices->getNVIRT( Iw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nvir_d * nact_w * nvir_w > 0 ){
            for ( int c = 0; c < nvir_w; c++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_wc = fock->get( Iw, nocc_w + w, n_oa_w + c );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_wc, FCF_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               if ( IR == 0 ){ // Ia == Ib  and  Ic == Id
                  for ( int d = 0; d < c; d++ ){
                     const int count = shift + d + ( c * ( c - 1 ) ) / 2;
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     double minus_one = -1.0; // ( c > d  --->  - delta_ad delta_bc )
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &minus_one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int count = shift + c + ( d * ( d - 1 ) ) / 2;
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0; // ( c < d  --->  + delta_ac delta_bd )
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = 0; d < c; d++ ){
                     const int count = shift + d + ( c * ( c - 1 ) ) / 2;
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     double minus_one = -1.0; // ( c > d  --->  - delta_ad delta_bc )
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &minus_one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int count = shift + c + ( d * ( d - 1 ) ) / 2;
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0; // ( c < d  --->  + delta_ac delta_bd )
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &one, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               } else {
                  double factor = (( Iw < IL ) ? 1.0 : -1.0 ); // ( c < d  --->  + delta_ac delta_bd ) and ( c > d  --->  - delta_ad delta_bc )
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( Iw < IL ) ? ( c + nvir_w * d ) : ( d + nvir_d * c ));
                     double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     char notrans = 'N';
                     int inc1 = 1;
                     dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int count = shift + (( Iw < IL ) ? ( c + nvir_w * d ) : ( d + nvir_d * c ));
                     double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * count;
                     double one = 1.0;
                     char trans = 'T';
                     int inc1 = 1;
                     dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                  }
               }
            }
         }
      }
   }

   // FBE singlet: < SB_xkyl E_wc SE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_singlet[ Ik x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iik x Ijl == Ix x Iy
      int SIZE_L = size_B_singlet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iik x Ijl x Iac
         int SIZE_R = size_E[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Iac
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         const int n_oa_w = nocc_w + nact_w;
         const int nvir_w = indices->getNVIRT( Iw );
         int linsize = 0;
         for ( int Iik = 0; Iik < num_irreps; Iik++ ){
            const int Ijl = Irreps::directProd( Iik, IL );
            if ( Iik <= Ijl ){
               const int nocc_ik = indices->getNOCC( Iik );
               linsize += (( IL == 0 ) ? ( nocc_ik * ( nocc_ik + 1 ) ) / 2 : nocc_ik * indices->getNOCC( Ijl ) );
            }
         }
         const int size_ij = linsize;
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nact_w * nvir_w * size_ij > 0 ){
            const int shift_E = shift_E_nonactive( indices, Iw, 0, IL, +1 );
            for ( int ac = 0; ac < nvir_w; ac++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_wc = fock->get( Iw, nocc_w + w, n_oa_w + ac );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_wc, FBE_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ij = 0; ij < size_ij; ij++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_L * ij;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE_R * ( shift_E + ac + nvir_w * ij );
                  char notrans = 'N';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ij = 0; ij < size_ij; ij++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_L * ij;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE_R * ( shift_E + ac + nvir_w * ij );
                  char trans = 'T';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }

   // FBE triplet: < TB_xkyl E_wc TE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_triplet[ Ik x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iik x Ijl == Ix x Iy
      int SIZE_L = size_B_triplet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iik x Ijl x Iac
         int SIZE_R = size_E[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Iac
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         const int n_oa_w = nocc_w + nact_w;
         const int nvir_w = indices->getNVIRT( Iw );
         int linsize = 0;
         for ( int Iik = 0; Iik < num_irreps; Iik++ ){
            const int Ijl = Irreps::directProd( Iik, IL );
            if ( Iik <= Ijl ){
               const int nocc_ik = indices->getNOCC( Iik );
               linsize += (( IL == 0 ) ? ( nocc_ik * ( nocc_ik - 1 ) ) / 2 : nocc_ik * indices->getNOCC( Ijl ) );
            }
         }
         const int size_ij = linsize;
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nact_w * nvir_w * size_ij > 0 ){
            const int shift_E = shift_E_nonactive( indices, Iw, 0, IL, -1 );
            for ( int ac = 0; ac < nvir_w; ac++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_wc = fock->get( Iw, nocc_w + w, n_oa_w + ac );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_wc, FBE_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ij = 0; ij < size_ij; ij++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_L * ij;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE_R * ( shift_E + ac + nvir_w * ij );
                  char notrans = 'N';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ij = 0; ij < size_ij; ij++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_L * ij;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE_R * ( shift_E + ac + nvir_w * ij );
                  char trans = 'T';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }

   // FFG singlet: < SF_cxdy E_kw SG_aibt > = 2 delta_ac delta_bd delta_ik FFG_singlet[ Ic x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iac x Ibd == Ix x Iy
      int SIZE_L = size_F_singlet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iac x Ibd x Iik
         int SIZE_R = size_G[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Iik
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         int linsize = 0;
         for ( int Iac = 0; Iac < num_irreps; Iac++ ){
            const int Ibd = Irreps::directProd( Iac, IL );
            if ( Iac <= Ibd ){
               const int nvir_ac = indices->getNVIRT( Iac );
               linsize += (( IL == 0 ) ? ( nvir_ac * ( nvir_ac + 1 ) ) / 2 : nvir_ac * indices->getNVIRT( Ibd ) );
            }
         }
         const int size_ab = linsize;
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nact_w * nocc_w * size_ab > 0 ){
            const int shift_G = shift_G_nonactive( indices, Iw, 0, IL, +1 );
            for ( int ik = 0; ik < nocc_w; ik++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_kw = fock->get( Iw, ik, nocc_w + w );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_kw, FFG_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ab = 0; ab < size_ab; ab++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_L * ab;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE_R * ( shift_G + ik + nocc_w * ab );
                  char notrans = 'N';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ab = 0; ab < size_ab; ab++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_L * ab;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE_R * ( shift_G + ik + nocc_w * ab );
                  char trans = 'T';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }
   
   // FFG triplet: < TF_cxdy E_kw TG_aibt > = 2 delta_ac delta_bd delta_ik FFG_triplet[ Ic x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iac x Ibd == Ix x Iy
      int SIZE_L = size_F_triplet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iac x Ibd x Iik
         int SIZE_R = size_G[ IR ];
         const int Iw = Irreps::directProd( IL, IR ); // Iw == Iik
         const int nocc_w = indices->getNOCC( Iw );
         const int nact_w = indices->getNDMRG( Iw );
         int linsize = 0;
         for ( int Iac = 0; Iac < num_irreps; Iac++ ){
            const int Ibd = Irreps::directProd( Iac, IL );
            if ( Iac <= Ibd ){
               const int nvir_ac = indices->getNVIRT( Iac );
               linsize += (( IL == 0 ) ? ( nvir_ac * ( nvir_ac - 1 ) ) / 2 : nvir_ac * indices->getNVIRT( Ibd ) );
            }
         }
         const int size_ab = linsize;
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nact_w * nocc_w * size_ab > 0 ){
            const int shift_G = shift_G_nonactive( indices, Iw, 0, IL, -1 );
            for ( int ik = 0; ik < nocc_w; ik++ ){
               for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
               for ( int w = 0; w < nact_w; w++ ){
                  double f_kw = fock->get( Iw, ik, nocc_w + w );
                  int inc1 = 1;
                  daxpy_( &total_size, &f_kw, FFG_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
               }
               for ( int ab = 0; ab < size_ab; ab++ ){
                  double * target = result + jump[ IL + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_L * ab;
                  double * origin = vector + jump[ IR + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE_R * ( shift_G + ik + nocc_w * ab );
                  char notrans = 'N';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &notrans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
               for ( int ab = 0; ab < size_ab; ab++ ){
                  double * origin = vector + jump[ IL + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_L * ab;
                  double * target = result + jump[ IR + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE_R * ( shift_G + ik + nocc_w * ab );
                  char trans = 'T';
                  double two = 2.0;
                  double one = 1.0;
                  int inc1 = 1;
                  dgemv_( &trans, &SIZE_L, &SIZE_R, &two, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
               }
            }
         }
      }
   }

   // FEH singlet: < SE_xkdl E_wc SH_aibj > = 2 delta_ik delta_jl ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FEH[ Ix ][ w ][ x ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ixw == Ic == Iik x Ijl x Id
      int SIZE = size_E[ IL ];
      const int nocc_w = indices->getNOCC( IL );
      const int nact_w = indices->getNDMRG( IL );
      const int n_oa_w = nocc_w + nact_w;
      const int nvir_w = indices->getNVIRT( IL );
      for ( int c = 0; c < nvir_w; c++ ){

         // workspace_c[ x ] = sum_w f_wc FEH[ Ixw == IL == Ic ][ w ][ x ]
         for ( int cnt = 0; cnt < SIZE; cnt++ ){ workspace[ cnt ] = 0.0; }
         for ( int w = 0; w < nact_w; w++ ){
            double f_wc = fock->get( IL, nocc_w + w, n_oa_w + c );
            int inc1 = 1;
            daxpy_( &SIZE, &f_wc, FEH[ IL ][ w ], &inc1, workspace, &inc1 );
         }

         for ( int Id = 0; Id < num_irreps; Id++ ){

            const int Icenter = Irreps::directProd( IL, Id );
            const int nvir_d = indices->getNVIRT( Id );

            if ( IL < Id ){ // Ic < Id
               for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                  const int Ij = Irreps::directProd( Ii, Icenter );
                  if ( Ii < Ij ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Ii, Ij, IL, Id, +1 );
                     const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           const double prefactor = 2 * vector[ jump_H + ij + size_ij * ( c + nvir_w * d ) ];
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                           }
                           result[ jump_H + ij + size_ij * ( c + nvir_w * d ) ] += 2 * value;
                        }
                     }
                  }
               }
            }

            if ( IL == Id ){ // Ic == Id == Ia == Ib --> Iik == Ijl
               for ( int Iij = 0; Iij < num_irreps; Iij++ ){
                  const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * shift_E_nonactive( indices, Id,  Iij, Iij, +1 );
                  const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Iij, Iij, IL,  Id, +1 );
                  const int size_ij = ( indices->getNOCC( Iij ) * ( indices->getNOCC( Iij ) + 1 ) ) / 2;
                  for ( int d = 0; d < c; d++ ){
                     const int count_cd = d + ( c * ( c + 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        const double prefactor = 2 * vector[ jump_H + ij + size_ij * count_cd ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int d = c; d < nvir_d; d++ ){
                     const int count_cd = c + ( d * ( d + 1 ) ) / 2;
                     const double factor = 2 * (( c == d ) ? SQRT2 : 1.0 );
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        const double prefactor = factor * vector[ jump_H + ij + size_ij * count_cd ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int d = 0; d < c; d++ ){
                     const int count_cd = d + ( c * ( c + 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                        }
                        result[ jump_H + ij + size_ij * count_cd ] += 2 * value;
                     }
                  }
                  for ( int d = c; d < nvir_d; d++ ){
                     const int count_cd = c + ( d * ( d + 1 ) ) / 2;
                     const double factor = 2 * (( c == d ) ? SQRT2 : 1.0 );
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                        }
                        result[ jump_H + ij + size_ij * count_cd ] += factor * value;
                     }
                  }
               }
            }

            if ( IL > Id ){ // Ic > Id
               for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                  const int Ij = Irreps::directProd( Ii, Icenter );
                  if ( Ii < Ij ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Ii, Ij, Id, IL, +1 );
                     const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           const double prefactor = 2 * vector[ jump_H + ij + size_ij * ( d + nvir_d * c ) ]; // ( d + nvir_d * c ) because Id < Ic
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                           }
                           result[ jump_H + ij + size_ij * ( d + nvir_d * c ) ] += 2 * value; // ( d + nvir_d * c ) because Id < Ic
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FEH triplet: < TE_xkdl E_wc TH_aibj > = 6 delta_ik delta_jl ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FEH[ Ix ][ w ][ x ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ixw == Ic == Iik x Ijl x Id
      int SIZE = size_E[ IL ];
      const int nocc_w = indices->getNOCC( IL );
      const int nact_w = indices->getNDMRG( IL );
      const int n_oa_w = nocc_w + nact_w;
      const int nvir_w = indices->getNVIRT( IL );
      for ( int c = 0; c < nvir_w; c++ ){

         // workspace_c[ x ] = sum_w f_wc FEH[ Ixw == IL == Ic ][ w ][ x ]
         for ( int cnt = 0; cnt < SIZE; cnt++ ){ workspace[ cnt ] = 0.0; }
         for ( int w = 0; w < nact_w; w++ ){
            double f_wc = fock->get( IL, nocc_w + w, n_oa_w + c );
            int inc1 = 1;
            daxpy_( &SIZE, &f_wc, FEH[ IL ][ w ], &inc1, workspace, &inc1 );
         }

         for ( int Id = 0; Id < num_irreps; Id++ ){

            const int Icenter = Irreps::directProd( IL, Id );
            const int nvir_d = indices->getNVIRT( Id );

            if ( IL < Id ){ // Ic < Id
               for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                  const int Ij = Irreps::directProd( Ii, Icenter );
                  if ( Ii < Ij ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Ii, Ij, IL, Id, -1 );
                     const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           const double prefactor = 6 * vector[ jump_H + ij + size_ij * ( c + nvir_w * d ) ];
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                           }
                           result[ jump_H + ij + size_ij * ( c + nvir_w * d ) ] += 6 * value;
                        }
                     }
                  }
               }
            }

            if ( IL == Id ){ // Ic == Id == Ia == Ib --> Iik == Ijl
               for ( int Iij = 0; Iij < num_irreps; Iij++ ){
                  const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * shift_E_nonactive( indices, Id,  Iij, Iij, -1 );
                  const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Iij, Iij, IL,  Id, -1 );
                  const int size_ij = ( indices->getNOCC( Iij ) * ( indices->getNOCC( Iij ) - 1 ) ) / 2;
                  for ( int d = 0; d < c; d++ ){
                     const int count_cd = d + ( c * ( c - 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        const double prefactor = 6 * vector[ jump_H + ij + size_ij * count_cd ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] -= prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int count_cd = c + ( d * ( d - 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        const double prefactor = 6 * vector[ jump_H + ij + size_ij * count_cd ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int d = 0; d < c; d++ ){
                     const int count_cd = d + ( c * ( c - 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                        }
                        result[ jump_H + ij + size_ij * count_cd ] -= 6 * value;
                     }
                  }
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int count_cd = c + ( d * ( d - 1 ) ) / 2;
                     for ( int ij = 0; ij < size_ij; ij++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                        }
                        result[ jump_H + ij + size_ij * count_cd ] += 6 * value;
                     }
                  }
               }
            }

            if ( IL > Id ){ // Ic > Id
               for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                  const int Ij = Irreps::directProd( Ii, Icenter );
                  if ( Ii < Ij ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Ii, Ij, Id, IL, -1 );
                     const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           const double prefactor = 6 * vector[ jump_H + ij + size_ij * ( d + nvir_d * c ) ]; // ( d + nvir_d * c ) because Id < Ic
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] -= prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int d = 0; d < nvir_d; d++ ){
                        for ( int ij = 0; ij < size_ij; ij++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_E + x + SIZE * ( d + nvir_d * ij ) ] * workspace[ x ];
                           }
                           result[ jump_H + ij + size_ij * ( d + nvir_d * c ) ] -= 6 * value; // ( d + nvir_d * c ) because Id < Ic
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FGH singlet: < SG_cldx E_kw SH_aibj > = 2 delta_ac delta_bd ( delta_il delta_jk + delta_ik delta_jl ) / sqrt( 1 + delta_ij ) FGH[ Ix ][ w ][ x ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ixw == Ik == Iac x Ibd x Il
      int SIZE = size_G[ IL ];
      const int nocc_w = indices->getNOCC( IL );
      const int nact_w = indices->getNDMRG( IL );
      for ( int k = 0; k < nocc_w; k++ ){

         // workspace_k[ x ] = sum_w f_kw FGH[ Ixw == IL == Ik ][ w ][ x ]
         for ( int cnt = 0; cnt < SIZE; cnt++ ){ workspace[ cnt ] = 0.0; }
         for ( int w = 0; w < nact_w; w++ ){
            double f_kw = fock->get( IL, k, nocc_w + w );
            int inc1 = 1;
            daxpy_( &SIZE, &f_kw, FGH[ IL ][ w ], &inc1, workspace, &inc1 );
         }

         for ( int Il = 0; Il < num_irreps; Il++ ){

            const int Icenter = Irreps::directProd( IL, Il );
            const int nocc_l = indices->getNOCC( Il );

            if ( IL < Il ){ // Ik < Il
               for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                  const int Ib = Irreps::directProd( Ia, Icenter );
                  if ( Ia < Ib ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, IL, Il, Ia, Ib, +1 );
                     const int size_ij = nocc_w * nocc_l;
                     const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const double prefactor = 2 * vector[ jump_H + k + nocc_w * l + size_ij * ab ];
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                           }
                           result[ jump_H + k + nocc_w * l + size_ij * ab ] += 2 * value;
                        }
                     }
                  }
               }
            }

            if ( IL == Il ){ // Ik == Il == Ii == Ij --> Iac == Ibd
               for ( int Iab = 0; Iab < num_irreps; Iab++ ){
                  const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Iab, Iab, +1 );
                  const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, IL, Il, Iab, Iab, +1 );
                  const int size_ij = ( nocc_w * ( nocc_w + 1 ) ) / 2;
                  const int size_ab = ( indices->getNVIRT( Iab ) * ( indices->getNVIRT( Iab ) + 1 ) ) / 2;
                  for ( int l = 0; l < k; l++ ){
                     const int count_kl = l + ( k * ( k + 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        const double prefactor = 2 * vector[ jump_H + count_kl + size_ij * ab ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int l = k; l < nocc_l; l++ ){
                     const int count_kl = k + ( l * ( l + 1 ) ) / 2;
                     const double factor = 2 * (( k == l ) ? SQRT2 : 1.0 );
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        const double prefactor = factor * vector[ jump_H + count_kl + size_ij * ab ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int l = 0; l < k; l++ ){
                     const int count_kl = l + ( k * ( k + 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                        }
                        result[ jump_H + count_kl + size_ij * ab ] += 2 * value;
                     }
                  }
                  for ( int l = k; l < nocc_l; l++ ){
                     const int count_kl = k + ( l * ( l + 1 ) ) / 2;
                     const double factor = 2 * (( k == l ) ? SQRT2 : 1.0 );
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                        }
                        result[ jump_H + count_kl + size_ij * ab ] += factor * value;
                     }
                  }
               }
            }

            if ( IL > Il ){ // Ik > Il
               for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                  const int Ib = Irreps::directProd( Ia, Icenter );
                  if ( Ia < Ib ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Il, IL, Ia, Ib, +1 );
                     const int size_ij = nocc_w * nocc_l;
                     const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const double prefactor = 2 * vector[ jump_H + l + nocc_l * k + size_ij * ab ]; // ( l + nocc_l * k ) because Il < Ik
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                           }
                           result[ jump_H + l + nocc_l * k + size_ij * ab ] += 2 * value; // ( l + nocc_l * k ) because Il < Ik
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FGH triplet: < TG_cldx E_kw TH_aibj > = 6 delta_ac delta_bd ( delta_il delta_jk - delta_ik delta_jl ) / sqrt( 1 + delta_ij ) FGH[ Ix ][ w ][ x ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ixw == Ik == Iac x Ibd x Il
      int SIZE = size_G[ IL ];
      const int nocc_w = indices->getNOCC( IL );
      const int nact_w = indices->getNDMRG( IL );
      for ( int k = 0; k < nocc_w; k++ ){

         // workspace_k[ x ] = sum_w f_kw FGH[ Ixw == IL == Ik ][ w ][ x ]
         for ( int cnt = 0; cnt < SIZE; cnt++ ){ workspace[ cnt ] = 0.0; }
         for ( int w = 0; w < nact_w; w++ ){
            double f_kw = fock->get( IL, k, nocc_w + w );
            int inc1 = 1;
            daxpy_( &SIZE, &f_kw, FGH[ IL ][ w ], &inc1, workspace, &inc1 );
         }

         for ( int Il = 0; Il < num_irreps; Il++ ){

            const int Icenter = Irreps::directProd( IL, Il );
            const int nocc_l = indices->getNOCC( Il );

            if ( IL < Il ){ // Ik < Il
               for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                  const int Ib = Irreps::directProd( Ia, Icenter );
                  if ( Ia < Ib ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, IL, Il, Ia, Ib, -1 );
                     const int size_ij = nocc_w * nocc_l;
                     const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const double prefactor = 6 * vector[ jump_H + k + nocc_w * l + size_ij * ab ];
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] -= prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                           }
                           result[ jump_H + k + nocc_w * l + size_ij * ab ] -= 6 * value;
                        }
                     }
                  }
               }
            }

            if ( IL == Il ){ // Ik == Il == Ii == Ij --> Iac == Ibd
               for ( int Iab = 0; Iab < num_irreps; Iab++ ){
                  const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Iab, Iab, -1 );
                  const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, IL, Il, Iab, Iab, -1 );
                  const int size_ij = ( nocc_w * ( nocc_w - 1 ) ) / 2;
                  const int size_ab = ( indices->getNVIRT( Iab ) * ( indices->getNVIRT( Iab ) - 1 ) ) / 2;
                  for ( int l = 0; l < k; l++ ){
                     const int count_kl = l + ( k * ( k - 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        const double prefactor = 6 * vector[ jump_H + count_kl + size_ij * ab ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int count_kl = k + ( l * ( l - 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        const double prefactor = 6 * vector[ jump_H + count_kl + size_ij * ab ];
                        for ( int x = 0; x < nact_w; x++ ){
                           result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] -= prefactor * workspace[ x ];
                        }
                     }
                  }
                  for ( int l = 0; l < k; l++ ){
                     const int count_kl = l + ( k * ( k - 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                        }
                        result[ jump_H + count_kl + size_ij * ab ] += 6 * value;
                     }
                  }
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int count_kl = k + ( l * ( l - 1 ) ) / 2;
                     for ( int ab = 0; ab < size_ab; ab++ ){
                        double value = 0.0;
                        for ( int x = 0; x < nact_w; x++ ){
                           value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                        }
                        result[ jump_H + count_kl + size_ij * ab ] -= 6 * value;
                     }
                  }
               }
            }

            if ( IL > Il ){ // Ik > Il
               for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                  const int Ib = Irreps::directProd( Ia, Icenter );
                  if ( Ia < Ib ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Il, IL, Ia, Ib, -1 );
                     const int size_ij = nocc_w * nocc_l;
                     const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const double prefactor = 6 * vector[ jump_H + l + nocc_l * k + size_ij * ab ]; // ( l + nocc_l * k ) because Il < Ik
                           for ( int x = 0; x < nact_w; x++ ){
                              result[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] += prefactor * workspace[ x ];
                           }
                        }
                     }
                     for ( int l = 0; l < nocc_l; l++ ){
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           double value = 0.0;
                           for ( int x = 0; x < nact_w; x++ ){
                              value += vector[ jump_G + x + SIZE * ( l + nocc_l * ab ) ] * workspace[ x ];
                           }
                           result[ jump_H + l + nocc_l * k + size_ij * ab ] += 6 * value; // ( l + nocc_l * k ) because Il < Ik
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FDE singlet: < D(blxy) E_kw SE_tiaj > = 1 delta_ab ( delta_ik delta_jl + delta_il delta_jk ) / sqrt( 1 + delta_ij ) FDE_singlet[ Ib x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Ib x Il
      int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ii x Ij
         int SIZE_R = size_E[ IR ];
         const int Ikw = Irreps::directProd( IL, IR ); // Ikw == Ik == Iw
         const int nocc_kw = indices->getNOCC( Ikw );
         const int nact_kw = indices->getNDMRG( Ikw );
         int total_size = SIZE_L * SIZE_R;
         for ( int k = 0; k < nocc_kw; k++ ){
            for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
            for ( int w = 0; w < nact_kw; w++ ){
               double f_kw = fock->get( Ikw, k, nocc_kw + w );
               int inc1 = 1;
               daxpy_( &total_size, &f_kw, FDE_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
            }
            for ( int Iab = 0; Iab < num_irreps; Iab++ ){
               const int nvir_ab = indices->getNVIRT( Iab );
               const int Il = Irreps::directProd( Iab, IL );
               const int nocc_l = indices->getNOCC( Il );
               const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D         ] + SIZE_L * shift_D_nonactive( indices, Il, Iab );
               const int jump_E = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE_R * (( Ikw <= Il ) ? shift_E_nonactive( indices, Iab, Ikw, Il,  +1 )
                                                                                                               : shift_E_nonactive( indices, Iab, Il,  Ikw, +1 ));
               if ( Ikw == Il ){ // irrep_k == irrep_l
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int cnt_kl = (( k < l ) ? ( k + ( l * ( l + 1 ) ) / 2 ) : ( l + ( k * ( k + 1 ) ) / 2 ));
                     double factor = (( k == l ) ? SQRT2 : 1.0 );
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * target = result + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * origin = vector + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * origin = vector + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * target = result + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               } else { // irrep_k != irrep_l
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int cnt_kl = (( Ikw < Il ) ? ( k + nocc_kw * l ) : ( l + nocc_l * k ));
                     double factor = 1.0;
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * target = result + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * origin = vector + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * origin = vector + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * target = result + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               }
            }
         }
      }
   }

   // FDE triplet: < D(blxy) E_kw TE_tiaj > = 3 delta_ab ( delta_ik delta_jl - delta_il delta_jk ) / sqrt( 1 + delta_ij ) FDE_triplet[ Ib x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Ib x Il
      int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ii x Ij
         int SIZE_R = size_E[ IR ];
         const int Ikw = Irreps::directProd( IL, IR ); // Ikw == Ik == Iw
         const int nocc_kw = indices->getNOCC( Ikw );
         const int nact_kw = indices->getNDMRG( Ikw );
         int total_size = SIZE_L * SIZE_R;
         for ( int k = 0; k < nocc_kw; k++ ){
            for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
            for ( int w = 0; w < nact_kw; w++ ){
               double f_kw = fock->get( Ikw, k, nocc_kw + w );
               int inc1 = 1;
               daxpy_( &total_size, &f_kw, FDE_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
            }
            for ( int Iab = 0; Iab < num_irreps; Iab++ ){
               const int nvir_ab = indices->getNVIRT( Iab );
               const int Il = Irreps::directProd( Iab, IL );
               const int nocc_l = indices->getNOCC( Il );
               const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D         ] + SIZE_L * shift_D_nonactive( indices, Il, Iab );
               const int jump_E = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE_R * (( Ikw <= Il ) ? shift_E_nonactive( indices, Iab, Ikw, Il,  -1 )
                                                                                                               : shift_E_nonactive( indices, Iab, Il,  Ikw, -1 ));
               if ( Ikw == Il ){ // irrep_k == irrep_l
                  for ( int l = 0; l < k; l++ ){
                     const int cnt_kl = l + ( k * ( k - 1 ) ) / 2;
                     double factor = -3.0;
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * target = result + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * origin = vector + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * origin = vector + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * target = result + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int cnt_kl = k + ( l * ( l - 1 ) ) / 2;
                     double factor = 3.0;
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * target = result + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * origin = vector + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * origin = vector + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * target = result + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               } else { // irrep_k != irrep_l
                  for ( int l = 0; l < nocc_l; l++ ){
                     const int cnt_kl = (( Ikw < Il ) ? ( k + nocc_kw * l ) : ( l + nocc_l * k ));
                     double factor = (( Ikw < Il ) ? 3.0 : -3.0 );
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * target = result + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * origin = vector + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ab = 0; ab < nvir_ab; ab++ ){
                        double * origin = vector + jump_D + SIZE_L * ( l + nocc_l * ab );
                        double * target = result + jump_E + SIZE_R * ( ab + nvir_ab * cnt_kl );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               }
            }
         }
      }
   }

   // FDG singlet: < D(djxy) E_wc SG_aibt > = 1 delta_ij ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FDG_singlet[ Ij x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Id x Ij
      int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ib x Ii
         int SIZE_R = size_E[ IR ];
         const int Iwc = Irreps::directProd( IL, IR ); // Iwc == Ic == Iw
         const int nocc_wc = indices->getNOCC( Iwc );
         const int nact_wc = indices->getNDMRG( Iwc );
         const int n_oa_wc = nocc_wc + nact_wc;
         const int nvir_wc = indices->getNVIRT( Iwc );
         int total_size = SIZE_L * SIZE_R;
         for ( int c = 0; c < nvir_wc; c++ ){
            for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
            for ( int w = 0; w < nact_wc; w++ ){
               double f_wc = fock->get( Iwc, nocc_wc + w, n_oa_wc + c );
               int inc1 = 1;
               daxpy_( &total_size, &f_wc, FDG_singlet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
            }
            for ( int Iij = 0; Iij < num_irreps; Iij++ ){
               const int nocc_ij = indices->getNOCC( Iij );
               const int Id = Irreps::directProd( Iij, IL );
               const int nvir_d = indices->getNVIRT( Id );
               const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_L * shift_D_nonactive( indices, Iij, Id );
               const int jump_G = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE_R * (( Iwc <= Id ) ? shift_G_nonactive( indices, Iij, Iwc, Id,  +1 )
                                                                                                               : shift_G_nonactive( indices, Iij, Id,  Iwc, +1 ));
               if ( Iwc == Id ){ // irrep_c == irrep_d
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int cnt_cd = (( c < d ) ? ( c + ( d * ( d + 1 ) ) / 2 ) : ( d + ( c * ( c + 1 ) ) / 2 ));
                     double factor = (( c == d ) ? SQRT2 : 1.0 );
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * target = result + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * origin = vector + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * origin = vector + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * target = result + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               } else { // irrep_c != irrep_d
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int cnt_cd = (( Iwc < Id ) ? ( c + nvir_wc * d ) : ( d + nvir_d * c ));
                     double factor = 1.0;
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * target = result + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * origin = vector + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * origin = vector + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * target = result + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               }
            }
         }
      }
   }
   
   // FDG triplet: < D(djxy) E_wc TG_aibt > = 3 delta_ij ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FDG_triplet[ Ij x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Id x Ij
      int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ib x Ii
         int SIZE_R = size_E[ IR ];
         const int Iwc = Irreps::directProd( IL, IR ); // Iwc == Ic == Iw
         const int nocc_wc = indices->getNOCC( Iwc );
         const int nact_wc = indices->getNDMRG( Iwc );
         const int n_oa_wc = nocc_wc + nact_wc;
         const int nvir_wc = indices->getNVIRT( Iwc );
         int total_size = SIZE_L * SIZE_R;
         for ( int c = 0; c < nvir_wc; c++ ){
            for ( int cnt = 0; cnt < total_size; cnt++ ){ workspace[ cnt ] = 0.0; }
            for ( int w = 0; w < nact_wc; w++ ){
               double f_wc = fock->get( Iwc, nocc_wc + w, n_oa_wc + c );
               int inc1 = 1;
               daxpy_( &total_size, &f_wc, FDG_triplet[ IL ][ IR ][ w ], &inc1, workspace, &inc1 );
            }
            for ( int Iij = 0; Iij < num_irreps; Iij++ ){
               const int nocc_ij = indices->getNOCC( Iij );
               const int Id = Irreps::directProd( Iij, IL );
               const int nvir_d = indices->getNVIRT( Id );
               const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_L * shift_D_nonactive( indices, Iij, Id );
               const int jump_G = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE_R * (( Iwc <= Id ) ? shift_G_nonactive( indices, Iij, Iwc, Id,  -1 )
                                                                                                               : shift_G_nonactive( indices, Iij, Id,  Iwc, -1 ));
               if ( Iwc == Id ){ // irrep_c == irrep_d
                  for ( int d = 0; d < c; d++ ){
                     const int cnt_cd = d + ( c * ( c - 1 ) ) / 2;
                     double factor = -3.0;
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * target = result + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * origin = vector + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * origin = vector + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * target = result + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int cnt_cd = c + ( d * ( d - 1 ) ) / 2;
                     double factor = 3.0;
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * target = result + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * origin = vector + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * origin = vector + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * target = result + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               } else { // irrep_c != irrep_d
                  for ( int d = 0; d < nvir_d; d++ ){
                     const int cnt_cd = (( Iwc < Id ) ? ( c + nvir_wc * d ) : ( d + nvir_d * c ));
                     double factor = (( Iwc < Id ) ? 3.0 : -3.0 );
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * target = result + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * origin = vector + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char notrans = 'N';
                        int inc1 = 1;
                        dgemv_( &notrans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                     for ( int ij = 0; ij < nocc_ij; ij++ ){
                        double * origin = vector + jump_D + SIZE_L * ( ij + nocc_ij * d );
                        double * target = result + jump_G + SIZE_R * ( ij + nocc_ij * cnt_cd );
                        double one = 1.0;
                        char trans = 'T';
                        int inc1 = 1;
                        dgemv_( &trans, &SIZE_L, &SIZE_R, &factor, workspace, &SIZE_L, origin, &inc1, &one, target, &inc1 );
                     }
                  }
               }
            }
         }
      }
   }

   delete [] workspace;

}

int CheMPS2::CASPT2::shift_H_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_j, const int irrep_a, const int irrep_b, const int ST ){

   assert( irrep_i <= irrep_j );
   assert( irrep_a <= irrep_b );
   const int irrep_prod = Irreps::directProd( irrep_i, irrep_j );
   const int irrep_virt = Irreps::directProd( irrep_a, irrep_b );
   assert( irrep_prod == irrep_virt );
   const int n_irreps = idx->getNirreps();

   int shift = 0;
   if ( irrep_prod == 0 ){
      for ( int Iij = 0; Iij < n_irreps; Iij++ ){
         for ( int Iab = 0; Iab < n_irreps; Iab++ ){
            if (( irrep_i == Iij ) && ( irrep_j == Iij ) && ( irrep_a == Iab ) && ( irrep_b == Iab )){
               Iij = n_irreps;
               Iab = n_irreps;
            } else {
               shift += ( idx->getNOCC( Iij ) * ( idx->getNOCC( Iij ) + ST ) * idx->getNVIRT( Iab ) * ( idx->getNVIRT( Iab ) + ST ) ) / 4;
            }
         }
      }
   } else {
      for ( int Ii = 0; Ii < n_irreps; Ii++ ){
         const int Ij = Irreps::directProd( irrep_prod, Ii );
         if ( Ii < Ij ){
            for ( int Ia = 0; Ia < n_irreps; Ia++ ){
               const int Ib = Irreps::directProd( irrep_prod, Ia );
               if ( Ia < Ib ){
                  if (( irrep_i == Ii ) && ( irrep_j == Ij ) && ( irrep_a == Ia ) && ( irrep_b == Ib )){
                     Ii = n_irreps;
                     Ia = n_irreps;
                  } else {
                     shift += idx->getNOCC( Ii ) * idx->getNOCC( Ij ) * idx->getNVIRT( Ia ) * idx->getNVIRT( Ib );
                  }
               }
            }
         }
      }
   }
   return shift;

}

int CheMPS2::CASPT2::shift_G_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_a, const int irrep_b, const int ST ){

   assert( irrep_a <= irrep_b );
   const int irrep_virt = Irreps::directProd( irrep_a,    irrep_b );
   const int irrep_prod = Irreps::directProd( irrep_virt, irrep_i );
   const int n_irreps = idx->getNirreps();

   int shift = 0;
   for ( int Ii = 0; Ii < n_irreps; Ii++ ){
      const int Ivirt = Irreps::directProd( Ii, irrep_prod );
      if ( Ivirt == 0 ){
         for ( int Iab = 0; Iab < n_irreps; Iab++ ){
            if (( irrep_i == Ii ) && ( irrep_a == Iab ) && ( irrep_b == Iab )){
               Ii  = n_irreps;
               Iab = n_irreps;
            } else {
               shift += ( idx->getNOCC( Ii ) * idx->getNVIRT( Iab ) * ( idx->getNVIRT( Iab ) + ST ) ) / 2;
            }
         }
      } else {
         for ( int Ia = 0; Ia < n_irreps; Ia++ ){
            const int Ib = Irreps::directProd( Ivirt, Ia );
            if ( Ia < Ib ){
               if (( irrep_i == Ii ) && ( irrep_a == Ia ) && ( irrep_b == Ib )){
                  Ii = n_irreps;
                  Ia = n_irreps;
               } else {
                  shift += idx->getNOCC( Ii ) * idx->getNVIRT( Ia ) * idx->getNVIRT( Ib );
               }
            }
         }
      }
   }
   return shift;

}

int CheMPS2::CASPT2::shift_E_nonactive( const DMRGSCFindices * idx, const int irrep_a, const int irrep_i, const int irrep_j, const int ST ){

   assert( irrep_i <= irrep_j );
   const int irrep_occ  = Irreps::directProd( irrep_i,   irrep_j );
   const int irrep_prod = Irreps::directProd( irrep_occ, irrep_a );
   const int n_irreps = idx->getNirreps();

   int shift = 0;
   for ( int Ia = 0; Ia < n_irreps; Ia++ ){
      const int Iocc = Irreps::directProd( Ia, irrep_prod );
      if ( Iocc == 0 ){
         for ( int Iij = 0; Iij < n_irreps; Iij++ ){
            if (( irrep_a == Ia ) && ( irrep_i == Iij ) && ( irrep_j == Iij )){
               Ia  = n_irreps;
               Iij = n_irreps;
            } else {
               shift += ( idx->getNVIRT( Ia ) * idx->getNOCC( Iij ) * ( idx->getNOCC( Iij ) + ST ) ) / 2;
            }
         }
      } else {
         for ( int Ii = 0; Ii < n_irreps; Ii++ ){
            const int Ij = Irreps::directProd( Iocc, Ii );
            if ( Ii < Ij ){
               if (( irrep_a == Ia ) && ( irrep_i == Ii ) && ( irrep_j == Ij )){
                  Ia = n_irreps;
                  Ii = n_irreps;
               } else {
                  shift += idx->getNVIRT( Ia ) * idx->getNOCC( Ii ) * idx->getNOCC( Ij );
               }
            }
         }
      }
   }
   return shift;

}

int CheMPS2::CASPT2::shift_F_nonactive( const DMRGSCFindices * idx, const int irrep_a, const int irrep_b, const int ST ){

   assert( irrep_a <= irrep_b );
   const int irr_prod = Irreps::directProd( irrep_a, irrep_b );
   const int n_irreps = idx->getNirreps();

   int shift = 0;
   if ( irr_prod == 0 ){
      for ( int Iab = 0; Iab < n_irreps; Iab++ ){
         if (( irrep_a == Iab ) && ( irrep_b == Iab )){
            Iab = n_irreps;
         } else {
            shift += ( idx->getNVIRT( Iab ) * ( idx->getNVIRT( Iab ) + ST ) ) / 2;
         }
      }
   } else {
      for ( int Ia = 0; Ia < n_irreps; Ia++ ){
         const int Ib = Irreps::directProd( Ia, irr_prod );
         if ( Ia < Ib ){
            if (( irrep_a == Ia ) && ( irrep_b == Ib )){
               Ia = n_irreps;
            } else {
               shift += idx->getNVIRT( Ia ) * idx->getNVIRT( Ib );
            }
         }
      }
   }
   return shift;

}

int CheMPS2::CASPT2::shift_B_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_j, const int ST ){

   assert( irrep_i <= irrep_j );
   const int irr_prod = Irreps::directProd( irrep_i, irrep_j );
   const int n_irreps = idx->getNirreps();

   int shift = 0;
   if ( irr_prod == 0 ){
      for ( int Iij = 0; Iij < n_irreps; Iij++ ){
         if (( Iij == irrep_i ) && ( Iij == irrep_j )){
            Iij = n_irreps;
         } else {
            shift += ( idx->getNOCC( Iij ) * ( idx->getNOCC( Iij ) + ST ) ) / 2;
         }
      }
   } else {
      for ( int Ii = 0; Ii < n_irreps; Ii++ ){
         const int Ij = Irreps::directProd( irr_prod, Ii );
         if ( Ii < Ij ){
            if (( Ii == irrep_i ) && ( Ij == irrep_j )){
               Ii = n_irreps;
            } else {
               shift += idx->getNOCC( Ii ) * idx->getNOCC( Ij );
            }
         }
      }
   }
   return shift;

}

int CheMPS2::CASPT2::shift_D_nonactive( const DMRGSCFindices * idx, const int irrep_i, const int irrep_a ){

   const int irrep_ia = Irreps::directProd( irrep_i, irrep_a );
   const int n_irreps = idx->getNirreps();
   
   int shift = 0;
   for ( int Ii = 0; Ii < n_irreps; Ii++ ){
      const int Ia = Irreps::directProd( irrep_ia, Ii );
      if (( Ii == irrep_i ) && ( Ia == irrep_a )){
         Ii = n_irreps;
      } else {
         shift += idx->getNOCC( Ii ) * idx->getNVIRT( Ia );
      }
   }
   return shift;

}

void CheMPS2::CASPT2::diagonal( double * result, const double ovlp_prefactor ) const{

   const double shifted_prefactor = ovlp_prefactor + 2 * sum_f_kk;

   // FAA: < E_zy E_jx ( f_pq E_pq ) E_ti E_uv > = delta_ij * ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_A[ irrep ];
      if ( SIZE > 0 ){
         const int NOCC = indices->getNOCC( irrep );
         for ( int count = 0; count < NOCC; count++ ){
            const double beta = shifted_prefactor - fock->get( irrep, count, count );
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
            for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FAA[ irrep ][ elem ] + beta; }
         }
      }
   }

   // FCC: < E_zy E_xb ( f_pq E_pq ) E_at E_uv > = delta_ab * ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_C[ irrep ];
      if ( SIZE > 0 ){
         const int NVIR = indices->getNVIRT( irrep );
         const int N_OA = indices->getNOCC( irrep ) + indices->getNDMRG( irrep );
         for ( int count = 0; count < NVIR; count++ ){
            const double beta = shifted_prefactor + fock->get( irrep, N_OA + count, N_OA + count );
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
            for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FCC[ irrep ][ elem ] + beta; }
         }
      }
   }

   // FDD: < E_yx E_jb ( f_pq E_pq ) E_ai E_tu > = delta_ab delta_ij ( FDD[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ xytu ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_D[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_a = Irreps::directProd( irrep_i, irrep );
            assert( shift == shift_D_nonactive( indices, irrep_i, irrep_a ) );
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            for ( int i = 0; i < NOCC_i; i++ ){
               const double f_ii = fock->get( irrep_i, i, i );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  const double beta = shifted_prefactor + f_aa - f_ii;
                  const int count = shift + i + NOCC_i * a;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FDD[ irrep ][ elem ] + beta; }
               }
            }
            shift += NOCC_i * NVIR_a;
         }
      }
   }

   // FBB singlet: < SB_xkyl | f_pq E_pq | SB_tiuj > = 2 delta_ik delta_jl ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )
   {
      const int SIZE = size_B_singlet[ 0 ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, +1 ) );
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               const double f_ii = fock->get( irrep_ij, i, i );
               for ( int j = i; j < nocc_ij; j++ ){
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = shifted_prefactor - f_ii - f_jj;
                  const int count = shift + i + ( j * ( j + 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ 0 ][ elem ] + beta ); }
               }
            }
            shift += ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_B_singlet[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               assert( shift == shift_B_nonactive( indices, irrep_i, irrep_j, +1 ) );
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  const double f_ii = fock->get( irrep_i, i, i );
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = shifted_prefactor - f_ii - f_jj;
                     const int count = shift + i + nocc_i * j;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ irrep ][ elem ] + beta ); }
                  }
               }
               shift += nocc_i * nocc_j;
            }
         }
      }
   }

   // FBB triplet: < TB_xkyl | f_pq E_pq | TB_tiuj > = 2 delta_ik delta_jl ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )
   {
      const int SIZE = size_B_triplet[ 0 ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, -1 ) );
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               const double f_ii = fock->get( irrep_ij, i, i );
               for ( int j = i+1; j < nocc_ij; j++ ){
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = shifted_prefactor - f_ii - f_jj;
                  const int count = shift + i + ( j * ( j - 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ 0 ][ elem ] + beta ); }
               }
            }
            shift += ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_B_triplet[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               assert( shift == shift_B_nonactive( indices, irrep_i, irrep_j, -1 ) );
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  const double f_ii = fock->get( irrep_i, i, i );
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = shifted_prefactor - f_ii - f_jj;
                     const int count = shift + i + nocc_i * j;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ irrep ][ elem ] + beta ); }
                  }
               }
               shift += nocc_i * nocc_j;
            }
         }
      }
   }

   // FFF singlet: < SF_cxdy | f_pq E_pq | SF_atbu > = 2 delta_ac delta_bd ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )
   {
      const int SIZE = size_F_singlet[ 0 ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, +1 ) );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a; b < NVIR_ab; b++ ){
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = shifted_prefactor + f_aa + f_bb;
                  const int count = shift + a + ( b * ( b + 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ 0 ][ elem ] + beta ); }
               }
            }
            shift += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_F_singlet[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int irrep_b = Irreps::directProd( irrep, irrep_a );
            if ( irrep_a < irrep_b ){
               assert( shift == shift_F_nonactive( indices, irrep_a, irrep_b, +1 ) );
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int NVIR_b = indices->getNVIRT( irrep_b );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = shifted_prefactor + f_aa + f_bb;
                     const int count = shift + a + NVIR_a * b;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ irrep ][ elem ] + beta ); }
                  }
               }
               shift += NVIR_a * NVIR_b;
            }
         }
      }
   }

   // FFF triplet: < TF_cxdy | f_pq E_pq | TF_atbu > = 2 delta_ac delta_bd ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )
   {
      const int SIZE = size_F_triplet[ 0 ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, -1 ) );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a+1; b < NVIR_ab; b++ ){
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = shifted_prefactor + f_aa + f_bb;
                  const int count = shift + a + ( b * ( b - 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ 0 ][ elem ] + beta ); }
               }
            }
            shift += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_F_triplet[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int irrep_b = Irreps::directProd( irrep, irrep_a );
            if ( irrep_a < irrep_b ){
               assert( shift == shift_F_nonactive( indices, irrep_a, irrep_b, -1 ) );
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int NVIR_b = indices->getNVIRT( irrep_b );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = shifted_prefactor + f_aa + f_bb;
                     const int count = shift + a + NVIR_a * b;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ irrep ][ elem ] + beta ); }
                  }
               }
               shift += NVIR_a * NVIR_b;
            }
         }
      }
   }

   // FEE singlet: < SE_ukbl | f_pq E_pq | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_E[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            const int irrep_occ = Irreps::directProd( irrep_a, irrep );
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
               if ( irrep_i == irrep_j ){
                  assert( shift == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, +1 ) );
                  const int NOCC_ij = indices->getNOCC( irrep_i );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int count = i + ( j * ( j + 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( shift + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += ( NVIR_a * NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  assert( shift == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, +1 ) );
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = i + NOCC_i * j;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( shift + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += NVIR_a * NOCC_i * NOCC_j;
               }
            }
         }
      }
   }

   // FEE triplet: < TE_ukbl | f_pq E_pq | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_E[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            const int irrep_occ = Irreps::directProd( irrep_a, irrep );
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
               if ( irrep_i == irrep_j ){
                  assert( shift == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, -1 ) );
                  const int NOCC_ij = indices->getNOCC( irrep_i );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = i+1; j < NOCC_ij; j++ ){
                        const int count = i + ( j * ( j - 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( shift + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += ( NVIR_a * NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  assert( shift == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, -1 ) );
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = i + NOCC_i * j;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( shift + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += NVIR_a * NOCC_i * NOCC_j;
               }
            }
         }
      }
   }
   
   // FGG singlet: < SG_cjdu | f_pq E_pq | SG_aibt > = 2 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_G[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int irrep_vir = Irreps::directProd( irrep_i, irrep );
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
               if ( irrep_a == irrep_b ){
                  assert( shift == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, +1 ) );
                  const int NVIR_ab = indices->getNVIRT( irrep_a );
                  const int N_OA_ab = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a );
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const int count = a + ( b * ( b + 1 ) ) / 2;
                        const double f_bb = fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( shift + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += ( NOCC_i * NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               }
               if ( irrep_a < irrep_b ){
                  assert( shift == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, +1 ) );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  for ( int a = 0; a < NVIR_a; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = a + NVIR_a * b;
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( shift + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += NOCC_i * NVIR_a * NVIR_b;
               }
            }
         }
      }
   }
   
   // FGG triplet: < TG_cjdu | f_pq E_pq | TG_aibt > = 6 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_G[ irrep ];
      if ( SIZE > 0 ){
         int shift = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int irrep_vir = Irreps::directProd( irrep_i, irrep );
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
               if ( irrep_a == irrep_b ){
                  assert( shift == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, -1 ) );
                  const int NVIR_ab = indices->getNVIRT( irrep_a );
                  const int N_OA_ab = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a );
                     for ( int b = a+1; b < NVIR_ab; b++ ){
                        const int count = a + ( b * ( b - 1 ) ) / 2;
                        const double f_bb = fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( shift + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += ( NOCC_i * NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
               }
               if ( irrep_a < irrep_b ){
                  assert( shift == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, -1 ) );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  for ( int a = 0; a < NVIR_a; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = a + NVIR_a * b;
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( shift + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  shift += NOCC_i * NVIR_a * NVIR_b;
               }
            }
         }
      }
   }

   // FHH singlet and triplet
   {
      int shift = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, +1 ) );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift;
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a; b < NVIR_ab; b++ ){
                  const int cnt_ab = a + ( b * ( b + 1 ) ) / 2;
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_ij, i, i );
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int cnt_ij = i + ( j * ( j + 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_ij, j, j );
                        const double term = shifted_prefactor + f_dot_1dm + f_aa + f_bb - f_ii - f_jj;
                        target[ cnt_ij + size_ij * cnt_ab ] = 4 * term;
                     }
                  }
               }
            }
            shift += size_ij * size_ab;
         }
      }
   }
   {
      int shift = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, -1 ) );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift;
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a+1; b < NVIR_ab; b++ ){
                  const int cnt_ab = a + ( b * ( b - 1 ) ) / 2;
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_ij, i, i );
                     for ( int j = i+1; j < NOCC_ij; j++ ){
                        const int cnt_ij = i + ( j * ( j - 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_ij, j, j );
                        const double term = shifted_prefactor + f_dot_1dm + f_aa + f_bb - f_ii - f_jj;
                        target[ cnt_ij + size_ij * cnt_ab ] = 12 * term;
                     }
                  }
               }
            }
            shift += size_ij * size_ab;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int shift = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_j = Irreps::directProd( irrep, irrep_i );
         if ( irrep_i < irrep_j ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NOCC_j = indices->getNOCC( irrep_j );
            const int size_ij = NOCC_i * NOCC_j;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep, irrep_a );
               if ( irrep_a < irrep_b ){
                  assert( shift == shift_H_nonactive( indices, irrep_i, irrep_j, irrep_a, irrep_b, 0 ) );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int size_ab = NVIR_a * NVIR_b;
                  double * target_singlet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift;
                  double * target_triplet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift;
                  for ( int a = 0; a < NVIR_a; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int cnt_ab = a + NVIR_a * b;
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           for ( int j = 0; j < NOCC_j; j++ ){
                              const int cnt_ij = i + NOCC_i * j;
                              const double f_jj = fock->get( irrep_j, j, j );
                              const double term = shifted_prefactor + f_dot_1dm + f_aa + f_bb - f_ii - f_jj;
                              target_singlet[ cnt_ij + size_ij * cnt_ab ] =  4 * term;
                              target_triplet[ cnt_ij + size_ij * cnt_ab ] = 12 * term;
                           }
                        }
                     }
                  }
                  shift += size_ij * size_ab;
               }
            }
         }
      }
   }

}

void CheMPS2::CASPT2::construct_rhs( const DMRGSCFmatrix * oei, const DMRGSCFintegrals * integrals ){

   /*
      VA:  < H E_ti E_uv > = sum_w ( t_iw + sum_k [ 2 (iw|kk) - (ik|kw) ] ) [ 2 delta_tw Gamma_uv - Gamma_tuwv - delta_wu Gamma_tv ]
                           + sum_xzy (ix|zy) SAA[ Ii ][ xyztuv ]

      VB:  < H E_ti E_uj >
           < S_tiuj | H > = sum_xy (ix|jy) SBB_singlet[ It x Iu ][ xytu ] / sqrt( 1 + delta_ij )
                          = sum_{x<=y} [ (ix|jy) + (iy|jx) ] SBB_singlet[ It x Iu ][ xytu ] * (( x==y ) ? 0.5 : 1.0 ) / sqrt( 1 + delta_ij )
           < T_tiuj | H > = sum_xy (ix|jy) SBB_triplet[ It x Iu ][ xytu ]
                          = sum_{x<y}  [ (ix|jy) - (iy|jx) ] SBB_triplet[ It x Iu ][ xytu ]

      VC:  < H E_at E_uv > = sum_w ( t_wa + sum_k [ 2 (wa|kk) - (wk|ka) ] - sum_y (wy|ya) ) < E_wt E_uv >
                           + sum_wxy (xy|wa) < E_xy E_wt E_uv >
           < H E_at E_uv > = sum_w ( t_wa + sum_k [ 2 (wa|kk) - (wk|ka) ] - sum_y (wy|ya) ) [ Gamma_wutv + delta_ut Gamma_wv ]
                           + sum_zxy (zy|xa) SCC[ Ia ][ xyztuv ]

      VD1: < H E_ai E_tu > = ( t_ia + sum_k [ 2 (ia|kk) - (ik|ka) ] ) [ 2 Gamma_tu ]
                           + sum_xy [ (ia|yx) - 0.5 * (ix|ya) ] SD1D1[ It x Iu ][ xytu ]
           < H E_ai E_tu > = ( t_ia + sum_k [ 2 (ia|kk) - (ik|ka) ] ) [ 2 Gamma_tu ]
                           + sum_xy (ia|yx) SD1D1[ It x Iu ][ xytu ]
                           + sum_xy (ix|ya) SD2D1[ It x Iu ][ xytu ]

      VD2: < H E_ti E_au > = ( t_ia + sum_k [ 2 (ia|kk) - (ik|ka) ] ) [ - Gamma_tu ]
                           + sum_xy (ia|xy) [ - Gamma_xtyu ]
                           + sum_xy (iy|xa) [ - Gamma_txyu ]
                           + sum_y [ 2 (it|ya) - (ia|yt) ] Gamma_yu
           < H E_ti E_au > = ( t_ia + sum_k [ 2 (ia|kk) - (ik|ka) ] ) [ - Gamma_tu ]
                           + sum_xy (ia|yx) SD1D2[ It x Iu ][ xytu ]
                           + sum_xy (ix|ya) SD2D2[ It x Iu ][ xytu ]

      VE:  < H E_ti E_aj >
           < S_tiaj | H > = sum_w [ (aj|wi) + (ai|wj) ] * 1 * SEE[ It ][ wt ] / sqrt( 1 + delta_ij )
           < T_tiaj | H > = sum_w [ (aj|wi) - (ai|wj) ] * 3 * SEE[ It ][ wt ]

      VF:  < H E_at E_bu >
           < S_atbu | H > = sum_xy (ax|by) SFF_singlet[ It x Iu ][ xytu ] / sqrt( 1 + delta_ab )
                          = sum_{x<=y} [ (ax|by) + (ay|bx) ] SFF_singlet[ It x Iu ][ xytu ] * (( x==y ) ? 0.5 : 1.0 ) / sqrt( 1 + delta_ab )
           < T_atbu | H > = sum_xy (ax|by) SFF_triplet[ It x Iu ][ xytu ]
                          = sum_{x<y}  [ (ax|by) - (ay|bx) ] SFF_triplet[ It x Iu ][ xytu ]

      VG:  < H E_ai E_bt >
           < S_aibt | H > = sum_w [ (ai|bw) + (bi|aw) ] * 1 * SGG[ It ][ wt ] / sqrt( 1 + delta_ab )
           < T_aibt | H > = sum_w [ (ai|bw) - (bi|aw) ] * 3 * SGG[ It ][ wt ]

      VH:  < H E_ai E_bj >
           < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ] / sqrt( ( 1 + delta_ij ) * ( 1 + delta_ab ) )
           < T_aibj | H > = 6 * [ (ai|bj) - (aj|bi) ]
   */

   const int LAS = indices->getDMRGcumulative( num_irreps );

   // First construct MAT[p,q] = ( t_pq + sum_k [ 2 (pq|kk) - (pk|kq) ] )
   DMRGSCFmatrix * MAT = new DMRGSCFmatrix( indices );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = indices->getNORB( irrep );
      for ( int row = 0; row < NORB; row++ ){
         for ( int col = row; col < NORB; col++ ){
            double value = oei->get( irrep, row, col );
            for ( int irrep_occ = 0; irrep_occ < num_irreps; irrep_occ++ ){
               const int NOCC = indices->getNOCC( irrep_occ );
               for ( int occ = 0; occ < NOCC; occ++ ){
                  value += ( 2 * integrals->FourIndexAPI( irrep, irrep_occ, irrep, irrep_occ, row, occ, col, occ )
                               - integrals->FourIndexAPI( irrep, irrep_occ, irrep_occ, irrep, row, occ, occ, col ) ); // Physics notation at 4-index
               }
            }
            MAT->set( irrep, row, col, value );
            MAT->set( irrep, col, row, value );
         }
      }
   }

   // and MAT2[p,q] = sum_y (py|yq)
   DMRGSCFmatrix * MAT2 = new DMRGSCFmatrix( indices );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NORB = indices->getNORB( irrep );
      for ( int row = 0; row < NORB; row++ ){
         for ( int col = row; col < NORB; col++ ){
            double value = 0.0;
            for ( int irrep_act = 0; irrep_act < num_irreps; irrep_act++ ){
               const int NOCC = indices->getNOCC( irrep_act );
               const int NACT = indices->getNDMRG( irrep_act );
               for ( int act = 0; act < NACT; act++ ){
                  value += integrals->FourIndexAPI( irrep, irrep_act, irrep_act, irrep, row, NOCC + act, NOCC + act, col ); // Physics notation at 4-index
               }
            }
            MAT2->set( irrep, row, col, value );
            MAT2->set( irrep, col, row, value );
         }
      }
   }

   vector_rhs = new double[ jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ] ];

   const int max_size = get_maxsize();
   double * workspace = new double[ max_size ];

   // VA
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      const int NACT = indices->getNDMRG( irrep );
      const int d_w  = indices->getDMRGcumulative( irrep );
      for ( int count_i = 0; count_i < NOCC; count_i++ ){

         double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + size_A[ irrep ] * count_i;

         // Fill workspace[ xyz ] with (ix|zy)
         // Fill target[ tuv ] with sum_w MAT[i,w] [ 2 delta_tw Gamma_uv - Gamma_tuwv - delta_wu Gamma_tv ]
         //                           = 2 MAT[i,t] Gamma_uv - MAT[i,u] Gamma_tv - sum_w MAT[i,w] Gamma_tuwv
         int jump_xyz = 0;
         for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
            const int occ_x = indices->getNOCC( irrep_x );
            const int num_x = indices->getNDMRG( irrep_x );
            const int d_x   = indices->getDMRGcumulative( irrep_x );
            for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
               const int irrep_z = Irreps::directProd( Irreps::directProd( irrep, irrep_x ), irrep_y );
               const int occ_y = indices->getNOCC( irrep_y );
               const int occ_z = indices->getNOCC( irrep_z );
               const int num_y = indices->getNDMRG( irrep_y );
               const int num_z = indices->getNDMRG( irrep_z );
               const int d_y   = indices->getDMRGcumulative( irrep_y );
               const int d_z   = indices->getDMRGcumulative( irrep_z );

               // workspace[ xyz ] = (ix|zy)
               for ( int z = 0; z < num_z; z++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double ix_zy = integrals->get_coulomb( irrep, irrep_x, irrep_z, irrep_y, count_i, occ_x + x, occ_z + z, occ_y + y );
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ] = ix_zy;
                     }
                  }
               }

               // target[ tuv ] = - sum_w MAT[i,w] Gamma_tuwv
               for ( int z = 0; z < num_z; z++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        double value = 0.0;
                        for ( int w = 0; w < NACT; w++ ){
                           value += MAT->get( irrep, count_i, NOCC + w ) * two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_w + w + LAS * ( d_z + z ) ) ) ];
                        }
                        target[ jump_xyz + x + num_x * ( y + num_y * z ) ] = - value;
                     }
                  }
               }

               // target[ tuv ] += 2 MAT[i,t] Gamma_uv
               if ( irrep_x == irrep ){
                  for ( int z = 0; z < num_z; z++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           target[ jump_xyz + x + num_x * ( y + num_y * z ) ] += 2 * MAT->get( irrep, count_i, occ_x + x ) * one_rdm[ d_y + y + LAS * ( d_z + z ) ];
                        }
                     }
                  }
               }

               // target[ tuv ] -= MAT[i,u] Gamma_tv
               if ( irrep_y == irrep ){
                  for ( int z = 0; z < num_z; z++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           target[ jump_xyz + x + num_x * ( y + num_y * z ) ] -= MAT->get( irrep, count_i, occ_y + y ) * one_rdm[ d_x + x + LAS * ( d_z + z ) ];
                        }
                     }
                  }
               }
               jump_xyz += num_x * num_y * num_z;
            }
         }
         assert( jump_xyz == size_A[ irrep ] );

         // Perform target[ tuv ] += sum_xzy (ix|zy) SAA[ It x Iu x Iv ][ xyztuv ]
         char notrans = 'N';
         int inc1 = 1;
         double one = 1.0;
         dgemv_( &notrans, &jump_xyz, &jump_xyz, &one, SAA[ irrep ], &jump_xyz, workspace, &inc1, &one, target, &inc1 );
      }
      assert( NOCC * size_A[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_A ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] );
   }

   // VC
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      const int NVIR = indices->getNVIRT( irrep );
      const int NACT = indices->getNDMRG( irrep );
      const int N_OA = NOCC + NACT;
      const int d_w  = indices->getDMRGcumulative( irrep );
      for ( int count_a = 0; count_a < NVIR; count_a++ ){

         double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + size_C[ irrep ] * count_a;

         // Fill workspace[ xyz ] with (zy|xa)
         // Fill target[ tuv ] with sum_w (MAT[w,a] - MAT2[w,a]) [ Gamma_wutv + delta_ut Gamma_wv ]
         int jump_xyz = 0;
         for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
            const int occ_x = indices->getNOCC( irrep_x );
            const int num_x = indices->getNDMRG( irrep_x );
            const int d_x   = indices->getDMRGcumulative( irrep_x );
            for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
               const int irrep_z = Irreps::directProd( Irreps::directProd( irrep, irrep_x ), irrep_y );
               const int occ_y = indices->getNOCC( irrep_y );
               const int occ_z = indices->getNOCC( irrep_z );
               const int num_y = indices->getNDMRG( irrep_y );
               const int num_z = indices->getNDMRG( irrep_z );
               const int d_y   = indices->getDMRGcumulative( irrep_y );
               const int d_z   = indices->getDMRGcumulative( irrep_z );

               // workspace[ xyz ] = (zy|xa)
               for ( int z = 0; z < num_z; z++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double zy_xa = integrals->get_coulomb( irrep_z, irrep_y, irrep_x, irrep, occ_z + z, occ_y + y, occ_x + x, N_OA + count_a );
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ] = zy_xa;
                     }
                  }
               }

               // target[ tuv ] = sum_w ( MAT[w,a] - MAT2[w,a] ) Gamma_wutv
               for ( int z = 0; z < num_z; z++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        double value = 0.0;
                        for ( int w = 0; w < NACT; w++ ){
                           value += ( ( MAT->get( irrep, NOCC + w, N_OA + count_a ) - MAT2->get( irrep, NOCC + w, N_OA + count_a ) )
                                    * two_rdm[ d_w + w + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_z + z ) ) ) ] );
                        }
                        target[ jump_xyz + x + num_x * ( y + num_y * z ) ] = value;
                     }
                  }
               }

               // target[ tuv ] += sum_w ( MAT[w,a] - MAT2[w,a] ) delta_ut Gamma_wv
               if (( irrep_z == irrep ) && ( irrep_x == irrep_y )){
                  for ( int z = 0; z < num_z; z++ ){ // v
                     double value = 0.0;
                     for ( int w = 0; w < NACT; w++ ){
                        value += ( ( MAT->get( irrep, NOCC + w, N_OA + count_a ) - MAT2->get( irrep, NOCC + w, N_OA + count_a ) )
                                 * one_rdm[ d_w + w + LAS * ( d_z + z ) ] );
                     }
                     for ( int xy = 0; xy < num_x; xy++ ){ // tu
                        target[ jump_xyz + xy + num_x * ( xy + num_x * z ) ] += value;
                     }
                  }
               }
               jump_xyz += num_x * num_y * num_z;
            }
         }
         assert( jump_xyz == size_C[ irrep ] );

         // Perform target[ tuv ] += sum_zxy (zy|xa) SCC[ It x Iu x Iv ][ xyztuv ]
         char notrans = 'N';
         int inc1 = 1;
         double one = 1.0;
         dgemv_( &notrans, &jump_xyz, &jump_xyz, &one, SCC[ irrep ], &jump_xyz, workspace, &inc1, &one, target, &inc1 );
      }
      assert( NVIR * size_C[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_C ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] );
   }
   delete MAT2;

   // VD1 and VD2
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int shift = 0;
      const int D2JUMP = size_D[ irrep ] / 2;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_a = Irreps::directProd( irrep_i, irrep );
         assert( shift == shift_D_nonactive( indices, irrep_i, irrep_a ) );
         const int NOCC_i  = indices->getNOCC( irrep_i );
         const int N_OA_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
         const int NVIR_a  = indices->getNVIRT( irrep_a );
         for ( int count_i = 0; count_i < NOCC_i; count_i++ ){
            for ( int count_a = 0; count_a < NVIR_a; count_a++ ){

               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + size_D[ irrep ] * ( shift + count_i + NOCC_i * count_a );
               const double MAT_ia = ( ( irrep_i == irrep_a ) ? MAT->get( irrep_i, count_i, N_OA_a + count_a ) : 0.0 );

               /* Fill workspace[          xy ] with (ia|yx)
                  Fill workspace[ D2JUMP + xy ] with (ix|ya)
                       then ( workspace * SDD )_D1 = sum_xy (ia|yx) SD1D1[ It x Iu ][ xytu ]
                                                   + sum_xy (ix|ya) SD2D1[ It x Iu ][ xytu ]
                        and ( workspace * SDD )_D2 = sum_xy (ia|yx) SD1D2[ It x Iu ][ xytu ]
                                                   + sum_xy (ix|ya) SD2D2[ It x Iu ][ xytu ]
                  Fill target[          tu ] = 2 MAT[i,a] Gamma_tu
                  Fill target[ D2JUMP + tu ] = - MAT[i,a] Gamma_tu */
               int jump_xy = 0;
               for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                  const int irrep_y = Irreps::directProd( irrep, irrep_x );
                  const int occ_x = indices->getNOCC( irrep_x );
                  const int occ_y = indices->getNOCC( irrep_y );
                  const int num_x = indices->getNDMRG( irrep_x );
                  const int num_y = indices->getNDMRG( irrep_y );

                  // workspace[ xy ] = (ia|yx)
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double ia_yx = integrals->get_coulomb( irrep_y, irrep_x, irrep_i, irrep_a, occ_y + y, occ_x + x, count_i, N_OA_a + count_a );
                        workspace[ jump_xy + x + num_x * y ] = ia_yx;
                     }
                  }

                  // workspace[ D2JUMP + xy ] = (ix|ya)
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double ix_ya = integrals->get_coulomb( irrep_i, irrep_x, irrep_y, irrep_a, count_i, occ_x + x, occ_y + y, N_OA_a + count_a );
                        workspace[ D2JUMP + jump_xy + x + num_x * y ] = ix_ya;
                     }
                  }

                  // target[ tu          ] = 2 MAT[i,a] Gamma_tu
                  // target[ D2JUMP + tu ] = - MAT[i,a] Gamma_tu
                  if ( irrep == 0 ){
                     const int d_xy = indices->getDMRGcumulative( irrep_x );
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           const double value = MAT_ia * one_rdm[ d_xy + x + LAS * ( d_xy + y ) ];
                           target[          jump_xy + x + num_x * y ] = 2 * value;
                           target[ D2JUMP + jump_xy + x + num_x * y ] =   - value;
                        }
                     }
                  } else {
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           target[          jump_xy + x + num_x * y ] = 0.0;
                           target[ D2JUMP + jump_xy + x + num_x * y ] = 0.0;
                        }
                     }
                  }
                  jump_xy += num_x * num_y;
               }
               jump_xy = 2 * jump_xy;
               assert( jump_xy == size_D[ irrep ] );

               // Perform target += workspace * SDD[ irrep ]
               char notrans = 'N';
               int inc1 = 1;
               double one = 1.0;
               dgemv_( &notrans, &jump_xy, &jump_xy, &one, SDD[ irrep ], &jump_xy, workspace, &inc1, &one, target, &inc1 );
            }
         }
         shift += NOCC_i * NVIR_a;
      }
      assert( shift * size_D[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_D ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] );
   }
   delete MAT;
   const double SQRT_0p5 = sqrt( 0.5 );

   // VB singlet and triplet
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      int shift = 0; // First do SINGLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, +1 ) );
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with [ (ix|jy) + (iy|jx) ] * (( x==y ) ? 0.5 : 1.0 )
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                        const double ix_jy = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y );
                        const double iy_jx = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + y, j, occ_xy + x );
                        workspace[ jump_xy + x + (y*(y+1))/2 ] = ( ix_jy + iy_jx ) * (( x==y ) ? 0.5 : 1.0 );
                     }
                  }
                  jump_xy += ( num_xy * ( num_xy + 1 ) ) / 2;
               }
               assert( jump_xy == size_B_singlet[ 0 ] );

               // Perform target[ tu ] = sum_{x<=y} [ (ix|jy) + (iy|jx) ] SBB_singlet[ It x Iu ][ xytu ] * (( x==y ) ? 0.5 : 1.0 ) / sqrt( 1 + delta_ij )
               char notrans = 'N';
               int inc1 = 1;
               double alpha = (( i == j ) ? SQRT_0p5 : 1.0 );
               double set = 0.0;
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ 0 ] * ( shift + i + (j*(j+1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         shift += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
      }
      assert( shift * size_B_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );

      shift = 0; // Then do TRIPLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, -1 ) );
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i+1; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with [ (ix|jy) - (iy|jx) ]
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                        const double ix_jy = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y );
                        const double iy_jx = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + y, j, occ_xy + x );
                        workspace[ jump_xy + x + (y*(y-1))/2 ] = ix_jy - iy_jx;
                     }
                  }
                  jump_xy += ( num_xy * ( num_xy - 1 ) ) / 2;
               }
               assert( jump_xy == size_B_triplet[ 0 ] );

               // Perform target[ tu ] = sum_{x<y} [ (ix|jy) - (iy|jx) ] SBB_triplet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int inc1 = 1;
               double alpha = 1.0;
               double set = 0.0;
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ 0 ] * ( shift + i + (j*(j-1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         shift += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
      }
      assert( shift * size_B_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int shift = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_j = Irreps::directProd( irrep, irrep_i );
         if ( irrep_i < irrep_j ){
            assert( shift == shift_B_nonactive( indices, irrep_i, irrep_j, 0 ) );
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NOCC_j = indices->getNOCC( irrep_j );
            for ( int i = 0; i < NOCC_i; i++ ){
               for ( int j = 0; j < NOCC_j; j++ ){

                  // Fill workspace[ xy ] with [ (ix|jy) + (iy|jx) ]
                  int jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int occ_x = indices->getNOCC( irrep_x );
                        const int occ_y = indices->getNOCC( irrep_y );
                        const int num_x = indices->getNDMRG( irrep_x );
                        const int num_y = indices->getNDMRG( irrep_y );

                        for ( int y = 0; y < num_y; y++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double ix_jy = integrals->get_coulomb( irrep_i, irrep_x, irrep_j, irrep_y, i, occ_x + x, j, occ_y + y );
                              const double iy_jx = integrals->get_coulomb( irrep_i, irrep_y, irrep_j, irrep_x, i, occ_y + y, j, occ_x + x );
                              workspace[ jump_xy + x + num_x * y ] = ix_jy + iy_jx;
                           }
                        }
                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_B_singlet[ irrep ] );

                  // Perform target[ tu ] = sum_{x<y} [ (ix|jy) + (iy|jx) ] SBB_singlet[ It x Iu ][ xytu ]
                  char notrans = 'N';
                  int inc1 = 1;
                  double alpha = 1.0;
                  double set = 0.0;
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ irrep ] * ( shift + i + NOCC_i * j );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_singlet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );

                  // Fill workspace[ xy ] with [ (ix|jy) - (iy|jx) ]
                  jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int occ_x = indices->getNOCC( irrep_x );
                        const int occ_y = indices->getNOCC( irrep_y );
                        const int num_x = indices->getNDMRG( irrep_x );
                        const int num_y = indices->getNDMRG( irrep_y );

                        for ( int y = 0; y < num_y; y++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double ix_jy = integrals->get_coulomb( irrep_i, irrep_x, irrep_j, irrep_y, i, occ_x + x, j, occ_y + y );
                              const double iy_jx = integrals->get_coulomb( irrep_i, irrep_y, irrep_j, irrep_x, i, occ_y + y, j, occ_x + x );
                              workspace[ jump_xy + x + num_x * y ] = ix_jy - iy_jx;
                           }
                        }
                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_B_triplet[ irrep ] );

                  // Perform target[ tu ] = sum_{x<y} [ (ix|jy) - (iy|jx) ] SBB_triplet[ It x Iu ][ xytu ]
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ irrep ] * ( shift + i + NOCC_i * j );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
               }
            }
            shift += NOCC_i * NOCC_j;
         }
      }
      assert( shift * size_B_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );
      assert( shift * size_B_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
   }

   // VF singlet and triplet
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      int shift = 0; // First do SINGLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, +1 ) );
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with [ (ax|by) + (ay|bx) ] * (( x==y ) ? 0.5 : 1.0 )
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                        const double ax_by = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b );
                        const double ay_bx = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + y, occ_xy + x, N_OA_ab + a, N_OA_ab + b );
                        workspace[ jump_xy + x + (y*(y+1))/2 ] = ( ax_by + ay_bx ) * (( x==y ) ? 0.5 : 1.0 );
                     }
                  }
                  jump_xy += ( num_xy * ( num_xy + 1 ) ) / 2;
               }
               assert( jump_xy == size_F_singlet[ 0 ] );

               // Perform target[ tu ] = sum_{x<=y} [ (ax|by) + (ay|bx) ] SFF_singlet[ It x Iu ][ xytu ] * (( x==y ) ? 0.5 : 1.0 ) / sqrt( 1 + delta_ab )
               char notrans = 'N';
               int inc1 = 1;
               double alpha = (( a == b ) ? SQRT_0p5 : 1.0 );
               double set = 0.0;
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ 0 ] * ( shift + a + (b*(b+1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         shift += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
      }
      assert( shift * size_F_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );

      shift = 0; // Then do TRIPLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, -1 ) );
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a+1; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with [ (ax|by) - (ay|bx) ]
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                        const double ax_by = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b );
                        const double ay_bx = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + y, occ_xy + x, N_OA_ab + a, N_OA_ab + b );
                        workspace[ jump_xy + x + (y*(y-1))/2 ] = ax_by - ay_bx;
                     }
                  }
                  jump_xy += ( num_xy * ( num_xy - 1 ) ) / 2;
               }
               assert( jump_xy == size_F_triplet[ 0 ] );

               // Perform target[ tu ] = sum_{x<y} [ (ax|by) - (ay|bx) ] SFF_triplet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int inc1 = 1;
               double alpha = 1.0;
               double set = 0.0;
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ 0 ] * ( shift + a + (b*(b-1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         shift += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
      }
      assert( shift * size_F_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int shift = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int irrep_b = Irreps::directProd( irrep, irrep_a );
         if ( irrep_a < irrep_b ){
            assert( shift == shift_F_nonactive( indices, irrep_a, irrep_b, 0 ) );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int NVIR_b = indices->getNVIRT( irrep_b );
            for ( int a = 0; a < NVIR_a; a++ ){
               for ( int b = 0; b < NVIR_b; b++ ){

                  // Fill workspace[ xy ] with [ (ax|by) + (ay|bx) ]
                  int jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int occ_x = indices->getNOCC( irrep_x );
                        const int occ_y = indices->getNOCC( irrep_y );
                        const int num_x = indices->getNDMRG( irrep_x );
                        const int num_y = indices->getNDMRG( irrep_y );

                        for ( int y = 0; y < num_y; y++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double ax_by = integrals->get_exchange( irrep_x, irrep_y, irrep_a, irrep_b, occ_x + x, occ_y + y, N_OA_a + a, N_OA_b + b );
                              const double ay_bx = integrals->get_exchange( irrep_y, irrep_x, irrep_a, irrep_b, occ_y + y, occ_x + x, N_OA_a + a, N_OA_b + b );
                              workspace[ jump_xy + x + num_x * y ] = ax_by + ay_bx;
                           }
                        }
                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_F_singlet[ irrep ] );

                  // Perform target[ tu ] = sum_{x<y} [ (ax|by) + (ay|bx) ] SFF_singlet[ It x Iu ][ xytu ]
                  char notrans = 'N';
                  int inc1 = 1;
                  double alpha = 1.0;
                  double set = 0.0;
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ irrep ] * ( shift + a + NVIR_a * b );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_singlet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );

                  // Fill workspace[ xy ] with [ (ax|by) - (ay|bx) ]
                  jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int occ_x = indices->getNOCC( irrep_x );
                        const int occ_y = indices->getNOCC( irrep_y );
                        const int num_x = indices->getNDMRG( irrep_x );
                        const int num_y = indices->getNDMRG( irrep_y );

                        for ( int y = 0; y < num_y; y++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double ax_by = integrals->get_exchange( irrep_x, irrep_y, irrep_a, irrep_b, occ_x + x, occ_y + y, N_OA_a + a, N_OA_b + b );
                              const double ay_bx = integrals->get_exchange( irrep_y, irrep_x, irrep_a, irrep_b, occ_y + y, occ_x + x, N_OA_a + a, N_OA_b + b );
                              workspace[ jump_xy + x + num_x * y ] = ax_by - ay_bx;
                           }
                        }
                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_F_triplet[ irrep ] );

                  // Perform target[ tu ] = sum_{x<y} [ (ax|by) - (ay|bx) ] SFF_triplet[ It x Iu ][ xytu ]
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ irrep ] * ( shift + a + NVIR_a * b );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
               }
            }
            shift += NVIR_a * NVIR_b;
         }
      }
      assert( shift * size_F_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );
      assert( shift * size_F_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
   }
   delete [] workspace;

   // VE singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int occ_t = indices->getNOCC( irrep );
      const int num_t = indices->getNDMRG( irrep );
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ];
      int shift_singlet = 0;
      int shift_triplet = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int NVIR_a = indices->getNVIRT( irrep_a );
         const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
         const int irrep_occ = Irreps::directProd( irrep_a, irrep );
         if ( irrep_occ == 0 ){
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               assert( shift_singlet == shift_E_nonactive( indices, irrep_a, irrep_ij, irrep_ij, +1 ) );
               assert( shift_triplet == shift_E_nonactive( indices, irrep_a, irrep_ij, irrep_ij, -1 ) );
               const int NOCC_ij = indices->getNOCC( irrep_ij );
               for ( int i = 0; i < NOCC_ij; i++ ){
                  for ( int j = i; j < NOCC_ij; j++ ){
                     const double ij_factor = (( i == j ) ? SQRT_0p5 : 1.0 );
                     for ( int a = 0; a < NVIR_a; a++ ){
                        const int count_aij_singlet = shift_singlet + a + NVIR_a * ( i + (j*(j+1))/2 );
                        const int count_aij_triplet = shift_triplet + a + NVIR_a * ( i + (j*(j-1))/2 );
                        for ( int t = 0; t < num_t; t++ ){
                           double value_singlet = 0.0;
                           double value_triplet = 0.0;
                           for ( int w = 0; w < num_t; w++ ){
                              const double SEE_wt = SEE[ irrep ][ w + num_t * t ];
                              const double aj_wi  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, i, occ_t + w, j, N_OA_a + a );
                              const double ai_wj  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, j, occ_t + w, i, N_OA_a + a );
                              value_singlet +=     SEE_wt * ( aj_wi + ai_wj );
                              value_triplet += 3 * SEE_wt * ( aj_wi - ai_wj );
                           }
                           target_singlet[ t + num_t * count_aij_singlet ] = value_singlet * ij_factor;
             if ( j > i ){ target_triplet[ t + num_t * count_aij_triplet ] = value_triplet; }
                        }
                     }
                  }
               }
               shift_singlet += ( NVIR_a * NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               shift_triplet += ( NVIR_a * NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
            }
         } else {
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_i, irrep_occ );
               if ( irrep_i < irrep_j ){
                  assert( shift_singlet == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, +1 ) );
                  assert( shift_triplet == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, -1 ) );
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const int count_aij_singlet = shift_singlet + a + NVIR_a * ( i + NOCC_i * j );
                           const int count_aij_triplet = shift_triplet + a + NVIR_a * ( i + NOCC_i * j );
                           for ( int t = 0; t < num_t; t++ ){
                              double value_singlet = 0.0;
                              double value_triplet = 0.0;
                              for ( int w = 0; w < num_t; w++ ){
                                 const double SEE_wt = SEE[ irrep ][ w + num_t * t ];
                                 const double aj_wi  = integrals->get_coulomb( irrep_i, irrep, irrep_j, irrep_a, i, occ_t + w, j, N_OA_a + a );
                                 const double ai_wj  = integrals->get_coulomb( irrep_j, irrep, irrep_i, irrep_a, j, occ_t + w, i, N_OA_a + a );
                                 value_singlet +=     SEE_wt * ( aj_wi + ai_wj );
                                 value_triplet += 3 * SEE_wt * ( aj_wi - ai_wj );
                              }
                              target_singlet[ t + num_t * count_aij_singlet ] = value_singlet;
                              target_triplet[ t + num_t * count_aij_triplet ] = value_triplet;
                           }
                        }
                     }
                  }
                  shift_singlet += NVIR_a * NOCC_i * NOCC_j;
                  shift_triplet += NVIR_a * NOCC_i * NOCC_j;
               }
            }
         }
      }
      assert( shift_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] );
      assert( shift_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] );
   }

   // VG singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int occ_t = indices->getNOCC( irrep );
      const int num_t = indices->getNDMRG( irrep );
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ];
      int shift_singlet = 0;
      int shift_triplet = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int NOCC_i = indices->getNOCC( irrep_i );
         const int irrep_virt = Irreps::directProd( irrep_i, irrep );
         if ( irrep_virt == 0 ){
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               assert( shift_singlet == shift_G_nonactive( indices, irrep_i, irrep_ab, irrep_ab, +1 ) );
               assert( shift_triplet == shift_G_nonactive( indices, irrep_i, irrep_ab, irrep_ab, -1 ) );
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               for ( int i = 0; i < NOCC_i; i++ ){
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const double ab_factor = (( a == b ) ? SQRT_0p5 : 1.0 );
                        const int count_abi_singlet = shift_singlet + i + NOCC_i * ( a + (b*(b+1))/2 );
                        const int count_abi_triplet = shift_triplet + i + NOCC_i * ( a + (b*(b-1))/2 );
                        for ( int t = 0; t < num_t; t++ ){
                           double value_singlet = 0.0;
                           double value_triplet = 0.0;
                           for ( int u = 0; u < num_t; u++ ){
                              const double SGG_ut = SGG[ irrep ][ u + num_t * t ];
                              const double ai_bu  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + a, N_OA_ab + b );
                              const double bi_au  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + b, N_OA_ab + a );
                              value_singlet +=     SGG_ut * ( ai_bu + bi_au );
                              value_triplet += 3 * SGG_ut * ( ai_bu - bi_au );
                           }
                           target_singlet[ t + num_t * count_abi_singlet ] = value_singlet * ab_factor;
             if ( b > a ){ target_triplet[ t + num_t * count_abi_triplet ] = value_triplet; }
                        }
                     }
                  }
               }
               shift_singlet += ( NOCC_i * NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               shift_triplet += ( NOCC_i * NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            }
         } else {
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_a, irrep_virt );
               if ( irrep_a < irrep_b ){
                  assert( shift_singlet == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, +1 ) );
                  assert( shift_triplet == shift_G_nonactive( indices, irrep_i, irrep_a, irrep_b, -1 ) );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int a = 0; a < NVIR_a; a++ ){
                        for ( int b = 0; b < NVIR_b; b++ ){
                           const int count_abi_singlet = shift_singlet + i + NOCC_i * ( a + NVIR_a * b );
                           const int count_abi_triplet = shift_triplet + i + NOCC_i * ( a + NVIR_a * b );
                           for ( int t = 0; t < num_t; t++ ){
                              double value_singlet = 0.0;
                              double value_triplet = 0.0;
                              for ( int u = 0; u < num_t; u++ ){
                                 const double SGG_ut = SGG[ irrep ][ u + num_t * t ];
                                 const double ai_bu  = integrals->get_exchange( irrep_i, irrep, irrep_a, irrep_b, i, occ_t + u, N_OA_a + a, N_OA_b + b );
                                 const double bi_au  = integrals->get_exchange( irrep_i, irrep, irrep_b, irrep_a, i, occ_t + u, N_OA_b + b, N_OA_a + a );
                                 value_singlet +=     SGG_ut * ( ai_bu + bi_au );
                                 value_triplet += 3 * SGG_ut * ( ai_bu - bi_au );
                              }
                              target_singlet[ t + num_t * count_abi_singlet ] = value_singlet;
                              target_triplet[ t + num_t * count_abi_triplet ] = value_triplet;
                           }
                        }
                     }
                  }
                  shift_singlet += NOCC_i * NVIR_a * NVIR_b;
                  shift_triplet += NOCC_i * NVIR_a * NVIR_b;
               }
            }
         }
      }
      assert( shift_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] );
      assert( shift_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] );
   }

   // VH singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int shift_singlet = 0;
      int shift_triplet = 0;
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ];
      if ( irrep == 0 ){ // irrep_i == irrep_j  and  irrep_a == irrep_b
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            const int linsize_singlet = ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
            const int linsize_triplet = ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               assert( shift_singlet == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, +1 ) );
               assert( shift_triplet == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, -1 ) );
               const int nvirt_ab = indices->getNVIRT( irrep_ab );
               const int noa_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               for ( int a = 0; a < nvirt_ab; a++ ){
                  for ( int b = a; b < nvirt_ab; b++ ){
                     const double ab_factor = (( a==b ) ? SQRT_0p5 : 1.0 );
                     for ( int i = 0; i < nocc_ij; i++ ){
                        for ( int j = i; j < nocc_ij; j++){
                           const double ij_factor = (( i==j ) ? SQRT_0p5 : 1.0 );
                           const int count_singlet = shift_singlet + i + (j*(j+1))/2 + linsize_singlet * ( a + (b*(b+1))/2 );
                           const int count_triplet = shift_triplet + i + (j*(j-1))/2 + linsize_triplet * ( a + (b*(b-1))/2 );
                           const double ai_bj = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, i, j, noa_ab + a, noa_ab + b );
                           const double aj_bi = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, j, i, noa_ab + a, noa_ab + b );
                           target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi ) * ij_factor * ab_factor;
   if ( (b-a)*(j-i) > 0 ){ target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi ); }
                        }
                     }
                  }
               }
               shift_singlet += linsize_singlet * ( ( nvirt_ab * ( nvirt_ab + 1 ) ) / 2 );
               shift_triplet += linsize_triplet * ( ( nvirt_ab * ( nvirt_ab - 1 ) ) / 2 );
            }
         }
      } else { // irrep_i < irrep_j = irrep_i x irrep   and   irrep_a < irrep_b = irrep_a x irrep
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
                  const int irrep_b = Irreps::directProd( irrep, irrep_a );
                  if ( irrep_a < irrep_b ){
                     assert( shift_singlet == shift_H_nonactive( indices, irrep_i, irrep_j, irrep_a, irrep_b, +1 ) );
                     assert( shift_triplet == shift_H_nonactive( indices, irrep_i, irrep_j, irrep_a, irrep_b, -1 ) );
                     const int nvir_a = indices->getNVIRT( irrep_a );
                     const int nvir_b = indices->getNVIRT( irrep_b );
                     const int noa_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                     const int noa_b  = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                     for ( int a = 0; a < nvir_a; a++ ){
                        for ( int b = 0; b < nvir_b; b++ ){
                           for ( int i = 0; i < nocc_i; i++ ){
                              for ( int j = 0; j < nocc_j; j++){
                                 const int count_singlet = shift_singlet + i + nocc_i * ( j + nocc_j * ( a + nvir_a * b ) );
                                 const int count_triplet = shift_triplet + i + nocc_i * ( j + nocc_j * ( a + nvir_a * b ) );
                                 const double ai_bj = integrals->get_exchange( irrep_i, irrep_j, irrep_a, irrep_b, i, j, noa_a + a, noa_b + b );
                                 const double aj_bi = integrals->get_exchange( irrep_j, irrep_i, irrep_a, irrep_b, j, i, noa_a + a, noa_b + b );
                                 target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi );
                                 target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi );
                              }
                           }
                        }
                     }
                     shift_singlet += nocc_i * nocc_j * nvir_a * nvir_b;
                     shift_triplet += nocc_i * nocc_j * nvir_a * nvir_b;
                  }
               }
            }
         }
      }
      assert( shift_singlet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] );
      assert( shift_triplet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] );
   }

}

int CheMPS2::CASPT2::jump_AC_active( const DMRGSCFindices * idx, const int irrep_t, const int irrep_u, const int irrep_v ){

   const int irrep_tuv = Irreps::directProd( irrep_t, Irreps::directProd( irrep_u, irrep_v ) );
   const int n_irreps = idx->getNirreps();

   int jump_ac = 0;
   for ( int It = 0; It < n_irreps; It++ ){
      for ( int Iu = 0; Iu < n_irreps; Iu++ ){
         const int Iv = Irreps::directProd( irrep_tuv, Irreps::directProd( It, Iu ) );
         if (( It == irrep_t ) && ( Iu == irrep_u ) && ( Iv == irrep_v )){
            It = n_irreps;
            Iu = n_irreps;
         } else {
            jump_ac += idx->getNDMRG( It ) * idx->getNDMRG( Iu ) * idx->getNDMRG( Iv );
         }
      }
   }
   return jump_ac;

}

int CheMPS2::CASPT2::jump_BF_active( const DMRGSCFindices * idx, const int irrep_t, const int irrep_u, const int ST ){

   assert( irrep_t <= irrep_u );
   const int irrep_tu = Irreps::directProd( irrep_t, irrep_u );
   const int n_irreps = idx->getNirreps();

   int jump_bf = 0;
   if ( irrep_tu == 0 ){
      for ( int Itu = 0; Itu < n_irreps; Itu++ ){
         if (( Itu == irrep_t ) && ( Itu == irrep_u )){
            Itu = n_irreps;
         } else {
            jump_bf += ( idx->getNDMRG( Itu ) * ( idx->getNDMRG( Itu ) + ST ) ) / 2;
         }
      }
   } else {
      for ( int It = 0; It < n_irreps; It++ ){
         const int Iu = Irreps::directProd( irrep_tu, It );
         if ( It < Iu ){
            if (( It == irrep_t ) && ( Iu == irrep_u )){
               It = n_irreps;
            } else {
               jump_bf += idx->getNDMRG( It ) * idx->getNDMRG( Iu );
            }
         }
      }
   }
   return jump_bf;

}

void CheMPS2::CASPT2::make_FDE_FDG(){

   /*
      FD1E singlet: < E_yx E_lb E_kw SE_tiaj > = 1 delta_ab ( delta_ik delta_jl + delta_il delta_jk ) / sqrt( 1 + delta_ij ) FD1E_singlet[ Ib x Il ][ It ][ w ][ xy, t ]
      FD2E singlet: < E_yb E_lx E_kw SE_tiaj > = 1 delta_ab ( delta_ik delta_jl + delta_il delta_jk ) / sqrt( 1 + delta_ij ) FD2E_singlet[ Ib x Il ][ It ][ w ][ xy, t ]
      FD1E triplet: < E_yx E_lb E_kw TE_tiaj > = 3 delta_ab ( delta_ik delta_jl - delta_il delta_jk ) / sqrt( 1 + delta_ij ) FD1E_triplet[ Ib x Il ][ It ][ w ][ xy, t ]
      FD2E triplet: < E_yb E_lx E_kw TE_tiaj > = 3 delta_ab ( delta_ik delta_jl - delta_il delta_jk ) / sqrt( 1 + delta_ij ) FD2E_triplet[ Ib x Il ][ It ][ w ][ xy, t ]

            FD1E_singlet[ Ib x Il ][ It ][ w ][ xy, t ] = ( + 2 delta_tw Gamma_yx
                                                            - delta_tx Gamma_yw
                                                            - Gamma_ytxw
                                                          )

            FD2E_singlet[ Ib x Il ][ It ][ w ][ xy, t ] = ( + Gamma_ytxw
                                                            + Gamma_ytwx
                                                            - delta_tx Gamma_yw
                                                            - delta_tw Gamma_yx
                                                          )

            FD1E_triplet[ Ib x Il ][ It ][ w ][ xy, t ] = ( + 2 delta_tw Gamma_yx
                                                            - delta_tx Gamma_yw
                                                            - Gamma_ytxw
                                                          )

            FD2E_triplet[ Ib x Il ][ It ][ w ][ xy, t ] = ( + Gamma_ytxw / 3
                                                            - Gamma_ytwx / 3
                                                            + delta_tx Gamma_yw
                                                            - delta_tw Gamma_yx
                                                          )

      FD1G singlet: < E_yx E_jd E_wc SG_aibt > = 1 delta_ij ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FD1G_singlet[ Ij x Id ][ It ][ w ][ xy, t ]
      FD2G singlet: < E_yd E_jx E_wc SG_aibt > = 1 delta_ij ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FD2G_singlet[ Ij x Id ][ It ][ w ][ xy, t ]
      FD1G triplet: < E_yx E_jd E_wc TG_aibt > = 3 delta_ij ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FD1G_triplet[ Ij x Id ][ It ][ w ][ xy, t ]
      FD2G triplet: < E_yd E_jx E_wc TG_aibt > = 3 delta_ij ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FD2G_triplet[ Ij x Id ][ It ][ w ][ xy, t ]

            FD1G_singlet[ Ij x Id ][ It ][ w ][ xy, t ] = ( + Gamma_ywxt
                                                            + delta_wx Gamma_yt
                                                          )

            FD2G_singlet[ Ij x Id ][ It ][ w ][ xy, t ] = ( + delta_xw Gamma_yt
                                                            - Gamma_ywtx
                                                            - Gamma_ywxt
                                                          )

            FD1G_triplet[ Ij x Id ][ It ][ w ][ xy, t ] = ( - Gamma_ywxt
                                                            - delta_wx Gamma_yt
                                                          )

            FD2G_triplet[ Ij x Id ][ It ][ w ][ xy, t ] = ( + delta_xw Gamma_yt
                                                            - Gamma ywtx / 3
                                                            + Gamma ywxt / 3
                                                          )
   */

   FDE_singlet = new double***[ num_irreps ];
   FDE_triplet = new double***[ num_irreps ];
   FDG_singlet = new double***[ num_irreps ];
   FDG_triplet = new double***[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      const int SIZE_left = size_D[ irrep_left ];
      const int D2JUMP    = SIZE_left / 2;
      FDE_singlet[ irrep_left ] = new double**[ num_irreps ];
      FDE_triplet[ irrep_left ] = new double**[ num_irreps ];
      FDG_singlet[ irrep_left ] = new double**[ num_irreps ];
      FDG_triplet[ irrep_left ] = new double**[ num_irreps ];

      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){

         const int SIZE_right = indices->getNDMRG( irrep_t );
         const int d_t     = indices->getDMRGcumulative( irrep_t );
         const int num_t   = indices->getNDMRG( irrep_t );
         const int irrep_w = Irreps::directProd( irrep_left, irrep_t );
         const int d_w     = indices->getDMRGcumulative( irrep_w );
         const int num_w   = indices->getNDMRG( irrep_w );
         FDE_singlet[ irrep_left ][ irrep_t ] = new double*[ num_w ];
         FDE_triplet[ irrep_left ][ irrep_t ] = new double*[ num_w ];
         FDG_singlet[ irrep_left ][ irrep_t ] = new double*[ num_w ];
         FDG_triplet[ irrep_left ][ irrep_t ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FDE_singlet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];
            FDE_triplet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];
            FDG_singlet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];
            FDG_triplet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * FDE_sing = FDE_singlet[ irrep_left ][ irrep_t ][ w ];
            double * FDE_trip = FDE_triplet[ irrep_left ][ irrep_t ][ w ];
            double * FDG_sing = FDG_singlet[ irrep_left ][ irrep_t ][ w ];
            double * FDG_trip = FDG_triplet[ irrep_left ][ irrep_t ][ w ];

            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int d_x     = indices->getDMRGcumulative( irrep_x );
               const int num_x   = indices->getNDMRG( irrep_x );
               const int irrep_y = Irreps::directProd( irrep_left, irrep_x );
               const int d_y     = indices->getDMRGcumulative( irrep_y );
               const int num_y   = indices->getNDMRG( irrep_y );

               for ( int t = 0; t < num_t; t++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double gamma_ytxw = two_rdm[ d_y + y + LAS * ( d_t + t + LAS * ( d_x + x + LAS * ( d_w + w ))) ];
                        const double gamma_ytwx = two_rdm[ d_y + y + LAS * ( d_t + t + LAS * ( d_w + w + LAS * ( d_x + x ))) ];
                        FDE_sing[          jump_row + x + num_x * y + SIZE_left * t ] = - gamma_ytxw;
                        FDE_sing[ D2JUMP + jump_row + x + num_x * y + SIZE_left * t ] = + gamma_ytxw + gamma_ytwx;
                        FDE_trip[          jump_row + x + num_x * y + SIZE_left * t ] = - gamma_ytxw;
                        FDE_trip[ D2JUMP + jump_row + x + num_x * y + SIZE_left * t ] = ( gamma_ytxw - gamma_ytwx ) / 3.0;
                        const double gamma_ywtx = two_rdm[ d_y + y + LAS * ( d_w + w + LAS * ( d_t + t + LAS * ( d_x + x ))) ];
                        const double gamma_ywxt = two_rdm[ d_y + y + LAS * ( d_w + w + LAS * ( d_x + x + LAS * ( d_t + t ))) ];
                        FDG_sing[          jump_row + x + num_x * y + SIZE_left * t ] = + gamma_ywxt;
                        FDG_sing[ D2JUMP + jump_row + x + num_x * y + SIZE_left * t ] = - gamma_ywxt - gamma_ywtx;
                        FDG_trip[          jump_row + x + num_x * y + SIZE_left * t ] = - gamma_ywxt;
                        FDG_trip[ D2JUMP + jump_row + x + num_x * y + SIZE_left * t ] = ( gamma_ywxt - gamma_ywtx ) / 3.0;
                     }
                  }
               }

               if (( irrep_t == irrep_w ) && ( irrep_x == irrep_y )){
                  for ( int x = 0; x < num_x; x++ ){
                     for ( int y = 0; y < num_x; y++ ){
                        const double gamma_yx = one_rdm[ d_y + y + LAS * ( d_x + x ) ];
                        FDE_sing[          jump_row + x + num_x * y + SIZE_left * w ] += 2 * gamma_yx;
                        FDE_sing[ D2JUMP + jump_row + x + num_x * y + SIZE_left * w ] -= gamma_yx;
                        FDE_trip[          jump_row + x + num_x * y + SIZE_left * w ] += 2 * gamma_yx;
                        FDE_trip[ D2JUMP + jump_row + x + num_x * y + SIZE_left * w ] -= gamma_yx;
                     }
                  }
               }

               if (( irrep_t == irrep_x ) && ( irrep_y == irrep_w )){
                  for ( int y = 0; y < num_y; y++ ){
                     const double gamma_yw = one_rdm[ d_y + y + LAS * ( d_w + w ) ];
                     for ( int tx = 0; tx < num_x; tx++ ){
                        FDE_sing[          jump_row + tx + num_x * y + SIZE_left * tx ] -= gamma_yw;
                        FDE_sing[ D2JUMP + jump_row + tx + num_x * y + SIZE_left * tx ] -= gamma_yw;
                        FDE_trip[          jump_row + tx + num_x * y + SIZE_left * tx ] -= gamma_yw;
                        FDE_trip[ D2JUMP + jump_row + tx + num_x * y + SIZE_left * tx ] += gamma_yw;
                     }
                  }
               }

               if (( irrep_t == irrep_y ) && ( irrep_w == irrep_x )){
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        const double gamma_yt = one_rdm[ d_y + y + LAS * ( d_t + t ) ];
                        FDG_sing[          jump_row + w + num_x * y + SIZE_left * t ] += gamma_yt;
                        FDG_sing[ D2JUMP + jump_row + w + num_x * y + SIZE_left * t ] += gamma_yt;
                        FDG_trip[          jump_row + w + num_x * y + SIZE_left * t ] -= gamma_yt;
                        FDG_trip[ D2JUMP + jump_row + w + num_x * y + SIZE_left * t ] += gamma_yt;
                     }
                  }
               }
               jump_row += num_x * num_y;
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FEH_FGH(){

   /*
      FEH singlet: < SE_xkdl E_wc SH_aibj > = 2 delta_ik delta_jl ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FEH[ Ix ][ w ][ x ]
      FEH triplet: < TE_xkdl E_wc TH_aibj > = 6 delta_ik delta_jl ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FEH[ Ix ][ w ][ x ]
      FGH singlet: < SG_cldx E_kw SH_aibj > = 2 delta_ac delta_bd ( delta_il delta_jk + delta_ik delta_jl ) / sqrt( 1 + delta_ij ) FGH[ Ix ][ w ][ x ]
      FGH triplet: < TG_cldx E_kw TH_aibj > = 6 delta_ac delta_bd ( delta_il delta_jk - delta_ik delta_jl ) / sqrt( 1 + delta_ij ) FGH[ Ix ][ w ][ x ]

            FEH[ Ix ][ w ][ x ] = + SEE[ Ix ][ xw ]
            FGH[ Ix ][ w ][ x ] = - SGG[ Ix ][ xw ]
   */

   FEH = new double**[ num_irreps ];
   FGH = new double**[ num_irreps ];

   for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){

      const int NACT = indices->getNDMRG( irrep_x );
      FEH[ irrep_x ] = new double*[ NACT ];
      FGH[ irrep_x ] = new double*[ NACT ];

      for ( int w = 0; w < NACT; w++ ){
         FEH[ irrep_x ][ w ] = new double[ NACT ];
         FGH[ irrep_x ][ w ] = new double[ NACT ];

         for ( int x = 0; x < NACT; x++ ){
            FEH[ irrep_x ][ w ][ x ] = + SEE[ irrep_x ][ x + NACT * w ];
            FGH[ irrep_x ][ w ][ x ] = - SGG[ irrep_x ][ x + NACT * w ];
         }
      }
   }

}

void CheMPS2::CASPT2::make_FBE_FFG_singlet(){

   /*
      FBE singlet : < SB_xkyl E_wc SE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_singlet[ Ik x Il ][ It ][ w ][ xy, t ]
      FFG singlet : < SF_cxdy E_kw SG_aibt > = 2 delta_ac delta_bd delta_ik FFG_singlet[ Ic x Id ][ It ][ w ][ xy, t ]

            FBE_singlet[ Ik x Il ][ It ][ w ][ xy, t ] = + SBB_singlet[ Ik x Il ][ xy, tw ]
            FFG_singlet[ Ic x Id ][ It ][ w ][ xy, t ] = - SFF_singlet[ Ic x Id ][ xy, tw ]
   */

   FBE_singlet = new double***[ num_irreps ];
   FFG_singlet = new double***[ num_irreps ];

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      assert( size_B_singlet[ irrep_left ] == size_F_singlet[ irrep_left ] ); // At construction
      const int SIZE_left = size_B_singlet[ irrep_left ];
      FBE_singlet[ irrep_left ] = new double**[ num_irreps ];
      FFG_singlet[ irrep_left ] = new double**[ num_irreps ];

      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){

         const int SIZE_right = indices->getNDMRG( irrep_t );
         const int num_t   = indices->getNDMRG( irrep_t );
         const int irrep_w = Irreps::directProd( irrep_left, irrep_t );
         const int num_w   = indices->getNDMRG( irrep_w );
         FBE_singlet[ irrep_left ][ irrep_t ] = new double*[ num_w ];
         FFG_singlet[ irrep_left ][ irrep_t ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FBE_singlet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];
            FFG_singlet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * BEptr = FBE_singlet[ irrep_left ][ irrep_t ][ w ];
            double * FGptr = FFG_singlet[ irrep_left ][ irrep_t ][ w ];

            if ( irrep_left == 0 ){ // irrep_t == irrep_w
               const int jump_singlet = jump_BF_active( indices, irrep_t, irrep_w, +1 );
               for ( int t = 0; t < w; t++ ){
                  for ( int xy = 0; xy < SIZE_left; xy++ ){
                     BEptr[ xy + SIZE_left * t ] = + SBB_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + t + ( w * ( w + 1 ) ) / 2 ) ];
                     FGptr[ xy + SIZE_left * t ] = - SFF_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + t + ( w * ( w + 1 ) ) / 2 ) ];
                  }
               }
               for ( int t = w; t < num_t; t++ ){
                  for ( int xy = 0; xy < SIZE_left; xy++ ){
                     BEptr[ xy + SIZE_left * t ] = + SBB_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + w + ( t * ( t + 1 ) ) / 2 ) ];
                     FGptr[ xy + SIZE_left * t ] = - SFF_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + w + ( t * ( t + 1 ) ) / 2 ) ];
                  }
               }
            } else { // irrep_t != irrep_w
               if ( irrep_t < irrep_w ){
                  const int jump_singlet = jump_BF_active( indices, irrep_t, irrep_w, +1 );
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int xy = 0; xy < SIZE_left; xy++ ){
                        BEptr[ xy + SIZE_left * t ] = + SBB_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + t + num_t * w ) ];
                        FGptr[ xy + SIZE_left * t ] = - SFF_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + t + num_t * w ) ];
                     }
                  }
               } else {
                  const int jump_singlet = jump_BF_active( indices, irrep_w, irrep_t, +1 );
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int xy = 0; xy < SIZE_left; xy++ ){
                        BEptr[ xy + SIZE_left * t ] = + SBB_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + w + num_w * t ) ];
                        FGptr[ xy + SIZE_left * t ] = - SFF_singlet[ irrep_left ][ xy + SIZE_left * ( jump_singlet + w + num_w * t ) ];
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FBE_FFG_triplet(){

   /*
      FBE triplet : < TB_xkyl E_wc TE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_triplet[ Ik x Il ][ It ][ w ][ xy, t ]
      FFG triplet : < TF_cxdy E_kw TG_aibt > = 2 delta_ac delta_bd delta_ik FFG_triplet[ Ic x Id ][ It ][ w ][ xy, t ]

            FBE_triplet[ Ik x Il ][ It ][ w ][ xy, t ] = + SBB_triplet[ Ik x Il ][ xy, tw ]
            FFG_triplet[ Ic x Id ][ It ][ w ][ xy, t ] = + SFF_triplet[ Ic x Id ][ xy, tw ]
   */

   FBE_triplet = new double***[ num_irreps ];
   FFG_triplet = new double***[ num_irreps ];

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      assert( size_B_triplet[ irrep_left ] == size_F_triplet[ irrep_left ] ); // At construction
      const int SIZE_left = size_B_triplet[ irrep_left ];
      FBE_triplet[ irrep_left ] = new double**[ num_irreps ];
      FFG_triplet[ irrep_left ] = new double**[ num_irreps ];

      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){

         const int SIZE_right = indices->getNDMRG( irrep_t );
         const int num_t   = indices->getNDMRG( irrep_t );
         const int irrep_w = Irreps::directProd( irrep_left, irrep_t );
         const int num_w   = indices->getNDMRG( irrep_w );
         FBE_triplet[ irrep_left ][ irrep_t ] = new double*[ num_w ];
         FFG_triplet[ irrep_left ][ irrep_t ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FBE_triplet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];
            FFG_triplet[ irrep_left ][ irrep_t ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * BEptr = FBE_triplet[ irrep_left ][ irrep_t ][ w ];
            double * FGptr = FFG_triplet[ irrep_left ][ irrep_t ][ w ];

            if ( irrep_left == 0 ){ // irrep_t == irrep_w
               const int jump_triplet = jump_BF_active( indices, irrep_t, irrep_w, -1 );
               for ( int t = 0; t < w; t++ ){
                  for ( int xy = 0; xy < SIZE_left; xy++ ){
                     BEptr[ xy + SIZE_left * t ] = + SBB_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + t + ( w * ( w - 1 ) ) / 2 ) ];
                     FGptr[ xy + SIZE_left * t ] = + SFF_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + t + ( w * ( w - 1 ) ) / 2 ) ];
                  }
               }
               for ( int xy = 0; xy < SIZE_left; xy++ ){
                  BEptr[ xy + SIZE_left * w ] = 0.0;
                  FGptr[ xy + SIZE_left * w ] = 0.0;
               }
               for ( int t = w+1; t < num_t; t++ ){
                  for ( int xy = 0; xy < SIZE_left; xy++ ){
                     BEptr[ xy + SIZE_left * t ] = - SBB_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + w + ( t * ( t - 1 ) ) / 2 ) ];
                     FGptr[ xy + SIZE_left * t ] = - SFF_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + w + ( t * ( t - 1 ) ) / 2 ) ];
                  }
               }
            } else { // irrep_t != irrep_w
               if ( irrep_t < irrep_w ){
                  const int jump_triplet = jump_BF_active( indices, irrep_t, irrep_w, -1 );
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int xy = 0; xy < SIZE_left; xy++ ){
                        BEptr[ xy + SIZE_left * t ] = + SBB_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + t + num_t * w ) ];
                        FGptr[ xy + SIZE_left * t ] = + SFF_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + t + num_t * w ) ];
                     }
                  }
               } else {
                  const int jump_triplet = jump_BF_active( indices, irrep_w, irrep_t, -1 );
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int xy = 0; xy < SIZE_left; xy++ ){
                        BEptr[ xy + SIZE_left * t ] = - SBB_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + w + num_w * t ) ];
                        FGptr[ xy + SIZE_left * t ] = - SFF_triplet[ irrep_left ][ xy + SIZE_left * ( jump_triplet + w + num_w * t ) ];
                     }
                  }
               }
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FAB_FCF_singlet(){

   /*
      FAB singlet : < E_zy E_lx | E_kw | SB_tiuj > = ( delta_ik delta_jl + delta_jk delta_il ) / sqrt( 1 + delta_ij ) * FAB_singlet[ Il ][ Ii x Ij ][ w ][ xyz, tu ]

      FCF singlet : < E_zy E_xd | E_wc | SF_atbu > = ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) * FCF_singlet[ Id ][ Ia x Ib ][ w ][ xyz, tu ]

            FAB_singlet[ Il ][ Ii x Ij ][ w ][ xyz, tu ] = ( + 2 delta_tw delta_ux Gamma_zy
                                                             + 2 delta_uw delta_tx Gamma_zy
                                                             - delta_tw Gamma_zuyx
                                                             - delta_tw delta_uy Gamma_zx
                                                             - delta_uw Gamma_ztyx
                                                             - delta_uw delta_ty Gamma_zx
                                                             - SAA[ Il ][ xyz, utw ]
                                                             - SAA[ Il ][ xyz, tuw ]
                                                           )

            FCF_singlet[ Id ][ Ia x Ib ][ w ][ xyz, tu ] = ( + SCC[ Id ][ xyz, uwt ]
                                                             + SCC[ Id ][ xyz, twu ]
                                                             - delta_uw Gamma_zxyt
                                                             - delta_uw delta_xy Gamma_zt
                                                             - delta_tw Gamma_zxyu
                                                             - delta_tw delta_xy Gamma_zu
                                                           )
   */

   FAB_singlet = new double***[ num_irreps ];
   FCF_singlet = new double***[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      assert( size_A[ irrep_left ] == size_C[ irrep_left ] ); // At construction
      const int SIZE_left = size_A[ irrep_left ];
      FAB_singlet[ irrep_left ] = new double**[ num_irreps ];
      FCF_singlet[ irrep_left ] = new double**[ num_irreps ];

      { // irrep_right == 0

         assert( size_B_singlet[ 0 ] == size_F_singlet[ 0 ] ); // At construction
         const int SIZE_right = size_B_singlet[ 0 ];
         const int num_w = indices->getNDMRG( irrep_left ); // irrep_w == irrep_left x irrep_right = irrep_left
         FAB_singlet[ irrep_left ][ 0 ] = new double*[ num_w ];
         FCF_singlet[ irrep_left ][ 0 ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FAB_singlet[ irrep_left ][ 0 ][ w ] = new double[ SIZE_left * SIZE_right ];
            FCF_singlet[ irrep_left ][ 0 ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * ABptr = FAB_singlet[ irrep_left ][ 0 ][ w ];
            double * CFptr = FCF_singlet[ irrep_left ][ 0 ][ w ];

            int jump_col = 0;
            for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
               const int d_ut   = indices->getDMRGcumulative( irrep_ut );
               const int num_ut = indices->getNDMRG( irrep_ut );
               assert( jump_col == jump_BF_active( indices, irrep_ut, irrep_ut, +1 ) );

               const int jump_AB = jump_AC_active( indices, irrep_ut, irrep_ut, irrep_left );
               const int jump_CF = jump_AC_active( indices, irrep_ut, irrep_left, irrep_ut );

               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                     for ( int xyz = 0; xyz < SIZE_left; xyz++ ){
                        ABptr[ xyz + SIZE_left * ( jump_col + t + ( u * ( u + 1 ) ) / 2 ) ] = ( - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB + u + num_ut * ( t + num_ut * w )) ]
                                                                                                - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB + t + num_ut * ( u + num_ut * w )) ] );
                        CFptr[ xyz + SIZE_left * ( jump_col + t + ( u * ( u + 1 ) ) / 2 ) ] = ( + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF + u + num_ut * ( w + num_w  * t )) ]
                                                                                                + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF + t + num_ut * ( w + num_w  * u )) ] );
                     }
                  }
               }

               int jump_row = 0;
               for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                  const int d_x   = indices->getDMRGcumulative( irrep_x );
                  const int num_x = indices->getNDMRG( irrep_x );
                  for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                     const int irrep_z = Irreps::directProd( Irreps::directProd( irrep_left, irrep_x ), irrep_y );
                     const int d_y   = indices->getDMRGcumulative( irrep_y );
                     const int num_y = indices->getNDMRG( irrep_y );
                     const int d_z   = indices->getDMRGcumulative( irrep_z );
                     const int num_z = indices->getNDMRG( irrep_z );
                     assert( jump_row == jump_AC_active( indices, irrep_x, irrep_y, irrep_z ) );

                     if ( irrep_ut == irrep_left ){

                        // FAB_singlet[ xyz,tu ] -= delta_tw Gamma_zuyx
                        for ( int u = w; u < num_ut; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zuyx = two_rdm[ d_z + z + LAS * ( d_ut + u + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                    ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( u * ( u + 1 ) ) / 2 ) ] -= gamma_zuyx;
                                 }
                              }
                           }
                        }

                        // FAB_singlet[ xyz,tu ] -= delta_tw delta_uy Gamma_zx
                        if ( irrep_ut == irrep_y ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                 for ( int uy = w; uy < num_y; uy++ ){
                                    ABptr[ jump_row + x + num_x * ( uy + num_y * z ) + SIZE_left * ( jump_col + w + ( uy * ( uy + 1 ) ) / 2 ) ] -= gamma_zx;
                                 }
                              }
                           }
                        }

                        // FAB_singlet[ xyz,tu ] += 2 delta_tw delta_ux Gamma_zy
                        if ( irrep_ut == irrep_x ){
                           for ( int y = 0; y < num_y; y++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                 for ( int ux = w; ux < num_x; ux++ ){
                                    ABptr[ jump_row + ux + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( ux * ( ux + 1 ) ) / 2 ) ] += 2 * gamma_zy;
                                 }
                              }
                           }
                        }

                        // FCF_singlet[ xyz,tu ] -= delta_tw Gamma_zxyu
                        for ( int u = w; u < num_ut; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyu = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_ut + u ))) ];
                                    CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( u * ( u + 1 ) ) / 2 ) ] -= gamma_zxyu;
                                 }
                              }
                           }
                        }

                        // FCF_singlet[ xyz,tu ] -= delta_tw delta_xy Gamma_zu
                        if ( irrep_x == irrep_y ){
                           for ( int u = w; u < num_ut; u++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zu = one_rdm[ d_z + z + LAS * ( d_ut + u ) ];
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + w + ( u * ( u + 1 ) ) / 2 ) ] -= gamma_zu;
                                 }
                              }
                           }
                        }

                        // FAB_singlet[ xyz,tu ] -= delta_uw Gamma_ztyx
                        for ( int t = 0; t <= w; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_ztyx = two_rdm[ d_z + z + LAS * ( d_ut + t + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                    ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + ( w * ( w + 1 ) ) / 2 ) ] -= gamma_ztyx;
                                 }
                              }
                           }
                        }

                        // FAB_singlet[ xyz,tu ] -= delta_uw delta_ty Gamma_zx
                        if ( irrep_ut == irrep_y ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                 for ( int ty = 0; ty <= w; ty++ ){
                                    ABptr[ jump_row + x + num_x * ( ty + num_y * z ) + SIZE_left * ( jump_col + ty + ( w * ( w + 1 ) ) / 2 ) ] -= gamma_zx;
                                 }
                              }
                           }
                        }

                        // FAB_singlet[ xyz,tu ] += 2 delta_uw delta_tx Gamma_zy
                        if ( irrep_ut == irrep_x ){
                           for ( int y = 0; y < num_y; y++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                 for ( int tx = 0; tx <= w; tx++ ){
                                    ABptr[ jump_row + tx + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + tx + ( w * ( w + 1 ) ) / 2 ) ] += 2 * gamma_zy;
                                 }
                              }
                           }
                        }

                        // FCF_singlet[ xyz,tu ] -= delta_uw Gamma_zxyt
                        for ( int t = 0; t <= w; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyt = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_ut + t ))) ];
                                    CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + ( w * ( w + 1 ) ) / 2 ) ] -= gamma_zxyt;
                                 }
                              }
                           }
                        }

                        // FCF_singlet[ xyz,tu ] -= delta_uw delta_xy Gamma_zt
                        if ( irrep_x == irrep_y ){
                           for ( int t = 0; t <= w; t++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zt = one_rdm[ d_z + z + LAS * ( d_ut + t ) ];
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + t + ( w * ( w + 1 ) ) / 2 ) ] -= gamma_zt;
                                 }
                              }
                           }
                        }
                     }
                     jump_row += num_x * num_y * num_z;
                  }
               }
               jump_col += ( num_ut * ( num_ut + 1 ) ) / 2;
            }
         }
      }

      for ( int irrep_right = 1; irrep_right < num_irreps; irrep_right++ ){

         assert( size_B_singlet[ irrep_right ] == size_F_singlet[ irrep_right ] ); // At construction
         const int SIZE_right = size_B_singlet[ irrep_right ];
         const int irrep_w = Irreps::directProd( irrep_left, irrep_right );
         const int num_w   = indices->getNDMRG( irrep_w );
         FAB_singlet[ irrep_left ][ irrep_right ] = new double*[ num_w ];
         FCF_singlet[ irrep_left ][ irrep_right ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FAB_singlet[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];
            FCF_singlet[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * ABptr = FAB_singlet[ irrep_left ][ irrep_right ][ w ];
            double * CFptr = FCF_singlet[ irrep_left ][ irrep_right ][ w ];

            int jump_col = 0;
            for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
               const int irrep_u = Irreps::directProd( irrep_right, irrep_t );
               if ( irrep_t < irrep_u ){
                  const int d_t   = indices->getDMRGcumulative( irrep_t );
                  const int num_t = indices->getNDMRG( irrep_t );
                  const int d_u   = indices->getDMRGcumulative( irrep_u );
                  const int num_u = indices->getNDMRG( irrep_u );
                  assert( jump_col == jump_BF_active( indices, irrep_t, irrep_u, +1 ) );

                  const int jump_AB1 = jump_AC_active( indices, irrep_u, irrep_t, irrep_w );
                  const int jump_AB2 = jump_AC_active( indices, irrep_t, irrep_u, irrep_w );
                  const int jump_CF1 = jump_AC_active( indices, irrep_u, irrep_w, irrep_t );
                  const int jump_CF2 = jump_AC_active( indices, irrep_t, irrep_w, irrep_u );

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int xyz = 0; xyz < SIZE_left; xyz++ ){
                           ABptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = ( - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB1 + u + num_u * ( t + num_t * w )) ]
                                                                                       - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB2 + t + num_t * ( u + num_u * w )) ] );
                           CFptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = ( + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF1 + u + num_u * ( w + num_w * t )) ]
                                                                                       + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF2 + t + num_t * ( w + num_w * u )) ] );
                        }
                     }
                  }

                  int jump_row = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int d_x   = indices->getDMRGcumulative( irrep_x );
                     const int num_x = indices->getNDMRG( irrep_x );
                     for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                        const int irrep_z = Irreps::directProd( Irreps::directProd( irrep_left, irrep_x ), irrep_y );
                        const int d_y   = indices->getDMRGcumulative( irrep_y );
                        const int num_y = indices->getNDMRG( irrep_y );
                        const int d_z   = indices->getDMRGcumulative( irrep_z );
                        const int num_z = indices->getNDMRG( irrep_z );
                        assert( jump_row == jump_AC_active( indices, irrep_x, irrep_y, irrep_z ) );

                        if ( irrep_t == irrep_w ){

                           // FAB_singlet[ xyz,tu ] -= delta_tw Gamma_zuyx
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zuyx = two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                       ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] -= gamma_zuyx;
                                    }
                                 }
                              }
                           }

                           // FAB_singlet[ xyz,tu ] -= delta_tw delta_uy Gamma_zx
                           if ( irrep_u == irrep_y ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                    for ( int uy = 0; uy < num_y; uy++ ){
                                       ABptr[ jump_row + x + num_x * ( uy + num_y * z ) + SIZE_left * ( jump_col + w + num_w * uy ) ] -= gamma_zx;
                                    }
                                 }
                              }
                           }

                           // FAB_singlet[ xyz,tu ] += 2 delta_tw delta_ux Gamma_zy
                           if ( irrep_u == irrep_x ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                    for ( int ux = 0; ux < num_x; ux++ ){
                                       ABptr[ jump_row + ux + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * ux ) ] += 2 * gamma_zy;
                                    }
                                 }
                              }
                           }

                           // FCF_singlet[ xyz,tu ] -= delta_tw Gamma_zxyu
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zxyu = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_u + u ))) ];
                                       CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] -= gamma_zxyu;
                                    }
                                 }
                              }
                           }

                           // FCF_singlet[ xyz,tu ] -= delta_tw delta_xy Gamma_zu
                           if ( irrep_x == irrep_y ){
                              for ( int u = 0; u < num_u; u++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zu = one_rdm[ d_z + z + LAS * ( d_u + u ) ];
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] -= gamma_zu;
                                    }
                                 }
                              }
                           }
                        }

                        if ( irrep_u == irrep_w ){

                           // FAB_singlet[ xyz,tu ] -= delta_uw Gamma_ztyx
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_ztyx = two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                       ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] -= gamma_ztyx;
                                    }
                                 }
                              }
                           }

                           // FAB_singlet[ xyz,tu ] -= delta_uw delta_ty Gamma_zx
                           if ( irrep_t == irrep_y ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                    for ( int ty = 0; ty < num_y; ty++ ){
                                       ABptr[ jump_row + x + num_x * ( ty + num_y * z ) + SIZE_left * ( jump_col + ty + num_y * w ) ] -= gamma_zx;
                                    }
                                 }
                              }
                           }

                           // FAB_singlet[ xyz,tu ] += 2 delta_uw delta_tx Gamma_zy
                           if ( irrep_t == irrep_x ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                    for ( int tx = 0; tx < num_x; tx++ ){
                                       ABptr[ jump_row + tx + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + tx + num_x * w ) ] += 2 * gamma_zy;
                                    }
                                 }
                              }
                           }

                           // FCF_singlet[ xyz,tu ] -= delta_uw Gamma_zxyt
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zxyt = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_t + t ))) ];
                                       CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] -= gamma_zxyt;
                                    }
                                 }
                              }
                           }

                           // FCF_singlet[ xyz,tu ] -= delta_uw delta_xy Gamma_zt
                           if ( irrep_x == irrep_y ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zt = one_rdm[ d_z + z + LAS * ( d_t + t ) ];
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] -= gamma_zt;
                                    }
                                 }
                              }
                           }
                        }
                        jump_row += num_x * num_y * num_z;
                     }
                  }
                  jump_col += num_t * num_u;
               }
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FAB_FCF_triplet(){

   /*
      FAB triplet : < E_zy E_lx | E_kw | TB_tiuj > = ( delta_ik delta_jl - delta_jk delta_il ) / sqrt( 1 + delta_ij ) * FAB_triplet[ Il ][ Ii x Ij ][ w ][ xyz, tu ]

      FCF triplet : < E_zy E_xd | E_wc | TF_atbu > = ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) * FCF_triplet[ Id ][ Ia x Ib ][ w ][ xyz, tu ]

            FAB_triplet[ Il ][ Ii x Ij ][ w ][ xyz, tu ] = ( + 6 delta_tw delta_ux Gamma_zy
                                                             - 6 delta_uw delta_tx Gamma_zy
                                                             - 3 delta_tw Gamma_zuyx
                                                             - 3 delta_tw delta_uy Gamma_zx
                                                             + 3 delta_uw Gamma_ztyx
                                                             + 3 delta_uw delta_ty Gamma_zx
                                                             - SAA[ Il ][ xyz, utw ]
                                                             + SAA[ Il ][ xyz, tuw ]
                                                           )

            FCF_triplet[ Id ][ Ia x Ib ][ w ][ xyz, tu ] = ( + SCC[ Id ][ xyz, uwt ]
                                                             - SCC[ Id ][ xyz, twu ]
                                                             - delta_uw Gamma_zxyt
                                                             - delta_uw delta_xy Gamma_zt
                                                             + delta_tw Gamma_zxyu
                                                             + delta_tw delta_xy Gamma_zu
                                                           )
   */

   FAB_triplet = new double***[ num_irreps ];
   FCF_triplet = new double***[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      assert( size_A[ irrep_left ] == size_C[ irrep_left ] ); // At construction
      const int SIZE_left = size_A[ irrep_left ];
      FAB_triplet[ irrep_left ] = new double**[ num_irreps ];
      FCF_triplet[ irrep_left ] = new double**[ num_irreps ];

      { // irrep_right == 0

         assert( size_B_triplet[ 0 ] == size_F_triplet[ 0 ] ); // At construction
         const int SIZE_right = size_B_triplet[ 0 ];
         const int num_w = indices->getNDMRG( irrep_left );
         FAB_triplet[ irrep_left ][ 0 ] = new double*[ num_w ];
         FCF_triplet[ irrep_left ][ 0 ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FAB_triplet[ irrep_left ][ 0 ][ w ] = new double[ SIZE_left * SIZE_right ];
            FCF_triplet[ irrep_left ][ 0 ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * ABptr = FAB_triplet[ irrep_left ][ 0 ][ w ];
            double * CFptr = FCF_triplet[ irrep_left ][ 0 ][ w ];

            int jump_col = 0;
            for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
               const int d_ut   = indices->getDMRGcumulative( irrep_ut );
               const int num_ut = indices->getNDMRG( irrep_ut );
               assert( jump_col == jump_BF_active( indices, irrep_ut, irrep_ut, -1 ) );

               const int jump_AB = jump_AC_active( indices, irrep_ut, irrep_ut, irrep_left );
               const int jump_CF = jump_AC_active( indices, irrep_ut, irrep_left, irrep_ut );

               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                     for ( int xyz = 0; xyz < SIZE_left; xyz++ ){
                        ABptr[ xyz + SIZE_left * ( jump_col + t + ( u * ( u - 1 ) ) / 2 ) ] = ( - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB + u + num_ut * ( t + num_ut * w )) ]
                                                                                                + SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB + t + num_ut * ( u + num_ut * w )) ] );
                        CFptr[ xyz + SIZE_left * ( jump_col + t + ( u * ( u - 1 ) ) / 2 ) ] = ( + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF + u + num_ut * ( w + num_w  * t )) ]
                                                                                                - SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF + t + num_ut * ( w + num_w  * u )) ] );
                     }
                  }
               }

               int jump_row = 0;
               for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                  const int d_x   = indices->getDMRGcumulative( irrep_x );
                  const int num_x = indices->getNDMRG( irrep_x );
                  for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                     const int irrep_z = Irreps::directProd( Irreps::directProd( irrep_left, irrep_x ), irrep_y );
                     const int d_y   = indices->getDMRGcumulative( irrep_y );
                     const int num_y = indices->getNDMRG( irrep_y );
                     const int d_z   = indices->getDMRGcumulative( irrep_z );
                     const int num_z = indices->getNDMRG( irrep_z );
                     assert( jump_row == jump_AC_active( indices, irrep_x, irrep_y, irrep_z ) );

                     if ( irrep_ut == irrep_left ){

                        // FAB_triplet[ xyz,tu ] -= 3 delta_tw Gamma_zuyx
                        for ( int u = w+1; u < num_ut; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zuyx = two_rdm[ d_z + z + LAS * ( d_ut + u + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                    ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * gamma_zuyx;
                                 }
                              }
                           }
                        }

                        // FAB_triplet[ xyz,tu ] -= 3 delta_tw delta_uy Gamma_zx
                        if ( irrep_ut == irrep_y ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                 for ( int uy = w+1; uy < num_y; uy++ ){
                                    ABptr[ jump_row + x + num_x * ( uy + num_y * z ) + SIZE_left * ( jump_col + w + ( uy * ( uy - 1 ) ) / 2 ) ] -= 3 * gamma_zx;
                                 }
                              }
                           }
                        }

                        // FAB_triplet[ xyz,tu ] += 6 delta_tw delta_ux Gamma_zy
                        if ( irrep_ut == irrep_x ){
                           for ( int y = 0; y < num_y; y++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                 for ( int ux = w+1; ux < num_x; ux++ ){
                                    ABptr[ jump_row + ux + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( ux * ( ux - 1 ) ) / 2 ) ] += 6 * gamma_zy;
                                 }
                              }
                           }
                        }

                        // FCF_triplet[ xyz,tu ] += delta_tw Gamma_zxyu
                        for ( int u = w+1; u < num_ut; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyu = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_ut + u ))) ];
                                    CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + ( u * ( u - 1 ) ) / 2 ) ] += gamma_zxyu;
                                 }
                              }
                           }
                        }

                        // FCF_triplet[ xyz,tu ] += delta_tw delta_xy Gamma_zu
                        if ( irrep_x == irrep_y ){
                           for ( int u = w+1; u < num_ut; u++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zu = one_rdm[ d_z + z + LAS * ( d_ut + u ) ];
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + w + ( u * ( u - 1 ) ) / 2 ) ] += gamma_zu;
                                 }
                              }
                           }
                        }

                        // FAB_triplet[ xyz,tu ] += 3 delta_uw Gamma_ztyx
                        for ( int t = 0; t < w; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_ztyx = two_rdm[ d_z + z + LAS * ( d_ut + t + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                    ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + ( w * ( w - 1 ) ) / 2 ) ] += 3 * gamma_ztyx;
                                 }
                              }
                           }
                        }

                        // FAB_triplet[ xyz,tu ] += 3 delta_uw delta_ty Gamma_zx
                        if ( irrep_ut == irrep_y ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                 for ( int ty = 0; ty < w; ty++ ){
                                    ABptr[ jump_row + x + num_x * ( ty + num_y * z ) + SIZE_left * ( jump_col + ty + ( w * ( w - 1 ) ) / 2 ) ] += 3 * gamma_zx;
                                 }
                              }
                           }
                        }

                        // FAB_triplet[ xyz,tu ] -= 6 delta_uw delta_tx Gamma_zy
                        if ( irrep_ut == irrep_x ){
                           for ( int y = 0; y < num_y; y++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                 for ( int tx = 0; tx < w; tx++ ){
                                    ABptr[ jump_row + tx + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + tx + ( w * ( w - 1 ) ) / 2 ) ] -= 6 * gamma_zy;
                                 }
                              }
                           }
                        }

                        // FCF_triplet[ xyz,tu ] -= delta_uw Gamma_zxyt
                        for ( int t = 0; t < w; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyt = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_ut + t ))) ];
                                    CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + ( w * ( w - 1 ) ) / 2 ) ] -= gamma_zxyt;
                                 }
                              }
                           }
                        }

                        // FCF_triplet[ xyz,tu ] -= delta_uw delta_xy Gamma_zt
                        if ( irrep_x == irrep_y ){
                           for ( int t = 0; t < w; t++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zt = one_rdm[ d_z + z + LAS * ( d_ut + t ) ];
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + t + ( w * ( w - 1 ) ) / 2 ) ] -= gamma_zt;
                                 }
                              }
                           }
                        }
                     }
                     jump_row += num_x * num_y * num_z;
                  }
               }
               jump_col += ( num_ut * ( num_ut - 1 ) ) / 2;
            }
         }
      }

      for ( int irrep_right = 1; irrep_right < num_irreps; irrep_right++ ){

         assert( size_B_triplet[ irrep_right ] == size_F_triplet[ irrep_right ] ); // At construction
         const int SIZE_right = size_B_triplet[ irrep_right ];
         const int irrep_w = Irreps::directProd( irrep_left, irrep_right );
         const int num_w   = indices->getNDMRG( irrep_w );
         FAB_triplet[ irrep_left ][ irrep_right ] = new double*[ num_w ];
         FCF_triplet[ irrep_left ][ irrep_right ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FAB_triplet[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];
            FCF_triplet[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * ABptr = FAB_triplet[ irrep_left ][ irrep_right ][ w ];
            double * CFptr = FCF_triplet[ irrep_left ][ irrep_right ][ w ];

            int jump_col = 0;
            for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
               const int irrep_u = Irreps::directProd( irrep_right, irrep_t );
               if ( irrep_t < irrep_u ){
                  const int d_t   = indices->getDMRGcumulative( irrep_t );
                  const int num_t = indices->getNDMRG( irrep_t );
                  const int d_u   = indices->getDMRGcumulative( irrep_u );
                  const int num_u = indices->getNDMRG( irrep_u );
                  assert( jump_col == jump_BF_active( indices, irrep_t, irrep_u, -1 ) );

                  const int jump_AB1 = jump_AC_active( indices, irrep_u, irrep_t, irrep_w ); 
                  const int jump_AB2 = jump_AC_active( indices, irrep_t, irrep_u, irrep_w );
                  const int jump_CF1 = jump_AC_active( indices, irrep_u, irrep_w, irrep_t ); 
                  const int jump_CF2 = jump_AC_active( indices, irrep_t, irrep_w, irrep_u );

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int xyz = 0; xyz < SIZE_left; xyz++ ){
                           ABptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = ( - SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB1 + u + num_u * ( t + num_t * w )) ]
                                                                                       + SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AB2 + t + num_t * ( u + num_u * w )) ] );
                           CFptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = ( + SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF1 + u + num_u * ( w + num_w * t )) ]
                                                                                       - SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CF2 + t + num_t * ( w + num_w * u )) ] );
                        }
                     }
                  }

                  int jump_row = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int d_x   = indices->getDMRGcumulative( irrep_x );
                     const int num_x = indices->getNDMRG( irrep_x );
                     for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                        const int irrep_z = Irreps::directProd( Irreps::directProd( irrep_left, irrep_x ), irrep_y );
                        const int d_y   = indices->getDMRGcumulative( irrep_y );
                        const int num_y = indices->getNDMRG( irrep_y );
                        const int d_z   = indices->getDMRGcumulative( irrep_z );
                        const int num_z = indices->getNDMRG( irrep_z );
                        assert( jump_row == jump_AC_active( indices, irrep_x, irrep_y, irrep_z ) );

                        if ( irrep_t == irrep_w ){

                           // FAB_triplet[ xyz,tu ] -= 3 delta_tw Gamma_zuyx
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zuyx = two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                       ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] -= 3 * gamma_zuyx;
                                    }
                                 }
                              }
                           }

                           // FAB_triplet[ xyz,tu ] -= 3 delta_tw delta_uy Gamma_zx
                           if ( irrep_u == irrep_y ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                    for ( int uy = 0; uy < num_y; uy++ ){
                                       ABptr[ jump_row + x + num_x * ( uy + num_y * z ) + SIZE_left * ( jump_col + w + num_w * uy ) ] -= 3 * gamma_zx;
                                    }
                                 }
                              }
                           }

                           // FAB_triplet[ xyz,tu ] += 6 delta_tw delta_ux Gamma_zy
                           if ( irrep_u == irrep_x ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                    for ( int ux = 0; ux < num_x; ux++ ){
                                       ABptr[ jump_row + ux + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * ux ) ] += 6 * gamma_zy;
                                    }
                                 }
                              }
                           }

                           // FCF_triplet[ xyz,tu ] += delta_tw Gamma_zxyu
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zxyu = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_u + u ))) ];
                                       CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] += gamma_zxyu;
                                    }
                                 }
                              }
                           }

                           // FCF_triplet[ xyz,tu ] += delta_tw delta_xy Gamma_zu
                           if ( irrep_x == irrep_y ){
                              for ( int u = 0; u < num_u; u++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zu = one_rdm[ d_z + z + LAS * ( d_u + u ) ];
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + w + num_w * u ) ] += gamma_zu;
                                    }
                                 }
                              }
                           }
                        }

                        if ( irrep_u == irrep_w ){

                           // FAB_triplet[ xyz,tu ] += 3 * delta_uw Gamma_ztyx
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_ztyx = two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_x + x ))) ];
                                       ABptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] += 3 * gamma_ztyx;
                                    }
                                 }
                              }
                           }

                           // FAB_triplet[ xyz,tu ] += 3 * delta_uw delta_ty Gamma_zx
                           if ( irrep_t == irrep_y ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zx = one_rdm[ d_z + z + LAS * ( d_x + x ) ];
                                    for ( int ty = 0; ty < num_y; ty++ ){
                                       ABptr[ jump_row + x + num_x * ( ty + num_y * z ) + SIZE_left * ( jump_col + ty + num_y * w ) ] += 3 * gamma_zx;
                                    }
                                 }
                              }
                           }

                           // FAB_triplet[ xyz,tu ] -= 6 delta_uw delta_tx Gamma_zy
                           if ( irrep_t == irrep_x ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zy = one_rdm[ d_z + z + LAS * ( d_y + y ) ];
                                    for ( int tx = 0; tx < num_x; tx++ ){
                                       ABptr[ jump_row + tx + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + tx + num_x * w ) ] -= 6 * gamma_zy;
                                    }
                                 }
                              }
                           }

                           // FCF_triplet[ xyz,tu ] -= delta_uw Gamma_zxyt
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       const double gamma_zxyt = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_t + t ))) ];
                                       CFptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] -= gamma_zxyt;
                                    }
                                 }
                              }
                           }

                           // FCF_triplet[ xyz,tu ] -= delta_uw delta_xy Gamma_zt
                           if ( irrep_x == irrep_y ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zt = one_rdm[ d_z + z + LAS * ( d_t + t ) ];
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       CFptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + t + num_t * w ) ] -= gamma_zt;
                                    }
                                 }
                              }
                           }
                        }
                        jump_row += num_x * num_y * num_z;
                     }
                  }
                  jump_col += num_t * num_u;
               }
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FAD_FCD(){

   /*
      FAD1: < E_zy E_jx E_wc E_ai E_tu > = delta_ac delta_ij FAD1[ Ij ][ Ii x Ia ][ w ][ xyz, tu ]
      FAD2: < E_zy E_jx E_wc E_ti E_au > = delta_ac delta_ij FAD2[ Ij ][ Ii x Ia ][ w ][ xyz, tu ]
      FCD1: < E_zy E_xb E_kw E_ai E_tu > = delta_ik delta_ab FCD1[ Ib ][ Ii x Ia ][ w ][ xyz, tu ]
      FCD2: < E_zy E_xb E_kw E_ti E_au > = delta_ik delta_ab FCD2[ Ib ][ Ii x Ia ][ w ][ xyz, tu ]

            FAD1[ Ij ][ Ii x Ia ][ w ][ xyz, tu ] = + SAA[ Ij ][ xyzwtu ]

            FAD2[ Ij ][ Ii x Ia ][ w ][ xyz, tu ] = + SAA[ Ij ][ xyztwu ]

            FCD1[ Ib ][ Ii x Ia ][ w ][ xyz, tu ] = - SCC[ Ib ][ xyzwtu ]

            FCD2[ Ib ][ Ii x Ia ][ w ][ xyz, tu ] = ( + 2 delta_tw Gamma_zxyu
                                                      + 2 delta_tw delta_xy Gamma_zu
                                                      + delta_tu Gamma_zxyw
                                                      + delta_tu delta_xy Gamma_zw
                                                      - SCC[ Ib ][ xyzutw ]
                                                    )
   */

   FAD = new double***[ num_irreps ];
   FCD = new double***[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_left = 0; irrep_left < num_irreps; irrep_left++ ){

      assert( size_A[ irrep_left ] == size_C[ irrep_left ] ); // At construction
      const int SIZE_left = size_A[ irrep_left ];
      FAD[ irrep_left ] = new double**[ num_irreps ];
      FCD[ irrep_left ] = new double**[ num_irreps ];

      for ( int irrep_right = 0; irrep_right < num_irreps; irrep_right++ ){

         const int SIZE_right = size_D[ irrep_right ];
         const int D2JUMP  = SIZE_right / 2;
         const int irrep_w = Irreps::directProd( irrep_left, irrep_right );
         const int d_w     = indices->getDMRGcumulative( irrep_w );
         const int num_w   = indices->getNDMRG( irrep_w );
         FAD[ irrep_left ][ irrep_right ] = new double*[ num_w ];
         FCD[ irrep_left ][ irrep_right ] = new double*[ num_w ];

         for ( int w = 0; w < num_w; w++ ){

            FAD[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];
            FCD[ irrep_left ][ irrep_right ][ w ] = new double[ SIZE_left * SIZE_right ];

            double * AD1ptr = FAD[ irrep_left ][ irrep_right ][ w ];
            double * AD2ptr = FAD[ irrep_left ][ irrep_right ][ w ] + SIZE_left * D2JUMP;
            double * CD1ptr = FCD[ irrep_left ][ irrep_right ][ w ];
            double * CD2ptr = FCD[ irrep_left ][ irrep_right ][ w ] + SIZE_left * D2JUMP;

            int jump_col = 0;
            for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
               const int irrep_u = Irreps::directProd( irrep_right, irrep_t );
               const int num_t = indices->getNDMRG( irrep_t );
               const int d_u   = indices->getDMRGcumulative( irrep_u );
               const int num_u = indices->getNDMRG( irrep_u );

               const int jump_AD1 = jump_AC_active( indices, irrep_w, irrep_t, irrep_u );
               const int jump_AD2 = jump_AC_active( indices, irrep_t, irrep_w, irrep_u );
               const int jump_CD2 = jump_AC_active( indices, irrep_u, irrep_t, irrep_w );

               for ( int t = 0; t < num_t; t++ ){
                  for ( int u = 0; u < num_u; u++ ){
                     for ( int xyz = 0; xyz < SIZE_left; xyz++ ){
                        AD1ptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = + SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AD1 + w + num_w * ( t + num_t * u )) ];
                        AD2ptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = + SAA[ irrep_left ][ xyz + SIZE_left * ( jump_AD2 + t + num_t * ( w + num_w * u )) ];
                        CD1ptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = - SCC[ irrep_left ][ xyz + SIZE_left * ( jump_AD1 + w + num_w * ( t + num_t * u )) ];
                        CD2ptr[ xyz + SIZE_left * ( jump_col + t + num_t * u ) ] = - SCC[ irrep_left ][ xyz + SIZE_left * ( jump_CD2 + u + num_u * ( t + num_t * w )) ];
                     }
                  }
               }

               int jump_row = 0;
               for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                  const int d_x   = indices->getDMRGcumulative( irrep_x );
                  const int num_x = indices->getNDMRG( irrep_x );
                  for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                     const int irrep_z = Irreps::directProd( Irreps::directProd( irrep_left, irrep_x ), irrep_y );
                     const int d_y   = indices->getDMRGcumulative( irrep_y );
                     const int num_y = indices->getNDMRG( irrep_y );
                     const int d_z   = indices->getDMRGcumulative( irrep_z );
                     const int num_z = indices->getNDMRG( irrep_z );
                     assert( jump_row == jump_AC_active( indices, irrep_x, irrep_y, irrep_z ) );

                     if ( irrep_t == irrep_w ){
                        // + 2 delta_tw Gamma_zxyu
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyu = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_u + u ))) ];
                                    CD2ptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + w + num_t * u ) ] += 2 * gamma_zxyu;
                                 }
                              }
                           }
                        }

                        // + 2 delta_tw delta_xy Gamma_zu
                        if ( irrep_x == irrep_y ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 const double gamma_zu = one_rdm[ d_z + z + LAS * ( d_u + u ) ];
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CD2ptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + w + num_t * u ) ] += 2 * gamma_zu;
                                 }
                              }
                           }
                        }
                     }

                     if ( irrep_t == irrep_u ){
                        // + delta_tu Gamma_zxyw
                        for ( int tu = 0; tu < num_u; tu++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    const double gamma_zxyw = two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_w + w ))) ];
                                    CD2ptr[ jump_row + x + num_x * ( y + num_y * z ) + SIZE_left * ( jump_col + tu + num_u * tu ) ] += gamma_zxyw;
                                 }
                              }
                           }
                        }

                        // + delta_tu delta_xy Gamma_zw
                        if ( irrep_x == irrep_y ){
                           for ( int z = 0; z < num_z; z++ ){
                              const double gamma_zw = one_rdm[ d_z + z + LAS * ( d_w + w ) ];
                              for ( int tu = 0; tu < num_u; tu++ ){
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    CD2ptr[ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE_left * ( jump_col + tu + num_u * tu ) ] += gamma_zw;
                                 }
                              }
                           }
                        }
                     }
                     jump_row += num_x * num_y * num_z;
                  }
               }
               jump_col += num_t * num_u;
            }
         }
      }
   }

}

void CheMPS2::CASPT2::make_FAA_FCC(){

   /*
      FAA: < E_zy E_jx ( f_pq E_pq ) E_ti E_uv > = delta_ij * ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )

            FAA[ Ii ][ xyztuv ] = ( - f_dot_4dm[ ztuyxv ]
                                    + 2 delta_tx f_dot_3dm[ zuyv ]
                                    - delta_uy f_dot_3dm[ ztvx ]
                                    - delta_ux f_dot_3dm[ ztyv ]
                                    - delta_ty f_dot_3dm[ zuxv ]
                                    + ( 2 delta_tx delta_uy - delta_ux delta_ty ) f_dot_2dm[ zv ]
                                    + sum_r SAA[ Ii ][ xyzruv ] f_rt
                                    + sum_r SAA[ Ii ][ xyztrv ] f_ru
                                    + sum_s SAA[ Ii ][ syztuv ] f_xs
                                    + sum_s SAA[ Ii ][ xsztuv ] f_ys
                                    - f_xt ( 2 * Gamma_{zuyv} + 2 * delta_yu Gamma_{zv} )
                                    - f_yu (   - Gamma_{ztvx} + 2 * delta_xt Gamma_{zv} )
                                    - f_yt (   - Gamma_{zuxv} -     delta_xu Gamma_{zv} )
                                    - f_xu (   - Gamma_{ztyv} -     delta_yt Gamma_{zv} )
                                  )

      FCC: < E_zy E_xb ( f_pq E_pq ) E_at E_uv > = delta_ab * ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )

            FCC[ Ia ][ xyztuv ] = ( + f_dot_4dm[ zxuytv ]
                                    + delta_uy f_dot_3dm[ xztv ]
                                    + delta_xy f_dot_3dm[ zutv ]
                                    + delta_ut f_dot_3dm[ zxyv ]
                                    + delta_ut delta_xy f_dot_2dm[ zv ]
                                    + sum_s f_ys SCC[ Ia ] [ xsztuv ]
                                    + sum_r f_ru SCC[ Ia ] [ xyztrv ]
                                    - f_yx ( Gamma_{zutv} + delta_ut Gamma_{zv} )
                                    - f_tu ( Gamma_{zxyv} + delta_yx Gamma_{zv} )
                                    - f_yu Gamma_{zxvt}
                                  )

   */

   FAA = new double*[ num_irreps ];
   FCC = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      assert( size_A[ irrep ] == size_C[ irrep ] ); // At construction
      const int SIZE = size_A[ irrep ];
      FAA[ irrep ] = new double[ SIZE * SIZE ];
      FCC[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int d_t    = indices->getDMRGcumulative( irrep_t );
         const int num_t  = indices->getNDMRG( irrep_t );
         const int nocc_t = indices->getNOCC( irrep_t );
         for ( int irrep_u = 0; irrep_u < num_irreps; irrep_u++ ){
            const int d_u     = indices->getDMRGcumulative( irrep_u );
            const int num_u   = indices->getNDMRG( irrep_u );
            const int nocc_u  = indices->getNOCC( irrep_u );
            const int irrep_v = Irreps::directProd( Irreps::directProd( irrep, irrep_t ), irrep_u );
            const int d_v     = indices->getDMRGcumulative( irrep_v );
            const int num_v   = indices->getNDMRG( irrep_v );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int d_x    = indices->getDMRGcumulative( irrep_x );
               const int num_x  = indices->getNDMRG( irrep_x );
               const int nocc_x = indices->getNOCC( irrep_x );
               for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                  const int d_y     = indices->getDMRGcumulative( irrep_y );
                  const int num_y   = indices->getNDMRG( irrep_y );
                  const int nocc_y  = indices->getNOCC( irrep_y );
                  const int irrep_z = Irreps::directProd( Irreps::directProd( irrep, irrep_x ), irrep_y );
                  const int d_z     = indices->getDMRGcumulative( irrep_z );
                  const int num_z   = indices->getNDMRG( irrep_z );

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){

                                    // FAA: - f_dot_4dm[ ztuyxv ]
                                    double val = - f_dot_4dm[ d_z + z + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_v + v ))))) ];

                                    // FAA: + sum_r SAA[ Ii ][ xyzruv ] f_rt
                                    for ( int r = 0; r < num_t; r++ ){
                                       val += ( fock->get( irrep_t, nocc_t + r, nocc_t + t )
                                              * SAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + r + num_t * ( u + num_u * v ) ) ] );
                                    }

                                    // FAA: + sum_r SAA[ Ii ][ xyztrv ] f_ru
                                    for ( int r = 0; r < num_u; r++ ){
                                       val += ( fock->get( irrep_u, nocc_u + r, nocc_u + u )
                                              * SAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( r + num_u * v ) ) ] );
                                    }

                                    // FAA: + sum_s SAA[ Ii ][ syztuv ] f_xs
                                    for ( int s = 0; s < num_x; s++ ){
                                       val += ( fock->get( irrep_x, nocc_x + x, nocc_x + s )
                                              * SAA[ irrep ][ jump_row + s + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ] );
                                    }

                                    // FAA: + sum_s SAA[ Ii ][ xsztuv ] f_ys
                                    for ( int s = 0; s < num_y; s++ ){
                                       val += ( fock->get( irrep_y, nocc_y + y, nocc_y + s )
                                              * SAA[ irrep ][ jump_row + x + num_x * ( s + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ] );
                                    }

                                    FAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ] = val;
                                 }
                              }
                           }
                        }
                     }
                  }

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){

                                    // FCC: + f_dot_4dm[ zxuytv ]
                                    double val = f_dot_4dm[ d_z + z + LAS * ( d_x + x + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_v + v ))))) ];

                                    // FCC: + sum_s f_ys SCC[ Ia ] [ xsztuv ]
                                    for ( int s = 0; s < num_y; s++ ){
                                       val += ( fock->get( irrep_y, nocc_y + y, nocc_y + s )
                                              * SCC[ irrep ][ jump_row + x + num_x * ( s + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ] );
                                    }

                                    // FCC: + sum_r f_ru SCC[ Ia ] [ xyztrv ]
                                    for ( int r = 0; r < num_u; r++ ){
                                       val += ( fock->get( irrep_u, nocc_u + r, nocc_u + u )
                                              * SCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( r + num_u * v ) ) ] );
                                    }

                                    FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ] = val;
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_t == irrep_x ){ // FAA: + 2 delta_tx f_dot_3dm[ zuyv ]
                     for ( int xt = 0; xt < num_t; xt++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + xt + num_x * ( y + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( u + num_u * v ) ) ]
                                       += 2 * f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_y ){ // FAA: - delta_uy f_dot_3dm[ tzxv ]
                     for ( int uy = 0; uy < num_u; uy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                       -= f_dot_3dm[ d_t + t + LAS * ( d_z + z + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_t == irrep_y ){ // FAA: - delta_ty f_dot_3dm[ zuxv ]
                     for ( int ty = 0; ty < num_t; ty++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + x + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( u + num_u * v ) ) ]
                                       -= f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_x ){ // FAA: - delta_ux f_dot_3dm[ ztyv ]
                     for ( int ux = 0; ux < num_u; ux++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + ux + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( ux + num_u * v ) ) ]
                                       -= f_dot_3dm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_y ){ // FCC: + delta_uy f_dot_3dm[ xztv ]
                     for ( int uy = 0; uy < num_y; uy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FCC[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_y * v ) ) ]
                                       += f_dot_3dm[ d_x + x + LAS * ( d_z + z + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_x == irrep_y ){ // FCC: + delta_xy f_dot_3dm[ zutv ]
                     for ( int xy = 0; xy < num_x; xy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                       += f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_t ){ // FCC: + delta_ut f_dot_3dm[ zxyv ]
                     for ( int ut = 0; ut < num_u; ut++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + ut + num_u * ( ut + num_u * v ) ) ]
                                       += f_dot_3dm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_x == irrep_t ){ // FAA: - 2 * f_xt * Gamma_{zuyv}
                     for ( int x = 0; x < num_x; x++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           const double f_xt = fock->get( irrep_t, nocc_t + x, nocc_t + t );
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          -= 2 * f_xt * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_y == irrep_u ){ // FAA: + f_yu * Gamma_{ztvx}
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_yu = fock->get( irrep_u, nocc_u + y, nocc_u + u );
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          += f_yu * two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_v + v + LAS * ( d_x + x ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_y == irrep_t ){ // FAA: + f_yt * Gamma_{zuxv}
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           const double f_yt = fock->get( irrep_t, nocc_t + y, nocc_t + t );
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          += f_yt * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_x == irrep_u ){ // FAA: + f_xu * Gamma_{ztyv}
                     for ( int x = 0; x < num_x; x++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_xu = fock->get( irrep_u, nocc_u + x, nocc_u + u );
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          += f_xu * two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_x == irrep_y ){ // FCC: - f_yx Gamma_{zutv}
                    for ( int x = 0; x < num_x; x++ ){
                        for ( int y = 0; y < num_x; y++ ){
                           const double f_yx = fock->get( irrep_x, nocc_x + y, nocc_x + x );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int u = 0; u < num_u; u++ ){
                                 for ( int t = 0; t < num_t; t++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_x * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          -= f_yx * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_t == irrep_u ){ // FCC: - f_tu Gamma_{zxyv}
                     for ( int t = 0; t < num_t; t++ ){
                        for ( int u = 0; u < num_t; u++ ){
                           const double f_tu = fock->get( irrep_t, nocc_t + t, nocc_t + u );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_t * v ) ) ]
                                          -= f_tu * two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_y == irrep_u ){ // FCC: - f_yu Gamma_{zxvt}
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int u = 0; u < num_y; u++ ){
                           const double f_yu = fock->get( irrep_y, nocc_y + y, nocc_y + u );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_y * v ) ) ]
                                          -= f_yu * two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_v + v + LAS * ( d_t + t ) ) ) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_t == irrep_x ) && ( irrep_u == irrep_y ) && ( irrep_z == irrep_v )){

                     // FAA: + 2 delta_tx delta_uy f_dot_2dm[ zv ]
                     for ( int xt = 0; xt < num_t; xt++ ){
                        for ( int uy = 0; uy < num_u; uy++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 FAA[ irrep ][ jump_row + xt + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( uy + num_u * v ) ) ]
                                    += 2 * f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }

                     // FAA: - 2 * f_xt * delta_yu Gamma_{zv}
                     for ( int x = 0; x < num_x; x++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           const double f_xt = fock->get( irrep_t, nocc_t + x, nocc_t + t );
                           for ( int uy = 0; uy < num_y; uy++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_y * v ) ) ]
                                       -= 2 * f_xt * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }

                     // FAA: - 2 * f_yu * delta_xt Gamma_{zv}
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_yu = fock->get( irrep_u, nocc_u + y, nocc_u + u );
                           for ( int xt = 0; xt < num_x; xt++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + xt + num_x * ( y + num_y * z ) + SIZE * ( jump_col + xt + num_x * ( u + num_u * v ) ) ]
                                       -= 2 * f_yu * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_x ) && ( irrep_t == irrep_y ) && ( irrep_z == irrep_v )){

                     // FAA: - delta_ux delta_ty f_dot_2dm[ zv ]
                     for ( int ty = 0; ty < num_t; ty++ ){
                        for ( int ux = 0; ux < num_u; ux++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 FAA[ irrep ][ jump_row + ux + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( ux + num_u * v ) ) ]
                                    -= f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }

                     // FAA: + f_xu * delta_yt * Gamma_{zv}
                     for ( int x = 0; x < num_x; x++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_xu = fock->get( irrep_u, nocc_u + x, nocc_u + u );
                           for ( int yt = 0; yt < num_y; yt++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + x + num_x * ( yt + num_y * z ) + SIZE * ( jump_col + yt + num_y * ( u + num_u * v ) ) ]
                                       += f_xu * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }

                     // FAA: + f_yt * delta_xu * Gamma_{zv}
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           const double f_yt = fock->get( irrep_t, nocc_t + y, nocc_t + t );
                           for ( int xu = 0; xu < num_x; xu++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FAA[ irrep ][ jump_row + xu + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( xu + num_x * v ) ) ]
                                       += f_yt * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_t ) && ( irrep_x == irrep_y ) && ( irrep_z == irrep_v )){

                     // FCC: + delta_ut delta_xy f_dot_2dm[ zv ]
                     for ( int xy = 0; xy < num_x; xy++ ){
                        for ( int tu = 0; tu < num_t; tu++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + tu + num_t * ( tu + num_t * v ) ) ]
                                    += f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }

                     // FCC: - f_tu delta_yx Gamma_{zv}
                     for ( int t = 0; t < num_t; t++ ){
                        for ( int u = 0; u < num_t; u++ ){
                           const double f_tu = fock->get( irrep_t, nocc_t + t, nocc_t + u );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + t + num_t * ( u + num_t * v ) ) ]
                                       -= f_tu * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }

                     // FCC: - f_yx delta_ut Gamma_{zv}
                     for ( int x = 0; x < num_x; x++ ){
                        for ( int y = 0; y < num_x; y++ ){
                           const double f_yx = fock->get( irrep_x, nocc_x + y, nocc_x + x );
                           for ( int ut = 0; ut < num_u; ut++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FCC[ irrep ][ jump_row + x + num_x * ( y + num_x * z ) + SIZE * ( jump_col + ut + num_u * ( ut + num_u * v ) ) ]
                                       -= f_yx * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y * num_z;
               }
            }
            jump_col += num_t * num_u * num_v;
         }
      }
   }

}

void CheMPS2::CASPT2::make_SAA_SCC(){

   /*
      SAA: < E_zy E_jx E_ti E_uv > = delta_ij SAA[ Ii ][ xyztuv ]

            SAA[ Ii ][ xyztuv ] = ( + 2 delta_tx Gamma_{zuyv}
                                    + 2 delta_tx delta_uy Gamma_{zv}
                                    - Gamma_{ztuyxv}
                                    - delta_uy Gamma_{tzxv}
                                    - delta_ty Gamma_{zuxv}
                                    - delta_ux Gamma_{ztyv}
                                    - delta_ux delta_ty Gamma_{zv}
                                  )

      SCC: < E_zy E_xb E_at E_uv > = delta_ab SCC[ Ia ][ xyztuv ]

            SCC[ Ia ][ xyztuv ] = ( + Gamma_{zxuytv}
                                    + delta_uy Gamma_{xztv}
                                    + delta_xy Gamma_{zutv}
                                    + delta_ut Gamma_{zxyv}
                                    + delta_ut delta_xy Gamma_{zv}
                                  )
   */

   SAA = new double*[ num_irreps ];
   SCC = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      assert( size_A[ irrep ] == size_C[ irrep ] ); // At construction
      const int SIZE = size_A[ irrep ];
      SAA[ irrep ] = new double[ SIZE * SIZE ];
      SCC[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int d_t   = indices->getDMRGcumulative( irrep_t );
         const int num_t = indices->getNDMRG( irrep_t );
         for ( int irrep_u = 0; irrep_u < num_irreps; irrep_u++ ){
            const int d_u     = indices->getDMRGcumulative( irrep_u );
            const int num_u   = indices->getNDMRG( irrep_u );
            const int irrep_v = Irreps::directProd( Irreps::directProd( irrep, irrep_t ), irrep_u );
            const int d_v     = indices->getDMRGcumulative( irrep_v );
            const int num_v   = indices->getNDMRG( irrep_v );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int d_x   = indices->getDMRGcumulative( irrep_x );
               const int num_x = indices->getNDMRG( irrep_x );
               for ( int irrep_y = 0; irrep_y < num_irreps; irrep_y++ ){
                  const int d_y     = indices->getDMRGcumulative( irrep_y );
                  const int num_y   = indices->getNDMRG( irrep_y );
                  const int irrep_z = Irreps::directProd( Irreps::directProd( irrep, irrep_x ), irrep_y );
                  const int d_z     = indices->getDMRGcumulative( irrep_z );
                  const int num_z   = indices->getNDMRG( irrep_z );

                  for ( int t = 0; t < num_t; t++ ){ // SAA: - Gamma_{ztuyxv}
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SAA[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                       = - three_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_v + v ))))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  for ( int t = 0; t < num_t; t++ ){ // SCC: + Gamma_{zxuytv}
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                       = three_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_v + v ))))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_t == irrep_x ){ // SAA: + 2 delta_tx Gamma_{zuyv}
                     for ( int xt = 0; xt < num_t; xt++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SAA[ irrep ][ jump_row + xt + num_x * ( y + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( u + num_u * v ) ) ]
                                       += 2 * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_y ){ // SAA: - delta_uy Gamma_{tzxv}
                     for ( int uy = 0; uy < num_u; uy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SAA[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                       -= two_rdm[ d_t + t + LAS * ( d_z + z + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_t == irrep_y ){ // SAA: - delta_ty Gamma_{zuxv}
                     for ( int ty = 0; ty < num_t; ty++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SAA[ irrep ][ jump_row + x + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( u + num_u * v ) ) ]
                                       -= two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_x ){ // SAA: - delta_ux Gamma_{ztyv}
                     for ( int ux = 0; ux < num_u; ux++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SAA[ irrep ][ jump_row + ux + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( ux + num_u * v ) ) ]
                                       -= two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_y ){ // SCC: + delta_uy Gamma_{xztv}
                     for ( int uy = 0; uy < num_u; uy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SCC[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                       += two_rdm[ d_x + x + LAS * ( d_z + z + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_x == irrep_y ){ // SCC: + delta_xy Gamma_{zutv}
                     for ( int xy = 0; xy < num_x; xy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int v = 0; v < num_v; v++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SCC[ irrep ][ jump_row + xy + num_x * ( xy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                       += two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if ( irrep_u == irrep_t ){ // SCC: + delta_ut Gamma_{zxyv}
                     for ( int ut = 0; ut < num_u; ut++ ){
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    SCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + ut + num_t * ( ut + num_u * v ) ) ]
                                       += two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                 }
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_t == irrep_x ) && ( irrep_u == irrep_y ) && ( irrep_z == irrep_v )){ // SAA: + 2 delta_tx delta_uy Gamma_{zv}
                     for ( int xt = 0; xt < num_t; xt++ ){
                        for ( int uy = 0; uy < num_u; uy++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 SAA[ irrep ][ jump_row + xt + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( uy + num_u * v ) ) ]
                                    += 2 * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_x ) && ( irrep_t == irrep_y ) && ( irrep_z == irrep_v )){ // SAA: - delta_ux delta_ty Gamma_{zv}
                     for ( int ty = 0; ty < num_t; ty++ ){
                        for ( int ux = 0; ux < num_u; ux++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 SAA[ irrep ][ jump_row + ux + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( ux + num_u * v ) ) ]
                                    -= one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_t ) && ( irrep_x == irrep_y ) && ( irrep_z == irrep_v )){ // SCC: + delta_ut delta_xy Gamma_{zv}
                     for ( int xy = 0; xy < num_x; xy++ ){
                        for ( int tu = 0; tu < num_t; tu++ ){
                           for ( int v = 0; v < num_z; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 SCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + tu + num_t * ( tu + num_t * v ) ) ]
                                    += one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y * num_z;
               }
            }
            jump_col += num_t * num_u * num_v;
         }
      }
   }

}

void CheMPS2::CASPT2::make_FDD(){

   /*
      FDD: < E_yx E_jb ( f_pq E_pq) E_ai E_tu > = delta_ab delta_ij ( FDD[ It x Iu ][ xytu ] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ It x Iu ][ xytu ] )

            FD1D1[ It x Iu ][ xytu ] = ( + sum_w f_xw SD1D1[ It x Iu ][ wytu ]
                                         + sum_v f_vt SD1D1[ Ix x Iy ][ xyvu ]
                                         + ( 2 * f_dot_3dm[ ytxu ] + 2 * delta_tx f_dot_2dm[ yu ] )
                                         - 2 f_xt Gamma_{yu}
                                       )

            FD2D2[ It x Iu ][ xytu ] = ( + sum_w f_xw SD2D2[ It x Iu ][ wytu ]
                                         + sum_v f_vt SD2D2[ Ix x Iy ][ xyvu ]
                                         + (   - f_dot_3dm[ ytux ] + 2 * delta_tx f_dot_2dm[ yu ] )
                                         - 2 f_xt Gamma_{yu}
                                       )

            FD1D2[ It x Iu ][ xytu ] = ( + sum_w f_xw SD1D2[ It x Iu ][ wytu ]
                                         + sum_v f_vt SD1D2[ Ix x Iy ][ xyvu ]
                                         + (   - f_dot_3dm[ ytxu ]     - delta_tx f_dot_2dm[ yu ] )
                                         + f_xt Gamma_{yu}
                                       )

            FD2D1[ It x Iu ][ xytu ] = ( + sum_w f_xw SD2D1[ It x Iu ][ wytu ]
                                         + sum_v f_vt SD2D1[ Ix x Iy ][ xyvu ]
                                         + (   - f_dot_3dm[ ytxu ]     - delta_tx f_dot_2dm[ yu ] )
                                         + f_xt Gamma_{yu}
                                       )
   */

   FDD = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int SIZE   = size_D[ irrep ];
      const int D2JUMP = SIZE / 2;
      FDD[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int d_t     = indices->getDMRGcumulative( irrep_t );
         const int num_t   = indices->getNDMRG( irrep_t );
         const int nocc_t  = indices->getNOCC( irrep_t );
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         const int d_u     = indices->getDMRGcumulative( irrep_u );
         const int num_u   = indices->getNDMRG( irrep_u );
         int jump_row = 0;
         for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
            const int d_x     = indices->getDMRGcumulative( irrep_x );
            const int num_x   = indices->getNDMRG( irrep_x );
            const int nocc_x  = indices->getNOCC( irrep_x );
            const int irrep_y = Irreps::directProd( irrep, irrep_x );
            const int d_y     = indices->getDMRGcumulative( irrep_y );
            const int num_y   = indices->getNDMRG( irrep_y );

            for ( int t = 0; t < num_t; t++ ){
               for ( int u = 0; u < num_u; u++ ){
                  for ( int x = 0; x < num_x; x++ ){
                     for ( int y = 0; y < num_y; y++ ){

                        const double f_3dm_ytxu = f_dot_3dm[ d_y + y + LAS * ( d_t + t + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                        const double f_3dm_ytux = f_dot_3dm[ d_y + y + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_x + x ))) ];
                        const int ptr = jump_row + x + num_x * y + SIZE * ( jump_col + t + num_t * u );

                        double valD1D1 = 2 * f_3dm_ytxu;
                        double valD1D2 =   - f_3dm_ytxu;
                        double valD2D1 =   - f_3dm_ytxu;
                        double valD2D2 =   - f_3dm_ytux;

                        for ( int v = 0; v < num_t; v++ ){
                           const double f_vt = fock->get( irrep_t, nocc_t + v, nocc_t + t );
                           valD1D1 += f_vt * SDD[ irrep ][ ptr + SIZE * ( v - t )                          ];
                           valD1D2 += f_vt * SDD[ irrep ][ ptr + SIZE * ( v - t ) +          SIZE * D2JUMP ];
                           valD2D1 += f_vt * SDD[ irrep ][ ptr + SIZE * ( v - t ) + D2JUMP                 ];
                           valD2D2 += f_vt * SDD[ irrep ][ ptr + SIZE * ( v - t ) + D2JUMP + SIZE * D2JUMP ];
                        }

                        for ( int w = 0; w < num_x; w++ ){
                           const double f_xw = fock->get( irrep_x, nocc_x + x, nocc_x + w );
                           valD1D1 += f_xw * SDD[ irrep ][ ptr + w - x                          ];
                           valD1D2 += f_xw * SDD[ irrep ][ ptr + w - x +          SIZE * D2JUMP ];
                           valD2D1 += f_xw * SDD[ irrep ][ ptr + w - x + D2JUMP                 ];
                           valD2D2 += f_xw * SDD[ irrep ][ ptr + w - x + D2JUMP + SIZE * D2JUMP ];
                        }

                        FDD[ irrep ][ ptr                          ] = valD1D1;
                        FDD[ irrep ][ ptr +          SIZE * D2JUMP ] = valD1D2;
                        FDD[ irrep ][ ptr + D2JUMP                 ] = valD2D1;
                        FDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ] = valD2D2;
                     }
                  }
               }
            }

            if (( irrep_x == irrep_t ) && ( irrep_y == irrep_u )){
               for ( int xt = 0; xt < num_x; xt++ ){
                  for ( int u = 0; u < num_u; u++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        const double f_2dm_yu = f_dot_2dm[ d_y + y + LAS * ( d_u + u ) ];
                        FDD[ irrep ][          jump_row + xt + num_x * y + SIZE * (          jump_col + xt + num_x * u ) ] += 2 * f_2dm_yu; // FD1D1
                        FDD[ irrep ][          jump_row + xt + num_x * y + SIZE * ( D2JUMP + jump_col + xt + num_x * u ) ] -=     f_2dm_yu; // FD1D2
                        FDD[ irrep ][ D2JUMP + jump_row + xt + num_x * y + SIZE * (          jump_col + xt + num_x * u ) ] -=     f_2dm_yu; // FD2D1
                        FDD[ irrep ][ D2JUMP + jump_row + xt + num_x * y + SIZE * ( D2JUMP + jump_col + xt + num_x * u ) ] += 2 * f_2dm_yu; // FD2D2
                     }
                  }
               }
               for ( int t = 0; t < num_t; t++ ){
                  for ( int u = 0; u < num_u; u++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           const double f_xt_gamma_yu = fock->get( irrep_t, nocc_t + x, nocc_t + t ) * one_rdm[ d_y + y + LAS * ( d_u + u ) ];
                           FDD[ irrep ][          jump_row + x + num_x * y + SIZE * (          jump_col + t + num_x * u ) ] -= 2 * f_xt_gamma_yu; // FD1D1
                           FDD[ irrep ][          jump_row + x + num_x * y + SIZE * ( D2JUMP + jump_col + t + num_x * u ) ] +=     f_xt_gamma_yu; // FD1D2
                           FDD[ irrep ][ D2JUMP + jump_row + x + num_x * y + SIZE * (          jump_col + t + num_x * u ) ] +=     f_xt_gamma_yu; // FD2D1
                           FDD[ irrep ][ D2JUMP + jump_row + x + num_x * y + SIZE * ( D2JUMP + jump_col + t + num_x * u ) ] -= 2 * f_xt_gamma_yu; // FD2D2
                        }
                     }
                  }
               }
            }
            jump_row += num_x * num_y;
         }
         jump_col += num_t * num_u;
      }
   }

}

void CheMPS2::CASPT2::make_SDD(){

   /*
      SD1D1: < E_yx E_jb E_ai E_tu > = delta_ab delta_ij SD1D1[ It x Iu ][ xytu ]

            SD1D1[ It x Iu ][ xytu ] = ( + 2 * Gamma_{ytxu}
                                         + 2 * delta_tx Gamma_{yu}
                                       )

      SD2D2: < E_yb E_jx E_ti E_au > = delta_ab delta_ij SD2D2[ xytu ]

            SD2D2[ It x Iu ][ xytu ] = ( + 2 * delta_tx Gamma_{yu}
                                         - Gamma_{ytux}
                                       )

      SD1D2: < E_yx E_jb E_ti E_au > = delta_ab delta_ij SD1D2[ xytu ]

            SD1D2[ It x Iu ][ xytu ] = ( - Gamma_{ytxu}
                                         - delta_tx Gamma_{yu}
                                       )

      SD2D1: < E_yb E_jx E_ai E_tu > = delta_ab delta_ij SD2D1[ xytu ]

            SD2D1[ It x Iu ][ xytu ] = ( - Gamma_{ytxu}
                                         - delta_tx Gamma_{yu}
                                       )
   */

   SDD = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_ai = 0; irrep_ai < num_irreps; irrep_ai++ ){

      const int SIZE   = size_D[ irrep_ai ];
      const int D2JUMP = SIZE / 2;
      SDD[ irrep_ai ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int d_t     = indices->getDMRGcumulative( irrep_t );
         const int num_t   = indices->getNDMRG( irrep_t );
         const int irrep_u = Irreps::directProd( irrep_ai, irrep_t );
         const int d_u     = indices->getDMRGcumulative( irrep_u );
         const int num_u   = indices->getNDMRG( irrep_u );
         int jump_row = 0;
         for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
            const int d_x     = indices->getDMRGcumulative( irrep_x );
            const int num_x   = indices->getNDMRG( irrep_x );
            const int irrep_y = Irreps::directProd( irrep_ai, irrep_x );
            const int d_y     = indices->getDMRGcumulative( irrep_y );
            const int num_y   = indices->getNDMRG( irrep_y );

            for ( int t = 0; t < num_t; t++ ){
               for ( int u = 0; u < num_u; u++ ){
                  for ( int x = 0; x < num_x; x++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        const double gamma_ytxu = two_rdm[ d_y + y + LAS * ( d_t + t + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                        const double gamma_ytux = two_rdm[ d_y + y + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_x + x ))) ];
                        SDD[ irrep_ai ][          jump_row + x + num_x * y + SIZE * (          jump_col + t + num_t * u ) ] = 2 * gamma_ytxu; // SD1D1
                        SDD[ irrep_ai ][          jump_row + x + num_x * y + SIZE * ( D2JUMP + jump_col + t + num_t * u ) ] =   - gamma_ytxu; // SD1D2
                        SDD[ irrep_ai ][ D2JUMP + jump_row + x + num_x * y + SIZE * (          jump_col + t + num_t * u ) ] =   - gamma_ytxu; // SD2D1
                        SDD[ irrep_ai ][ D2JUMP + jump_row + x + num_x * y + SIZE * ( D2JUMP + jump_col + t + num_t * u ) ] =   - gamma_ytux; // SD2D2
                     }
                  }
               }
            }

            if (( irrep_x == irrep_t ) && ( irrep_y == irrep_u )){
               for ( int xt = 0; xt < num_x; xt++ ){
                  for ( int u = 0; u < num_u; u++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        const double gamma_yu = one_rdm[ d_y + y + LAS * ( d_u + u ) ];
                        SDD[ irrep_ai ][          jump_row + xt + num_x * y + SIZE * (          jump_col + xt + num_x * u ) ] += 2 * gamma_yu; // SD1D1
                        SDD[ irrep_ai ][          jump_row + xt + num_x * y + SIZE * ( D2JUMP + jump_col + xt + num_x * u ) ] -=     gamma_yu; // SD1D2
                        SDD[ irrep_ai ][ D2JUMP + jump_row + xt + num_x * y + SIZE * (          jump_col + xt + num_x * u ) ] -=     gamma_yu; // SD2D1
                        SDD[ irrep_ai ][ D2JUMP + jump_row + xt + num_x * y + SIZE * ( D2JUMP + jump_col + xt + num_x * u ) ] += 2 * gamma_yu; // SD2D2
                     }
                  }
               }
            }
            jump_row += num_x * num_y;
         }
         jump_col += num_t * num_u;
      }
   }

}

void CheMPS2::CASPT2::make_FBB_FFF_singlet(){

   /*
      FBB singlet: < SB_xkyl | ( f_pq E_pq ) | SB_tiuj >
                      = 2 delta_ik delta_jl ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )

      FFF singlet: < SF_cxdy | ( f_pq E_pq ) | SF_atbu >
                      = 2 delta_ac delta_bd ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )

         FBB_singlet[ Iij ][ xytu ] = ( + f_dot_3dm[ utyx ]
                                        + f_dot_3dm[ tuyx ]
                                        + 2 delta_uy delta_tx f_dot_1dm
                                        + 2 delta_ux delta_ty f_dot_1dm
                                        -   delta_uy f_dot_2dm[ tx ]
                                        -   delta_tx f_dot_2dm[ uy ]
                                        -   delta_ux f_dot_2dm[ ty ]
                                        -   delta_ty f_dot_2dm[ ux ]
                                        + sum_v f_vu SBB_singlet[ xytv ]
                                        + sum_v f_vt SBB_singlet[ xyvu ]
                                        + sum_w f_yw SBB_singlet[ xwtu ]
                                        + sum_w f_xw SBB_singlet[ wytu ]
                                        - f_xt ( 2 delta_uy - Gamma_{uy} )
                                        - f_yu ( 2 delta_tx - Gamma_{tx} )
                                        - f_xu ( 2 delta_ty - Gamma_{ty} )
                                        - f_yt ( 2 delta_ux - Gamma_{ux} )
                                      )

         FFF_singlet[ Iab ][ xytu ] = ( + f_dot_3dm[ yxut ]
                                        + f_dot_3dm[ yxtu ]
                                      )

   */

   FBB_singlet = new double*[ num_irreps ];
   FFF_singlet = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_singlet[ 0 ] == size_F_singlet[ 0 ] ); // At construction
      const int SIZE = size_B_singlet[ 0 ];
      FBB_singlet[ 0 ] = new double[ SIZE * SIZE ];
      FFF_singlet[ 0 ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
         const int d_ut    = indices->getDMRGcumulative( irrep_ut );
         const int num_ut  = indices->getNDMRG( irrep_ut );
         const int nocc_ut = indices->getNOCC( irrep_ut );
         int jump_row = 0;
         for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
            const int d_xy    = indices->getDMRGcumulative( irrep_xy );
            const int num_xy  = indices->getNDMRG( irrep_xy );
            const int nocc_xy = indices->getNOCC( irrep_xy );
            const int shift   = jump_row + SIZE * jump_col;

            for ( int t = 0; t < num_ut; t++ ){
               for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy

                        // + f_dot_3dm_{utyx} + f_dot_3dm_{utxy}
                        const double fdotsum = ( f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                               + f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );

                        double sbb_value = fdotsum;

                        // FBB: + sum_v f_vu SBB_singlet[ xytv ] = + sum_(v>=t) f_vu SBB_singlet[ xytv ] + sum_(v<t) f_vu SBB_singlet[ xyvt ]
                        for ( int v = 0; v < t; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( v + ( t * ( t + 1 ) ) / 2 ) ] );
                        }
                        for ( int v = t; v < num_ut; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( v * ( v + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_v f_vt SBB_singlet[ xyvu ] = + sum_(v<u) f_vt SBB_singlet[ xyvu ] + sum_(v>=u) f_vt SBB_singlet[ xyuv ]
                        for ( int v = 0; v < u; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( v + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int v = u; v < num_ut; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( u + ( v * ( v + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_yw SBB_singlet[ xwtu ] = + sum_(w>=x) f_yw SBB_singlet[ xwtu ] + sum_(w<x) f_yw SBB_singlet[ wxtu ]
                        for ( int w = 0; w < x; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_singlet[ 0 ][ shift + w + ( x * ( x + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int w = x; w < num_xy; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_singlet[ 0 ][ shift + x + ( w * ( w + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_xw SBB_singlet[ wytu ] = + sum_(w<y) f_xw SBB_singlet[ wytu ] + sum_(w>=y) f_xw SBB_singlet[ ywtu ]
                        for ( int w = 0; w < y; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_singlet[ 0 ][ shift + w + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int w = y; w < num_xy; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_singlet[ 0 ][ shift + y + ( w * ( w + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }

                        FBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] = sbb_value;
                        FFF_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] = fdotsum;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy and d_ut = d_xy

               // FBB: + 2 ( delta_uy delta_tx + delta_ux delta_ty ) f_dot_1dm --> last term only if x <= y == t <= u or hence all equal
               for ( int t = 0; t < num_ut; t++ ){
                  FBB_singlet[ 0 ][ shift + t + ( t * ( t + 1 ) ) / 2 + SIZE * ( t + ( t * ( t + 1 ) ) / 2 ) ] += 4 * f_dot_1dm;
                  for ( int u = t+1; u < num_ut; u++ ){
                     FBB_singlet[ 0 ][ shift + t + ( u * ( u + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] += 2 * f_dot_1dm;
                  }
               }

               // FBB: - delta_uy f_dot_2dm[ tx ]
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                     for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                        const double val_tx = f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ];
                        FBB_singlet[ 0 ][ shift + x + ( uy * ( uy + 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy + 1 ) ) / 2 ) ] -= val_tx;
                     }
                  }
               }

               // FBB: - delta_tx f_dot_2dm[ uy ]
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                     for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                        const double val_uy = f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ];
                        FBB_singlet[ 0 ][ shift + tx + ( y * ( y + 1 ) ) / 2 + SIZE * ( tx + ( u * ( u + 1 ) ) / 2 ) ] -= val_uy;
                     }
                  }
               }

               // FBB: - delta_ux f_dot_2dm[ ty ]
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                     for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                        const double val_ty = f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ];
                        FBB_singlet[ 0 ][ shift + ux + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux + 1 ) ) / 2 ) ] -= val_ty;
                     }
                  }
               }

               // FBB: - delta_ty f_dot_2dm[ ux ]
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                     for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                        const double val_ux = f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ];
                        FBB_singlet[ 0 ][ shift + x + ( ty * ( ty + 1 ) ) / 2 + SIZE * ( ty + ( u * ( u + 1 ) ) / 2 ) ] -= val_ux;
                     }
                  }
               }

               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                     for ( int x = 0; x < num_ut; x++ ){
                        for ( int y = x; y < num_ut; y++ ){ // 0 <= x <= y < num_xy == num_ut
                           const double value = ( fock->get( irrep_ut, nocc_ut + x, nocc_ut + t ) * SEE[ irrep_ut ][ u + num_ut * y ]    // - f_xt ( 2 delta_uy - Gamma_{uy} )
                                                + fock->get( irrep_ut, nocc_ut + y, nocc_ut + u ) * SEE[ irrep_ut ][ t + num_ut * x ]    // - f_yu ( 2 delta_tx - Gamma_{tx} )
                                                + fock->get( irrep_ut, nocc_ut + x, nocc_ut + u ) * SEE[ irrep_ut ][ t + num_ut * y ]    // - f_xu ( 2 delta_ty - Gamma_{ty} )
                                                + fock->get( irrep_ut, nocc_ut + y, nocc_ut + t ) * SEE[ irrep_ut ][ u + num_ut * x ] ); // - f_yt ( 2 delta_ux - Gamma_{ux} )
                           FBB_singlet[ 0 ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] -= value;
                        }
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy + 1 ) ) / 2;
         }
         jump_col += ( num_ut * ( num_ut + 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_singlet[ irrep ] == size_F_singlet[ irrep ] ); // At construction
      const int SIZE = size_B_singlet[ irrep ];
      FBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
      FFF_singlet[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         if ( irrep_t < irrep_u ){
            const int d_t    = indices->getDMRGcumulative( irrep_t );
            const int num_t  = indices->getNDMRG( irrep_t );
            const int nocc_t = indices->getNOCC( irrep_t );
            const int d_u    = indices->getDMRGcumulative( irrep_u );
            const int num_u  = indices->getNDMRG( irrep_u );
            const int nocc_u = indices->getNOCC( irrep_u );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int irrep_y = Irreps::directProd( irrep, irrep_x );
               if ( irrep_x < irrep_y ){
                  const int d_x    = indices->getDMRGcumulative( irrep_x );
                  const int num_x  = indices->getNDMRG( irrep_x );
                  const int nocc_x = indices->getNOCC( irrep_x );
                  const int d_y    = indices->getDMRGcumulative( irrep_y );
                  const int num_y  = indices->getNDMRG( irrep_y );
                  const int nocc_y = indices->getNOCC( irrep_y );
                  const int shift  = jump_row + SIZE * jump_col;

                  /*
                      + 2 delta_ux delta_ty f_dot_1dm
                      - delta_ux f_dot_2dm[ ty ]
                      - delta_ty f_dot_2dm[ ux ]
                      - f_xu ( 2 delta_ty - Gamma_{ty} )
                      - f_yt ( 2 delta_ux - Gamma_{ux} )

                      Because irrep_x < irrep_y = irrep_t < irrep_u is incompatible with irrep_x == irrep_u, these terms are zero
                  */

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           for ( int y = 0; y < num_y; y++ ){

                              // + f_dot_3dm_{utyx} + f_dot_3dm_{utxy}
                              const double fdotsum = ( f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                     + f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );

                              double sbb_value = fdotsum;

                              // FBB: + sum_v f_vu SBB_singlet[ xytv ]
                              for ( int v = 0; v < num_u; v++ ){
                                 sbb_value += fock->get( irrep_u, nocc_u + v, nocc_u + u ) * SBB_singlet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * v ) ];
                              }

                              // FBB: + sum_v f_vt SBB_singlet[ xyvu ]
                              for ( int v = 0; v < num_t; v++ ){
                                 sbb_value += fock->get( irrep_t, nocc_t + v, nocc_t + t ) * SBB_singlet[ irrep ][ shift + x + num_x * y + SIZE * ( v + num_t * u ) ];
                              }

                              // FBB: + sum_w f_yw SBB_singlet[ xwtu ]
                              for ( int w = 0; w < num_y; w++ ){
                                 sbb_value += fock->get( irrep_y, nocc_y + y, nocc_y + w ) * SBB_singlet[ irrep ][ shift + x + num_x * w + SIZE * ( t + num_t * u ) ];
                              }

                              // FBB: + sum_w f_xw SBB_singlet[ wytu ]
                              for ( int w = 0; w < num_x; w++ ){
                                 sbb_value += fock->get( irrep_x, nocc_x + x, nocc_x + w ) * SBB_singlet[ irrep ][ shift + w + num_x * y + SIZE * ( t + num_t * u ) ];
                              }

                              FBB_singlet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] = sbb_value;
                              FFF_singlet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] = fdotsum;
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // num_t == num_x  and  num_u == num_y

                     // 2 delta_uy delta_tx f_dot_1dm
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           FBB_singlet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 2 * f_dot_1dm;
                        }
                     }

                     // - delta_tx f_dot_2dm[ uy ]
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int u = 0; u < num_y; u++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double val_uy = f_dot_2dm[ d_u + u + LAS * ( d_u + y ) ];
                              FBB_singlet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= val_uy;
                           }
                        }
                     }

                     // - delta_uy f_dot_2dm[ tx ]
                     for ( int yu = 0; yu < num_y; yu++ ){
                        for ( int t = 0; t < num_x; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double val_tx = f_dot_2dm[ d_t + t + LAS * ( d_t + x ) ];
                              FBB_singlet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_t * yu ) ] -= val_tx;
                           }
                        }
                     }

                     for ( int t = 0; t < num_x; t++ ){ // num_t == num_x
                        for ( int u = 0; u < num_y; u++ ){ // num_u == num_y
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double value = ( fock->get( irrep_x, nocc_x + x, nocc_x + t ) * SEE[ irrep_y ][ u + num_y * y ]    // - f_xt ( 2 delta_uy - Gamma_{uy} )
                                                      + fock->get( irrep_y, nocc_y + y, nocc_y + u ) * SEE[ irrep_x ][ t + num_x * x ] ); // - f_yu ( 2 delta_tx - Gamma_{tx} )
                                 FBB_singlet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] -= value;
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_SBB_SFF_singlet(){

   /*
      SBB singlet: | SB_tiuj > = ( E_ti E_uj + E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 >  with  i <= j and t <= u

            < SB_xkyl | SB_tiuj > = 2 delta_ik delta_jl SBB_singlet[ It x Iu ][ xytu ]

            SBB_singlet[ It x Iu ][ xytu ] = ( + Gamma_{utyx}
                                               + Gamma_{utxy}
                                               + 2 delta_uy delta_tx
                                               + 2 delta_ux delta_ty
                                               -   delta_uy Gamma_{tx}
                                               -   delta_tx Gamma_{uy}
                                               -   delta_ux Gamma_{ty}
                                               -   delta_ty Gamma_{ux}
                                             )

      SFF singlet: | SF_atbu > = ( E_at E_bu + E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 >  with  a <= b and t <= u

            < SF_cxdy | SF_atbu > = 2 delta_ac delta_bd SFF_singlet[ It x Iu ][ xytu ]

            SFF_singlet[ It x Iu ][ xytu ] = ( + Gamma_{yxut}
                                               + Gamma_{yxtu}
                                             )
   */

   SBB_singlet = new double*[ num_irreps ];
   SFF_singlet = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_singlet[ 0 ] == size_F_singlet[ 0 ] ); // At construction
      const int SIZE = size_B_singlet[ 0 ];
      SBB_singlet[ 0 ] = new double[ SIZE * SIZE ];
      SFF_singlet[ 0 ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
         const int d_ut   = indices->getDMRGcumulative( irrep_ut );
         const int num_ut = indices->getNDMRG( irrep_ut );
         int jump_row = 0;
         for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
            const int d_xy   = indices->getDMRGcumulative( irrep_xy );
            const int num_xy = indices->getNDMRG( irrep_xy );

            // Gamma_{utyx} + Gamma_{utxy}
            for ( int t = 0; t < num_ut; t++ ){
               for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                        const double value = ( two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                             + two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                        const int pointer = jump_row + x + (y*(y+1))/2 + SIZE * ( jump_col + t + (u*(u+1))/2 );
                        SBB_singlet[ 0 ][ pointer ] = value;
                        SFF_singlet[ 0 ][ pointer ] = value;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy

               // + 2 ( delta_uy delta_tx + delta_ux delta_ty ) --> last term only if x <= y = t <= u or hence all equal
               for ( int t = 0; t < num_ut; t++ ){
                  SBB_singlet[ 0 ][ jump_row + t + (t*(t+1))/2 + SIZE * ( jump_col + t + (t*(t+1))/2 ) ] += 4.0;
                  for ( int u = t+1; u < num_ut; u++ ){
                     SBB_singlet[ 0 ][ jump_row + t + (u*(u+1))/2 + SIZE * ( jump_col + t + (u*(u+1))/2 ) ] += 2.0;
                  }
               }

               // - delta_uy Gamma_{tx}
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                     for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                        const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_xy + x ) ];
                        SBB_singlet[ 0 ][ jump_row + x + (uy*(uy+1))/2 + SIZE * ( jump_col + t + (uy*(uy+1))/2 ) ] -= gamma_tx;
                     }
                  }
               }

               // - delta_tx Gamma_{uy}
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                     for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                        const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_xy + y ) ];
                        SBB_singlet[ 0 ][ jump_row + tx + (y*(y+1))/2 + SIZE * ( jump_col + tx + (u*(u+1))/2 ) ] -= gamma_uy;
                     }
                  }
               }

               // - delta_ux Gamma_{ty}
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                     for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                        const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_xy + y ) ];
                        SBB_singlet[ 0 ][ jump_row + ux + (y*(y+1))/2 + SIZE * ( jump_col + t + (ux*(ux+1))/2 ) ] -= gamma_ty;
                     }
                  }
               }

               // - delta_ty Gamma_{ux}
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                     for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                        const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_xy + x ) ];
                        SBB_singlet[ 0 ][ jump_row + x + (ty*(ty+1))/2 + SIZE * ( jump_col + ty + (u*(u+1))/2 ) ] -= gamma_ux;
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy + 1 ) ) / 2;
         }
         jump_col += ( num_ut * ( num_ut + 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_singlet[ irrep ] == size_F_singlet[ irrep ] ); // At construction
      const int SIZE = size_B_singlet[ irrep ];
      SBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
      SFF_singlet[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         if ( irrep_t < irrep_u ){
            const int d_t   = indices->getDMRGcumulative( irrep_t );
            const int num_t = indices->getNDMRG( irrep_t );
            const int d_u   = indices->getDMRGcumulative( irrep_u );
            const int num_u = indices->getNDMRG( irrep_u );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int irrep_y = Irreps::directProd( irrep, irrep_x );
               if ( irrep_x < irrep_y ){
                  const int d_x   = indices->getDMRGcumulative( irrep_x );
                  const int num_x = indices->getNDMRG( irrep_x );
                  const int d_y   = indices->getDMRGcumulative( irrep_y );
                  const int num_y = indices->getNDMRG( irrep_y );

                  /*
                      2 delta_ux delta_ty - delta_ux Gamma_{ty} - delta_ty Gamma_{ux}
                      Because irrep_x < irrep_y = irrep_t < irrep_u is incompatible with irrep_x == irrep_u, these terms are zero
                  */

                  // Gamma_{utyx} + Gamma_{utxy}
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double value = ( two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                   + two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                              const int pointer = jump_row + x + num_x * y + SIZE * ( jump_col + t + num_t * u );
                              SBB_singlet[ irrep ][ pointer ] = value;
                              SFF_singlet[ irrep ][ pointer ] = value;
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // num_t == num_x  and  num_u == num_y

                     // 2 delta_uy delta_tx
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           SBB_singlet[ irrep ][ jump_row + xt + num_x * yu + SIZE * ( jump_col + xt + num_x * yu ) ] += 2.0;
                        }
                     }

                     // - delta_tx Gamma_{uy}
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int u = 0; u < num_y; u++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double gamma_uy = one_rdm[ d_u + u + LAS * ( d_y + y ) ];
                              SBB_singlet[ irrep ][ jump_row + xt + num_x * y + SIZE * ( jump_col + xt + num_x * u ) ] -= gamma_uy;
                           }
                        }
                     }

                     // - delta_uy Gamma_{tx}
                     for ( int yu = 0; yu < num_y; yu++ ){
                        for ( int t = 0; t < num_x; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double gamma_tx = one_rdm[ d_t + t + LAS * ( d_x + x ) ];
                              SBB_singlet[ irrep ][ jump_row + x + num_x * yu + SIZE * ( jump_col + t + num_t * yu ) ] -= gamma_tx;
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_FBB_FFF_triplet(){

   /*
      FBB triplet: < TB_xkyl | ( f_pq E_pq ) | TB_tiuj >
                      = 2 delta_ik delta_jl ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )

      FFF triplet: < TF_cxdy | ( f_pq E_pq ) | TF_atbu >
                      = 2 delta_ac delta_bd ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )

         FBB_triplet[ Iij ][ xytu ] = ( + f_dot_3dm[ utyx ]
                                        - f_dot_3dm[ tuyx ]
                                        + 6 delta_uy delta_tx f_dot_1dm
                                        - 6 delta_ux delta_ty f_dot_1dm
                                        - 3 delta_uy f_dot_2dm[ tx ]
                                        - 3 delta_tx f_dot_2dm[ uy ]
                                        + 3 delta_ux f_dot_2dm[ ty ]
                                        + 3 delta_ty f_dot_2dm[ ux ]
                                        + sum_v f_vu SBB_triplet[ xytv ]
                                        + sum_v f_vt SBB_triplet[ xyvu ]
                                        + sum_w f_yw SBB_triplet[ xwtu ]
                                        + sum_w f_xw SBB_triplet[ wytu ]
                                        - 3 f_xt ( 2 delta_uy - Gamma_{uy} )
                                        - 3 f_yu ( 2 delta_tx - Gamma_{tx} )
                                        + 3 f_xu ( 2 delta_ty - Gamma_{ty} )
                                        + 3 f_yt ( 2 delta_ux - Gamma_{ux} )
                                      )

         FFF_triplet[ Iab ][ xytu ] = ( + f_dot_3dm[ yxut ]
                                        - f_dot_3dm[ yxtu ]
                                      )
   */

   FBB_triplet = new double*[ num_irreps ];
   FFF_triplet = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_triplet[ 0 ] == size_F_triplet[ 0 ] ); // At construction
      const int SIZE = size_B_triplet[ 0 ];
      FBB_triplet[ 0 ] = new double[ SIZE * SIZE ];
      FFF_triplet[ 0 ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
         const int d_ut    = indices->getDMRGcumulative( irrep_ut );
         const int num_ut  = indices->getNDMRG( irrep_ut );
         const int nocc_ut = indices->getNOCC( irrep_ut );
         int jump_row = 0;
         for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
            const int d_xy    = indices->getDMRGcumulative( irrep_xy );
            const int num_xy  = indices->getNDMRG( irrep_xy );
            const int nocc_xy = indices->getNOCC( irrep_xy );
            const int shift   = jump_row + SIZE * jump_col;

            for ( int t = 0; t < num_ut; t++ ){
               for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy

                        // + f_dot_3dm_{utyx} - f_dot_3dm_{utxy}
                        const double fdotdiff = ( f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                                - f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );

                        double sbb_value = fdotdiff;

                        // FBB: + sum_v f_vu SBB_triplet[ xytv ] = + sum_(v>=t) f_vu SBB_triplet[ xytv ] - sum_(v<t) f_vu SBB_triplet[ xyvt ]
                        for ( int v = 0; v < t; v++ ){ // Swap ( t <--> u ) --> minus sign
                           sbb_value -= ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( v + ( t * ( t - 1 ) ) / 2 ) ] );
                        }
                        for ( int v = t+1; v < num_ut; v++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( v * ( v - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_v f_vt SBB_triplet[ xyvu ] = + sum_(v<u) f_vt SBB_triplet[ xyvu ] - sum_(v>=u) f_vt SBB_triplet[ xyuv ]
                        for ( int v = 0; v < u; v++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( v + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int v = u+1; v < num_ut; v++ ){ // Swap ( t <--> u ) --> minus sign
                           sbb_value -= ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( u + ( v * ( v - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_yw SBB_triplet[ xwtu ] = + sum_(w>=x) f_yw SBB_triplet[ xwtu ] - sum_(w<x) f_yw SBB_triplet[ wxtu ]
                        for ( int w = 0; w < x; w++ ){ // Swap ( x <--> y ) --> minus sign
                           sbb_value -= ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_triplet[ 0 ][ shift + w + ( x * ( x - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int w = x+1; w < num_xy; w++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_triplet[ 0 ][ shift + x + ( w * ( w - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_xw SBB_triplet[ wytu ] = + sum_(w<y) f_xw SBB_triplet[ wytu ] - sum_(w>=y) f_xw SBB_triplet[ ywtu ]
                        for ( int w = 0; w < y; w++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_triplet[ 0 ][ shift + w + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int w = y+1; w < num_xy; w++ ){ // Swap ( x <--> y ) --> minus sign
                           sbb_value -= ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_triplet[ 0 ][ shift + y + ( w * ( w - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }

                        FBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] = sbb_value;
                        FFF_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] = fdotdiff;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy and d_ut = d_xy

               // FBB: + 6 ( delta_uy delta_tx - delta_ux delta_ty ) f_dot_1dm --> last term only if x < y = t < u or NEVER
               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t+1; u < num_ut; u++ ){
                     FBB_triplet[ 0 ][ shift + t + ( u * ( u - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] += 6 * f_dot_1dm;
                  }
               }

               // FBB: - 3 delta_uy f_dot_2dm[ tx ]
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                     for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                        const double val_tx = f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ];
                        FBB_triplet[ 0 ][ shift + x + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy - 1 ) ) / 2 ) ] -= 3 * val_tx;
                     }
                  }
               }

               // FBB: - 3 delta_tx f_dot_2dm[ uy ]
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                     for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                        const double val_uy = f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ];
                        FBB_triplet[ 0 ][ shift + tx + ( y * ( y - 1 ) ) / 2 + SIZE * ( tx + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * val_uy;
                     }
                  }
               }

               // FBB: + 3 delta_ux f_dot_2dm[ ty ]
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                     for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                        const double val_ty = f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ];
                        FBB_triplet[ 0 ][ shift + ux + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux - 1 ) ) / 2 ) ] += 3 * val_ty;
                     }
                  }
               }

               // FBB: + 3 delta_ty f_dot_2dm[ ux ]
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                     for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                        const double val_ux = f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ];
                        FBB_triplet[ 0 ][ shift + x + ( ty * ( ty - 1 ) ) / 2 + SIZE * ( ty + ( u * ( u - 1 ) ) / 2 ) ] += 3 * val_ux;
                     }
                  }
               }

               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                     for ( int x = 0; x < num_ut; x++ ){
                        for ( int y = x+1; y < num_ut; y++ ){ // 0 <= x < y < num_xy == num_ut
                           const double value = ( fock->get( irrep_ut, nocc_ut + x, nocc_ut + t ) * SEE[ irrep_ut ][ u + num_ut * y ]    // - 3 f_xt ( 2 delta_uy - Gamma_{uy} )
                                                + fock->get( irrep_ut, nocc_ut + y, nocc_ut + u ) * SEE[ irrep_ut ][ t + num_ut * x ]    // - 3 f_yu ( 2 delta_tx - Gamma_{tx} )
                                                - fock->get( irrep_ut, nocc_ut + x, nocc_ut + u ) * SEE[ irrep_ut ][ t + num_ut * y ]    // + 3 f_xu ( 2 delta_ty - Gamma_{ty} )
                                                - fock->get( irrep_ut, nocc_ut + y, nocc_ut + t ) * SEE[ irrep_ut ][ u + num_ut * x ] ); // + 3 f_yt ( 2 delta_ux - Gamma_{ux} )
                           FBB_triplet[ 0 ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * value;
                        }
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy - 1 ) ) / 2;
         }
         jump_col += ( num_ut * ( num_ut - 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_triplet[ irrep ] == size_F_triplet[ irrep ] ); // At construction
      const int SIZE = size_B_triplet[ irrep ];
      FBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
      FFF_triplet[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         if ( irrep_t < irrep_u ){
            const int d_t    = indices->getDMRGcumulative( irrep_t );
            const int num_t  = indices->getNDMRG( irrep_t );
            const int nocc_t = indices->getNOCC( irrep_t );
            const int d_u    = indices->getDMRGcumulative( irrep_u );
            const int num_u  = indices->getNDMRG( irrep_u );
            const int nocc_u = indices->getNOCC( irrep_u );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int irrep_y = Irreps::directProd( irrep, irrep_x );
               if ( irrep_x < irrep_y ){
                  const int d_x    = indices->getDMRGcumulative( irrep_x );
                  const int num_x  = indices->getNDMRG( irrep_x );
                  const int nocc_x = indices->getNOCC( irrep_x );
                  const int d_y    = indices->getDMRGcumulative( irrep_y );
                  const int num_y  = indices->getNDMRG( irrep_y );
                  const int nocc_y = indices->getNOCC( irrep_y );
                  const int shift  = jump_row + SIZE * jump_col;

                  /*
                      - 6 delta_ux delta_ty f_dot_1dm
                      + 3 delta_ux f_dot_2dm[ ty ]
                      + 3 delta_ty f_dot_2dm[ ux ]
                      + 3 f_xu ( 2 delta_ty - Gamma_{ty} )
                      + 3 f_yt ( 2 delta_ux - Gamma_{ux} )

                      Because irrep_x < irrep_y = irrep_t < irrep_u is incompatible with irrep_x == irrep_u, these terms are zero
                  */

                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           for ( int y = 0; y < num_y; y++ ){

                              // + f_dot_3dm_{utyx} - f_dot_3dm_{utxy}
                              const double fdotdiff = ( f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                      - f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );

                              double sbb_value = fdotdiff;

                              // FBB: + sum_v f_vu SBB_triplet[ xytv ]
                              for ( int v = 0; v < num_u; v++ ){
                                 sbb_value += fock->get( irrep_u, nocc_u + v, nocc_u + u ) * SBB_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * v ) ];
                              }

                              // FBB: + sum_v f_vt SBB_triplet[ xyvu ]
                              for ( int v = 0; v < num_t; v++ ){
                                 sbb_value += fock->get( irrep_t, nocc_t + v, nocc_t + t ) * SBB_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( v + num_t * u ) ];
                              }

                              // FBB: + sum_w f_yw SBB_triplet[ xwtu ]
                              for ( int w = 0; w < num_y; w++ ){
                                 sbb_value += fock->get( irrep_y, nocc_y + y, nocc_y + w ) * SBB_triplet[ irrep ][ shift + x + num_x * w + SIZE * ( t + num_t * u ) ];
                              }

                              // FBB: + sum_w f_xw SBB_triplet[ wytu ]
                              for ( int w = 0; w < num_x; w++ ){
                                 sbb_value += fock->get( irrep_x, nocc_x + x, nocc_x + w ) * SBB_triplet[ irrep ][ shift + w + num_x * y + SIZE * ( t + num_t * u ) ];
                              }

                              FBB_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] = sbb_value;
                              FFF_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] = fdotdiff;
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // num_t == num_x  and  num_u == num_y

                     // + 6 delta_uy delta_tx f_dot_1dm
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           FBB_triplet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 6 * f_dot_1dm;
                        }
                     }

                     // - 3 delta_tx f_dot_2dm[ uy ]
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int u = 0; u < num_y; u++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double val_uy = f_dot_2dm[ d_u + u + LAS * ( d_u + y ) ];
                              FBB_triplet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= 3 * val_uy;
                           }
                        }
                     }

                     // - 3 delta_uy f_dot_2dm[ tx ]
                     for ( int yu = 0; yu < num_y; yu++ ){
                        for ( int t = 0; t < num_x; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double val_tx = f_dot_2dm[ d_t + t + LAS * ( d_t + x ) ];
                              FBB_triplet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_x * yu ) ] -= 3 * val_tx;
                           }
                        }
                     }

                     for ( int t = 0; t < num_x; t++ ){ // num_t == num_x
                        for ( int u = 0; u < num_y; u++ ){ // num_u == num_y
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double value = ( fock->get( irrep_x, nocc_x + x, nocc_x + t ) * SEE[ irrep_y ][ u + num_y * y ]    // - 3 f_xt ( 2 delta_uy - Gamma_{uy} )
                                                      + fock->get( irrep_y, nocc_y + y, nocc_y + u ) * SEE[ irrep_x ][ t + num_x * x ] ); // - 3 f_yu ( 2 delta_tx - Gamma_{tx} )
                                 FBB_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] -= 3 * value;
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_SBB_SFF_triplet(){

   /*
      SBB triplet: | TB_tiuj > = ( E_ti E_uj - E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 >  with  i < j and t < u

            < TB_xkyl | TB_tiuj > = 2 delta_ik delta_jl SBB_triplet[ It x Iu ][ xytu ]

            SBB_triplet[ It x Iu ][ xytu ] = ( + Gamma_{utyx}
                                               - Gamma_{utxy}
                                               + 6 delta_uy delta_tx
                                               - 6 delta_ux delta_ty
                                               - 3 delta_uy Gamma_{tx}
                                               - 3 delta_tx Gamma_{uy}
                                               + 3 delta_ux Gamma_{ty}
                                               + 3 delta_ty Gamma_{ux}
                                             )

      SFF triplet: | TF_atbu > = ( E_at E_bu - E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 >  with  a < b and t < u

            < TF_cxdy | TF_atbu > = 2 delta_ac delta_bd SFF_triplet[ It x Iu ][ xytu ]

            SFF_triplet[ It x Iu ][ xytu ] = ( + Gamma_{yxut}
                                               - Gamma_{yxtu}
                                             )
   */

   SBB_triplet = new double*[ num_irreps ];
   SFF_triplet = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_triplet[ 0 ] == size_F_triplet[ 0 ] ); // At construction
      const int SIZE = size_B_triplet[ 0 ];
      SBB_triplet[ 0 ] = new double[ SIZE * SIZE ];
      SFF_triplet[ 0 ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
         const int d_ut   = indices->getDMRGcumulative( irrep_ut );
         const int num_ut = indices->getNDMRG( irrep_ut );
         int jump_row = 0;
         for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
            const int d_xy   = indices->getDMRGcumulative( irrep_xy );
            const int num_xy = indices->getNDMRG( irrep_xy );

            // Gamma_{utyx} - Gamma_{utxy}
            for ( int t = 0; t < num_ut; t++ ){
               for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                        const double value = ( two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                             - two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                        const int pointer = jump_row + x + (y*(y-1))/2 + SIZE * ( jump_col + t + (u*(u-1))/2 );
                        SBB_triplet[ 0 ][ pointer ] = value;
                        SFF_triplet[ 0 ][ pointer ] = value;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal --> num_ut == num_xy

               // + 6 ( delta_uy delta_tx - delta_ux delta_ty ) --> last term only if x < y = t < u or NEVER
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int uy = tx+1; uy < num_ut; uy++ ){
                     SBB_triplet[ 0 ][ jump_row + tx + (uy*(uy-1))/2 + SIZE * ( jump_col + tx + (uy*(uy-1))/2 ) ] += 6.0;
                  }
               }

               // - 3 delta_uy Gamma_{tx}
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                     for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                        const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_xy + x ) ];
                        SBB_triplet[ 0 ][ jump_row + x + (uy*(uy-1))/2 + SIZE * ( jump_col + t + (uy*(uy-1))/2 ) ] -= 3 * gamma_tx;
                     }
                  }
               }

               // - 3 delta_tx Gamma_{uy}
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                     for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                        const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_xy + y ) ];
                        SBB_triplet[ 0 ][ jump_row + tx + (y*(y-1))/2 + SIZE * ( jump_col + tx + (u*(u-1))/2 ) ] -= 3 * gamma_uy;
                     }
                  }
               }

               // + 3 delta_ux Gamma_{ty}
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                     for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                        const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_xy + y ) ];
                        SBB_triplet[ 0 ][ jump_row + ux + (y*(y-1))/2 + SIZE * ( jump_col + t + (ux*(ux-1))/2 ) ] += 3 * gamma_ty;
                     }
                  }
               }

               // + 3 delta_ty Gamma_{ux}
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                     for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                        const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_xy + x ) ];
                        SBB_triplet[ 0 ][ jump_row + x + (ty*(ty-1))/2 + SIZE * ( jump_col + ty + (u*(u-1))/2 ) ] += 3 * gamma_ux;
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy - 1 ) ) / 2;
         }
         jump_col += ( num_ut * ( num_ut - 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_triplet[ irrep ] == size_F_triplet[ irrep ] ); // At construction
      const int SIZE = size_B_triplet[ irrep ];
      SBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
      SFF_triplet[ irrep ] = new double[ SIZE * SIZE ];

      int jump_col = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         const int irrep_u = Irreps::directProd( irrep, irrep_t );
         if ( irrep_t < irrep_u ){
            const int d_t   = indices->getDMRGcumulative( irrep_t );
            const int num_t = indices->getNDMRG( irrep_t );
            const int d_u   = indices->getDMRGcumulative( irrep_u );
            const int num_u = indices->getNDMRG( irrep_u );
            int jump_row = 0;
            for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
               const int irrep_y = Irreps::directProd( irrep, irrep_x );
               if ( irrep_x < irrep_y ){
                  const int d_x   = indices->getDMRGcumulative( irrep_x );
                  const int num_x = indices->getNDMRG( irrep_x );
                  const int d_y   = indices->getDMRGcumulative( irrep_y );
                  const int num_y = indices->getNDMRG( irrep_y );

                  /*
                      - 6 delta_ux delta_ty + 3 delta_ux Gamma_{ty} + 3 delta_ty Gamma_{ux}
                      Because irrep_x < irrep_y = irrep_t < irrep_u is incompatible with irrep_x == irrep_u, these terms are zero
                  */

                  // Gamma_{utyx} - Gamma_{utxy}
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int u = 0; u < num_u; u++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double value = ( two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                   - two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                              const int pointer = jump_row + x + num_x * y + SIZE * ( jump_col + t + num_t * u );
                              SBB_triplet[ irrep ][ pointer ] = value;
                              SFF_triplet[ irrep ][ pointer ] = value;
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // --> num_t == num_x   and   num_u == num_y

                     // + 6 delta_uy delta_tx
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           SBB_triplet[ irrep ][ jump_row + xt + num_x * yu + SIZE * ( jump_col + xt + num_x * yu ) ] += 6.0;
                        }
                     }

                     // - 3 delta_tx Gamma_{uy}
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int u = 0; u < num_y; u++ ){
                           for ( int y = 0; y < num_y; y++ ){
                              const double gamma_uy = one_rdm[ d_u + u + LAS * ( d_y + y ) ];
                              SBB_triplet[ irrep ][ jump_row + xt + num_x * y + SIZE * ( jump_col + xt + num_x * u ) ] -= 3 * gamma_uy;
                           }
                        }
                     }

                     // - 3 delta_uy Gamma_{tx}
                     for ( int yu = 0; yu < num_y; yu++ ){
                        for ( int t = 0; t < num_x; t++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double gamma_tx = one_rdm[ d_t + t + LAS * ( d_x + x ) ];
                              SBB_triplet[ irrep ][ jump_row + x + num_x * yu + SIZE * ( jump_col + t + num_x * yu ) ] -= 3 * gamma_tx;
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_FEE_FGG(){

   /*

      FEE singlet: | SE_tiaj > = ( E_ti E_aj + E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <= j
      FEE triplet: | TE_tiaj > = ( E_ti E_aj - E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <  j

         < SE_ukbl | f_pq E_pq | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
         < TE_ukbl | f_pq E_pq | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )

         FEE[ It ][ ut ] = ( + 2 * delta_ut * f_dot_1dm
                             - f_dot_2dm[ It ][ tu ]
                             + 2 * f_ut
                             - sum_v Gamma_uv f_vt
                             - sum_w f_uw Gamma_wt
                           )

      FGG singlet: | SG_aibt > = ( E_ai E_bt + E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <= b
      FGG triplet: | TG_aibt > = ( E_ai E_bt - E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <  b

         < SG_cjdu | f_pq E_pq | SG_aibt > = 2 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
         < TG_cjdu | f_pq E_pq | TG_aibt > = 6 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )

         FGG[ It ][ ut ] = ( + f_dot_2dm[ It ][ ut ]
                           )
   */

   FEE = new double*[ num_irreps ];
   FGG = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int SIZE = indices->getNDMRG( irrep );
      FEE[ irrep ] = new double[ SIZE * SIZE ];
      FGG[ irrep ] = new double[ SIZE * SIZE ];
      const int jumpx = indices->getDMRGcumulative( irrep );
      const int NOCC  = indices->getNOCC( irrep );

      //FGG
      for ( int t = 0; t < SIZE; t++ ){
         for ( int u = 0; u < SIZE; u++ ){
            FGG[ irrep ][ u + SIZE * t ] = f_dot_2dm[ jumpx + u + LAS * ( jumpx + t ) ];
         }
      }

      //FEE
      for ( int t = 0; t < SIZE; t++ ){
         for ( int u = 0; u < SIZE; u++ ){

            double value = 2 * fock->get( irrep, NOCC + u, NOCC + t ) - f_dot_2dm[ jumpx + u + LAS * ( jumpx + t ) ];
            for ( int vw = 0; vw < SIZE; vw++ ){
               value -= ( one_rdm[ jumpx + u + LAS * ( jumpx + vw ) ] * fock->get( irrep, NOCC + vw, NOCC + t )
                        + fock->get( irrep, NOCC + u, NOCC + vw ) * one_rdm[ jumpx + vw + LAS * ( jumpx + t ) ] );
            }
            FEE[ irrep ][ u + SIZE * t ] = value;
         }
         FEE[ irrep ][ t + SIZE * t ] += 2 * f_dot_1dm;
      }
   }

}

void CheMPS2::CASPT2::make_SEE_SGG(){

   /*
      SEE singlet: | SE_tiaj > = ( E_ti E_aj + E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <= j
      SEE triplet: | TE_tiaj > = ( E_ti E_aj - E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <  j

            < SE_ukbl | SE_tiaj > = 2 delta_ab delta_ik delta_jl SEE[ It ][ ut ]
            < TE_ukbl | TE_tiaj > = 6 delta_ab delta_ik delta_jl SEE[ It ][ ut ]

            SEE[ It ][ ut ] = ( + 2 delta_tu
                                - Gamma_tu
                              )

      SGG singlet: | SG_aibt > = ( E_ai E_bt + E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <= b
      SGG triplet: | TG_aibt > = ( E_ai E_bt - E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <  b

            < SG_cjdu | SG_aibt > = 2 delta_ij delta_ac delta_bd SGG[ It ][ ut ]
            < TG_cjdu | TG_aibt > = 6 delta_ij delta_ac delta_bd SGG[ It ][ ut ]

            SGG[ It ][ ut ] = ( + Gamma_ut
                              )
   */

   SEE = new double*[ num_irreps ];
   SGG = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){

      const int SIZE = indices->getNDMRG( irrep_ut );
      SEE[ irrep_ut ] = new double[ SIZE * SIZE ];
      SGG[ irrep_ut ] = new double[ SIZE * SIZE ];
      const int d_ut = indices->getDMRGcumulative( irrep_ut );

      for ( int t = 0; t < SIZE; t++ ){
         for ( int u = 0; u < SIZE; u++ ){
            const double gamma_ut = one_rdm[ d_ut + u + LAS * ( d_ut + t ) ];
            SEE[ irrep_ut ][ u + SIZE * t ] = - gamma_ut;
            SGG[ irrep_ut ][ u + SIZE * t ] =   gamma_ut;
         }
         SEE[ irrep_ut ][ t + SIZE * t ] += 2.0;
      }
   }

}



