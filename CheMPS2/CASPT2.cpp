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

   make_FAA_FCC(); // Needs to be called AFTER make_S**()!
   make_FDD(); // Needs to be called AFTER make_S**()!
   make_FEE_FGG(); // Needs to be called AFTER make_S**()!
   make_FBB_FFF_singlet(); // Needs to be called AFTER make_S**()!
   make_FBB_FFF_triplet(); // Needs to be called AFTER make_S**()!
   
   make_FEH_FGH(); // Needs to be called AFTER make_S**()!

   delete [] f_dot_3dm;
   delete [] f_dot_2dm;

   recreate(); // Remove the overlap matrices

   {

      int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];

      double * diag_fock = new double[ total_size ];
      diagonal( diag_fock, -E_FOCK );
      
      double min_eig = 1e18;
      for ( int elem = 0; elem < total_size; elem++ ){ min_eig = min( min_eig, diag_fock[ elem ] ); }
      cout << "Minimum eigenvalue diag(FOCK) = " << min_eig << endl;

      // Calculate E(CASPT2-D) = - < Psi0 | H P_SD [ blockdiag(F) - E_FOCK * S ]^{-1} P_SD H | Psi0 >
      double energy_caspt2d = 0.0;
      for ( int elem = 0; elem < total_size; elem++ ){ energy_caspt2d -= vector_rhs[ elem ] * vector_rhs[ elem ] / diag_fock[ elem ]; }
      
      delete [] diag_fock;
      
      cout.precision(8);
      cout << std::fixed;
      cout << "E(CASPT2-D)                = " << energy_caspt2d << endl;
      cout << "Test 8 according to molcas = " << -0.1596306078 << endl;
      cout.unsetf( std::ios::floatfield );
      cout.precision(15);

   }

}

CheMPS2::CASPT2::~CASPT2(){

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      delete [] FAA[ irrep ];
      delete [] FCC[ irrep ];
      delete [] FDD[ irrep ];
      delete [] FEE[ irrep ];
      delete [] FGG[ irrep ];
      delete [] FEH[ irrep ];
      delete [] FGH[ irrep ];
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
   delete [] FEH;
   delete [] FGH;
   delete [] FBB_singlet;
   delete [] FBB_triplet;
   delete [] FFF_singlet;
   delete [] FFF_triplet;

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

   /*** Type H singlet : c_aibj ( E_ai E_bj + E_bi E_aj ) / sqrt( 2 - ( 1 - delta_ij ) * ( 1 - delta_ab ) ) | 0 > with a <= b and i <= j
                         c_aibj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + count_aibj ]
        Type H triplet : c_aibj ( E_ai E_bj - E_bi E_aj ) / sqrt( 2 - ( 1 - delta_ij ) * ( 1 - delta_ab ) ) | 0 > with a <  b and i <  j
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
         const long long nvir2 = indices->getNVIRT( i2 );
         for ( int i3 = 0; i3 < num_irreps; i3++ ){
            const int i4 = Irreps::directProd( Irreps::directProd( i1, i2 ), i3 );
            const long long nocc3 = indices->getNOCC( i3 );
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

int CheMPS2::CASPT2::recreatehelper( double * FOCK, double * OVLP, int SIZE, double * work, double * eigs, int lwork ){

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

void CheMPS2::CASPT2::recreatehelper2( double * OVLP, int OLDSIZE, int NEWSIZE, double * rhs_old, double * rhs_new, const int num_rhs ){

   // rhs <-- eigs^{-0.5} V^T rhs
   int inc1 = 1;
   double set = 0.0;
   double one = 1.0;
   char trans = 'T';
   for ( int sector = 0; sector < num_rhs; sector++ ){
      dgemv_( &trans, &OLDSIZE, &NEWSIZE, &one, OVLP, &OLDSIZE, rhs_old + OLDSIZE * sector, &inc1, &set, rhs_new + NEWSIZE * sector, &inc1 );
   }

}

void CheMPS2::CASPT2::recreatehelper_left( double * ROT, int OLDSIZE, int NEWSIZE, double * matrix, int RIGHTSIZE, double * work ){

   if ( NEWSIZE * OLDSIZE * RIGHTSIZE == 0 ){ return; }

   // work = ROT^T * matrix
   double set = 0.0;
   double one = 1.0;
   char trans   = 'T';
   char notrans = 'N';
   dgemm_( &trans, &notrans, &NEWSIZE, &RIGHTSIZE, &OLDSIZE, &one, ROT, &OLDSIZE, matrix, &OLDSIZE, &set, work, &NEWSIZE );

   // matrix = work
   int size_copy = NEWSIZE * RIGHTSIZE;
   int inc1 = 1;
   dcopy_( &size_copy, work, &inc1, matrix, &inc1 );

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
      newsize_A[ irrep ] = recreatehelper( FAA[ irrep ], SAA[ irrep ], size_A[ irrep ], work, eigs, lwork );
      newsize_C[ irrep ] = recreatehelper( FCC[ irrep ], SCC[ irrep ], size_C[ irrep ], work, eigs, lwork );
      newsize_D[ irrep ] = recreatehelper( FDD[ irrep ], SDD[ irrep ], size_D[ irrep ], work, eigs, lwork );
      newsize_E[ irrep ] = recreatehelper( FEE[ irrep ], SEE[ irrep ], size_E[ irrep ], work, eigs, lwork );
      newsize_G[ irrep ] = recreatehelper( FGG[ irrep ], SGG[ irrep ], size_G[ irrep ], work, eigs, lwork );
      newsize_B_singlet[ irrep ] = recreatehelper( FBB_singlet[ irrep ], SBB_singlet[ irrep ], size_B_singlet[ irrep ], work, eigs, lwork );
      newsize_B_triplet[ irrep ] = recreatehelper( FBB_triplet[ irrep ], SBB_triplet[ irrep ], size_B_triplet[ irrep ], work, eigs, lwork );
      newsize_F_singlet[ irrep ] = recreatehelper( FFF_singlet[ irrep ], SFF_singlet[ irrep ], size_F_singlet[ irrep ], work, eigs, lwork );
      newsize_F_triplet[ irrep ] = recreatehelper( FFF_triplet[ irrep ], SFF_triplet[ irrep ], size_F_triplet[ irrep ], work, eigs, lwork );
   }

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      recreatehelper_left( SEE[ irrep ], size_E[ irrep ], newsize_E[ irrep ], FEH[ irrep ], indices->getNVIRT( irrep ), work );
      recreatehelper_left( SGG[ irrep ], size_G[ irrep ], newsize_G[ irrep ], FGH[ irrep ], indices->getNOCC( irrep ),  work );
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
         recreatehelper2( SAA[ irrep ], size_A[ irrep ], newsize_A[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_B_singlet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_B_singlet[ irrep ];
         assert( num_rhs * size_B_singlet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_B_singlet[ irrep ];
         recreatehelper2( SBB_singlet[ irrep ], size_B_singlet[ irrep ], newsize_B_singlet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_B_triplet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_B_triplet[ irrep ];
         assert( num_rhs * size_B_triplet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_B_triplet[ irrep ];
         recreatehelper2( SBB_triplet[ irrep ], size_B_triplet[ irrep ], newsize_B_triplet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_C[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_C;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_C[ irrep ];
         assert( num_rhs * size_C[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_C[ irrep ];
         recreatehelper2( SCC[ irrep ], size_C[ irrep ], newsize_C[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_D[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_D;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_D[ irrep ];
         assert( num_rhs * size_D[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_D[ irrep ];
         recreatehelper2( SDD[ irrep ], size_D[ irrep ], newsize_D[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_E[ irrep ] > 0 ){
         const int ptr1 = irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET;
         const int num_rhs1 = ( jump[ ptr1 + 1 ] - jump[ ptr1 ] ) / size_E[ irrep ];
         assert( num_rhs1 * size_E[ irrep ] == jump[ ptr1 + 1 ] - jump[ ptr1 ] );
         helper[ ptr1 ] = num_rhs1 * newsize_E[ irrep ];
         recreatehelper2( SEE[ irrep ], size_E[ irrep ], newsize_E[ irrep ], vector_rhs + jump[ ptr1 ], tempvector_rhs + jump[ ptr1 ], num_rhs1 );
         const int ptr2 = irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET;
         const int num_rhs2 = ( jump[ ptr2 + 1 ] - jump[ ptr2 ] ) / size_E[ irrep ];
         assert( num_rhs2 * size_E[ irrep ] == jump[ ptr2 + 1 ] - jump[ ptr2 ] );
         helper[ ptr2 ] = num_rhs2 * newsize_E[ irrep ];
         recreatehelper2( SEE[ irrep ], size_E[ irrep ], newsize_E[ irrep ], vector_rhs + jump[ ptr2 ], tempvector_rhs + jump[ ptr2 ], num_rhs2 );
      }

      if ( newsize_F_singlet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_F_singlet[ irrep ];
         assert( num_rhs * size_F_singlet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_F_singlet[ irrep ];
         recreatehelper2( SFF_singlet[ irrep ], size_F_singlet[ irrep ], newsize_F_singlet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_F_triplet[ irrep ] > 0 ){
         const int ptr = irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET;
         const int num_rhs = ( jump[ ptr + 1 ] - jump[ ptr ] ) / size_F_triplet[ irrep ];
         assert( num_rhs * size_F_triplet[ irrep ] == jump[ ptr + 1 ] - jump[ ptr ] );
         helper[ ptr ] = num_rhs * newsize_F_triplet[ irrep ];
         recreatehelper2( SFF_triplet[ irrep ], size_F_triplet[ irrep ], newsize_F_triplet[ irrep ], vector_rhs + jump[ ptr ], tempvector_rhs + jump[ ptr ], num_rhs );
      }

      if ( newsize_G[ irrep ] > 0 ){
         const int ptr1 = irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET;
         const int num_rhs1 = ( jump[ ptr1 + 1 ] - jump[ ptr1 ] ) / size_G[ irrep ];
         assert( num_rhs1 * size_G[ irrep ] == jump[ ptr1 + 1 ] - jump[ ptr1 ] );
         helper[ ptr1 ] = num_rhs1 * newsize_G[ irrep ];
         recreatehelper2( SGG[ irrep ], size_G[ irrep ], newsize_G[ irrep ], vector_rhs + jump[ ptr1 ], tempvector_rhs + jump[ ptr1 ], num_rhs1 );
         const int ptr2 = irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET;
         const int num_rhs2 = ( jump[ ptr2 + 1 ] - jump[ ptr2 ] ) / size_G[ irrep ];
         assert( num_rhs2 * size_G[ irrep ] == jump[ ptr2 + 1 ] - jump[ ptr2 ] );
         helper[ ptr2 ] = num_rhs2 * newsize_G[ irrep ];
         recreatehelper2( SGG[ irrep ], size_G[ irrep ], newsize_G[ irrep ], vector_rhs + jump[ ptr2 ], tempvector_rhs + jump[ ptr2 ], num_rhs2 );
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
      maxsize = max( max( max( max( max( max( max( max( max( max( max(
                                        size_A[irrep],
                                        size_C[irrep] ),
                                        size_D[irrep] ),
                                        size_E[irrep] ),
                                        size_G[irrep] ),
                                        size_B_singlet[irrep] ),
                                        size_B_triplet[irrep] ),
                                        size_F_singlet[irrep] ),
                                        size_F_triplet[irrep] ),
                                        indices->getNOCC( irrep ) ),
                                        indices->getNVIRT( irrep ) ),
                                        maxsize );
   }
   if ( maxsize <= 2 ){ maxsize = 3; }
   return maxsize;

}

void CheMPS2::CASPT2::matvec( double * vector, double * result, double * diag_fock ) const{

   /*
         TODO  | A  Bsinglet  Btriplet  C     D1     D2    Esinglet  Etriplet  Fsinglet  Ftriplet  Gsinglet  Gtriplet  Hsinglet  Htriplet
      ---------+-------------------------------------------------------------------------------------------------------------------------
      A        | OK    x       x        0     x      x     GRAD      GRAD      0         0         0         0         0         0
      Bsinglet | x     OK      0        0     0      0     x         0         0         0         0         0         0         0
      Btriplet | x     0       OK       0     0      0     0         x         0         0         0         0         0         0
      C        | 0     0       0        OK    x      x     0         0         x         x         GRAD      GRAD      0         0
      D1       | x     0       0        x     OK     OK    x         x         0         0         x         x         GRAD      GRAD
      D2       | x     0       0        x     OK     OK    x         x         0         0         x         x         GRAD      GRAD
      Esinglet | GRAD  x       0        0     x      x     OK        0         0         0         0         0         xy        0
      Etriplet | GRAD  0       x        0     x      x     0         OK        0         0         0         0         0         xy
      Fsinglet | 0     0       0        x     0      0     0         0         OK        0         x         0         0         0
      Ftriplet | 0     0       0        x     0      0     0         0         0         OK        0         x         0         0
      Gsinglet | 0     0       0        GRAD  x      x     0         0         x         0         OK        0         xy        0
      Gtriplet | 0     0       0        GRAD  x      x     0         0         0         x         0         OK        0         xy
      Hsinglet | 0     0       0        0     GRAD   GRAD  xy        0         0         0         xy        0         OK        0
      Htriplet | 0     0       0        0     GRAD   GRAD  0         xy        0         0         0         xy        0         OK
      
   */

   const int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
   for ( int elem = 0; elem < total_size; elem++ ){ result[ elem ] = diag_fock[ elem ] * vector[ elem ]; }

   /*
      FEH singlet: < SE_xkcl | F | SH_aibj > = +2 delta_ik delta_jl ( delta_bc FEH[ Ix ][ xa ] + delta_ac FEH[ Ix ][ xb ] )
      FEH triplet: < TE_xkcl | F | TH_aibj > = +6 delta_ik delta_jl ( delta_bc FEH[ Ix ][ xa ] - delta_ac FEH[ Ix ][ xb ] )
      FGH singlet: < SG_ckdx | F | SH_aibj > = -2 delta_ac delta_bd ( delta_ik FGH[ Ix ][ xj ] + delta_jk FGH[ Ix ][ xi ] )
      FGH triplet: < TG_ckdx | F | TH_aibj > = -6 delta_ac delta_bd ( delta_ik FGH[ Ix ][ xj ] - delta_jk FGH[ Ix ][ xi ] )
   */

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
         int jump_ai = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_a = Irreps::directProd( irrep_i, irrep );
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            for ( int i = 0; i < NOCC_i; i++ ){
               const double f_ii = fock->get( irrep_i, i, i );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  const double beta = shifted_prefactor + f_aa - f_ii;
                  const int count = jump_ai + i + NOCC_i * a;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FDD[ irrep ][ elem ] + beta; }
               }
            }
            jump_ai += NOCC_i * NVIR_a;
         }
      }
   }

   // FBB singlet: < SB_xkyl | f_pq E_pq | SB_tiuj > = 2 delta_ik delta_jl ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )
   {
      const int SIZE = size_B_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               const double f_ii = fock->get( irrep_ij, i, i );
               for ( int j = i; j < nocc_ij; j++ ){
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = shifted_prefactor - f_ii - f_jj;
                  const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ 0 ][ elem ] + beta ); }
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_B_singlet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  const double f_ii = fock->get( irrep_i, i, i );
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = shifted_prefactor - f_ii - f_jj;
                     const int count = jump_ij + i + nocc_i * j;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ irrep ][ elem ] + beta ); }
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }

   // FBB triplet: < TB_xkyl | f_pq E_pq | TB_tiuj > = 2 delta_ik delta_jl ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )
   {
      const int SIZE = size_B_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               const double f_ii = fock->get( irrep_ij, i, i );
               for ( int j = i+1; j < nocc_ij; j++ ){
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = shifted_prefactor - f_ii - f_jj;
                  const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ 0 ][ elem ] + beta ); }
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_B_triplet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  const double f_ii = fock->get( irrep_i, i, i );
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = shifted_prefactor - f_ii - f_jj;
                     const int count = jump_ij + i + nocc_i * j;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ irrep ][ elem ] + beta ); }
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }

   // FFF singlet: < SF_cxdy | f_pq E_pq | SF_atbu > = 2 delta_ac delta_bd ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )
   {
      const int SIZE = size_F_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a; b < NVIR_ab; b++ ){
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = shifted_prefactor + f_aa + f_bb;
                  const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ 0 ][ elem ] + beta ); }
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_F_singlet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int irrep_b = Irreps::directProd( irrep, irrep_a );
            if ( irrep_a < irrep_b ){
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int NVIR_b = indices->getNVIRT( irrep_b );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = shifted_prefactor + f_aa + f_bb;
                     const int count = jump_ab + a + NVIR_a * b;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ irrep ][ elem ] + beta ); }
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }

   // FFF triplet: < TF_cxdy | f_pq E_pq | TF_atbu > = 2 delta_ac delta_bd ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )
   {
      const int SIZE = size_F_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a+1; b < NVIR_ab; b++ ){
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = shifted_prefactor + f_aa + f_bb;
                  const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ 0 ][ elem ] + beta ); }
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      const int SIZE = size_F_triplet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int irrep_b = Irreps::directProd( irrep, irrep_a );
            if ( irrep_a < irrep_b ){
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int NVIR_b = indices->getNVIRT( irrep_b );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = shifted_prefactor + f_aa + f_bb;
                     const int count = jump_ab + a + NVIR_a * b;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ irrep ][ elem ] + beta ); }
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }

   // FEE singlet: < SE_ukbl | f_pq E_pq | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_E[ irrep ];
      if ( SIZE > 0 ){
         int jump_aij = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            const int irrep_occ = Irreps::directProd( irrep_a, irrep );
            int jump_ij = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
               if ( irrep_i == irrep_j ){
                  const int NOCC_ij = indices->getNOCC( irrep_i );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ij += NOCC_i * NOCC_j;
               }
            }
            jump_aij += NVIR_a * jump_ij;
         }
      }
   }

   // FEE triplet: < TE_ukbl | f_pq E_pq | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_E[ irrep ];
      if ( SIZE > 0 ){
         int jump_aij = 0;
         for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            const int irrep_occ = Irreps::directProd( irrep_a, irrep );
            int jump_ij = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
               if ( irrep_i == irrep_j ){
                  const int NOCC_ij = indices->getNOCC( irrep_i );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = i+1; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     const double f_ii = fock->get( irrep_i, i, i );
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        const double f_jj = fock->get( irrep_j, j, j );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                           const double beta = shifted_prefactor + f_aa - f_ii - f_jj;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ij += NOCC_i * NOCC_j;
               }
            }
            jump_aij += NVIR_a * jump_ij;
         }
      }
   }
   
   // FGG singlet: < SG_cjdu | f_pq E_pq | SG_aibt > = 2 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_G[ irrep ];
      if ( SIZE > 0 ){
         int jump_abi = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int irrep_vir = Irreps::directProd( irrep_i, irrep );
            int jump_ab = 0;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
               if ( irrep_a == irrep_b ){
                  const int NVIR_ab = indices->getNVIRT( irrep_a );
                  const int N_OA_ab = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a );
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                        const double f_bb = fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               }
               if ( irrep_a < irrep_b ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  for ( int a = 0; a < NVIR_a; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ab += NVIR_a * NVIR_b;
               }
            }
            jump_abi += NOCC_i * jump_ab;
         }
      }
   }
   
   // FGG triplet: < TG_cjdu | f_pq E_pq | TG_aibt > = 6 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int SIZE = size_G[ irrep ];
      if ( SIZE > 0 ){
         int jump_abi = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int irrep_vir = Irreps::directProd( irrep_i, irrep );
            int jump_ab = 0;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
               if ( irrep_a == irrep_b ){
                  const int NVIR_ab = indices->getNVIRT( irrep_a );
                  const int N_OA_ab = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a );
                     for ( int b = a+1; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                        const double f_bb = fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
               }
               if ( irrep_a < irrep_b ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  for ( int a = 0; a < NVIR_a; a++ ){
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double f_ii = fock->get( irrep_i, i, i );
                           const double beta = shifted_prefactor + f_aa + f_bb - f_ii;
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                        }
                     }
                  }
                  jump_ab += NVIR_a * NVIR_b;
               }
            }
            jump_abi += NOCC_i * jump_ab;
         }
      }
   }

   // FHH singlet and triplet
   {
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
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
                        target[ cnt_ij + size_ij * cnt_ab ] = ((( a==b ) && ( i==j )) ? 8 : 4 ) * term;
                     }
                  }
               }
            }
            jump_aibj += size_ij * size_ab;
         }
      }
   }
   {
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
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
            jump_aibj += size_ij * size_ab;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int jump_aibj = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_j = Irreps::directProd( irrep, irrep_i );
         if ( irrep_i < irrep_j ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NOCC_j = indices->getNOCC( irrep_j );
            const int size_ij = NOCC_i * NOCC_j;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep, irrep_a );
               if ( irrep_a < irrep_b ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int size_ab = NVIR_a * NVIR_b;
                  double * target_singlet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
                  double * target_triplet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
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
                  jump_aibj += size_ij * size_ab;
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
           < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ] / sqrt( 2 - ( 1 - delta_ij ) * ( 1 - delta_ab ) )
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
      int jump_ai = 0;
      const int D2JUMP = size_D[ irrep ] / 2;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_a = Irreps::directProd( irrep_i, irrep );
         const int NOCC_i  = indices->getNOCC( irrep_i );
         const int N_OA_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
         const int NVIR_a  = indices->getNVIRT( irrep_a );
         for ( int count_i = 0; count_i < NOCC_i; count_i++ ){
            for ( int count_a = 0; count_a < NVIR_a; count_a++ ){

               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + size_D[ irrep ] * ( jump_ai + count_i + NOCC_i * count_a );
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
         jump_ai += NOCC_i * NVIR_a;
      }
      assert( jump_ai * size_D[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_D ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] );
   }
   delete MAT;

   // VB singlet and triplet
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      int jump_ij = 0; // First do SINGLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with [ (ix|jy) + (iy|jx) ] * (( x==y ) ? 0.5 : 1.0 )
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
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
               double alpha = (( i == j ) ? sqrt( 0.5 ) : 1.0 );
               double set = 0.0;
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ 0 ] * ( jump_ij + i + (j*(j+1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         jump_ij += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
      }
      assert( jump_ij * size_B_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );

      jump_ij = 0; // Then do TRIPLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i+1; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with [ (ix|jy) - (iy|jx) ]
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
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
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ 0 ] * ( jump_ij + i + (j*(j-1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         jump_ij += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
      }
      assert( jump_ij * size_B_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int jump_ij = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int irrep_j = Irreps::directProd( irrep, irrep_i );
         if ( irrep_i < irrep_j ){
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NOCC_j = indices->getNOCC( irrep_j );
            for ( int i = 0; i < NOCC_i; i++ ){
               for ( int j = 0; j < NOCC_j; j++ ){

                  // Fill workspace[ xy ] with [ (ix|jy) + (iy|jx) ]
                  int jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int d_x   = indices->getDMRGcumulative( irrep_x );
                        const int d_y   = indices->getDMRGcumulative( irrep_y );
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
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ irrep ] * ( jump_ij + i + NOCC_i * j );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_singlet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );

                  // Fill workspace[ xy ] with [ (ix|jy) - (iy|jx) ]
                  jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int d_x   = indices->getDMRGcumulative( irrep_x );
                        const int d_y   = indices->getDMRGcumulative( irrep_y );
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
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ irrep ] * ( jump_ij + i + NOCC_i * j );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
               }
            }
            jump_ij += NOCC_i * NOCC_j;
         }
      }
      assert( jump_ij * size_B_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );
      assert( jump_ij * size_B_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
   }

   // VF singlet and triplet
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      int jump_ab = 0; // First do SINGLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with [ (ax|by) + (ay|bx) ] * (( x==y ) ? 0.5 : 1.0 )
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
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
               double alpha = (( a == b ) ? sqrt( 0.5 ) : 1.0 );
               double set = 0.0;
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ 0 ] * ( jump_ab + a + (b*(b+1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
      }
      assert( jump_ab * size_F_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );

      jump_ab = 0; // Then do TRIPLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a+1; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with [ (ax|by) - (ay|bx) ]
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
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
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ 0 ] * ( jump_ab + a + (b*(b-1))/2 );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
         }
         jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
      }
      assert( jump_ab * size_F_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int jump_ab = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int irrep_b = Irreps::directProd( irrep, irrep_a );
         if ( irrep_a < irrep_b ){
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
                        const int d_x = indices->getDMRGcumulative( irrep_x );
                        const int d_y = indices->getDMRGcumulative( irrep_y );
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
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ irrep ] * ( jump_ab + a + NVIR_a * b );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_singlet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );

                  // Fill workspace[ xy ] with [ (ax|by) - (ay|bx) ]
                  jump_xy = 0;
                  for ( int irrep_x = 0; irrep_x < num_irreps; irrep_x++ ){
                     const int irrep_y = Irreps::directProd( irrep, irrep_x );
                     if ( irrep_x < irrep_y ){
                        const int d_x = indices->getDMRGcumulative( irrep_x );
                        const int d_y = indices->getDMRGcumulative( irrep_y );
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
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ irrep ] * ( jump_ab + a + NVIR_a * b );
                  dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
               }
            }
            jump_ab += NVIR_a * NVIR_b;
         }
      }
      assert( jump_ab * size_F_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );
      assert( jump_ab * size_F_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
   }
   delete [] workspace;

   // VE singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int occ_t = indices->getNOCC( irrep );
      const int num_t = indices->getNDMRG( irrep );
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ];
      int jump_aij_singlet = 0;
      int jump_aij_triplet = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int NVIR_a = indices->getNVIRT( irrep_a );
         const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
         const int irrep_occ = Irreps::directProd( irrep_a, irrep );
         if ( irrep_occ == 0 ){
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               const int NOCC_ij = indices->getNOCC( irrep_ij );
               for ( int i = 0; i < NOCC_ij; i++ ){
                  for ( int j = i; j < NOCC_ij; j++ ){
                     for ( int a = 0; a < NVIR_a; a++ ){
                        const int count_aij_singlet = jump_aij_singlet + a + NVIR_a * ( i + (j*(j+1))/2 );
                        const int count_aij_triplet = jump_aij_triplet + a + NVIR_a * ( i + (j*(j-1))/2 );
                        for ( int t = 0; t < num_t; t++ ){
                           double value_singlet = 0.0;
                           double value_triplet = 0.0;
                           for ( int w = 0; w < num_t; w++ ){
                              const double SEE_wt = SEE[ irrep ][ w + num_t * t ];
                              const double aj_wi  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, i, occ_t + w, j, N_OA_a + a );
                              const double ai_wj  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, j, occ_t + w, i, N_OA_a + a );
                              value_singlet +=     SEE_wt * ( aj_wi + ai_wj ) * (( i == j ) ? sqrt( 0.5 ) : 1.0 );
                              value_triplet += 3 * SEE_wt * ( aj_wi - ai_wj );
                           }
                           target_singlet[ t + num_t * count_aij_singlet ] = value_singlet;
             if ( j > i ){ target_triplet[ t + num_t * count_aij_triplet ] = value_triplet; }
                        }
                     }
                  }
               }
               jump_aij_singlet += ( NVIR_a * NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               jump_aij_triplet += ( NVIR_a * NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
            }
         } else {
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep_i, irrep_occ );
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const int count_aij_singlet = jump_aij_singlet + a + NVIR_a * ( i + NOCC_i * j );
                           const int count_aij_triplet = jump_aij_triplet + a + NVIR_a * ( i + NOCC_i * j );
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
                  jump_aij_singlet += NVIR_a * NOCC_i * NOCC_j;
                  jump_aij_triplet += NVIR_a * NOCC_i * NOCC_j;
               }
            }
         }
      }
      assert( jump_aij_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] );
      assert( jump_aij_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] );
   }

   // VG singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int occ_t = indices->getNOCC( irrep );
      const int num_t = indices->getNDMRG( irrep );
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ];
      int jump_abi_singlet = 0;
      int jump_abi_triplet = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int NOCC_i = indices->getNOCC( irrep_i );
         const int irrep_virt = Irreps::directProd( irrep_i, irrep );
         if ( irrep_virt == 0 ){
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               for ( int i = 0; i < NOCC_i; i++ ){
                  for ( int a = 0; a < NVIR_ab; a++ ){
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const int count_abi_singlet = jump_abi_singlet + i + NOCC_i * ( a + (b*(b+1))/2 );
                        const int count_abi_triplet = jump_abi_triplet + i + NOCC_i * ( a + (b*(b-1))/2 );
                        for ( int t = 0; t < num_t; t++ ){
                           double value_singlet = 0.0;
                           double value_triplet = 0.0;
                           for ( int u = 0; u < num_t; u++ ){
                              const double SGG_ut = SGG[ irrep ][ u + num_t * t ];
                              const double ai_bu  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + a, N_OA_ab + b );
                              const double bi_au  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + b, N_OA_ab + a );
                              value_singlet +=     SGG_ut * ( ai_bu + bi_au ) * (( a == b ) ? sqrt( 0.5 ) : 1.0 );
                              value_triplet += 3 * SGG_ut * ( ai_bu - bi_au );
                           }
                           target_singlet[ t + num_t * count_abi_singlet ] = value_singlet;
             if ( b > a ){ target_triplet[ t + num_t * count_abi_triplet ] = value_triplet; }
                        }
                     }
                  }
               }
               jump_abi_singlet += ( NOCC_i * NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               jump_abi_triplet += ( NOCC_i * NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            }
         } else {
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep_a, irrep_virt );
               if ( irrep_a < irrep_b ){
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int a = 0; a < NVIR_a; a++ ){
                        for ( int b = 0; b < NVIR_b; b++ ){
                           const int count_abi_singlet = jump_abi_singlet + i + NOCC_i * ( a + NVIR_a * b );
                           const int count_abi_triplet = jump_abi_triplet + i + NOCC_i * ( a + NVIR_a * b );
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
                  jump_abi_singlet += NOCC_i * NVIR_a * NVIR_b;
                  jump_abi_triplet += NOCC_i * NVIR_a * NVIR_b;
               }
            }
         }
      }
      assert( jump_abi_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] );
      assert( jump_abi_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] );
   }

   // VH singlet and triplet
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int jump_aibj_singlet = 0;
      int jump_aibj_triplet = 0;
      double * target_singlet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ];
      double * target_triplet = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ];
      if ( irrep == 0 ){ // irrep_i == irrep_j  and  irrep_a == irrep_b
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            const int linsize_singlet = ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
            const int linsize_triplet = ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int nvirt_ab = indices->getNVIRT( irrep_ab );
               const int noa_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               for ( int a = 0; a < nvirt_ab; a++ ){
                  for ( int b = a; b < nvirt_ab; b++ ){
                     for ( int i = 0; i < nocc_ij; i++ ){
                        for ( int j = i; j < nocc_ij; j++){
                           const int count_singlet = jump_aibj_singlet + i + (j*(j+1))/2 + linsize_singlet * ( a + (b*(b+1))/2 );
                           const int count_triplet = jump_aibj_triplet + i + (j*(j-1))/2 + linsize_triplet * ( a + (b*(b-1))/2 );
                           const double ai_bj = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, i, j, noa_ab + a, noa_ab + b );
                           const double aj_bi = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, j, i, noa_ab + a, noa_ab + b );
                           target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi ) * ((( i==j ) || ( a==b )) ? sqrt( 0.5 ) : 1.0 ); 
   if ( (b-a)*(j-i) > 0 ){ target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi ); }
                        }
                     }
                  }
               }
               jump_aibj_singlet += linsize_singlet * ( ( nvirt_ab * ( nvirt_ab + 1 ) ) / 2 );
               jump_aibj_triplet += linsize_triplet * ( ( nvirt_ab * ( nvirt_ab - 1 ) ) / 2 );
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
                     const int nvir_a = indices->getNVIRT( irrep_a );
                     const int nvir_b = indices->getNVIRT( irrep_b );
                     const int noa_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                     const int noa_b  = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                     for ( int a = 0; a < nvir_a; a++ ){
                        for ( int b = 0; b < nvir_b; b++ ){
                           for ( int i = 0; i < nocc_i; i++ ){
                              for ( int j = 0; j < nocc_j; j++){
                                 const int count_singlet = jump_aibj_singlet + i + nocc_i * ( j + nocc_j * ( a + nvir_a * b ) );
                                 const int count_triplet = jump_aibj_triplet + i + nocc_i * ( j + nocc_j * ( a + nvir_a * b ) );
                                 const double ai_bj = integrals->get_exchange( irrep_i, irrep_j, irrep_a, irrep_b, i, j, noa_a + a, noa_b + b );
                                 const double aj_bi = integrals->get_exchange( irrep_j, irrep_i, irrep_a, irrep_b, j, i, noa_a + a, noa_b + b );
                                 target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi );
                                 target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi );
                              }
                           }
                        }
                     }
                     jump_aibj_singlet += nocc_i * nocc_j * nvir_a * nvir_b;
                     jump_aibj_triplet += nocc_i * nocc_j * nvir_a * nvir_b;
                  }
               }
            }
         }
      }
      assert( jump_aibj_singlet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] );
      assert( jump_aibj_triplet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] );
   }

}

void CheMPS2::CASPT2::make_FEH_FGH(){

   /*
      FEH singlet: < SE_xkcl | F | SH_aibj > = +2 delta_ik delta_jl ( delta_bc FEH[ Ix ][ xa ] + delta_ac FEH[ Ix ][ xb ] )
      FEH triplet: < TE_xkcl | F | TH_aibj > = +6 delta_ik delta_jl ( delta_bc FEH[ Ix ][ xa ] - delta_ac FEH[ Ix ][ xb ] )
      FGH singlet: < SG_ckdx | F | SH_aibj > = -2 delta_ac delta_bd ( delta_ik FGH[ Ix ][ xj ] + delta_jk FGH[ Ix ][ xi ] )
      FGH triplet: < TG_ckdx | F | TH_aibj > = -6 delta_ac delta_bd ( delta_ik FGH[ Ix ][ xj ] - delta_jk FGH[ Ix ][ xi ] )

         FEH[ Ix ][ xc ] = sum_w SEE[ Ix ][ xw ] fock[ wc ]
         FGH[ Ix ][ xk ] = sum_w SGG[ Ix ][ xw ] fock[ wk ]

   */

   FEH = new double*[ num_irreps ];
   FGH = new double*[ num_irreps ];

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int NOCC = indices->getNOCC( irrep );
      const int NACT = indices->getNDMRG( irrep );
      const int N_OA = NOCC + NACT;
      const int NVIR = indices->getNVIRT( irrep );

      FEH[ irrep ] = new double[ NACT * NVIR ];
      FGH[ irrep ] = new double[ NACT * NOCC ];

      for ( int c = 0; c < NVIR; c++ ){
         for ( int x = 0; x < NACT; x++ ){
            double value = 0.0;
            for ( int w = 0; w < NACT; w++ ){ value += SEE[ irrep ][ x + NACT * w ] * fock->get( irrep, NOCC + w, N_OA + c ); }
            FEH[ irrep ][ x + NACT * c ] = value;
         }
      }

      for ( int k = 0; k < NOCC; k++ ){
         for ( int x = 0; x < NACT; x++ ){
            double value = 0.0;
            for ( int w = 0; w < NACT; w++ ){ value += SGG[ irrep ][ x + NACT * w ] * fock->get( irrep, NOCC + w, k ); }
            FEH[ irrep ][ x + NACT * k ] = value;
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



