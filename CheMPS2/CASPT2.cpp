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
#include <iostream>
#include <math.h>
#include <assert.h>
#include <algorithm>

#include "CASPT2.h"
#include "Lapack.h"

using std::cout;
using std::endl;
using std::max;

CheMPS2::CASPT2::CASPT2(DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei_in, DMRGSCFmatrix * fock_in, double * one_dm, double * two_dm, double * three_dm, double * contract){

   indices    = idx;
   integrals  = ints;
   oei        = oei_in;
   fock       = fock_in;
   one_rdm    = one_dm;
   two_rdm    = two_dm;
   three_rdm  = three_dm;
   contracted = contract;
   num_irreps = indices->getNirreps();
   E_FOCK     = calc_fock_expectation(); // Uses indices, num_irreps, fock, one_rdm
   
   vector_helper(); // Needs to be called BEFORE make_S**()!
   
   make_SAA_SCC();
   make_SDD();
   make_SEE_SGG();
   make_SBB_SFF_singlet();
   make_SBB_SFF_triplet();
   
   construct_rhs(); // Needs to be called AFTER make_S**()!

}

CheMPS2::CASPT2::~CASPT2(){

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
   
   delete [] size_AC;
   delete [] size_D;
   delete [] size_BF_singlet;
   delete [] size_BF_triplet;
   delete [] jump;
   delete [] vector_rhs;

}

double CheMPS2::CASPT2::calc_fock_expectation() const{

   double value = 0.0;
   int passed = 0;
   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
   
      const int NOCC  = indices->getNOCC(  irrep );
      const int NDMRG = indices->getNDMRG( irrep );
      
      for ( int orb = 0; orb < NOCC; orb++ ){ value += 2 * fock->get( irrep, orb, orb ); }
      for ( int row = 0; row < NDMRG; row++ ){
         for ( int col = 0; col < NDMRG; col++){
            value += one_rdm[ passed + row + LAS * ( passed + col ) ] * fock->get( irrep, row, col );
         }
      }
      
      passed += NDMRG;
      
   }
   
   cout << "CASPT2::calc_fock_expectation : E_FOCK = " << value << endl;
   
   return value;

}

int CheMPS2::CASPT2::vector_helper(){

   int * helper = new int[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];

   /*** Type A : c_tiuv E_ti E_uv | 0 >
                 c_tiuv = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + count_tuv + size_AC[ irrep ] * count_i ]
        Type C : c_atuv E_at E_uv | 0 >
                 c_atuv = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + count_tuv + size_AC[ irrep ] * count_a ]

        1/ (TYPE A) count_i = 0 .. NOCC[  irrep ]
           (TYPE C) count_a = 0 .. NVIRT[ irrep ]
        2/ jump_tuv = 0
           irrep_t = 0 .. num_irreps
              irrep_u = 0 .. num_irreps
                 irrep_v = irrep x irrep_t x irrep_u
                    ---> count_tuv = jump_tuv + t + NDMRG[ irrep_t ] * ( u + NDMRG[ irrep_u ] * v )
                    jump_tuv += NDMRG[ irrep_t ] * NDMRG[ irrep_u ] * NDMRG[ irrep_v ]
   */

   size_AC = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int linsize_AC = 0;
      for ( int irrep_t = 0; irrep_t < num_irreps; irrep_t++ ){
         for ( int irrep_u = 0; irrep_u < num_irreps; irrep_u++ ){
            const int irrep_v = Irreps::directProd( Irreps::directProd( irrep, irrep_t ), irrep_u );
            linsize_AC += indices->getNDMRG( irrep_t ) * indices->getNDMRG( irrep_u ) * indices->getNDMRG( irrep_v );
         }
      }
      size_AC[ irrep ] = linsize_AC;
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_A ] = indices->getNOCC(  irrep ) * size_AC[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_C ] = indices->getNVIRT( irrep ) * size_AC[ irrep ];
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
         const int nocc_i  = indices->getNOCC(  irrep_i );
         const int nvirt_a = indices->getNVIRT( irrep_a );
         jump_ai += nocc_i * nvirt_a;
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_D ] = jump_ai * size_D[ irrep ];
   }

   /*** Type B singlet : c_tiuj ( E_ti E_uj + E_tj E_ui ) | 0 > with i <= j and t <= u
                         c_tiuj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + count_tu + size_BF_singlet[ irrep ] * count_ij ]
        Type B triplet : c_tiuj ( E_ti E_uj - E_tj E_ui ) | 0 > with i <  j and t <  u
                         c_tiuj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + count_tu + size_BF_triplet[ irrep ] * count_ij ]
        Type F singlet : c_atbu ( E_at E_bu + E_bt E_au ) | 0 > with a <= b and t <= u
                         c_atbu = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + count_tu + size_BF_singlet[ irrep ] * count_ab ]
        Type F triplet : c_atbu ( E_at E_bu - E_bt E_au ) | 0 > with a <  b and t <  u
                         c_atbu = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + count_tu + size_BF_triplet[ irrep ] * count_ab ]

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
   
   
   size_BF_singlet = new int[ num_irreps ];
   size_BF_triplet = new int[ num_irreps ];
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int jump_tu_singlet = 0;
      int jump_tu_triplet = 0;
      if ( irrep == 0 ){
         for ( int irrep_tu = 0; irrep_tu < num_irreps; irrep_tu++ ){ // irrep_u == irrep_t
            const int nact_tu = indices->getNDMRG( irrep_tu );
            jump_tu_singlet += ( nact_tu * ( nact_tu + 1 ) ) / 2;
            jump_tu_triplet += ( nact_tu * ( nact_tu - 1 ) ) / 2;
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
      size_BF_singlet[ irrep ] = jump_tu_singlet;
      size_BF_triplet[ irrep ] = jump_tu_triplet;
   
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
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] = linsize_B_singlet * size_BF_singlet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] = linsize_B_triplet * size_BF_triplet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] = linsize_F_singlet * size_BF_singlet[ irrep ];
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] = linsize_F_triplet * size_BF_triplet[ irrep ];
   }
   
   /*** Type E singlet : c_tiaj ( E_ti E_aj + E_tj E_ai ) | 0 > with i <= j
                         c_tiaj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + count_t + NDMRG[ irrep ] * count_aij ]
        Type E triplet : c_tiaj ( E_ti E_aj - E_tj E_ai ) | 0 > with i <  j
                         c_tiaj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + count_t + NDMRG[ irrep ] * count_aij ]

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
   
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int linsize_E_singlet = 0;
      int linsize_E_triplet = 0;
      for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
         const int nvirt_a = indices->getNVIRT( irrep_a );
         const int irrep_occ = Irreps::directProd( irrep, irrep_a );
         if ( irrep_occ == 0 ){ // irrep_i == irrep_j
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               const int nocc_ij = indices->getNOCC( irrep_ij );
               linsize_E_singlet += ( nvirt_a * nocc_ij * ( nocc_ij + 1 ))/2;
               linsize_E_triplet += ( nvirt_a * nocc_ij * ( nocc_ij - 1 ))/2;
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
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] = linsize_E_singlet * indices->getNDMRG( irrep );
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] = linsize_E_triplet * indices->getNDMRG( irrep );
   }
   
   /*** Type G singlet : c_aibt ( E_ai E_bt + E_bi E_at ) | 0 > with a <= b
                         c_aibt = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + count_t + NDMRG[ irrep ] * count_abi ]
        Type G triplet : c_aibt ( E_ai E_bt - E_bi E_at ) | 0 > with a <  b
                         c_aibt = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + count_t + NDMRG[ irrep ] * count_abi ]

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
   
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int linsize_G_singlet = 0;
      int linsize_G_triplet = 0;
      for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
         const int nocc_i = indices->getNOCC( irrep_i );
         const int irrep_virt = Irreps::directProd( irrep, irrep_i );
         if ( irrep_virt == 0 ){ // irrep_a == irrep_b
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int nvirt_ab = indices->getNVIRT( irrep_ab );
               linsize_G_singlet += ( nocc_i * nvirt_ab * ( nvirt_ab + 1 ))/2;
               linsize_G_triplet += ( nocc_i * nvirt_ab * ( nvirt_ab - 1 ))/2;
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
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] = linsize_G_singlet * indices->getNDMRG( irrep );
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] = linsize_G_triplet * indices->getNDMRG( irrep );
   }
   
   /*** Type H singlet : c_aibj ( E_ai E_bj + E_bi E_aj ) | 0 > with a <= b and i <= j
                         c_aibj = vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + count_aibj ]
        Type H triplet : c_aibj ( E_ai E_bj - E_bi E_aj ) | 0 > with a <  b and i <  j
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
         linsize_H_triplet = linsize_H_singlet;
      }
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] = linsize_H_singlet;
      helper[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] = linsize_H_triplet;
   }
   
   jump = new int[ CHEMPS2_CASPT2_NUM_CASES * num_irreps + 1 ];
   jump[ 0 ] = 0;
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES * num_irreps; cnt++ ){ jump[ cnt+1 ] = jump[ cnt ] + helper[ cnt ]; }
   delete [] helper;
   
   const int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
   cout << "CASPT2::vector_helper : Total size of the V_SD space = " << total_size << endl;
   
   /*int check = 0;
   for ( int i1 = 0; i1 < num_irreps; i1++ ){
      for ( int i2 = 0; i2 < num_irreps; i2++ ){
         for ( int i3 = 0; i3 < num_irreps; i3++ ){
            const int i4 = Irreps::directProd( Irreps::directProd( i1, i2 ), i3 ); // V_SD: E_12 E_34 | 0 >
            check +=     indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNDMRG( i4 ); // A:  E_ti E_uv | 0 >
            check +=     indices->getNDMRG( i1 ) * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // C:  E_at E_uv | 0 >
            check += 2 * indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // D:  E_ai E_tu | 0 >  and  E_ti E_au | 0 >
            check +=     indices->getNOCC( i1 )  * indices->getNOCC( i2 )  * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // E:  E_ti E_aj | 0 >
            check +=     indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNVIRT( i3 ) * indices->getNVIRT( i4 ); // G:  E_ai E_bt | 0 >
            if ( i2 < i4 ){
               // 2 < 4 and irrep_2 < irrep_4
               check += indices->getNDMRG( i1 ) * indices->getNDMRG( i3 ) * indices->getNOCC( i2 )  * indices->getNOCC( i4 );  // B:  E_ti E_uj | 0 >
               check += indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNOCC( i2 )  * indices->getNOCC( i4 );  // H:  E_ai E_bj | 0 >
               check += indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNDMRG( i2 ) * indices->getNDMRG( i4 ); // F:  E_at E_bu | 0 >
            }
            if ( i2 == i4 ){ // i2 == i4 implies i1 == i3
               // 2 < 4 and irrep_2 == irrep_4
               check += ( indices->getNDMRG( i1 ) * indices->getNDMRG( i3 ) * indices->getNOCC( i2 )  * ( indices->getNOCC( i2 )  - 1 ) ) / 2; // B:  E_ti E_uj | 0 >
               check += ( indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNOCC( i2 )  * ( indices->getNOCC( i2 )  - 1 ) ) / 2; // H:  E_ai E_bj | 0 >
               check += ( indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNDMRG( i2 ) * ( indices->getNDMRG( i2 ) - 1 ) ) / 2; // F:  E_at E_bu | 0 >
               // 2 == 4 and 1 <= 3
               check += ( indices->getNDMRG( i1 ) * ( indices->getNDMRG( i3 ) + 1 ) * indices->getNOCC( i2 )  ) / 2; // B:  E_ti E_uj | 0 >
               check += ( indices->getNVIRT( i1 ) * ( indices->getNVIRT( i3 ) + 1 ) * indices->getNDMRG( i2 ) ) / 2; // F:  E_at E_bu | 0 >
               check += ( indices->getNVIRT( i1 ) * ( indices->getNVIRT( i3 ) + 1 ) * indices->getNOCC( i2 )  ) / 2; // H:  E_ai E_bj | 0 >
            }
         }
      }
   }
   assert( check == total_size );*/

   return total_size;

}

void CheMPS2::CASPT2::construct_rhs(){

   /*
      VA:  < H E_ti E_uv > = sum_w ( t_iw + sum_k [ 2 (iw|kk) - (ik|kw) ] ) [ 2 delta_tw Gamma_uv - Gamma_tuwv - delta_wu Gamma_tv ]
                           + sum_xzy (ix|zy) SAA[ It x Iu x Iv ][ xyztuv ]

      VB:  < H E_ti E_uj >
           < S_tiuj | H > = sum_xy (ix|jy) SBB_singlet[ It x Iu ][ xytu ]
           < T_tiuj | H > = sum_xy (ix|jy) SBB_triplet[ It x Iu ][ xytu ]

      VC:  < H E_at E_uv > = sum_w ( t_wa + sum_k [ 2 (wa|kk) - (wk|ka) ] ) < E_wt E_uv >
                           + sum_wxy (xy|wa) < E_xy E_wt E_uv >
           < H E_at E_uv > = sum_w ( t_wa + sum_k [ 2 (wa|kk) - (wk|ka) ] ) [ Gamma_wutv + delta_ut Gamma_wv ]
                           + sum_zxy (zy|xa) SCC[ It x Iu x Iv ][ xyztuv ]

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
           < S_tiaj | H > = sum_w [ (aj|wi) + (ai|wj) ] * 1 * SEE[ It ][ wt ]  -----------> OK
           < T_tiaj | H > = sum_w [ (aj|wi) - (ai|wj) ] * 3 * SEE[ It ][ wt ]  -----------> OK

      VF:  < H E_at E_bu >
           < S_atbu | H > = sum_xy (ax|by) SFF_singlet[ It x Iu ][ xytu ]
           < T_atbu | H > = sum_xy (ax|by) SFF_triplet[ It x Iu ][ xytu ]

      VG:  < H E_ai E_bt >
           < S_aibt | H > = sum_u [ (ai|bu) + (bi|au) ] * 1 * SGG[ It ][ ut ]  -----------> OK
           < T_aibt | H > = sum_u [ (ai|bu) - (bi|au) ] * 3 * SGG[ It ][ ut ]  -----------> OK

      VH:  < H E_ai E_bj >
           < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ]  -----------> OK
           < T_aibj | H > = 6 * [ (ai|bj) - (aj|bi) ]  -----------> OK
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

   vector_rhs = new double[ jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ] ];

   int max_size = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      max_size = max( max( max( max( size_AC[ irrep ], size_D[ irrep ] ), size_BF_singlet[ irrep ] ), size_BF_triplet[ irrep ] ), max_size );
   }
   double * workspace = new double[ max_size ];

   // VA: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + count_tuv + size_AC[ irrep ] * count_i ]
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      const int NACT = indices->getNDMRG( irrep );
      const int d_w  = indices->getDMRGcumulative( irrep );
      for ( int count_i = 0; count_i < NOCC; count_i++ ){

         double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + size_AC[ irrep ] * count_i;

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
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ]
                           = integrals->get_coulomb( irrep, irrep_x, irrep_z, irrep_y, count_i, occ_x + x, occ_z + z, occ_y + y );
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
         assert( jump_xyz == size_AC[ irrep ] );

         // Perform target[ tuv ] += sum_xzy (ix|zy) SAA[ It x Iu x Iv ][ xyztuv ]
         char notrans = 'N';
         int int1 = 1;
         double alpha = 1.0;
         double beta = 1.0; //ADD
         dgemm_( &notrans, &notrans, &int1, &jump_xyz, &jump_xyz, &alpha, workspace, &int1, SAA[ irrep ], &jump_xyz, &beta, target, &int1 );
      }
      assert( NOCC * size_AC[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_A ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] );
   }
   
   // VC: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + count_tuv + size_AC[ irrep ] * count_a ]
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC  = indices->getNOCC( irrep );
      const int NVIR  = indices->getNVIRT( irrep );
      const int NACT  = indices->getNDMRG( irrep );
      const int N_OA  = NOCC + NACT;
      const int d_w   = indices->getDMRGcumulative( irrep );
      for ( int count_a = 0; count_a < NVIR; count_a++ ){

         double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + size_AC[ irrep ] * count_a;

         // Fill workspace[ xyz ] with (zy|xa)
         // Fill target[ tuv ] with sum_w MAT[w,a] [ Gamma_wutv + delta_ut Gamma_wv ]
         //                       = sum_w MAT[w,a] Gamma_wutv + sum_w MAT[w,a] delta_ut Gamma_wv
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
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ]
                           = integrals->get_coulomb( irrep_z, irrep_y, irrep_x, irrep, occ_z + z, occ_y + y, occ_x + x, N_OA + count_a );
                     }
                  }
               }

               // target[ tuv ] = sum_w MAT[w,a] Gamma_wutv
               for ( int z = 0; z < num_z; z++ ){
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        double value = 0.0;
                        for ( int w = 0; w < NACT; w++ ){
                           value += MAT->get( irrep, NOCC + w, N_OA + count_a ) * two_rdm[ d_w + w + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_z + z ) ) ) ];
                        }
                        target[ jump_xyz + x + num_x * ( y + num_y * z ) ] = value;
                     }
                  }
               }

               // target[ tuv ] += sum_w MAT[w,a] delta_ut Gamma_wv
               if (( irrep_z == irrep ) && ( irrep_x == irrep_y )){
                  for ( int z = 0; z < num_z; z++ ){ // v
                     double value = 0.0;
                     for ( int w = 0; w < NACT; w++ ){
                        value += MAT->get( irrep, NOCC + w, N_OA + count_a ) * one_rdm[ d_w + w + LAS * ( d_z + z ) ];
                     }
                     for ( int xy = 0; xy < num_x; xy++ ){ // tu
                        target[ jump_xyz + xy + num_x * ( xy + num_y * z ) ] += value;
                     }
                  }
               }

               jump_xyz += num_x * num_y * num_z;
            }
         }
         assert( jump_xyz == size_AC[ irrep ] );

         // Perform target[ tuv ] += sum_zxy (zy|xa) SCC[ It x Iu x Iv ][ xyztuv ]
         char notrans = 'N';
         int int1 = 1;
         double alpha = 1.0;
         double beta = 1.0; //ADD
         dgemm_( &notrans, &notrans, &int1, &jump_xyz, &jump_xyz, &alpha, workspace, &int1, SCC[ irrep ], &jump_xyz, &beta, target, &int1 );
      }
      assert( NVIR * size_AC[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_C ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] );
   }

   // VD1: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + count_tu +          size_D[ irrep ] * count_ai ]
   // VD2: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + count_tu + D2JUMP + size_D[ irrep ] * count_ai ]
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
               const double MAT_ia = ( ( irrep == 0 ) ? MAT->get( irrep, count_i, N_OA_a + count_a ) : 0.0 );

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
                  const int occ_x   = indices->getNOCC( irrep_x );
                  const int occ_y   = indices->getNOCC( irrep_y );
                  const int num_x   = indices->getNDMRG( irrep_x );
                  const int num_y   = indices->getNDMRG( irrep_y );

                  // workspace[          xy ] = (ia|yx)
                  // workspace[ D2JUMP + xy ] = (ix|ya)
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        workspace[          jump_xy + x + num_x * y ] =
                           integrals->get_coulomb( irrep_y, irrep_x, irrep_i, irrep_a, occ_y + y, occ_x + x, count_i, N_OA_a + count_a ); // (ia|yx)
                        workspace[ D2JUMP + jump_xy + x + num_x * y ] =
                           integrals->get_coulomb( irrep_i, irrep_x, irrep_y, irrep_a, count_i, occ_x + x, occ_y + y, N_OA_a + count_a ); // (ix|ya)
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
               int int1 = 1;
               double alpha = 1.0;
               double beta = 1.0; //ADD
               dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SDD[ irrep ], &jump_xy, &beta, target, &int1 );

            }
         }
         jump_ai += NOCC_i * NVIR_a;
      }
      assert( jump_ai * size_D[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_D ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] );
   }

   // VB_singlet: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + count_tu + size_BF_singlet[ irrep ] * count_ij ]
   // VB_triplet: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + count_tu + size_BF_triplet[ irrep ] * count_ij ]
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      const int irrep = 0;

      int jump_ij = 0; // First do SINGLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with (ix|jy)
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                        workspace[ jump_xy + x + (y*(y+1))/2 ] =
                           integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y ); // (ix|jy)
                     }
                  }

                  jump_xy += ( num_xy * ( num_xy + 1 ) ) / 2;
               }
               assert( jump_xy == size_BF_singlet[ irrep ] );

               // Perform target[ tu ] = sum_xy (ix|jy) SBB_singlet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int int1 = 1;
               double alpha = 1.0;
               double beta = 0.0; //SET
               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_BF_singlet[ irrep ] * ( jump_ij + i + (j*(j+1))/2 );
               dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SBB_singlet[ irrep ], &jump_xy, &beta, target, &int1 );

            }
         }
         jump_ij += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
      }
      assert( jump_ij * size_BF_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );

      jump_ij = 0; // Then do TRIPLET
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         for ( int i = 0; i < NOCC_ij; i++ ){
            for ( int j = i+1; j < NOCC_ij; j++ ){

               // Fill workspace[ xy ] with (ix|jy)
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                        workspace[ jump_xy + x + (y*(y-1))/2 ] =
                           integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y ); // (ix|jy)
                     }
                  }

                  jump_xy += ( num_xy * ( num_xy - 1 ) ) / 2;
               }
               assert( jump_xy == size_BF_triplet[ irrep ] );

               // Perform target[ tu ] = sum_xy (ix|jy) SBB_triplet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int int1 = 1;
               double alpha = 1.0;
               double beta = 0.0; //SET
               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_BF_triplet[ irrep ] * ( jump_ij + i + (j*(j-1))/2 );
               dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SBB_triplet[ irrep ], &jump_xy, &beta, target, &int1 );

            }
         }
         jump_ij += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
      }
      assert( jump_ij * size_BF_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
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

                  // Fill workspace[ xy ] with (ix|jy)
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
                              workspace[ jump_xy + x + num_x * y ] =
                                 integrals->get_coulomb( irrep_i, irrep_x, irrep_j, irrep_y, i, occ_x + x, j, occ_y + y ); // (ix|jy)
                           }
                        }

                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_BF_singlet[ irrep ] );
                  assert( jump_xy == size_BF_triplet[ irrep ] );

                  // Perform target[ tu ] = sum_xy (ix|jy) SBB_singlet[ It x Iu ][ xytu ]
                  char notrans = 'N';
                  int int1 = 1;
                  double alpha = 1.0;
                  double beta = 0.0; //SET
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_BF_singlet[ irrep ] * ( jump_ij + i + NOCC_i * j );
                  dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SBB_singlet[ irrep ], &jump_xy, &beta, target, &int1 );

                  // Perform target[ tu ] = sum_xy (ix|jy) SBB_triplet[ It x Iu ][ xytu ]
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_BF_triplet[ irrep ] * ( jump_ij + i + NOCC_i * j );
                  dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SBB_triplet[ irrep ], &jump_xy, &beta, target, &int1 );

               }
            }
            jump_ij += NOCC_i * NOCC_j;
         }
      }
      assert( jump_ij * size_BF_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );
      assert( jump_ij * size_BF_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
   }

   // VF_singlet: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + count_tu + size_BF_singlet[ irrep ] * count_ab ]
   // VF_triplet: vector[ jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + count_tu + size_BF_triplet[ irrep ] * count_ab ]
   { // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
      const int irrep = 0;

      int jump_ab = 0; // First do SINGLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with (ax|by)
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                        workspace[ jump_xy + x + (y*(y+1))/2 ] =
                           integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b ); // (ax|by)
                     }
                  }

                  jump_xy += ( num_xy * ( num_xy + 1 ) ) / 2;
               }
               assert( jump_xy == size_BF_singlet[ irrep ] );

               // Perform target[ tu ] = sum_xy (ax|by) SFF_singlet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int int1 = 1;
               double alpha = 1.0;
               double beta = 0.0; //SET
               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_BF_singlet[ irrep ] * ( jump_ab + a + (b*(b+1))/2 );
               dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SFF_singlet[ irrep ], &jump_xy, &beta, target, &int1 );

            }
         }
         jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
      }
      assert( jump_ab * size_BF_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );

      jump_ab = 0; // Then do TRIPLET
      for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
         const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
         const int NVIR_ab = indices->getNVIRT( irrep_ab );
         for ( int a = 0; a < NVIR_ab; a++ ){
            for ( int b = a+1; b < NVIR_ab; b++ ){

               // Fill workspace[ xy ] with (ax|by)
               int jump_xy = 0;
               for ( int irrep_xy = 0; irrep_xy < num_irreps; irrep_xy++ ){
                  const int d_xy   = indices->getDMRGcumulative( irrep_xy );
                  const int occ_xy = indices->getNOCC( irrep_xy );
                  const int num_xy = indices->getNDMRG( irrep_xy );

                  for ( int x = 0; x < num_xy; x++ ){
                     for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                        workspace[ jump_xy + x + (y*(y-1))/2 ] =
                           integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b ); // (ax|by)
                     }
                  }

                  jump_xy += ( num_xy * ( num_xy - 1 ) ) / 2;
               }
               assert( jump_xy == size_BF_triplet[ irrep ] );

               // Perform target[ tu ] = sum_xy (ax|by) SFF_triplet[ It x Iu ][ xytu ]
               char notrans = 'N';
               int int1 = 1;
               double alpha = 1.0;
               double beta = 0.0; //SET
               double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_BF_triplet[ irrep ] * ( jump_ab + a + (b*(b-1))/2 );
               dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SFF_triplet[ irrep ], &jump_xy, &beta, target, &int1 );

            }
         }
         jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
      }
      assert( jump_ab * size_BF_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
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

                  // Fill workspace[ xy ] with (ax|by)
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
                              workspace[ jump_xy + x + num_x * y ] =
                                 integrals->get_exchange( irrep_x, irrep_y, irrep_a, irrep_b, occ_x + x, occ_y + y, N_OA_a + a, N_OA_b + b ); // (ax|by)
                           }
                        }

                        jump_xy += num_x * num_y;
                     }
                  }
                  assert( jump_xy == size_BF_singlet[ irrep ] );
                  assert( jump_xy == size_BF_triplet[ irrep ] );

                  // Perform target[ tu ] = sum_xy (ix|jy) SFF_singlet[ It x Iu ][ xytu ]
                  char notrans = 'N';
                  int int1 = 1;
                  double alpha = 1.0;
                  double beta = 0.0; //SET
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_BF_singlet[ irrep ] * ( jump_ab + a + NVIR_a * b );
                  dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SFF_singlet[ irrep ], &jump_xy, &beta, target, &int1 );

                  // Perform target[ tu ] = sum_xy (ix|jy) SFF_triplet[ It x Iu ][ xytu ]
                  target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_BF_triplet[ irrep ] * ( jump_ab + a + NVIR_a * b );
                  dgemm_( &notrans, &notrans, &int1, &jump_xy, &jump_xy, &alpha, workspace, &int1, SFF_triplet[ irrep ], &jump_xy, &beta, target, &int1 );

               }
            }
            jump_ab += NVIR_a * NVIR_b;
         }
      }
      assert( jump_ab * size_BF_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );
      assert( jump_ab * size_BF_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
   }

   delete [] workspace;
   delete MAT;

}

void CheMPS2::CASPT2::make_SAA_SCC(){

   /*
      SAA: < E_zy E_jx E_ti E_uv >
           = delta_ji ( 2 delta_tx Gamma_{zuyv}
                      + 2 delta_tx delta_uy Gamma_{zv}
                      - Gamma_{ztuyxv}
                      - delta_uy Gamma_{tzxv}
                      - delta_ty Gamma_{zuxv}
                      - delta_ux Gamma_{ztyv}
                      - delta_ux delta_ty Gamma_{zv} )

           with It x Ii x Iu x Iv = Iz x Iy x Ij x Ix = Itrivial

      SCC: < E_zy E_xb E_at E_uv >
           = delta_ab ( Gamma_{zxuytv}
                      + delta_uy Gamma_{xztv}
                      + delta_xy Gamma_{zutv}
                      + delta_ut Gamma_{zxyv}
                      + delta_ut delta_xy Gamma_{zv} )

           with It x Ii x Iu x Iv = Iz x Iy x Ij x Ix = Itrivial
   */
   
   SAA = new double*[ num_irreps ];
   SCC = new double*[ num_irreps ];
   
   const int LAS = indices->getDMRGcumulative( num_irreps );
   
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int SIZE = size_AC[ irrep ];
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
                           for ( int v = 0; v < num_v; v++ ){
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

void CheMPS2::CASPT2::make_SDD(){

   /*
      SD1D1: < E_yx E_jb E_ai E_tu >
             = delta_ab delta_ji ( 2 * Gamma_{ytxu}
                                 + 2 * delta_tx Gamma_{yu} )
                                 
      SD2D2: < E_yb E_jx E_ti E_au >
             = delta_ab delta_ji ( 2 * delta_tx Gamma_{yu}
                                 - Gamma_{ytux} )
                                 
      SD1D2: < E_yx E_jb E_ti E_au >
             = delta_ab delta_ji ( - Gamma_{ytxu}
                                   - delta_tx Gamma_{yu} )
                                 
      SD2D1: < E_yb E_jx E_ai E_tu >
             = delta_ab delta_ji ( - Gamma_{ytxu}
                                   - delta_tx Gamma_{yu} )
                      
      with Iu x It x Ii x Ia = Iy x Ix x Ij x Ib = Itrivial
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

void CheMPS2::CASPT2::make_SBB_SFF_singlet(){

   /*
      SBB singlet: | S_tiuj > = ( E_ti E_uj + E_tj E_ui ) | 0 >  with  i <= j and t <= u

                   < S_xkyl | S_tiuj > = 2 ( delta_ki delta_lj + delta_kj delta_li ) ( Gamma_{utyx}
                                                                                     + Gamma_{utxy}
                                                                                     + 2 delta_uy delta_tx
                                                                                     + 2 delta_ux delta_ty
                                                                                     -   delta_uy Gamma_{tx}
                                                                                     -   delta_tx Gamma_{uy}
                                                                                     -   delta_ux Gamma_{ty}
                                                                                     -   delta_ty Gamma_{ux} )
                                                                                     
      SFF singlet: | S_atbu > = ( E_at E_bu + E_bt E_au ) | 0 >  with  a <= b and t <= u

                   < S_cxdy | S_atbu > = 2 ( delta_ac delta_bd + delta_ad delta_bc ) ( Gamma_{yxut} + Gamma_{yxtu} )
      
      with Iia x Ijb x It x Iu = Ikc x Ild x Ix x Iy = Itrivial
      or hence Iia x Ijb = Ikc x Ild = It x Iu = Ix x Iy
      therefore:
         - if ( Iia x Ijb == 0 ) : It == Iu  and  Ix == Iy
         - if ( Iia x Ijb != 0 ) : It != Iu  and  Ix != Iy
   */

   SBB_singlet = new double*[ num_irreps ];
   SFF_singlet = new double*[ num_irreps ];
   
   const int LAS = indices->getDMRGcumulative( num_irreps );
   
   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      const int irrep = 0;

      const int SIZE = size_BF_singlet[ irrep ];
      SBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
      SFF_singlet[ irrep ] = new double[ SIZE * SIZE ];
      
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
                        SBB_singlet[ irrep ][ pointer ] = value;
                        SFF_singlet[ irrep ][ pointer ] = value;
                     }
                  }
               }
            }
            
            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy
            
               // + 2 ( delta_uy delta_tx + delta_ux delta_ty ) --> last term only if x <= y = t <= u or hence all equal
               for ( int t = 0; t < num_ut; t++ ){
                  SBB_singlet[ irrep ][ jump_row + t + (t*(t+1))/2 + SIZE * ( jump_col + t + (t*(t+1))/2 ) ] += 4.0;
                  for ( int u = t+1; u < num_ut; u++ ){
                     SBB_singlet[ irrep ][ jump_row + t + (u*(u+1))/2 + SIZE * ( jump_col + t + (u*(u+1))/2 ) ] += 2.0;
                  }
               }
               
               // - delta_uy Gamma_{tx}
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                     for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                        const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_xy + x ) ];
                        SBB_singlet[ irrep ][ jump_row + x + (uy*(uy+1))/2 + SIZE * ( jump_col + t + (uy*(uy+1))/2 ) ] -= gamma_tx;
                     }
                  }
               }
               
               // - delta_tx Gamma_{uy}
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                     for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                        const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_xy + y ) ];
                        SBB_singlet[ irrep ][ jump_row + tx + (y*(y+1))/2 + SIZE * ( jump_col + tx + (u*(u+1))/2 ) ] -= gamma_uy;
                     }
                  }
               }
               
               // - delta_ux Gamma_{ty}
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                     for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                        const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_xy + y ) ];
                        SBB_singlet[ irrep ][ jump_row + ux + (y*(y+1))/2 + SIZE * ( jump_col + t + (ux*(ux+1))/2 ) ] -= gamma_ty;
                     }
                  }
               }
               
               // - delta_ty Gamma_{ux}
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                     for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                        const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_xy + x ) ];
                        SBB_singlet[ irrep ][ jump_row + x + (ty*(ty+1))/2 + SIZE * ( jump_col + ty + (u*(u+1))/2 ) ] -= gamma_ux;
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

      const int SIZE = size_BF_singlet[ irrep ];
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

void CheMPS2::CASPT2::make_SBB_SFF_triplet(){

   /*
      SBB triplet: | T_tiuj > = ( E_ti E_uj + E_tj E_ui ) | 0 >  with  i < j and t < u

                   < T_xkyl | T_tiuj > = 2 ( delta_ki delta_lj - delta_kj delta_li ) ( Gamma_{utyx}
                                                                                     - Gamma_{utxy}
                                                                                     + 6 delta_uy delta_tx
                                                                                     - 6 delta_ux delta_ty
                                                                                     - 3 delta_uy Gamma_{tx}
                                                                                     - 3 delta_tx Gamma_{uy}
                                                                                     + 3 delta_ux Gamma_{ty}
                                                                                     + 3 delta_ty Gamma_{ux} )
                                                                                     
      SFF triplet: | T_atbu > = ( E_at E_bu - E_bt E_au ) | 0 >  with  a <= b and t <= u

                   < T_cxdy | T_atbu > = 2 ( delta_ac delta_bd - delta_ad delta_bc ) ( Gamma_{yxut} - Gamma_{yxtu} )
      
      with Iia x Ijb x It x Iu = Ikc x Ild x Ix x Iy = Itrivial
      or hence Iia x Ijb = Ikc x Ild = It x Iu = Ix x Iy
      therefore:
         - if ( Iia x Ijb == 0 ) : It == Iu  and  Ix == Iy
         - if ( Iia x Ijb != 0 ) : It != Iu  and  Ix != Iy
   */
   
   SBB_triplet = new double*[ num_irreps ];
   SFF_triplet = new double*[ num_irreps ];
   
   const int LAS = indices->getDMRGcumulative( num_irreps );
   
   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      const int irrep = 0;

      const int SIZE = size_BF_triplet[ irrep ];
      SBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
      SFF_triplet[ irrep ] = new double[ SIZE * SIZE ];
      
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
                        const double gamma_utyx = two_rdm[ d_ut + u + LAS * ( d_ut + t + LAS * ( d_xy + y + LAS * ( d_xy + x ))) ];
                        const double gamma_utxy = two_rdm[ d_ut + u + LAS * ( d_ut + t + LAS * ( d_xy + x + LAS * ( d_xy + y ))) ];
                        const double value = ( two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                             - two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                        const int pointer = jump_row + x + (y*(y-1))/2 + SIZE * ( jump_col + t + (u*(u-1))/2 );
                        SBB_triplet[ irrep ][ pointer ] = value;
                        SFF_triplet[ irrep ][ pointer ] = value;
                     }
                  }
               }
            }
            
            if ( irrep_ut == irrep_xy ){ // All four irreps are equal --> num_ut == num_xy

               // + 6 ( delta_uy delta_tx - delta_ux delta_ty ) --> last term only if x < y = t < u or NEVER
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int uy = tx+1; uy < num_ut; uy++ ){
                     SBB_triplet[ irrep ][ jump_row + tx + (uy*(uy-1))/2 + SIZE * ( jump_col + tx + (uy*(uy-1))/2 ) ] += 6.0;
                  }
               }
               
               // - 3 delta_uy Gamma_{tx}
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                     for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                        const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_xy + x ) ];
                        SBB_triplet[ irrep ][ jump_row + x + (uy*(uy-1))/2 + SIZE * ( jump_col + t + (uy*(uy-1))/2 ) ] -= 3 * gamma_tx;
                     }
                  }
               }
               
               // - 3 delta_tx Gamma_{uy}
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                     for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                        const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_xy + y ) ];
                        SBB_triplet[ irrep ][ jump_row + tx + (y*(y-1))/2 + SIZE * ( jump_col + tx + (u*(u-1))/2 ) ] -= 3 * gamma_uy;
                     }
                  }
               }
               
               // + 3 delta_ux Gamma_{ty}
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                     for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                        const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_xy + y ) ];
                        SBB_triplet[ irrep ][ jump_row + ux + (y*(y-1))/2 + SIZE * ( jump_col + t + (ux*(ux-1))/2 ) ] += 3 * gamma_ty;
                     }
                  }
               }
               
               // + 3 delta_ty Gamma_{ux}
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                     for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                        const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_xy + x ) ];
                        SBB_triplet[ irrep ][ jump_row + x + (ty*(ty-1))/2 + SIZE * ( jump_col + ty + (u*(u-1))/2 ) ] += 3 * gamma_ux;
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

      // Fill the triplet overlap matrix
      {
         const int SIZE = size_BF_triplet[ irrep ];
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

                        // 6 delta_uy delta_tx
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

}

void CheMPS2::CASPT2::make_SEE_SGG(){

   /*
      SEE singlet: | S_tiaj > = ( E_ti E_aj + E_tj E_ai ) | 0 >  with  i <= j
      SEE triplet: | T_tiaj > = ( E_ti E_aj - E_tj E_ai ) | 0 >  with  i <  j

                   < S_ukbl | S_tiaj > = 2 delta_ab ( 2 delta_tu - Gamma_tu ) ( delta_ki delta_lj + delta_kj delta_li )
                   < T_ukbl | T_tiaj > = 6 delta_ab ( 2 delta_tu - Gamma_tu ) ( delta_ki delta_lj - delta_kj delta_li )
                   
                   with It x Ii x Ia x Ij = Iu x Ik x Ib x Il = Itrivial
                   
      SGG singlet: | S_aibt > = ( E_ai E_bt + E_bi E_at ) | 0 >  with  a <= b
      SGG triplet: | T_aibt > = ( E_ai E_bt - E_bi E_at ) | 0 >  with  a <  b

                   < S_cjdu | S_aibt > = 2 delta_ji Gamma_ut ( delta_ca delta_db + delta_cb delta_da )
                   < T_cjdu | T_aibt > = 6 delta_ji Gamma_ut ( delta_ca delta_db - delta_cb delta_da )
      
                   with Ia x Ii x Ib x It = Ic x Ij x Id x Iu = Itrivial
   */
   
   SEE = new double*[ num_irreps ];
   SGG = new double*[ num_irreps ];
   
   const int LAS = indices->getDMRGcumulative( num_irreps );
   
   for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){

      const int SIZE = indices->getNDMRG( irrep_ut );
      SEE[ irrep_ut ] = new double[ SIZE * SIZE ];
      SGG[ irrep_ut ] = new double[ SIZE * SIZE ];
      const int d_ut   = indices->getDMRGcumulative( irrep_ut );
      const int num_ut = indices->getNDMRG( irrep_ut );
      
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



