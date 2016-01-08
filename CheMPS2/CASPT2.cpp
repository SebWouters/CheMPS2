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

using std::cout;
using std::endl;
using std::max;

CheMPS2::CASPT2::CASPT2(DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei_in, DMRGSCFmatrix * fock_in, double * one_dm, double * two_dm, double * three_dm, double * contract_4dm){

   indices    = idx;
   oei        = oei_in;
   fock       = fock_in;
   one_rdm    = one_dm;
   two_rdm    = two_dm;
   three_rdm  = three_dm;
   f_dot_4dm  = contract_4dm;
   num_irreps = indices->getNirreps();
   E_FOCK     = create_f_dots();
   
   vector_helper(); // Needs to be called BEFORE make_S**()!
   
   make_SAA_SCC();
   make_SDD();
   make_SEE_SGG();
   make_SBB_SFF_singlet();
   make_SBB_SFF_triplet();
   
   construct_rhs( ints ); // Needs to be called AFTER make_S**()!
   
   make_FAA_FCC();
   make_FDD();
   make_FEE_FGG();
   make_FBB_FFF_singlet();
   make_FBB_FFF_triplet();

   {
      
      int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
      
      double * dummy  = new double[ total_size ];
      double * dummy2 = new double[ total_size ];
      
      /*cout << "Start A         : " << jump[ num_irreps * CHEMPS2_CASPT2_A ] << endl;
      cout << "Start B singlet : " << jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] << endl;
      cout << "Start B triplet : " << jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] << endl;
      cout << "Start C         : " << jump[ num_irreps * CHEMPS2_CASPT2_C ] << endl;
      cout << "Start D         : " << jump[ num_irreps * CHEMPS2_CASPT2_D ] << endl;
      cout << "Start E singlet : " << jump[ num_irreps * CHEMPS2_CASPT2_E_SINGLET ] << endl;
      cout << "Start E triplet : " << jump[ num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] << endl;
      cout << "Start F singlet : " << jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] << endl;
      cout << "Start F triplet : " << jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] << endl;
      cout << "Start G singlet : " << jump[ num_irreps * CHEMPS2_CASPT2_G_SINGLET ] << endl;
      cout << "Start G triplet : " << jump[ num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] << endl;
      cout << "Start H singlet : " << jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] << endl;
      cout << "Start H triplet : " << jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] << endl;*/
      
      // Calculate P_SD [ blockdiag(F) - E_FOCK * S ]^{-1} P_SD H | Psi0 >
      inverse_blockdiag( vector_rhs, dummy, 1.0, -E_FOCK );
      
      // Calculate E(CASPT2-D) = - < Psi0 | H P_SD [ blockdiag(F) - E_FOCK * S ]^{-1} P_SD H | Psi0 >
      int inc = 1;
      const double energy_caspt2d = - ddot_( &total_size, dummy, &inc, vector_rhs, &inc );
      cout << "E(CASPT2-D) = " << energy_caspt2d << endl;
      
      // Calculate P_SD [ blockdiag(F) - E_FOCK * S ] P_SD [ blockdiag(F) - E_FOCK * S ]^{-1} P_SD H | Psi0 > which should equal P_SD H | Psi0 >
      matvec_blockdiag( dummy, dummy2, 1.0, -E_FOCK );
      
      // Test whether dummy2 == vector_rhs
      double rms = 0.0;
      for ( int cnt = 0; cnt < total_size; cnt++ ){
         const double diff = dummy2[ cnt ] - vector_rhs[ cnt ];
         rms += diff * diff;
      }
      rms = sqrt( rms );
      cout << "Norm of P_SD ( [ blockdiag(F) - E * S ] P_SD [ blockdiag(F) - E * S ]^{-1} P_SD - I ) H | Psi0 > = " << rms << endl;

      delete [] dummy;
      delete [] dummy2;
   }

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

   delete [] SAA;
   delete [] SCC;
   delete [] SDD;
   delete [] SEE;
   delete [] SGG;
   delete [] SBB_singlet;
   delete [] SBB_triplet;
   delete [] SFF_singlet;
   delete [] SFF_triplet;
   delete [] FAA;
   delete [] FCC;
   delete [] FDD;
   delete [] FEE;
   delete [] FGG;
   delete [] FBB_singlet;
   delete [] FBB_triplet;
   delete [] FFF_singlet;
   delete [] FFF_triplet;
   
   delete [] size_AC;
   delete [] size_D;
   delete [] size_BF_singlet;
   delete [] size_BF_triplet;
   delete [] jump;
   delete [] vector_rhs;
   delete [] f_dot_3dm;
   delete [] f_dot_2dm;

}

double CheMPS2::CASPT2::create_f_dots(){

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

      const int NDMRG1 = indices->getNDMRG( irrep1 );
      const int jump1  = indices->getDMRGcumulative( irrep1 );

      for ( int row1 = 0; row1 < NDMRG1; row1++ ){
         for ( int col1 = 0; col1 < NDMRG1; col1++ ){

            double value = 0.0;
            for ( int irrep2 = 0; irrep2 < num_irreps; irrep2++ ){

               const int NOCC2  = indices->getNOCC(  irrep2 );
               const int NDMRG2 = indices->getNDMRG( irrep2 );
               const int jump2  = indices->getDMRGcumulative( irrep2 );

               for ( int row2 = 0; row2 < NDMRG2; row2++ ){
                  for ( int col2 = 0; col2 < NDMRG2; col2++ ){
                     value += ( fock->get( irrep2, NOCC2 + row2, NOCC2 + col2 )
                              * two_rdm[ jump1 + row1 + LAS * ( jump2 + row2 + LAS * ( jump1 + col1 + LAS * ( jump2 + col2 ) ) ) ] );
                  }
               }
            }
            f_dot_2dm[ jump1 + row1 + LAS * ( jump1 + col1 ) ] = value;
         }
      }
   }

   for ( int irrep1 = 0; irrep1 < num_irreps; irrep1++ ){
      const int jump1  = indices->getDMRGcumulative( irrep1 );
      const int NDMRG1 = indices->getNDMRG( irrep1 );
      for ( int irrep2 = 0; irrep2 < num_irreps; irrep2++ ){
         const int jump2  = indices->getDMRGcumulative( irrep2 );
         const int NDMRG2 = indices->getNDMRG( irrep2 );
         const int irr_12 = Irreps::directProd( irrep1, irrep2 );
         for ( int irrep3 = 0; irrep3 < num_irreps; irrep3++ ){
            const int irrep4 = Irreps::directProd( irr_12, irrep3 );
            const int jump3  = indices->getDMRGcumulative( irrep3 );
            const int jump4  = indices->getDMRGcumulative( irrep4 );
            const int NDMRG3 = indices->getNDMRG( irrep3 );
            const int NDMRG4 = indices->getNDMRG( irrep4 );

            for ( int i1 = 0; i1 < NDMRG1; i1++ ){
               for ( int i2 = 0; i2 < NDMRG2; i2++ ){
                  for ( int i3 = 0; i3 < NDMRG3; i3++ ){
                     for ( int i4 = 0; i4 < NDMRG4; i4++ ){

                        double value = 0.0;
                        int jumpx = 0;
                        for ( int irrepx = 0; irrepx < num_irreps; irrepx++ ){

                           const int NOCCx  = indices->getNOCC(  irrepx );
                           const int NDMRGx = indices->getNDMRG( irrepx );

                           for ( int rowx = 0; rowx < NDMRGx; rowx++ ){
                              for ( int colx = 0; colx < NDMRGx; colx++ ){
                                 value += ( fock->get( irrepx, NOCCx + rowx, NOCCx + colx )
                                          * three_rdm[ jump1 + i1 + LAS*( jump2 + i2 + LAS*( jumpx + rowx + LAS*( jump3 + i3 + LAS*( jump4 + i4 + LAS*( jumpx + colx ))))) ] );
                              }
                           }
                           jumpx += NDMRGx;
                        }
                        f_dot_3dm[ jump1 + i1 + LAS * ( jump2 + i2 + LAS * ( jump3 + i3 + LAS * ( jump4 + i4 ) ) ) ] = value;
                     }
                  }
               }
            }
         }
      }
   }

   const double expectation_value = 2 * sum_f_kk() + f_dot_1dm;
   cout << "CASPT2::create_f_dots : < F > = " << expectation_value << endl;
   return expectation_value;

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
   
   const long long check = total_vector_length();
   assert( check == total_size );

   return total_size;

}

long long CheMPS2::CASPT2::total_vector_length() const{

   long long length = 0;
   for ( int i1 = 0; i1 < num_irreps; i1++ ){
      for ( int i2 = 0; i2 < num_irreps; i2++ ){
         for ( int i3 = 0; i3 < num_irreps; i3++ ){
            const int i4 = Irreps::directProd( Irreps::directProd( i1, i2 ), i3 );
            length +=     indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNDMRG( i4 ); // A:  E_ti E_uv | 0 >
            length +=     indices->getNDMRG( i1 ) * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // C:  E_at E_uv | 0 >
            length += 2 * indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // D:  E_ai E_tu | 0 >  and  E_ti E_au | 0 >
            length +=     indices->getNOCC( i1 )  * indices->getNOCC( i2 )  * indices->getNDMRG( i3 ) * indices->getNVIRT( i4 ); // E:  E_ti E_aj | 0 >
            length +=     indices->getNOCC( i1 )  * indices->getNDMRG( i2 ) * indices->getNVIRT( i3 ) * indices->getNVIRT( i4 ); // G:  E_ai E_bt | 0 >
            if ( i2 < i4 ){
               // 2 < 4 and irrep_2 < irrep_4
               length += indices->getNDMRG( i1 ) * indices->getNDMRG( i3 ) * indices->getNOCC( i2 )  * indices->getNOCC( i4 );  // B:  E_ti E_uj | 0 >
               length += indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNOCC( i2 )  * indices->getNOCC( i4 );  // H:  E_ai E_bj | 0 >
               length += indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNDMRG( i2 ) * indices->getNDMRG( i4 ); // F:  E_at E_bu | 0 >
            }
            if ( i2 == i4 ){ // i2 == i4 implies i1 == i3
               // 2 < 4 and irrep_2 == irrep_4
               length += ( indices->getNDMRG( i1 ) * indices->getNDMRG( i3 ) * indices->getNOCC( i2 )  * ( indices->getNOCC( i2 )  - 1 ) ) / 2; // B:  E_ti E_uj | 0 >
               length += ( indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNOCC( i2 )  * ( indices->getNOCC( i2 )  - 1 ) ) / 2; // H:  E_ai E_bj | 0 >
               length += ( indices->getNVIRT( i1 ) * indices->getNVIRT( i3 ) * indices->getNDMRG( i2 ) * ( indices->getNDMRG( i2 ) - 1 ) ) / 2; // F:  E_at E_bu | 0 >
               // 2 == 4 and 1 <= 3
               length += ( indices->getNDMRG( i1 ) * ( indices->getNDMRG( i3 ) + 1 ) * indices->getNOCC( i2 )  ) / 2; // B:  E_ti E_uj | 0 >
               length += ( indices->getNVIRT( i1 ) * ( indices->getNVIRT( i3 ) + 1 ) * indices->getNDMRG( i2 ) ) / 2; // F:  E_at E_bu | 0 >
               length += ( indices->getNVIRT( i1 ) * ( indices->getNVIRT( i3 ) + 1 ) * indices->getNOCC( i2 )  ) / 2; // H:  E_ai E_bj | 0 >
            }
         }
      }
   }

   return length;

}

void CheMPS2::CASPT2::matvec_blockdiag( double * vector, double * result, const double fock_prefactor, const double ovlp_prefactor ) const{

   /*
         TODO  | A  Bsinglet  Btriplet  C     D1     D2    Esinglet  Etriplet  Fsinglet  Ftriplet  Gsinglet  Gtriplet  Hsinglet  Htriplet
      ---------+-------------------------------------------------------------------------------------------------------------------------
      A        | OK    x       x        0     x      x     GRAD      GRAD      0         0         0         0         0         0
      Bsinglet | x     OK      0        0     0      0     x         x         0         0         0         0         0         0
      Btriplet | x     0       OK       0     0      0     x         x         0         0         0         0         0         0
      C        | 0     0       0        OK    x      x     0         0         x         x         GRAD      GRAD      0         0
      D1       | x     0       0        x     OK     OK    x         x         0         0         x         x         GRAD      GRAD
      D2       | x     0       0        x     OK     OK    x         x         0         0         x         x         GRAD      GRAD
      Esinglet | GRAD  x       x        0     x      x     OK        0         0         0         0         0         x         x
      Etriplet | GRAD  x       x        0     x      x     0         OK        0         0         0         0         x         x
      Fsinglet | 0     0       0        x     0      0     0         0         OK        0         x         x         0         0
      Ftriplet | 0     0       0        x     0      0     0         0         0         OK        x         x         0         0
      Gsinglet | 0     0       0        GRAD  x      x     0         0         x         x         OK        0         x         x
      Gtriplet | 0     0       0        GRAD  x      x     0         0         x         x         0         OK        x         x
      Hsinglet | 0     0       0        0     GRAD   GRAD  x         x         0         0         x         x         OK        0
      Htriplet | 0     0       0        0     GRAD   GRAD  x         x         0         0         x         x         0         OK
      
   */

   char notrans = 'N';
   int one      = 1;
   double set   = 0.0;
   double add   = 1.0;

   const double shifted_prefactor = ovlp_prefactor + 2 * fock_prefactor * sum_f_kk();

   double timings[ CHEMPS2_CASPT2_NUM_CASES ];
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES; cnt++ ){ timings[ cnt ] = 0.0; }
   struct timeval start, end;

   // FAA: < E_zy E_jx ( f_pq E_pq ) E_ti E_uv > = delta_ji * ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_AC[ irrep ];
      if ( SIZE > 0 ){
         const int NOCC = indices->getNOCC( irrep );
         for ( int count = 0; count < NOCC; count++ ){
            double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
            double alpha = fock_prefactor;
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FAA[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
            alpha = shifted_prefactor - fock_prefactor * fock->get( irrep, count, count );
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SAA[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_A ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FCC: < E_zy E_xb ( f_pq E_pq ) E_at E_uv > = delta_ba * ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_AC[ irrep ];
      if ( SIZE > 0 ){
         const int NVIR = indices->getNVIRT( irrep );
         const int N_OA = indices->getNOCC( irrep ) + indices->getNDMRG( irrep );
         for ( int count = 0; count < NVIR; count++ ){
            double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
            double alpha = fock_prefactor;
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FCC[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
            alpha = shifted_prefactor + fock_prefactor * fock->get( irrep, N_OA + count, N_OA + count );
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SCC[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_C ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FDD: < E_yx E_jb ( f_pq E_pq) E_ai E_tu > = delta_ab delta_ji ( FDD[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ xytu ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_D[ irrep ];
      if ( SIZE > 0 ){
         int jump_ai = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_a = Irreps::directProd( irrep_i, irrep );
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            for ( int i = 0; i < NOCC_i; i++ ){
               for ( int a = 0; a < NVIR_a; a++ ){
                  const int count = jump_ai + i + NOCC_i * a;
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  double alpha = fock_prefactor;
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FDD[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                  alpha = shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a ) - fock->get( irrep_i, i, i ) );
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SDD[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
               }
            }
            jump_ai += NOCC_i * NVIR_a;
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_D ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FBB singlet: < S_xkyl | f_pq E_pq | S_tiuj > = 2 ( delta_ki delta_lj + delta_kj delta_li ) ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               for ( int j = i; j < nocc_ij; j++ ){
                  const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                  const int ij_factor = (( i == j ) ? 2 : 1 );
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count; // irrep = 0
                  double alpha = 2 * ij_factor * fock_prefactor;
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FBB_singlet[ 0 ], &SIZE, origin, &SIZE, &set, target, &SIZE );
                  alpha = 2 * ij_factor * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_ij, i, i ) + fock->get( irrep_ij, j, j ) ) );
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SBB_singlet[ 0 ], &SIZE, origin, &SIZE, &add, target, &SIZE );
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_singlet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  for ( int j = 0; j < nocc_j; j++ ){
                     const int count = jump_ij + i + nocc_i * j;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     double alpha = 2 * fock_prefactor;
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FBB_singlet[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                     alpha = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_i, i, i ) + fock->get( irrep_j, j, j ) ) );
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SBB_singlet[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_B_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FBB triplet: < T_xkyl | f_pq E_pq | T_tiuj > = 2 ( delta_ki delta_lj - delta_kj delta_li ) ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               for ( int j = i+1; j < nocc_ij; j++ ){
                  const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count; // irrep = 0
                  double alpha = 2 * fock_prefactor;
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FBB_triplet[ 0 ], &SIZE, origin, &SIZE, &set, target, &SIZE );
                  alpha = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_ij, i, i ) + fock->get( irrep_ij, j, j ) ) );
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SBB_triplet[ 0 ], &SIZE, origin, &SIZE, &add, target, &SIZE );
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_triplet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  for ( int j = 0; j < nocc_j; j++ ){
                     const int count = jump_ij + i + nocc_i * j;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     double alpha = 2 * fock_prefactor;
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FBB_triplet[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                     alpha = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_i, i, i ) + fock->get( irrep_j, j, j ) ) );
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SBB_triplet[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_B_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FFF singlet: < S_cxdy | f_pq E_pq | S_atbu > = 2 ( delta_ac delta_bd + delta_ad delta_bc ) ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               for ( int b = a; b < NVIR_ab; b++ ){
                  const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                  const int ab_factor = (( a == b ) ? 2 : 1 );
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count; // irrep = 0
                  double alpha = 2 * ab_factor * fock_prefactor;
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FFF_singlet[ 0 ], &SIZE, origin, &SIZE, &set, target, &SIZE );
                  alpha = 2 * ab_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a )
                                                                                 + fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b ) ) );
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SFF_singlet[ 0 ], &SIZE, origin, &SIZE, &add, target, &SIZE );
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_singlet[ irrep ];
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
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const int count = jump_ab + a + NVIR_a * b;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     double alpha = 2 * fock_prefactor;
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FFF_singlet[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                     alpha = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                        + fock->get( irrep_b, N_OA_b + b, N_OA_b + b ) ) );
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SFF_singlet[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_F_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FFF triplet: < T_cxdy | f_pq E_pq | T_atbu > = 2 ( delta_ac delta_bd - delta_ad delta_bc ) ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               for ( int b = a+1; b < NVIR_ab; b++ ){
                  const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count; // irrep = 0
                  double alpha = 2 * fock_prefactor;
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FFF_triplet[ 0 ], &SIZE, origin, &SIZE, &set, target, &SIZE );
                  alpha = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a ) + fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b ) ) );
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SFF_triplet[ 0 ], &SIZE, origin, &SIZE, &add, target, &SIZE );
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_triplet[ irrep ];
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
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const int count = jump_ab + a + NVIR_a * b;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     double alpha = 2 * fock_prefactor;
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FFF_triplet[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                     alpha = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                        + fock->get( irrep_b, N_OA_b + b, N_OA_b + b ) ) );
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SFF_triplet[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_F_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FEE singlet: < S_ukbl | f_pq E_pq | S_tiaj > = 2 delta_ab ( delta_ki delta_lj + delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                        const int ij_factor = (( i == j ) ? 2 : 1 );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 2 * ij_factor * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FEE[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 2 * ij_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                          - fock->get( irrep_i,          i,          i )
                                                                                          - fock->get( irrep_j,          j,          j ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SEE[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 2 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FEE[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                              - fock->get( irrep_i,          i,          i )
                                                                              - fock->get( irrep_j,          j,          j ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SEE[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_E_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FEE triplet: < T_ukbl | f_pq E_pq | T_tiaj > = 6 delta_ab ( delta_ki delta_lj - delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int j = i+1; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 6 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FEE[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                              - fock->get( irrep_i,          i,          i )
                                                                              - fock->get( irrep_j,          j,          j ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SEE[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 6 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FEE[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                              - fock->get( irrep_i,          i,          i )
                                                                              - fock->get( irrep_j,          j,          j ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SEE[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_E_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FGG singlet: < S_cjdu | f_pq E_pq | S_aibt > = 2 delta_ji ( delta_ca delta_db + delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                        const int ab_factor = (( a == b ) ? 2 : 1 );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 2 * ab_factor * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FGG[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 2 * ab_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a )
                                                                                          + fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b )
                                                                                          - fock->get( irrep_i,           i,           i ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SGG[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 2 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FGG[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                              + fock->get( irrep_b, N_OA_b + b, N_OA_b + b )
                                                                              - fock->get( irrep_i,          i,          i ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SGG[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_G_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FGG triplet: < T_cjdu | f_pq E_pq | T_aibt > = 6 delta_ji ( delta_ca delta_db - delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int b = a+1; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 6 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FGG[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a )
                                                                              + fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b )
                                                                              - fock->get( irrep_i,           i,           i ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SGG[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 6 * fock_prefactor;
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, FGG[irrep], &SIZE, origin, &SIZE, &set, target, &SIZE );
                           alpha = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                              + fock->get( irrep_b, N_OA_b + b, N_OA_b + b )
                                                                              - fock->get( irrep_i,          i,          i ) ) );
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, SGG[irrep], &SIZE, origin, &SIZE, &add, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_G_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FHH singlet
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
            double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj; // irrep = 0
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj; // irrep = 0
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a; b < NVIR_ab; b++ ){
                  const int cnt_ab = a + ( b * ( b + 1 ) ) / 2;
                  const int ab_factor = (( a == b ) ? 2 : 1 );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_ij, i, i );
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int cnt_ij = i + ( j * ( j + 1 ) ) / 2;
                        const int ij_factor = (( i == j ) ? 2 : 1 );
                        const double f_jj = fock->get( irrep_ij, j, j );
                        const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                        target[ cnt_ij + size_ij * cnt_ab ] = 4 * alpha * ij_factor * ab_factor * origin[ cnt_ij + size_ij * cnt_ab ];
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
            const int NOCC_i  = indices->getNOCC( irrep_i );
            const int NOCC_j  = indices->getNOCC( irrep_j );
            const int size_ij = NOCC_i * NOCC_j;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep, irrep_a );
               if ( irrep_a < irrep_b ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int size_ab = NVIR_a * NVIR_b;
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
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
                              const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                              target[ cnt_ij + size_ij * cnt_ab ] = 4 * alpha * origin[ cnt_ij + size_ij * cnt_ab ];
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_H_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FHH triplet
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj; // irrep = 0
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj; // irrep = 0
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
                        const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                        target[ cnt_ij + size_ij * cnt_ab ] = 12 * alpha * origin[ cnt_ij + size_ij * cnt_ab ];
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
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
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
                              const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                              target[ cnt_ij + size_ij * cnt_ab ] = 12 * alpha * origin[ cnt_ij + size_ij * cnt_ab ];
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_H_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   double total = 0.0;
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES; cnt++ ){ total += timings[ cnt ]; }
   cout << "Application of [ blockdiag(F) - E * S ] : Wall time = " << total << " seconds" << endl;
   /*cout << "      A         : " << timings[ CHEMPS2_CASPT2_A         ] << endl;
   cout << "      B singlet : " << timings[ CHEMPS2_CASPT2_B_SINGLET ] << endl;
   cout << "      B triplet : " << timings[ CHEMPS2_CASPT2_B_TRIPLET ] << endl;
   cout << "      C         : " << timings[ CHEMPS2_CASPT2_C         ] << endl;
   cout << "      D         : " << timings[ CHEMPS2_CASPT2_D         ] << endl;
   cout << "      E singlet : " << timings[ CHEMPS2_CASPT2_E_SINGLET ] << endl;
   cout << "      E triplet : " << timings[ CHEMPS2_CASPT2_E_TRIPLET ] << endl;
   cout << "      F singlet : " << timings[ CHEMPS2_CASPT2_F_SINGLET ] << endl;
   cout << "      F triplet : " << timings[ CHEMPS2_CASPT2_F_TRIPLET ] << endl;
   cout << "      G singlet : " << timings[ CHEMPS2_CASPT2_G_SINGLET ] << endl;
   cout << "      G triplet : " << timings[ CHEMPS2_CASPT2_G_TRIPLET ] << endl;
   cout << "      H singlet : " << timings[ CHEMPS2_CASPT2_H_SINGLET ] << endl;
   cout << "      H triplet : " << timings[ CHEMPS2_CASPT2_H_TRIPLET ] << endl;*/

}

void CheMPS2::CASPT2::inverse_blockdiag( double * vector, double * result, const double fock_prefactor, const double ovlp_prefactor ) const{
   
   char notrans = 'N';
   char trans   = 'T';
   char uplo    = 'U';
   char jobz    = 'V';
   int one      = 1;
   double set   = 0.0;
   
   int maxsize = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      maxsize = max( max( max( max( max( size_AC[irrep],
                                         size_D[irrep] ),
                                         size_BF_singlet[irrep] ),
                                         size_BF_triplet[irrep] ),
                                         indices->getNDMRG( irrep ) ),
                                         maxsize );
   }
   if ( maxsize <= 2 ){ maxsize = 3; }
   
   double * work1 = new double[ maxsize * maxsize ];
   double * work2 = new double[ maxsize * maxsize ];
   double * eigs  = new double[ maxsize * maxsize ];
   int lwork = maxsize * maxsize;

   const double shifted_prefactor = ovlp_prefactor + 2 * fock_prefactor * sum_f_kk();

   double timings[ CHEMPS2_CASPT2_NUM_CASES ];
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES; cnt++ ){ timings[ cnt ] = 0.0; }
   struct timeval start, end;

   // FAA: < E_zy E_jx ( f_pq E_pq ) E_ti E_uv > = delta_ji * ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_AC[ irrep ];
      if ( SIZE > 0 ){
         const int NOCC = indices->getNOCC( irrep );
         for ( int count = 0; count < NOCC; count++ ){
            const double alpha1 = fock_prefactor;
            const double alpha2 = shifted_prefactor - fock_prefactor * fock->get( irrep, count, count );
            for ( int col = 0; col < SIZE; col++ ){
               for ( int row = 0; row < SIZE; row++ ){
                  work1[ row + SIZE * col ] = alpha1 * FAA[ irrep ][ row + SIZE * col ] + alpha2 * SAA[ irrep ][ row + SIZE * col ];
               }
            }
            int info;
            dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
            /*for ( int col = 0; col < SIZE; col++ ){
               cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
            }*/
            double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
            double alpha = 1.0;
            dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
            for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_A ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FCC: < E_zy E_xb ( f_pq E_pq ) E_at E_uv > = delta_ba * ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_AC[ irrep ];
      if ( SIZE > 0 ){
         const int NVIR = indices->getNVIRT( irrep );
         const int N_OA = indices->getNOCC( irrep ) + indices->getNDMRG( irrep );
         for ( int count = 0; count < NVIR; count++ ){
            const double alpha1 = fock_prefactor;
            const double alpha2 = shifted_prefactor + fock_prefactor * fock->get( irrep, N_OA + count, N_OA + count );
            for ( int col = 0; col < SIZE; col++ ){
               for ( int row = 0; row < SIZE; row++ ){
                  work1[ row + SIZE * col ] = alpha1 * FCC[ irrep ][ row + SIZE * col ] + alpha2 * SCC[ irrep ][ row + SIZE * col ];
               }
            }
            int info;
            dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
            /*for ( int col = 0; col < SIZE; col++ ){
               cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
            }*/
            double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
            double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
            double alpha = 1.0;
            dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
            for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
            dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_C ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FDD: < E_yx E_jb ( f_pq E_pq) E_ai E_tu > = delta_ab delta_ji ( FDD[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ xytu ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = size_D[ irrep ];
      if ( SIZE > 0 ){
         int jump_ai = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_a = Irreps::directProd( irrep_i, irrep );
            const int NOCC_i = indices->getNOCC( irrep_i );
            const int NVIR_a = indices->getNVIRT( irrep_a );
            const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
            for ( int i = 0; i < NOCC_i; i++ ){
               for ( int a = 0; a < NVIR_a; a++ ){
                  const double alpha1 = fock_prefactor;
                  const double alpha2 = shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a ) - fock->get( irrep_i, i, i ) );
                  for ( int col = 0; col < SIZE; col++ ){
                     for ( int row = 0; row < SIZE; row++ ){
                        work1[ row + SIZE * col ] = alpha1 * FDD[ irrep ][ row + SIZE * col ] + alpha2 * SDD[ irrep ][ row + SIZE * col ];
                     }
                  }
                  int info;
                  dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                  /*for ( int col = 0; col < SIZE; col++ ){
                     cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                  }*/
                  const int count = jump_ai + i + NOCC_i * a;
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * count;
                  double alpha = 1.0;
                  dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                  for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
               }
            }
            jump_ai += NOCC_i * NVIR_a;
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_D ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FBB singlet: < S_xkyl | f_pq E_pq | S_tiuj > = 2 ( delta_ki delta_lj + delta_kj delta_li ) ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               for ( int j = i; j < nocc_ij; j++ ){
                  const int ij_factor = (( i == j ) ? 2 : 1 );
                  const double alpha1 = 2 * ij_factor * fock_prefactor;
                  const double alpha2 = 2 * ij_factor * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_ij, i, i ) + fock->get( irrep_ij, j, j ) ) );
                  for ( int col = 0; col < SIZE; col++ ){
                     for ( int row = 0; row < SIZE; row++ ){
                        work1[ row + SIZE * col ] = alpha1 * FBB_singlet[ 0 ][ row + SIZE * col ] + alpha2 * SBB_singlet[ 0 ][ row + SIZE * col ];
                     }
                  }
                  int info;
                  dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                  /*for ( int col = 0; col < SIZE; col++ ){
                     cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                  }*/
                  const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count; // irrep = 0
                  double alpha = 1.0;
                  dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                  for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_singlet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double alpha1 = 2 * fock_prefactor;
                     const double alpha2 = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_i, i, i ) + fock->get( irrep_j, j, j ) ) );
                     for ( int col = 0; col < SIZE; col++ ){
                        for ( int row = 0; row < SIZE; row++ ){
                           work1[ row + SIZE * col ] = alpha1 * FBB_singlet[ irrep ][ row + SIZE * col ] + alpha2 * SBB_singlet[ irrep ][ row + SIZE * col ];
                        }
                     }
                     int info;
                     dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                     /*for ( int col = 0; col < SIZE; col++ ){
                        cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                     }*/
                     const int count = jump_ij + i + nocc_i * j;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * count;
                     double alpha = 1.0;
                     dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                     for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_B_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FBB triplet: < T_xkyl | f_pq E_pq | T_tiuj > = 2 ( delta_ki delta_lj - delta_kj delta_li ) ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int nocc_ij = indices->getNOCC( irrep_ij );
            for ( int i = 0; i < nocc_ij; i++ ){
               for ( int j = i+1; j < nocc_ij; j++ ){
                  const double alpha1 = 2 * fock_prefactor;
                  const double alpha2 = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_ij, i, i ) + fock->get( irrep_ij, j, j ) ) );
                  for ( int col = 0; col < SIZE; col++ ){
                     for ( int row = 0; row < SIZE; row++ ){
                        work1[ row + SIZE * col ] = alpha1 * FBB_triplet[ 0 ][ row + SIZE * col ] + alpha2 * SBB_triplet[ 0 ][ row + SIZE * col ];
                     }
                  }
                  int info;
                  dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                  /*for ( int col = 0; col < SIZE; col++ ){
                     cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                  }*/
                  const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count; // irrep = 0
                  double alpha = 1.0;
                  dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                  for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
               }
            }
            jump_ij += ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_triplet[ irrep ];
      if ( SIZE > 0 ){
         int jump_ij = 0;
         for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
            const int irrep_j = Irreps::directProd( irrep, irrep_i );
            if ( irrep_i < irrep_j ){
               const int nocc_i = indices->getNOCC( irrep_i );
               const int nocc_j = indices->getNOCC( irrep_j );
               for ( int i = 0; i < nocc_i; i++ ){
                  for ( int j = 0; j < nocc_j; j++ ){
                     const double alpha1 = 2 * fock_prefactor;
                     const double alpha2 = 2 * ( shifted_prefactor - fock_prefactor * ( fock->get( irrep_i, i, i ) + fock->get( irrep_j, j, j ) ) );
                     for ( int col = 0; col < SIZE; col++ ){
                        for ( int row = 0; row < SIZE; row++ ){
                           work1[ row + SIZE * col ] = alpha1 * FBB_triplet[ irrep ][ row + SIZE * col ] + alpha2 * SBB_triplet[ irrep ][ row + SIZE * col ];
                        }
                     }
                     int info;
                     dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                     /*for ( int col = 0; col < SIZE; col++ ){
                        cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                     }*/
                     const int count = jump_ij + i + nocc_i * j;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * count;
                     double alpha = 1.0;
                     dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                     for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                  }
               }
               jump_ij += nocc_i * nocc_j;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_B_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FFF singlet: < S_cxdy | f_pq E_pq | S_atbu > = 2 ( delta_ac delta_bd + delta_ad delta_bc ) ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_singlet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               for ( int b = a; b < NVIR_ab; b++ ){
                  const int ab_factor = (( a == b ) ? 2 : 1 );
                  const double alpha1 = 2 * ab_factor * fock_prefactor;
                  const double alpha2 = 2 * ab_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a )
                                                                                               + fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b ) ) );
                  for ( int col = 0; col < SIZE; col++ ){
                     for ( int row = 0; row < SIZE; row++ ){
                        work1[ row + SIZE * col ] = alpha1 * FFF_singlet[ 0 ][ row + SIZE * col ] + alpha2 * SFF_singlet[ 0 ][ row + SIZE * col ];
                     }
                  }
                  int info;
                  dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                  /*for ( int col = 0; col < SIZE; col++ ){
                     cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                  }*/
                  const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count; // irrep = 0
                  double alpha = 1.0;
                  dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                  for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_singlet[ irrep ];
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
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double alpha1 = 2 * fock_prefactor;
                     const double alpha2 = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                      + fock->get( irrep_b, N_OA_b + b, N_OA_b + b ) ) );
                     for ( int col = 0; col < SIZE; col++ ){
                        for ( int row = 0; row < SIZE; row++ ){
                           work1[ row + SIZE * col ] = alpha1 * FFF_singlet[ irrep ][ row + SIZE * col ] + alpha2 * SFF_singlet[ irrep ][ row + SIZE * col ];
                        }
                     }
                     int info;
                     dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                     /*for ( int col = 0; col < SIZE; col++ ){
                        cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                     }*/
                     const int count = jump_ab + a + NVIR_a * b;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * count;
                     double alpha = 1.0;
                     dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                     for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_F_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FFF triplet: < T_cxdy | f_pq E_pq | T_atbu > = 2 ( delta_ac delta_bd - delta_ad delta_bc ) ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int SIZE = size_BF_triplet[ 0 ];
      if ( SIZE > 0 ){
         int jump_ab = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            for ( int a = 0; a < NVIR_ab; a++ ){
               for ( int b = a+1; b < NVIR_ab; b++ ){
                  const double alpha1 = 2 * fock_prefactor;
                  const double alpha2 = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a )
                                                                                   + fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b ) ) );
                  for ( int col = 0; col < SIZE; col++ ){
                     for ( int row = 0; row < SIZE; row++ ){
                        work1[ row + SIZE * col ] = alpha1 * FFF_triplet[ 0 ][ row + SIZE * col ] + alpha2 * SFF_triplet[ 0 ][ row + SIZE * col ];
                     }
                  }
                  int info;
                  dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                  /*for ( int col = 0; col < SIZE; col++ ){
                     cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                  }*/
                  const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                  double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count; // irrep = 0
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count; // irrep = 0
                  double alpha = 1.0;
                  dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                  for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                  dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
               }
            }
            jump_ab += ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
         }
      }
   }
   for ( int irrep = 1; irrep < num_irreps; irrep++ ){
      int SIZE = size_BF_triplet[ irrep ];
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
                  for ( int b = 0; b < NVIR_b; b++ ){
                     const double alpha1 = 2 * fock_prefactor;
                     const double alpha2 = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                      + fock->get( irrep_b, N_OA_b + b, N_OA_b + b ) ) );
                     for ( int col = 0; col < SIZE; col++ ){
                        for ( int row = 0; row < SIZE; row++ ){
                           work1[ row + SIZE * col ] = alpha1 * FFF_triplet[ irrep ][ row + SIZE * col ] + alpha2 * SFF_triplet[ irrep ][ row + SIZE * col ];
                        }
                     }
                     int info;
                     dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                     /*for ( int col = 0; col < SIZE; col++ ){
                        cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                     }*/
                     const int count = jump_ab + a + NVIR_a * b;
                     double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * count;
                     double alpha = 1.0;
                     dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                     for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                     dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                  }
               }
               jump_ab += NVIR_a * NVIR_b;
            }
         }
      }
   }
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_F_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FEE singlet: < S_ukbl | f_pq E_pq | S_tiaj > = 2 delta_ab ( delta_ki delta_lj + delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j + 1 ) ) / 2;
                        const int ij_factor = (( i == j ) ? 2 : 1 );
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double alpha1 = 2 * ij_factor * fock_prefactor;
                           const double alpha2 = 2 * ij_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                                        - fock->get( irrep_i,          i,          i )
                                                                                                        - fock->get( irrep_j,          j,          j ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FEE[ irrep ][ row + SIZE * col ] + alpha2 * SEE[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double alpha1 = 2 * fock_prefactor;
                           const double alpha2 = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                            - fock->get( irrep_i,          i,          i )
                                                                                            - fock->get( irrep_j,          j,          j ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FEE[ irrep ][ row + SIZE * col ] + alpha2 * SEE[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_E_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   // FEE triplet: < T_ukbl | f_pq E_pq | T_tiaj > = 6 delta_ab ( delta_ki delta_lj - delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int j = i+1; j < NOCC_ij; j++ ){
                        const int count = jump_ij + i + ( j * ( j - 1 ) ) / 2;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double alpha1 = 6 * fock_prefactor;
                           const double alpha2 = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                            - fock->get( irrep_i,          i,          i )
                                                                                            - fock->get( irrep_j,          j,          j ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FEE[ irrep ][ row + SIZE * col ] + alpha2 * SEE[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
                        }
                     }
                  }
                  jump_ij += ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
               }
               if ( irrep_i < irrep_j ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  for ( int i = 0; i < NOCC_i; i++ ){
                     for ( int j = 0; j < NOCC_j; j++ ){
                        const int count = jump_ij + i + NOCC_i * j;
                        for ( int a = 0; a < NVIR_a; a++ ){
                           const double alpha1 = 6 * fock_prefactor;
                           const double alpha2 = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                            - fock->get( irrep_i,          i,          i )
                                                                                            - fock->get( irrep_j,          j,          j ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FEE[ irrep ][ row + SIZE * col ] + alpha2 * SEE[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( jump_aij + a + NVIR_a * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_E_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FGG singlet: < S_cjdu | f_pq E_pq | S_aibt > = 2 delta_ji ( delta_ca delta_db + delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int b = a; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b + 1 ) ) / 2;
                        const int ab_factor = (( a == b ) ? 2 : 1 );
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double alpha1 = 2 * ab_factor * fock_prefactor;
                           const double alpha2 = 2 * ab_factor * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a )
                                                                                                        + fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b )
                                                                                                        - fock->get( irrep_i,           i,           i ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FGG[ irrep ][ row + SIZE * col ] + alpha2 * SGG[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double alpha1 = 2 * fock_prefactor;
                           const double alpha2 = 2 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                            + fock->get( irrep_b, N_OA_b + b, N_OA_b + b )
                                                                                            - fock->get( irrep_i,          i,          i ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FGG[ irrep ][ row + SIZE * col ] + alpha2 * SGG[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_G_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FGG triplet: < T_cjdu | f_pq E_pq | T_aibt > = 6 delta_ji ( delta_ca delta_db - delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
   gettimeofday( &start, NULL );
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      int SIZE = indices->getNDMRG( irrep );
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
                     for ( int b = a+1; b < NVIR_ab; b++ ){
                        const int count = jump_ab + a + ( b * ( b - 1 ) ) / 2;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double alpha1 = 6 * fock_prefactor;
                           const double alpha2 = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_ab + a, N_OA_ab + a )
                                                                                            + fock->get( irrep_b, N_OA_ab + b, N_OA_ab + b )
                                                                                            - fock->get( irrep_i,           i,           i ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FGG[ irrep ][ row + SIZE * col ] + alpha2 * SGG[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
                     for ( int b = 0; b < NVIR_b; b++ ){
                        const int count = jump_ab + a + NVIR_a * b;
                        for ( int i = 0; i < NOCC_i; i++ ){
                           const double alpha1 = 6 * fock_prefactor;
                           const double alpha2 = 6 * ( shifted_prefactor + fock_prefactor * ( fock->get( irrep_a, N_OA_a + a, N_OA_a + a )
                                                                                            + fock->get( irrep_b, N_OA_b + b, N_OA_b + b )
                                                                                            - fock->get( irrep_i,          i,          i ) ) );
                           for ( int col = 0; col < SIZE; col++ ){
                              for ( int row = 0; row < SIZE; row++ ){
                                 work1[ row + SIZE * col ] = alpha1 * FGG[ irrep ][ row + SIZE * col ] + alpha2 * SGG[ irrep ][ row + SIZE * col ];
                              }
                           }
                           int info;
                           dsyev_( &jobz, &uplo, &SIZE, work1, &SIZE, eigs, work2, &lwork, &info );
                           /*for ( int col = 0; col < SIZE; col++ ){
                              cout << "eigs[ " << col << " ] = " << eigs[ col ] << endl;
                           }*/
                           double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( jump_abi + i + NOCC_i * count );
                           double alpha = 1.0;
                           dgemm_( &trans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, origin, &SIZE, &set, work2, &SIZE );
                           for ( int col = 0; col < SIZE; col++ ){ work2[ col ] = work2[ col ] / eigs[ col ]; }
                           dgemm_( &notrans, &notrans, &SIZE, &one, &SIZE, &alpha, work1, &SIZE, work2, &SIZE, &set, target, &SIZE );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_G_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   delete [] work1;
   delete [] work2;
   delete [] eigs;

   // FHH singlet
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
            double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj; // irrep = 0
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj; // irrep = 0
            for ( int a = 0; a < NVIR_ab; a++ ){
               const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
               for ( int b = a; b < NVIR_ab; b++ ){
                  const int cnt_ab = a + ( b * ( b + 1 ) ) / 2;
                  const int ab_factor = (( a == b ) ? 2 : 1 );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  for ( int i = 0; i < NOCC_ij; i++ ){
                     const double f_ii = fock->get( irrep_ij, i, i );
                     for ( int j = i; j < NOCC_ij; j++ ){
                        const int cnt_ij = i + ( j * ( j + 1 ) ) / 2;
                        const int ij_factor = (( i == j ) ? 2 : 1 );
                        const double f_jj = fock->get( irrep_ij, j, j );
                        const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                        target[ cnt_ij + size_ij * cnt_ab ] = origin[ cnt_ij + size_ij * cnt_ab ] / ( 4 * alpha * ij_factor * ab_factor );
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
            const int NOCC_i  = indices->getNOCC( irrep_i );
            const int NOCC_j  = indices->getNOCC( irrep_j );
            const int size_ij = NOCC_i * NOCC_j;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep, irrep_a );
               if ( irrep_a < irrep_b ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int size_ab = NVIR_a * NVIR_b;
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + jump_aibj;
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
                              const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                              target[ cnt_ij + size_ij * cnt_ab ] = origin[ cnt_ij + size_ij * cnt_ab ] / ( 4 * alpha );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_H_SINGLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   
   // FHH triplet
   gettimeofday( &start, NULL );
   { // First do irrep == 0
      int jump_aibj = 0;
      for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
         const int NOCC_ij = indices->getNOCC( irrep_ij );
         const int size_ij = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int size_ab = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            double * origin = vector + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj; // irrep = 0
            double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj; // irrep = 0
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
                        const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                        target[ cnt_ij + size_ij * cnt_ab ] = origin[ cnt_ij + size_ij * cnt_ab ] / ( 12 * alpha );
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
                  double * origin = vector + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + jump_aibj;
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
                              const double alpha = shifted_prefactor + fock_prefactor * ( f_dot_1dm + f_aa + f_bb - f_ii - f_jj );
                              target[ cnt_ij + size_ij * cnt_ab ] = origin[ cnt_ij + size_ij * cnt_ab ] / ( 12 * alpha );
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
   gettimeofday( &end, NULL );
   timings[ CHEMPS2_CASPT2_H_TRIPLET ] = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );

   double total = 0.0;
   for ( int cnt = 0; cnt < CHEMPS2_CASPT2_NUM_CASES; cnt++ ){ total += timings[ cnt ]; }
   cout << "Application of [ blockdiag(F) - E * S ]^{-1} : Wall time = " << total << " seconds" << endl;
   /*cout << "      A         : " << timings[ CHEMPS2_CASPT2_A         ] << endl;
   cout << "      B singlet : " << timings[ CHEMPS2_CASPT2_B_SINGLET ] << endl;
   cout << "      B triplet : " << timings[ CHEMPS2_CASPT2_B_TRIPLET ] << endl;
   cout << "      C         : " << timings[ CHEMPS2_CASPT2_C         ] << endl;
   cout << "      D         : " << timings[ CHEMPS2_CASPT2_D         ] << endl;
   cout << "      E singlet : " << timings[ CHEMPS2_CASPT2_E_SINGLET ] << endl;
   cout << "      E triplet : " << timings[ CHEMPS2_CASPT2_E_TRIPLET ] << endl;
   cout << "      F singlet : " << timings[ CHEMPS2_CASPT2_F_SINGLET ] << endl;
   cout << "      F triplet : " << timings[ CHEMPS2_CASPT2_F_TRIPLET ] << endl;
   cout << "      G singlet : " << timings[ CHEMPS2_CASPT2_G_SINGLET ] << endl;
   cout << "      G triplet : " << timings[ CHEMPS2_CASPT2_G_TRIPLET ] << endl;
   cout << "      H singlet : " << timings[ CHEMPS2_CASPT2_H_SINGLET ] << endl;
   cout << "      H triplet : " << timings[ CHEMPS2_CASPT2_H_TRIPLET ] << endl;*/

}

void CheMPS2::CASPT2::construct_rhs( const DMRGSCFintegrals * integrals ){

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
           < S_tiaj | H > = sum_w [ (aj|wi) + (ai|wj) ] * 1 * SEE[ It ][ wt ]
           < T_tiaj | H > = sum_w [ (aj|wi) - (ai|wj) ] * 3 * SEE[ It ][ wt ]

      VF:  < H E_at E_bu >
           < S_atbu | H > = sum_xy (ax|by) SFF_singlet[ It x Iu ][ xytu ]
           < T_atbu | H > = sum_xy (ax|by) SFF_triplet[ It x Iu ][ xytu ]

      VG:  < H E_ai E_bt >
           < S_aibt | H > = sum_u [ (ai|bu) + (bi|au) ] * 1 * SGG[ It ][ ut ]
           < T_aibt | H > = sum_u [ (ai|bu) - (bi|au) ] * 3 * SGG[ It ][ ut ]

      VH:  < H E_ai E_bj >
           < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ]
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

   vector_rhs = new double[ jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ] ];

   int max_size = 0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      max_size = max( max( max( max( size_AC[ irrep ], size_D[ irrep ] ), size_BF_singlet[ irrep ] ), size_BF_triplet[ irrep ] ), max_size );
   }
   double * workspace = new double[ max_size ];

   // VA
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
                        const double ix_zy = integrals->get_coulomb( irrep, irrep_x, irrep_z, irrep_y, count_i, occ_x + x, occ_z + z, occ_y + y );
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ] = ix_zy;
                        //assert( fabs( ix_zy - integrals->FourIndexAPI( irrep, irrep_z, irrep_x, irrep_y, count_i, occ_z + z, occ_x + x, occ_y + y ) ) < 1e-30 );
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

   // VC
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
                        const double zy_xa = integrals->get_coulomb( irrep_z, irrep_y, irrep_x, irrep, occ_z + z, occ_y + y, occ_x + x, N_OA + count_a );
                        workspace[ jump_xyz + x + num_x * ( y + num_y * z ) ] = zy_xa;
                        //assert( fabs( zy_xa - integrals->FourIndexAPI( irrep_z, irrep_x, irrep_y, irrep, occ_z + z, occ_x + x, occ_y + y, N_OA + count_a ) ) < 1e-30 );
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
                  const int occ_x = indices->getNOCC( irrep_x );
                  const int occ_y = indices->getNOCC( irrep_y );
                  const int num_x = indices->getNDMRG( irrep_x );
                  const int num_y = indices->getNDMRG( irrep_y );

                  // workspace[ xy ] = (ia|yx)
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double ia_yx = integrals->get_coulomb( irrep_y, irrep_x, irrep_i, irrep_a, occ_y + y, occ_x + x, count_i, N_OA_a + count_a );
                        workspace[ jump_xy + x + num_x * y ] = ia_yx;
                        //assert( fabs( ia_yx - integrals->FourIndexAPI( irrep_i, irrep_y, irrep_a, irrep_x, count_i, occ_y + y, N_OA_a + count_a, occ_x + x ) ) < 1e-30 );
                     }
                  }

                  // workspace[ D2JUMP + xy ] = (ix|ya)
                  for ( int y = 0; y < num_y; y++ ){
                     for ( int x = 0; x < num_x; x++ ){
                        const double ix_ya = integrals->get_coulomb( irrep_i, irrep_x, irrep_y, irrep_a, count_i, occ_x + x, occ_y + y, N_OA_a + count_a );
                        workspace[ D2JUMP + jump_xy + x + num_x * y ] = ix_ya;
                        //assert( fabs( ix_ya - integrals->FourIndexAPI( irrep_i, irrep_y, irrep_x, irrep_a, count_i, occ_y + y, occ_x + x, N_OA_a + count_a ) ) < 1e-30 );
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

   delete MAT;

   // VB singlet and triplet
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
                        const double ix_jy = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y );
                        workspace[ jump_xy + x + (y*(y+1))/2 ] = ix_jy;
                        //assert( fabs( ix_jy - integrals->FourIndexAPI( irrep_ij, irrep_ij, irrep_xy, irrep_xy, i, j, occ_xy + x, occ_xy + y ) ) < 1e-30 );
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
                        const double ix_jy = integrals->get_coulomb( irrep_ij, irrep_xy, irrep_ij, irrep_xy, i, occ_xy + x, j, occ_xy + y );
                        workspace[ jump_xy + x + (y*(y-1))/2 ] = ix_jy;
                        //assert( fabs( ix_jy - integrals->FourIndexAPI( irrep_ij, irrep_ij, irrep_xy, irrep_xy, i, j, occ_xy + x, occ_xy + y ) ) < 1e-30 );
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
                              const double ix_jy = integrals->get_coulomb( irrep_i, irrep_x, irrep_j, irrep_y, i, occ_x + x, j, occ_y + y );
                              workspace[ jump_xy + x + num_x * y ] = ix_jy;
                              //assert( fabs( ix_jy - integrals->FourIndexAPI( irrep_i, irrep_j, irrep_x, irrep_y, i, j, occ_x + x, occ_y + y ) ) < 1e-30 );
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

   // VF singlet and triplet
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
                        const double ax_by = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b );
                        workspace[ jump_xy + x + (y*(y+1))/2 ] = ax_by;
                        //assert( fabs( ax_by - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_xy, irrep_xy, N_OA_ab + a, N_OA_ab + b, occ_xy + x, occ_xy + y ) ) < 1e-30 );
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
                        const double ax_by = integrals->get_exchange( irrep_xy, irrep_xy, irrep_ab, irrep_ab, occ_xy + x, occ_xy + y, N_OA_ab + a, N_OA_ab + b );
                        workspace[ jump_xy + x + (y*(y-1))/2 ] = ax_by;
                        //assert( fabs( ax_by - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_xy, irrep_xy, N_OA_ab + a, N_OA_ab + b, occ_xy + x, occ_xy + y ) ) < 1e-30 );
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
                        const int d_x = indices->getDMRGcumulative( irrep_x );
                        const int d_y = indices->getDMRGcumulative( irrep_y );
                        const int occ_x = indices->getNOCC( irrep_x );
                        const int occ_y = indices->getNOCC( irrep_y );
                        const int num_x = indices->getNDMRG( irrep_x );
                        const int num_y = indices->getNDMRG( irrep_y );

                        for ( int y = 0; y < num_y; y++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              const double ax_by = integrals->get_exchange( irrep_x, irrep_y, irrep_a, irrep_b, occ_x + x, occ_y + y, N_OA_a + a, N_OA_b + b );
                              workspace[ jump_xy + x + num_x * y ] = ax_by;
                              //assert( fabs( ax_by - integrals->FourIndexAPI( irrep_a, irrep_b, irrep_x, irrep_y, N_OA_a + a, N_OA_b + b, occ_x + x, occ_y + y ) ) < 1e-30 );
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
                              // < S_tiaj | H > = sum_w [ (aj|wi) + (ai|wj) ] * 1 * SEE[ It ][ wt ]
                              // < T_tiaj | H > = sum_w [ (aj|wi) - (ai|wj) ] * 3 * SEE[ It ][ wt ]
                              const double SEE_wt = SEE[ irrep ][ w + num_t * t ];
                              const double aj_wi  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, i, occ_t + w, j, N_OA_a + a );
                              const double ai_wj  = integrals->get_coulomb( irrep_ij, irrep, irrep_ij, irrep_a, j, occ_t + w, i, N_OA_a + a );
                              value_singlet +=     SEE_wt * ( aj_wi + ai_wj );
                              value_triplet += 3 * SEE_wt * ( aj_wi - ai_wj );
                              //assert( fabs( aj_wi - integrals->FourIndexAPI( irrep_a, irrep, irrep_ij, irrep_ij, N_OA_a + a, occ_t + w, j, i ) ) < 1e-30 );
                              //assert( fabs( ai_wj - integrals->FourIndexAPI( irrep_a, irrep, irrep_ij, irrep_ij, N_OA_a + a, occ_t + w, i, j ) ) < 1e-30 );
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
                                 // < S_tiaj | H > = sum_w [ (aj|wi) + (ai|wj) ] * 1 * SEE[ It ][ wt ]
                                 // < T_tiaj | H > = sum_w [ (aj|wi) - (ai|wj) ] * 3 * SEE[ It ][ wt ]
                                 const double SEE_wt = SEE[ irrep ][ w + num_t * t ];
                                 const double aj_wi  = integrals->get_coulomb( irrep_i, irrep, irrep_j, irrep_a, i, occ_t + w, j, N_OA_a + a );
                                 const double ai_wj  = integrals->get_coulomb( irrep_j, irrep, irrep_i, irrep_a, j, occ_t + w, i, N_OA_a + a );
                                 value_singlet +=     SEE_wt * ( aj_wi + ai_wj );
                                 value_triplet += 3 * SEE_wt * ( aj_wi - ai_wj );
                                 //assert( fabs( aj_wi - integrals->FourIndexAPI( irrep_a, irrep, irrep_j, irrep_i, N_OA_a + a, occ_t + w, j, i ) ) < 1e-30 );
                                 //assert( fabs( ai_wj - integrals->FourIndexAPI( irrep_a, irrep, irrep_i, irrep_j, N_OA_a + a, occ_t + w, i, j ) ) < 1e-30 );
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
                              // < S_aibt | H > = sum_u [ (ai|bu) + (bi|au) ] * 1 * SGG[ It ][ ut ]
                              // < T_aibt | H > = sum_u [ (ai|bu) - (bi|au) ] * 3 * SGG[ It ][ ut ]
                              const double SGG_ut = SGG[ irrep ][ u + num_t * t ];
                              const double ai_bu  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + a, N_OA_ab + b );
                              const double bi_au  = integrals->get_exchange( irrep_i, irrep, irrep_ab, irrep_ab, i, occ_t + u, N_OA_ab + b, N_OA_ab + a );
                              value_singlet +=     SGG_ut * ( ai_bu + bi_au );
                              value_triplet += 3 * SGG_ut * ( ai_bu - bi_au );
                              //assert( fabs( ai_bu - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_i, irrep, N_OA_ab + a, N_OA_ab + b, i, occ_t + u ) ) < 1e-30 );
                              //assert( fabs( bi_au - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_i, irrep, N_OA_ab + b, N_OA_ab + a, i, occ_t + u ) ) < 1e-30 );
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
                                 // < S_aibt | H > = sum_u [ (ai|bu) + (bi|au) ] * 1 * SGG[ It ][ ut ]
                                 // < T_aibt | H > = sum_u [ (ai|bu) - (bi|au) ] * 3 * SGG[ It ][ ut ]
                                 const double SGG_ut = SGG[ irrep ][ u + num_t * t ];
                                 const double ai_bu  = integrals->get_exchange( irrep_i, irrep, irrep_a, irrep_b, i, occ_t + u, N_OA_a + a, N_OA_b + b );
                                 const double bi_au  = integrals->get_exchange( irrep_i, irrep, irrep_b, irrep_a, i, occ_t + u, N_OA_b + b, N_OA_a + a );
                                 value_singlet +=     SGG_ut * ( ai_bu + bi_au );
                                 value_triplet += 3 * SGG_ut * ( ai_bu - bi_au );
                                 //assert( fabs( ai_bu - integrals->FourIndexAPI( irrep_a, irrep_b, irrep_i, irrep, N_OA_a + a, N_OA_b + b, i, occ_t + u ) ) < 1e-30 );
                                 //assert( fabs( bi_au - integrals->FourIndexAPI( irrep_b, irrep_a, irrep_i, irrep, N_OA_b + b, N_OA_a + a, i, occ_t + u ) ) < 1e-30 );
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
                           // < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ]
                           // < T_aibj | H > = 6 * [ (ai|bj) - (aj|bi) ]
                           const double ai_bj = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, i, j, noa_ab + a, noa_ab + b );
                           const double aj_bi = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, j, i, noa_ab + a, noa_ab + b );
                                target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi );
   if (( b > a ) && ( j > i )){ target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi ); }
                           //assert( fabs( ai_bj - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_ij, irrep_ij, noa_ab + a, noa_ab + b, i, j ) ) < 1e-30 );
                           //assert( fabs( aj_bi - integrals->FourIndexAPI( irrep_ab, irrep_ab, irrep_ij, irrep_ij, noa_ab + a, noa_ab + b, j, i ) ) < 1e-30 );
                        }
                     }
                  }
               }
               jump_aibj_singlet += ( nocc_ij * ( nocc_ij + 1 ) * nvirt_ab * ( nvirt_ab + 1 )) / 4;
               jump_aibj_triplet += ( nocc_ij * ( nocc_ij - 1 ) * nvirt_ab * ( nvirt_ab - 1 )) / 4;
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
                                 // < S_aibj | H > = 2 * [ (ai|bj) + (aj|bi) ]
                                 // < T_aibj | H > = 6 * [ (ai|bj) - (aj|bi) ]
                                 const double ai_bj = integrals->get_exchange( irrep_i, irrep_j, irrep_a, irrep_b, i, j, noa_a + a, noa_b + b );
                                 const double aj_bi = integrals->get_exchange( irrep_j, irrep_i, irrep_a, irrep_b, j, i, noa_a + a, noa_b + b );
                                 target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi );
                                 target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi );
                                 //assert( fabs( ai_bj - integrals->FourIndexAPI( irrep_a, irrep_b, irrep_i, irrep_j, noa_a + a, noa_b + b, i, j ) ) < 1e-30 );
                                 //assert( fabs( aj_bi - integrals->FourIndexAPI( irrep_a, irrep_b, irrep_j, irrep_i, noa_a + a, noa_b + b, j, i ) ) < 1e-30 );
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

double CheMPS2::CASPT2::sum_f_kk() const{

   double total = 0.0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      for ( int orb = 0; orb < NOCC; orb++ ){
         total += fock->get( irrep, orb, orb );
      }
   }
   return total;

}

void CheMPS2::CASPT2::make_FAA_FCC(){

   /*
      FAA: sum_cd f_cd < E_zy E_jx E_cd E_ti E_uv > = 0
      
           sum_kl f_kl < E_zy E_jx E_kl E_ti E_uv > = delta_ji ( 2 * ( sum_k f_kk ) - f_ij ) SAA[ Ii ][ xyztuv ]
                                                      ( because in the occupied-occupied block the fock operator is diagonal! )
           
           sum_rs f_rs < E_zy E_jx E_rs E_ti E_uv > ( with rs active indices )
           
              = delta_ji [ - f_dot_4dm[ ztuyxv ]
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
                         ]
                         
      FCC: sum_cd f_cd < E_zy E_xb E_cd E_at E_uv > = delta_ba f_ba SCC[ Ia ][ xyztuv ]
                                                      ( because in the virtual-virtual block the fock operator is diagonal! )
           
           sum_kl f_kl < E_zy E_xb E_kl E_at E_uv > = delta_ba ( 2 sum_k f_kk ) SCC[ Ia ][ xyztuv ]
           
           sum_rs f_rs < E_zy E_xb E_rs E_at E_uv > ( with rs active indices )
           
              = delta_ba [ + f_dot_4dm[ zxuytv ]
                           + delta_uy f_dot_3dm[ xztv ]
                           + delta_xy f_dot_3dm[ zutv ]
                           + delta_ut f_dot_3dm[ zxyv ]
                           + delta_ut delta_xy f_dot_2dm[ zv ]
                           + sum_s f_ys SCC[ Ia ] [ xsztuv ]
                           + sum_r f_ru SCC[ Ia ] [ xyztrv ]
                           + f_yx ( Gamma_{zutv} + delta_ut Gamma_{zv} )
                           + f_tu ( Gamma_{zxyv} + delta_yx Gamma_{zv} )
                           - f_yu Gamma_{zxvt}
                         ]

   */
   
   FAA = new double*[ num_irreps ];
   FCC = new double*[ num_irreps ];
   
   const int LAS = indices->getDMRGcumulative( num_irreps );
   
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int SIZE = size_AC[ irrep ];
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

                  for ( int t = 0; t < num_t; t++ ){ // FCC: + f_dot_4dm[ zxuytv ]
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
                     for ( int uy = 0; uy < num_u; uy++ ){
                        for ( int t = 0; t < num_t; t++ ){
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    FCC[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
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
                                    FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
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
                                    FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + ut + num_t * ( ut + num_u * v ) ) ]
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
                        for ( int y = 0; y < num_y; y++ ){
                           const double f_yx = fock->get( irrep_x, nocc_x + y, nocc_x + x );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int u = 0; u < num_u; u++ ){
                                 for ( int t = 0; t < num_t; t++ ){
                                    for ( int z = 0; z < num_z; z++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
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
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_tu = fock->get( irrep_t, nocc_t + t, nocc_t + u );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
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
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_yu = fock->get( irrep_y, nocc_y + y, nocc_y + u );
                           for ( int v = 0; v < num_v; v++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
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

void CheMPS2::CASPT2::make_FDD(){

   /*
      FD1D1: < E_yx E_jb ( f_pq E_pq) E_ai E_tu > = delta_ab delta_ji ( FD1D1[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SD1D1[ xytu ] )
      FD2D2: < E_yb E_jx ( f_pq E_pq) E_ti E_au > = delta_ab delta_ji ( FD2D2[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SD2D2[ xytu ] )
      FD1D2: < E_yx E_jb ( f_pq E_pq) E_ti E_au > = delta_ab delta_ji ( FD1D2[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SD1D2[ xytu ] )
      FD2D1: < E_yb E_jx ( f_pq E_pq) E_ai E_tu > = delta_ab delta_ji ( FD2D1[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SD2D1[ xytu ] )

            FD1D1[ xytu ] = + sum_w f_xw SD1D1[ wytu ]
                            + sum_v f_vt SD1D1[ xyvu ]
                            + ( 2 * f_dot_3dm[ ytxu ] + 2 * delta_tx f_dot_2dm[ yu ] )
                            - 2 f_xt Gamma_{yu}

            FD2D2[ xytu ] = + sum_w f_xw SD2D2[ wytu ]
                            + sum_v f_vt SD2D2[ xyvu ]
                            + (   - f_dot_3dm[ ytux ] + 2 * delta_tx f_dot_2dm[ yu ] )
                            - 2 f_xt Gamma_{yu}
                           
            FD1D2[ xytu ] = + sum_w f_xw SD1D2[ wytu ]
                            + sum_v f_vt SD1D2[ xyvu ]
                            + (   - f_dot_3dm[ ytxu ]     - delta_tx f_dot_2dm[ yu ] )
                            + f_xt Gamma_{yu}
                           
            FD2D1[ xytu ] = + sum_w f_xw SD2D1[ wytu ]
                            + sum_v f_vt SD2D1[ xyvu ]
                            + (   - f_dot_3dm[ ytxu ]     - delta_tx f_dot_2dm[ yu ] )
                            + f_xt Gamma_{yu}
      
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

void CheMPS2::CASPT2::make_FBB_FFF_singlet(){

   /*
      FBB singlet: < S_xkyl | ( f_pq E_pq ) | S_tiuj >
                      = 2 ( delta_ki delta_lj + delta_kj delta_li ) ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )

      FFF singlet: < S_cxdy | ( f_pq E_pq ) | S_atbu >
                      = 2 ( delta_ac delta_bd + delta_ad delta_bc ) ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )

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
      const int irrep = 0;

      const int SIZE = size_BF_singlet[ irrep ];
      FBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
      FFF_singlet[ irrep ] = new double[ SIZE * SIZE ];

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
                                        * SBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( v + ( t * ( t + 1 ) ) / 2 ) ] );
                        }
                        for ( int v = t; v < num_ut; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( v * ( v + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_v f_vt SBB_singlet[ xyvu ] = + sum_(v<u) f_vt SBB_singlet[ xyvu ] + sum_(v>=u) f_vt SBB_singlet[ xyuv ]
                        for ( int v = 0; v < u; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( v + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int v = u; v < num_ut; v++ ){
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( u + ( v * ( v + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_yw SBB_singlet[ xwtu ] = + sum_(w>=x) f_yw SBB_singlet[ xwtu ] + sum_(w<x) f_yw SBB_singlet[ wxtu ]
                        for ( int w = 0; w < x; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_singlet[ irrep ][ shift + w + ( x * ( x + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int w = x; w < num_xy; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_singlet[ irrep ][ shift + x + ( w * ( w + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_xw SBB_singlet[ wytu ] = + sum_(w<y) f_xw SBB_singlet[ wytu ] + sum_(w>=y) f_xw SBB_singlet[ ywtu ]
                        for ( int w = 0; w < y; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_singlet[ irrep ][ shift + w + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }
                        for ( int w = y; w < num_xy; w++ ){
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_singlet[ irrep ][ shift + y + ( w * ( w + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] );
                        }

                        FBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] = sbb_value;
                        FFF_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] = fdotsum;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy and d_ut = d_xy

               // FBB: + 2 ( delta_uy delta_tx + delta_ux delta_ty ) f_dot_1dm --> last term only if x <= y == t <= u or hence all equal
               for ( int t = 0; t < num_ut; t++ ){
                  FBB_singlet[ irrep ][ shift + t + ( t * ( t + 1 ) ) / 2 + SIZE * ( t + ( t * ( t + 1 ) ) / 2 ) ] += 4 * f_dot_1dm;
                  for ( int u = t+1; u < num_ut; u++ ){
                     FBB_singlet[ irrep ][ shift + t + ( u * ( u + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] += 2 * f_dot_1dm;
                  }
               }

               // FBB: - delta_uy f_dot_2dm[ tx ]
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                     for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                        const double val_tx = f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ];
                        FBB_singlet[ irrep ][ shift + x + ( uy * ( uy + 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy + 1 ) ) / 2 ) ] -= val_tx;
                     }
                  }
               }

               // FBB: - delta_tx f_dot_2dm[ uy ]
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                     for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                        const double val_uy = f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ];
                        FBB_singlet[ irrep ][ shift + tx + ( y * ( y + 1 ) ) / 2 + SIZE * ( tx + ( u * ( u + 1 ) ) / 2 ) ] -= val_uy;
                     }
                  }
               }

               // FBB: - delta_ux f_dot_2dm[ ty ]
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                     for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                        const double val_ty = f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ];
                        FBB_singlet[ irrep ][ shift + ux + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux + 1 ) ) / 2 ) ] -= val_ty;
                     }
                  }
               }

               // FBB: - delta_ty f_dot_2dm[ ux ]
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                     for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                        const double val_ux = f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ];
                        FBB_singlet[ irrep ][ shift + x + ( ty * ( ty + 1 ) ) / 2 + SIZE * ( ty + ( u * ( u + 1 ) ) / 2 ) ] -= val_ux;
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
                           FBB_singlet[ irrep ][ shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] -= value;
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

      const int SIZE = size_BF_singlet[ irrep ];
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

void CheMPS2::CASPT2::make_FBB_FFF_triplet(){

   /*
      FBB triplet: < T_xkyl | ( f_pq E_pq ) | T_tiuj >
                      = 2 ( delta_ki delta_lj - delta_kj delta_li ) ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )

      FFF triplet: < T_cxdy | ( f_pq E_pq ) | T_atbu >
                      = 2 ( delta_ac delta_bd - delta_ad delta_bc ) ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )

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
      const int irrep = 0;

      const int SIZE = size_BF_triplet[ irrep ];
      FBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
      FFF_triplet[ irrep ] = new double[ SIZE * SIZE ];

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

                        // FBB: + sum_v f_vu SBB_triplet[ xytv ] = + sum_(v>=t) f_vu SBB_triplet[ xytv ] + sum_(v<t) f_vu SBB_triplet[ xyvt ]
                        for ( int v = 0; v < t; v++ ){ // Swap ( t <--> u ) --> minus sign
                           sbb_value -= ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( v + ( t * ( t - 1 ) ) / 2 ) ] );
                        }
                        for ( int v = t+1; v < num_ut; v++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + u )
                                        * SBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( v * ( v - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_v f_vt SBB_triplet[ xyvu ] = + sum_(v<u) f_vt SBB_triplet[ xyvu ] + sum_(v>=u) f_vt SBB_triplet[ xyuv ]
                        for ( int v = 0; v < u; v++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( v + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int v = u+1; v < num_ut; v++ ){ // Swap ( t <--> u ) --> minus sign
                           sbb_value -= ( fock->get( irrep_ut, nocc_ut + v, nocc_ut + t )
                                        * SBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( u + ( v * ( v - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_yw SBB_triplet[ xwtu ] = + sum_(w>=x) f_yw SBB_triplet[ xwtu ] + sum_(w<x) f_yw SBB_triplet[ wxtu ]
                        for ( int w = 0; w < x; w++ ){ // Swap ( x <--> y ) --> minus sign
                           sbb_value -= ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_triplet[ irrep ][ shift + w + ( x * ( x - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int w = x+1; w < num_xy; w++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + y, nocc_xy + w )
                                        * SBB_triplet[ irrep ][ shift + x + ( w * ( w - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }

                        // FBB: + sum_w f_xw SBB_triplet[ wytu ] = + sum_(w<y) f_xw SBB_triplet[ wytu ] + sum_(w>=y) f_xw SBB_triplet[ ywtu ]
                        for ( int w = 0; w < y; w++ ){ // Natural replacement --> plus sign
                           sbb_value += ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_triplet[ irrep ][ shift + w + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }
                        for ( int w = y+1; w < num_xy; w++ ){ // Swap ( x <--> y ) --> minus sign
                           sbb_value -= ( fock->get( irrep_xy, nocc_xy + x, nocc_xy + w )
                                        * SBB_triplet[ irrep ][ shift + y + ( w * ( w - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] );
                        }

                        FBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] = sbb_value;
                        FFF_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] = fdotdiff;
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy and d_ut = d_xy

               // FBB: + 6 ( delta_uy delta_tx - delta_ux delta_ty ) f_dot_1dm --> last term only if x < y = t < u or NEVER
               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t+1; u < num_ut; u++ ){
                     FBB_triplet[ irrep ][ shift + t + ( u * ( u - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] += 6 * f_dot_1dm;
                  }
               }

               // FBB: - 3 delta_uy f_dot_2dm[ tx ]
               for ( int uy = 0; uy < num_ut; uy++ ){
                  for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                     for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                        const double val_tx = f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ];
                        FBB_triplet[ irrep ][ shift + x + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy - 1 ) ) / 2 ) ] -= 3 * val_tx;
                     }
                  }
               }

               // FBB: - 3 delta_tx f_dot_2dm[ uy ]
               for ( int tx = 0; tx < num_ut; tx++ ){
                  for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                     for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                        const double val_uy = f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ];
                        FBB_triplet[ irrep ][ shift + tx + ( y * ( y - 1 ) ) / 2 + SIZE * ( tx + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * val_uy;
                     }
                  }
               }

               // FBB: + 3 delta_ux f_dot_2dm[ ty ]
               for ( int ux = 0; ux < num_ut; ux++ ){
                  for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                     for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                        const double val_ty = f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ];
                        FBB_triplet[ irrep ][ shift + ux + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux - 1 ) ) / 2 ) ] += 3 * val_ty;
                     }
                  }
               }

               // FBB: + 3 delta_ty f_dot_2dm[ ux ]
               for ( int ty = 0; ty < num_ut; ty++ ){
                  for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                     for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                        const double val_ux = f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ];
                        FBB_triplet[ irrep ][ shift + x + ( ty * ( ty - 1 ) ) / 2 + SIZE * ( ty + ( u * ( u - 1 ) ) / 2 ) ] += 3 * val_ux;
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
                           FBB_triplet[ irrep ][ shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * value;
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

      const int SIZE = size_BF_triplet[ irrep ];
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

                              // + f_dot_3dm_{utyx} + f_dot_3dm_{utxy}
                              const double fdotsum = ( f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                     + f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );

                              double sbb_value = fdotsum;

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
                              FFF_triplet[ irrep ][ shift + x + num_x * y + SIZE * ( t + num_t * u ) ] = fdotsum;
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

void CheMPS2::CASPT2::make_FEE_FGG(){

   /*

      E singlet: | SE_tiaj > = ( E_ti E_aj + E_tj E_ai ) | 0 >  with  i <= j
      E triplet: | TE_tiaj > = ( E_ti E_aj - E_tj E_ai ) | 0 >  with  i <  j

         < SE_ukbl | ( sum_pq f_pq E_pq ) | SE_tiaj > = 2 delta_ab ( delta_ki delta_lj + delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
         < TE_ukbl | ( sum_pq f_pq E_pq ) | TE_tiaj > = 6 delta_ab ( delta_ki delta_lj - delta_kj delta_li ) ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )

         FEE[ It ][ ut ] = ( + 2 * delta_ut * f_dot_1dm
                             - f_dot_2dm[ It ][ tu ]
                             + 2 * f_ut
                             - sum_v Gamma_uv f_vt
                             - sum_w f_uw Gamma_wt
                           )

      G singlet: | SG_aibt > = ( E_ai E_bt + E_bi E_at ) | 0 >  with  a <= b
      G triplet: | TG_aibt > = ( E_ai E_bt - E_bi E_at ) | 0 >  with  a <  b

         < SG_cjdu | ( sum_pq f_pq E_pq ) | SG_aibt > = 2 delta_ji ( delta_ca delta_db + delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
         < TG_cjdu | ( sum_pq f_pq E_pq ) | TG_aibt > = 6 delta_ji ( delta_ca delta_db - delta_cb delta_da ) ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )

         FGG[ It ][ ut ] = f_dot_2dm[ It ][ ut ]

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



