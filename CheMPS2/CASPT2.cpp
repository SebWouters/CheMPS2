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
#include <iostream>
#include <math.h>
#include <assert.h>
#include <algorithm>
#include <sys/time.h>

#include "CASPT2.h"
#include "Lapack.h"
#include "Options.h"
#include "ConjugateGradient.h"
#include "Davidson.h"
#include "Special.h"

using std::cout;
using std::endl;
using std::min;
using std::max;

CheMPS2::CASPT2::CASPT2( DMRGSCFindices * idx, DMRGSCFintegrals * ints, DMRGSCFmatrix * oei, DMRGSCFmatrix * fock_in, double * one_dm, double * two_dm, double * three_dm, double * contract_4dm, const double IPEA ){

   indices    = idx;
   fock       = fock_in;
   one_rdm    = one_dm;
   two_rdm    = two_dm;
   three_rdm  = three_dm;
   f_dot_4dm  = contract_4dm;
   num_irreps = indices->getNirreps();

   struct timeval start, end;
   gettimeofday( &start, NULL );

   create_f_dots();
   vector_helper();

   make_AA_CC( true, 0.0 );
   make_DD( true, 0.0 );
   make_EE_GG( true, 0.0 );
   make_BB_FF_singlet( true, 0.0 );
   make_BB_FF_triplet( true, 0.0 );

   construct_rhs( oei, ints );

   make_AA_CC( false, IPEA );
   make_DD( false, IPEA );
   make_EE_GG( false, IPEA );
   make_BB_FF_singlet( false, IPEA );
   make_BB_FF_triplet( false, IPEA );
   make_FAD_FCD();
   make_FEH_FGH();
   make_FAB_FCF_singlet();
   make_FAB_FCF_triplet();
   make_FBE_FFG_singlet();
   make_FBE_FFG_triplet();
   make_FDE_FDG();

   delete [] f_dot_3dm;
   delete [] f_dot_2dm;

   gettimeofday( &end, NULL );
   double elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   cout << "CASPT2 : Wall time tensors    = " << elapsed << " seconds" << endl;
   gettimeofday( &start, NULL );

   recreate();

   gettimeofday( &end, NULL );
   elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   cout << "CASPT2 : Wall time diag(ovlp) = " << elapsed << " seconds" << endl;

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

double CheMPS2::CASPT2::solve( const double imag_shift, const bool CONJUGATE_GRADIENT ) const{

   struct timeval start, end;
   gettimeofday( &start, NULL );

   // Normalizations of sectors   A  Bs Bt C  D  Es Et Fs Ft Gs Gt Hs Ht
   const int normalizations[] = { 1, 2, 2, 1, 1, 2, 6, 2, 2, 2, 6, 4, 12 };
   const bool apply_shift = (( fabs( imag_shift ) > 0.0 ) ? true : false );

   int total_size = jump[ CHEMPS2_CASPT2_NUM_CASES * num_irreps ];
   double * diag_fock = new double[ total_size ];
   diagonal( diag_fock );
   double min_eig = diag_fock[ 0 ];
   for ( int elem = 1; elem < total_size; elem++ ){ min_eig = min( min_eig, diag_fock[ elem ] ); }
   cout << "CASPT2 : Solution algorithm   = " << (( CONJUGATE_GRADIENT ) ? "Conjugate Gradient" : "Davidson" ) << endl;
   cout << "CASPT2 : Minimum(diagonal)    = " << min_eig << endl;

   ConjugateGradient * CG = (( CONJUGATE_GRADIENT ) ? new ConjugateGradient( total_size, CheMPS2::CONJ_GRADIENT_RTOL, CheMPS2::CONJ_GRADIENT_PRECOND_CUTOFF, false ) : NULL );
   Davidson * DAVID = (( CONJUGATE_GRADIENT ) ? NULL : new Davidson( total_size,
                                                                     CheMPS2::DAVIDSON_NUM_VEC,
                                                                     CheMPS2::DAVIDSON_NUM_VEC_KEEP,
                                                                     CheMPS2::CONJ_GRADIENT_RTOL,
                                                                     CheMPS2::CONJ_GRADIENT_PRECOND_CUTOFF,
                                                                     true, // debug_print
                                                                     'L' )); // Linear problem
   double ** pointers = new double*[ 3 ];
   char instruction = (( CONJUGATE_GRADIENT ) ? CG->step( pointers ) : DAVID->FetchInstruction( pointers ));
   assert( instruction == 'A' );
   for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 0 ][ elem ] = vector_rhs[ elem ] / diag_fock[ elem ]; } // Initial guess of F * x = V
   for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 1 ][ elem ] =  diag_fock[ elem ]; } // Diagonal of the operator F
   for ( int elem = 0; elem < total_size; elem++ ){ pointers[ 2 ][ elem ] = vector_rhs[ elem ]; } // RHS of the linear problem F * x = V
   int inc1 = 1;
   const double E2_DIAGONAL = - ddot_( &total_size, pointers[ 0 ], &inc1, pointers[ 2 ], &inc1 );
   instruction = (( CONJUGATE_GRADIENT ) ? CG->step( pointers ) : DAVID->FetchInstruction( pointers ));
   assert( instruction == 'B' );
   while ( instruction == 'B' ){
      matvec( pointers[ 0 ], pointers[ 1 ], diag_fock );
      if ( apply_shift ){ add_shift( pointers[ 0 ], pointers[ 1 ], diag_fock, imag_shift, normalizations ); }
      instruction = (( CONJUGATE_GRADIENT ) ? CG->step( pointers ) : DAVID->FetchInstruction( pointers ));
   }
   assert( instruction == 'C' );
   const double E2_NONVARIATIONAL = - ddot_( &total_size, pointers[ 0 ], &inc1, vector_rhs, &inc1 );
   const double rnorm = pointers[ 1 ][ 0 ];
   cout << "CASPT2 : Number of iterations = " << (( CONJUGATE_GRADIENT ) ? CG->get_num_matvec() : DAVID->GetNumMultiplications() ) << endl;
   cout << "CASPT2 : Residual norm        = " << rnorm << endl;
   matvec( pointers[ 0 ], pointers[ 1 ], diag_fock ); // pointers[ 1 ] is a WORK array when instruction == 'C'
   const double E2_VARIATIONAL = 2 * E2_NONVARIATIONAL + ddot_( &total_size, pointers[ 0 ], &inc1, pointers[ 1 ], &inc1 );
   delete [] diag_fock;

   const double inproduct = inproduct_vectors( pointers[ 0 ], pointers[ 0 ], normalizations );
   const double reference_weight = 1.0 / ( 1.0 + inproduct );
   cout << "CASPT2 : Reference weight     = " << reference_weight << endl;
   //energy_per_sector( pointers[ 0 ] );
   delete [] pointers;
   if ( CG    != NULL ){ delete CG;    }
   if ( DAVID != NULL ){ delete DAVID; }

   gettimeofday( &end, NULL );
   double elapsed = ( end.tv_sec - start.tv_sec ) + 1e-6 * ( end.tv_usec - start.tv_usec );
   cout << "CASPT2 : Wall time solution   = " << elapsed << " seconds" << endl;

   cout << "CASPT2 : E2 [DIAGONAL]        = " << E2_DIAGONAL << endl;
   cout << "CASPT2 : E2 [NON-VARIATIONAL] = " << E2_NONVARIATIONAL << endl;
   cout << "CASPT2 : E2 [VARIATIONAL]     = " << E2_VARIATIONAL << endl;
   return E2_VARIATIONAL;

}

void CheMPS2::CASPT2::add_shift( double * vector, double * result, double * diag_fock, const double imag_shift, const int * normalizations ) const{

   for ( int sector = 0; sector < CHEMPS2_CASPT2_NUM_CASES; sector++ ){
      const int start = jump[ num_irreps * sector         ];
      const int stop  = jump[ num_irreps * ( sector + 1 ) ];
      const double factor = imag_shift * imag_shift * normalizations[ sector ] * normalizations[ sector ];
      for ( int elem = start; elem < stop; elem ++ ){
         result[ elem ] += factor * vector[ elem ] / diag_fock[ elem ];
      }
   }

}

double CheMPS2::CASPT2::inproduct_vectors( double * first, double * second, const int * normalizations ) const{

   int inc1 = 1;
   double value = 0.0;
   for ( int sector = 0; sector < CHEMPS2_CASPT2_NUM_CASES; sector++ ){
      int pointer = jump[ num_irreps * sector         ];
      int size    = jump[ num_irreps * ( sector + 1 ) ] - pointer;
      value += normalizations[ sector ] * ddot_( &size, first + pointer, &inc1, second + pointer, &inc1 );
   }
   return value;

}

void CheMPS2::CASPT2::energy_per_sector( double * solution ) const{

   int inc1 = 1;
   double energies[ CHEMPS2_CASPT2_NUM_CASES ];
   for ( int sector = 0; sector < CHEMPS2_CASPT2_NUM_CASES; sector++ ){
      int pointer = jump[ num_irreps * sector         ];
      int size    = jump[ num_irreps * ( sector + 1 ) ] - pointer;
      energies[ sector ] = - ddot_( &size, solution + pointer, &inc1, vector_rhs + pointer, &inc1 );
   }
   cout << "************************************************" << endl;
   cout << "*   CASPT2 non-variational energy per sector   *" << endl;
   cout << "************************************************" << endl;
   cout << "   A or VJTU = " << energies[ CHEMPS2_CASPT2_A ] << endl;
   cout << "   B or VJTI = " << energies[ CHEMPS2_CASPT2_B_SINGLET ] + energies[ CHEMPS2_CASPT2_B_TRIPLET ] << endl;
   cout << "   C or ATVX = " << energies[ CHEMPS2_CASPT2_C ] << endl;
   cout << "   D or AIVX = " << energies[ CHEMPS2_CASPT2_D ] << endl;
   cout << "   E or VJAI = " << energies[ CHEMPS2_CASPT2_E_SINGLET ] + energies[ CHEMPS2_CASPT2_E_TRIPLET ] << endl;
   cout << "   F or BVAT = " << energies[ CHEMPS2_CASPT2_F_SINGLET ] + energies[ CHEMPS2_CASPT2_F_TRIPLET ] << endl;
   cout << "   G or BJAT = " << energies[ CHEMPS2_CASPT2_G_SINGLET ] + energies[ CHEMPS2_CASPT2_G_TRIPLET ] << endl;
   cout << "   H or BJAI = " << energies[ CHEMPS2_CASPT2_H_SINGLET ] + energies[ CHEMPS2_CASPT2_H_TRIPLET ] << endl;
   cout << "************************************************" << endl;

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
      for ( int orb = 0; orb < NDMRG; orb++ ){
         value += fock->get( irrep, NOCC + orb, NOCC + orb ) * one_rdm[ jumpx + orb + LAS * ( jumpx + orb ) ];
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

               for ( int orbx = 0; orbx < NDMRGx; orbx++ ){
                  value += fock->get( irrepx, NOCCx + orbx, NOCCx + orbx ) * two_rdm[ i1 + LAS * ( jumpx + orbx + LAS * ( i2 + LAS * ( jumpx + orbx ))) ];
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

                           for ( int orbx = 0; orbx < NDMRGx; orbx++ ){
                              value += ( fock->get( irrepx, NOCCx + orbx, NOCCx + orbx )
                                       * three_rdm[ i1 + LAS*( i2 + LAS*( jumpx + orbx + LAS*( i3 + LAS*( i4 + LAS*( jumpx + orbx ))))) ] );
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

   /*
   double sum_f_kk = 0.0;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      const int NOCC = indices->getNOCC( irrep );
      for ( int orb = 0; orb < NOCC; orb++ ){
         sum_f_kk += fock->get( irrep, orb, orb );
      }
   }

   const double E_FOCK = 2 * sum_f_kk + f_dot_1dm;
   cout << "CASPT2 : < 0 | F | 0 >        = " << E_FOCK << endl;
   */

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
   assert( total_size == vector_length( indices ) );
   cout << "CASPT2 : Old size V_SD space  = " << total_size << endl;
   return total_size;

}

long long CheMPS2::CASPT2::vector_length( const DMRGSCFindices * idx ){

   const int NUM_IRREPS = idx->getNirreps();
   long long length = 0;
   for ( int i1 = 0; i1 < NUM_IRREPS; i1++ ){
      const long long nocc1 = idx->getNOCC( i1 );
      const long long nact1 = idx->getNDMRG( i1 );
      const long long nvir1 = idx->getNVIRT( i1 );
      for ( int i2 = 0; i2 < NUM_IRREPS; i2++ ){
         const long long nocc2 = idx->getNOCC( i2 );
         const long long nact2 = idx->getNDMRG( i2 );
         for ( int i3 = 0; i3 < NUM_IRREPS; i3++ ){
            const int i4 = Irreps::directProd( Irreps::directProd( i1, i2 ), i3 );
            const long long nact3 = idx->getNDMRG( i3 );
            const long long nvir3 = idx->getNVIRT( i3 );
            const long long nocc4 = idx->getNOCC( i4 );
            const long long nact4 = idx->getNDMRG( i4 );
            const long long nvir4 = idx->getNVIRT( i4 );

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

   if ( NEWSIZE == 0 ){ return NEWSIZE; }

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
   double * work  = new double[ lwork ];
   double * eigs  = new double[ maxsize ];

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

   double * temp = NULL;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
       temp = new double[ size_A[ irrep ] ]; dcopy_( size_A + irrep, FAA[ irrep ], &inc1, temp, &inc1 ); delete [] FAA[ irrep ]; FAA[ irrep ] = temp;
       temp = new double[ size_C[ irrep ] ]; dcopy_( size_C + irrep, FCC[ irrep ], &inc1, temp, &inc1 ); delete [] FCC[ irrep ]; FCC[ irrep ] = temp;
       temp = new double[ size_D[ irrep ] ]; dcopy_( size_D + irrep, FDD[ irrep ], &inc1, temp, &inc1 ); delete [] FDD[ irrep ]; FDD[ irrep ] = temp;
       temp = new double[ size_E[ irrep ] ]; dcopy_( size_E + irrep, FEE[ irrep ], &inc1, temp, &inc1 ); delete [] FEE[ irrep ]; FEE[ irrep ] = temp;
       temp = new double[ size_G[ irrep ] ]; dcopy_( size_G + irrep, FGG[ irrep ], &inc1, temp, &inc1 ); delete [] FGG[ irrep ]; FGG[ irrep ] = temp;
       temp = new double[ size_B_singlet[ irrep ] ]; dcopy_( size_B_singlet + irrep, FBB_singlet[ irrep ], &inc1, temp, &inc1 ); delete [] FBB_singlet[ irrep ]; FBB_singlet[ irrep ] = temp;
       temp = new double[ size_B_triplet[ irrep ] ]; dcopy_( size_B_triplet + irrep, FBB_triplet[ irrep ], &inc1, temp, &inc1 ); delete [] FBB_triplet[ irrep ]; FBB_triplet[ irrep ] = temp;
       temp = new double[ size_F_singlet[ irrep ] ]; dcopy_( size_F_singlet + irrep, FFF_singlet[ irrep ], &inc1, temp, &inc1 ); delete [] FFF_singlet[ irrep ]; FFF_singlet[ irrep ] = temp;
       temp = new double[ size_F_triplet[ irrep ] ]; dcopy_( size_F_triplet + irrep, FFF_triplet[ irrep ], &inc1, temp, &inc1 ); delete [] FFF_triplet[ irrep ]; FFF_triplet[ irrep ] = temp;
   }

   const int total_size = jump[ num_irreps * CHEMPS2_CASPT2_NUM_CASES ];
   cout << "CASPT2 : New size V_SD space  = " << total_size << endl;

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

void CheMPS2::CASPT2::matmat( char totrans, int rowdim, int coldim, int sumdim, double alpha, double * matrix, int ldaM, double * origin, int ldaO, double * target, int ldaT ){

   double add = 1.0;
   char notrans = 'N';
   dgemm_( &totrans, &notrans, &rowdim, &coldim, &sumdim, &alpha, matrix, &ldaM, origin, &ldaO, &add, target, &ldaT );

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
   #pragma omp simd
   for ( int elem = 0; elem < vectorlength; elem++ ){ result[ elem ] = diag_fock[ elem ] * vector[ elem ]; }
   const int maxlinsize = get_maxsize();
   double * workspace = new double[ maxlinsize * maxlinsize ];
   const double SQRT2 = sqrt( 2.0 );

   // FAD: < A(xjyz) E_wc D(aitu) > = delta_ac delta_ij FAD[ Ij ][ Ii x Ia ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ii == Ij == Ix x Iy x Iz
      const int SIZE_L = size_A[ IL ];
      const int nocc_ij = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ia
         const int SIZE_R = size_D[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + nocc_ij * ac );
               matmat( 'N', SIZE_L, nocc_ij, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, nocc_ij, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
            }
         }
      }
   }

   // FCD: < C(bxyz) E_kw D(aitu) > = delta_ik delta_ab FCD[ Ib ][ Ii x Ia ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ia == Ib
      const int SIZE_L = size_C[ IL ];
      const int nvir_ab = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ia
         const int SIZE_R = size_D[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_R * ( shift + ik );
               const int LDA_R = SIZE_R * nocc_w;
               matmat( 'N', SIZE_L, nvir_ab, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, nvir_ab, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
            }
         }
      }
   }

   // FAB singlet: < A(xlyz) E_kw SB_tiuj > = ( delta_ik delta_jl + delta_jk delta_il ) / sqrt( 1 + delta_ij ) * FAB_singlet[ Il ][ Ii x Ij ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Il == Ix x Iy x Iz
      const int SIZE_L = size_A[ IL ];
      const int nocc_l = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ij
         const int SIZE_R = size_B_singlet[ IR ];
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
                  if ( k > 0 ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ];
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * ( shift + ( k * ( k + 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, k, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, k, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
                  #pragma omp parallel for schedule(static)
                  for ( int l = k; l < nocc_l; l++ ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * ( shift + k + ( l * ( l + 1 ) ) / 2 );
                     const double factor = (( k == l ) ? SQRT2 : 1.0 );
                     matmat( 'N', SIZE_L, 1, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, 1, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
               } else {
                  const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ];
                  const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE_R * ( shift + (( Iw < IL ) ? k : nocc_l * k ));
                  const int LDA_R = (( Iw < IL ) ? SIZE_R * nocc_w : SIZE_R );
                  matmat( 'N', SIZE_L, nocc_l, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
                  matmat( 'T', SIZE_R, nocc_l, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
               }
            }
         }
      }
   }

   // FAB triplet: < A(xlyz) E_kw TB_tiuj > = ( delta_ik delta_jl - delta_jk delta_il ) * FAB_triplet[ Il ][ Ii x Ij ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Il == Ix x Iy x Iz
      const int SIZE_L = size_A[ IL ];
      const int nocc_l = indices->getNOCC( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ii x Ij
         const int SIZE_R = size_B_triplet[ IR ];
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
                  if ( k > 0 ){ // ( k > l  --->  - delta_jk delta_il )
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ];
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * ( shift + ( k * ( k - 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, k, SIZE_R, -1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, k, SIZE_L, -1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
                  #pragma omp parallel for schedule(static)
                  for ( int l = k+1; l < nocc_l; l++ ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ] + SIZE_L * l;
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * ( shift + k + ( l * ( l - 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, 1, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L ); // ( k < l  --->  + delta_ik delta_jl )
                     matmat( 'T', SIZE_R, 1, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
               } else {
                  const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_A         ];
                  const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE_R * ( shift + (( Iw < IL ) ? k : nocc_l * k ));
                  const int LDA_R = (( Iw < IL ) ? SIZE_R * nocc_w : SIZE_R );
                  const double factor = (( Iw < IL ) ? 1.0 : -1.0 ); // ( k < l  --->  + delta_ik delta_jl ) and ( k > l  --->  - delta_jk delta_il )
                  matmat( 'N', SIZE_L, nocc_l, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
                  matmat( 'T', SIZE_R, nocc_l, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
               }
            }
         }
      }
   }

   // FCF singlet: < C(dxyz) E_wc SF_atbu > = ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) * FCF_singlet[ Id ][ Ia x Ib ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Id == Ix x Iy x Iz
      const int SIZE_L = size_C[ IL ];
      const int nvir_d = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ia x Ib
         const int SIZE_R = size_F_singlet[ IR ];
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
                  if ( c > 0 ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ];
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * ( shift + ( c * ( c + 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, c, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, c, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
                  #pragma omp parallel for schedule(static)
                  for ( int d = c; d < nvir_d; d++ ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * ( shift + c + ( d * ( d + 1 ) ) / 2 );
                     const double factor = (( c == d ) ? SQRT2 : 1.0 );
                     matmat( 'N', SIZE_L, 1, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, 1, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
               } else {
                  const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ];
                  const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE_R * ( shift + (( Iw < IL ) ? c : nvir_d * c ));
                  const int LDA_R = (( Iw < IL ) ? SIZE_R * nvir_w : SIZE_R );
                  matmat( 'N', SIZE_L, nvir_d, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
                  matmat( 'T', SIZE_R, nvir_d, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
               }
            }
         }
      }
   }

   // FCF triplet: < C(dxyz) E_wc TF_atbu > = ( delta_ac delta_bd - delta_ad delta_bc ) * FCF_triplet[ Id ][ Ia x Ib ][ w ][ (xyz),(tu) ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Id == Ix x Iy x Iz
      const int SIZE_L = size_C[ IL ];
      const int nvir_d = indices->getNVIRT( IL );
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It x Iu == Ia x Ib
         const int SIZE_R = size_F_triplet[ IR ];
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
                  if ( c > 0 ){ // ( c > d  --->  - delta_ad delta_bc )
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ];
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * ( shift + ( c * ( c - 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, c, SIZE_R, -1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                     matmat( 'T', SIZE_R, c, SIZE_L, -1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
                  #pragma omp parallel for schedule(static)
                  for ( int d = c+1; d < nvir_d; d++ ){
                     const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ] + SIZE_L * d;
                     const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * ( shift + c + ( d * ( d - 1 ) ) / 2 );
                     matmat( 'N', SIZE_L, 1, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L ); // ( c < d  --->  + delta_ac delta_bd )
                     matmat( 'T', SIZE_R, 1, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                  }
               } else {
                  const double factor = (( Iw < IL ) ? 1.0 : -1.0 ); // ( c < d  --->  + delta_ac delta_bd ) and ( c > d  --->  - delta_ad delta_bc )
                  const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_C         ];
                  const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE_R * ( shift + (( Iw < IL ) ? c : nvir_d * c ));
                  const int LDA_R = (( Iw < IL ) ? SIZE_R * nvir_w : SIZE_R );
                  matmat( 'N', SIZE_L, nvir_d, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
                  matmat( 'T', SIZE_R, nvir_d, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
               }
            }
         }
      }
   }

   // FBE singlet: < SB_xkyl E_wc SE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_singlet[ Ik x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iik x Ijl == Ix x Iy
      const int SIZE_L = size_B_singlet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iik x Ijl x Iac
         const int SIZE_R = size_E[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_B_SINGLET ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE_R * ( shift_E + ac );
               const int LDA_R = SIZE_R * nvir_w;
               matmat( 'N', SIZE_L, size_ij, SIZE_R, 2.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, size_ij, SIZE_L, 2.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
            }
         }
      }
   }

   // FBE triplet: < TB_xkyl E_wc TE_tiaj > = 2 delta_ac delta_ik delta_jl FBE_triplet[ Ik x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iik x Ijl == Ix x Iy
      const int SIZE_L = size_B_triplet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iik x Ijl x Iac
         const int SIZE_R = size_E[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE_R * ( shift_E + ac );
               const int LDA_R = SIZE_R * nvir_w;
               matmat( 'N', SIZE_L, size_ij, SIZE_R, 2.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, size_ij, SIZE_L, 2.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
            }
         }
      }
   }

   // FFG singlet: < SF_cxdy E_kw SG_aibt > = 2 delta_ac delta_bd delta_ik FFG_singlet[ Ic x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iac x Ibd == Ix x Iy
      const int SIZE_L = size_F_singlet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iac x Ibd x Iik
         const int SIZE_R = size_G[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_F_SINGLET ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE_R * ( shift_G + ik );
               const int LDA_R = SIZE_R * nocc_w;
               matmat( 'N', SIZE_L, size_ab, SIZE_R, 2.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, size_ab, SIZE_L, 2.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
            }
         }
      }
   }

   // FFG triplet: < TF_cxdy E_kw TG_aibt > = 2 delta_ac delta_bd delta_ik FFG_triplet[ Ic x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Iac x Ibd == Ix x Iy
      const int SIZE_L = size_F_triplet[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR == It == Iac x Ibd x Iik
         const int SIZE_R = size_G[ IR ];
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
               const int ptr_L = jump[ IL + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ];
               const int ptr_R = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE_R * ( shift_G + ik );
               const int LDA_R = SIZE_R * nocc_w;
               matmat( 'N', SIZE_L, size_ab, SIZE_R, 2.0, workspace, SIZE_L, vector + ptr_R, LDA_R,  result + ptr_L, SIZE_L );
               matmat( 'T', SIZE_R, size_ab, SIZE_L, 2.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, LDA_R  );
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
      if ( SIZE * nact_w * nvir_w > 0 ){
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
               const int LDA_E = SIZE * nvir_d;

               if ( IL < Id ){ // Ic < Id
                  for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                     const int Ij = Irreps::directProd( Ii, Icenter );
                     if ( Ii < Ij ){
                        const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, +1 );
                        const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Ii, Ij, IL, Id, +1 );
                        const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                        #pragma omp parallel for schedule(static)
                        for ( int d = 0; d < nvir_d; d++ ){
                           const int ptr_H = jump_H + size_ij * ( c + nvir_w * d );
                           const int ptr_E = jump_E + SIZE * d;
                           matmat( 'N', SIZE, size_ij, 1,    2.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                           matmat( 'T', 1,    size_ij, SIZE, 2.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
                        }
                     }
                  }
               }

               if ( IL == Id ){ // Ic == Id == Ia == Ib --> Iik == Ijl
                  for ( int Iij = 0; Iij < num_irreps; Iij++ ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * shift_E_nonactive( indices, Id,  Iij, Iij, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Iij, Iij, IL,  Id, +1 );
                     const int size_ij = ( indices->getNOCC( Iij ) * ( indices->getNOCC( Iij ) + 1 ) ) / 2;
                     #pragma omp parallel for schedule(static)
                     for ( int d = 0; d < c; d++ ){
                        const int ptr_H = jump_H + size_ij * ( d + ( c * ( c + 1 ) ) / 2 );
                        const int ptr_E = jump_E + SIZE * d;
                        matmat( 'N', SIZE, size_ij, 1,    2.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                        matmat( 'T', 1,    size_ij, SIZE, 2.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
                     }
                     #pragma omp parallel for schedule(static)
                     for ( int d = c; d < nvir_d; d++ ){
                        const int ptr_H = jump_H + size_ij * ( c + ( d * ( d + 1 ) ) / 2 );
                        const int ptr_E = jump_E + SIZE * d;
                        const double factor = 2 * (( c == d ) ? SQRT2 : 1.0 );
                        matmat( 'N', SIZE, size_ij, 1,    factor, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                        matmat( 'T', 1,    size_ij, SIZE, factor, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
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
                        #pragma omp parallel for schedule(static)
                        for ( int d = 0; d < nvir_d; d++ ){
                           const int ptr_H = jump_H + size_ij * ( d + nvir_d * c );
                           const int ptr_E = jump_E + SIZE * d;
                           matmat( 'N', SIZE, size_ij, 1,    2.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                           matmat( 'T', 1,    size_ij, SIZE, 2.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
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
      if ( SIZE * nact_w * nvir_w > 0 ){
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
               const int LDA_E = SIZE * nvir_d;

               if ( IL < Id ){ // Ic < Id
                  for ( int Ii = 0; Ii < num_irreps; Ii++ ){
                     const int Ij = Irreps::directProd( Ii, Icenter );
                     if ( Ii < Ij ){
                        const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * shift_E_nonactive( indices, Id, Ii, Ij, -1 );
                        const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Ii, Ij, IL, Id, -1 );
                        const int size_ij = indices->getNOCC( Ii ) * indices->getNOCC( Ij );
                        #pragma omp parallel for schedule(static)
                        for ( int d = 0; d < nvir_d; d++ ){
                           const int ptr_H = jump_H + size_ij * ( c + nvir_w * d );
                           const int ptr_E = jump_E + SIZE * d;
                           matmat( 'N', SIZE, size_ij, 1,    6.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                           matmat( 'T', 1,    size_ij, SIZE, 6.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
                        }
                     }
                  }
               }

               if ( IL == Id ){ // Ic == Id == Ia == Ib --> Iik == Ijl
                  for ( int Iij = 0; Iij < num_irreps; Iij++ ){
                     const int jump_E = jump[ IL + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * shift_E_nonactive( indices, Id,  Iij, Iij, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Iij, Iij, IL,  Id, -1 );
                     const int size_ij = ( indices->getNOCC( Iij ) * ( indices->getNOCC( Iij ) - 1 ) ) / 2;
                     #pragma omp parallel for schedule(static)
                     for ( int d = 0; d < c; d++ ){
                        const int ptr_H = jump_H + size_ij * ( d + ( c * ( c - 1 ) ) / 2 );
                        const int ptr_E = jump_E + SIZE * d;
                        matmat( 'N', SIZE, size_ij, 1,    -6.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                        matmat( 'T', 1,    size_ij, SIZE, -6.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
                     }
                     #pragma omp parallel for schedule(static)
                     for ( int d = c+1; d < nvir_d; d++ ){
                        const int ptr_H = jump_H + size_ij * ( c + ( d * ( d - 1 ) ) / 2 );
                        const int ptr_E = jump_E + SIZE * d;
                        matmat( 'N', SIZE, size_ij, 1,    6.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                        matmat( 'T', 1,    size_ij, SIZE, 6.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
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
                        #pragma omp parallel for schedule(static)
                        for ( int d = 0; d < nvir_d; d++ ){
                           const int ptr_H = jump_H + size_ij * ( d + nvir_d * c );
                           const int ptr_E = jump_E + SIZE * d;
                           matmat( 'N', SIZE, size_ij, 1,    -6.0, workspace, SIZE, vector + ptr_H, 1,     result + ptr_E, LDA_E );
                           matmat( 'T', 1,    size_ij, SIZE, -6.0, workspace, SIZE, vector + ptr_E, LDA_E, result + ptr_H, 1     );
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
      if ( SIZE * nocc_w * nact_w > 0 ){
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
                        const int colsize = nocc_l * indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                        if ( colsize > 0 ){
                           const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, +1 );
                           const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, IL, Il, Ia, Ib, +1 ) + k;
                           matmat( 'N', SIZE, colsize, 1,    2.0, workspace, SIZE, vector + jump_H, nocc_w, result + jump_G, SIZE   );
                           matmat( 'T', 1,    colsize, SIZE, 2.0, workspace, SIZE, vector + jump_G, SIZE,   result + jump_H, nocc_w );
                        }
                     }
                  }
               }

               if ( IL == Il ){ // Ik == Il == Ii == Ij --> Iac == Ibd   and   nocc_l == nocc_w
                  for ( int Iab = 0; Iab < num_irreps; Iab++ ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Iab, Iab, +1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, IL, Il, Iab, Iab, +1 );
                     const int size_ij = ( nocc_w * ( nocc_w + 1 ) ) / 2;
                     const int size_ab = ( indices->getNVIRT( Iab ) * ( indices->getNVIRT( Iab ) + 1 ) ) / 2;
                     const int LDA_G = SIZE * nocc_l;
                     #pragma omp parallel for schedule(static)
                     for ( int l = 0; l < k; l++ ){
                        const int ptr_H = jump_H + ( l + ( k * ( k + 1 ) ) / 2 );
                        const int ptr_G = jump_G + SIZE * l;
                        matmat( 'N', SIZE, size_ab, 1,    2.0, workspace, SIZE, vector + ptr_H, size_ij, result + ptr_G, LDA_G   );
                        matmat( 'T', 1,    size_ab, SIZE, 2.0, workspace, SIZE, vector + ptr_G, LDA_G,   result + ptr_H, size_ij );
                     }
                     #pragma omp parallel for schedule(static)
                     for ( int l = k; l < nocc_l; l++ ){
                        const int ptr_H = jump_H + ( k + ( l * ( l + 1 ) ) / 2 );
                        const int ptr_G = jump_G + SIZE * l;
                        const double factor = 2 * (( k == l ) ? SQRT2 : 1.0 );
                        matmat( 'N', SIZE, size_ab, 1,    factor, workspace, SIZE, vector + ptr_H, size_ij, result + ptr_G, LDA_G   );
                        matmat( 'T', 1,    size_ab, SIZE, factor, workspace, SIZE, vector + ptr_G, LDA_G,   result + ptr_H, size_ij );
                     }
                  }
               }

               if ( IL > Il ){ // Ik > Il
                  for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                     const int Ib = Irreps::directProd( Ia, Icenter );
                     if ( Ia < Ib ){
                        const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, +1 );
                        const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift_H_nonactive( indices, Il, IL, Ia, Ib, +1 );
                        const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                        #pragma omp parallel for schedule(static)
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const int ptr_H = jump_H + nocc_l * ( k + nocc_w * ab );
                           const int ptr_G = jump_G + SIZE * nocc_l * ab;
                           matmat( 'N', SIZE, nocc_l, 1,    2.0, workspace, SIZE, vector + ptr_H, 1,    result + ptr_G, SIZE );
                           matmat( 'T', 1,    nocc_l, SIZE, 2.0, workspace, SIZE, vector + ptr_G, SIZE, result + ptr_H, 1    );
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
      if ( SIZE * nocc_w * nact_w > 0 ){
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
                        const int colsize = nocc_l * indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                        if ( colsize > 0 ){
                           const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, -1 );
                           const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, IL, Il, Ia, Ib, -1 ) + k;
                           matmat( 'N', SIZE, colsize, 1,    -6.0, workspace, SIZE, vector + jump_H, nocc_w, result + jump_G, SIZE   );
                           matmat( 'T', 1,    colsize, SIZE, -6.0, workspace, SIZE, vector + jump_G, SIZE,   result + jump_H, nocc_w );
                        }
                     }
                  }
               }

               if ( IL == Il ){ // Ik == Il == Ii == Ij --> Iac == Ibd   and   nocc_l == nocc_w
                  for ( int Iab = 0; Iab < num_irreps; Iab++ ){
                     const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Iab, Iab, -1 );
                     const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, IL, Il, Iab, Iab, -1 );
                     const int size_ij = ( nocc_w * ( nocc_w - 1 ) ) / 2;
                     const int size_ab = ( indices->getNVIRT( Iab ) * ( indices->getNVIRT( Iab ) - 1 ) ) / 2;
                     const int LDA_G = SIZE * nocc_l;
                     #pragma omp parallel for schedule(static)
                     for ( int l = 0; l < k; l++ ){
                        const int ptr_H = jump_H + ( l + ( k * ( k - 1 ) ) / 2 );
                        const int ptr_G = jump_G + SIZE * l;
                        matmat( 'N', SIZE, size_ab, 1,    6.0, workspace, SIZE, vector + ptr_H, size_ij, result + ptr_G, LDA_G   );
                        matmat( 'T', 1,    size_ab, SIZE, 6.0, workspace, SIZE, vector + ptr_G, LDA_G,   result + ptr_H, size_ij );
                     }
                     #pragma omp parallel for schedule(static)
                     for ( int l = k+1; l < nocc_l; l++ ){
                        const int ptr_H = jump_H + ( k + ( l * ( l - 1 ) ) / 2 );
                        const int ptr_G = jump_G + SIZE * l;
                        matmat( 'N', SIZE, size_ab, 1,    -6.0, workspace, SIZE, vector + ptr_H, size_ij, result + ptr_G, LDA_G   );
                        matmat( 'T', 1,    size_ab, SIZE, -6.0, workspace, SIZE, vector + ptr_G, LDA_G,   result + ptr_H, size_ij );
                     }
                  }
               }

               if ( IL > Il ){ // Ik > Il
                  for ( int Ia = 0; Ia < num_irreps; Ia++ ){
                     const int Ib = Irreps::directProd( Ia, Icenter );
                     if ( Ia < Ib ){
                        const int jump_G = jump[ IL + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * shift_G_nonactive( indices, Il, Ia, Ib, -1 );
                        const int jump_H = jump[ Icenter + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift_H_nonactive( indices, Il, IL, Ia, Ib, -1 );
                        const int size_ab = indices->getNVIRT( Ia ) * indices->getNVIRT( Ib );
                        #pragma omp parallel for schedule(static)
                        for ( int ab = 0; ab < size_ab; ab++ ){
                           const int ptr_H = jump_H + nocc_l * ( k + nocc_w * ab );
                           const int ptr_G = jump_G + SIZE * nocc_l * ab;
                           matmat( 'N', SIZE, nocc_l, 1,    6.0, workspace, SIZE, vector + ptr_H, 1,    result + ptr_G, SIZE );
                           matmat( 'T', 1,    nocc_l, SIZE, 6.0, workspace, SIZE, vector + ptr_G, SIZE, result + ptr_H, 1    );
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
      const int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ii x Ij
         const int SIZE_R = size_E[ IR ];
         const int Ikw = Irreps::directProd( IL, IR ); // Ikw == Ik == Iw
         const int nocc_kw = indices->getNOCC( Ikw );
         const int nact_kw = indices->getNDMRG( Ikw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nocc_kw * nact_kw > 0 ){
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
                  if ( nvir_ab * nocc_l > 0 ){
                     const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D         ] + SIZE_L * shift_D_nonactive( indices, Il, Iab );
                     const int jump_E = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE_R * (( Ikw <= Il ) ? shift_E_nonactive( indices, Iab, Ikw, Il,  +1 )
                                                                                                                     : shift_E_nonactive( indices, Iab, Il,  Ikw, +1 ));
                     if ( Ikw == Il ){ // irrep_k == irrep_l
                        #pragma omp parallel for schedule(static)
                        for ( int l = 0; l < nocc_l; l++ ){
                           const double factor = (( k == l ) ? SQRT2 : 1.0 );
                           const int ptr_L = jump_D + SIZE_L * l;
                           const int ptr_R = jump_E + SIZE_R * nvir_ab * (( k < l ) ? ( k + ( l * ( l + 1 ) ) / 2 ) : ( l + ( k * ( k + 1 ) ) / 2 ));
                           const int LDA_L = SIZE_L * nocc_l;
                           matmat( 'N', SIZE_L, nvir_ab, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, LDA_L  );
                           matmat( 'T', SIZE_R, nvir_ab, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, LDA_L,  result + ptr_R, SIZE_R );
                        }
                     } else { // irrep_k != irrep_l
                        #pragma omp parallel for schedule(static)
                        for ( int l = 0; l < nocc_l; l++ ){
                           const int ptr_L = jump_D + SIZE_L * l;
                           const int ptr_R = jump_E + SIZE_R * nvir_ab * (( Ikw < Il ) ? ( k + nocc_kw * l ) : ( l + nocc_l * k ));
                           const int LDA_L = SIZE_L * nocc_l;
                           matmat( 'N', SIZE_L, nvir_ab, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, LDA_L  );
                           matmat( 'T', SIZE_R, nvir_ab, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, LDA_L,  result + ptr_R, SIZE_R );
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FDE triplet: < D(blxy) E_kw TE_tiaj > = 3 delta_ab ( delta_ik delta_jl - delta_il delta_jk ) / sqrt( 1 + delta_ij ) FDE_triplet[ Ib x Il ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Ib x Il
      const int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ii x Ij
         const int SIZE_R = size_E[ IR ];
         const int Ikw = Irreps::directProd( IL, IR ); // Ikw == Ik == Iw
         const int nocc_kw = indices->getNOCC( Ikw );
         const int nact_kw = indices->getNDMRG( Ikw );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nocc_kw * nact_kw > 0 ){
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
                  if ( nvir_ab * nocc_l > 0 ){
                     const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D         ] + SIZE_L * shift_D_nonactive( indices, Il, Iab );
                     const int jump_E = jump[ IR + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE_R * (( Ikw <= Il ) ? shift_E_nonactive( indices, Iab, Ikw, Il,  -1 )
                                                                                                                     : shift_E_nonactive( indices, Iab, Il,  Ikw, -1 ));
                     if ( Ikw == Il ){ // irrep_k == irrep_l
                        #pragma omp parallel for schedule(static)
                        for ( int l = 0; l < k; l++ ){
                           const int ptr_L = jump_D + SIZE_L * l;
                           const int ptr_R = jump_E + SIZE_R * nvir_ab * ( l + ( k * ( k - 1 ) ) / 2 );
                           const int LDA_L = SIZE_L * nocc_l;
                           matmat( 'N', SIZE_L, nvir_ab, SIZE_R, -3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, LDA_L  );
                           matmat( 'T', SIZE_R, nvir_ab, SIZE_L, -3.0, workspace, SIZE_L, vector + ptr_L, LDA_L,  result + ptr_R, SIZE_R );
                        }
                        #pragma omp parallel for schedule(static)
                        for ( int l = k+1; l < nocc_l; l++ ){
                           const int ptr_L = jump_D + SIZE_L * l;
                           const int ptr_R = jump_E + SIZE_R * nvir_ab * ( k + ( l * ( l - 1 ) ) / 2 );
                           const int LDA_L = SIZE_L * nocc_l;
                           matmat( 'N', SIZE_L, nvir_ab, SIZE_R, 3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, LDA_L  );
                           matmat( 'T', SIZE_R, nvir_ab, SIZE_L, 3.0, workspace, SIZE_L, vector + ptr_L, LDA_L,  result + ptr_R, SIZE_R );
                        }
                     } else { // irrep_k != irrep_l
                        #pragma omp parallel for schedule(static)
                        for ( int l = 0; l < nocc_l; l++ ){
                           const double factor = (( Ikw < Il ) ? 3.0 : -3.0 );
                           const int ptr_L = jump_D + SIZE_L * l;
                           const int ptr_R = jump_E + SIZE_R * nvir_ab * (( Ikw < Il ) ? ( k + nocc_kw * l ) : ( l + nocc_l * k ));
                           const int LDA_L = SIZE_L * nocc_l;
                           matmat( 'N', SIZE_L, nvir_ab, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, LDA_L  );
                           matmat( 'T', SIZE_R, nvir_ab, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, LDA_L,  result + ptr_R, SIZE_R );
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FDG singlet: < D(djxy) E_wc SG_aibt > = 1 delta_ij ( delta_ac delta_bd + delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FDG_singlet[ Ij x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Id x Ij
      const int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ib x Ii
         const int SIZE_R = size_E[ IR ];
         const int Iwc = Irreps::directProd( IL, IR ); // Iwc == Ic == Iw
         const int nocc_wc = indices->getNOCC( Iwc );
         const int nact_wc = indices->getNDMRG( Iwc );
         const int n_oa_wc = nocc_wc + nact_wc;
         const int nvir_wc = indices->getNVIRT( Iwc );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nvir_wc * nact_wc > 0 ){
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
                  if ( nvir_d * nocc_ij > 0 ){
                     const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_L * shift_D_nonactive( indices, Iij, Id );
                     const int jump_G = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE_R * (( Iwc <= Id ) ? shift_G_nonactive( indices, Iij, Iwc, Id,  +1 )
                                                                                                                     : shift_G_nonactive( indices, Iij, Id,  Iwc, +1 ));
                     if ( Iwc == Id ){ // irrep_c == irrep_d
                        if ( c > 0 ){
                           const int ptr_L = jump_D;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c * ( c + 1 ) ) / 2;
                           matmat( 'N', SIZE_L, c * nocc_ij, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, c * nocc_ij, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
                        #pragma omp parallel for schedule(static)
                        for ( int d = c; d < nvir_d; d++ ){
                           const double factor = (( c == d ) ? SQRT2 : 1.0 );
                           const int ptr_L = jump_D + SIZE_L * nocc_ij * d;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c + ( d * ( d + 1 ) ) / 2 );
                           matmat( 'N', SIZE_L, nocc_ij, SIZE_R, factor, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, nocc_ij, SIZE_L, factor, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
                     } else { // irrep_c != irrep_d
                        if ( Iwc < Id ){
                           #pragma omp parallel for schedule(static)
                           for ( int d = 0; d < nvir_d; d++ ){
                              const int ptr_L = jump_D + SIZE_L * nocc_ij * d;
                              const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c + nvir_wc * d );
                              matmat( 'N', SIZE_L, nocc_ij, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                              matmat( 'T', SIZE_R, nocc_ij, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                           }
                        } else {
                           const int ptr_L = jump_D;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * nvir_d * c;
                           matmat( 'N', SIZE_L, nvir_d * nocc_ij, SIZE_R, 1.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, nvir_d * nocc_ij, SIZE_L, 1.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // FDG triplet: < D(djxy) E_wc TG_aibt > = 3 delta_ij ( delta_ac delta_bd - delta_ad delta_bc ) / sqrt( 1 + delta_ab ) FDG_triplet[ Ij x Id ][ It ][ w ][ xy, t ]
   for ( int IL = 0; IL < num_irreps; IL++ ){ // IL == Ix x Iy == Id x Ij
      const int SIZE_L = size_D[ IL ];
      for ( int IR = 0; IR < num_irreps; IR++ ){ // IR = It = Ia x Ib x Ii
         const int SIZE_R = size_E[ IR ];
         const int Iwc = Irreps::directProd( IL, IR ); // Iwc == Ic == Iw
         const int nocc_wc = indices->getNOCC( Iwc );
         const int nact_wc = indices->getNDMRG( Iwc );
         const int n_oa_wc = nocc_wc + nact_wc;
         const int nvir_wc = indices->getNVIRT( Iwc );
         int total_size = SIZE_L * SIZE_R;
         if ( total_size * nvir_wc * nact_wc > 0 ){
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
                  if ( nvir_d * nocc_ij > 0 ){
                     const int jump_D = jump[ IL + num_irreps * CHEMPS2_CASPT2_D ] + SIZE_L * shift_D_nonactive( indices, Iij, Id );
                     const int jump_G = jump[ IR + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE_R * (( Iwc <= Id ) ? shift_G_nonactive( indices, Iij, Iwc, Id,  -1 )
                                                                                                                     : shift_G_nonactive( indices, Iij, Id,  Iwc, -1 ));
                     if ( Iwc == Id ){ // irrep_c == irrep_d
                        if ( c > 0 ){
                           const int ptr_L = jump_D;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c * ( c - 1 ) ) / 2;
                           matmat( 'N', SIZE_L, c * nocc_ij, SIZE_R, -3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, c * nocc_ij, SIZE_L, -3.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
                        #pragma omp parallel for schedule(static)
                        for ( int d = c+1; d < nvir_d; d++ ){
                           const int ptr_L = jump_D + SIZE_L * nocc_ij * d;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c + ( d * ( d - 1 ) ) / 2 );
                           matmat( 'N', SIZE_L, nocc_ij, SIZE_R, 3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, nocc_ij, SIZE_L, 3.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
                     } else { // irrep_c != irrep_d
                        if ( Iwc < Id ){
                           #pragma omp parallel for schedule(static)
                           for ( int d = 0; d < nvir_d; d++ ){
                              const int ptr_L = jump_D + SIZE_L * nocc_ij * d;
                              const int ptr_R = jump_G + SIZE_R * nocc_ij * ( c + nvir_wc * d );
                              matmat( 'N', SIZE_L, nocc_ij, SIZE_R, 3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                              matmat( 'T', SIZE_R, nocc_ij, SIZE_L, 3.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                           }
                        } else {
                           const int ptr_L = jump_D;
                           const int ptr_R = jump_G + SIZE_R * nocc_ij * nvir_d * c;
                           matmat( 'N', SIZE_L, nvir_d * nocc_ij, SIZE_R, -3.0, workspace, SIZE_L, vector + ptr_R, SIZE_R, result + ptr_L, SIZE_L );
                           matmat( 'T', SIZE_R, nvir_d * nocc_ij, SIZE_L, -3.0, workspace, SIZE_L, vector + ptr_L, SIZE_L, result + ptr_R, SIZE_R );
                        }
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

void CheMPS2::CASPT2::diagonal( double * result ) const{

   #pragma omp parallel
   {

      // FAA: < E_zy E_jx | F | E_ti E_uv > = delta_ij * ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_A[ irrep ];
         if ( SIZE > 0 ){
            const int NOCC = indices->getNOCC( irrep );
            #pragma omp for schedule(static)
            for ( int count = 0; count < NOCC; count++ ){
               const double beta = - f_dot_1dm - fock->get( irrep, count, count );
               double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_A ] + SIZE * count;
               #pragma omp simd
               for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FAA[ irrep ][ elem ] + beta; }
            }
         }
      }

      // FCC: < E_zy E_xb | F | E_at E_uv > = delta_ab * ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_C[ irrep ];
         if ( SIZE > 0 ){
            const int NVIR = indices->getNVIRT( irrep );
            const int N_OA = indices->getNOCC( irrep ) + indices->getNDMRG( irrep );
            #pragma omp for schedule(static)
            for ( int count = 0; count < NVIR; count++ ){
               const double beta = - f_dot_1dm + fock->get( irrep, N_OA + count, N_OA + count );
               double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_C ] + SIZE * count;
               #pragma omp simd
               for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FCC[ irrep ][ elem ] + beta; }
            }
         }
      }

      // FDD: < E_yx E_jb | F | E_ai E_tu > = delta_ab delta_ij ( FDD[ xytu] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ xytu ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_D[ irrep ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_a = Irreps::directProd( irrep_i, irrep );
               const int NOCC_i  = indices->getNOCC( irrep_i );
               const int N_OA_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int SIZE_ia = NOCC_i * indices->getNVIRT( irrep_a );
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < SIZE_ia; combined++ ){
                  const int i = combined % NOCC_i;
                  const int a = combined / NOCC_i;
                  const double f_ii = fock->get( irrep_i, i, i );
                  const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                  const double beta = - f_dot_1dm + f_aa - f_ii;
                  double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + SIZE * ( shift + combined );
                  #pragma omp simd
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = FDD[ irrep ][ elem ] + beta; }
               }
               shift += SIZE_ia;
            }
         }
      }

      int triangle_idx[] = { 0, 0 };

      // FBB singlet: < SB_xkyl | F | SB_tiuj > = 2 delta_ik delta_jl ( FBB_singlet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Iij ][ xytu ] )
      {
         const int SIZE = size_B_singlet[ 0 ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               const int nocc_ij = indices->getNOCC( irrep_ij );
               const int size_ij = ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < size_ij; combined++ ){
                  Special::invert_triangle_two( combined, triangle_idx );
                  const int i = triangle_idx[ 0 ];
                  const int j = triangle_idx[ 1 ];
                  const double f_ii = fock->get( irrep_ij, i, i );
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = - f_dot_1dm - f_ii - f_jj;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * ( shift + combined );
                  #pragma omp simd
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ 0 ][ elem ] + beta ); }
               }
               shift += size_ij;
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
                  const int nocc_i  = indices->getNOCC( irrep_i );
                  const int size_ij = nocc_i * indices->getNOCC( irrep_j );
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < size_ij; combined++ ){
                     const int i = combined % nocc_i;
                     const int j = combined / nocc_i;
                     const double f_ii = fock->get( irrep_i, i, i );
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = - f_dot_1dm - f_ii - f_jj;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + SIZE * ( shift + combined );
                     #pragma omp simd
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_singlet[ irrep ][ elem ] + beta ); }
                  }
                  shift += size_ij;
               }
            }
         }
      }

      // FBB triplet: < TB_xkyl | F | TB_tiuj > = 2 delta_ik delta_jl ( FBB_triplet[ Iij ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Iij ][ xytu ] )
      {
         const int SIZE = size_B_triplet[ 0 ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
               const int nocc_ij = indices->getNOCC( irrep_ij );
               const int size_ij = ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < size_ij; combined++ ){
                  Special::invert_lower_triangle_two( combined, triangle_idx );
                  const int i = triangle_idx[ 0 ];
                  const int j = triangle_idx[ 1 ];
                  const double f_ii = fock->get( irrep_ij, i, i );
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double beta = - f_dot_1dm - f_ii - f_jj;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * ( shift + combined );
                  #pragma omp simd
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ 0 ][ elem ] + beta ); }
               }
               shift += size_ij;
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
                  const int nocc_i  = indices->getNOCC( irrep_i );
                  const int size_ij = nocc_i * indices->getNOCC( irrep_j );
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < size_ij; combined++ ){
                     const int i = combined % nocc_i;
                     const int j = combined / nocc_i;
                     const double f_ii = fock->get( irrep_i, i, i );
                     const double f_jj = fock->get( irrep_j, j, j );
                     const double beta = - f_dot_1dm - f_ii - f_jj;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + SIZE * ( shift + combined );
                     #pragma omp simd
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FBB_triplet[ irrep ][ elem ] + beta ); }
                  }
                  shift += size_ij;
               }
            }
         }
      }

      // FFF singlet: < SF_cxdy | F | SF_atbu > = 2 delta_ac delta_bd ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )
      {
         const int SIZE = size_F_singlet[ 0 ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               const int SIZE_ab = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < SIZE_ab; combined++ ){
                  Special::invert_triangle_two( combined, triangle_idx );
                  const int a = triangle_idx[ 0 ];
                  const int b = triangle_idx[ 1 ];
                  const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = - f_dot_1dm + f_aa + f_bb;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * ( shift + combined );
                  #pragma omp simd
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ 0 ][ elem ] + beta ); }
               }
               shift += SIZE_ab;
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
                  const int NVIR_a  = indices->getNVIRT( irrep_a );
                  const int SIZE_ab = NVIR_a * indices->getNVIRT( irrep_b );
                  const int N_OA_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b  = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < SIZE_ab; combined++ ){
                     const int a = combined % NVIR_a;
                     const int b = combined / NVIR_a;
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = - f_dot_1dm + f_aa + f_bb;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + SIZE * ( shift + combined );
                     #pragma omp simd
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_singlet[ irrep ][ elem ] + beta ); }
                  }
                  shift += SIZE_ab;
               }
            }
         }
      }

      // FFF triplet: < TF_cxdy | F | TF_atbu > = 2 delta_ac delta_bd ( FFF_triplet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Iab ][ xytu ] )
      {
         const int SIZE = size_F_triplet[ 0 ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               const int SIZE_ab = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < SIZE_ab; combined++ ){
                  Special::invert_lower_triangle_two( combined, triangle_idx );
                  const int a = triangle_idx[ 0 ];
                  const int b = triangle_idx[ 1 ];
                  const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double beta = - f_dot_1dm + f_aa + f_bb;
                  double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * ( shift + combined );
                  #pragma omp simd
                  for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ 0 ][ elem ] + beta ); }
               }
               shift += SIZE_ab;
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
                  const int NVIR_a  = indices->getNVIRT( irrep_a );
                  const int SIZE_ab = NVIR_a * indices->getNVIRT( irrep_b );
                  const int N_OA_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b  = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < SIZE_ab; combined++ ){
                     const int a = combined % NVIR_a;
                     const int b = combined / NVIR_a;
                     const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                     const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                     const double beta = - f_dot_1dm + f_aa + f_bb;
                     double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + SIZE * ( shift + combined );
                     #pragma omp simd
                     for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FFF_triplet[ irrep ][ elem ] + beta ); }
                  }
                  shift += SIZE_ab;
               }
            }
         }
      }

      // FEE singlet: < SE_ukbl | F | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_E[ irrep ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int irrep_occ = Irreps::directProd( irrep_a, irrep );
               for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
                  if ( irrep_i == irrep_j ){
                     const int SIZE_aij = ( NVIR_a * NOCC_i * ( NOCC_i + 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_aij; combined++ ){
                        const int a = combined % NVIR_a;
                        Special::invert_triangle_two( combined / NVIR_a, triangle_idx );
                        const int i = triangle_idx[ 0 ];
                        const int j = triangle_idx[ 1 ];
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_jj = fock->get( irrep_j, j, j );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double beta = - f_dot_1dm + f_aa - f_ii - f_jj;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_aij;
                  }
                  if ( irrep_i < irrep_j ){
                     const int SIZE_aij = NVIR_a * NOCC_i * indices->getNOCC( irrep_j );
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_aij; combined++ ){
                        const int a = combined % NVIR_a;
                        const int temp = combined / NVIR_a;
                        const int i = temp % NOCC_i;
                        const int j = temp / NOCC_i;
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_jj = fock->get( irrep_j, j, j );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double beta = - f_dot_1dm + f_aa - f_ii - f_jj;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FEE[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_aij;
                  }
               }
            }
         }
      }

      // FEE triplet: < TE_ukbl | F | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_E[ irrep ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int irrep_occ = Irreps::directProd( irrep_a, irrep );
               for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int irrep_j = Irreps::directProd( irrep_occ, irrep_i );
                  if ( irrep_i == irrep_j ){
                     const int SIZE_aij = ( NVIR_a * NOCC_i * ( NOCC_i - 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_aij; combined++ ){
                        const int a = combined % NVIR_a;
                        Special::invert_lower_triangle_two( combined / NVIR_a, triangle_idx );
                        const int i = triangle_idx[ 0 ];
                        const int j = triangle_idx[ 1 ];
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_jj = fock->get( irrep_j, j, j );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double beta = - f_dot_1dm + f_aa - f_ii - f_jj;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_aij;
                  }
                  if ( irrep_i < irrep_j ){
                     const int SIZE_aij = NVIR_a * NOCC_i * indices->getNOCC( irrep_j );
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_aij; combined++ ){
                        const int a = combined % NVIR_a;
                        const int temp = combined / NVIR_a;
                        const int i = temp % NOCC_i;
                        const int j = temp / NOCC_i;
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_jj = fock->get( irrep_j, j, j );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double beta = - f_dot_1dm + f_aa - f_ii - f_jj;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FEE[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_aij;
                  }
               }
            }
         }
      }

      // FGG singlet: < SG_cjdu | F | SG_aibt > = 2 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_G[ irrep ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int NOCC_i = indices->getNOCC( irrep_i );
               const int irrep_vir = Irreps::directProd( irrep_i, irrep );
               for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
                  if ( irrep_a == irrep_b ){
                     const int SIZE_iab = ( NOCC_i * NVIR_a * ( NVIR_a + 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_iab; combined++ ){
                        const int i = combined % NOCC_i;
                        Special::invert_triangle_two( combined / NOCC_i, triangle_idx );
                        const int a = triangle_idx[ 0 ];
                        const int b = triangle_idx[ 1 ];
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double f_bb = fock->get( irrep_b, N_OA_a + b, N_OA_a + b );
                        const double beta = - f_dot_1dm + f_aa + f_bb - f_ii;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_iab;
                  }
                  if ( irrep_a < irrep_b ){
                     const int SIZE_iab = NOCC_i * NVIR_a * indices->getNVIRT( irrep_b );
                     const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_iab; combined++ ){
                        const int i = combined % NOCC_i;
                        const int temp = combined / NOCC_i;
                        const int a = temp % NVIR_a;
                        const int b = temp / NVIR_a;
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        const double beta = - f_dot_1dm + f_aa + f_bb - f_ii;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 2 * ( FGG[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_iab;
                  }
               }
            }
         }
      }

      // FGG triplet: < TG_cjdu | F | TG_aibt > = 6 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int SIZE = size_G[ irrep ];
         if ( SIZE > 0 ){
            int shift = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int NOCC_i = indices->getNOCC( irrep_i );
               const int irrep_vir = Irreps::directProd( irrep_i, irrep );
               for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int irrep_b = Irreps::directProd( irrep_vir, irrep_a );
                  if ( irrep_a == irrep_b ){
                     const int SIZE_iab = ( NOCC_i * NVIR_a * ( NVIR_a - 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_iab; combined++ ){
                        const int i = combined % NOCC_i;
                        Special::invert_lower_triangle_two( combined / NOCC_i, triangle_idx );
                        const int a = triangle_idx[ 0 ];
                        const int b = triangle_idx[ 1 ];
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double f_bb = fock->get( irrep_b, N_OA_a + b, N_OA_a + b );
                        const double beta = - f_dot_1dm + f_aa + f_bb - f_ii;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_iab;
                  }
                  if ( irrep_a < irrep_b ){
                     const int SIZE_iab = NOCC_i * NVIR_a * indices->getNVIRT( irrep_b );
                     const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < SIZE_iab; combined++ ){
                        const int i = combined % NOCC_i;
                        const int temp = combined / NOCC_i;
                        const int a = temp % NVIR_a;
                        const int b = temp / NVIR_a;
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        const double beta = - f_dot_1dm + f_aa + f_bb - f_ii;
                        double * target = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] + SIZE * ( shift + combined );
                        #pragma omp simd
                        for ( int elem = 0; elem < SIZE; elem++ ){ target[ elem ] = 6 * ( FGG[ irrep ][ elem ] + beta ); }
                     }
                     shift += SIZE_iab;
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
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift;
               const int size_ijab = ( size_ij * NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < size_ijab; combined++ ){
                  Special::invert_triangle_two( combined % size_ij, triangle_idx );
                  const int i = triangle_idx[ 0 ];
                  const int j = triangle_idx[ 1 ];
                  Special::invert_triangle_two( combined / size_ij, triangle_idx );
                  const int a = triangle_idx[ 0 ];
                  const int b = triangle_idx[ 1 ];
                  const double f_ii = fock->get( irrep_ij, i, i );
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double term = f_aa + f_bb - f_ii - f_jj;
                  target[ combined ] = 4 * term;
               }
               shift += size_ijab;
            }
         }
      }
      {
         int shift = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            const int NOCC_ij = indices->getNOCC( irrep_ij );
            const int size_ij = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
            for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
               const int NVIR_ab = indices->getNVIRT( irrep_ab );
               const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
               double * target = result + jump[ num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift;
               const int size_ijab = ( size_ij * NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < size_ijab; combined++ ){
                  Special::invert_lower_triangle_two( combined % size_ij, triangle_idx );
                  const int i = triangle_idx[ 0 ];
                  const int j = triangle_idx[ 1 ];
                  Special::invert_lower_triangle_two( combined / size_ij, triangle_idx );
                  const int a = triangle_idx[ 0 ];
                  const int b = triangle_idx[ 1 ];
                  const double f_ii = fock->get( irrep_ij, i, i );
                  const double f_jj = fock->get( irrep_ij, j, j );
                  const double f_aa = fock->get( irrep_ab, N_OA_ab + a, N_OA_ab + a );
                  const double f_bb = fock->get( irrep_ab, N_OA_ab + b, N_OA_ab + b );
                  const double term = f_aa + f_bb - f_ii - f_jj;
                  target[ combined ] = 12 * term;
               }
               shift += size_ijab;
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
               for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
                  const int irrep_b = Irreps::directProd( irrep, irrep_a );
                  if ( irrep_a < irrep_b ){
                     const int NVIR_a = indices->getNVIRT( irrep_a );
                     const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                     const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                     double * target_singlet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] + shift;
                     double * target_triplet = result + jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] + shift;
                     const int size_ijab = NOCC_i * NOCC_j * NVIR_a * indices->getNVIRT( irrep_b );
                     #pragma omp for schedule(static)
                     for ( int combined = 0; combined < size_ijab; combined++ ){
                        const int i = combined % NOCC_i;
                        const int temp1 = combined / NOCC_i;
                        const int j = temp1 % NOCC_j;
                        const int temp2 = temp1 / NOCC_j;
                        const int a = temp2 % NVIR_a;
                        const int b = temp2 / NVIR_a;
                        const double f_ii = fock->get( irrep_i, i, i );
                        const double f_jj = fock->get( irrep_j, j, j );
                        const double f_aa = fock->get( irrep_a, N_OA_a + a, N_OA_a + a );
                        const double f_bb = fock->get( irrep_b, N_OA_b + b, N_OA_b + b );
                        const double term = f_aa + f_bb - f_ii - f_jj;
                        target_singlet[ combined ] =  4 * term;
                        target_triplet[ combined ] = 12 * term;
                     }
                     shift += size_ijab;
                  }
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

   #pragma omp parallel
   {
      const int max_size = get_maxsize();
      double * workspace = new double[ max_size ];

      // VA
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         if ( size_A[ irrep ] > 0 ){
            const int NOCC = indices->getNOCC( irrep );
            const int NACT = indices->getNDMRG( irrep );
            const int d_w  = indices->getDMRGcumulative( irrep );
            #pragma omp for schedule(static)
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
      }

      // VC
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         if ( size_C[ irrep ] > 0 ){
            const int NOCC = indices->getNOCC( irrep );
            const int NVIR = indices->getNVIRT( irrep );
            const int NACT = indices->getNDMRG( irrep );
            const int N_OA = NOCC + NACT;
            const int d_w  = indices->getDMRGcumulative( irrep );
            #pragma omp for schedule(static)
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
      }
      #pragma omp single
      {
         delete MAT2;
      }

      // VD1 and VD2
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         if ( size_D[ irrep ] > 0 ){
            int shift = 0;
            const int D2JUMP = size_D[ irrep ] / 2;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_a = Irreps::directProd( irrep_i, irrep );
               assert( shift == shift_D_nonactive( indices, irrep_i, irrep_a ) );
               const int NOCC_i = indices->getNOCC( irrep_i );
               const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
               const int NVIR_a = indices->getNVIRT( irrep_a );
               const int loopsize = NOCC_i * NVIR_a;
               #pragma omp for schedule(static)
               for ( int combined = 0; combined < loopsize; combined++ ){

                  const int count_i = combined % NOCC_i;
                  const int count_a = combined / NOCC_i;
                  double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] + size_D[ irrep ] * ( shift + combined );
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
               shift += loopsize;
            }
            assert( shift * size_D[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_D ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_D ] );
         }
      }
      #pragma omp single
      {
         delete MAT;
      }
      const double SQRT_0p5 = sqrt( 0.5 );
      int triangle_idx[] = { 0, 0 };

      // VB singlet and triplet
      if ( size_B_singlet[ 0 ] > 0 ){ // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
         int shift = 0; // First do SINGLET
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, +1 ) );
            const int NOCC_ij = indices->getNOCC( irrep_ij );
            const int loopsize = ( NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
            #pragma omp for schedule(static)
            for ( int combined = 0; combined < loopsize; combined++ ){

               Special::invert_triangle_two( combined, triangle_idx );
               const int i = triangle_idx[ 0 ];
               const int j = triangle_idx[ 1 ];
               assert( combined == i + ( j * ( j + 1 ) ) / 2 );

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
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ 0 ] * ( shift + combined );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
            shift += loopsize;
         }
         assert( shift * size_B_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );
      }
      if ( size_B_triplet[ 0 ] > 0 ){ // Then do TRIPLET
         int shift = 0;
         for ( int irrep_ij = 0; irrep_ij < num_irreps; irrep_ij++ ){
            assert( shift == shift_B_nonactive( indices, irrep_ij, irrep_ij, -1 ) );
            const int NOCC_ij = indices->getNOCC( irrep_ij );
            const int loopsize = ( NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
            #pragma omp for schedule(static)
            for ( int combined = 0; combined < loopsize; combined++ ){

               Special::invert_lower_triangle_two( combined, triangle_idx );
               const int i = triangle_idx[ 0 ];
               const int j = triangle_idx[ 1 ];
               assert( combined == i + ( j * ( j - 1 ) ) / 2 );

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
               double * target = vector_rhs + jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ 0 ] * ( shift + combined );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
            shift += loopsize;
         }
         assert( shift * size_B_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
      }
      for ( int irrep = 1; irrep < num_irreps; irrep++ ){
         assert( size_B_singlet[ irrep ] == size_B_triplet[ irrep ] );
         if ( size_B_singlet[ irrep ] > 0 ){
            int shift = 0;
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep, irrep_i );
               if ( irrep_i < irrep_j ){
                  assert( shift == shift_B_nonactive( indices, irrep_i, irrep_j, 0 ) );
                  const int NOCC_i = indices->getNOCC( irrep_i );
                  const int NOCC_j = indices->getNOCC( irrep_j );
                  const int loopsize = NOCC_i * NOCC_j;
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < loopsize; combined++ ){

                     const int i = combined % NOCC_i;
                     const int j = combined / NOCC_i;

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
                     double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] + size_B_singlet[ irrep ] * ( shift + combined );
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
                     target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] + size_B_triplet[ irrep ] * ( shift + combined );
                     dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SBB_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
                  }
                  shift += loopsize;
               }
            }
            assert( shift * size_B_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_SINGLET ] );
            assert( shift * size_B_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_B_TRIPLET ] );
         }
      }

      // VF singlet and triplet
      if ( size_F_singlet[ 0 ] > 0 ){ // First do irrep == Ii x Ij == Ix x Iy == It x Iu == 0
         int shift = 0; // First do SINGLET
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, +1 ) );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int loopsize = ( NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
            #pragma omp for schedule(static)
            for ( int combined = 0; combined < loopsize; combined++ ){

               Special::invert_triangle_two( combined, triangle_idx );
               const int a = triangle_idx[ 0 ];
               const int b = triangle_idx[ 1 ];
               assert( combined == a + ( b * ( b + 1 ) ) / 2 );

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
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ 0 ] * ( shift + combined );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_singlet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
            shift += loopsize;
         }
         assert( shift * size_F_singlet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );
      }
      if ( size_F_triplet[ 0 ] > 0 ){ // Then do TRIPLET
         int shift = 0;
         for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
            assert( shift == shift_F_nonactive( indices, irrep_ab, irrep_ab, -1 ) );
            const int N_OA_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
            const int NVIR_ab = indices->getNVIRT( irrep_ab );
            const int loopsize = ( NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
            #pragma omp for schedule(static)
            for ( int combined = 0; combined < loopsize; combined++ ){

               Special::invert_lower_triangle_two( combined, triangle_idx );
               const int a = triangle_idx[ 0 ];
               const int b = triangle_idx[ 1 ];
               assert( combined == a + ( b * ( b - 1 ) ) / 2 );

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
               double * target = vector_rhs + jump[ 0 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ 0 ] * ( shift + combined );
               dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ 0 ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
            }
            shift += loopsize;
         }
         assert( shift * size_F_triplet[ 0 ] == jump[ 1 + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
      }
      for ( int irrep = 1; irrep < num_irreps; irrep++ ){
         assert( size_F_singlet[ irrep ] == size_F_triplet[ irrep ] );
         if ( size_F_singlet[ irrep ] > 0 ){
            int shift = 0;
            for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
               const int irrep_b = Irreps::directProd( irrep, irrep_a );
               if ( irrep_a < irrep_b ){
                  assert( shift == shift_F_nonactive( indices, irrep_a, irrep_b, 0 ) );
                  const int N_OA_a = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                  const int N_OA_b = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                  const int NVIR_a = indices->getNVIRT( irrep_a );
                  const int NVIR_b = indices->getNVIRT( irrep_b );
                  const int loopsize = NVIR_a * NVIR_b;
                  #pragma omp for schedule(static)
                  for ( int combined = 0; combined < loopsize; combined++ ){

                     const int a = combined % NVIR_a;
                     const int b = combined / NVIR_a;

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
                     double * target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] + size_F_singlet[ irrep ] * ( shift + combined );
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
                     target = vector_rhs + jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] + size_F_triplet[ irrep ] * ( shift + combined );
                     dgemv_( &notrans, &jump_xy, &jump_xy, &alpha, SFF_triplet[ irrep ], &jump_xy, workspace, &inc1, &set, target, &inc1 );
                  }
                  shift += loopsize;
               }
            }
            assert( shift * size_F_singlet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_SINGLET ] );
            assert( shift * size_F_triplet[ irrep ] == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_F_TRIPLET ] );
         }
      }
      delete [] workspace;

      // VE singlet and triplet
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int occ_t = indices->getNOCC( irrep );
         const int num_t = indices->getNDMRG( irrep );
         if ( num_t > 0 ){
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
                     const int loopsize_singlet = ( NVIR_a * NOCC_ij * ( NOCC_ij + 1 ) ) / 2;
                     const int loopsize_triplet = ( NVIR_a * NOCC_ij * ( NOCC_ij - 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined_singlet = 0; combined_singlet < loopsize_singlet; combined_singlet++ ){

                        const int a = combined_singlet % NVIR_a;
                        const int remainder = combined_singlet / NVIR_a;
                        Special::invert_triangle_two( remainder, triangle_idx );
                        const int i = triangle_idx[ 0 ];
                        const int j = triangle_idx[ 1 ];
                        const int combined_triplet = a + NVIR_a * ( i + ( j * ( j - 1 ) ) / 2 );
                        const double ij_factor = (( i == j ) ? SQRT_0p5 : 1.0 );

                        const int count_aij_singlet = shift_singlet + combined_singlet;
                        const int count_aij_triplet = shift_triplet + combined_triplet;
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
                     shift_singlet += loopsize_singlet;
                     shift_triplet += loopsize_triplet;
                  }
               } else {
                  for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
                     const int irrep_j = Irreps::directProd( irrep_i, irrep_occ );
                     if ( irrep_i < irrep_j ){
                        assert( shift_singlet == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, +1 ) );
                        assert( shift_triplet == shift_E_nonactive( indices, irrep_a, irrep_i, irrep_j, -1 ) );
                        const int NOCC_i = indices->getNOCC( irrep_i );
                        const int NOCC_j = indices->getNOCC( irrep_j );
                        const int loopsize = NOCC_i * NOCC_j * NVIR_a;
                        #pragma omp for schedule(static)
                        for ( int combined = 0; combined < loopsize; combined++ ){

                           const int a = combined % NVIR_a;
                           const int remainder = combined / NVIR_a;
                           const int i = remainder % NOCC_i;
                           const int j = remainder / NOCC_i;

                           const int count_aij_singlet = shift_singlet + combined;
                           const int count_aij_triplet = shift_triplet + combined;
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
                        shift_singlet += loopsize;
                        shift_triplet += loopsize;
                     }
                  }
               }
            }
            assert( shift_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_SINGLET ] );
            assert( shift_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_E_TRIPLET ] );
         }
      }

      // VG singlet and triplet
      for ( int irrep = 0; irrep < num_irreps; irrep++ ){
         const int occ_t = indices->getNOCC( irrep );
         const int num_t = indices->getNDMRG( irrep );
         if ( num_t > 0 ){
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
                     const int loopsize_singlet = ( NOCC_i * NVIR_ab * ( NVIR_ab + 1 ) ) / 2;
                     const int loopsize_triplet = ( NOCC_i * NVIR_ab * ( NVIR_ab - 1 ) ) / 2;
                     #pragma omp for schedule(static)
                     for ( int combined_singlet = 0; combined_singlet < loopsize_singlet; combined_singlet++ ){

                        const int i = combined_singlet % NOCC_i;
                        const int remainder = combined_singlet / NOCC_i;
                        Special::invert_triangle_two( remainder, triangle_idx );
                        const int a = triangle_idx[ 0 ];
                        const int b = triangle_idx[ 1 ];
                        const int combined_triplet = i + NOCC_i * ( a + ( b * ( b - 1 ) ) / 2 );
                        const double ab_factor = (( a == b ) ? SQRT_0p5 : 1.0 );

                        const int count_abi_singlet = shift_singlet + combined_singlet;
                        const int count_abi_triplet = shift_triplet + combined_triplet;
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
                     shift_singlet += loopsize_singlet;
                     shift_triplet += loopsize_triplet;
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
                        const int loopsize = NOCC_i * NVIR_a * NVIR_b;
                        #pragma omp for schedule(static)
                        for ( int combined = 0; combined < loopsize; combined++ ){

                           const int i = combined % NOCC_i;
                           const int remainder = combined / NOCC_i;
                           const int a = remainder % NVIR_a;
                           const int b = remainder / NVIR_a;

                           const int count_abi_singlet = shift_singlet + combined;
                           const int count_abi_triplet = shift_triplet + combined;
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
                        shift_singlet += loopsize;
                        shift_triplet += loopsize;
                     }
                  }
               }
            }
            assert( shift_singlet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_SINGLET ] );
            assert( shift_triplet * num_t == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_G_TRIPLET ] );
         }
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
               const int linsize_ij_singlet = ( nocc_ij * ( nocc_ij + 1 ) ) / 2;
               const int linsize_ij_triplet = ( nocc_ij * ( nocc_ij - 1 ) ) / 2;
               for ( int irrep_ab = 0; irrep_ab < num_irreps; irrep_ab++ ){
                  assert( shift_singlet == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, +1 ) );
                  assert( shift_triplet == shift_H_nonactive( indices, irrep_ij, irrep_ij, irrep_ab, irrep_ab, -1 ) );
                  const int nvirt_ab = indices->getNVIRT( irrep_ab );
                  const int noa_ab = indices->getNOCC( irrep_ab ) + indices->getNDMRG( irrep_ab );
                  const int linsize_ab_singlet = ( nvirt_ab * ( nvirt_ab + 1 ) ) / 2;
                  const int linsize_ab_triplet = ( nvirt_ab * ( nvirt_ab - 1 ) ) / 2;
                  #pragma omp for schedule(static)
                  for ( int combined_ab_singlet = 0; combined_ab_singlet < linsize_ab_singlet; combined_ab_singlet++ ){

                     Special::invert_triangle_two( combined_ab_singlet, triangle_idx );
                     const int a = triangle_idx[ 0 ];
                     const int b = triangle_idx[ 1 ];
                     const double ab_factor = (( a == b ) ? SQRT_0p5 : 1.0 );
                     const int combined_ab_triplet = a + ( b * ( b - 1 ) ) / 2;

                     for ( int i = 0; i < nocc_ij; i++ ){
                        for ( int j = i; j < nocc_ij; j++ ){
                           const double ij_factor = (( i == j ) ? SQRT_0p5 : 1.0 );
                           const int counter_singlet = shift_singlet + i + ( j * ( j + 1 ) ) / 2 + linsize_ij_singlet * combined_ab_singlet;
                           const int counter_triplet = shift_triplet + i + ( j * ( j - 1 ) ) / 2 + linsize_ij_triplet * combined_ab_triplet;

                           const double ai_bj = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, i, j, noa_ab + a, noa_ab + b );
                           const double aj_bi = integrals->get_exchange( irrep_ij, irrep_ij, irrep_ab, irrep_ab, j, i, noa_ab + a, noa_ab + b );
                           target_singlet[ counter_singlet ] = 2 * ( ai_bj + aj_bi ) * ij_factor * ab_factor;
      if ((a<b) && (i<j)){ target_triplet[ counter_triplet ] = 6 * ( ai_bj - aj_bi ); }
                        }
                     }
                  }
                  shift_singlet += linsize_ij_singlet * linsize_ab_singlet;
                  shift_triplet += linsize_ij_triplet * linsize_ab_triplet;
               }
            }
         } else { // irrep_i < irrep_j = irrep_i x irrep   and   irrep_a < irrep_b = irrep_a x irrep
            for ( int irrep_i = 0; irrep_i < num_irreps; irrep_i++ ){
               const int irrep_j = Irreps::directProd( irrep, irrep_i );
               if ( irrep_i < irrep_j ){
                  const int nocc_i = indices->getNOCC( irrep_i );
                  const int nocc_j = indices->getNOCC( irrep_j );
                  const int linsize_ij = nocc_i * nocc_j;
                  for ( int irrep_a = 0; irrep_a < num_irreps; irrep_a++ ){
                     const int irrep_b = Irreps::directProd( irrep, irrep_a );
                     if ( irrep_a < irrep_b ){
                        assert( shift_singlet == shift_H_nonactive( indices, irrep_i, irrep_j, irrep_a, irrep_b, +1 ) );
                        assert( shift_triplet == shift_H_nonactive( indices, irrep_i, irrep_j, irrep_a, irrep_b, -1 ) );
                        const int nvir_a = indices->getNVIRT( irrep_a );
                        const int nvir_b = indices->getNVIRT( irrep_b );
                        const int noa_a  = indices->getNOCC( irrep_a ) + indices->getNDMRG( irrep_a );
                        const int noa_b  = indices->getNOCC( irrep_b ) + indices->getNDMRG( irrep_b );
                        const int linsize_ab = nvir_a * nvir_b;
                        #pragma omp for schedule(static)
                        for ( int combined_ab = 0; combined_ab < linsize_ab; combined_ab++ ){

                           const int a = combined_ab % nvir_a;
                           const int b = combined_ab / nvir_a;

                           for ( int j = 0; j < nocc_j; j++ ){
                              for ( int i = 0; i < nocc_i; i++ ){
                                 const int count_singlet = shift_singlet + i + nocc_i * ( j + nocc_j * combined_ab );
                                 const int count_triplet = shift_triplet + i + nocc_i * ( j + nocc_j * combined_ab );
                                 const double ai_bj = integrals->get_exchange( irrep_i, irrep_j, irrep_a, irrep_b, i, j, noa_a + a, noa_b + b );
                                 const double aj_bi = integrals->get_exchange( irrep_j, irrep_i, irrep_a, irrep_b, j, i, noa_a + a, noa_b + b );
                                 target_singlet[ count_singlet ] = 2 * ( ai_bj + aj_bi );
                                 target_triplet[ count_triplet ] = 6 * ( ai_bj - aj_bi );
                              }
                           }
                        }
                        shift_singlet += linsize_ij * linsize_ab;
                        shift_triplet += linsize_ij * linsize_ab;
                     }
                  }
               }
            }
         }
         assert( shift_singlet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_SINGLET ] );
         assert( shift_triplet == jump[ 1 + irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] - jump[ irrep + num_irreps * CHEMPS2_CASPT2_H_TRIPLET ] );
      }
   }// End of #pragma omp parallel

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

void CheMPS2::CASPT2::make_AA_CC( const bool OVLP, const double IPEA ){

   /*
      SAA: < E_zy E_jx | 1 | E_ti E_uv > = delta_ij ( SAA[ Ii ][ xyztuv ] )
      FAA: < E_zy E_jx | F | E_ti E_uv > = delta_ij ( FAA[ Ii ][ xyztuv ] + ( 2 sum_k f_kk - f_ii ) SAA[ Ii ][ xyztuv ] )

         SAA[ Ii ][ xyztuv ] = ( - Gamma_{ztuyxv}
                                 + 2 delta_tx Gamma_{zuyv}
                                 -   delta_uy Gamma_{tzxv}
                                 -   delta_ty Gamma_{zuxv}
                                 -   delta_ux Gamma_{ztyv}
                                 + ( 2 delta_tx delta_uy - delta_ux delta_ty ) Gamma_{zv}
                               )

         FAA[ Ii ][ xyztuv ] = ( - f_dot_4dm[ ztuyxv ] + ( f_tt + f_uu + f_xx + f_yy ) SAA[ Ii ][ xyztuv ]
                                 + 2 delta_tx ( f_dot_3dm[ zuyv ] - f_tt Gamma_{zuyv} )
                                 -   delta_uy ( f_dot_3dm[ ztvx ] - f_uu Gamma_{ztvx} )
                                 -   delta_ty ( f_dot_3dm[ zuxv ] - f_tt Gamma_{zuxv} )
                                 -   delta_ux ( f_dot_3dm[ ztyv ] - f_uu Gamma_{ztyv} )
                                 + ( 2 delta_tx delta_uy - delta_ux delta_ty ) ( f_dot_2dm[ zv ] - ( f_tt + f_uu ) Gamma_{zv} )
                               )

      SCC: < E_zy E_xb | 1 | E_at E_uv > = delta_ab ( SCC[ Ia ][ xyztuv ] )
      FCC: < E_zy E_xb | F | E_at E_uv > = delta_ab ( FCC[ Ia ][ xyztuv ] + ( 2 sum_k f_kk + f_aa ) SCC[ Ia ][ xyztuv ] )

         SCC[ Ia ][ xyztuv ] = ( + Gamma_{zxuytv}
                                 + delta_uy Gamma_{xztv}
                                 + delta_xy Gamma_{zutv}
                                 + delta_ut Gamma_{zxyv}
                                 + delta_ut delta_xy Gamma_{zv}
                               )

         FCC[ Ia ][ xyztuv ] = ( + f_dot_4dm[ zxuytv ] + ( f_uu + f_yy ) SCC[ Ia ][ xyztuv ]
                                 + delta_uy ( f_dot_3dm[ xztv ] - f_uu Gamma_{xztv} )
                                 + delta_xy ( f_dot_3dm[ zutv ] - f_xx Gamma_{zutv} )
                                 + delta_ut ( f_dot_3dm[ zxyv ] - f_tt Gamma_{zxyv} )
                                 + delta_ut delta_xy ( f_dot_2dm[ zv ] - ( f_tt + f_xx ) Gamma_{zv} )
                               )
   */

   if ( OVLP ){ SAA = new double*[ num_irreps ];
                SCC = new double*[ num_irreps ]; }
   else {       FAA = new double*[ num_irreps ];
                FCC = new double*[ num_irreps ]; }

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      assert( size_A[ irrep ] == size_C[ irrep ] ); // At construction
      const int SIZE = size_A[ irrep ];
      if ( OVLP ){ SAA[ irrep ] = new double[ SIZE * SIZE ];
                   SCC[ irrep ] = new double[ SIZE * SIZE ]; }
      else {       FAA[ irrep ] = new double[ SIZE * SIZE ];
                   FCC[ irrep ] = new double[ SIZE * SIZE ]; }

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

                  // SAA: - Gamma_{ztuyxv}
                  // FAA: - f_dot_4dm[ ztuyxv ] + ( f_tt + f_uu + f_xx + f_yy ) SAA[ Ii ][ xyztuv ]
                  if ( OVLP ){
                     #pragma omp parallel for schedule(static)
                     for ( int v = 0; v < num_v; v++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       const int ptr = jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) );
                                       SAA[ irrep ][ ptr ] = - three_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_v + v ))))) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     #pragma omp parallel for schedule(static)
                     for ( int v = 0; v < num_v; v++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_uu = fock->get( irrep_u, nocc_u + u, nocc_u + u );
                           for ( int t = 0; t < num_t; t++ ){
                              const double f_tt = fock->get( irrep_t, nocc_t + t, nocc_t + t );
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    const double f_yy = fock->get( irrep_y, nocc_y + y, nocc_y + y );
                                    for ( int x = 0; x < num_x; x++ ){
                                       const double f_xx = fock->get( irrep_x, nocc_x + x, nocc_x + x );
                                       const int ptr = jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) );
                                       FAA[ irrep ][ ptr ] = - f_dot_4dm[ d_z + z + LAS * ( d_t + t + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_v + v ))))) ]
                                                             + ( f_tt + f_uu + f_xx + f_yy ) * SAA[ irrep ][ ptr ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SCC: + Gamma_{zxuytv}
                  // FCC: + f_dot_4dm[ zxuytv ] + ( f_uu + f_yy ) SCC[ Ia ][ xyztuv ]
                  if ( OVLP ){
                     #pragma omp parallel for schedule(static)
                     for ( int v = 0; v < num_v; v++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       const int ptr = jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) );
                                       SCC[ irrep ][ ptr ] = three_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_v + v ))))) ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  } else {
                     #pragma omp parallel for schedule(static)
                     for ( int v = 0; v < num_v; v++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_uu = fock->get( irrep_u, nocc_u + u, nocc_u + u );
                           for ( int t = 0; t < num_t; t++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    const double f_yy = fock->get( irrep_y, nocc_y + y, nocc_y + y );
                                    for ( int x = 0; x < num_x; x++ ){
                                       const int ptr = jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) );
                                       FCC[ irrep ][ ptr ] = f_dot_4dm[ d_z + z + LAS * ( d_x + x + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_v + v ))))) ]
                                                           + ( f_yy + f_uu ) * SCC[ irrep ][ ptr ];
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: + 2 delta_tx Gamma_{zuyv}
                  // FAA: + 2 delta_tx ( f_dot_3dm[ zuyv ] - f_tt Gamma_{zuyv} )
                  if ( irrep_t == irrep_x ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int xt = 0; xt < num_t; xt++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       SAA[ irrep ][ jump_row + xt + num_x * ( y + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( u + num_u * v ) ) ]
                                          += 2 * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int tx = 0; tx < num_t; tx++ ){
                                 const double f_tt = fock->get( irrep_t, nocc_t + tx, nocc_t + tx );
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       FAA[ irrep ][ jump_row + tx + num_x * ( y + num_y * z ) + SIZE * ( jump_col + tx + num_t * ( u + num_u * v ) ) ]
                                          += 2 * ( f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ))) ]
                                            - f_tt * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_y + y + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: - delta_uy Gamma_{ztvx}
                  // FAA: - delta_uy ( f_dot_3dm[ ztvx ] - f_uu Gamma_{ztvx} )
                  if ( irrep_u == irrep_y ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int uy = 0; uy < num_u; uy++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       SAA[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                          -= two_rdm[ d_t + t + LAS * ( d_z + z + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int uy = 0; uy < num_u; uy++ ){
                              const double f_uu = fock->get( irrep_u, nocc_u + uy, nocc_u + uy );
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                          -= ( f_dot_3dm[ d_t + t + LAS * ( d_z + z + LAS * ( d_x + x + LAS * ( d_v + v ))) ]
                                        - f_uu * two_rdm[ d_t + t + LAS * ( d_z + z + LAS * ( d_x + x + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: - delta_ty Gamma_{zuxv}
                  // FAA: - delta_ty ( f_dot_3dm[ zuxv ] - f_tt Gamma_{zuxv} )
                  if ( irrep_t == irrep_y ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int ty = 0; ty < num_t; ty++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       SAA[ irrep ][ jump_row + x + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( u + num_u * v ) ) ]
                                          -= two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int ty = 0; ty < num_t; ty++ ){
                                 const double f_tt = fock->get( irrep_t, nocc_t + ty, nocc_t + ty );
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FAA[ irrep ][ jump_row + x + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( u + num_u * v ) ) ]
                                          -= ( f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ))) ]
                                        - f_tt * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_x + x + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: - delta_ux Gamma_{ztyv}
                  // FAA: - delta_ux ( f_dot_3dm[ ztyv ] - f_uu Gamma_{ztyv} )
                  if ( irrep_u == irrep_x ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int ux = 0; ux < num_u; ux++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       SAA[ irrep ][ jump_row + ux + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( ux + num_u * v ) ) ]
                                          -= two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int ux = 0; ux < num_u; ux++ ){
                              const double f_uu = fock->get( irrep_u, nocc_u + ux, nocc_u + ux );
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int y = 0; y < num_y; y++ ){
                                       FAA[ irrep ][ jump_row + ux + num_x * ( y + num_y * z ) + SIZE * ( jump_col + t + num_t * ( ux + num_u * v ) ) ]
                                          -= ( f_dot_3dm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ))) ]
                                        - f_uu * two_rdm[ d_z + z + LAS * ( d_t + t + LAS * ( d_y + y + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SCC: + delta_uy Gamma_{xztv}
                  // FCC: + delta_uy ( f_dot_3dm[ xztv ] - f_uu Gamma_{xztv} )
                  if ( irrep_u == irrep_y ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int uy = 0; uy < num_u; uy++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       SCC[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_u * v ) ) ]
                                          += two_rdm[ d_x + x + LAS * ( d_z + z + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int uy = 0; uy < num_y; uy++ ){
                              const double f_uu = fock->get( irrep_y, nocc_y + uy, nocc_y + uy );
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( uy + num_y * v ) ) ]
                                          += ( f_dot_3dm[ d_x + x + LAS * ( d_z + z + LAS * ( d_t + t + LAS * ( d_v + v ))) ]
                                        - f_uu * two_rdm[ d_x + x + LAS * ( d_z + z + LAS * ( d_t + t + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SCC: + delta_xy Gamma_{zutv}
                  // FCC: + delta_xy ( f_dot_3dm[ zutv ] - f_xx Gamma_{zutv} )
                  if ( irrep_x == irrep_y ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       SCC[ irrep ][ jump_row + xy + num_x * ( xy + num_y * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          += two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int u = 0; u < num_u; u++ ){
                              for ( int t = 0; t < num_t; t++ ){
                                 for ( int z = 0; z < num_z; z++ ){
                                    for ( int xy = 0; xy < num_x; xy++ ){
                                       const double f_xx = fock->get( irrep_x, nocc_x + xy, nocc_x + xy );
                                       FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + t + num_t * ( u + num_u * v ) ) ]
                                          += ( f_dot_3dm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ))) ]
                                        - f_xx * two_rdm[ d_z + z + LAS * ( d_u + u + LAS * ( d_t + t + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SCC: + delta_ut Gamma_{zxyv}
                  // FCC: + delta_ut ( f_dot_3dm[ zxyv ] - f_tt Gamma_{zxyv} )
                  if ( irrep_u == irrep_t ){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int ut = 0; ut < num_u; ut++ ){
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       SCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + ut + num_t * ( ut + num_u * v ) ) ]
                                          += two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ))) ];
                                    }
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int ut = 0; ut < num_u; ut++ ){
                              const double f_tt = fock->get( irrep_t, nocc_t + ut, nocc_t + ut );
                              for ( int z = 0; z < num_z; z++ ){
                                 for ( int y = 0; y < num_y; y++ ){
                                    for ( int x = 0; x < num_x; x++ ){
                                       FCC[ irrep ][ jump_row + x + num_x * ( y + num_y * z ) + SIZE * ( jump_col + ut + num_u * ( ut + num_u * v ) ) ]
                                          += ( f_dot_3dm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ))) ]
                                        - f_tt * two_rdm[ d_z + z + LAS * ( d_x + x + LAS * ( d_y + y + LAS * ( d_v + v ))) ] );
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: + 2 delta_tx delta_uy Gamma_{zv}
                  // FAA: + 2 delta_tx delta_uy ( f_dot_2dm[ zv ] - ( f_tt + f_uu ) Gamma_{zv} )
                  if (( irrep_t == irrep_x ) && ( irrep_u == irrep_y ) && ( irrep_z == irrep_v )){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int uy = 0; uy < num_u; uy++ ){
                                 for ( int xt = 0; xt < num_t; xt++ ){
                                    SAA[ irrep ][ jump_row + xt + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( uy + num_u * v ) ) ]
                                       += 2 * one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int uy = 0; uy < num_u; uy++ ){
                                 const double f_uu = fock->get( irrep_u, nocc_u + uy, nocc_u + uy );
                                 for ( int xt = 0; xt < num_t; xt++ ){
                                    const double f_tt = fock->get( irrep_t, nocc_t + xt, nocc_t + xt );
                                    FAA[ irrep ][ jump_row + xt + num_x * ( uy + num_y * z ) + SIZE * ( jump_col + xt + num_t * ( uy + num_u * v ) ) ]
                                       += 2 * ( f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ] - ( f_tt + f_uu ) * one_rdm[ d_z + z + LAS * ( d_v + v ) ] );
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SAA: - delta_ux delta_ty Gamma_{zv}
                  // FAA: - delta_ux delta_ty ( f_dot_2dm[ zv ] - ( f_tt + f_uu ) Gamma_{zv} )
                  if (( irrep_u == irrep_x ) && ( irrep_t == irrep_y ) && ( irrep_z == irrep_v )){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int ux = 0; ux < num_u; ux++ ){
                                 for ( int ty = 0; ty < num_t; ty++ ){
                                    SAA[ irrep ][ jump_row + ux + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( ux + num_u * v ) ) ]
                                       -= one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int ux = 0; ux < num_u; ux++ ){
                                 const double f_uu = fock->get( irrep_u, nocc_u + ux, nocc_u + ux );
                                 for ( int ty = 0; ty < num_t; ty++ ){
                                    const double f_tt = fock->get( irrep_t, nocc_t + ty, nocc_t + ty );
                                    FAA[ irrep ][ jump_row + ux + num_x * ( ty + num_y * z ) + SIZE * ( jump_col + ty + num_t * ( ux + num_u * v ) ) ]
                                       -= ( f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ] - ( f_tt + f_uu ) * one_rdm[ d_z + z + LAS * ( d_v + v ) ] );
                                 }
                              }
                           }
                        }
                     }
                  }

                  // SCC: + delta_ut delta_xy Gamma_{zv}
                  // FCC: + delta_ut delta_xy ( f_dot_2dm[ zv ] - ( f_tt + f_xx ) Gamma_{zv} )
                  if (( irrep_u == irrep_t ) && ( irrep_x == irrep_y ) && ( irrep_z == irrep_v )){
                     if ( OVLP ){
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_z; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int tu = 0; tu < num_t; tu++ ){
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    SCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + tu + num_t * ( tu + num_t * v ) ) ]
                                       += one_rdm[ d_z + z + LAS * ( d_v + v ) ];
                                 }
                              }
                           }
                        }
                     } else {
                        #pragma omp parallel for schedule(static)
                        for ( int v = 0; v < num_v; v++ ){
                           for ( int z = 0; z < num_z; z++ ){
                              for ( int tu = 0; tu < num_t; tu++ ){
                                 const double f_tt = fock->get( irrep_t, nocc_t + tu, nocc_t + tu );
                                 for ( int xy = 0; xy < num_x; xy++ ){
                                    const double f_xx = fock->get( irrep_x, nocc_x + xy, nocc_x + xy );
                                    FCC[ irrep ][ jump_row + xy + num_x * ( xy + num_x * z ) + SIZE * ( jump_col + tu + num_t * ( tu + num_t * v ) ) ]
                                       += ( f_dot_2dm[ d_z + z + LAS * ( d_v + v ) ] - ( f_tt + f_xx ) * one_rdm[ d_z + z + LAS * ( d_v + v ) ] );
                                 }
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y * num_z;
               }
            }
            if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
               // A: E_ti E_uv | 0 >   --->   t: excitation into,   u: excitation into, v: excitation out of
               // C: E_at E_uv | 0 >   --->   t: excitation out of, u: excitation into, v: excitation out of
               #pragma omp parallel for schedule(static)
               for ( int v = 0; v < num_v; v++ ){
                  const double gamma_vv = one_rdm[ ( d_v + v ) * ( 1 + LAS ) ];
                  for ( int u = 0; u < num_u; u++ ){
                     const double gamma_uu = one_rdm[ ( d_u + u ) * ( 1 + LAS ) ];
                     for ( int t = 0; t < num_t; t++ ){
                        const double gamma_tt = one_rdm[ ( d_t + t ) * ( 1 + LAS ) ];
                        const double ipea_A_tuv = 0.5 * IPEA * ( 2.0 + gamma_tt + gamma_uu - gamma_vv );
                        const double ipea_C_tuv = 0.5 * IPEA * ( 4.0 - gamma_tt + gamma_uu - gamma_vv );
                        const int ptr = ( jump_col + t + num_t * ( u + num_u * v ) ) * ( 1 + SIZE );
                        FAA[ irrep ][ ptr ] += ipea_A_tuv * SAA[ irrep ][ ptr ];
                        FCC[ irrep ][ ptr ] += ipea_C_tuv * SCC[ irrep ][ ptr ];
                     }
                  }
               }
            }
            jump_col += num_t * num_u * num_v;
         }
      }
   }

}

void CheMPS2::CASPT2::make_DD( const bool OVLP, const double IPEA ){

   /*
      SDD: < D(bjxy) | 1 | D(aitu) > = delta_ab delta_ij ( SDD[ It x Iu ][ xytu ] )
      FDD: < D(bjxy) | F | D(aitu) > = delta_ab delta_ij ( FDD[ It x Iu ][ xytu ] + ( 2 sum_k f_kk + f_aa - f_ii ) SDD[ It x Iu ][ xytu ] )

            SD1D1[ It x Iu ][ xytu ] = ( + 2 * Gamma_{ytxu}
                                         + 2 * delta_tx Gamma_{yu}
                                       )

            FD1D1[ It x Iu ][ xytu ] = ( + 2 * f_dot_3dm[ ytxu ] + ( f_xx + f_tt ) SD1D1[ It x Iu ][ xytu ]
                                         + 2 * delta_tx ( f_dot_2dm[ yu ] - f_tt Gamma_{yu} )
                                       )

            SD2D2[ It x Iu ][ xytu ] = ( - Gamma_{ytux}
                                         + 2 * delta_tx Gamma_{yu}
                                       )

            FD2D2[ It x Iu ][ xytu ] = ( - f_dot_3dm[ ytux ] + ( f_xx + f_tt ) SD2D2[ It x Iu ][ xytu ]
                                         + 2 * delta_tx ( f_dot_2dm[ yu ] - f_tt Gamma_{yu} )
                                       )

            SD1D2[ It x Iu ][ xytu ] = ( - Gamma_{ytxu}
                                         - delta_tx Gamma_{yu}
                                       )

            FD1D2[ It x Iu ][ xytu ] = ( - f_dot_3dm[ ytxu ] + ( f_xx + f_tt ) SD1D2[ It x Iu ][ xytu ]
                                         - delta_tx ( f_dot_2dm[ yu ] - f_tt Gamma_{yu} )
                                       )

            SD2D1[ It x Iu ][ xytu ] = ( - Gamma_{ytxu}
                                         - delta_tx Gamma_{yu}
                                       )

            FD2D1[ It x Iu ][ xytu ] = ( - f_dot_3dm[ ytxu ] + ( f_xx + f_tt ) SD2D1[ It x Iu ][ xytu ]
                                         - delta_tx ( f_dot_2dm[ yu ] - f_tt Gamma_{yu} )
                                       )
   */

   if ( OVLP ){ SDD = new double*[ num_irreps ]; }
   else {       FDD = new double*[ num_irreps ]; }

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep = 0; irrep < num_irreps; irrep++ ){

      const int SIZE   = size_D[ irrep ];
      const int D2JUMP = SIZE / 2;
      if ( OVLP ){ SDD[ irrep ] = new double[ SIZE * SIZE ]; }
      else {       FDD[ irrep ] = new double[ SIZE * SIZE ]; }

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

            if ( OVLP ){
               #pragma omp parallel for schedule(static)
               for ( int u = 0; u < num_u; u++ ){
                  for ( int t = 0; t < num_t; t++ ){
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           const double gamma_ytxu = two_rdm[ d_y + y + LAS * ( d_t + t + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                           const double gamma_ytux = two_rdm[ d_t + t + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                           const int ptr = jump_row + x + num_x * y + SIZE * ( jump_col + t + num_t * u );
                           SDD[ irrep ][ ptr                          ] = 2 * gamma_ytxu;
                           SDD[ irrep ][ ptr +          SIZE * D2JUMP ] =   - gamma_ytxu;
                           SDD[ irrep ][ ptr + D2JUMP                 ] =   - gamma_ytxu;
                           SDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ] =   - gamma_ytux;
                        }
                     }
                  }
               }
            } else {
               #pragma omp parallel for schedule(static)
               for ( int u = 0; u < num_u; u++ ){
                  for ( int t = 0; t < num_t; t++ ){
                     const double f_tt = fock->get( irrep_t, nocc_t + t, nocc_t + t );
                     for ( int y = 0; y < num_y; y++ ){
                        for ( int x = 0; x < num_x; x++ ){
                           const double f_xx = fock->get( irrep_x, nocc_x + x, nocc_x + x );
                           const double f_3dm_ytxu = f_dot_3dm[ d_y + y + LAS * ( d_t + t + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                           const double f_3dm_ytux = f_dot_3dm[ d_t + t + LAS * ( d_y + y + LAS * ( d_x + x + LAS * ( d_u + u ))) ];
                           const int ptr = jump_row + x + num_x * y + SIZE * ( jump_col + t + num_t * u );
                           FDD[ irrep ][ ptr                          ] = 2 * f_3dm_ytxu + ( f_xx + f_tt ) * SDD[ irrep ][ ptr                          ];
                           FDD[ irrep ][ ptr +          SIZE * D2JUMP ] =   - f_3dm_ytxu + ( f_xx + f_tt ) * SDD[ irrep ][ ptr +          SIZE * D2JUMP ];
                           FDD[ irrep ][ ptr + D2JUMP                 ] =   - f_3dm_ytxu + ( f_xx + f_tt ) * SDD[ irrep ][ ptr + D2JUMP                 ];
                           FDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ] =   - f_3dm_ytux + ( f_xx + f_tt ) * SDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ];
                        }
                     }
                  }
               }
            }

            if (( irrep_x == irrep_t ) && ( irrep_y == irrep_u )){
               if ( OVLP ){
                  #pragma omp parallel for schedule(static)
                  for ( int u = 0; u < num_y; u++ ){
                     for ( int xt = 0; xt < num_x; xt++ ){
                        for ( int y = 0; y < num_y; y++ ){
                           const double gamma_yu = one_rdm[ d_y + y + LAS * ( d_y + u ) ];
                           const int ptr = jump_row + xt + num_x * y + SIZE * ( jump_col + xt + num_x * u );
                           SDD[ irrep ][ ptr                          ] += 2 * gamma_yu;
                           SDD[ irrep ][ ptr +          SIZE * D2JUMP ] -=     gamma_yu;
                           SDD[ irrep ][ ptr + D2JUMP                 ] -=     gamma_yu;
                           SDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ] += 2 * gamma_yu;
                        }
                     }
                  }
               } else {
                  #pragma omp parallel for schedule(static)
                  for ( int u = 0; u < num_y; u++ ){
                     for ( int xt = 0; xt < num_x; xt++ ){
                        const double f_tt = fock->get( irrep_t, nocc_t + xt, nocc_t + xt );
                        for ( int y = 0; y < num_y; y++ ){
                           const double val_yu = ( f_dot_2dm[ d_y + y + LAS * ( d_y + u ) ]
                                            - f_tt * one_rdm[ d_y + y + LAS * ( d_y + u ) ] );
                           const int ptr = jump_row + xt + num_x * y + SIZE * ( jump_col + xt + num_x * u );
                           FDD[ irrep ][ ptr                          ] += 2 * val_yu;
                           FDD[ irrep ][ ptr +          SIZE * D2JUMP ] -=     val_yu;
                           FDD[ irrep ][ ptr + D2JUMP                 ] -=     val_yu;
                           FDD[ irrep ][ ptr + D2JUMP + SIZE * D2JUMP ] += 2 * val_yu;
                        }
                     }
                  }
               }
            }
            jump_row += num_x * num_y;
         }
         if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
            // D1: E_ai E_tu | 0 >   --->   t: excitation into, u excitation out of
            // D2: E_ti E_au | 0 >   --->   t: excitation into, u excitation out of
            #pragma omp parallel for schedule(static)
            for ( int u = 0; u < num_u; u++ ){
               const double gamma_uu = one_rdm[ ( d_u + u ) * ( 1 + LAS ) ];
               for ( int t = 0; t < num_t; t++ ){
                  const double gamma_tt = one_rdm[ ( d_t + t ) * ( 1 + LAS ) ];
                  const double ipea_tu = 0.5 * IPEA * ( 2.0 + gamma_tt - gamma_uu );
                  const int ptr1 = (          jump_col + t + num_t * u ) * ( 1 + SIZE );
                  const int ptr2 = ( D2JUMP + jump_col + t + num_t * u ) * ( 1 + SIZE );
                  FDD[ irrep ][ ptr1 ] += ipea_tu * SDD[ irrep ][ ptr1 ];
                  FDD[ irrep ][ ptr2 ] += ipea_tu * SDD[ irrep ][ ptr2 ];
               }
            }
         }
         jump_col += num_t * num_u;
      }
   }

}

void CheMPS2::CASPT2::make_BB_FF_singlet( const bool OVLP, const double IPEA ){

   /*
      | SB_tiuj > = ( E_ti E_uj + E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 >  with  i <= j and t <= u

      SBB singlet: < SB_xkyl | 1 | SB_tiuj > = 2 delta_ik delta_jl ( SBB_singlet[ Itu ][ xytu ] )
      FBB singlet: < SB_xkyl | F | SB_tiuj > = 2 delta_ik delta_jl ( FBB_singlet[ Itu ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_singlet[ Itu ][ xytu ] )

            SBB_singlet[ Itu ][ xytu ] = ( + Gamma_{utyx}
                                           + Gamma_{utxy}
                                           + 2 ( delta_uy delta_tx + delta_ux delta_ty )
                                           - delta_uy Gamma_{tx}
                                           - delta_tx Gamma_{uy}
                                           - delta_ux Gamma_{ty}
                                           - delta_ty Gamma_{ux}
                                         )

            FBB_singlet[ Itu ][ xytu ] = ( + f_dot_3dm[ utyx ]
                                           + f_dot_3dm[ tuyx ]
                                           + ( f_tt + f_uu + f_xx + f_yy ) SBB_singlet[ xytu ]
                                           + 2 ( delta_uy delta_tx + delta_ux delta_ty ) ( f_dot_1dm - f_tt - f_uu )
                                           -   delta_uy ( f_dot_2dm[ tx ] - f_uu Gamma_{tx} )
                                           -   delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
                                           -   delta_ux ( f_dot_2dm[ ty ] - f_uu Gamma_{ty} )
                                           -   delta_ty ( f_dot_2dm[ ux ] - f_tt Gamma_{ux} )
                                         )

      | SF_atbu > = ( E_at E_bu + E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 >  with  a <= b and t <= u

      SFF singlet: < SF_cxdy | 1 | SF_atbu > = 2 delta_ac delta_bd ( SFF_singlet[ Itu ][ xytu ] )
      FFF singlet: < SF_cxdy | F | SF_atbu > = 2 delta_ac delta_bd ( FFF_singlet[ Iab ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_singlet[ Iab ][ xytu ] )

            SFF_singlet[ Itu ][ xytu ] = ( + Gamma_{yxut}
                                           + Gamma_{yxtu}
                                         )

            FFF_singlet[ Itu ][ xytu ] = ( + f_dot_3dm[ yxut ]
                                           + f_dot_3dm[ yxtu ]
                                         )
   */

   if ( OVLP ){ SBB_singlet = new double*[ num_irreps ];
                SFF_singlet = new double*[ num_irreps ]; }
   else {       FBB_singlet = new double*[ num_irreps ];
                FFF_singlet = new double*[ num_irreps ]; }

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_singlet[ 0 ] == size_F_singlet[ 0 ] ); // At construction
      const int SIZE = size_B_singlet[ 0 ];
      if ( OVLP ){ SBB_singlet[ 0 ] = new double[ SIZE * SIZE ];
                   SFF_singlet[ 0 ] = new double[ SIZE * SIZE ]; }
      else {       FBB_singlet[ 0 ] = new double[ SIZE * SIZE ];
                   FFF_singlet[ 0 ] = new double[ SIZE * SIZE ]; }

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

            // SBB singlet: + Gamma_{utyx} + Gamma_{utxy}
            // SFF singlet: + Gamma_{yxut} + Gamma_{yxtu}
            // FBB singlet: + f_dot_3dm[ utyx ] + f_dot_3dm[ tuyx ] + ( f_tt + f_uu + f_xx + f_yy ) SBB_singlet[ xytu ]
            // FFF singlet: + f_dot_3dm[ yxut ] + f_dot_3dm[ yxtu ]
            if ( OVLP ){
               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                     for ( int x = 0; x < num_xy; x++ ){
                        for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                           const double value = ( two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                                + two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                           const int ptr = shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 );
                           SBB_singlet[ 0 ][ ptr ] = value;
                           SFF_singlet[ 0 ][ ptr ] = value;
                        }
                     }
                  }
               }
            } else {
               for ( int t = 0; t < num_ut; t++ ){
                  const double f_tt = fock->get( irrep_ut, nocc_ut + t, nocc_ut + t );
                  for ( int u = t; u < num_ut; u++ ){ // 0 <= t <= u < num_ut
                     const double f_uu = fock->get( irrep_ut, nocc_ut + u, nocc_ut + u );
                     for ( int x = 0; x < num_xy; x++ ){
                        const double f_xx = fock->get( irrep_xy, nocc_xy + x, nocc_xy + x );
                        for ( int y = x; y < num_xy; y++ ){ // 0 <= x <= y < num_xy
                           const double f_yy = fock->get( irrep_xy, nocc_xy + y, nocc_xy + y );
                           const int ptr = shift + x + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 );
                           const double fdotsum = ( f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                                  + f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                           FBB_singlet[ 0 ][ ptr ] = fdotsum + ( f_tt + f_uu + f_xx + f_yy ) * SBB_singlet[ 0 ][ ptr ];
                           FFF_singlet[ 0 ][ ptr ] = fdotsum;
                        }
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal: num_ut = num_xy and d_ut = d_xy

               // SBB singlet: + 2 ( delta_uy delta_tx + delta_ux delta_ty )
               // FBB singlet: + 2 ( delta_uy delta_tx + delta_ux delta_ty ) ( f_dot_1dm - ( f_tt + f_uu ) )
               if ( OVLP ){
                  for ( int t = 0; t < num_ut; t++ ){
                     SBB_singlet[ 0 ][ shift + t + ( t * ( t + 1 ) ) / 2 + SIZE * ( t + ( t * ( t + 1 ) ) / 2 ) ] += 4.0;
                     for ( int u = t+1; u < num_ut; u++ ){
                        SBB_singlet[ 0 ][ shift + t + ( u * ( u + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] += 2.0;
                     }
                  }
               } else {
                  for ( int t = 0; t < num_ut; t++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + t, nocc_ut + t );
                     FBB_singlet[ 0 ][ shift + t + ( t * ( t + 1 ) ) / 2 + SIZE * ( t + ( t * ( t + 1 ) ) / 2 ) ] += 4 * ( f_dot_1dm - 2 * f_tt );
                     for ( int u = t+1; u < num_ut; u++ ){
                        const double f_uu = fock->get( irrep_ut, nocc_ut + u, nocc_ut + u );
                        FBB_singlet[ 0 ][ shift + t + ( u * ( u + 1 ) ) / 2 + SIZE * ( t + ( u * ( u + 1 ) ) / 2 ) ] += 2 * ( f_dot_1dm - f_tt - f_uu );
                     }
                  }
               }

               // SBB singlet: - delta_uy Gamma_{tx}
               // FBB singlet: - delta_uy ( f_dot_2dm[ tx ] - f_uu Gamma_{tx} )
               if ( OVLP ){
                  for ( int uy = 0; uy < num_ut; uy++ ){
                     for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                        for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                           const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_ut + x ) ];
                           SBB_singlet[ 0 ][ shift + x + ( uy * ( uy + 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy + 1 ) ) / 2 ) ] -= gamma_tx;
                        }
                     }
                  }
               } else {
                  for ( int uy = 0; uy < num_ut; uy++ ){
                     const double f_uu = fock->get( irrep_ut, nocc_ut + uy, nocc_ut + uy );
                     for ( int t = 0; t <= uy; t++ ){ // 0 <= t <= uy < num_ut
                        for ( int x = 0; x <= uy; x++ ){ // 0 <= x <= uy < num_ut
                           const double val_tx = ( f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ]
                                            - f_uu * one_rdm[ d_ut + t + LAS * ( d_ut + x ) ] );
                           FBB_singlet[ 0 ][ shift + x + ( uy * ( uy + 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy + 1 ) ) / 2 ) ] -= val_tx;
                        }
                     }
                  }
               }

               // SBB singlet: - delta_tx Gamma_{uy}
               // FBB singlet: - delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
               if ( OVLP ){
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                        for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                           const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_ut + y ) ];
                           SBB_singlet[ 0 ][ shift + tx + ( y * ( y + 1 ) ) / 2 + SIZE * ( tx + ( u * ( u + 1 ) ) / 2 ) ] -= gamma_uy;
                        }
                     }
                  }
               } else {
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + tx, nocc_ut + tx );
                     for ( int u = tx; u < num_ut; u++ ){ // 0 <= tx <= u < num_ut
                        for ( int y = tx; y < num_ut; y++ ){ // 0 <= tx <= y < num_xy = num_ut
                           const double val_uy = ( f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ]
                                            - f_tt * one_rdm[ d_ut + u + LAS * ( d_ut + y ) ] );
                           FBB_singlet[ 0 ][ shift + tx + ( y * ( y + 1 ) ) / 2 + SIZE * ( tx + ( u * ( u + 1 ) ) / 2 ) ] -= val_uy;
                        }
                     }
                  }
               }

               // SBB singlet: - delta_ux Gamma_{ty}
               // FBB singlet: - delta_ux ( f_dot_2dm[ ty ] - f_uu Gamma_{ty} )
               if ( OVLP ){
                  for ( int ux = 0; ux < num_ut; ux++ ){
                     for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                        for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                           const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_ut + y ) ];
                           SBB_singlet[ 0 ][ shift + ux + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux + 1 ) ) / 2 ) ] -= gamma_ty;
                        }
                     }
                  }
               } else {
                  for ( int ux = 0; ux < num_ut; ux++ ){
                     const double f_uu = fock->get( irrep_ut, nocc_ut + ux, nocc_ut + ux );
                     for ( int t = 0; t <= ux; t++ ){ // 0 <= t <= ux < num_ut
                        for ( int y = ux; y < num_ut; y++ ){ // 0 <= ux <= y < num_xy = num_ut
                           const double val_ty = ( f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ]
                                            - f_uu * one_rdm[ d_ut + t + LAS * ( d_ut + y ) ] );
                           FBB_singlet[ 0 ][ shift + ux + ( y * ( y + 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux + 1 ) ) / 2 ) ] -= val_ty;
                        }
                     }
                  }
               }

               // SBB singlet: - delta_ty Gamma_{ux}
               // FBB singlet: - delta_ty ( f_dot_2dm[ ux ] - f_tt Gamma_{ux} )
               if ( OVLP ){
                  for ( int ty = 0; ty < num_ut; ty++ ){
                     for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                        for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                           const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_ut + x ) ];
                           SBB_singlet[ 0 ][ shift + x + ( ty * ( ty + 1 ) ) / 2 + SIZE * ( ty + ( u * ( u + 1 ) ) / 2 ) ] -= gamma_ux;
                        }
                     }
                  }
               } else {
                  for ( int ty = 0; ty < num_ut; ty++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + ty, nocc_ut + ty );
                     for ( int u = ty; u < num_ut; u++ ){ // 0 <= ty <= u < num_ut
                        for ( int x = 0; x <= ty; x++ ){ // 0 <= x <= ty < num_ut
                           const double val_ux = ( f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ]
                                            - f_tt * one_rdm[ d_ut + u + LAS * ( d_ut + x ) ] );
                           FBB_singlet[ 0 ][ shift + x + ( ty * ( ty + 1 ) ) / 2 + SIZE * ( ty + ( u * ( u + 1 ) ) / 2 ) ] -= val_ux;
                        }
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy + 1 ) ) / 2;
         }
         if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
            // B: E_ti E_uj | 0 >   --->   tu: excitation into
            // F: E_at E_bu | 0 >   --->   tu: excitation out of
            for ( int u = 0; u < num_ut; u++ ){
               const double gamma_uu = one_rdm[ ( d_ut + u ) * ( 1 + LAS ) ];
               for ( int t = 0; t <= u; t++ ){ // 0 <= t <= u < num_ut
                  const double gamma_tt = one_rdm[ ( d_ut + t ) * ( 1 + LAS ) ];
                  const double ipea_B_tu = 0.5 * IPEA * ( gamma_tt + gamma_uu );
                  const double ipea_F_tu = 0.5 * IPEA * ( 4.0 - gamma_tt - gamma_uu );
                  const int ptr = ( jump_col + t + ( u * ( u + 1 ) ) / 2 ) * ( 1 + SIZE );
                  FBB_singlet[ 0 ][ ptr ] += ipea_B_tu * SBB_singlet[ 0 ][ ptr ];
                  FFF_singlet[ 0 ][ ptr ] += ipea_F_tu * SFF_singlet[ 0 ][ ptr ];
               }
            }
         }
         jump_col += ( num_ut * ( num_ut + 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_singlet[ irrep ] == size_F_singlet[ irrep ] ); // At construction
      const int SIZE = size_B_singlet[ irrep ];
      if ( OVLP ){ SBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
                   SFF_singlet[ irrep ] = new double[ SIZE * SIZE ]; }
      else {       FBB_singlet[ irrep ] = new double[ SIZE * SIZE ];
                   FFF_singlet[ irrep ] = new double[ SIZE * SIZE ]; }

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

                  // SBB singlet: + Gamma_{utyx} + Gamma_{utxy}
                  // SFF singlet: + Gamma_{yxut} + Gamma_{yxtu}
                  // FBB singlet: + f_dot_3dm[ utyx ] + f_dot_3dm[ tuyx ] + ( f_tt + f_uu + f_xx + f_yy ) SBB_singlet[ xytu ]
                  // FFF singlet: + f_dot_3dm[ yxut ] + f_dot_3dm[ yxtu ]
                  if ( OVLP ){
                     for ( int t = 0; t < num_t; t++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double value = ( two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                      + two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                                 const int ptr = shift + x + num_x * y + SIZE * ( t + num_t * u );
                                 SBB_singlet[ irrep ][ ptr ] = value;
                                 SFF_singlet[ irrep ][ ptr ] = value;
                              }
                           }
                        }
                     }
                  } else {
                     for ( int t = 0; t < num_t; t++ ){
                        const double f_tt = fock->get( irrep_t, nocc_t + t, nocc_t + t );
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_uu = fock->get( irrep_u, nocc_u + u, nocc_u + u );
                           for ( int x = 0; x < num_x; x++ ){
                              const double f_xx = fock->get( irrep_x, nocc_x + x, nocc_x + x );
                              for ( int y = 0; y < num_y; y++ ){
                                 const double f_yy = fock->get( irrep_y, nocc_y + y, nocc_y + y );
                                 const int ptr = shift + x + num_x * y + SIZE * ( t + num_t * u );
                                 const double fdotsum = ( f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                        + f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                                 FBB_singlet[ irrep ][ ptr ] = fdotsum + ( f_xx + f_yy + f_tt + f_uu ) * SBB_singlet[ irrep ][ ptr ];
                                 FFF_singlet[ irrep ][ ptr ] = fdotsum;
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // num_t == num_x  and  num_u == num_y

                     // SBB singlet: + 2 delta_uy delta_tx
                     // FBB singlet: + 2 delta_uy delta_tx ( f_dot_1dm - f_tt - f_uu )
                     if ( OVLP ){
                        for ( int xt = 0; xt < num_x; xt++ ){
                           for ( int yu = 0; yu < num_y; yu++ ){
                              SBB_singlet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 2.0;
                           }
                        }
                     } else {
                        for ( int xt = 0; xt < num_x; xt++ ){
                           const double f_tt = fock->get( irrep_x, nocc_x + xt, nocc_x + xt );
                           for ( int yu = 0; yu < num_y; yu++ ){
                              const double f_uu = fock->get( irrep_y, nocc_y + yu, nocc_y + yu );
                              FBB_singlet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 2 * ( f_dot_1dm - f_tt - f_uu );
                           }
                        }
                     }

                     // SBB singlet: - delta_tx Gamma_{uy}
                     // FBB singlet: - delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
                     if ( OVLP ){
                        for ( int xt = 0; xt < num_x; xt++ ){
                           for ( int u = 0; u < num_y; u++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double gamma_uy = one_rdm[ d_u + u + LAS * ( d_u + y ) ];
                                 SBB_singlet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= gamma_uy;
                              }
                           }
                        }
                     } else {
                        for ( int xt = 0; xt < num_x; xt++ ){
                           const double f_tt = fock->get( irrep_x, nocc_x + xt, nocc_x + xt );
                           for ( int u = 0; u < num_y; u++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double val_uy = ( f_dot_2dm[ d_u + u + LAS * ( d_u + y ) ]
                                                  - f_tt * one_rdm[ d_u + u + LAS * ( d_u + y ) ] );
                                 FBB_singlet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= val_uy;
                              }
                           }
                        }
                     }

                     // SBB singlet: - delta_uy Gamma_{tx}
                     // FBB singlet: - delta_uy ( f_dot_2dm[ tx ] - f_uu Gamma_{tx} )
                     if ( OVLP ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           for ( int t = 0; t < num_x; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 const double gamma_tx = one_rdm[ d_t + t + LAS * ( d_t + x ) ];
                                 SBB_singlet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_t * yu ) ] -= gamma_tx;
                              }
                           }
                        }
                     } else {
                        for ( int yu = 0; yu < num_y; yu++ ){
                           const double f_uu = fock->get( irrep_y, nocc_y + yu, nocc_y + yu );
                           for ( int t = 0; t < num_x; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 const double val_tx = ( f_dot_2dm[ d_t + t + LAS * ( d_t + x ) ]
                                                  - f_uu * one_rdm[ d_t + t + LAS * ( d_t + x ) ] );
                                 FBB_singlet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_t * yu ) ] -= val_tx;
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
               // B: E_ti E_uj | 0 >   --->   tu: excitation into
               // F: E_at E_bu | 0 >   --->   tu: excitation out of
               for ( int u = 0; u < num_u; u++ ){
                  const double gamma_uu = one_rdm[ ( d_u + u ) * ( 1 + LAS ) ];
                  for ( int t = 0; t < num_t; t++ ){
                     const double gamma_tt = one_rdm[ ( d_t + t ) * ( 1 + LAS ) ];
                     const double ipea_B_tu = 0.5 * IPEA * ( gamma_tt + gamma_uu );
                     const double ipea_F_tu = 0.5 * IPEA * ( 4.0 - gamma_tt - gamma_uu );
                     const int ptr = ( jump_col + t + num_t * u ) * ( 1 + SIZE );
                     FBB_singlet[ irrep ][ ptr ] += ipea_B_tu * SBB_singlet[ irrep ][ ptr ];
                     FFF_singlet[ irrep ][ ptr ] += ipea_F_tu * SFF_singlet[ irrep ][ ptr ];
                  }
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_BB_FF_triplet( const bool OVLP, const double IPEA ){

   /*
      | TB_tiuj > = ( E_ti E_uj - E_tj E_ui ) / sqrt( 1 + delta_ij ) | 0 >  with  i < j and t < u

      SBB triplet: < TB_xkyl | 1 | TB_tiuj > = 2 delta_ik delta_jl ( SBB_triplet[ Itu ][ xytu ] )
      FBB triplet: < TB_xkyl | F | TB_tiuj > = 2 delta_ik delta_jl ( FBB_triplet[ Itu ][ xytu ] + ( 2 sum_n f_nn - f_ii - f_jj ) * SBB_triplet[ Itu ][ xytu ] )

            SBB_triplet[ Itu ][ xytu ] = ( + Gamma_{utyx}
                                           - Gamma_{utxy}
                                           + 6 delta_uy delta_tx
                                           - 6 delta_ux delta_ty
                                           - 3 delta_uy Gamma_{tx}
                                           - 3 delta_tx Gamma_{uy}
                                           + 3 delta_ux Gamma_{ty}
                                           + 3 delta_ty Gamma_{ux}
                                         )

            FBB_triplet[ Itu ][ xytu ] = ( + f_dot_3dm[ utyx ]
                                           - f_dot_3dm[ tuyx ]
                                           + ( f_tt + f_uu + f_xx + f_yy ) SBB_triplet[ xytu ]
                                           + 6 delta_uy delta_tx ( f_dot_1dm - f_tt - f_uu )
                                           - 6 delta_ux delta_ty ( f_dot_1dm - f_tt - f_uu )
                                           - 3 delta_uy ( f_dot_2dm[ tx ] - f_uu Gamma_{tx} )
                                           - 3 delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
                                           + 3 delta_ux ( f_dot_2dm[ ty ] - f_uu Gamma_{ty} )
                                           + 3 delta_ty ( f_dot_2dm[ ux ] - f_tt Gamma_{ux} )
                                         )

      | TF_atbu > = ( E_at E_bu - E_bt E_au ) / sqrt( 1 + delta_ab ) | 0 >  with  a < b and t < u

      SFF triplet: < TF_cxdy | 1 | TF_atbu > = 2 delta_ac delta_bd ( SFF_triplet[ Itu ][ xytu ] )
      FFF triplet: < TF_cxdy | F | TF_atbu > = 2 delta_ac delta_bd ( FFF_triplet[ Itu ][ xytu ] + ( 2 sum_n f_nn + f_aa + f_bb ) * SFF_triplet[ Itu ][ xytu ] )

            SFF_triplet[ Itu ][ xytu ] = ( + Gamma_{yxut}
                                           - Gamma_{yxtu}
                                         )

            FFF_triplet[ Itu ][ xytu ] = ( + f_dot_3dm[ yxut ]
                                           - f_dot_3dm[ yxtu ]
                                         )
   */

   if ( OVLP ){ SBB_triplet = new double*[ num_irreps ];
                SFF_triplet = new double*[ num_irreps ]; }
   else {       FBB_triplet = new double*[ num_irreps ];
                FFF_triplet = new double*[ num_irreps ]; }

   const int LAS = indices->getDMRGcumulative( num_irreps );

   { // First do irrep == Iia x Ijb == 0 -->  It == Iu and Ix == Iy
      assert( size_B_triplet[ 0 ] == size_F_triplet[ 0 ] ); // At construction
      const int SIZE = size_B_triplet[ 0 ];
      if ( OVLP ){ SBB_triplet[ 0 ] = new double[ SIZE * SIZE ];
                   SFF_triplet[ 0 ] = new double[ SIZE * SIZE ]; }
      else {       FBB_triplet[ 0 ] = new double[ SIZE * SIZE ];
                   FFF_triplet[ 0 ] = new double[ SIZE * SIZE ]; }

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

            // SBB triplet: + Gamma_{utyx} - Gamma_{utxy}
            // SFF triplet: + Gamma_{yxut} - Gamma_{yxtu}
            // FBB triplet: + f_dot_3dm[ utyx ] - f_dot_3dm[ tuyx ] + ( f_tt + f_uu + f_xx + f_yy ) SBB_triplet[ xytu ]
            // FFF triplet: + f_dot_3dm[ yxut ] - f_dot_3dm[ yxtu ]
            if ( OVLP ){
               for ( int t = 0; t < num_ut; t++ ){
                  for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                     for ( int x = 0; x < num_xy; x++ ){
                        for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                           const int ptr = shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 );
                           const double value = ( two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                                - two_rdm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                           SBB_triplet[ 0 ][ ptr ] = value;
                           SFF_triplet[ 0 ][ ptr ] = value;
                        }
                     }
                  }
               }
            } else {
               for ( int t = 0; t < num_ut; t++ ){
                  const double f_tt = fock->get( irrep_ut, nocc_ut + t, nocc_ut + t );
                  for ( int u = t+1; u < num_ut; u++ ){ // 0 <= t < u < num_ut
                     const double f_uu = fock->get( irrep_ut, nocc_ut + u, nocc_ut + u );
                     for ( int x = 0; x < num_xy; x++ ){
                        const double f_xx = fock->get( irrep_xy, nocc_xy + x, nocc_xy + x );
                        for ( int y = x+1; y < num_xy; y++ ){ // 0 <= x < y < num_xy
                           const double f_yy = fock->get( irrep_xy, nocc_xy + y, nocc_xy + y );
                           const int ptr = shift + x + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( u * ( u - 1 ) ) / 2 );
                           const double fdotdiff = ( f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + t + LAS * ( d_ut + u ))) ]
                                                   - f_dot_3dm[ d_xy + x + LAS * ( d_xy + y + LAS * ( d_ut + u + LAS * ( d_ut + t ))) ] );
                           FBB_triplet[ 0 ][ ptr ] = fdotdiff + ( f_tt + f_uu + f_xx + f_yy ) * SBB_triplet[ 0 ][ ptr ];
                           FFF_triplet[ 0 ][ ptr ] = fdotdiff;
                        }
                     }
                  }
               }
            }

            if ( irrep_ut == irrep_xy ){ // All four irreps are equal --> num_ut == num_xy

               // SBB triplet: + 6 ( delta_uy delta_tx - delta_ux delta_ty )
               // FBB triplet: + 6 ( delta_uy delta_tx - delta_ux delta_ty ) ( f_dot_1dm - f_tt - f_uu )
               if ( OVLP ){
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     for ( int uy = tx+1; uy < num_ut; uy++ ){
                        SBB_triplet[ 0 ][ shift + tx + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( tx + ( uy * ( uy - 1 ) ) / 2 ) ] += 6.0;
                     }
                  }
               } else {
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + tx, nocc_ut + tx );
                     for ( int uy = tx+1; uy < num_ut; uy++ ){
                        const double f_uu = fock->get( irrep_ut, nocc_ut + uy, nocc_ut + uy );
                        FBB_triplet[ 0 ][ shift + tx + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( tx + ( uy * ( uy - 1 ) ) / 2 ) ] += 6 * ( f_dot_1dm - f_tt - f_uu );
                     }
                  }
               }

               // SBB triplet: - 3 delta_uy Gamma_{tx}
               // FBB triplet: - 3 delta_uy ( f_dot_2dm[ tx ] - f_yu Gamma_{tx} )
               if ( OVLP ){
                  for ( int uy = 0; uy < num_ut; uy++ ){
                     for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                        for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                           const double gamma_tx = one_rdm[ d_ut + t + LAS * ( d_xy + x ) ];
                           SBB_triplet[ 0 ][ shift + x + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy - 1 ) ) / 2 ) ] -= 3 * gamma_tx;
                        }
                     }
                  }
               } else {
                  for ( int uy = 0; uy < num_ut; uy++ ){
                     const double f_uu = fock->get( irrep_ut, nocc_ut + uy, nocc_ut + uy );
                     for ( int t = 0; t < uy; t++ ){ // 0 <= t < uy < num_ut
                        for ( int x = 0; x < uy; x++ ){ // 0 <= x < uy < num_ut
                           const double val_tx = ( f_dot_2dm[ d_ut + t + LAS * ( d_ut + x ) ]
                                            - f_uu * one_rdm[ d_ut + t + LAS * ( d_ut + x ) ] );
                           FBB_triplet[ 0 ][ shift + x + ( uy * ( uy - 1 ) ) / 2 + SIZE * ( t + ( uy * ( uy - 1 ) ) / 2 ) ] -= 3 * val_tx;
                        }
                     }
                  }
               }

               // SBB triplet: - 3 delta_tx Gamma_{uy}
               // FBB triplet: - 3 delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
               if ( OVLP ){
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                        for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                           const double gamma_uy = one_rdm[ d_ut + u + LAS * ( d_xy + y ) ];
                           SBB_triplet[ 0 ][ shift + tx + ( y * ( y - 1 ) ) / 2 + SIZE * ( tx + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * gamma_uy;
                        }
                     }
                  }
               } else {
                  for ( int tx = 0; tx < num_ut; tx++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + tx, nocc_ut + tx );
                     for ( int u = tx+1; u < num_ut; u++ ){ // 0 <= tx < u < num_ut
                        for ( int y = tx+1; y < num_ut; y++ ){ // 0 <= tx < y < num_xy = num_ut
                           const double val_uy = ( f_dot_2dm[ d_ut + u + LAS * ( d_ut + y ) ]
                                            - f_tt * one_rdm[ d_ut + u + LAS * ( d_ut + y ) ] );
                           FBB_triplet[ 0 ][ shift + tx + ( y * ( y - 1 ) ) / 2 + SIZE * ( tx + ( u * ( u - 1 ) ) / 2 ) ] -= 3 * val_uy;
                        }
                     }
                  }
               }

               // SBB triplet: + 3 delta_ux Gamma_{ty}
               // FBB triplet: + 3 delta_ux ( f_dot_2dm[ ty ] - f_uu Gamma_{ty} )
               if ( OVLP ){
                  for ( int ux = 0; ux < num_ut; ux++ ){
                     for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                        for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                           const double gamma_ty = one_rdm[ d_ut + t + LAS * ( d_xy + y ) ];
                           SBB_triplet[ 0 ][ shift + ux + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux - 1 ) ) / 2 ) ] += 3 * gamma_ty;
                        }
                     }
                  }
               } else {
                  for ( int ux = 0; ux < num_ut; ux++ ){
                     const double f_uu = fock->get( irrep_ut, nocc_ut + ux, nocc_ut + ux );
                     for ( int t = 0; t < ux; t++ ){ // 0 <= t < ux < num_ut
                        for ( int y = ux+1; y < num_ut; y++ ){ // 0 <= ux < y < num_xy = num_ut
                           const double val_ty = ( f_dot_2dm[ d_ut + t + LAS * ( d_ut + y ) ]
                                            - f_uu * one_rdm[ d_ut + t + LAS * ( d_ut + y ) ] );
                           FBB_triplet[ 0 ][ shift + ux + ( y * ( y - 1 ) ) / 2 + SIZE * ( t + ( ux * ( ux - 1 ) ) / 2 ) ] += 3 * val_ty;
                        }
                     }
                  }
               }

               // SBB triplet: + 3 delta_ty Gamma_{ux}
               // FBB triplet: + 3 delta_ty ( f_dot_2dm[ ux ] - f_tt Gamma_{ux} )
               if ( OVLP ){
                  for ( int ty = 0; ty < num_ut; ty++ ){
                     for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                        for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                           const double gamma_ux = one_rdm[ d_ut + u + LAS * ( d_xy + x ) ];
                           SBB_triplet[ 0 ][ shift + x + ( ty * ( ty - 1 ) ) / 2 + SIZE * ( ty + ( u * ( u - 1 ) ) / 2 ) ] += 3 * gamma_ux;
                        }
                     }
                  }
               } else {
                  for ( int ty = 0; ty < num_ut; ty++ ){
                     const double f_tt = fock->get( irrep_ut, nocc_ut + ty, nocc_ut + ty );
                     for ( int u = ty+1; u < num_ut; u++ ){ // 0 <= ty < u < num_ut
                        for ( int x = 0; x < ty; x++ ){ // 0 <= x < ty < num_ut
                           const double val_ux = ( f_dot_2dm[ d_ut + u + LAS * ( d_ut + x ) ]
                                            - f_tt * one_rdm[ d_ut + u + LAS * ( d_ut + x ) ] );
                           FBB_triplet[ 0 ][ shift + x + ( ty * ( ty - 1 ) ) / 2 + SIZE * ( ty + ( u * ( u - 1 ) ) / 2 ) ] += 3 * val_ux;
                        }
                     }
                  }
               }
            }
            jump_row += ( num_xy * ( num_xy - 1 ) ) / 2;
         }
         if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
            // B: E_ti E_uj | 0 >   --->   tu: excitation into
            // F: E_at E_bu | 0 >   --->   tu: excitation out of
            for ( int u = 0; u < num_ut; u++ ){
               const double gamma_uu = one_rdm[ ( d_ut + u ) * ( 1 + LAS ) ];
               for ( int t = 0; t < u; t++ ){ // 0 <= t < u < num_ut
                  const double gamma_tt = one_rdm[ ( d_ut + t ) * ( 1 + LAS ) ];
                  const double ipea_B_tu = 0.5 * IPEA * ( gamma_tt + gamma_uu );
                  const double ipea_F_tu = 0.5 * IPEA * ( 4.0 - gamma_tt - gamma_uu );
                  const int ptr = ( jump_col + t + ( u * ( u - 1 ) ) / 2 ) * ( 1 + SIZE );
                  FBB_triplet[ 0 ][ ptr ] += ipea_B_tu * SBB_triplet[ 0 ][ ptr ];
                  FFF_triplet[ 0 ][ ptr ] += ipea_F_tu * SFF_triplet[ 0 ][ ptr ];
               }
            }
         }
         jump_col += ( num_ut * ( num_ut - 1 ) ) / 2;
      }
   }

   for ( int irrep = 1; irrep < num_irreps; irrep++ ){ // Then do irrep == Iia x Ijb != 0 -->  It != Iu and Ix != Iy

      assert( size_B_triplet[ irrep ] == size_F_triplet[ irrep ] ); // At construction
      const int SIZE = size_B_triplet[ irrep ];
      if ( OVLP ){ SBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
                   SFF_triplet[ irrep ] = new double[ SIZE * SIZE ]; }
      else {       FBB_triplet[ irrep ] = new double[ SIZE * SIZE ];
                   FFF_triplet[ irrep ] = new double[ SIZE * SIZE ]; }

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

                  // SBB triplet: + Gamma_{utyx} - Gamma_{utxy}
                  // SFF triplet: + Gamma_{yxut} - Gamma_{yxtu}
                  // FBB triplet: + f_dot_3dm[ utyx ] - f_dot_3dm[ tuyx ] + ( f_tt + f_uu + f_xx + f_yy ) SBB_triplet[ xytu ]
                  // FFF triplet: + f_dot_3dm[ yxut ] - f_dot_3dm[ yxtu ]
                  if ( OVLP ){
                     for ( int t = 0; t < num_t; t++ ){
                        for ( int u = 0; u < num_u; u++ ){
                           for ( int x = 0; x < num_x; x++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const int ptr = shift + x + num_x * y + SIZE * ( t + num_t * u );
                                 const double value = ( two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                      - two_rdm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                                 SBB_triplet[ irrep ][ ptr ] = value;
                                 SFF_triplet[ irrep ][ ptr ] = value;
                              }
                           }
                        }
                     }
                  } else {
                     for ( int t = 0; t < num_t; t++ ){
                        const double f_tt = fock->get( irrep_t, nocc_t + t, nocc_t + t );
                        for ( int u = 0; u < num_u; u++ ){
                           const double f_uu = fock->get( irrep_u, nocc_u + u, nocc_u + u );
                           for ( int x = 0; x < num_x; x++ ){
                              const double f_xx = fock->get( irrep_x, nocc_x + x, nocc_x + x );
                              for ( int y = 0; y < num_y; y++ ){
                                 const double f_yy = fock->get( irrep_y, nocc_y + y, nocc_y + y );
                                 const int ptr = shift + x + num_x * y + SIZE * ( t + num_t * u );
                                 const double fdotdiff = ( f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_t + t + LAS * ( d_u + u ))) ]
                                                         - f_dot_3dm[ d_x + x + LAS * ( d_y + y + LAS * ( d_u + u + LAS * ( d_t + t ))) ] );
                                 FBB_triplet[ irrep ][ ptr ] = fdotdiff + ( f_tt + f_uu + f_xx + f_yy ) * SBB_triplet[ irrep ][ ptr ];
                                 FFF_triplet[ irrep ][ ptr ] = fdotdiff;
                              }
                           }
                        }
                     }
                  }

                  if (( irrep_u == irrep_y ) && ( irrep_t == irrep_x )){ // --> num_t == num_x   and   num_u == num_y

                     // SBB triplet: + 6 delta_uy delta_tx
                     // FBB triplet: + 6 delta_uy delta_tx ( f_dot_1dm - f_tt - f_uu )
                     if ( OVLP ){
                        for ( int xt = 0; xt < num_x; xt++ ){
                           for ( int yu = 0; yu < num_y; yu++ ){
                              SBB_triplet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 6.0;
                           }
                        }
                     } else {
                        for ( int xt = 0; xt < num_x; xt++ ){
                           const double f_tt = fock->get( irrep_x, nocc_x + xt, nocc_x + xt );
                           for ( int yu = 0; yu < num_y; yu++ ){
                              const double f_uu = fock->get( irrep_y, nocc_y + yu, nocc_y + yu );
                              FBB_triplet[ irrep ][ shift + xt + num_x * yu + SIZE * ( xt + num_x * yu ) ] += 6 * ( f_dot_1dm - f_tt - f_uu );
                           }
                        }
                     }

                     // SBB triplet: - 3 delta_tx Gamma_{uy}
                     // FBB triplet: - 3 delta_tx ( f_dot_2dm[ uy ] - f_tt Gamma_{uy} )
                     if ( OVLP ){
                        for ( int xt = 0; xt < num_x; xt++ ){
                           for ( int u = 0; u < num_y; u++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double gamma_uy = one_rdm[ d_u + u + LAS * ( d_y + y ) ];
                                 SBB_triplet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= 3 * gamma_uy;
                              }
                           }
                        }
                     } else {
                        for ( int xt = 0; xt < num_x; xt++ ){
                           const double f_tt = fock->get( irrep_x, nocc_x + xt, nocc_x + xt );
                           for ( int u = 0; u < num_y; u++ ){
                              for ( int y = 0; y < num_y; y++ ){
                                 const double val_uy = ( f_dot_2dm[ d_u + u + LAS * ( d_u + y ) ]
                                                  - f_tt * one_rdm[ d_u + u + LAS * ( d_u + y ) ] );
                                 FBB_triplet[ irrep ][ shift + xt + num_x * y + SIZE * ( xt + num_x * u ) ] -= 3 * val_uy;
                              }
                           }
                        }
                     }

                     // SBB triplet: - 3 delta_uy Gamma_{tx}
                     // FBB triplet: - 3 delta_uy ( f_dot_2dm[ tx ] - f_uu Gamma_{tx} )
                     if ( OVLP ){
                        for ( int yu = 0; yu < num_y; yu++ ){
                           for ( int t = 0; t < num_x; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 const double gamma_tx = one_rdm[ d_t + t + LAS * ( d_x + x ) ];
                                 SBB_triplet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_x * yu ) ] -= 3 * gamma_tx;
                              }
                           }
                        }
                     } else {
                        for ( int yu = 0; yu < num_y; yu++ ){
                           const double f_uu = fock->get( irrep_y, nocc_y + yu, nocc_y + yu );
                           for ( int t = 0; t < num_x; t++ ){
                              for ( int x = 0; x < num_x; x++ ){
                                 const double val_tx = ( f_dot_2dm[ d_t + t + LAS * ( d_t + x ) ]
                                                  - f_uu * one_rdm[ d_t + t + LAS * ( d_t + x ) ] );
                                 FBB_triplet[ irrep ][ shift + x + num_x * yu + SIZE * ( t + num_x * yu ) ] -= 3 * val_tx;
                              }
                           }
                        }
                     }
                  }
                  jump_row += num_x * num_y;
               }
            }
            if (( OVLP == false ) && ( fabs( IPEA ) > 0.0 )){
               // B: E_ti E_uj | 0 >   --->   tu: excitation into
               // F: E_at E_bu | 0 >   --->   tu: excitation out of
               for ( int u = 0; u < num_u; u++ ){
                  const double gamma_uu = one_rdm[ ( d_u + u ) * ( 1 + LAS ) ];
                  for ( int t = 0; t < num_t; t++ ){
                     const double gamma_tt = one_rdm[ ( d_t + t ) * ( 1 + LAS ) ];
                     const double ipea_B_tu = 0.5 * IPEA * ( gamma_tt + gamma_uu );
                     const double ipea_F_tu = 0.5 * IPEA * ( 4.0 - gamma_tt - gamma_uu );
                     const int ptr = ( jump_col + t + num_t * u ) * ( 1 + SIZE );
                     FBB_triplet[ irrep ][ ptr ] += ipea_B_tu * SBB_triplet[ irrep ][ ptr ];
                     FFF_triplet[ irrep ][ ptr ] += ipea_F_tu * SFF_triplet[ irrep ][ ptr ];
                  }
               }
            }
            jump_col += num_t * num_u;
         }
      }
   }

}

void CheMPS2::CASPT2::make_EE_GG( const bool OVLP, const double IPEA ){

   /*
      | SE_tiaj > = ( E_ti E_aj + E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <= j
      | TE_tiaj > = ( E_ti E_aj - E_tj E_ai ) / sqrt( 1 + delta_ij ) | 0 >  with  i <  j

         < SE_ukbl | 1 | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( SEE[ It ][ ut ] )
         < SE_ukbl | F | SE_tiaj > = 2 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )
         < TE_ukbl | 1 | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( SEE[ It ][ ut ] )
         < TE_ukbl | F | TE_tiaj > = 6 delta_ab delta_ik delta_jl ( FEE[ It ][ ut ] + ( 2 sum_k f_kk + f_aa - f_ii - f_jj ) SEE[ It ][ ut ] )

            SEE[ It ][ ut ] = ( + 2 delta_tu
                                - Gamma_tu
                              )

            FEE[ It ][ ut ] = ( + 2 * delta_ut * ( f_dot_1dm + f_tt )
                                - f_dot_2dm[ It ][ tu ] - ( f_tt + f_uu ) Gamma_ut
                              )

      | SG_aibt > = ( E_ai E_bt + E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <= b
      | TG_aibt > = ( E_ai E_bt - E_bi E_at ) / sqrt( 1 + delta_ab ) | 0 >  with  a <  b

         < SG_cjdu | 1 | SG_aibt > = 2 delta_ij delta_ac delta_bd ( SGG[ It ][ ut ] )
         < SG_cjdu | F | SG_aibt > = 2 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )
         < TG_cjdu | 1 | TG_aibt > = 6 delta_ij delta_ac delta_bd ( SGG[ It ][ ut ] )
         < TG_cjdu | F | TG_aibt > = 6 delta_ij delta_ac delta_bd ( FGG[ It ][ ut ] + ( 2 sum_k f_kk + f_aa + f_bb - f_ii ) SGG[ It ][ ut ] )

            SGG[ It ][ ut ] = ( + Gamma_ut
                              )

            FGG[ It ][ ut ] = ( + f_dot_2dm[ It ][ ut ]
                              )
   */

   if ( OVLP ){ SEE = new double*[ num_irreps ];
                SGG = new double*[ num_irreps ]; }
   else {       FEE = new double*[ num_irreps ];
                FGG = new double*[ num_irreps ]; }

   const int LAS = indices->getDMRGcumulative( num_irreps );

   for ( int irrep_ut = 0; irrep_ut < num_irreps; irrep_ut++ ){
      const int d_ut  = indices->getDMRGcumulative( irrep_ut );
      const int SIZE  = indices->getNDMRG( irrep_ut );
      const int NOCC  = indices->getNOCC( irrep_ut );
      if ( OVLP ){
         SEE[ irrep_ut ] = new double[ SIZE * SIZE ];
         SGG[ irrep_ut ] = new double[ SIZE * SIZE ];
         for ( int t = 0; t < SIZE; t++ ){
            for ( int u = 0; u < SIZE; u++ ){
               const double gamma_ut = one_rdm[ d_ut + u + LAS * ( d_ut + t ) ];
               SEE[ irrep_ut ][ u + SIZE * t ] = - gamma_ut;
               SGG[ irrep_ut ][ u + SIZE * t ] =   gamma_ut;
            }
            SEE[ irrep_ut ][ t + SIZE * t ] += 2.0;
         }
      } else {
         FEE[ irrep_ut ] = new double[ SIZE * SIZE ];
         FGG[ irrep_ut ] = new double[ SIZE * SIZE ];
         for ( int t = 0; t < SIZE; t++ ){
            const double f_tt = fock->get( irrep_ut, NOCC + t, NOCC + t );
            for ( int u = 0; u < SIZE; u++ ){
               const double f_uu = fock->get( irrep_ut, NOCC + u, NOCC + u );
               const double gamma_ut =   one_rdm[ d_ut + u + LAS * ( d_ut + t ) ];
               const double fdot2_ut = f_dot_2dm[ d_ut + u + LAS * ( d_ut + t ) ];
               FEE[ irrep_ut ][ u + SIZE * t ] = - fdot2_ut - ( f_tt + f_uu ) * gamma_ut;
               FGG[ irrep_ut ][ u + SIZE * t ] =   fdot2_ut;
            }
            FEE[ irrep_ut ][ t + SIZE * t ] += 2 * ( f_dot_1dm + f_tt );
         }
         if ( fabs( IPEA ) > 0.0 ){
            // E: E_ti E_aj | 0 >   --->   t: excitation into
            // G: E_ai E_bt | 0 >   --->   t: excitation out of
            for ( int t = 0; t < SIZE; t++ ){
               const double gamma_tt = one_rdm[ ( d_ut + t ) * ( 1 + LAS ) ];
               const double ipea_E_t = 0.5 * IPEA * ( gamma_tt );
               const double ipea_G_t = 0.5 * IPEA * ( 2.0 - gamma_tt );
               const int ptr = t * ( 1 + SIZE );
               FEE[ irrep_ut ][ ptr ] += ipea_E_t * SEE[ irrep_ut ][ ptr ];
               FGG[ irrep_ut ][ ptr ] += ipea_G_t * SGG[ irrep_ut ][ ptr ];
            }
         }
      }
   }

}


