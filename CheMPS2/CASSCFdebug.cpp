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

#include "CASSCF.h"

using std::cout;
using std::endl;

void CheMPS2::CASSCF::checkHF( int * docc, int * socc ){

   double EnergyHF = NUCL_ORIG;

   cout << "Single particle energy levels : " << endl;
   for ( int irrep = 0; irrep < num_irreps; irrep++ ){
      for ( int orb = 0; orb < iHandler->getNORB( irrep ); orb++ ){

         double SPenergy = TMAT_ORIG->get( irrep, orb, orb );

         const int num_beta  = (( orb < docc[ irrep ] ) ? 1 : 0 );
         const int num_alpha = (( orb < docc[ irrep ] + socc[ irrep ] ) ? 1 : 0 );
         const int num_total = num_alpha + num_beta;

         EnergyHF += num_total * TMAT_ORIG->get( irrep, orb, orb );

         for ( int irrep2 = 0; irrep2 < num_irreps; irrep2++ ){
            for ( int orb2 = 0; orb2 < iHandler->getNORB( irrep2 ); orb2++ ){

               const int num_beta2  = (( orb2 < docc[ irrep2 ] ) ? 1 : 0 );
               const int num_alpha2 = (( orb2 < docc[ irrep2 ] + socc[ irrep2 ] ) ? 1 : 0 );
               const int num_total2 = num_alpha2 + num_beta2;

               SPenergy += ( num_total2 * VMAT_ORIG->get( irrep, irrep2, irrep, irrep2, orb, orb2, orb, orb2 )
                           - num_alpha2 * VMAT_ORIG->get( irrep, irrep, irrep2, irrep2, orb, orb, orb2, orb2 ) );

               EnergyHF += 0.5 * num_total * num_total2 * VMAT_ORIG->get( irrep, irrep2, irrep, irrep2, orb, orb2, orb, orb2 );
               EnergyHF -= 0.5 * ( num_alpha * num_alpha2 + num_beta * num_beta2 ) * VMAT_ORIG->get( irrep, irrep, irrep2, irrep2, orb, orb, orb2, orb2 );

            }
         }
         cout << "   Orb " << iHandler->getOrigNOCCstart( irrep ) + orb << " : "<< orb + 1 << SymmInfo.getIrrepName(irrep) << " = " << SPenergy << endl;
      }
   }
   cout << "HF energy = " << EnergyHF << endl;

}

void CheMPS2::CASSCF::coeff_fe2( DMRG * theDMRG ){

   /*
      For the iron dimer:
         int NOCC[]  = {  5,  0,  2,  2,  0,  5,  2,  2 }; // 36 core electrons
         int NDMRG[] = {  6,  2,  3,  3,  2,  6,  3,  3 }; // 16 active space electrons, 28 active space orbitals
   */

   assert( nOrbDMRG == 28 );

   //               Ag                   B1g      B2g         B3g         Au       B1u                  B2u         B3u
   int coeff0[] = { 2, 2, 1, 0, 0, 0,    1, 0,    1, 0, 0,    1, 0, 0,    1, 0,    1, 1, 1, 0, 0, 0,    2, 0, 0,    2, 0, 0 };
   int coeff1[] = { 2, 1, 1, 0, 0, 0,    1, 0,    2, 0, 0,    1, 0, 0,    1, 0,    2, 1, 1, 0, 0, 0,    2, 0, 0,    1, 0, 0 };
   int coeff2[] = { 2, 1, 1, 0, 0, 0,    1, 0,    1, 0, 0,    2, 0, 0,    1, 0,    2, 1, 1, 0, 0, 0,    1, 0, 0,    2, 0, 0 };

   int coeff3[] = { 2, 2, 1, 0, 0, 0,    2, 0,    1, 0, 0,    1, 0, 0,    1, 0,    1, 1, 0, 0, 0, 0,    2, 0, 0,    2, 0, 0 };
   int coeff4[] = { 2, 1, 1, 0, 0, 0,    2, 0,    2, 0, 0,    1, 0, 0,    1, 0,    2, 1, 0, 0, 0, 0,    2, 0, 0,    1, 0, 0 };
   int coeff5[] = { 2, 1, 1, 0, 0, 0,    2, 0,    1, 0, 0,    2, 0, 0,    1, 0,    2, 1, 0, 0, 0, 0,    1, 0, 0,    2, 0, 0 };

   int coeff6[] = { 2, 2, 1, 0, 0, 0,    1, 0,    1, 0, 0,    1, 0, 0,    1, 0,    2, 1, 1, 0, 0, 0,    2, 0, 0,    2, 0, 0 };
   int coeff7[] = { 2, 2, 1, 0, 0, 0,    1, 0,    1, 0, 0,    1, 0, 0,    1, 0,    1, 1, 0, 0, 0, 0,    2, 0, 0,    2, 0, 0 };

   const double value0 = theDMRG->getSpecificCoefficient( coeff0 );
   const double value1 = theDMRG->getSpecificCoefficient( coeff1 );
   const double value2 = theDMRG->getSpecificCoefficient( coeff2 );
   cout << "Coeff of main contribution   ^9 Sigma_g^- = " << value0 << endl;
   cout << "Coeff of | pi_x > excitation ^9 Sigma_g^- = " << value1 << endl;
   cout << "Coeff of | pi_y > excitation ^9 Sigma_g^- = " << value2 << endl;
   const double value3 = theDMRG->getSpecificCoefficient( coeff3 );
   const double value4 = theDMRG->getSpecificCoefficient( coeff4 );
   const double value5 = theDMRG->getSpecificCoefficient( coeff5 );
   cout << "Coeff of main contribution   ^7 Delta_u   = " << value3 << endl;
   cout << "Coeff of | pi_x > excitation ^7 Delta_u   = " << value4 << endl;
   cout << "Coeff of | pi_y > excitation ^7 Delta_u   = " << value5 << endl;
   const double value6 = theDMRG->getSpecificCoefficient( coeff6 );
   const double value7 = theDMRG->getSpecificCoefficient( coeff7 );
   cout << "Coeff of main contrib  anion ^8 Sigma_u^- = " << value6 << endl;
   cout << "Coeff of main contrib cation ^8 Sigma_u^- = " << value7 << endl;

}

