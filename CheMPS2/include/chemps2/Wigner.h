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

#ifndef WIGNER_CHEMPS2_H
#define WIGNER_CHEMPS2_H

#define CHEMPS2_WIGNER_FACTORIAL_MAX 191
#define CHEMPS2_WIGNER_MAX_2J        95   // Maximum factorial = (4j+1)!   <=>   2j = 95

namespace CheMPS2{
/** Wigner class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date May 23, 2016

    The Wigner class allows to calculate Wigner-nj symbols.
*/
   class Wigner{

      public:

         //! Two times the maximum value of the angular momentum which is allowed
         static int max_2j();

         //! Wigner-3j symbol (gsl api)
         /** \param two_ja Two times the first spin
             \param two_jb Two times the second spin
             \param two_jc Two times the third spin
             \param two_ma Two times the first spin projection
             \param two_mb Two times the second spin projection
             \param two_mc Two times the thirs spin projection
             \return The corresponding wigner-3j value */
         static double wigner3j( const int two_ja, const int two_jb, const int two_jc, const int two_ma, const int two_mb, const int two_mc );

         //! Wigner-6j symbol (gsl api)
         /** \param two_ja Two times the first spin
             \param two_jb Two times the second spin
             \param two_jc Two times the third spin
             \param two_jd Two times the fourth spin
             \param two_je Two times the fifth spin
             \param two_jf Two times the sixth spin
             \return The corresponding wigner-6j value */
         static double wigner6j( const int two_ja, const int two_jb, const int two_jc, const int two_jd, const int two_je, const int two_jf );

         //! Wigner-9j symbol (gsl api)
         /** \param two_ja Two times the first spin
             \param two_jb Two times the second spin
             \param two_jc Two times the third spin
             \param two_jd Two times the fourth spin
             \param two_je Two times the fifth spin
             \param two_jf Two times the sixth spin
             \param two_jg Two times the seventh spin
             \param two_jh Two times the eigth spin
             \param two_ji Two times the nineth spin
             \return The corresponding wigner-9j value */
         static double wigner9j( const int two_ja, const int two_jb, const int two_jc, const int two_jd, const int two_je, const int two_jf, const int two_jg, const int two_jh, const int two_ji );

      private:

         // List of square roots of factorials
         static const long double sqrt_fact[ CHEMPS2_WIGNER_FACTORIAL_MAX + 1 ];

         // Test triangle conditions
         static bool triangle_fails( const int two_ja, const int two_jb, const int two_jc );

         // Delta function for the Wigner-6j terms
         static long double sqrt_delta( const int two_ja, const int two_jb, const int two_jc );

   };
}

#endif
