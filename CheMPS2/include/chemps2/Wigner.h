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

#ifndef WIGNER_CHEMPS2_H
#define WIGNER_CHEMPS2_H

#define CHEMPS2_WIGNER_MAX_FACTORIAL 312  // prime number no. 64
#define CHEMPS2_WIGNER_NUM_PRIMENUMS 64   // prime numbers 2 to 311 (313 = next one)
#define CHEMPS2_WIGNER_MAX_2J        144  // maximum factorial 6j symbols = (4j+1)! ==> max j = 77

namespace CheMPS2{

   struct prime_powers{
      short pow[ CHEMPS2_WIGNER_NUM_PRIMENUMS ];
   };

/** Wigner class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date May 23, 2016

    The Wigner class allows to calculate Wigner-nj symbols based on prime factorization of the factorials [WIGNER1].

    \section biblio_wigner References

    [WIGNER1]  Johansson and Forssen, Siam J. Sci. Comput. 38 (1), A376-A384 (2016). http://dx.doi.org/10.1137/15M1021908 \n
*/
   class Wigner{

      public:

         //! Get the maximum value of 2J
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

         // List of prime numbers
         static const short prime_numbers[ CHEMPS2_WIGNER_NUM_PRIMENUMS ];

         // List of factorial prime factorizations
         static const prime_powers factorials[ CHEMPS2_WIGNER_MAX_FACTORIAL + 1 ];

         // Test triangle conditions
         static bool triangle_fails( const int two_ja, const int two_jb, const int two_jc );

         // Delta function for the Wigner-6j terms
         static void delta( const int two_ja, const int two_jb, const int two_jc, short * pow );

         // Additional factor under the square root for the Wigner-3j symbols
         static void delta3j( const int two_ja, const int two_jb, const int two_jc, const int two_ma, const int two_mb, const int two_mc, short * pow );

         // Compute one of the terms in the k-summation for the Wigner-3j term (without phase)
         static void term3j( const int alpha1, const int alpha2, const int beta1, const int beta2, const int beta3, const int k, short * pow );

         // Compute one of the terms in the k-summation for the Wigner-6j term (without phase)
         static void term6j( const int alpha1, const int alpha2, const int alpha3, const int alpha4, const int beta1, const int beta2, const int beta3, const int k, short * pow );

         // Extract powers from pow_delta (of which a sqrt should be taken) so that all its entries are 0 or 1 and adjust the terms accordingly
         static void shift_sqrt( short * target, prime_powers * terms, const int num_terms );

         // Extract the denominator and adjust the terms accordingly
         static void extract_denominator( short * denom, prime_powers * terms, const int num_terms );

         // Convert the prime number factorization into a double
         static long double form_number_simple( short * pow );

   };
}

#endif
