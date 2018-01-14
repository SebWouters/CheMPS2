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

#ifndef SPECIAL_CHEMPS2_H
#define SPECIAL_CHEMPS2_H

namespace CheMPS2{
/** Special class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date April 7, 2016

    The special class contains static special functions used at various places throughout the code. */
   class Special{

      public:

         //! Phase function
         /** \param power_times_two Two times the power of the phase (-1)^{ power }
             \return The phase (-1)^{ power_times_two / 2 } */
         static int phase( const int power_times_two ){ return (((( power_times_two / 2 ) % 2 ) != 0 ) ? -1 : 1 ); }

         //! Triangle function for two variables
         /** \param global The triangle counter global = i + ( j * ( j + 1 ) ) / 2 with i <= j
             \param result Integer array of size 2 where to store i <= j */
         static void invert_triangle_two( const int global, int * result ){

            int j = 0;
            while ((( j + 1 ) * ( j + 2 )) <= ( 2 * global )){ j++; }
            result[ 1 ] = j;
            result[ 0 ] = global - ( j * ( j + 1 ) ) / 2;

         }

         //! LOWER triangle function for two variables
         /** \param global The lower triangle counter global = i + ( j * ( j - 1 ) ) / 2 with i < j
             \param result Integer array of size 2 where to store i < j */
         static void invert_lower_triangle_two( const int global, int * result ){

            int j = 1;
            while (( j * ( j + 1 )) <= ( 2 * global )){ j++; }
            result[ 1 ] = j;
            result[ 0 ] = global - ( j * ( j - 1 ) ) / 2;

         }

         //! Triangle function for three variables
         /** \param global The triangle counter global = i + ( j * ( j + 1 ) ) / 2 + ( k * ( k + 1 ) * ( k + 2 ) ) / 6 with i <= j <= k
             \param result Integer array of size 3 where to store i <= j <= k */
         static void invert_triangle_three( const int global, int * result ){

            int k = 0;
            while ((( k + 1 ) * ( k + 2 ) * ( k + 3 )) <= ( 6 * global )){ k++; }
            const int remainder = global - ( k * ( k + 1 ) * ( k + 2 ) ) / 6;
            int j = 0;
            while ((( j + 1 ) * ( j + 2 )) <= ( 2 * remainder )){ j++; }
            result[ 2 ] = k;
            result[ 1 ] = j;
            result[ 0 ] = remainder - ( j * ( j + 1 ) ) / 2;

         }

   };
}

#endif
