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

#ifndef EXCITATION_CHEMPS2_H
#define EXCITATION_CHEMPS2_H

#include "TensorL.h"
#include "TensorO.h"
#include "SyBookkeeper.h"
#include "Sobject.h"
#include "Options.h"

namespace CheMPS2{
/** Excitation class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date April 8, 2016

    The Excitation class contains multipliction routines to calculate ( alpha * E_tl + beta * E_lt + gamma ) x vector, where the orbital irreps of l and t are the same. */
   class Excitation{

      public:

         //! Matrix-vector multiplication S_up = ( alpha * E_{orb1,orb2} + beta * E_{orb2,orb1} + gamma ) x S_down
         /** \param book_up   The upper symmetry bookkeeper (result vector)
             \param book_down The lower symmetry bookkeeper (origin vector)
             \param orb1      The first  orbital with orb1 < orb2
             \param orb2      The second orbital with orb1 < orb2
             \param alpha     Prefactor of E_{orb1,orb2}
             \param beta      Prefactor of E_{orb2,orb1}
             \param gamma     Prefactor of operator 1
             \param S_up      Where the result should be stored
             \param S_down    The original S-object
             \param overlaps  Overlap tensors
             \param regular   Regular    renormalized single quantized operator
             \param trans     Transposed renormalized single quantized operator */
         static double matvec( const SyBookkeeper * book_up, const SyBookkeeper * book_down, const int orb1, const int orb2, const double alpha, const double beta, const double gamma, Sobject * S_up, Sobject * S_down, TensorO ** overlaps, TensorL ** regular, TensorL ** trans );

      private:

         static void clear( const int ikappa, Sobject * S_up );
         static double neighbours( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, const double beta, const double gamma, Sobject * S_up, Sobject * S_down );
         static void first_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Rtrans );
         static void second_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Rregular );
         static double third_left( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Rovlp, double * workmem );
         static void first_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Ltrans );
         static void second_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Lregular );
         static double third_right( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Lovlp, double * workmem );
         static void first_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double alpha, Sobject * S_up, Sobject * S_down, TensorL * Ltrans, TensorL * Rtrans, double * workmem );
         static void second_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double beta, Sobject * S_up, Sobject * S_down, TensorL * Lregular, TensorL * Rregular, double * workmem );
         static double third_middle( const int ikappa, const SyBookkeeper * book_up, const SyBookkeeper * book_down, const double gamma, Sobject * S_up, Sobject * S_down, TensorO * Lovlp, TensorO * Rovlp, double * workmem1, double * workmem2 );

   };
}

#endif
