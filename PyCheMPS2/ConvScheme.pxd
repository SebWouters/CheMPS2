#
#   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
#   Copyright (C) 2013-2018 Sebastian Wouters
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

cdef extern from "chemps2/ConvergenceScheme.h" namespace "CheMPS2":
    cdef cppclass ConvergenceScheme:
        ConvergenceScheme(const int) except +
        int getNInstructions()
        void setInstruction(const int, const int, const double, const int, const double)
        void set_instruction(const int, const int, const double, const int, const double, const double)
        int get_D(const int)
        double get_energy_conv(const int)
        int get_max_sweeps(const int)
        double get_noise_prefactor(const int)
        double get_dvdson_rtol(const int)


