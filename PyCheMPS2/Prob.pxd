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

cimport Ham

cdef extern from "chemps2/Problem.h" namespace "CheMPS2":
    cdef cppclass Problem:
        Problem(const Ham.Hamiltonian *, const int, const int, const int) except +
        int gL()
        int gSy()
        int gIrrep(const int)
        int gTwoS()
        int gN()
        int gIrrep()
        double gEconst()
        double gMxElement(const int, const int, const int, const int)
        void setMxElement(const int, const int, const int, const int, const double)
        void SetupReorderD2h()

