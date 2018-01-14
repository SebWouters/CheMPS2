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
cimport ConvScheme
cimport DMRGSCFopt
from libcpp cimport bool

cdef extern from "chemps2/FCI.h" namespace "CheMPS2":
    cdef cppclass FCI:
        FCI(Ham.Hamiltonian *, const unsigned int, const unsigned int, const int, const double, const int) except +
        unsigned int getVecLength(const int)
        unsigned int LowestEnergyDeterminant()
        void FillRandom(const unsigned int, double *)
        double GSDavidson(double *)
        double Fill2RDM(double *, double *)
        void Fill3RDM(double *, double *)
        void Fill4RDM(double *, double *)
        void Diag4RDM(double * vector, double * three_rdm, const unsigned int orbz, double * output)
        double CalcSpinSquared(double *)
        double getFCIcoeff(int *, int *, double *)
        void RetardedGF(const double, const double, const unsigned int, const unsigned int, const bool, const double, double *, Ham.Hamiltonian *, double *, double *)
        void RetardedGF_addition(const double, const double, const unsigned int, const unsigned int, const bool, const double, double *, Ham.Hamiltonian *, double *, double *, double *, double *, double *)
        void RetardedGF_removal(const double, const double, const unsigned int, const unsigned int, const bool, const double, double *, Ham.Hamiltonian *, double *, double *, double *, double *, double *)
        void GFmatrix_addition(const double, const double, const double, int *, const unsigned int, int *, const unsigned int, const bool, double *, Ham.Hamiltonian *, double *, double *)
        void GFmatrix_removal(const double, const double, const double, int *, const unsigned int, int *, const unsigned int, const bool, double *, Ham.Hamiltonian *, double *, double *)
        void DensityResponseGF(const double, const double, const unsigned int, const unsigned int, const double, double *, double *, double *)
        void DensityResponseGF_forward(const double, const double, const unsigned int, const unsigned int, const double, double *, double *, double *, double *, double *, double *)
        void DensityResponseGF_backward(const double, const double, const unsigned int, const unsigned int, const double, double *, double *, double *, double *, double *, double *)

