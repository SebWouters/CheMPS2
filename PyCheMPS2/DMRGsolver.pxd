#
#   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
#   Copyright (C) 2013-2016 Sebastian Wouters
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

cimport ConvScheme
cimport Prob
cimport Corr
cimport TwoRDM
cimport ThreeRDM
from libcpp.string cimport string

cdef extern from "chemps2/DMRG.h" namespace "CheMPS2":
    cdef cppclass DMRG:
        DMRG(const Prob.Problem *, const ConvScheme.ConvergenceScheme *, const bint makechkpt, const string tmpfolder) except +
        double Solve()
        void calc2DMandCorrelations()
        void calc_rdms_and_correlations(const bint do_3rdm)
        void Diag4RDM(double * output, const int ham_orbz, const bint last_case)
        TwoRDM.TwoDM * get2DM()
        ThreeRDM.ThreeDM * get3DM()
        Corr.Correlations * getCorrelations()
        void deleteStoredMPS()
        void deleteStoredOperators()
        void activateExcitations(const int)
        void newExcitation(const double)
        double getFCIcoefficient(int *, int *)

