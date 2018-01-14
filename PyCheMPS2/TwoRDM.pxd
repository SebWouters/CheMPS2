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

cdef extern from "chemps2/TwoDM.h" namespace "CheMPS2":
    cdef cppclass TwoDM:
        double getTwoDMA_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4)
        double getTwoDMB_HAM(const int cnt1, const int cnt2, const int cnt3, const int cnt4)
        double trace()
        double energy()
        void save()
        void read()

