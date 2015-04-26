#
#   CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
#   Copyright (C) 2013-2015 Sebastian Wouters
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

cdef extern from "chemps2/DMRGSCFoptions.h" namespace "CheMPS2":
    cdef cppclass DMRGSCFoptions:
        DMRGSCFoptions() except +
        bint getDoDIIS()
        double getDIISGradientBranch()
        int getNumDIISVecs()
        bint getStoreDIIS()
        int getMaxIterations()
        double getGradientThreshold()
        bint getStoreUnitary()
        int getWhichActiveSpace()
        bint getDumpCorrelations()
        bint getStateAveraging()
        void setDoDIIS(const bint)
        void setDIISGradientBranch(const double)
        void setNumDIISVecs(const int)
        void setStoreDIIS(const bint)
        void setMaxIterations(const int)
        void setGradientThreshold(const double)
        void setStoreUnitary(const bint)
        void setWhichActiveSpace(const int)
        void setDumpCorrelations(const bint)
        void setStateAveraging(const bint)

