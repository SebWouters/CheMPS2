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

import numpy as np
import sys
import PyCheMPS2
import ctypes

# Set the seed of the random number generator and cout.precision
Initializer = PyCheMPS2.PyInitialize()
Initializer.Init()

# Read in the FCIDUMP
psi4group = 7 # d2h: see chemps2/Irreps.h
filename  = b'../../tests/matrixelements/N2.CCPVDZ.FCIDUMP'
orbirreps = np.array([-1, -1], dtype=ctypes.c_int) # CheMPS2 reads it in from FCIDUMP
Ham = PyCheMPS2.PyHamiltonian( -1, psi4group, orbirreps, filename )

# Define the symmetry sector
TwoS = 0     # Two times the targeted spin
N = 14       # The number of electrons
Irrep = 0    # The targeted irrep

# Define the CASSCF
DOCC  = np.array([ 3, 0, 0, 0, 0, 2, 1, 1 ], dtype=ctypes.c_int) # see N2.ccpvdz.out
SOCC  = np.array([ 0, 0, 0, 0, 0, 0, 0, 0 ], dtype=ctypes.c_int)
NOCC  = np.array([ 1, 0, 0, 0, 0, 1, 0, 0 ], dtype=ctypes.c_int)
NDMRG = np.array([ 2, 0, 1, 1, 0, 2, 1, 1 ], dtype=ctypes.c_int)
NVIRT = np.array([ 4, 1, 2, 2, 1, 4, 2, 2 ], dtype=ctypes.c_int)
theDMRGSCF = PyCheMPS2.PyCASSCF(Ham, DOCC, SOCC, NOCC, NDMRG, NVIRT)

# Setting up the ConvergenceScheme
OptScheme = PyCheMPS2.PyConvergenceScheme( 1 ) # 1 instruction
OptScheme.set_instruction( 0, 1000, 1e-8, 20, 0.0, 1e-8 )

# Setting the DMRGSCFoptions and run DMRGSCF
root_num = 1 # Ground state only
theDMRGSCFoptions = PyCheMPS2.PyDMRGSCFoptions()
theDMRGSCFoptions.setDoDIIS(True)
theDMRGSCFoptions.setWhichActiveSpace(1) # 1 means natural orbitals
theDMRGSCFoptions.setDumpCorrelations(True)
Energy = theDMRGSCF.solve(N, TwoS, Irrep, OptScheme, root_num, theDMRGSCFoptions)

# Clean-up
if theDMRGSCFoptions.getStoreUnitary():
    theDMRGSCF.deleteStoredUnitary()
if theDMRGSCFoptions.getStoreDIIS():
    theDMRGSCF.deleteStoredDIIS()

# The order of deallocation matters!
del theDMRGSCFoptions
del OptScheme
del theDMRGSCF
del Ham
del Initializer

# Check whether the test succeeded
if (np.fabs(Energy + 109.103502335253) < 1e-8):
    print("================> Did test 8 succeed : yes")
else:
    print("================> Did test 8 succeed : no")

