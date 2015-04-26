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

import numpy as np
import sys
import PyCheMPS2
import ReadinHamiltonianFCIDUMP
import ctypes

# Set the seed of the random number generator and cout.precision
Initializer = PyCheMPS2.PyInitialize()
Initializer.Init()

# Read in the FCIDUMP
Ham = ReadinHamiltonianFCIDUMP.Read('../../tests/matrixelements/N2.CCPVDZ.FCIDUMP', 'd2h')
DOCC = np.array([ 3, 0, 0, 0, 0, 2, 1, 1 ], dtype=ctypes.c_int) # see N2.ccpvdz.out
SOCC = np.zeros([ 8 ], dtype=ctypes.c_int)
L = Ham.getL()

# Define the symmetry sector
TwoS = 0     # Two times the targeted spin
N = 14       # The number of electrons
Irrep = 0    # The targeted irrep

# Define the CASSCF
theDMRGSCF = PyCheMPS2.PyCASSCF(Ham, DOCC, SOCC)

# Define the active space
Nocc  = np.zeros([L], dtype=ctypes.c_int)
NDMRG = np.zeros([L], dtype=ctypes.c_int)
Nvirt = np.zeros([L], dtype=ctypes.c_int)
Nocc[0] = Nocc[5] = 1
Nocc[1] = Nocc[2] = Nocc[3] = Nocc[4] = Nocc[6] = Nocc[7] = 0
NDMRG[0] = NDMRG[5] = 4
NDMRG[1] = NDMRG[4] = 0
NDMRG[2] = NDMRG[3] = NDMRG[6] = NDMRG[7] = 1
Nvirt[0] = Nvirt[5] = 2
Nvirt[2] = Nvirt[3] = Nvirt[6] = Nvirt[7] = 2
Nvirt[1] = Nvirt[4] = 1
theDMRGSCF.setupStart(Nocc,NDMRG,Nvirt)

# Setting up the ConvergenceScheme
# setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
OptScheme = PyCheMPS2.PyConvergenceScheme(1) # 1 instruction
OptScheme.setInstruction(0, 1000, 1e-8, 20, 0.0)

# Setting the DMRGSCFoptions and run DMRGSCF
rootNum = 1 # Ground state only
theDMRGSCFoptions = PyCheMPS2.PyDMRGSCFoptions()
theDMRGSCFoptions.setDoDIIS(True)
theDMRGSCFoptions.setWhichActiveSpace(2) # 2 means localized orbitals
Energy = theDMRGSCF.doCASSCFnewtonraphson(N, TwoS, Irrep, OptScheme, rootNum, theDMRGSCFoptions)

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
if (np.fabs(Energy + 109.15104350802) < 1e-10):
    print "================> Did test 9 succeed : yes"
else:
    print "================> Did test 9 succeed : no"

