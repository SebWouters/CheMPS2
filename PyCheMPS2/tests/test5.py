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
filename  = b'../../tests/matrixelements/N2.STO3G.FCIDUMP'
orbirreps = np.array([-1, -1], dtype=ctypes.c_int) # CheMPS2 reads it in from FCIDUMP
Ham = PyCheMPS2.PyHamiltonian( -1, psi4group, orbirreps, filename )

# Define the symmetry sector
TwoS = 0     # Two times the targeted spin
N = 14       # The number of electrons
Irrep = 0    # The targeted irrep

# Setting up the Problem
Prob = PyCheMPS2.PyProblem(Ham, TwoS, N, Irrep)
Prob.SetupReorderD2h()

# Setting up the ConvergenceScheme
# setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
OptScheme = PyCheMPS2.PyConvergenceScheme(1) # 1 instruction
OptScheme.setInstruction(0, 1000, 1e-12, 100, 0.0)

# Do DMRG calculation for the ground state
theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
Energy0 = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()

# Calculate the lowest two excited states
theDMRG.activateExcitations(10)
theDMRG.newExcitation(20.0)
Energy1 = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()
theDMRG.newExcitation(20.0)
Energy2 = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()

# Clean-up
# theDMRG.deleteStoredMPS()
theDMRG.deleteStoredOperators()

# The order of deallocation matters!
del theDMRG
del OptScheme
del Prob
del Ham
del Initializer

# Check whether the test succeeded
OK0 = (np.fabs(Energy0 + 107.648250974014) < 1e-8)
OK1 = (np.fabs(Energy1 + 106.944757308768) < 1e-8)
OK2 = (np.fabs(Energy2 + 106.92314213886 ) < 1e-8)
if (OK0 and OK1 and OK2):
    print("================> Did test 5 succeed : yes")
else:
    print("================> Did test 5 succeed : no")

