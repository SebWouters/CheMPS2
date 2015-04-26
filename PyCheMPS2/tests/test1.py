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
Ham = ReadinHamiltonianFCIDUMP.Read('../../tests/matrixelements/N2.STO3G.FCIDUMP', 'd2h')

# Define the symmetry sector
Nelec = 14 # The number of electrons
TwoS  = np.array([0, 2, 4, 4, 2, 2], dtype=ctypes.c_int) # Two times the targeted spin
Irrep = np.array([0, 5, 0, 5, 2, 6], dtype=ctypes.c_int) # The targeted irreps (Ag, B1u, Ag, ... )
FCIenergies  = np.zeros([len(TwoS)] ,dtype=ctypes.c_double)
DMRGenergies = np.ones( [len(TwoS)] ,dtype=ctypes.c_double)

for cnt in range(0, len(TwoS)):

    # Do FCI calculation
    Nel_up   = ( Nelec + TwoS[cnt] ) / 2
    Nel_down = ( Nelec - TwoS[cnt] ) / 2
    maxMemWorkMB = 10.0
    FCIverbose = 1
    theFCI = PyCheMPS2.PyFCI(Ham, Nel_up, Nel_down, Irrep[cnt], maxMemWorkMB, FCIverbose)
    GSvector = np.zeros([ theFCI.getVecLength() ], dtype=ctypes.c_double)
    GSvector[ theFCI.LowestEnergyDeterminant() ] = 1.0
    FCIenergies[cnt] = theFCI.GSDavidson(GSvector)
    theFCI.CalcSpinSquared(GSvector)
    del theFCI
    
    # Setting up the Problem
    Prob = PyCheMPS2.PyProblem(Ham, TwoS[cnt], Nelec, Irrep[cnt])
    # Prob.SetupReorderD2h() # Not used in cpp file
    
    # To perform DMRG, a set of convergence instructions should be added as well (normally more than 1 instruction should be used)
    OptScheme = PyCheMPS2.PyConvergenceScheme(1) # 1 instruction
    # setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
    OptScheme.setInstruction(0, 2000, 1e-10, 100, 0.0)
    
    # Do DMRG calculation
    theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
    DMRGenergies[cnt] = theDMRG.Solve()
    theDMRG.calc2DMandCorrelations()

    # Clean-up
    # theDMRG.deleteStoredMPS()
    theDMRG.deleteStoredOperators()
    del theDMRG
    del OptScheme
    del Prob
    
# Clean up
del Ham
del Initializer

# Check whether the test succeeded
success = True
for cnt in range(0, len(TwoS)):
    success = success and (np.fabs(FCIenergies[cnt] - DMRGenergies[cnt]) < 1e-10)
if (success):
    print "================> Did test 1 succeed : yes"
else:
    print "================> Did test 1 succeed : no"

