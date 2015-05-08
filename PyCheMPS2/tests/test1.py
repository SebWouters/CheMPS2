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
RMS_detcoeff = np.ones( [len(TwoS)] ,dtype=ctypes.c_double)

# Obtain a list of single species Slater determinants and their irrep
def stringlist(theHam, theNel):
    arrays = []
    irreps = []
    for counter in range(2**theHam.getL()):
        thestring = np.array(list(bin(counter)[2:]))
        theints   = map(int, thestring)
        if (np.sum(theints) == theNel):
            thearray = np.zeros([theHam.getL()], dtype=ctypes.c_int)
            thearray[theHam.getL() - len(theints):] = theints
            theirrep = 0
            for orb in range(theHam.getL()):
                if (thearray[ orb ] == 1):
                    theirrep = theirrep ^ theHam.getOrbitalIrrep(orb) # XOR
            arrays.append(thearray)
            irreps.append(theirrep)
    return (arrays, irreps)

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
    
    # Setting up the Problem
    Prob = PyCheMPS2.PyProblem(Ham, TwoS[cnt], Nelec, Irrep[cnt])
    Prob.SetupReorderD2h() # Determinant coefficient comparison OK both with option ON and OFF
    
    # To perform DMRG, a set of convergence instructions should be added as well (normally more than 1 instruction should be used)
    OptScheme = PyCheMPS2.PyConvergenceScheme(1) # 1 instruction
    # setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
    OptScheme.setInstruction(0, 2000, 1e-10, 100, 0.0)
    
    # Do DMRG calculation
    theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
    DMRGenergies[cnt] = theDMRG.Solve()
    theDMRG.calc2DMandCorrelations()
    
    # Compare the FCI and DMRG determinant coefficients
    list_alpha, irrep_alpha = stringlist(Ham, Nel_up  )
    list_beta,  irrep_beta  = stringlist(Ham, Nel_down)
    rms_abs = 0.0
    for count_alpha in range(len(irrep_alpha)):
        for count_beta in range(len(irrep_beta)):
            if ((irrep_alpha[count_alpha] ^ irrep_beta[count_beta]) == Irrep[cnt]):
                dmrg_coeff = theDMRG.getFCIcoefficient(list_alpha[count_alpha], list_beta[count_beta])
                fci_coeff  =  theFCI.getFCIcoefficient(list_alpha[count_alpha], list_beta[count_beta], GSvector)
                temp       = abs(dmrg_coeff) - abs(fci_coeff)
                rms_abs   += temp * temp
    RMS_detcoeff[cnt] = np.sqrt(rms_abs)
    print "RMS difference FCI and DMRG determinant coefficients =", RMS_detcoeff[cnt]

    # Clean-up
    # theDMRG.deleteStoredMPS()
    theDMRG.deleteStoredOperators()
    del theFCI
    del theDMRG
    del OptScheme
    del Prob
    
# Clean up
del Ham
del Initializer

# Check whether the test succeeded
success = True
for cnt in range(len(TwoS)):
    success = success and (np.fabs(FCIenergies[cnt] - DMRGenergies[cnt]) < 1e-10)
    success = success and (RMS_detcoeff[cnt] < 1e-6) # Energy converges quadratically in wfn error, cfr. EPJD 68 (9), 272 (2014)
if (success):
    print "================> Did test 1 succeed : yes"
else:
    print "================> Did test 1 succeed : no"

