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

########################################
### Square 2D Hubbard model with PBC ###
########################################

L_linear = 3                    # Linear size
L_square = L_linear * L_linear  # Number of orbitals
group = 0                       # C1 symmetry
U =  5.0                        # On-site repulsion
T = -1.0                        # Hopping term

Nelec = 9                       # Number of electrons
TwoS = 1                        # Two times the spin
Irrep = 0                       # Irrep = A (C1 symmetry)

# The Hamiltonian initializes all its matrix elements to 0.0
orbirreps = np.zeros([L_square], dtype=ctypes.c_int)
Ham = PyCheMPS2.PyHamiltonian(L_square, group, orbirreps)

# Fill with the site-basis matrix elements
for orb in range(L_square):
    Ham.setVmat(orb,orb,orb,orb,U)
for ix in range(L_linear):
    for iy in range(L_linear):
        idx1 = ix + L_linear * iy                        # This site
        idx2 = (( ix + 1 ) % L_linear) + L_linear * iy   # Right neighbour (PBC)
        idx3 = ix + L_linear * ((( iy + 1 ) % L_linear)) # Upper neighbour (PBC)
        Ham.setTmat(idx1,idx2,T)
        Ham.setTmat(idx1,idx3,T)

# Setting up the Problem
Prob = PyCheMPS2.PyProblem(Ham, TwoS, Nelec, Irrep)

# Setting up the ConvergenceScheme
# setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
OptScheme = PyCheMPS2.PyConvergenceScheme(2) # 2 instructions
OptScheme.setInstruction(0,  500, 1e-10,  3, 0.05)
OptScheme.setInstruction(1, 1000, 1e-10, 10, 0.0 )

# Run ground state calculation
theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
EnergySite = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()

# Clean-up DMRG
# theDMRG.deleteStoredMPS()
theDMRG.deleteStoredOperators()
del theDMRG

#################################################################################################################
### Hack: overwrite the matrix elements in momentum space (4-fold symmetry!!!) directly in the Problem object ###
#################################################################################################################
theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme) # Prob->construct_mxelem() is called in DMRG constructor
for orb1 in range(L_square):
    k1x = orb1 % L_linear
    k1y = orb1 / L_linear
    Telem1 = 2*T*(np.cos((2*np.pi*k1x)/L_linear) + np.cos((2*np.pi*k1y)/L_linear))
    for orb2 in range(L_square):
        k2x = orb2 % L_linear
        k2y = orb2 / L_linear
        Telem2 = 2*T*(np.cos((2*np.pi*k2x)/L_linear) + np.cos((2*np.pi*k2y)/L_linear))
        for orb3 in range(L_square):
            k3x = orb3 % L_linear
            k3y = orb3 / L_linear
            for orb4 in range(L_square):
                k4x = orb4 % L_linear
                k4y = orb4 / L_linear
                kx_conservation = False
                if (((k1x+k2x) % L_linear) == ((k3x+k4x) % L_linear)):
                    kx_conservation = True
                ky_conservation = False
                if (((k1y+k2y) % L_linear) == ((k3y+k4y) % L_linear)):
                    ky_conservation = True
                temp = 0.0
                if ( kx_conservation and ky_conservation ):
                    temp += U/L_square
                if (( orb1 == orb3 ) and ( orb2 == orb4 )):
                    temp += (Telem1+Telem2)/(Nelec-1)
                Prob.setMxElement(orb1,orb2,orb3,orb4,temp)
theDMRG.PreSolve() # New matrix elements require reconstruction of complementary renormalized operators
EnergyMomentum = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()

# Clean-up
# theDMRG.deleteStoredMPS()
theDMRG.deleteStoredOperators()
del theDMRG
del OptScheme
del Prob
del Ham
del Initializer

# Check whether the test succeeded
if (np.fabs(EnergySite - EnergyMomentum) < 1e-8):
    print("================> Did test 9 succeed : yes")
else:
    print("================> Did test 9 succeed : no")

