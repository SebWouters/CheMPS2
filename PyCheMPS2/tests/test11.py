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

########################
### 1D Hubbard model ###
########################

L = 10       # Number of lattice sites
Group = 0    # C1 symmetry
U = 2.0      # On-site repulsion
T = -1.0     # Hopping term

TwoS  = 5    # Two times the targeted spin
Nelec = 9    # The number of electrons
Irrep = 0    # The targeted irrep

# The Hamiltonian initializes all its matrix elements to 0.0
orbirreps = np.zeros([L], dtype=ctypes.c_int)
Ham = PyCheMPS2.PyHamiltonian(L, Group, orbirreps)
for cnt in range(0, L):
    Ham.setVmat(cnt, cnt, cnt, cnt, U)
for cnt in range(0, L-1):
    Ham.setTmat(cnt, cnt+1, T)

# Setting up the Problem
Prob = PyCheMPS2.PyProblem(Ham, TwoS, Nelec, Irrep)

# Setting up the ConvergenceScheme
# setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
OptScheme = PyCheMPS2.PyConvergenceScheme(2) # 2 instructions
OptScheme.setInstruction(0,   30, 1e-10,  3, 0.1)
OptScheme.setInstruction(1, 1000, 1e-10, 10, 0.0)

# Do DMRG calculation and print the correlations
theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
EnergyDMRG = theDMRG.Solve()
do_3rdm = True
theDMRG.calc_rdms_and_correlations( do_3rdm )
theDMRG.printCorrelations()

# Get a diagonal part of the 4-RDM from DMRG
ham_orbz = 3
dmrg_diag_4rdm = np.zeros([ L**6 ], dtype=ctypes.c_double)
theDMRG.Symm4RDM( dmrg_diag_4rdm, ham_orbz, ham_orbz, True )

# Do FCI calculation
Nel_up   = ( Nelec + TwoS ) / 2
Nel_down = ( Nelec - TwoS ) / 2
maxMemWorkMB = 10.0
FCIverbose = 1
theFCI = PyCheMPS2.PyFCI(Ham, Nel_up, Nel_down, Irrep, maxMemWorkMB, FCIverbose)
GSvector = np.zeros([ theFCI.getVecLength() ], dtype=ctypes.c_double)
theFCI.FillRandom( theFCI.getVecLength() , GSvector )
EnergyFCI = theFCI.GSDavidson(GSvector)
theFCI.CalcSpinSquared(GSvector)
L = Ham.getL()
TwoRDM = np.zeros([ L**4 ], dtype=ctypes.c_double)
theFCI.Fill2RDM(GSvector, TwoRDM)
RMSerror2DM = 0.0
for i in range(L):
    for j in range(L):
        for k in range(L):
            for l in range(L):
                temp = TwoRDM[i + L*(j + L*(k + L*l))] - theDMRG.get2DMA(i,j,k,l)
                RMSerror2DM += temp*temp
ThreeRDM = np.zeros([ L**6 ], dtype=ctypes.c_double)
theFCI.Fill3RDM(GSvector, ThreeRDM)
RMSerror3DM = 0.0
for i in range(L):
    for j in range(L):
        for k in range(L):
            for l in range(L):
                for m in range(L):
                    for n in range(L):
                        temp = ThreeRDM[i + L*(j + L*(k + L*(l + L*(m + L*n))))] - theDMRG.get3DM(i,j,k,l,m,n)
                        RMSerror3DM += temp*temp
fci_diag_4rdm = np.zeros([ L**6 ], dtype=ctypes.c_double)
theFCI.Diag4RDM( GSvector, ThreeRDM, ham_orbz, fci_diag_4rdm )
RMSerror4DM = np.linalg.norm( 0.5 * dmrg_diag_4rdm - fci_diag_4rdm )
RMSerror2DM = np.sqrt(RMSerror2DM)
RMSerror3DM = np.sqrt(RMSerror3DM)
print("Frobenius norm of the difference of the DMRG and FCI 2-RDM =", RMSerror2DM)
print("Frobenius norm of the difference of the DMRG and FCI 3-RDM =", RMSerror3DM)
print("Frobenius norm of the difference of the DMRG and FCI diag(4-RDM) for fixed orbital", ham_orbz, "=", RMSerror4DM)

# Clean-up
# theDMRG.deleteStoredMPS()
theDMRG.deleteStoredOperators()
del theFCI
del theDMRG
del OptScheme
del Prob
del Ham
del Initializer

# Check whether the test succeeded
if ((np.fabs(EnergyDMRG - EnergyFCI) < 1e-8) and (RMSerror2DM < 1e-3) and (RMSerror3DM < 1e-3) and (RMSerror4DM < 1e-3)):
    print("================> Did test 11 succeed : yes")
else:
    print("================> Did test 11 succeed : no")

