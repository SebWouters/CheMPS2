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
import math as m
import sys
import PyCheMPS2
import ctypes
import os

# Set the seed of the random number generator and cout.precision
Initializer = PyCheMPS2.PyInitialize()
Initializer.Init()

# Read in the FCIDUMP
psi4group = 7 # d2h: see chemps2/Irreps.h
filename  = b'../../tests/matrixelements/O2.CCPVDZ.FCIDUMP'
orbirreps = np.array([-1, -1], dtype=ctypes.c_int) # CheMPS2 reads it in from FCIDUMP
Ham = PyCheMPS2.PyHamiltonian( -1, psi4group, orbirreps, filename )
L = Ham.getL()

# Dump the Hamiltonian HDF5 file
Ham.save()

# Create a second Hamiltonian to read it back in
groupNum = Ham.getNGroup()
orbirreps = np.zeros([L], dtype=ctypes.c_int)
for orb in range(0, L):
    orbirreps[orb] = Ham.getOrbitalIrrep(orb)
HamLoad = PyCheMPS2.PyHamiltonian(L, groupNum, orbirreps)

# Read it back in
HamLoad.read()

# Compare both
RMSabs = 0.0
temp = Ham.getEconst()
RMSabs += temp * temp
RMS = 0.0
temp = Ham.getEconst() - HamLoad.getEconst()
RMS += temp * temp
for i1 in range(0,L):
    for i2 in range(0,L):
        temp = Ham.getTmat(i1,i2)
        RMSabs += temp * temp
        temp = Ham.getTmat(i1,i2) - HamLoad.getTmat(i1,i2)
        RMS += temp * temp
        for i3 in range(0,L):
            for i4 in range(0,L):
                temp = Ham.getVmat(i1,i2,i3,i4)
                RMSabs += temp * temp
                temp = Ham.getVmat(i1,i2,i3,i4) - HamLoad.getVmat(i1,i2,i3,i4)
                RMS += temp * temp
RMS = m.sqrt(RMS)
RMSabs = m.sqrt(RMSabs)
print("The 2-norm of all Hamiltonian matrix elements is", RMSabs)
print("The RMS difference between Ham and HamLoad is", RMS)

# Clean-up
os.system('ls -alh CheMPS2_Ham*.h5')
os.system('rm CheMPS2_Ham*.h5')
os.system('ls -alh CheMPS2_Ham*.h5')
del HamLoad
del Ham
del Initializer

# Check whether the test succeeded
if (RMS < 1e-10):
    print("================> Did test 7 succeed : yes")
else:
    print("================> Did test 7 succeed : no")

