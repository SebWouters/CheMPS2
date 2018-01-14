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

#######################
### BCS Hamiltonian ###
#######################

eps = np.array([ -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5 ], dtype=ctypes.c_double)
L = len( eps )
g = -1.0
power = 0.0

Nelec = L # Number of fermions in the model = Number of single-particle states
TwoS = 0  # Twice the total spin
Irrep = 0 # No point group is used, Irrep should ALWAYS be zero.

'''
   Model: h_ij = delta_ij eps[i]
          v_ijkl = delta_ij delta_kl g ( eps[i] * eps[k] ) ^ {power}
          h_ijkl = v_ijkl + ( delta_ik h_jl + delta_jl h_ik ) / ( N - 1 )
          Ham = 0.5 sum_ijkl h_ijkl sum_sigma,tau a^+_{i,sigma} a^+_{j,tau} a_{l,tau} a_{k,sigma}
'''

# The Hamiltonian initializes all its matrix elements to 0.0
orbirreps = np.zeros( [ L ], dtype=ctypes.c_int )
group = 0
Ham = PyCheMPS2.PyHamiltonian( L, group, orbirreps )

# Setting up the Problem
Prob = PyCheMPS2.PyProblem( Ham, TwoS, Nelec, Irrep )

# Setting up the ConvergenceScheme
# setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
OptScheme = PyCheMPS2.PyConvergenceScheme( 2 )
OptScheme.setInstruction( 0,  100, 1e-10, 10, 0.5 )
OptScheme.setInstruction( 1, 1000, 1e-10, 10, 0.0 )

# Run ground state calculation
theDMRG = PyCheMPS2.PyDMRG( Prob, OptScheme )
###############################################################################################
### Hack: overwrite the matrix elements with 4-fold symmetry directly in the Problem object ###
###############################################################################################
for orb1 in range( L ):
   for orb2 in range( L ):
      eri = g * ( abs( eps[ orb1 ] * eps[ orb2 ] )**power )
      oei = ( eps[ orb1 ] + eps[ orb2 ] ) / ( Nelec - 1 )
      if ( orb1 == orb2 ):
         Prob.setMxElement( orb1, orb1, orb2, orb2, eri + oei )
      else:
         Prob.setMxElement( orb1, orb1, orb2, orb2, eri )
         Prob.setMxElement( orb1, orb2, orb1, orb2, oei )
theDMRG.PreSolve() # New matrix elements require reconstruction of complementary renormalized operators
Energy = theDMRG.Solve()
theDMRG.calc2DMandCorrelations()
theDMRG.printCorrelations()

# Clean-up
# theDMRG.deleteStoredMPS()
theDMRG.deleteStoredOperators()
del theDMRG
del OptScheme
del Prob
del Ham
del Initializer

# Check whether the test succeeded
if ( np.fabs( Energy + 25.5134137600604 ) < 1e-8 ):
    print("================> Did test 12 succeed : yes")
else:
    print("================> Did test 12 succeed : no")

