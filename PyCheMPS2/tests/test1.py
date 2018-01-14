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
orbirreps = np.array( [ -1, -1 ], dtype=ctypes.c_int ) # CheMPS2 reads it in from FCIDUMP
Ham = PyCheMPS2.PyHamiltonian( -1, psi4group, orbirreps, filename )

# Define the symmetry sector
Nelec  = 14 # The number of electrons
TwoS   = np.array( [ 0, 2, 4, 4, 2, 2 ], dtype=ctypes.c_int ) # Two times the targeted spin
Irrep  = np.array( [ 0, 5, 0, 5, 2, 6 ], dtype=ctypes.c_int ) # The targeted irreps (Ag, B1u, Ag, ... )
E_FCI  = np.zeros( [ len( TwoS ) ] ,dtype=ctypes.c_double )
E_DMRG = np.zeros( [ len( TwoS ) ] ,dtype=ctypes.c_double )
C_dev  = np.zeros( [ len( TwoS ) ] ,dtype=ctypes.c_double )
S_FCI  = np.zeros( [ len( TwoS ) ] ,dtype=ctypes.c_double )

# Obtain a list of single species Slater determinants and their irrep
def stringlist( local_ham, num_elec ):
    arrays = []
    irreps = []
    for counter in range( 2**local_ham.getL() ):
        thestring = np.array( list( bin( counter )[ 2 : ] ) )
        theints   = map( int, thestring )
        if ( np.sum( theints ) == num_elec ):
            thearray = np.zeros( [ local_ham.getL() ], dtype=ctypes.c_int )
            thearray[ local_ham.getL() - len( theints ): ] = theints
            theirrep = 0
            for orb in range( local_ham.getL() ):
                if ( thearray[ orb ] == 1 ):
                    theirrep = theirrep ^ local_ham.getOrbitalIrrep( orb )
            arrays.append( thearray )
            irreps.append( theirrep )
    return ( arrays, irreps )

def relative_phase( L, string_up, string_down ):

    phase = 1
    for orb_down in range( L - 1 ):
        if ( string_down[ orb_down ] == 1 ):
            sum_up = np.sum( string_up[ orb_down + 1 : ] )
            if (( sum_up % 2 ) == 1 ):
                phase *= -1
    return phase

for sector in range( len( TwoS ) ):

    # Setting up the Problem
    prob = PyCheMPS2.PyProblem( Ham, TwoS[ sector ], Nelec, Irrep[ sector ] )

    # To perform DMRG, a set of convergence instructions should be provided
    opt_scheme = PyCheMPS2.PyConvergenceScheme( 2 ) # 2 instructions
    # setInstruction( counter, virtual_dimension, energy_convergence, max_sweeps, noise_prefactor, dvdson_rtol )
    opt_scheme.set_instruction( 0,  500, 1e-10,  2, 0.0, 1e-5  )
    opt_scheme.set_instruction( 1, 1000, 1e-10, 30, 0.0, 1e-10 ) # Tight convergence for accurate FCI coefficients

    # Do DMRG calculation
    dmrg_solver = PyCheMPS2.PyDMRG( prob, opt_scheme )
    E_DMRG[ sector ] = dmrg_solver.Solve()
    dmrg_solver.calc2DMandCorrelations()

    # Do FCI calculation
    nelec_up   = ( Nelec + TwoS[ sector ] ) / 2
    nelec_down = ( Nelec - TwoS[ sector ] ) / 2
    workmem_mb = 10.0
    verbose    = 1
    fci_solver = PyCheMPS2.PyFCI( Ham, nelec_up, nelec_down, Irrep[ sector ], workmem_mb, verbose )
    GSvector   = np.zeros( [ fci_solver.getVecLength() ], dtype=ctypes.c_double )
    GSvector[ fci_solver.LowestEnergyDeterminant() ] = 1.0
    E_FCI[ sector ] = fci_solver.GSDavidson( GSvector )
    S_FCI[ sector ] = fci_solver.CalcSpinSquared( GSvector )

    # Compare the FCI and DMRG determinant coefficients
    list_up,   irrep_up   = stringlist( Ham, nelec_up   )
    list_down, irrep_down = stringlist( Ham, nelec_down )
    rms_error1 = 0.0
    rms_error2 = 0.0
    for count_up in range( len( irrep_up ) ):
        for count_down in range( len( irrep_down ) ):
            if (( irrep_up[ count_up ] ^ irrep_down[ count_down ] ) == Irrep[ sector ] ):
                dmrg_coeff = dmrg_solver.getFCIcoefficient( list_up[ count_up ], list_down[ count_down ] )
                fci_coeff  =  fci_solver.getFCIcoefficient( list_up[ count_up ], list_down[ count_down ], GSvector )
                phase_diff = relative_phase( Ham.getL(), list_up[ count_up ], list_down[ count_down ] )
                temp1      = dmrg_coeff - phase_diff * fci_coeff
                temp2      = dmrg_coeff + phase_diff * fci_coeff
                rms_error1 += temp1 * temp1
                rms_error2 += temp2 * temp2
    C_dev[ sector ] = np.sqrt( min( rms_error1, rms_error2 ) ) # The global phase of the wavefunction is arbitrary, hence.
    print("RMS difference FCI and DMRG determinant coefficients =", C_dev[ sector ])

    # Clean-up
    # dmrg_solver.deleteStoredMPS()
    dmrg_solver.deleteStoredOperators()
    del fci_solver
    del dmrg_solver
    del opt_scheme
    del prob

# Clean up
del Ham
del Initializer

# Check whether the test succeeded
success = True
for sector in range( len( TwoS ) ):
    success = success and ( np.fabs( E_FCI[ sector ] - E_DMRG[ sector ] ) < 1e-8 )
    success = success and ( C_dev[ sector ] < 1e-5 )
    success = success and ( np.fabs( S_FCI[ sector ] - 0.25 * TwoS[ sector ] * ( TwoS[ sector ] + 2 ) ) < 1e-8 )
if ( success ):
    print("================> Did test 1 succeed : yes")
else:
    print("================> Did test 1 succeed : no")

