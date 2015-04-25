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
import ctypes

def Read(filename, groupName):

    groupNumber = -1
    if ('c1' in groupName) and (len(groupName)==2):
        groupNumber = 0
        irrep_molpro = np.array([ 1 ], dtype=ctypes.c_int)
    if ('ci' in groupName) and (len(groupName)==2):
        groupNumber = 1
        irrep_molpro = np.array([ 1, 2 ], dtype=ctypes.c_int)
    if ('c2' in groupName) and (len(groupName)==2):
        groupNumber = 2
        irrep_molpro = np.array([ 1, 2 ], dtype=ctypes.c_int)
    if ('cs' in groupName) and (len(groupName)==2):
        groupNumber = 3
        irrep_molpro = np.array([ 1, 2 ], dtype=ctypes.c_int)
    if ('d2' in groupName) and (len(groupName)==2):
        groupNumber = 4
        irrep_molpro = np.array([ 1, 4, 3, 2 ], dtype=ctypes.c_int)
    if ('c2v' in groupName) and (len(groupName)==3):
        groupNumber = 5
        irrep_molpro = np.array([ 1, 4, 2, 3 ], dtype=ctypes.c_int)
    if ('c2h' in groupName) and (len(groupName)==3):
        groupNumber = 6
        irrep_molpro = np.array([ 1, 4, 2, 3 ], dtype=ctypes.c_int)
    if ('d2h' in groupName) and (len(groupName)==3):
        groupNumber = 7
        irrep_molpro = np.array([ 1, 4, 6, 7, 8, 5, 3, 2 ], dtype=ctypes.c_int)
    if (groupNumber==-1):
        print "ERROR: Group ", groupName, " not recognized."
    else:
        print "The Hamiltonian has symmetry group", groupName, "with psi4-number", groupNumber

    thefile = open(filename, 'r')
    line = thefile.readline() #Number of orbitals
    pos1 = line.find('NORB')
    pos1 = line.find('=', pos1)
    pos2 = line.find(',', pos1)
    part = line[pos1+1:pos2]
    L = int(part)
    
    line = thefile.readline() #Orbital irreps
    orbirreps = np.zeros([L], dtype=ctypes.c_int)
    pos1 = line.find('=')
    for orb in range(L):
        pos2 = line.find(',', pos1+1)
        part = line[pos1+1:pos2]
        molproirrep = int(part)
        orbirreps[orb] = -1
        for cnt in range( len( irrep_molpro ) ):
            if ( molproirrep == irrep_molpro[ cnt ] ):
                orbirreps[ orb ] = cnt
        pos1 = pos2
    print "orb2irrep =", orbirreps
    Ham = PyCheMPS2.PyHamiltonian(L, groupNumber, orbirreps)
    Ham.setEconst(0.0)
    for orb1 in range(L):
        for orb2 in range(L):
            if ( orbirreps[orb1] == orbirreps[orb2] ):
                Ham.setTmat(orb1, orb2, 0.0)
            for orb3 in range(L):
                for orb4 in range(L):
                    if ( orbirreps[orb1] ^ orbirreps[orb2] == orbirreps[orb3] ^ orbirreps[orb4] ):
                        Ham.setVmat(orb1, orb2, orb3, orb4, 0.0)
    
    line = thefile.readline() #Skip
    line = thefile.readline() #Skip
    
    stop = False
    while (not stop):
        line  = thefile.readline()
        parts = line.split()
        value = float(parts[0])
        ind1  = int(parts[1])
        ind2  = int(parts[2])
        ind3  = int(parts[3])
        ind4  = int(parts[4])
        
        if ( ind4 != 0 ):
            #CheMPS2 uses physics notation; FCIDUMP chemistry notation
            Ham.setVmat( ind1-1, ind3-1, ind2-1, ind4-1, value )
        else:
            if ( ind2 != 0 ):
                Ham.setTmat( ind1-1, ind2-1, value )
            else:
                Ham.setEconst( value )
                stop = True
                
    thefile.close()
    return Ham

