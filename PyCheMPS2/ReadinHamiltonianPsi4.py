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

def Read(filename):
    thefile = open(filename, 'r')
    stop = False
    while not stop:
        line = thefile.readline()
        stop = 'CheMPS' in line
    line = thefile.readline()
    groupName = line.split()[3]
    groupNumber = -1
    if ('c1' in groupName) and (len(groupName)==2):
        groupNumber = 0
    if ('ci' in groupName) and (len(groupName)==2):
        groupNumber = 1
    if ('c2' in groupName) and (len(groupName)==2):
        groupNumber = 2
    if ('cs' in groupName) and (len(groupName)==2):
        groupNumber = 3
    if ('d2' in groupName) and (len(groupName)==2):
        groupNumber = 4
    if ('c2v' in groupName) and (len(groupName)==3):
        groupNumber = 5
    if ('c2h' in groupName) and (len(groupName)==3):
        groupNumber = 6
    if ('d2h' in groupName) and (len(groupName)==3):
        groupNumber = 7
    if (groupNumber==-1):
        print "ERROR: Group ", groupName, " not recognized."
    else:
        print "The Hamiltonian has symmetry group", groupName, "with Psi4-number", groupNumber
    line = thefile.readline() #Number of irreps: skip
    line = thefile.readline()
    Econst = float(line.split()[4])
    print "The nuclear repulsion energy =",Econst
    line = thefile.readline()
    L = int(line.split()[5])
    print "The number of MOs =",L
    line = thefile.readline() #Only text: skip
    line = thefile.readline()
    linesplit = line.split()
    orbirreps = np.zeros([L], dtype=ctypes.c_int)
    for cnt in range(0, L):
        orbirreps[cnt] = int(linesplit[cnt])
    print "The orbital irreps are", orbirreps
    Ham = PyCheMPS2.PyHamiltonian(L, groupNumber, orbirreps)
    line = thefile.readline() #DOCC: skip
    line = thefile.readline() #SOCC: skip
    line = thefile.readline() #Header OEI: skip
    Ham.setEconst(Econst)
    stop = False
    while not stop:
        line = thefile.readline()
        stop = 'TEI' in line
        if not stop:
            linesplit = line.split()
            row = int(linesplit[0])
            col = int(linesplit[1])
            value = float(linesplit[2])
            Ham.setTmat(row, col, value)
    stop = False
    while not stop:
        line = thefile.readline()
        stop = '***' in line
        if not stop:
            linesplit = line.split()
            i1 = int(linesplit[0])
            i2 = int(linesplit[1])
            i3 = int(linesplit[2])
            i4 = int(linesplit[3])
            value = float(linesplit[4])
            Ham.setVmat(i1, i3, i2, i4, value) #CheMPS2 uses physics notation; Psi4 chemistry notation
    thefile.close()
    return Ham
    
def DOCCandSOCC(filename):
    thefile = open(filename, 'r')
    stop = False
    while not stop:
        line = thefile.readline()
        stop = 'CheMPS' in line
    line = thefile.readline() #Point group: skip
    line = thefile.readline()
    numIrreps = int(line.split()[2])
    DOCC = np.zeros([numIrreps], dtype=ctypes.c_int)
    SOCC = np.zeros([numIrreps], dtype=ctypes.c_int)
    line = thefile.readline() #Nuclear repulsion: skip
    line = thefile.readline() #Number of orbitals: skip
    line = thefile.readline() #Text line: skip
    line = thefile.readline() #Orbital irreps: skip
    line = thefile.readline()
    linesplit = line.split()[2:]
    for cnt in range(0, len(linesplit)):
        DOCC[cnt] = int(linesplit[cnt])
    line = thefile.readline()
    linesplit = line.split()[2:]
    for cnt in range(0, len(linesplit)):
        SOCC[cnt] = int(linesplit[cnt])
    return (DOCC, SOCC)

