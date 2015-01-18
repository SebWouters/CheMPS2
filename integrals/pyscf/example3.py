import sys
sys.path.append('/home/seba') #Folder in which PySCF is installed
from pyscf import gto, scf, ao2mo
from dmrgci import dmrgci
import numpy as np

##############################################
#   Standard PySCF RHF or ROHF calculation   #
##############################################
mol = gto.Mole()
mol.atom = '''
  N       0.0000   0.0000   0.000000000
  N       0.0000   0.0000   1.120797413
           '''
mol.basis = 'cc-pvdz'
mol.symmetry = 1
mol.charge = 0
mol.spin = 0 #2*S; multiplicity-1
mol.build()
mf = scf.RHF(mol)
mf.verbose = 0
mf.scf()

#################################
#   Creating input for dmrgci   #
#################################
TwoS  = 0             # 2*S
Nelec = mol.nelectron # Total number of electrons
Irrep = 0             # Wavefunction irrep
DSU2        = np.array([   500,  1000])
Econv       = np.array([ 1e-10, 1e-10])
MaxSweeps   = np.array([     3,    10])
NoisePrefac = np.array([   0.1,   0.0])

# In the following two lines, the orbitals are specified by their R(O)HF energy ordering
frozen = np.array([ 0, 1 ])                                  # Keep 1s core orbitals of N2 frozen
active = np.array([ 2, 4, 10, 13, 7, 8, 3, 9, 16, 17, 5, 6]) # 4 Ag, 1 B2g, 1 B3g, 4 B1u, 1 B2u, 1 B3u in active space

##########################################
#   Do dmrgci calculation with CheMPS2   #
##########################################
E_dmrgci, TwoRDM = dmrgci( mf, TwoS, Nelec, Irrep, DSU2, Econv, MaxSweeps, NoisePrefac, frozen, active )
print "Energy( dmrgci ) =", E_dmrgci

