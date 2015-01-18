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
   O  0.000000000  0.00000000  0.000000000;
   H  0.790689766  0.00000000  0.612217330;
   H -0.790689766  0.00000000  0.612217330
           '''
mol.basis = '6-31g'
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
DSU2        = np.array([   100,  3000])
Econv       = np.array([ 1e-10, 1e-10])
MaxSweeps   = np.array([     3,    10])
NoisePrefac = np.array([   0.1,   0.0])

# In the following two lines, the orbitals are specified by their R(O)HF energy ordering
frozen = np.array([ ])                      # FCI all electron all orbital calculation
active = np.array(range( mf.mol.nao_nr() )) # FCI all electron all orbital calculation

##########################################
#   Do dmrgci calculation with CheMPS2   #
##########################################
E_dmrgci, TwoRDM = dmrgci( mf, TwoS, Nelec, Irrep, DSU2, Econv, MaxSweeps, NoisePrefac, frozen, active )

##############################
#   Example usage of 2-RDM   #
##############################
Norbs = mf.mol.nao_nr()
OEINT = np.dot(np.dot( mf.mo_coeff.T, mf.mol.intor('cint1e_kin_sph') + mf.mol.intor('cint1e_nuc_sph') ), mf.mo_coeff )
TEINT = ao2mo.outcore.full_iofree(mf.mol, mf.mo_coeff, compact=False).reshape(Norbs,Norbs,Norbs,Norbs)
TwoRDMenergy = mf.mol.energy_nuc()
TwoRDM = np.array( TwoRDM )
OneRDM = np.einsum('ijkk->ij', TwoRDM) / ( Nelec - 1.0 )
TwoRDMenergy += np.sum( np.multiply( OneRDM , OEINT ) )
TwoRDMenergy += 0.5 * np.sum( np.multiply( TwoRDM , TEINT ) )
print "Energy( 2-RDM * Ham ) =", TwoRDMenergy
print "Energy( dmrgci )      =", E_dmrgci

