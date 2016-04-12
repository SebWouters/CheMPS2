import sys
sys.path.append('/usr/lib/python2.7/site-packages')
import PyCheMPS2
from pyscf import gto, scf, ao2mo, symm, fci
import numpy as np
import ctypes
import time

##################
#   Molecule 1   #
##################
mol1 = gto.Mole()
mol1.atom = '''
   N  0.0000   0.0000    0.000000
   N  0.0000   0.0000    1.090105
           '''
mol1.basis = 'sto-3g'
mol1.symmetry = 1
mol1.charge = 0
mol1.spin = 0 #2*S; multiplicity-1
mol1.build()

##################
#   Molecule 2   #
##################
mol2 = gto.Mole()
mol2.atom = '''
   O  0.000000000  0.00000000  0.000000000
   H  0.790689766  0.00000000  0.612217330
   H -0.790689766  0.00000000  0.612217330
           '''
mol2.basis = '6-31g'
mol2.symmetry = 1
mol2.charge = 0
mol2.spin = 0 #2*S; multiplicity-1
mol2.build()

for geval in range(1):

    if ( geval == 0 ):
        mol = mol1
        psi4group = 7 #d2h
        print "####################"
        print "#   N2 in sto-3g   #"
        print "####################"
    if ( geval == 1 ):
        mol = mol2
        psi4group = 5 #c2v
        print "####################"
        print "#   H2O in 6-31g   #"
        print "####################"
    mf = scf.RHF(mol)
    mf.verbose = 3
    mf.scf()

    #############################
    #   Build the Hamiltonian   #
    #############################
    L         = mf.mol.nao_nr()
    N         = mf.mol.nelectron
    CONST     = mf.mol.energy_nuc()
    OEImo     = np.dot( np.dot( mf.mo_coeff.T, mf.mol.intor('cint1e_kin_sph') + mf.mol.intor('cint1e_nuc_sph') ), mf.mo_coeff )
    TEImo     = ao2mo.outcore.full_iofree(mf.mol, mf.mo_coeff, compact=False).reshape(L,L,L,L)
    orb2irrep = np.array(symm.label_orb_symm(mf.mol, mf.mol.irrep_id, mf.mol.symm_orb, mf.mo_coeff), dtype=ctypes.c_int)
    for cnt in range( len( orb2irrep ) ):
        orb2irrep[ cnt ] = orb2irrep[ cnt ] % 10

    #######################################
    #   PySCF FCI calculation and 4-RDM   #
    #######################################
    fci_solver = fci.FCI(mol, mf.mo_coeff)
    fci_energy, fci_vector = fci_solver.kernel()
    print "PySCF     FCI energy =", fci_energy + CONST
    pyscf_start = time.time()
    dm1, dm2, dm3, dm4 = fci.rdm.make_dm1234('FCI4pdm_kern_spin0', fci_vector, fci_vector, L, N)
    dm1, dm2, dm3, dm4 = fci.rdm.reorder_dm1234(dm1, dm2, dm3, dm4, inplace=True)
    pyscf_trace = np.einsum( 'ii->', np.einsum( 'jjkl->kl', np.einsum( 'iijklm->jklm', np.einsum( 'iijklmno->jklmno', dm4 ) ) ) )
    pyscf_dm4 = np.array( dm4, copy=True )
    pyscf_dm4 = pyscf_dm4.swapaxes( 1, 2 ) # dm4[ 1, 5, 2, 6, 3, 7, 4, 8 ] -> dm4[ 1, 2, 5, 6, 3, 7, 4, 8 ]
    pyscf_dm4 = pyscf_dm4.swapaxes( 2, 4 ) # dm4[ 1, 2, 5, 6, 3, 7, 4, 8 ] -> dm4[ 1, 2, 3, 6, 5, 7, 4, 8 ]
    pyscf_dm4 = pyscf_dm4.swapaxes( 3, 6 ) # dm4[ 1, 2, 3, 6, 5, 7, 4, 8 ] -> dm4[ 1, 2, 3, 4, 5, 7, 6, 8 ]
    pyscf_dm4 = pyscf_dm4.swapaxes( 5, 6 ) # dm4[ 1, 2, 3, 4, 5, 7, 6, 8 ] -> dm4[ 1, 2, 3, 4, 5, 6, 7, 8 ]
    pyscf_end = time.time()

    #############################
    #   PyCheMPS2 Hamiltonian   #
    #############################
    Initializer = PyCheMPS2.PyInitialize()
    Initializer.Init()
    Ham = PyCheMPS2.PyHamiltonian( L, psi4group, orb2irrep, 'none' )
    Ham.setEconst( CONST )
    for cnt1 in range( L ):
        for cnt2 in range( L ):
            irrep12 = Ham.getOrbitalIrrep( cnt1 ) ^ Ham.getOrbitalIrrep( cnt2 )
            if ( irrep12 == 0 ):
                Ham.setTmat(cnt1, cnt2, OEImo[cnt1, cnt2])
            for cnt3 in range( L ):
                for cnt4 in range( L ):
                    irrep34 = Ham.getOrbitalIrrep( cnt3 ) ^ Ham.getOrbitalIrrep( cnt4 )
                    if ( irrep12 == irrep34 ):
                        Ham.setVmat(cnt1, cnt2, cnt3, cnt4, TEImo[cnt1, cnt3, cnt2, cnt4]) #From chemist to physics notation

    ###########################################
    #   PyCheMPS2 FCI calculation and 3-RDM   #
    ###########################################
    Irrep     = 0
    Nel_up    = ( N + mol.spin ) / 2
    Nel_down  = ( N - mol.spin ) / 2
    work_mem  = 25.0
    verbosity = 0
    theFCI    = PyCheMPS2.PyFCI(Ham, Nel_up, Nel_down, Irrep, work_mem, verbosity)
    GSvector  = np.zeros([ theFCI.getVecLength() ], dtype=ctypes.c_double)
    GSvector[ theFCI.LowestEnergyDeterminant() ] = 1.0
    EnergyFCI = theFCI.GSDavidson(GSvector)
    print "PyCheMPS2 FCI energy =", EnergyFCI
    start2    = time.time()
    FourRDM   = np.zeros([ L**8 ], dtype=ctypes.c_double)
    theFCI.Fill4RDM( GSvector, FourRDM )
    FourRDM   = FourRDM.reshape(L,L,L,L,L,L,L,L)
    print "PySCF     Tr(\'dm4\')  =", pyscf_trace
    print "PyCheMPS2 Tr(4-RDM)  =", np.einsum( 'ii->', np.einsum( 'jkjl->kl', np.einsum( 'ijkilm->jklm', np.einsum( 'ijkliqrs->jklqrs', FourRDM ) ) ) )
    end2      = time.time()
    print "Time 4-RDM PySCF     =", pyscf_end - pyscf_start, "seconds."
    print "Time 4-RDM PyCheMPS2 =", end2      - start2,      "seconds."

    ##############################
    #   Compare the two 4-RDMs   #
    ##############################
    for orb1 in range(L):
        for orb2 in range(L):
            for orb3 in range(L):
                for orb4 in range(L):
                    for orb5 in range(L):
                        for orb6 in range(L):
                            for orb7 in range(L):
                                for orb8 in range(L):
                                    temp = FourRDM[ orb1, orb2, orb3, orb4, orb5, orb6, orb7, orb8 ] - pyscf_dm4[ orb1, orb2, orb3, orb4, orb5, orb6, orb7, orb8 ]
                                    if abs( temp ) > 1e-5:
                                        print "4-RDM[",orb1,",",orb2,",",orb3,",",orb4,",",orb5,",",orb6,",",orb7,",",orb8,"] diff =", temp
    print "RMS difference PySCF and PyCheMPS2 4-RDM =", np.linalg.norm( FourRDM - pyscf_dm4 )



