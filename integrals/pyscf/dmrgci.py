import sys
sys.path.append('/home/seba') #Folder in which PySCF is installed
from pyscf import gto, scf, symm, ao2mo
from call_chemps2 import call_chemps2
import numpy as np

def dmrgci( mf , TwoS , Nelec , Irrep , DSU2 , Econv , MaxSweeps , NoisePrefac , frozen , active ):

    ###  Check whether the frozen and active arrays make sense  ###
    Norbs = mf.mol.nao_nr()
    assert( np.sum( frozen < Norbs ) == len(frozen) )
    assert( np.sum( active < Norbs ) == len(active) )
    assert( np.sum( frozen >= 0 ) == len(frozen) )
    assert( np.sum( active >= 0 ) == len(active) )
    mo_occ = mf.mo_occ
    for item in frozen:
        assert( mo_occ[item] == 2 )
        for item2 in active:
            assert( item != item2 )
    frozen = frozen[ frozen.argsort() ]
    active = active[ active.argsort() ]

    ###  Get orbital information  ###
    if ( len(frozen) == 0 ):
        torotate = active
    else:
        torotate = np.concatenate((frozen, active))
    mo_coeff = mf.mo_coeff[:,torotate]
    Orbsym   = np.array(symm.label_orb_symm(mf.mol, mf.mol.irrep_id, mf.mol.symm_orb, mf.mo_coeff))[torotate] # Same conventions in PySCF and CheMPS2

    ###  Build the rotated MO matrix elements  ###
    Lrot  = len(torotate)
    CONST = mf.mol.energy_nuc()
    OEINT = np.dot(np.dot( mo_coeff.T, mf.mol.intor('cint1e_kin_sph') + mf.mol.intor('cint1e_nuc_sph') ), mo_coeff )
    #TEINT in chemical notation; also possible to pass different shapes; see http://www.sunqm.net/pyscf/tutorial.html#hf-mp2-mcscf
    TEINT = ao2mo.outcore.full_iofree(mf.mol, mo_coeff, compact=False).reshape(Lrot,Lrot,Lrot,Lrot)
    
    ###  Fold in the frozen core electrons  ###
    for orb in range(len(frozen)):
        CONST += 2 * OEINT[orb, orb]
        for orb2 in range(len(frozen)):
            CONST += 2 * TEINT[orb, orb, orb2, orb2] - TEINT[orb, orb2, orb, orb2]
        for orb2 in range(len(frozen), len(torotate)):
            for orb3 in range(len(frozen), len(torotate)):
                OEINT[orb2, orb3] += 2 * TEINT[orb2, orb3, orb, orb] - TEINT[orb2, orb, orb3, orb]
            
    ###  Cut out the relevant integrals  ###
    OEINT  = OEINT[len(frozen):,len(frozen):]
    TEINT  = TEINT[len(frozen):,len(frozen):,len(frozen):,len(frozen):]
    Orbsym = Orbsym[len(frozen):]

    ###  Reorder per irrep  ### --> not strictly necessary for CheMPS2, but converges faster
    idx    = Orbsym.argsort()
    OEINT  = OEINT[:,idx]
    OEINT  = OEINT[idx,:]
    TEINT  = TEINT[:,:,:,idx]
    TEINT  = TEINT[:,:,idx,:]
    TEINT  = TEINT[:,idx,:,:]
    TEINT  = TEINT[idx,:,:,:]
    Orbsym = Orbsym[idx]

    ###  Do DMRG calculation with CheMPS2  ###
    EnergyDMRG, TwoRDM = call_chemps2( len(Orbsym), mf.mol.groupname, Orbsym, CONST, OEINT, TEINT, TwoS, Nelec - 2*len(frozen), Irrep, DSU2, Econv, MaxSweeps, NoisePrefac )

    ###  Reorder back to original ordering  ###
    idx2 = idx.argsort()
    assert( np.linalg.norm( idx[idx2] - idx2[idx] ) < 1e-10 )
    TwoRDM = TwoRDM[:,:,:,idx2]
    TwoRDM = TwoRDM[:,:,idx2,:]
    TwoRDM = TwoRDM[:,idx2,:,:]
    TwoRDM = TwoRDM[idx2,:,:,:]

    return ( EnergyDMRG , TwoRDM )

