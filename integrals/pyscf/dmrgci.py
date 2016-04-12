import sys
sys.path.append('/home/seba') #Folder in which PySCF is installed
from pyscf import gto, scf, symm, ao2mo
from call_chemps2 import call_chemps2
import numpy as np

def fetchJK_mo( mf , DM_mo , Cmat ):

   DM_ao = np.dot( np.dot( Cmat, DM_mo ), Cmat.T )
   JK_ao = scf.hf.get_veff( mf.mol, DM_ao, 0, 0, 1 ) #Last 3 numbers: dm_last, vhf_last, hermi
   JK_mo = np.dot( np.dot( Cmat.T, JK_ao ), Cmat )
   return JK_mo

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
    if (len(frozen) > 0):
        mo_fro = mf.mo_coeff[:, frozen]
    mo_act = mf.mo_coeff[:, active]
    Orbsym = np.array(symm.label_orb_symm(mf.mol, mf.mol.irrep_id, mf.mol.symm_orb, mf.mo_coeff))[active] # Same conventions in PySCF and CheMPS2
    for cnt in range( len( Orbsym ) ):
        # https://github.com/sunqm/pyscf/blob/master/symm/basis.py#L128
        Orbsym[ cnt ] = Orbsym[ cnt ] % 10
    
    ###  Get the JK matrix for the doubly occupied frozen orbitals  ###
    if (len(frozen) > 0):
        DM_docc = np.zeros([ Norbs ], dtype=float)
        DM_docc[ frozen ] = 2.0
        JK_docc = fetchJK_mo( mf, np.diag(DM_docc), mf.mo_coeff )
        JK_fro  = ( JK_docc[:, frozen] )[frozen, :]
        JK_act  = ( JK_docc[:, active] )[active, :]

    ###  Build the rotated MO matrix elements  ###
    Lrot  = len(active)
    OEIao = mf.mol.intor('cint1e_kin_sph') + mf.mol.intor('cint1e_nuc_sph')
    CONST = mf.mol.energy_nuc()
    OEImo = np.dot( np.dot( mo_act.T, OEIao ), mo_act )
    if (len(frozen) > 0):
        CONST += np.einsum( 'ii->', 2 * np.dot( np.dot( mo_fro.T, OEIao ), mo_fro ) + JK_fro )
        OEImo += JK_act
    #TEINT in chemical notation; also possible to pass different shapes; see http://www.sunqm.net/pyscf/tutorial.html#hf-mp2-mcscf
    TEImo = ao2mo.outcore.full_iofree(mf.mol, mo_act, compact=False).reshape(Lrot,Lrot,Lrot,Lrot)

    ###  Reorder per irrep  ### --> not strictly necessary for CheMPS2, but converges faster
    idx    = Orbsym.argsort()
    OEImo  = OEImo[:,idx]
    OEImo  = OEImo[idx,:]
    TEImo  = TEImo[:,:,:,idx]
    TEImo  = TEImo[:,:,idx,:]
    TEImo  = TEImo[:,idx,:,:]
    TEImo  = TEImo[idx,:,:,:]
    Orbsym = Orbsym[idx]

    ###  Do DMRG calculation with CheMPS2  ###
    EnergyDMRG, TwoRDM = call_chemps2( Lrot, mf.mol.groupname, Orbsym, CONST, OEImo, TEImo, TwoS, Nelec - 2*len(frozen), Irrep, DSU2, Econv, MaxSweeps, NoisePrefac )

    ###  Reorder back to original ordering  ###
    idx2 = idx.argsort()
    assert( np.linalg.norm( idx[idx2] - idx2[idx] ) < 1e-10 )
    TwoRDM = TwoRDM[:,:,:,idx2]
    TwoRDM = TwoRDM[:,:,idx2,:]
    TwoRDM = TwoRDM[:,idx2,:,:]
    TwoRDM = TwoRDM[idx2,:,:,:]

    return ( EnergyDMRG , TwoRDM )

