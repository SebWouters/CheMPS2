import sys
sys.path.append('/usr/lib/python2.7/site-packages') # Folder which contains PyCheMPS2.so
import PyCheMPS2
import ctypes
import numpy as np

def call_chemps2( Norbs, GroupName, Orbsym, CONST, OEINT, TEINT, TwoS, Nelec, Irrep, DSU2, Econv, MaxSweeps, NoisePrefac ):

    assert( len( DSU2 ) == len( Econv ) )
    assert( len( DSU2 ) == len( MaxSweeps) )
    assert( len( DSU2 ) == len( NoisePrefac ) )
    
    CheMPS2_groupnamemapping = {'C1' : 0, 'Ci' : 1, 'C2' : 2, 'Cs' : 3, 'D2' : 4, 'C2v': 5, 'C2h': 6, 'D2h': 7, 'Dooh': 7}
    Group = CheMPS2_groupnamemapping[ GroupName ]
    print "Point group of the molecule =",GroupName

    Initializer = PyCheMPS2.PyInitialize()
    Initializer.Init()

    Ham = PyCheMPS2.PyHamiltonian( Norbs, Group, np.array(Orbsym, dtype=ctypes.c_int) )
    Ham.setEconst( CONST )
    for orb1 in range(Norbs):
        for orb2 in range(Norbs):
            if ( Orbsym[orb1] ^ Orbsym[orb2] == 0 ):
                Ham.setTmat(orb1, orb2, OEINT[orb1, orb2])
            for orb3 in range(Norbs):
                for orb4 in range(Norbs):
                    if ( Orbsym[orb1] ^ Orbsym[orb2] ^ Orbsym[orb3] ^ Orbsym[orb4] == 0 ):
                        # CheMPS2 uses physics notation !
                        Ham.setVmat(orb1, orb2, orb3, orb4, TEINT[orb1, orb3, orb2, orb4])

    Prob = PyCheMPS2.PyProblem(Ham, TwoS, Nelec, Irrep)
    Prob.SetupReorderD2h() # Does not matter if the group is not d2h

    OptScheme = PyCheMPS2.PyConvergenceScheme( len( DSU2 ) )
    for instruction in range( len( DSU2 ) ):
        OptScheme.setInstruction(instruction, DSU2[instruction], Econv[instruction], MaxSweeps[instruction], NoisePrefac[instruction])

    theDMRG = PyCheMPS2.PyDMRG(Prob, OptScheme)
    EnergyDMRG = theDMRG.Solve()
    theDMRG.calc2DMandCorrelations()

    TwoRDM = np.zeros([Norbs, Norbs, Norbs, Norbs], dtype=ctypes.c_double)
    for orb1 in range(Norbs):
        for orb2 in range(Norbs):
            for orb3 in range(Norbs):
                for orb4 in range(Norbs):
                    # CheMPS2 uses physics notation !
                    TwoRDM[orb1, orb3, orb2, orb4] = theDMRG.get2DMA(orb1, orb2, orb3, orb4)
    
    theDMRG.deleteStoredMPS()
    theDMRG.deleteStoredOperators()
    del theDMRG
    del OptScheme
    del Prob
    del Ham
    del Initializer
    
    return ( EnergyDMRG, TwoRDM )

