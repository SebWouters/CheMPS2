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
cimport numpy as np
np.import_array()

cimport Init
cimport ConvScheme
cimport Ham
cimport Prob
cimport Corr
cimport TwoRDM
cimport DMRGsolver
cimport DMRGSCFopt
cimport DMRGSCF
cimport FullCI

cdef class PyInitialize:
    cdef Init.Initialize * thisptr
    def __cinit__(self):
        self.thisptr = new Init.Initialize()
    def __dealloc__(self):
        del self.thisptr
    def Init(self):
        self.thisptr.Init()

cdef class PyConvergenceScheme:
    cdef ConvScheme.ConvergenceScheme * thisptr
    def __cinit__(self, int nInstructions):
        self.thisptr = new ConvScheme.ConvergenceScheme(nInstructions)
    def __dealloc__(self):
        del self.thisptr
    def setInstruction(self, int instruction, int D, double Econv, int nMax, double noisePrefactor):
        self.thisptr.setInstruction(instruction, D, Econv, nMax, noisePrefactor)
    def getD(self, int instruction):
        return self.thisptr.getD(instruction)
    def getEconv(self, int instruction):
        return self.thisptr.getEconv(instruction)
    def getMaxSweeps(self, int instruction):
        return self.thisptr.getMaxSweeps(instruction)
    def getNoisePrefactor(self, int instruction):
        return self.thisptr.getNoisePrefactor(instruction)

cdef class PyHamiltonian:
    cdef Ham.Hamiltonian * thisptr
    def __cinit__(self, int Norbitals, int nGroup, np.ndarray[int, ndim=1, mode="c"] OrbIrreps not None):
        assert OrbIrreps.flags['C_CONTIGUOUS']
        assert OrbIrreps.shape[0] == Norbitals
        self.thisptr = new Ham.Hamiltonian(Norbitals, nGroup, &OrbIrreps[0])
    def __dealloc__(self):
        del self.thisptr
    def getL(self):
        return self.thisptr.getL()
    def getNGroup(self):
        return self.thisptr.getNGroup()
    def getOrbitalIrrep(self, int orb):
        return self.thisptr.getOrbitalIrrep(orb)
    def setEconst(self, double value):
        self.thisptr.setEconst(value)
    def setTmat(self, int index1, int index2, double value):
        self.thisptr.setTmat(index1, index2, value)
    def setVmat(self, int index1, int index2, int index3, int index4, double value):
        self.thisptr.setVmat(index1, index2, index3, index4, value)
    def getEconst(self):
        return self.thisptr.getEconst()
    def getTmat(self, int index1, int index2):
        return self.thisptr.getTmat(index1, index2)
    def getVmat(self, int index1, int index2, int index3, int index4):
        return self.thisptr.getVmat(index1, index2, index3, index4)
    def save(self):
        self.thisptr.save()
    def read(self):
        self.thisptr.read()
        
cdef class PyProblem:
    cdef Prob.Problem * thisptr
    def __cinit__(self, PyHamiltonian Hami, int TwoS, int N, int Irrep):
        self.thisptr = new Prob.Problem(Hami.thisptr, TwoS, N, Irrep)
    def __dealloc__(self):
        del self.thisptr
    def gL(self):
        return self.thisptr.gL()
    def gSy(self):
        return self.thisptr.gSy()
    def gIrrep(self, int orb):
        return self.thisptr.gIrrep(orb)
    def gTwoS(self):
        return self.thisptr.gTwoS()
    def gN(self):
        return self.thisptr.gN()
    def gIrrep(self):
        return self.thisptr.gIrrep()
    def gEconst(self):
        return self.thisptr.gEconst()
    def gMxElement(self, int index1, int index2, int index3, int index4):
        return self.thisptr.gMxElement(index1, index2, index3, index4)
    def SetupReorderD2h(self):
        self.thisptr.SetupReorderD2h()
        
cdef class PyDMRG:
    cdef DMRGsolver.DMRG * thisptr
    def __cinit__(self, PyProblem Probl, PyConvergenceScheme OptScheme):
        self.thisptr = new DMRGsolver.DMRG(Probl.thisptr, OptScheme.thisptr)
    def __dealloc__(self):
        del self.thisptr
    def Solve(self):
        return self.thisptr.Solve()
    def calc2DMandCorrelations(self):
        self.thisptr.calc2DMandCorrelations()
    def deleteStoredMPS(self):
        self.thisptr.deleteStoredMPS()
    def deleteStoredOperators(self):
        self.thisptr.deleteStoredOperators()
    def activateExcitations(self, int nExcitations):
        self.thisptr.activateExcitations(nExcitations)
    def newExcitation(self, const double Eshift):
        self.thisptr.newExcitation(Eshift)
    #Access functions of the Corr.Correlations class
    def getCspin(self, int row, int col):
        return self.thisptr.getCorrelations().getCspin_HAM(row, col)
    def getCdens(self, int row, int col):
        return self.thisptr.getCorrelations().getCdens_HAM(row, col)
    def getCspinflip(self, int row, int col):
        return self.thisptr.getCorrelations().getCspinflip_HAM(row, col)
    def getCdirad(self, int row, int col):
        return self.thisptr.getCorrelations().getCdirad_HAM(row, col)
    def getMutInfo(self, int row, int col):
        return self.thisptr.getCorrelations().getMutualInformation_HAM(row, col)
    def getSingleOrbEntropy(self, int index):
        return self.thisptr.getCorrelations().SingleOrbitalEntropy_HAM(index)
    def getMutInfoDistance(self, int power):
        return self.thisptr.getCorrelations().MutualInformationDistance(power)
    def printCorrelations(self):
        self.thisptr.getCorrelations().Print()
    #Access functions of the TwoRDM.TwoDM class
    def get2DMA(self, int i1, int i2, int i3, int i4):
        return self.thisptr.get2DM().getTwoDMA_HAM(i1, i2, i3, i4)
    def get2DMB(self, int i1, int i2, int i3, int i4):
        return self.thisptr.get2DM().getTwoDMB_HAM(i1, i2, i3, i4)
    def get2DMenergy(self):
        return self.thisptr.get2DM().calcEnergy()
    def getDoubleTrace2DMA(self):
        return self.thisptr.get2DM().doubletrace2DMA()
    def save(self):
        self.thisptr.get2DM().save()
    def read(self):
        self.thisptr.get2DM().read()
    def getFCIcoefficient(self, np.ndarray[int, ndim=1, mode="c"] alpha not None, np.ndarray[int, ndim=1, mode="c"] beta not None):
        assert alpha.flags['C_CONTIGUOUS']
        assert  beta.flags['C_CONTIGUOUS']
        return self.thisptr.getFCIcoefficient(&alpha[0],&beta[0])

cdef class PyDMRGSCFoptions:
    cdef DMRGSCFopt.DMRGSCFoptions * thisptr
    def __cinit__(self):
        self.thisptr = new DMRGSCFopt.DMRGSCFoptions()
    def __dealloc__(self):
        del self.thisptr
    def getDoDIIS(self):
        return self.thisptr.getDoDIIS()
    def getDIISGradientBranch(self):
        return self.thisptr.getDIISGradientBranch()
    def getNumDIISVecs(self):
        return self.thisptr.getNumDIISVecs()
    def getStoreDIIS(self):
        return self.thisptr.getStoreDIIS()
    def getMaxIterations(self):
        return self.thisptr.getMaxIterations()
    def getGradientThreshold(self):
        return self.thisptr.getGradientThreshold()
    def getStoreUnitary(self):
        return self.thisptr.getStoreUnitary()
    def getWhichActiveSpace(self):
        return self.thisptr.getWhichActiveSpace()
    def getDumpCorrelations(self):
        return self.thisptr.getDumpCorrelations()
    def getStateAveraging(self):
        return self.thisptr.getStateAveraging()
    def setDoDIIS(self, bint val):
        self.thisptr.setDoDIIS(val)
    def setDIISGradientBranch(self, double val):
        self.thisptr.setDIISGradientBranch(val)
    def setNumDIISVecs(self, int val):
        self.thisptr.setNumDIISVecs(val)
    def setStoreDIIS(self, bint val):
        self.thisptr.setStoreDIIS(val)
    def setMaxIterations(self, int val):
        self.thisptr.setMaxIterations(val)
    def setGradientThreshold(self, double val):
        self.thisptr.setGradientThreshold(val)
    def setStoreUnitary(self, bint val):
        self.thisptr.setStoreUnitary(val)
    def setWhichActiveSpace(self, int val):
        self.thisptr.setWhichActiveSpace(val)
    def setDumpCorrelations(self, bint val):
        self.thisptr.setDumpCorrelations(val)
    def setStateAveraging(self, bint val):
        self.thisptr.setStateAveraging(val)
        
cdef class PyCASSCF:
    cdef DMRGSCF.CASSCF * thisptr
    def __cinit__(self, PyHamiltonian theHam, np.ndarray[int, ndim=1, mode="c"] DOCC not None, np.ndarray[int, ndim=1, mode="c"] SOCC not None):
        assert DOCC.flags['C_CONTIGUOUS']
        assert SOCC.flags['C_CONTIGUOUS']
        self.thisptr = new DMRGSCF.CASSCF(theHam.thisptr, &DOCC[0], &SOCC[0])
    def __dealloc__(self):
        del self.thisptr
    def setupStart(self, np.ndarray[int, ndim=1, mode="c"] Nocc not None, np.ndarray[int, ndim=1, mode="c"] NDMRG not None, np.ndarray[int, ndim=1, mode="c"] Nvirt not None):
        assert  Nocc.flags['C_CONTIGUOUS']
        assert NDMRG.flags['C_CONTIGUOUS']
        assert Nvirt.flags['C_CONTIGUOUS']
        self.thisptr.setupStart(&Nocc[0], &NDMRG[0], &Nvirt[0])
    def doCASSCFnewtonraphson(self, int Nel, int TwoS, int Irrep, PyConvergenceScheme OptScheme, int rootNum, PyDMRGSCFoptions theDMRGSCFopts):
        return self.thisptr.doCASSCFnewtonraphson(Nel, TwoS, Irrep, OptScheme.thisptr, rootNum, theDMRGSCFopts.thisptr)
    def deleteStoredUnitary(self):
        self.thisptr.deleteStoredUnitary()
    def deleteStoredDIIS(self):
        self.thisptr.deleteStoredDIIS()
        
cdef class PyFCI:
    cdef FullCI.FCI * thisptr
    def __cinit__(self, PyHamiltonian theHam, unsigned int Nel_up, unsigned int Nel_down, int TargetIrrep, double maxMemWorkMB=100.0, int FCIverbose=1):
        self.thisptr = new FullCI.FCI(theHam.thisptr, Nel_up, Nel_down, TargetIrrep, maxMemWorkMB, FCIverbose)
    def __dealloc__(self):
        del self.thisptr
    def getVecLength(self):
        return self.thisptr.getVecLength(0)
    def LowestEnergyDeterminant(self):
        return self.thisptr.LowestEnergyDeterminant()
    def GSDavidson(self, np.ndarray[double, ndim=1, mode="c"] inoutput not None):
        assert inoutput.flags['C_CONTIGUOUS']
        Energy = self.thisptr.GSDavidson(&inoutput[0])
        return Energy
    def CalcSpinSquared(self, np.ndarray[double, ndim=1, mode="c"] GSvector not None):
        assert GSvector.flags['C_CONTIGUOUS']
        SpinSquared = self.thisptr.CalcSpinSquared(&GSvector[0])
        return SpinSquared
    def Fill2RDM(self, np.ndarray[double, ndim=1, mode="c"] GSvector not None, np.ndarray[double, ndim=1, mode="c"] TwoRDM not None):
        assert GSvector.flags['C_CONTIGUOUS']
        assert   TwoRDM.flags['C_CONTIGUOUS']
        EnergyByContraction = self.thisptr.Fill2RDM(&GSvector[0], &TwoRDM[0])
        return EnergyByContraction
    def FillRandom(self, unsigned long long vecLength, np.ndarray[double, ndim=1, mode="c"] vector not None):
        assert vector.flags['C_CONTIGUOUS']
        self.thisptr.FillRandom(vecLength, &vector[0])
    def getFCIcoefficient(self, np.ndarray[int, ndim=1, mode="c"] alpha not None, np.ndarray[int, ndim=1, mode="c"] beta not None, np.ndarray[double, ndim=1, mode="c"] GSvector not None):
        assert    alpha.flags['C_CONTIGUOUS']
        assert     beta.flags['C_CONTIGUOUS']
        assert GSvector.flags['C_CONTIGUOUS']
        return self.thisptr.getFCIcoeff(&alpha[0], &beta[0], &GSvector[0])
    def RetardedGF(self, double omega, double eta, int orb_alpha, int orb_beta, bint isUp, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None, PyHamiltonian Hami):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        self.thisptr.RetardedGF(omega, eta, orb_alpha, orb_beta, isUp, GSenergy, &GSvector[0], Hami.thisptr, &RePart[0], &ImPart[0])
        return (RePart[0], ImPart[0])
    def RetardedGF_addition(self, double omega, double eta, int orb_alpha, int orb_beta, bint isUp, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None, PyHamiltonian Hami, np.ndarray[double, ndim=1, mode="c"] Re2RDM not None, np.ndarray[double, ndim=1, mode="c"] Im2RDM not None, np.ndarray[double, ndim=1, mode="c"] Add2RDM not None):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        assert   Re2RDM.flags['C_CONTIGUOUS']
        assert   Im2RDM.flags['C_CONTIGUOUS']
        assert  Add2RDM.flags['C_CONTIGUOUS']
        self.thisptr.RetardedGF_addition(omega, eta, orb_alpha, orb_beta, isUp, GSenergy, &GSvector[0], Hami.thisptr, &RePart[0], &ImPart[0], &Re2RDM[0], &Im2RDM[0], &Add2RDM[0])
        return (RePart[0], ImPart[0])
    def RetardedGF_removal(self, double omega, double eta, int orb_alpha, int orb_beta, bint isUp, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None, PyHamiltonian Hami, np.ndarray[double, ndim=1, mode="c"] Re2RDM not None, np.ndarray[double, ndim=1, mode="c"] Im2RDM not None, np.ndarray[double, ndim=1, mode="c"] Rem2RDM not None):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        assert   Re2RDM.flags['C_CONTIGUOUS']
        assert   Im2RDM.flags['C_CONTIGUOUS']
        assert  Rem2RDM.flags['C_CONTIGUOUS']
        self.thisptr.RetardedGF_removal(omega, eta, orb_alpha, orb_beta, isUp, GSenergy, &GSvector[0], Hami.thisptr, &RePart[0], &ImPart[0], &Re2RDM[0], &Im2RDM[0], &Rem2RDM[0])
        return (RePart[0], ImPart[0])
    def GFmatrix_add(self, double alpha, double beta, double eta, np.ndarray[int, ndim=1, mode="c"] orbsLeft not None, np.ndarray[int, ndim=1, mode="c"] orbsRight not None, bint isUp, np.ndarray[double, ndim=1, mode="c"] GSvector not None, PyHamiltonian Hami):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([len(orbsLeft)*len(orbsRight)])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([len(orbsLeft)*len(orbsRight)])
        assert  GSvector.flags['C_CONTIGUOUS']
        assert    RePart.flags['C_CONTIGUOUS']
        assert    ImPart.flags['C_CONTIGUOUS']
        assert  orbsLeft.flags['C_CONTIGUOUS']
        assert orbsRight.flags['C_CONTIGUOUS']
        self.thisptr.GFmatrix_addition(alpha, beta, eta, &orbsLeft[0], len(orbsLeft), &orbsRight[0], len(orbsRight), isUp, &GSvector[0], Hami.thisptr, &RePart[0], &ImPart[0])
        return ( RePart, ImPart )
    def GFmatrix_rem(self, double alpha, double beta, double eta, np.ndarray[int, ndim=1, mode="c"] orbsLeft not None, np.ndarray[int, ndim=1, mode="c"] orbsRight not None, bint isUp, np.ndarray[double, ndim=1, mode="c"] GSvector not None, PyHamiltonian Hami):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([len(orbsLeft)*len(orbsRight)])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([len(orbsLeft)*len(orbsRight)])
        assert  GSvector.flags['C_CONTIGUOUS']
        assert    RePart.flags['C_CONTIGUOUS']
        assert    ImPart.flags['C_CONTIGUOUS']
        assert  orbsLeft.flags['C_CONTIGUOUS']
        assert orbsRight.flags['C_CONTIGUOUS']
        self.thisptr.GFmatrix_removal(alpha, beta, eta, &orbsLeft[0], len(orbsLeft), &orbsRight[0], len(orbsRight), isUp, &GSvector[0], Hami.thisptr, &RePart[0], &ImPart[0])
        return ( RePart, ImPart )
    def DensityResponseGF(self, double omega, double eta, int orb_alpha, int orb_beta, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        self.thisptr.DensityResponseGF(omega, eta, orb_alpha, orb_beta, GSenergy, &GSvector[0], &RePart[0], &ImPart[0])
        return (RePart[0], ImPart[0])
    def DensityResponseGF_forward(self, double omega, double eta, int orb_alpha, int orb_beta, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None, np.ndarray[double, ndim=1, mode="c"] Re2RDM not None, np.ndarray[double, ndim=1, mode="c"] Im2RDM not None, np.ndarray[double, ndim=1, mode="c"] Dens2RDM not None):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        assert   Re2RDM.flags['C_CONTIGUOUS']
        assert   Im2RDM.flags['C_CONTIGUOUS']
        assert Dens2RDM.flags['C_CONTIGUOUS']
        self.thisptr.DensityResponseGF_forward(omega, eta, orb_alpha, orb_beta, GSenergy, &GSvector[0], &RePart[0], &ImPart[0], &Re2RDM[0], &Im2RDM[0], &Dens2RDM[0])
        return (RePart[0], ImPart[0])
    def DensityResponseGF_backward(self, double omega, double eta, int orb_alpha, int orb_beta, double GSenergy, np.ndarray[double, ndim=1, mode="c"] GSvector not None, np.ndarray[double, ndim=1, mode="c"] Re2RDM not None, np.ndarray[double, ndim=1, mode="c"] Im2RDM not None, np.ndarray[double, ndim=1, mode="c"] Dens2RDM not None):
        cdef np.ndarray[double, ndim=1, mode="c"] RePart = np.zeros([1])
        cdef np.ndarray[double, ndim=1, mode="c"] ImPart = np.zeros([1])
        assert GSvector.flags['C_CONTIGUOUS']
        assert   RePart.flags['C_CONTIGUOUS']
        assert   ImPart.flags['C_CONTIGUOUS']
        assert   Re2RDM.flags['C_CONTIGUOUS']
        assert   Im2RDM.flags['C_CONTIGUOUS']
        assert Dens2RDM.flags['C_CONTIGUOUS']
        self.thisptr.DensityResponseGF_backward(omega, eta, orb_alpha, orb_beta, GSenergy, &GSvector[0], &RePart[0], &ImPart[0], &Re2RDM[0], &Im2RDM[0], &Dens2RDM[0])
        return (RePart[0], ImPart[0])


