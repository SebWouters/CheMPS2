#### Version 1.8.9 LTS update (2018-10-29):
* Bug fix long double literal Wigner.cpp

#### Version 1.8.8 LTS update (2018-09-26):
* Print the number of MPS variables

#### Version 1.8.7 LTS update (2018-03-25):
* Bump up SO version from 2 to 3

#### Version 1.8.6 LTS update (2018-03-16):
* Set initial guess to an ROHF determinant (for @quanp)

#### Version 1.8.5 LTS update (2018-01-18):
* Fix build; bump up required CMake (thanks @quanp)

#### Version 1.8.4 LTS update (2017-05-23):
* MOLCAS state average CASPT2 interface (thanks @quanp)
* Update build system for psi4 (thanks @loriab)

#### Version 1.8.3 LTS update (2016-11-15):
* Custom orbital reordering for MOLCAS interface
* Fix building PyCheMPS2 (thanks @loriab)

#### Version 1.8.2 LTS update (2016-10-26):
* Fix PyCheMPS2 build

#### Version 1.8.1 LTS update (2016-10-25):
* MPS checkpoint for MOLCAS interface
* Fiedler order checkpoint for MOLCAS interface

#### Version 1.8 LTS (2016-08-24):
* Fix slow convergence linear Davidson algorithm
* Fix MPI bug binary
* MOLCAS interface binary (2-RDM, 3-RDM, and F.4-RDM dump)
* Fortran example to read in CASPT2 intermediates for MOLCAS

#### Version 1.7.3 (2016-07-14):
* DMRG workshop in sphinx documentation

#### Version 1.7.2 (2016-07-07):
* Rotate molden files with DMRG-SCF unitary from executable

#### Version 1.7.1 (2016-06-29):
* Global reordering of the active space orbitals (mixing
  irreps) based on Fiedler vector of the exchange matrix
* Davidson algorithm for linear systems (for CASPT2)

#### Version 1.7 (2016-05-24):
* Class ConjugateGradient
* Class TensorOperator
* Class ThreeDM: DMRG 3-RDM in O(L^4 D^3 + L^6 D^2)
* Class Cumulant: Cumulant reconstruction of the DMRG 4-RDM
* Class CASPT2: Internally contracted CASPT2
* Class Molden: Rotate molpro/psi4 molden to CAS orbitals
* Class Excitation: Build |1> = ( b + a * ( E_ij + E_ji )) |0>
* Part of DMRG 4-RDM via DMRG::Symm4RDM
* DMRG sweeps after RDM calculation allowed
* FCI 3-RDM in FCI::Fill3RDM
* FCI 4-RDM in FCI::Fill4RDM
* FCI diagonal 4-RDM in FCI::Diag4RDM
* FCI 4-RDM contraction with Fock operator in FCI::Fock4RDM
* Blockwise ERI rotations with DMRGSCFrotations using disk
* Update binary for DMRG-CASPT2
* Break API: bump up SO version
* Deprecate TwoDMstorage class
* Deprecate TensorDiag class
* Deprecate TensorSwap class
* Deprecate TensorS0Abase class
* Deprecate TensorS1Bbase class
* Deprecate TensorF0Cbase class
* Deprecate TensorF1Dbase class
* Deprecate TensorA class
* Deprecate TensorB class
* Deprecate TensorC class
* Deprecate TensorD class

#### Version 1.6 (2015-08-26):
* Disk i/o improvements with HDF5's hyperslab
* Performance counters in DMRG class
* Faster preconditioner FCI Green's function solver
* Bug fix FCIDUMP read-in
* chemps2 binary
* manpage for binary

#### Version 1.5 (2015-06-18):
* DMRG-CI plugin pyscf
* DMRG-SCF plugin psi4 (official release for psi4)
* Fix bug small electron number FCI class
* FCIDUMP file support
* Sphinx documentation
* DMRG class supports 4-fold permutation symmetry (i.o. 8-fold)
* Hybrid MPI & OpenMP for DMRG + 2-RDM (not for DMRG-SCF yet)

#### Version 1.4 (2014-11-23):
* 2-RDM storage class
* Optimization OpenMP over symmetry blocks Heff
* DIIS acceleration DMRG-SCF
* Two-orbital mutual information and correlation functions
* Augmented Hessian NR Edmiston-Ruedenberg orbital localization
* PyCheMPS2: python interface to libchemps2
* FCI Green's function solver
* State-averaged DMRG-SCF

#### Version 1.0 (2014-04-08):
* Initial release

