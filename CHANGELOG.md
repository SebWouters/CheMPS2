#### Current HEAD
* Break API: bump up SO version (psi4 interface does NOT need to change!)
* Class TensorOperator: replaces TensorDiag, TensorSwap, TensorS0Abase, TensorS1Bbase, TensorF0Cbase, TensorF1Dbase, TensorA, TensorB, TensorC, and TensorD
* Deprecate TwoDMstorage class (storage handling occurs in TwoDM class)
* Further optimization allowed after RDM calculation
* 3-RDM in FCI class & two comparisons with pyscf in integrals/pyscf/debug_3rdm.py

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

