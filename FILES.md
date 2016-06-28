List of files in CheMPS2
------------------------

[CheMPS2/CASPT2.cpp](CheMPS2/CASPT2.cpp) contains an implementation of
internally contracted CASPT2. The user can specify an IPEA shift and/or
an imaginary level shift to mitigate possible intruder state problems.
The linear CASPT2 equations can be solved with either Davidson's
algorithm, or with the conjugate gradient method. Note that the overlap
matrix is always diagonalized, and that a cumulant approximation of the
4-RDM should therefore be avoided.

[CheMPS2/CASSCF.cpp](CheMPS2/CASSCF.cpp) contains the functionality
to construct the active space Hamiltonian for DMRG-SCF.

[CheMPS2/CASSCFdebug.cpp](CheMPS2/CASSCFdebug.cpp) contains two
functions: one for calculating the ROHF energy; and one for fetching FCI
coefficients to determine the point group symmetry of certain electronic
states of the iron dimer.

[CheMPS2/CASSCFnewtonraphson.cpp](CheMPS2/CASSCFnewtonraphson.cpp)
contains all DMRG-SCF functions which are specific to the augmented Hessian
Newton-Raphson update scheme, including the functions to calculate the
gradient and Hessian.

[CheMPS2/CASSCFpt2.cpp](CheMPS2/CASSCFpt2.cpp) provides the interface between
the CASSCF and CASPT2 classes. The routines for the 3-RDM and the Fock
operator contraction with the 4-RDM are called here.

[CheMPS2/ConjugateGradient.cpp](CheMPS2/ConjugateGradient.cpp) is an
implementation of the conjugate gradient algorithm, in the style of the
Davidson class.

[CheMPS2/ConvergenceScheme.cpp](CheMPS2/ConvergenceScheme.cpp) contains
all functions of the ConvergenceScheme class, which contains the instructions
for the subsequent DMRG sweeps.

[CheMPS2/Correlations.cpp](CheMPS2/Correlations.cpp) contains all the
functionality to calculate the spin, density, and spin-flip correlation
functions as well as the two-orbital mutual information.

[CheMPS2/Cumulant.cpp](CheMPS2/Cumulant.cpp) contains static member functions
to reconstruct the 4-RDM from lower order reduced density matrices. There is
also a fast contraction of the Fock operator with the cumulant-reconstructed
4-RDM.

[CheMPS2/Davidson.cpp](CheMPS2/Davidson.cpp) is an implementation of
Davidson's algorithm, both for eigenvalue problems and linear equations.

[CheMPS2/DIIS.cpp](CheMPS2/DIIS.cpp) contains a DIIS convergence
speed-up for DMRG-SCF.

[CheMPS2/DMRG.cpp](CheMPS2/DMRG.cpp) contains the constructor and
destructor of the DMRG class, as well as the top-level sweep functions.

[CheMPS2/DMRGfock.cpp](CheMPS2/DMRGfock.cpp) contains the functionality to
express a symmetry (spin, particle number, and point group) conserving
single-particle excitation on top of an MPS as a new MPS.

[CheMPS2/DMRGmpsio.cpp](CheMPS2/DMRGmpsio.cpp) contains the store and
load functions for the DMRG checkpoint file (the MPS and the SyBookkeeper).

[CheMPS2/DMRGoperators3RDM.cpp](CheMPS2/DMRGoperators3RDM.cpp) contains all
update functions for the renormalized operators specific for the ThreeDM and
the Correlations.

[CheMPS2/DMRGoperators.cpp](CheMPS2/DMRGoperators.cpp) contains all
functions related to the DMRG renormalized operators: saving to disk,
loading from disk, and updating.

[CheMPS2/DMRGSCFindices.cpp](CheMPS2/DMRGSCFindices.cpp) contains the
index conversions for the DMRG-SCF algorithm.

[CheMPS2/DMRGSCFintegrals.cpp](CheMPS2/DMRGSCFintegrals.cpp) is a container
class for two-body matrix elements with at most two virtual indices.

[CheMPS2/DMRGSCFmatrix.cpp](CheMPS2/DMRGSCFmatrix.cpp) is a container
class for orbital matrices which are blockdiagonal in the irreps.

[CheMPS2/DMRGSCFoptions.cpp](CheMPS2/DMRGSCFoptions.cpp) is a container
class to pass the DMRG-SCF options to the augmented Hessian Newton-Raphson
and CASPT2 routines of the CASSCF class.

[CheMPS2/DMRGSCFrotations.cpp](CheMPS2/DMRGSCFrotations.cpp)
contains static member functions for the two-body matrix element rotations
for the CASSCF and Edmiston-Ruedenberg classes.

[CheMPS2/DMRGSCFunitary.cpp](CheMPS2/DMRGSCFunitary.cpp) contains the
storage and handling of the unitary matrix and its nonredundant
skew-symmetric parametrization required for the DMRG-SCF algorithm.

[CheMPS2/DMRGSCFwtilde.cpp](CheMPS2/DMRGSCFwtilde.cpp) is a container
class to store an intermediate for the DMRG-SCF Hessian.

[CheMPS2/DMRGtechnics.cpp](CheMPS2/DMRGtechnics.cpp) contains the
functions related to the RDM and excited-state calculations.

[CheMPS2/EdmistonRuedenberg.cpp](CheMPS2/EdmistonRuedenberg.cpp) contains
an orbital localization function based on the Edmiston-Ruedenberg cost function
and an augmented Hessian Newton-Raphson optimizer. It also contains the
functionality to compute the Fiedler vector of the exchange matrix, to reorder
the active space orbitals in a black-box fashion.

[CheMPS2/Excitation.cpp](CheMPS2/Excitation.cpp) contains matrix-vector
multiplication routines for spin-conserving single-particle excitations.

[CheMPS2/FCI.cpp](CheMPS2/FCI.cpp) contains a fast determinant-based full
configuration interaction (FCI) solver. The eigenvalue problem is solved with
Davidson's algorithm. Green's functions can also be computed.

[CheMPS2/FourIndex.cpp](CheMPS2/FourIndex.cpp) contains all functions of
the FourIndex container class for the two-body matrix elements.

[CheMPS2/Hamiltonian.cpp](CheMPS2/Hamiltonian.cpp) contains all functions
of the Hamiltonian class, including functions to get or set specific variables,
as well as to save and load the Hamiltonian on disk.

[CheMPS2/Heff.cpp](CheMPS2/Heff.cpp) contains top-level functions to perform
the DMRG effective Hamiltonian times vector multiplication for Davidson's
algorithm.

[CheMPS2/HeffDiagonal.cpp](CheMPS2/HeffDiagonal.cpp) contains the
functions to calculate the diagonal elements of the effective Hamiltonian.
These are required as preconditioner in Davidson's algorithm.
    
[CheMPS2/HeffDiagrams1.cpp](CheMPS2/HeffDiagrams1.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess vector
multiplication.

[CheMPS2/HeffDiagrams2.cpp](CheMPS2/HeffDiagrams2.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess vector
multiplication.

[CheMPS2/HeffDiagrams3.cpp](CheMPS2/HeffDiagrams3.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess vector
multiplication.

[CheMPS2/HeffDiagrams4.cpp](CheMPS2/HeffDiagrams4.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess vector
multiplication.

[CheMPS2/HeffDiagrams5.cpp](CheMPS2/HeffDiagrams5.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess vector
multiplication.

[CheMPS2/Initialize.cpp](CheMPS2/Initialize.cpp) sets the seed
of the random number generator and cout.precision (added for PyCheMPS2).

[CheMPS2/Irreps.cpp](CheMPS2/Irreps.cpp) contains the psi4 symmetry
conventions.

[CheMPS2/Molden.cpp](CheMPS2/Molden.cpp) contains the functionality to
rotate an R(O)HF molden file generated by molpro or psi4 to the new CAS space
defined by the DMRGSCFunitary HDF5 checkpoint file.

[CheMPS2/PrintLicense.cpp](CheMPS2/PrintLicense.cpp) contains a function
which prints the license disclaimer.

[CheMPS2/Problem.cpp](CheMPS2/Problem.cpp) contains all Problem class
functions. This wrapper class allows to set the desired symmetry sector for
the DMRG algorithm. It allows to compute fourfold permutation symmetric
Hamiltonians with DMRG (see tests 9 and 12).

[CheMPS2/Sobject.cpp](CheMPS2/Sobject.cpp) contains all Sobject class
functions. This class constructs, stores, and decomposes the reduced two-site
object.

[CheMPS2/SyBookkeeper.cpp](CheMPS2/SyBookkeeper.cpp) contains all
SyBookkeeper functions. This class keeps track of the FCI and DMRG virtual
dimensions of all symmetry sectors at all boundaries.

[CheMPS2/Tensor3RDM.cpp](CheMPS2/Tensor3RDM.cpp) contains all
initialization functions for the spin-reduced renormalized
operators of three second quantized operators.

[CheMPS2/TensorF0.cpp](CheMPS2/TensorF0.cpp) contains all TensorF0
functions. This class stores and handles the reduced spin-0 part of
two sandwiched second quantized operators, of which the particle
symmetry sectors are equal.

[CheMPS2/TensorF1.cpp](CheMPS2/TensorF1.cpp) contains all TensorF1
functions. This class stores and handles the reduced spin-1 part of
two sandwiched second quantized operators, of which the particle
symmetry sectors are equal.

[CheMPS2/TensorGYZ.cpp](CheMPS2/TensorGYZ.cpp) contains the contruct
and update functions for the G-, Y-, and Z-tensors. They are required
for the two-orbital mutual information.

[CheMPS2/TensorKM.cpp](CheMPS2/TensorKM.cpp) contains the contruct and
update functions for the K- and M-tensors. It is required for the
two-orbital mutual information.

[CheMPS2/TensorL.cpp](CheMPS2/TensorL.cpp) contains all TensorL functions.
This class stores and handles the reduced spin-1/2 part of a single
sandwiched second quantized operator. The class has been updated to allow
for a different bra and ket wavefunction, which is needed for
[CheMPS2/DMRGfock.cpp](CheMPS2/DMRGfock.cpp).

[CheMPS2/TensorO.cpp](CheMPS2/TensorO.cpp) handles the tensors required
to calculate the overlap between two MPSs.

[CheMPS2/TensorOperator.cpp](CheMPS2/TensorOperator.cpp) implements the
storage and handling of tensor operators with a given spin, particle
number, and point group irrep. It replaces the deprecated TensorDiag,
TensorSwap, TensorS0Abase, TensorS1Bbase, TensorF0Cbase, TensorF1Dbase,
TensorA, TensorB, TensorC, and TensorD classes.

[CheMPS2/TensorQ.cpp](CheMPS2/TensorQ.cpp) contains all TensorQ functions.
This class stores and handles the complementary reduced spin-1/2 part of
three sandwiched second quantized operators.

[CheMPS2/TensorS0.cpp](CheMPS2/TensorS0.cpp) contains all TensorS0
functions. This class stores and handles the reduced spin-0 part of
two sandwiched second quantized operators, of which the particle symmetry
sectors differ by 2.

[CheMPS2/TensorS1.cpp](CheMPS2/TensorS1.cpp) contains all TensorS1
functions. This class stores and handles the reduced spin-1 part of
two sandwiched second quantized operators, of which the particle symmetry
sectors differ by 2.

[CheMPS2/TensorT.cpp](CheMPS2/TensorT.cpp) contains all TensorT functions.
This class stores and handles the reduced part of an MPS site-tensor. It
also contains the functionality for QR- and LQ-decomposition of MPS tensors.

[CheMPS2/TensorX.cpp](CheMPS2/TensorX.cpp) contains all TensorX functions.
This class stores and handles the complementary reduced spin-0 part of four
sandwiched second quantized operators, which is of course diagonal in the
symmetry sectors.

[CheMPS2/ThreeDM.cpp](CheMPS2/ThreeDM.cpp) contains all functions to calculate
and store the 3-RDM from the DMRG-optimized MPS.

[CheMPS2/TwoDM.cpp](CheMPS2/TwoDM.cpp) contains all functions to calculate
and store the 2-RDM from the DMRG-optimized MPS.

[CheMPS2/TwoIndex.cpp](CheMPS2/TwoIndex.cpp) contains all functions of the
TwoIndex container class for the one-body matrix elements.

[CheMPS2/Wigner.cpp](CheMPS2/Wigner.cpp) contains static member functions
to compute Wigner 3j, 6j, and 9j symbols. The API has been chosen to match
GSL's gsl_sf_coupling_3j, gsl_sf_coupling_6j, and gsl_sf_coupling_9j.

[CheMPS2/executable.cpp](CheMPS2/executable.cpp) builds to the chemps2
executable, which allows to use libchemps2 from the command line.

[CheMPS2/include/chemps2/CASPT2.h](CheMPS2/include/chemps2/CASPT2.h) contains the definitions of the CASPT2 class.

[CheMPS2/include/chemps2/CASSCF.h](CheMPS2/include/chemps2/CASSCF.h) contains the definitions of the CASSCF class.

[CheMPS2/include/chemps2/ConjugateGradient.h](CheMPS2/include/chemps2/ConjugateGradient.h) contains the definitions of the ConjugateGradient class.

[CheMPS2/include/chemps2/ConvergenceScheme.h](CheMPS2/include/chemps2/ConvergenceScheme.h) contains the definitions of the ConvergenceScheme class.

[CheMPS2/include/chemps2/Correlations.h](CheMPS2/include/chemps2/Correlations.h) contains the definitions of the Correlations class.

[CheMPS2/include/chemps2/Cumulant.h](CheMPS2/include/chemps2/Cumulant.h) contains the definitions of the Cumulant class.

[CheMPS2/include/chemps2/Davidson.h](CheMPS2/include/chemps2/Davidson.h) contains the definitions of the Davidson class.

[CheMPS2/include/chemps2/DIIS.h](CheMPS2/include/chemps2/DIIS.h) contains the definitions of the DIIS class.

[CheMPS2/include/chemps2/DMRG.h](CheMPS2/include/chemps2/DMRG.h) contains the definitions of the DMRG class.

[CheMPS2/include/chemps2/DMRGSCFindices.h](CheMPS2/include/chemps2/DMRGSCFindices.h) contains the definitions of the DMRGSCFindices class.

[CheMPS2/include/chemps2/DMRGSCFintegrals.h](CheMPS2/include/chemps2/DMRGSCFintegrals.h) contains the definitions of the DMRGSCFintegrals class.

[CheMPS2/include/chemps2/DMRGSCFmatrix.h](CheMPS2/include/chemps2/DMRGSCFmatrix.h) contains the definitions of the DMRGSCFmatrix class.

[CheMPS2/include/chemps2/DMRGSCFoptions.h](CheMPS2/include/chemps2/DMRGSCFoptions.h) contains the definitions of the DMRGSCFoptions container class.

[CheMPS2/include/chemps2/DMRGSCFrotations.h](CheMPS2/include/chemps2/DMRGSCFrotations.h) contains the definitions of the DMRGSCFrotations class.

[CheMPS2/include/chemps2/DMRGSCFunitary.h](CheMPS2/include/chemps2/DMRGSCFunitary.h) contains the definitions of the DMRGSCFunitary class.

[CheMPS2/include/chemps2/DMRGSCFwtilde.h](CheMPS2/include/chemps2/DMRGSCFwtilde.h) contains the definitions of the DMRGSCFwtilde class.

[CheMPS2/include/chemps2/EdmistonRuedenberg.h](CheMPS2/include/chemps2/EdmistonRuedenberg.h) contains the definitions of the EdmistonRuedenberg class.

[CheMPS2/include/chemps2/Excitation.h](CheMPS2/include/chemps2/Excitation.h) contains the definitions of the Excitation class.

[CheMPS2/include/chemps2/FCI.h](CheMPS2/include/chemps2/FCI.h) contains the definitions of the FCI class.

[CheMPS2/include/chemps2/FourIndex.h](CheMPS2/include/chemps2/FourIndex.h) contains the definitions of the FourIndex class.

[CheMPS2/include/chemps2/Hamiltonian.h](CheMPS2/include/chemps2/Hamiltonian.h) contains the definitions of the Hamiltonian class.

[CheMPS2/include/chemps2/Heff.h](CheMPS2/include/chemps2/Heff.h) contains the definitions of the Heff class.

[CheMPS2/include/chemps2/Initialize.h](CheMPS2/include/chemps2/Initialize.h) contains the definitions of the Initialize class.

[CheMPS2/include/chemps2/Irreps.h](CheMPS2/include/chemps2/Irreps.h) contains the definitions of the Irreps class.

[CheMPS2/include/chemps2/Lapack.h](CheMPS2/include/chemps2/Lapack.h) contains the definitions of the external BLAS and LAPACK routines.

[CheMPS2/include/chemps2/Molden.h](CheMPS2/include/chemps2/Molden.h) contains the definitions of the Molden class.

[CheMPS2/include/chemps2/MPIchemps2.h](CheMPS2/include/chemps2/MPIchemps2.h)
contains the distribution of (complementary) renormalized operators over MPI
processes, as well as wrappers for the MPI communication routines in the C API.

[CheMPS2/include/chemps2/MyHDF5.h](CheMPS2/include/chemps2/MyHDF5.h) forces
the use of the HDF5 1.8 API, e.g. H5Gcreate2 instead of H5Gcreate1, a known
issue in Ubuntu 12.04.

[CheMPS2/include/chemps2/Options.h](CheMPS2/include/chemps2/Options.h)
contains all the options of the CheMPS2 namespace. Here the checkpoint
storage names and folders can be set, as well as parameters related to
memory usage and convergence.

[CheMPS2/include/chemps2/Problem.h](CheMPS2/include/chemps2/Problem.h) contains the definitions of the Problem class.

[CheMPS2/include/chemps2/Sobject.h](CheMPS2/include/chemps2/Sobject.h) contains the definitions of the Sobject class.

[CheMPS2/include/chemps2/Special.h](CheMPS2/include/chemps2/Special.h) contains special functions needed in various parts of libchemps2.

[CheMPS2/include/chemps2/SyBookkeeper.h](CheMPS2/include/chemps2/SyBookkeeper.h) contains the definitions of the SyBookkeeper class.

[CheMPS2/include/chemps2/Tensor3RDM.h](CheMPS2/include/chemps2/Tensor3RDM.h) contains the definitions of the Tensor3RDM class.

[CheMPS2/include/chemps2/TensorF0.h](CheMPS2/include/chemps2/TensorF0.h) contains the definitions of the TensorF0 class.

[CheMPS2/include/chemps2/TensorF1.h](CheMPS2/include/chemps2/TensorF1.h) contains the definitions of the TensorF1 class.

[CheMPS2/include/chemps2/TensorGYZ.h](CheMPS2/include/chemps2/TensorGYZ.h) contains the definitions of the TensorGYZ class.

[CheMPS2/include/chemps2/Tensor.h](CheMPS2/include/chemps2/Tensor.h) contains the definitions of the virtual Tensor class.

[CheMPS2/include/chemps2/TensorKM.h](CheMPS2/include/chemps2/TensorKM.h) contains the definitions of the TensorKM class.

[CheMPS2/include/chemps2/TensorL.h](CheMPS2/include/chemps2/TensorL.h) contains the definitions of the TensorL class.

[CheMPS2/include/chemps2/TensorO.h](CheMPS2/include/chemps2/TensorO.h) contains the definitions of the TensorO class.

[CheMPS2/include/chemps2/TensorOperator.h](CheMPS2/include/chemps2/TensorOperator.h) contains the definitions of the TensorOperator class.

[CheMPS2/include/chemps2/TensorQ.h](CheMPS2/include/chemps2/TensorQ.h) contains the definitions of the TensorQ class.

[CheMPS2/include/chemps2/TensorS0.h](CheMPS2/include/chemps2/TensorS0.h) contains the definitions of the TensorS0 class.

[CheMPS2/include/chemps2/TensorS1.h](CheMPS2/include/chemps2/TensorS1.h) contains the definitions of the TensorS1 class.

[CheMPS2/include/chemps2/TensorT.h](CheMPS2/include/chemps2/TensorT.h) contains the definitions of the TensorT class.

[CheMPS2/include/chemps2/TensorX.h](CheMPS2/include/chemps2/TensorX.h) contains the definitions of the TensorX class.

[CheMPS2/include/chemps2/ThreeDM.h](CheMPS2/include/chemps2/ThreeDM.h) contains the definitions of the ThreeDM class.

[CheMPS2/include/chemps2/TwoDM.h](CheMPS2/include/chemps2/TwoDM.h) contains the definitions of the TwoDM class.

[CheMPS2/include/chemps2/TwoIndex.h](CheMPS2/include/chemps2/TwoIndex.h) contains the definitions of the TwoIndex class.

[CheMPS2/include/chemps2/Wigner.h](CheMPS2/include/chemps2/Wigner.h) contains the definitions of the Wigner class.

Please note that these files are documented with doxygen comments. The
[doxygen html output](http://sebwouters.github.io/CheMPS2/doxygen/index.html)
can be consulted online.


List of files to perform test runs
----------------------------------

[tests/test1.cpp.in](tests/test1.cpp.in) contains several DMRG ground
state calculations in different symmetry sectors for the N2 molecule (d2h
symmetry) in the minimal STO-3G basis set.

[tests/test2.cpp.in](tests/test2.cpp.in) contains a ground state DMRG
calculation of the ^1A1 state of H2O (c2v symmetry) in the 6-31G basis set.

[tests/test3.cpp.in](tests/test3.cpp.in) contains a ground state DMRG
calculation of the ^1A1 state of CH4 (c2v symmetry) in the STO-3G basis set.

[tests/test4.cpp.in](tests/test4.cpp.in) contains a ground state DMRG
calculation of the ^6A state of a linear Hubbard chain (forced c1 symmetry)
with 10 sites and open boundary conditions, containing 9 fermions (just below
half-filling).

[tests/test5.cpp.in](tests/test5.cpp.in) contains an excited state DMRG
calculation in the ^1Ag symmetry sector of N2 (d2h symmetry) in the minimal
STO-3G basis set. The ground and two lowest excited states are determined.

[tests/test6.cpp.in](tests/test6.cpp.in) contains a state-averaged
DMRG-SCF calculation of the first excited state of the ^1Ag sector of O2 (d2h
symmetry) in the CC-pVDZ basis set. The two 1s core orbitals are kept doubly
occupied, and two Ag, B2g, B3g, B1u, B2u, and B3u orbitals are chosen as
active space. A significant speedup is obtained with DIIS.

[tests/test7.cpp.in](tests/test7.cpp.in) reads in
[tests/matrixelements/O2.CCPVDZ.FCIDUMP](tests/matrixelements/O2.CCPVDZ.FCIDUMP),
stores these matrix elements to disk, reads them back in from disk, and
compares the two versions.

[tests/test8.cpp.in](tests/test8.cpp.in) contains a DMRG-SCF ground state
calculation of the ^1Ag state of N2 (d2h symmetry) in the CC-pVDZ basis set.
The two 1s core orbitals are kept doubly occupied. The next two Ag and B1u
orbitals (sigma bonding and antibonding), as well as one B2g, B3g, B2u, and
B3u orbital (pi bonding and antibonding) are chosen as active space. A
significant speedup is obtained with DIIS. This test is smaller than test6,
and is included for debugging with valgrind.

[tests/test9.cpp.in](tests/test9.cpp.in) contains a ground state DMRG
calculation of a half-filled square 3 by 3 Hubbard lattice, both in the site
basis and in the momentum basis. For the latter, the matrix elements only have
fourfold permutation symmetry.

[tests/test10.cpp.in](tests/test10.cpp.in) is a copy of
[tests/test3.cpp.in](tests/test3.cpp.in), in which the FCI and DMRG
2- and 3-RDM are compared. This test also shows that after calculating the
2- and/or 3-RDM, it is possible to continue sweeping.

[tests/test11.cpp.in](tests/test11.cpp.in) is a copy of
[tests/test4.cpp.in](tests/test4.cpp.in), in which the FCI and DMRG
2- and 3-RDM are compared for a wavefunction with higher multiplicity.

[tests/test12.cpp.in](tests/test12.cpp.in) contains a ground state DMRG
calculation of a BCS Hamiltonian. The matrix elements only have fourfold
permutation symmetry.

[tests/test13.cpp.in](tests/test13.cpp.in) is a copy of the CASSCF
calculation in [tests/test8.cpp.in](tests/test8.cpp.in), but with full
configuration interaction (FCI) as active space solver. In addition,
the CASPT2 variational second order perturbation correction energy
is calculated.

[tests/test14.cpp.in](tests/test14.cpp.in) is a copy of the CASSCF
calculation in [tests/test8.cpp.in](tests/test8.cpp.in) with a slightly
larger active space, and which works with ordered localized orbitals instead
of natural orbitals. The localization occurs by means of Edmiston-Ruedenberg,
and the ordering based on the Fiedler vector of the exchange matrix.
In addition a calculation of the CASPT2 variational second order
perturbation correction energy in the localized (i.e. not pseudocanonical)
basis is performed.

[tests/matrixelements/CH4.STO3G.FCIDUMP](tests/matrixelements/CH4.STO3G.FCIDUMP)
contains the matrix elements for test3 and test10.

[tests/matrixelements/H2O.631G.FCIDUMP](tests/matrixelements/H2O.631G.FCIDUMP)
contains the matrix elements for test2.

[tests/matrixelements/N2.STO3G.FCIDUMP](tests/matrixelements/N2.STO3G.FCIDUMP)
contains the matrix elements for test1 and test5.

[tests/matrixelements/O2.CCPVDZ.FCIDUMP](tests/matrixelements/O2.CCPVDZ.FCIDUMP)
contains the matrix elements for test6 and test7.

[tests/matrixelements/N2.CCPVDZ.FCIDUMP](tests/matrixelements/N2.CCPVDZ.FCIDUMP)
contains the matrix elements for test8, test13, and test14.

The python tests in [PyCheMPS2/tests/](PyCheMPS2/tests/) are an identical
conversion of the c++ tests.

These test files illustrate how to use libchemps2. Note that the
tests are too small to see (near) linear scaling with the number of cores,
although improvement should still be noticeable.

