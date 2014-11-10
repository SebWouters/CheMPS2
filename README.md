![CheMPS2 Logo](https://raw.githubusercontent.com/sebwouters/CheMPS2/master/CheMPS2/CheMPS2logo.png?raw=true)

CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
==============================================================================

Copyright (C) 2013, 2014 Sebastian Wouters <sebastianwouters@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


How to acknowledge CheMPS2
--------------------------

To acknowledge CheMPS2, please cite

1. S. Wouters, W. Poelmans, P.W. Ayers and D. Van Neck,
   CheMPS2: a free open-source spin-adapted implementation of the density
            matrix renormalization group for ab initio quantum chemistry,
   Computer Physics Communications 185, 1501 (2014),
   [doi:10.1016/j.cpc.2014.01.019](http://dx.doi.org/10.1016/j.cpc.2014.01.019)

        @article{CheMPS2cite1,
            author = {Sebastian Wouters and Ward Poelmans and Paul W. Ayers and Dimitri {Van Neck}},
            title = {CheMPS2: a free open-source spin-adapted implementation of the density matrix renormalization group for ab initio quantum chemistry},
            journal = {Computer Physics Communications},
            year = {2014},
            volume = {185},
            number = {6},
            pages = {1501},
            doi = {10.1016/j.cpc.2014.01.019}
        }

2. S. Wouters and D. Van Neck,
   The density matrix renormalization group for ab initio quantum chemistry,
   The European Physical Journal D 68, 272 (2014),
   [doi:10.1140/epjd/e2014-50500-1](http://dx.doi.org/10.1140/epjd/e2014-50500-1)

        @article{CheMPS2cite2,
            author = {Sebastian Wouters and Dimitri {Van Neck}},
            title = {The density matrix renormalization group for ab initio quantum chemistry},
            journal = {The European Physical Journal D},
            year = {2014},
            volume = {68},
            number = {9},
            pages = {272},
            doi = {10.1140/epjd/e2014-50500-1}
        }

Build
-----

### 1. Build CheMPS2 with CMake

CheMPS2 can be built with CMake. The files

    ./CMakeLists.txt
    ./CheMPS2/CMakeLists.txt
    ./tests/CMakeLists.txt
    ./PyCheMPS2/CMakeLists.txt

provide a minimal compilation. Start in ```./``` and run:

    > mkdir build
    > cd build
    
CMake generates makefiles based on the user's specifications:

    > CXX=option1 cmake .. -DMKL=option2 -DBUILD_DOCUMENTATION=option3
    -DCMAKE_INSTALL_PREFIX=option4
    
Option1 is the c++ compiler; typically ```g++``` or ```icpc``` on Linux.
Option2 can be ```ON``` or ```OFF``` and is used to switch on the
intel math kernel library.
Option3 can be ```ON``` or ```OFF``` and is used to switch on doxygen
documentation.
Option4 is the prefix of the installation directory; typically
```/usr``` on Linux. The CheMPS2 library is then installed in
```prefix/lib``` and the headers in ```prefix/include/CheMPS2```.
If one or more of the required libraries are not found, please use the
command

    > CMAKE_INCLUDE_PATH=option5 CMAKE_LIBRARY_PATH=option6 CXX=option1
    cmake .. -DMKL=option2 -DBUILD_DOCUMENTATION=option3
    -DCMAKE_INSTALL_PREFIX=option4
    
instead, where option5 and option6 are respectively the missing
colon-separated include and library paths, e.g.

    CMAKE_INCLUDE_PATH=/my_libs/lib1/include/:/my_libs/lib2/include/

and

    CMAKE_LIBRARY_PATH=/my_libs/lib1/lib/:/my_libs/lib2/lib/

To compile, run:

    > make

To install, run:

    > make install

For non-standard installation directories, please remember to append
the library path to ```LD_LIBRARY_PATH``` in your .bashrc.


### 2. Testing CheMPS2

To test CheMPS2, start in ```./build```, and run:

    > cd tests/
    > ./test1
    > ./test2
    > ./test3
    > ./test4
    > ./test5
    > ./test6
    > ./test7
    > ./test8
    > ./test9

These cpp tests should end with a line stating whether or not they succeeded.
They only require a very limited amount of memory (order 10-120 MB).


### 3. Build PyCheMPS2 with Cython

PyCheMPS2, a python inferface to CheMPS2, can be built with Cython. The
installation above generated the file ```./build/setup.py```, in which
the CheMPS2 library and include paths have been assigned. This section
assumes that CheMPS2 has been installed as described above! Start in
```./build```, and run:

    > python setup.py build_ext -i


### 4. Testing PyCheMPS2

To test PyCheMPS2, start in ```./build```, and run:

    > cd PyTests/
    > python test1.py
    > python test2.py
    > python test3.py
    > python test4.py
    > python test5.py
    > python test6.py
    > python test7.py
    > python test8.py
    > python test9.py

These tests do exactly the same thing as the cpp tests above, and illustrate
the usage of the python interface to CheMPS2. The tests should end with a
line stating whether or not they succeeded. They only require a very limited
amount of memory (order 10-120 MB).


### 5. Doxygen documentation

To build and view the Doxygen manual, the documentation flag should have
been on: ```-DBUILD_DOCUMENTATION=ON```. Start in ```./build``` and run:

    > make doc
    > cd LaTeX-documents
    > make
    > evince refman.pdf &
    > cd ../html
    > firefox index.html &


User manual
-----------

Doxygen output can be generated, see the section "Build" in README.md, or
the [Doxygen html output](http://sebwouters.github.io/CheMPS2/index.html).

For a more concise overview of the concepts and ideas used in CheMPS2,
please read the code release paper "[CheMPS2: a free open-source spin-adapted
implementation of the density matrix renormalization group for ab initio
quantum chemistry](http://doi.dx.org/10.1016/j.cpc.2014.01.019)", one of the
references in CITATIONS.


Matrix elements from Psi4
-------------------------

CheMPS2 has a Hamiltonian object which is able to read in matrix elements
from a plugin to [Psi4, Ab initio quantum chemistry](http://www.psicode.org),
which works on version psi4.0b5 and higher.

To make use of this feature, build Psi4 with the plugin option, and then run:

    > psi4 --new-plugin mointegrals
    > cd mointegrals

Now, replace the file ```mointegrals.cc``` with either:

1. [```mointegrals/mointegrals.cc_PRINT```](mointegrals/mointegrals.cc_PRINT)
   to print the matrix elements as text. Examples of output generated with this
   plugin can be found in [```tests/matrixelements```](tests/matrixelements)

2. [```mointegrals/mointegrals.cc_SAVEHAM```](mointegrals/mointegrals.cc_SAVEHAM)
   to store all unique matrix elements (remember that there is eightfold
   permutation symmetry) in binary form with HDF5. See the Doxygen manual for
   more information on CheMPS2::Hamiltonian.

For case 2, the ```Makefile``` should be adjusted. Replace
    ```$(PSILIBS)```
with
    ```$(PSILIBS) -L${CheMPS2_BINARY_DIR}/CheMPS2/ -lCheMPS2```
and replace
    ```$(CXXINCLUDE)```
with
    ```$(CXXINCLUDE) -I${CheMPS2_SOURCE_DIR}/CheMPS2/include/```

To compile the Psi4 plugin, run:

    > make


List of files in the CheMPS2 library
------------------------------------

[```CheMPS2/CASSCF.cpp```](CheMPS2/CASSCF.cpp) contains the constructor and
destructor of the CASSCF class, as well as the functions which allow to
update, get or set its variables.
    
[```CheMPS2/CASSCFdebug.cpp```](CheMPS2/CASSCFdebug.cpp) contains two
functions: one for calculating the ROHF energy; and one for fetching FCI
coefficients to determine the point group symmetry of certain electronic
states of the carbon dimer.
    
[```CheMPS2/CASSCFnewtonraphson.cpp```](CheMPS2/CASSCFnewtonraphson.cpp)
contains all DMRG-SCF functions which are specific to the augmented Hessian
Newton-Raphson update scheme, including the functions to calculate the
gradient and Hessian.
    
[```CheMPS2/ConvergenceScheme.cpp```](CheMPS2/ConvergenceScheme.cpp) contains
all functions of the ConvergenceScheme class, which contains the instructions
for the subsequent DMRG sweeps.

[```CheMPS2/Correlations.cpp```](CheMPS2/Correlations.cpp) contains all the
functionality to calculate the spin, density, and spin-flip correlation
functions as well as the two-orbital mutual information.

[```CheMPS2/DIIS.cpp```](CheMPS2/DIIS.cpp) contains a DIIS convergence
speed-up for DMRG-SCF.

[```CheMPS2/DMRG.cpp```](CheMPS2/DMRG.cpp) contains the constructor and
destructor of the DMRG class, as well as the top-level sweep functions.

[```CheMPS2/DMRGmpsio.cpp```](CheMPS2/DMRGmpsio.cpp) contains the store and
load functions for the DMRG checkpoint file (the MPS).

[```CheMPS2/DMRGoperators.cpp```](CheMPS2/DMRGoperators.cpp) contains all
functions related to the DMRG renormalized operators: saving to disk, loading
from disk, and updating.

[```CheMPS2/DMRGSCFindices.cpp```](CheMPS2/DMRGSCFindices.cpp) contains the
index conversions for the DMRG-SCF algorithm.

[```CheMPS2/DMRGSCFoptions.cpp```](CheMPS2/DMRGSCFoptions.cpp) is a container
class to pass the DMRGSCF options to the augmented Hessian Newton-Raphson
routine.

[```CheMPS2/DMRGSCFunitary.cpp```](CheMPS2/DMRGSCFunitary.cpp) contains the
storage and handling of the unitary matrix and its nonredundant skew-symmetric
parametrization required for the DMRG-SCF algorithm.

[```CheMPS2/DMRGSCFVmatRotations.cpp```](CheMPS2/DMRGSCFVmatRotations.cpp)
performs the two-body matrix element rotations for the DMRGSCF and
Edmiston-Ruedenberg classes.

[```CheMPS2/DMRGtechnics.cpp```](CheMPS2/DMRGtechnics.cpp) contains the
functions related to the 2-RDM and the excited-state calculations.

[```CheMPS2/EdmistonRuedenberg.cpp```](CheMPS2/EdmistonRuedenberg.cpp) contains
an orbital localization function based on the Edmiston-Ruedenberg cost function
and an augmented Hessian Newton-Raphson optimizer.

[```CheMPS2/FCI.cpp```](CheMPS2/FCI.cpp) contains a full configuration
interaction solver based on Davidson's algorithm.

[```CheMPS2/FourIndex.cpp```](CheMPS2/FourIndex.cpp) contains all functions of
the FourIndex container class for the two-body matrix elements.

[```CheMPS2/Hamiltonian.cpp```](CheMPS2/Hamiltonian.cpp) contains all functions
of the Hamiltonian class, including functions to get or set specific variables,
as well as the save and load functions to store the Hamiltonian on disk.

[```CheMPS2/Heff.cpp```](CheMPS2/Heff.cpp) contains the constructor and
destructor of the Heff class, as well as the top-level functions to perform
the effective Hamiltonian times guess-vector multiplication and Davidson's
algorithm.

[```CheMPS2/HeffDiagonal.cpp```](CheMPS2/HeffDiagonal.cpp) contains the
functions to calculate the diagonal elements of the effective Hamiltonian.
These are required as preconditioner in Davidson's algorithm.
    
[```CheMPS2/HeffDiagrams1.cpp```](CheMPS2/HeffDiagrams1.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess-vector
multiplication.

[```CheMPS2/HeffDiagrams2.cpp```](CheMPS2/HeffDiagrams2.cpp) contains a subset
 of functions to perform the effective Hamiltonian times guess-vector
 multiplication.

[```CheMPS2/HeffDiagrams3.cpp```](CheMPS2/HeffDiagrams3.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess-vector
multiplication.

[```CheMPS2/HeffDiagrams4.cpp```](CheMPS2/HeffDiagrams4.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess-vector
multiplication.

[```CheMPS2/HeffDiagrams5.cpp```](CheMPS2/HeffDiagrams5.cpp) contains a subset
of functions to perform the effective Hamiltonian times guess-vector
multiplication.

[```CheMPS2/Initialize.cpp```](CheMPS2/Initialize.cpp) allows to set the seed
of the random number generator and cout.precision (added for PyCheMPS2).

[```CheMPS2/Irreps.cpp```](CheMPS2/Irreps.cpp) contains the Psi4 symmetry
conventions.

[```CheMPS2/PrintLicense.cpp```](CheMPS2/PrintLicense.cpp) contains a function
which prints the license disclaimer.

[```CheMPS2/Problem.cpp```](CheMPS2/Problem.cpp) contains all Problem class
functions. This wrapper class allows to set the desired symmetry sector for
the DMRG algorithm.

[```CheMPS2/Sobject.cpp```](CheMPS2/Sobject.cpp) contains all Sobject class
functions. This class constructs, stores, and decomposes the reduced two-site
object.

[```CheMPS2/SyBookkeeper.cpp```](CheMPS2/SyBookkeeper.cpp) contains all
SyBookkeeper functions. This class keeps track of the FCI and DMRG virtual
dimensions of all symmetry sectors at all boundaries.

[```CheMPS2/TensorA.cpp```](CheMPS2/TensorA.cpp) contains the TensorA
functions. This class stores and handles the complementary reduced spin-0
S0-tensors.

[```CheMPS2/TensorB.cpp```](CheMPS2/TensorB.cpp) contains the TensorB
functions. This class stores and handles the complementary reduced spin-1
S1-tensors.

[```CheMPS2/TensorC.cpp```](CheMPS2/TensorC.cpp) contains the TensorC
functions. This class stores and handles the complementary reduced spin-0
F0-tensors.

[```CheMPS2/TensorD.cpp```](CheMPS2/TensorD.cpp) contains the TensorD
functions. This class stores and handles the complementary reduced spin-1
F1-tensors.

[```CheMPS2/TensorDiag.cpp```](CheMPS2/TensorDiag.cpp) contains the TensorDiag functions. This storage
class handles reduced tensors which are diagonal in the symmetry sectors.

[```CheMPS2/TensorF0Cbase.cpp```](CheMPS2/TensorF0Cbase.cpp) contains the common storage and handling
functions of the TensorF0 and TensorC classes.

[```CheMPS2/TensorF0.cpp```](CheMPS2/TensorF0.cpp) contains all TensorF0 functions. This class
stores and handles the reduced spin-0 part of two sandwiched second
quantized operators, of which the particle symmetry sectors are equal.

[```CheMPS2/TensorF1.cpp```](CheMPS2/TensorF1.cpp) contains all TensorF1 functions. This class
stores and handles the reduced spin-1 part of two sandwiched second
quantized operators, of which the particle symmetry sectors are equal.

[```CheMPS2/TensorF1Dbase.cpp```](CheMPS2/TensorF1Dbase.cpp) contains the common storage and handling
functions of the TensorF1 and TensorD classes.

[```CheMPS2/TensorGYZ.cpp```](CheMPS2/TensorGYZ.cpp) contains the contruct and update functions
for the G-, Y-, and Z-tensors. They are required for the two-orbital mutual
information.

[```CheMPS2/TensorK.cpp```](CheMPS2/TensorK.cpp) contains the contruct and update functions
for the K-tensor. It is required for the two-orbital mutual information.

[```CheMPS2/TensorL.cpp```](CheMPS2/TensorL.cpp) contains all TensorL functions. This class stores
and handles the reduced spin-1/2 part of a single sandwiched second quantized
operator.

[```CheMPS2/TensorM.cpp```](CheMPS2/TensorM.cpp) contains the contruct and update functions
for the M-tensor. It is required for the two-orbital mutual information.

[```CheMPS2/TensorO.cpp```](CheMPS2/TensorO.cpp) implements the storage and handling of the partial
terms which are required to calculate the overlap between two MPSs.

[```CheMPS2/TensorQ.cpp```](CheMPS2/TensorQ.cpp) contains all TensorQ functions. This class stores
and handles the complementary reduced spin-1/2 part of three sandwiched second
quantized operators.

[```CheMPS2/TensorS0Abase.cpp```](CheMPS2/TensorS0Abase.cpp) contains the common storage and handling
functions of the TensorS0 and TensorA classes.

[```CheMPS2/TensorS0.cpp```](CheMPS2/TensorS0.cpp) contains all TensorS0 functions. This class
stores and handles the reduced spin-0 part of two sandwiched second
quantized operators, of which the particle symmetry sectors differ by 2.

[```CheMPS2/TensorS1Bbase.cpp```](CheMPS2/TensorS1Bbase.cpp) contains the common storage and handling
functions of the TensorS1 and TensorB classes.

[```CheMPS2/TensorS1.cpp```](CheMPS2/TensorS1.cpp) contains all TensorS1 functions. This class
stores and handles the reduced spin-1 part of two sandwiched second
quantized operators, of which the particle symmetry sectors differ by 2.

[```CheMPS2/TensorSwap.cpp```](CheMPS2/TensorSwap.cpp) contains the common storage and handling
functions of the TensorL and TensorQ classes.

[```CheMPS2/TensorT.cpp```](CheMPS2/TensorT.cpp) contains all TensorT functions. This class
stores and handles the reduced part of an MPS site-tensor.

[```CheMPS2/TensorX.cpp```](CheMPS2/TensorX.cpp) contains all TensorX functions. This class
stores and handles the complementary reduced spin-0 part of four sandwiched
second quantized operators, which is of course diagonal in the symmetry
sectors.

[```CheMPS2/TwoDM.cpp```](CheMPS2/TwoDM.cpp) contains all functions to calculate the
2-RDM from the DMRG-optimized MPS.

[```CheMPS2/TwoDMstorage.cpp```](CheMPS2/TwoDMstorage.cpp) contains all functions to store the
2-RDM from the DMRG-optimized MPS.

[```CheMPS2/TwoIndex.cpp```](CheMPS2/TwoIndex.cpp) contains all functions of the TwoIndex container
class for the one-body matrix elements.
    
[```CheMPS2/include/CASSCF.h```](CheMPS2/include/CASSCF.h) contains the definitions of the CASSCF class.

[```CheMPS2/include/ConvergenceScheme.h```](CheMPS2/include/ConvergenceScheme.h) contains the definitions of the
ConvergenceScheme class.

[```CheMPS2/include/Correlations.h```](CheMPS2/include/Correlations.h) contains the definitions of the
Correlations class.

[```CheMPS2/include/DIIS.h```](CheMPS2/include/DIIS.h) contains the definitions of the DIIS class.

[```CheMPS2/include/DMRG.h```](CheMPS2/include/DMRG.h) contains the definitions of the DMRG class.

[```CheMPS2/include/DMRGSCFindices.h```](CheMPS2/include/DMRGSCFindices.h) contains the definitions of the
DMRGSCFindices class.

[```CheMPS2/include/DMRGSCFoptions.h```](CheMPS2/include/DMRGSCFoptions.h) contains the definitions of the
DMRGSCFoptions container class.

[```CheMPS2/include/DMRGSCFunitary.h```](CheMPS2/include/DMRGSCFunitary.h) contains the definitions of the
DMRGSCFunitary class.

[```CheMPS2/include/DMRGSCFVmatRotations.h```](CheMPS2/include/DMRGSCFVmatRotations.h) contains the definitions of the
DMRGSCFVmatRotations class.

[```CheMPS2/include/EdmistonRuedenberg.h```](CheMPS2/include/EdmistonRuedenberg.h) contains the definitions of the
EdmistonRuedenberg class.

[```CheMPS2/include/FCI.h```](CheMPS2/include/FCI.h) contains the definitions of the FCI class.

[```CheMPS2/include/FourIndex.h```](CheMPS2/include/FourIndex.h) contains the definitions of the FourIndex
class.

[```CheMPS2/include/Gsl.h```](CheMPS2/include/Gsl.h) contains the definitions of the external GSL
routines.

[```CheMPS2/include/Hamiltonian.h```](CheMPS2/include/Hamiltonian.h) contains the definitions of the Hamiltonian
class.

[```CheMPS2/include/Heff.h```](CheMPS2/include/Heff.h) contains the definitions of the Heff class.

[```CheMPS2/include/Initialize.h```](CheMPS2/include/Initialize.h) contains the definitions of the Initialize
class.

[```CheMPS2/include/Irreps.h```](CheMPS2/include/Irreps.h) contains the definitions of the Irrep class.

[```CheMPS2/include/Lapack.h```](CheMPS2/include/Lapack.h) contains the definitions of the external BLAS
and LAPACK routines.

[```CheMPS2/include/MyHDF5.h```](CheMPS2/include/MyHDF5.h) forces the use of the HDF5 1.8 API, e.g. 
H5Gcreate2 instead of H5Gcreate1, a known issue in Ubuntu 12.04.

[```CheMPS2/include/Options.h```](CheMPS2/include/Options.h) contains all the options of the CheMPS2
namespace. Here the checkpoint storage names and folders can be set, as well
as parameters related to memory usage and convergence.

[```CheMPS2/include/Problem.h```](CheMPS2/include/Problem.h) contains the definitions of the Problem class.

[```CheMPS2/include/Sobject.h```](CheMPS2/include/Sobject.h) contains the definitions of the Sobject class.

[```CheMPS2/include/SyBookkeeper.h```](CheMPS2/include/SyBookkeeper.h) contains the definitions of the
SyBookkeeper class.

[```CheMPS2/include/TensorA.h```](CheMPS2/include/TensorA.h) contains the definitions of the TensorA class.

[```CheMPS2/include/TensorB.h```](CheMPS2/include/TensorB.h) contains the definitions of the TensorB class.

[```CheMPS2/include/TensorC.h```](CheMPS2/include/TensorC.h) contains the definitions of the TensorC class.

[```CheMPS2/include/TensorD.h```](CheMPS2/include/TensorD.h) contains the definitions of the TensorD class.

[```CheMPS2/include/TensorDiag.h```](CheMPS2/include/TensorDiag.h) contains the definitions of the
TensorDiag class.

[```CheMPS2/include/TensorF0Cbase.h```](CheMPS2/include/TensorF0Cbase.h) contains the definitions of the
TensorF0Cbase class.

[```CheMPS2/include/TensorF0.h```](CheMPS2/include/TensorF0.h) contains the definitions of the TensorF0 class.

[```CheMPS2/include/TensorF1Dbase.h```](CheMPS2/include/TensorF1Dbase.h) contains the definitions of the
TensorF1Dbase class.

[```CheMPS2/include/TensorF1.h```](CheMPS2/include/TensorF1.h) contains the definitions of the TensorF1 class.

[```CheMPS2/include/TensorGYZ.h```](CheMPS2/include/TensorGYZ.h) contains the definitions of the TensorGYZ
class.

[```CheMPS2/include/Tensor.h```](CheMPS2/include/Tensor.h) contains the definitions of the virtual Tensor
class.

[```CheMPS2/include/TensorK.h```](CheMPS2/include/TensorK.h) contains the definitions of the TensorK class.

[```CheMPS2/include/TensorL.h```](CheMPS2/include/TensorL.h) contains the definitions of the TensorL class.

[```CheMPS2/include/TensorM.h```](CheMPS2/include/TensorM.h) contains the definitions of the TensorM class.

[```CheMPS2/include/TensorO.h```](CheMPS2/include/TensorO.h) contains the definitions of the TensorO class.

[```CheMPS2/include/TensorQ.h```](CheMPS2/include/TensorQ.h) contains the definitions of the TensorQ class.

[```CheMPS2/include/TensorS0Abase.h```](CheMPS2/include/TensorS0Abase.h) contains the definitions of the
TensorS0Abase class.

[```CheMPS2/include/TensorS0.h```](CheMPS2/include/TensorS0.h) contains the definitions of the TensorS0 class.

[```CheMPS2/include/TensorS1Bbase.h```](CheMPS2/include/TensorS1Bbase.h) contains the definitions of the
TensorS1Bbase class.

[```CheMPS2/include/TensorS1.h```](CheMPS2/include/TensorS1.h) contains the definitions of the TensorS1 class.

[```CheMPS2/include/TensorSwap.h```](CheMPS2/include/TensorSwap.h) contains the definitions of the TensorSwap
class.

[```CheMPS2/include/TensorT.h```](CheMPS2/include/TensorT.h) contains the definitions of the TensorT class.

[```CheMPS2/include/TensorX.h```](CheMPS2/include/TensorX.h) contains the definitions of the TensorX class.

[```CheMPS2/include/TwoDM.h```](CheMPS2/include/TwoDM.h) contains the definitions of the TwoDM class.

[```CheMPS2/include/TwoDMstorage.h```](CheMPS2/include/TwoDMstorage.h) contains the definitions of the
TwoDMstorage class.

[```CheMPS2/include/TwoIndex.h```](CheMPS2/include/TwoIndex.h) contains the definitions of the TwoIndex
class.

Please note that these files are documented with Doxygen-readable comments.
Search for the section "Build" in [```README.md```](README.md) to see how a
manual can be generated from these comments.


List of files to perform test runs
----------------------------------

[```tests/test1.cpp.in```](tests/test1.cpp.in) contains several DMRG ground
state calculations in different symmetry sectors for the N2 molecule (d2h
symmetry) in the minimal STO-3G basis set.

[```tests/test2.cpp.in```](tests/test2.cpp.in) contains a ground state DMRG
calculation of the ^1A1 state of H2O (c2v symmetry) in the 6-31G basis set.

[```tests/test3.cpp.in```](tests/test3.cpp.in) contains a ground state DMRG
calculation of the ^1A1 state of CH4 (c2v symmetry) in the STO-3G basis set.

[```tests/test4.cpp.in```](tests/test4.cpp.in) contains a ground state DMRG
calculation of the ^6A state of a linear Hubbard chain (forced c1 symmetry)
with 10 sites and open boundary conditions, containing 9 fermions (just below
half-filling).

[```tests/test5.cpp.in```](tests/test5.cpp.in) contains an excited state DMRG
calculation in the ^1Ag symmetry sector of N2 (d2h symmetry) in the minimal
STO-3G basis set. The ground and two lowest excited states are determined.

[```tests/test6.cpp.in```](tests/test6.cpp.in) contains a DMRG-SCF ground state
calculation of the ^1Ag state of O2 (d2h symmetry) in the CC-pVDZ basis set.
The two 1s core orbitals are kept frozen, and two Ag, B2g, B3g, B1u, B2u, and
B3u orbitals are chosen as active space. A significant speedup is obtained with
DIIS.

[```tests/test7.cpp.in```](tests/test7.cpp.in) reads in
[```tests/matrixelements/O2_CCPVDZ.dat```](tests/matrixelements/O2_CCPVDZ.dat),
stores these matrix elements to disk, reads them back in from disk, and
compares the two versions.

[```tests/test8.cpp.in```](tests/test8.cpp.in) contains a DMRG-SCF ground state
calculation of the ^1Ag state of N2 (d2h symmetry) in the CC-pVDZ basis set.
The two 1s core orbitals are kept frozen. The next two Ag and B1u orbitals
(sigma bonding and antibonding), as well as one B2g, B3g, B2u, and B3u orbital
(pi bonding and antibonding) are chosen as active space. A significant speedup
is obtained with DIIS. This test is smaller than test6, and is included for
debugging with valgrind.

[```tests/test9.cpp.in```](tests/test9.cpp.in) is a copy of
[```tests/test8.cpp.in```](tests/test8.cpp.in), with a slightly larger active
space and which works with ordered localized orbitals instead of natural
orbitals. The localization occurs by means of Edmiston-Ruedenberg, and the
ordering based on the Fiedler vector with the exchange matrix as cost function.

[```tests/matrixelements/CH4_N10_S0_c2v_I0.dat```](tests/matrixelements/CH4_N10_S0_c2v_I0.dat)
contains the matrix elements for test3.

[```tests/matrixelements/H2O_631G_eq.dat```](tests/matrixelements/H2O_631G_eq.dat)
contains the matrix elements for test2.

[```tests/matrixelements/N2_N14_S0_d2h_I0.dat```](tests/matrixelements/N2_N14_S0_d2h_I0.dat)
contains the matrix elements for test1 and test5.

[```tests/matrixelements/O2_CCPVDZ.dat```](tests/matrixelements/O2_CCPVDZ.dat)
contains the matrix elements for test6 and test7.

[```tests/matrixelements/N2_CCPVDZ.dat```](tests/matrixelements/N2_CCPVDZ.dat)
contains the matrix elements for test8 and test9.

[```PyCheMPS2/test1.py.in```](PyCheMPS2/test1.py.in) is a Python implementation
of a part of [```tests/test1.cpp.in```](tests/test1.cpp.in)

[```PyCheMPS2/test2.py.in```](PyCheMPS2/test2.py.in) is the Python version of
[```tests/test2.cpp.in```](tests/test2.cpp.in)

[```PyCheMPS2/test3.py.in```](PyCheMPS2/test3.py.in) is the Python version of
[```tests/test3.cpp.in```](tests/test3.cpp.in)

[```PyCheMPS2/test4.py.in```](PyCheMPS2/test4.py.in) is the Python version of
[```tests/test4.cpp.in```](tests/test4.cpp.in)

[```PyCheMPS2/test5.py.in```](PyCheMPS2/test5.py.in) is the Python version of
[```tests/test5.cpp.in```](tests/test5.cpp.in)

[```PyCheMPS2/test6.py.in```](PyCheMPS2/test6.py.in) is the Python version of
[```tests/test6.cpp.in```](tests/test6.cpp.in)

[```PyCheMPS2/test7.py.in```](PyCheMPS2/test7.py.in) is the Python version of
[```tests/test7.cpp.in```](tests/test7.cpp.in)

[```PyCheMPS2/test8.py.in```](PyCheMPS2/test8.py.in) is the Python version of
[```tests/test8.cpp.in```](tests/test8.cpp.in)

[```PyCheMPS2/test9.py.in```](PyCheMPS2/test9.py.in) is the Python version of
[```tests/test9.cpp.in```](tests/test9.cpp.in)

These test files illustrate how to use the CheMPS2 library. The tests only
require a very limited amount of memory (order 10-120 MB).

