![CheMPS2 Logo](../../blob/master/CheMPS2/CheMPS2logo.png?raw=true)

CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
==============================================================================

Copyright (C) 2013 Sebastian Wouters <sebastianwouters@gmail.com>

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
   Submitted to Computer Physics Communications,
   [arXiv:1312.2415](http://arxiv.org/abs/1312.2415)

        @article{CheMPS2cite1,
            author = {Sebastian Wouters and Ward Poelmans and Paul W. Ayers and Dimitri {Van Neck}},
            title = {CheMPS2: a free open-source spin-adapted implementation of the density matrix renormalization group for ab initio quantum chemistry},
            journal = {ArXiv e-prints},
            archivePrefix = "arXiv",
            eprint = {1312.2415},
            primaryClass = "cond-mat.str-el",
            year = {2013}
        }

2. S. Wouters, P.A. Limacher, D. Van Neck and P.W. Ayers,
   Longitudinal static optical properties of hydrogen chains: Finite field
                extrapolations of matrix product state calculations,
   The Journal of Chemical Physics 136, 134110 (2012),
   [doi:10.1063/1.3700087](http://dx.doi.org/10.1063/1.3700087)

        @article{CheMPS2cite2,
            author = {Sebastian Wouters and Peter A. Limacher and Dimitri {Van Neck} and Paul W. Ayers},
            title = {Longitudinal static optical properties of hydrogen chains: Finite field extrapolations of matrix product state calculations},
            journal = {The Journal of Chemical Physics},
            year = {2012},
            volume = {136},
            number = {13}, 
            pages = {134110},
            doi = {10.1063/1.3700087} 
        }


List of files in the CheMPS2 library
------------------------------------

    CheMPS2/CASSCF.cpp
    CheMPS2/CASSCFdebug.cpp
    CheMPS2/CASSCFhamiltonianrotation.cpp
    CheMPS2/CASSCFnewtonraphson.cpp
    CheMPS2/ConvergenceScheme.cpp
    CheMPS2/DMRG.cpp
    CheMPS2/DMRGmpsio.cpp
    CheMPS2/DMRGoperators.cpp
    CheMPS2/DMRGtechnics.cpp
    CheMPS2/FourIndex.cpp
    CheMPS2/Hamiltonian.cpp
    CheMPS2/Heff.cpp
    CheMPS2/HeffDiagonal.cpp
    CheMPS2/HeffDiagrams1.cpp
    CheMPS2/HeffDiagrams2.cpp
    CheMPS2/HeffDiagrams3.cpp
    CheMPS2/HeffDiagrams4.cpp
    CheMPS2/HeffDiagrams5.cpp
    CheMPS2/Irreps.cpp
    CheMPS2/PrintLicense.cpp
    CheMPS2/Problem.cpp
    CheMPS2/Sobject.cpp
    CheMPS2/SyBookkeeper.cpp
    CheMPS2/TensorA.cpp
    CheMPS2/TensorB.cpp
    CheMPS2/TensorC.cpp
    CheMPS2/TensorD.cpp
    CheMPS2/TensorDiag.cpp
    CheMPS2/TensorF0Cbase.cpp
    CheMPS2/TensorF0.cpp
    CheMPS2/TensorF1.cpp
    CheMPS2/TensorF1Dbase.cpp
    CheMPS2/TensorL.cpp
    CheMPS2/TensorO.cpp
    CheMPS2/TensorQ.cpp
    CheMPS2/TensorS0Abase.cpp
    CheMPS2/TensorS0.cpp
    CheMPS2/TensorS1Bbase.cpp
    CheMPS2/TensorS1.cpp
    CheMPS2/TensorSwap.cpp
    CheMPS2/TensorT.cpp
    CheMPS2/TensorX.cpp
    CheMPS2/TwoDM.cpp
    CheMPS2/TwoIndex.cpp
    CheMPS2/include/CASSCF.h
    CheMPS2/include/ConvergenceScheme.h
    CheMPS2/include/DMRG.h
    CheMPS2/include/FourIndex.h
    CheMPS2/include/Gsl.h
    CheMPS2/include/Hamiltonian.h
    CheMPS2/include/Heff.h
    CheMPS2/include/Irreps.h
    CheMPS2/include/Lapack.h
    CheMPS2/include/Options.h
    CheMPS2/include/Problem.h
    CheMPS2/include/Sobject.h
    CheMPS2/include/SyBookkeeper.h
    CheMPS2/include/TensorA.h
    CheMPS2/include/TensorB.h
    CheMPS2/include/TensorC.h
    CheMPS2/include/TensorD.h
    CheMPS2/include/TensorDiag.h
    CheMPS2/include/TensorF0Cbase.h
    CheMPS2/include/TensorF0.h
    CheMPS2/include/TensorF1Dbase.h
    CheMPS2/include/TensorF1.h
    CheMPS2/include/Tensor.h
    CheMPS2/include/TensorL.h
    CheMPS2/include/TensorO.h
    CheMPS2/include/TensorQ.h
    CheMPS2/include/TensorS0Abase.h
    CheMPS2/include/TensorS0.h
    CheMPS2/include/TensorS1Bbase.h
    CheMPS2/include/TensorS1.h
    CheMPS2/include/TensorSwap.h
    CheMPS2/include/TensorT.h
    CheMPS2/include/TensorX.h
    CheMPS2/include/TwoDM.h
    CheMPS2/include/TwoIndex.h

Please note that these files are documented with Doxygen-readable comments.
Search for the section "Build" in README to see how a manual can be
generated from these comments.


List of files to perform test runs
----------------------------------

    tests/test1.cpp
    tests/test2.cpp
    tests/test3.cpp
    tests/test4.cpp
    tests/test5.cpp
    tests/test6.cpp
    tests/matrixelements/CH4_N10_S0_c2v_I0.dat
    tests/matrixelements/H6_N6_S0_d2h_I0.dat
    tests/matrixelements/N2_N14_S0_d2h_I0.dat
    tests/matrixelements/O2_CCPVDZ.dat

These test files illustrate how to use the CheMPS2 library.
They only require a very limited amount of memory (order 10-100 MB).


Matrix elements from Psi4
-------------------------

CheMPS2 has a Hamiltonian object which is able to read in matrix elements
from a plugin to Psi4
[Psi4, Ab initio quantum chemistry](http://www.psicode.org),
which works on version psi4.0b5 and higher.

To make use of this feature, build Psi4 with the plugin option, and then run:

    > psi4 --new-plugin mointegrals
    > cd mointegrals

Now, replace the file ```mointegrals.cc``` with either:

1. ```./mointegrals/mointegrals.cc_PRINT``` to print the matrix elements as
   text. Examples of output generated with this plugin can be found in
   ```./tests/matrixelements```

2. ```./mointegrals/mointegrals.cc_SAVEHAM``` to store all unique matrix
   elements (remember that there is eightfold permutation symmetry) in
   binary form with HDF5. See the Doxygen manual for more information on
   CheMPS2::Hamiltonian.

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


Build
-----

### 1. Build CheMPS2 with CMake

CheMPS2 can be build with CMake. The files

    ./CMakeLists.txt
    ./CheMPS2/CMakeLists.txt
    ./tests/CMakeLists.txt

provide a minimal compilation. Start in ```./``` and run:

    > mkdir build
    > cd build
    
CMake generates makefiles based on the user's specifications:

    > CXX=option1 cmake .. -DMKL=option2 -DBUILD_DOCUMENTATION=option3
    
Option1 is the c++ compiler; typically ```g++``` or ```icpc``` on Linux.
Option2 can be ```ON``` or ```OFF``` and is used to switch on the
intel math kernel library.
Option3 can be ```ON``` or ```OFF``` and is used to switch on doxygen
documentation.

To compile, run:

    > make


### 2. Testing CheMPS2

To test CheMPS2, start in ```./build```, and run:

    > cd tests/
    > ./test1
    > ./test2
    > ./test3
    > ./test4
    > ./test5
    > ./test6

The tests should end with a line stating whether or not they succeeded.
They only require a very limited amount of memory (order 10-100 MB).

### 3. Doxygen documentation

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

Doxygen output can be generated, see the section "Build" in README, or
the [Doxygen html output](http://sebwouters.github.io/CheMPS2/index.html).

For a more concise overview of the concepts and ideas used in CheMPS2,
please read the code release paper "CheMPS2: a free open-source spin-adapted
implementation of the density matrix renormalization group for ab initio
quantum chemistry", one of the references in CITATIONS.

