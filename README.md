[![Build Status](https://travis-ci.org/SebWouters/CheMPS2.svg?branch=master)](https://travis-ci.org/SebWouters/CheMPS2)

![CheMPS2 Logo](https://raw.githubusercontent.com/sebwouters/CheMPS2/master/CheMPS2/CheMPS2logo.png?raw=true)

CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
==============================================================================

Copyright (C) 2013-2016 Sebastian Wouters <sebastianwouters@gmail.com>

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

Information
-----------

CheMPS2 is a scientific library which contains a spin-adapted implementation 
of the density matrix renormalization group (DMRG) for ab initio quantum 
chemistry. This wavefunction method allows to obtain numerical accuracy in
active spaces beyond the capabilities of full configuration interaction (FCI),
and allows to extract the 2-, 3-, and 4-particle reduced density matrices
(2-, 3- and 4-RDM) of the active space.

For general active spaces up to 40 electrons in 40 orbitals can be handled
with DMRG, and for one-dimensional active spaces up to 100 electrons in 100
orbitals. The 2-RDM of these active spaces can also be easily extracted,
while the 3- and 4-RDM are limited to about 28 orbitals.

When the active space size becomes prohibitively expensive for FCI, DMRG can
be used to replace the FCI solver in the complete active space self consistent
field (CASSCF) method and the corresponding complete active space second order
perturbation theory (CASPT2). The corresponding methods are called DMRG-SCF
and DMRG-CASPT2, respectively. For DMRG-SCF the active space 2-RDM is required,
and for DMRG-CASPT2 the active space 4-RDM.

CheMPS2 is designed for high-performance computers, with a hybrid
parallelization for mixed distributed and shared memory architectures,
realized with the Message Passing Interface ([MPI](http://www.mpi-forum.org/))
and the Open Multi-Processing ([OpenMP](http://openmp.org/wp/)) API.

The CheMPS2 library can be interfaced with quantum chemistry codes which
can handle R(O)HF calculations and molecular orbital matrix elements. This
has been done for [psi4](http://www.psicode.org) and
[pyscf](https://github.com/sunqm/pyscf), as described in the
[user manual](http://sebwouters.github.io/CheMPS2/index.html). Usage of the
library is illustrated in the [c++](#4-test-libchemps2) and
[python](#7-test-pychemps2) tests.

The CheMPS2 binary allows to perform DMRG-SCF and DMRG-CASPT2 calculations
based on input files. Molecular orbital matrix elements should then be
provided in FCIDUMP format. Usage of the binary is illustrated in the
[binary](#5-test-the-chemps2-binary) example.


User manual
-----------

Information on CheMPS2 can be found in the following documents:

1. [publications](#how-to-acknowledge-chemps2)
2. [user manual](http://sebwouters.github.io/CheMPS2/index.html)
3. [doxygen html output](http://sebwouters.github.io/CheMPS2/doxygen/index.html)

The [user manual](http://sebwouters.github.io/CheMPS2/index.html) contains
elaborate information on

* the installation of CheMPS2
* the DMRG, DMRG-SCF, and DMRG-CASPT2 algorithms
* the symmetries which are exploited in CheMPS2
* how to generate matrix elements with plugins to [psi4](http://www.psicode.org)
* how to perform DMRG, DMRG-SCF, and DMRG-CASPT2 calculations
* the interfaces to [psi4](http://www.psicode.org) and [pyscf](https://github.com/sunqm/pyscf)


How to acknowledge CheMPS2
--------------------------

To acknowledge CheMPS2, please cite

1. S. Wouters, W. Poelmans, P.W. Ayers and D. Van Neck,
   CheMPS2: a free open-source spin-adapted implementation of the density
   matrix renormalization group for ab initio quantum chemistry,
   *Computer Physics Communications* **185** (6), 1501-1514 (2014),
   [doi:10.1016/j.cpc.2014.01.019](http://dx.doi.org/10.1016/j.cpc.2014.01.019)

        @article{CheMPS2cite1,
            author = {Sebastian Wouters and Ward Poelmans and Paul W. Ayers and Dimitri {Van Neck}},
            title = {CheMPS2: a free open-source spin-adapted implementation of the density matrix renormalization group for ab initio quantum chemistry},
            journal = {Computer Physics Communications},
            year = {2014},
            volume = {185},
            number = {6},
            pages = {1501-1514},
            doi = {10.1016/j.cpc.2014.01.019}
        }

2. S. Wouters and D. Van Neck,
   The density matrix renormalization group for ab initio quantum chemistry,
   *European Physical Journal D* **68** (9), 272 (2014),
   [doi:10.1140/epjd/e2014-50500-1](http://dx.doi.org/10.1140/epjd/e2014-50500-1)

        @article{CheMPS2cite2,
            author = {Sebastian Wouters and Dimitri {Van Neck}},
            title = {The density matrix renormalization group for ab initio quantum chemistry},
            journal = {European Physical Journal D},
            year = {2014},
            volume = {68},
            number = {9},
            pages = {272},
            doi = {10.1140/epjd/e2014-50500-1}
        }

3. S. Wouters, T. Bogaerts, P. Van Der Voort, V. Van Speybroeck and D. Van Neck,
   Communication: DMRG-SCF study of the singlet, triplet, and quintet states of oxo-Mn(Salen),
   *Journal of Chemical Physics* **140** (24), 241103 (2014),
   [doi:10.1063/1.4885815](http://dx.doi.org/10.1063/1.4885815)
   
        @article{CheMPS2cite3,
            author = {Sebastian Wouters and Thomas Bogaerts and Pascal {Van Der Voort} and Veronique {Van Speybroeck} and Dimitri {Van Neck}},
            title = {Communication: DMRG-SCF study of the singlet, triplet, and quintet states of oxo-Mn(Salen)},
            journal = {Journal of Chemical Physics},
            year = {2014},
            volume = {140},
            number = {24},
            pages = {241103},
            doi = {10.1063/1.4885815}
        }


Installation
------------

### 1. Dependencies

CheMPS2 can be built with [CMake](http://www.cmake.org/) and depends on

-   BLAS
-   LAPACK
-   GSL ([GNU Scientific library](http://www.gnu.org/software/gsl/))
-   HDF5 ([Hierarchical Data Format Release 5](http://www.hdfgroup.org/HDF5/))

It is parallelized for shared memory architectures with the Open
Multi-Processing ([OpenMP](http://openmp.org/wp/)) API and for
distributed memory architectures with the Message Passing Interface
([MPI](http://www.mpi-forum.org/)). A hybrid combination of both
parallelization strategies is supported.

### 2. Download

It is advised to clone the CheMPS2 git repository from github. In your
terminal, run:

    > cd /sourcefolder
    > git clone 'https://github.com/sebwouters/chemps2'
    > cd chemps2

That way, future updates and bug fixes can be easily pulled in:

    > cd /sourcefolder/chemps2
    > git pull

### 3. Build the chemps2 library and binary

The files

    /sourcefolder/chemps2/CMakeLists.txt
    /sourcefolder/chemps2/CheMPS2/CMakeLists.txt
    /sourcefolder/chemps2/tests/CMakeLists.txt
    /sourcefolder/chemps2/sphinx/CMakeLists.txt

provide a minimal compilation. In your terminal, run:

    > cd /sourcefolder/chemps2
    > mkdir build
    > cd build

CMake generates makefiles based on the user’s specifications:

    > CXX=option1 cmake .. -DMKL=option2 -DCMAKE_INSTALL_PREFIX=/option3 -DWITH_MPI=option4

1.  Option1 is the `c++` compiler; typically `g++`, `icpc`, or `clang++`
    on Linux. It is advised to use the intel compiler if available.
2.  Option2 can be `ON` or `OFF` and is used to switch on the intel math
    kernel library.
3.  /option3 is the prefix of the installation directory; typically `/usr`
    or `/usr/local` on Linux. On my computer, libchemps2 is then installed
    in `/option3/lib/x86_64-linux-gnu`, the headers in
    `/option3/include/chemps2`, and the binary in
    `/option3/bin/chemps2`.
4.  Option4 can be `ON` or `OFF` and is used to switch on the possibility
    to compile with MPI. Please note that the compiler should then provide
    `mpi.h`. Option1 should hence be the `mpic++` compiler; typically
    `mpic++` or `mpiicpc` on Linux. It is advised to use the intel compiler
    if available.

If one or more of the required libraries are not found, use the command

    > CMAKE_INCLUDE_PATH=option5 CMAKE_LIBRARY_PATH=option6 CXX=option1 cmake .. -DMKL=option2 -DCMAKE_INSTALL_PREFIX=/option3 -DWITH_MPI=option4

instead, where option5 and option6 are respectively the missing
colon-separated include and library paths:

    CMAKE_INCLUDE_PATH=/my_libs/lib1/include:/my_libs/lib2/include
    CMAKE_LIBRARY_PATH=/my_libs/lib1/lib:/my_libs/lib2/lib
    
For operating systems based on debian, the HDF5 headers are located in
the folder `/usr/include/hdf5/serial`. If CMake complains about the HDF5
headers, try to pass it with the option
`-DHDF5_INCLUDE_DIRS=/usr/include/hdf5/serial`.

To compile, run:

    > make

To install, run:

    > make install

For non-standard installation directories, please remember to append the
library path to `LD_LIBRARY_PATH` in your `.bashrc`.

### 4. Test libchemps2

To test libchemps2 for compilation **without MPI**, run:

    > cd /sourcefolder/chemps2/build
    > make test

To test libchemps2 for compilation **with MPI**, run:
    
    > cd /sourcefolder/chemps2/build/tests
    > OMP_NUM_THREADS=YYY mpirun -np ZZZ ./test1
    > OMP_NUM_THREADS=YYY mpirun -np ZZZ ./test2
    ...

`YYY` specifies the number of threads per process and `ZZZ` the number of
processes. Note that the tests are too small to see (near) linear scaling
with the number of cores, although improvement should still be noticeable.

### 5. Test the chemps2 binary

To test the chemps2 binary for compilation **without MPI**, run:

    > man /sourcefolder/chemps2/chemps2.1
    > cd /sourcefolder/chemps2/build/CheMPS2
    > ./chemps2 --help
    > cp /sourcefolder/chemps2/tests/test14.input .
    > sed -i "s/path\/to/sourcefolder\/chemps2\/tests\/matrixelements/" test14.input
    > cat test14.input
    > ./chemps2 --file=test14.input
    
Note that when you use the CASPT2 checkpoint, and want to restart a
calculation at a later point, you should

 1. switch the option `SCF_ACTIVE_SPACE` to `I`
 2. remove the `CheMPS2_DIIS.h5` checkpoint

in order to ensure that **exactly** the same orbitals are used in the
different runs.
    
To test the chemps2 binary for compilation **with MPI**, prepend the binary with:

    > OMP_NUM_THREADS=YYY mpirun -np ZZZ ./chemps2 --file=test14.input

### 6. Build PyCheMPS2

PyCheMPS2 is a python interface to libchemps2, for
compilation **without MPI**. It can be built with
[Cython](http://cython.org/). The installation is independent of
CMake and assumes that you have installed the CheMPS2 library with
`make install`. For non-standard installation directories of CheMPS2,
please remember to append the library path to `LD_LIBRARY_PATH` in
your `.bashrc`. In addition, the include path should be appended
to `CPATH`:

    > export CPATH=${CPATH}:/option3/include

where `/option3` is the option provided to CMake with
`-DCMAKE_INSTALL_PREFIX=/option3` above. For operating systems based
on debian, the HDF5 headers are located in the folder
`/usr/include/hdf5/serial`. If it was explicitly passed to CMake, it
should also be appended to `CPATH`:

    > export CPATH=${CPATH}:/option3/include:/usr/include/hdf5/serial

The python wrapper can be installed with:

    > cd /sourcefolder/chemps2/PyCheMPS2
    > python setup.py build_ext -L ${LD_LIBRARY_PATH}
    > python setup.py install --prefix=/option3

On my machine, the python wrapper is installed to the folder
`/option3/lib/python2.7/site-packages`, but the folder `lib` and
the distribution of python can vary.

Compilation of PyCheMPS2 occurs by linking to the `c++` library in the
installation directory. The installation of PyCheMPS2 will fail if that
library is not properly installed. If you have pulled a newer version of
CheMPS2, please remember to reinstall the `c++` library first, before
reinstalling PyCheMPS2!


### 7. Test PyCheMPS2

When libchemps2 has been compiled **without MPI**, PyCheMPS2 can be tested
by running (remember that the python site-packages folder can vary):

    > cd /sourcefolder/chemps2/PyCheMPS2/tests
    > export PYTHONPATH=${PYTHONPATH}:/option3/lib/python2.7/site-packages
    > python test1.py
    > python test2.py
    ...

If you compiled the `c++` library with `-DMKL=ON`, you might get the error

    Intel MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so.

This issue of using Intel’s MKL inside python is known and reported. To
get the python tests to run, you can set the variable `LD_PRELOAD` in
order to preload libmkl\_rt.so. On my system, this is done with

    > export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so

The python tests do exactly the same thing as the `c++` tests above, and
illustrate the usage of the python interface to libchemps2. The tests
should end with a line stating whether or not they succeeded. Note that the
tests are too small to see (near) linear scaling with the number of cores,
although improvement should still be noticeable.


List of files in CheMPS2
------------------------

The files in the CheMPS2 library, as well as the files required to perform
test runs, can be found in [FILES.md](FILES.md).


