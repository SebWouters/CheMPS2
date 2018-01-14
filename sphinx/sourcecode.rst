.. CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
   Copyright (C) 2013-2018 Sebastian Wouters

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

.. index:: Source code
.. index:: Download
.. index:: Installation
.. index:: Build
.. index:: Dependencies
.. index:: Tests

Installation
============

Dependencies
------------

CheMPS2 can be built with `CMake <http://www.cmake.org/>`_ and depends on

* BLAS
* LAPACK
* HDF5 (`Hierarchical Data Format Release 5 <http://www.hdfgroup.org/HDF5/>`_)

It is parallelized for shared memory architectures with the Open Multi-Processing (`OpenMP <http://openmp.org/wp/>`_) API and for distributed memory architectures with the Message Passing Interface (`MPI <http://www.mpi-forum.org/>`_). A hybrid combination of both parallelization strategies is supported.

Download
--------

It is advised to clone the CheMPS2 git repository from github. In your terminal, run:

.. code-block:: bash

    $ cd /sourcefolder
    $ git clone 'https://github.com/sebwouters/chemps2'
    $ cd chemps2
    
That way, future updates and bug fixes can be easily pulled in:

.. code-block:: bash

    $ cd /sourcefolder/chemps2
    $ git pull

Build the chemps2 library and binary
------------------------------------

The files

.. code-block:: bash

    /sourcefolder/chemps2/CMakeLists.txt
    /sourcefolder/chemps2/CheMPS2/CMakeLists.txt
    /sourcefolder/chemps2/tests/CMakeLists.txt
    /sourcefolder/chemps2/sphinx/CMakeLists.txt

provide a minimal compilation. In your terminal, run:

.. code-block:: bash

    $ cd /sourcefolder/chemps2
    $ mkdir build
    $ cd build

CMake generates makefiles based on the user's specifications:

.. code-block:: bash

    $ CXX=option1 cmake .. -DMKL=option2 -DCMAKE_INSTALL_PREFIX=/option3 -DWITH_MPI=option4

#.  Option1 is the ``c++`` compiler; typically ``g++``, ``icpc``, or ``clang++`` on Linux. It is advised to use the intel compiler if available.
#.  Option2 can be ``ON`` or ``OFF`` and is used to switch on the intel math kernel library.
#.  /option3 is the prefix of the installation directory; typically ``/usr`` or ``/usr/local`` on Linux. On my computer, libchemps2 is then installed in ``/option3/lib/x86_64-linux-gnu/``, the headers in ``/option3/include/chemps2/``, and the binary in ``/option3/bin/chemps2``.
#.  Option4 can be ``ON`` or ``OFF`` and is used to switch on the possibility to compile with MPI. Please note that the compiler should then provide ``mpi.h``. Option1 should hence be the ``mpic++`` compiler; typically ``mpic++`` or ``mpiicpc`` on Linux. It is advised to use the intel compiler if available.

If one or more of the required libraries are not found, use the command

.. code-block:: bash

    $ CMAKE_INCLUDE_PATH=option5 CMAKE_LIBRARY_PATH=option6 CXX=option1 cmake .. -DMKL=option2 -DCMAKE_INSTALL_PREFIX=/option3 -DWITH_MPI=option4

instead, where option5 and option6 are respectively the missing colon-separated include and library paths:

.. code-block:: bash
    
    CMAKE_INCLUDE_PATH=/my_libs/lib1/include:/my_libs/lib2/include
    CMAKE_LIBRARY_PATH=/my_libs/lib1/lib:/my_libs/lib2/lib

Remarks:

#.  For operating systems based on debian, the HDF5 headers are located in the folder ``/usr/include/hdf5/serial/``. If CMake complains about the HDF5 headers, try to pass it with the option ``-DHDF5_INCLUDE_DIRS=/usr/include/hdf5/serial``.
#.  Sometimes it might be necessary to specify the MKL libraries, e.g. for mixed-type GCC and single-threaded MKL compilation with the option ``-DLAPACK_LIBRARIES="/opt/intel/mkl/lib/intel64/libmkl_intel_lp64.so;/opt/intel/mkl/lib/intel64/libmkl_sequential.so;/opt/intel/mkl/lib/intel64/libmkl_core.so"``.
#.  For building with GCC, errors involving unresolved symbols or a message ``plugin needed to handle lto object`` may indicate a failure of the interprocedural optimization. This can be resolved by passing full locations to gcc toolchain utilites to the ``setup`` command above: ``-DCMAKE_RANLIB=/path/to/gcc-ranlib -DCMAKE_AR=/path/to/gcc-ar`` .

To compile, run:

.. code-block:: bash
    
    $ make

To install, run:

.. code-block:: bash
    
    $ make install

For non-standard installation directories, please remember to append the library path to ``LD_LIBRARY_PATH`` in your ``.bashrc``.

Test libchemps2
---------------

To test libchemps2 for compilation **without MPI**, run:

.. code-block:: bash
    
    $ cd /sourcefolder/chemps2/build
    $ make test
    
To test libchemps2 for compilation **with MPI**, run:

.. code-block:: bash
    
    $ cd /sourcefolder/chemps2/build/tests
    $ OMP_NUM_THREADS=YYY mpirun -np ZZZ ./test1
    $ OMP_NUM_THREADS=YYY mpirun -np ZZZ ./test2
    ...

``YYY`` specifies the number of threads per process and ``ZZZ`` the number of processes. Note that the tests are too small to see (near) linear scaling with the number of cores, although improvement should still be noticeable.

Test the chemps2 binary
-----------------------

To test the chemps2 binary for compilation **without MPI**, run:

.. code-block:: bash

    $ man /sourcefolder/chemps2/chemps2.1
    $ cd /sourcefolder/chemps2/build/CheMPS2
    $ ./chemps2 --help
    $ cp /sourcefolder/chemps2/tests/test14.input .
    $ sed -i "s/path\/to/sourcefolder\/chemps2\/tests\/matrixelements/" test14.input
    $ cat test14.input
    $ ./chemps2 --file=test14.input

Note that when you use the CASPT2 checkpoint, and want to restart a
calculation at a later point, you should

    1. switch the option ``SCF_ACTIVE_SPACE`` to ``I``
    2. remove the ``CheMPS2_DIIS.h5`` checkpoint

in order to ensure that **exactly** the same orbitals are used in the different runs.
    
To test the chemps2 binary for compilation **with MPI**, prepend the binary with:

.. code-block:: bash

    $ OMP_NUM_THREADS=YYY mpirun -np ZZZ ./chemps2 --file=test14.input

Build PyCheMPS2
---------------

PyCheMPS2 is a python interface to libchemps2, for compilation **without MPI**. It can be built with `Cython <http://cython.org/>`_. The installation is independent of CMake and assumes that you have installed the CheMPS2 library with ``make install``. For non-standard installation directories of CheMPS2, please remember to append the library path to ``LD_LIBRARY_PATH`` in your ``.bashrc``. In addition, the include path should be appended to ``CPATH``:

.. code-block:: bash

    $ export CPATH=${CPATH}:/option3/include
    
where ``/option3`` is the option provided to CMake with ``-DCMAKE_INSTALL_PREFIX=/option3`` above. For operating systems based on debian, the HDF5 headers are located in the folder ``/usr/include/hdf5/serial/``. If it was explicitly passed to CMake, it should also be appended to ``CPATH``:

.. code-block:: bash

    $ export CPATH=${CPATH}:/option3/include:/usr/include/hdf5/serial

The python wrapper can be installed with:

.. code-block:: bash

    $ cd /sourcefolder/chemps2/PyCheMPS2
    $ python setup.py build_ext -L ${LD_LIBRARY_PATH}
    $ python setup.py install --prefix=/option3

On my machine, the python wrapper is installed to the folder ``/option3/lib/python2.7/site-packages``, but the folder ``lib`` and the distribution of python can vary.

Compilation of PyCheMPS2 occurs by linking to the ``c++`` library in the installation directory. The installation of PyCheMPS2 will fail if that library is not properly installed. If you have pulled a newer version of CheMPS2, please remember to reinstall the ``c++`` library first, before reinstalling PyCheMPS2!

Test PyCheMPS2
--------------

When libchemps2 has been compiled **without MPI**, PyCheMPS2 can be tested by running (remember that the python site-packages folder can vary):

.. code-block:: bash

    $ cd /sourcefolder/chemps2/PyCheMPS2/tests
    $ export PYTHONPATH=${PYTHONPATH}:/option3/lib/python2.7/site-packages
    $ python test1.py
    $ python test2.py
    ...


If you compiled the ``c++`` library with ``-DMKL=ON``, you might get the error

.. code-block:: bash

    Intel MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so.

This issue of using Intel's MKL inside python is known and reported. To get the python tests to run, you can set the variable ``LD_PRELOAD`` in order to preload libmkl_rt.so. On my system, this is done with

.. code-block:: bash

    $ export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so

The python tests do exactly the same thing as the ``c++`` tests above, and illustrate the usage of the python interface to libchemps2. The tests should end with a line stating whether or not they succeeded. Note that the tests are too small to see (near) linear scaling with the number of cores, although improvement should still be noticeable.


