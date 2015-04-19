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
* GSL (`GNU Scientific library <http://www.gnu.org/software/gsl/>`_)
* HDF5 (`Hierarchical Data Format Release 5 <http://www.hdfgroup.org/HDF5/>`_)

It is parallelized for shared memory architectures with the Open Multi-Processing (`OpenMP <http://openmp.org/wp/>`_) API. In the future, the parallelization of CheMPS2 for shared memory architectures will be extended to a hybrid scheme with both shared (OpenMP) and distributed (MPI) memory parallelization.

Download
--------

It is advised to clone the CheMPS2 git repository from github. In your terminal, run:

.. code-block:: bash

    $ cd /myfolder
    $ git clone 'https://github.com/sebwouters/chemps2'
    $ cd chemps2
    
That way, future updates and bug fixes can be easily pulled in:

.. code-block:: bash

    $ cd /myfolder/chemps2
    $ git pull

Build libchemps2
----------------

The files

.. code-block:: bash

    /myfolder/chemps2/CMakeLists.txt
    /myfolder/chemps2/CheMPS2/CMakeLists.txt
    /myfolder/chemps2/tests/CMakeLists.txt
    /myfolder/chemps2/PyCheMPS2/CMakeLists.txt
    /myfolder/chemps2/sphinx/CMakeLists.txt

provide a minimal compilation. In your terminal, run:

.. code-block:: bash

    $ cd /myfolder/chemps2
    $ mkdir build
    $ cd build

CMake generates makefiles based on the user's specifications:

.. code-block:: bash

    $ CXX=option1 cmake .. -DMKL=option2 -DBUILD_DOCUMENTATION=option3 -DCMAKE_INSTALL_PREFIX=option4 -DSTATIC_ONLY=option5 -DBUILD_SPHINX=option6

#. Option1 is the ``c++`` compiler; typically ``g++``, ``icpc``, or ``clang++`` on Linux.
#. Option2 can be ``ON`` or ``OFF`` and is used to switch on the intel math kernel library.
#. Option3 can be ``ON`` or ``OFF`` and is used to switch on the possibility to compile the doxygen documentation.
#. Option4 is the prefix of the installation directory; typically ``/usr`` or ``/usr/local`` on Linux. On my computer, libchemps2 is then installed in ``/prefix/lib/x86_64-linux-gnu`` and the headers in ``/prefix/include/chemps2``.
#. Option5 can be ``ON`` or ``OFF`` and is used to avoid building the shared library (please use the default here).
#. Option6 can be ``ON`` or ``OFF`` and is used to switch on the possibility to compile the user manual with sphinx.

If one or more of the required libraries are not found, use the command

.. code-block:: bash

    $ CMAKE_INCLUDE_PATH=option7 CMAKE_LIBRARY_PATH=option8 CXX=option1 cmake .. -DMKL=option2 -DBUILD_DOCUMENTATION=option3 -DCMAKE_INSTALL_PREFIX=option4 -DSTATIC_ONLY=option5 -DBUILD_SPHINX=option6

instead, where option7 and option8 are respectively the missing colon-separated include and library paths:

.. code-block:: bash
    
    CMAKE_INCLUDE_PATH=/my_libs/lib1/include:/my_libs/lib2/include
    CMAKE_LIBRARY_PATH=/my_libs/lib1/lib:/my_libs/lib2/lib

To compile, run:

.. code-block:: bash
    
    $ make

To install, run:

.. code-block:: bash
    
    $ make install

For non-standard installation directories, please remember to append the library path to ``LD_LIBRARY_PATH`` in your ``.bashrc``.

Test libchemps2
---------------

To test libchemps2, run:

.. code-block:: bash
    
    $ cd /myfolder/chemps2/build
    $ make test

The tests only require a very limited amount of memory (order 10-120 MB).

Build PyCheMPS2
---------------

PyCheMPS2, a python inferface to libchemps2, can be built with `Cython <http://cython.org/>`_. The installation above generated the file ``/myfolder/chemps2/build/setup.py``, in which the library and include paths have been set to the ones in the build directory. In your terminal, run:

.. code-block:: bash
    
    $ cd /myfolder/chemps2/build
    $ python setup.py build_ext -i

Compilation of PyCheMPS2 occurs by linking to the ``c++`` library in the build directory. However, when you want to run the tests in the next section, you should first install the ``c++`` library! If you have pulled a newer version of CheMPS2, please remember to (re)install the ``c++`` library first, before using PyCheMPS2!

Test PyCheMPS2
--------------

To test PyCheMPS2, run:

.. code-block:: bash
    
    $ cd /myfolder/chemps2/build/PyTests
    $ python test1.py
    $ python test2.py
    $ python test3.py
    $ python test4.py
    $ python test5.py
    $ python test6.py
    $ python test7.py
    $ python test8.py
    $ python test9.py

If you compiled the ``c++`` library with ``-DMKL=ON``, you might get the error

.. code-block:: bash

    Intel MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so.

This issue of using Intel's MKL inside python is known and reported. To get the python tests to run, you can set the variable ``LD_PRELOAD`` in order to preload libmkl_rt.so. On my system, this is done with

.. code-block:: bash

    $ export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so

The python tests do exactly the same thing as the ``c++`` tests above, and illustrate the usage of the python interface to libchemps2. The tests should end with a line stating whether or not they succeeded. They only require a very limited amount of memory (order 10-120 MB).

Doxygen
-------

To build the doxygen manual, the ``BUILD_DOCUMENTATION`` flag should have been on: ``-DBUILD_DOCUMENTATION=ON``. In your terminal, run:

.. code-block:: bash
    
    $ cd /myfolder/chemps2/build
    $ make doc
    $ cd LaTeX-documents
    $ make
    $ evince refman.pdf &
    $ cd ../html
    $ firefox index.html &
    
The `doxygen html output <http://sebwouters.github.io/CheMPS2/doxygen/index.html>`_ can also be consulted online.

Sphinx user manual
------------------

To build the sphinx user manual, the ``BUILD_SPHINX`` flag should have been on: ``-DBUILD_SPHINX=ON``. In your terminal, run:

.. code-block:: bash

    $ cd /myfolder/chemps2/build
    $ make sphinx
    $ cd sphinx/html
    $ firefox index.html &

The `sphinx user manual <http://sebwouters.github.io/CheMPS2/index.html>`_ can also be consulted online.


