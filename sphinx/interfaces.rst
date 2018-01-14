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

.. index:: Matrix elements
.. index:: psi4
.. index:: pyscf

Interfaces to psi4 and pyscf
============================

CheMPS2 is currently interfaced with two ab initio quantum chemistry packages:

#. The ``dmrg`` plugin to `psi4 <http://www.psicode.org/>`_ allows to perform DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations with CheMPS2. Since April 2015, CheMPS2 is also an integral part of psi4.

#. There are also DMRG-CI and DMRG-SCF interfaces between `pyscf <http://sunqm.github.io/pyscf/>`_ and CheMPS2.


psi4 ``dmrg`` plugin
--------------------

DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations can be performed directly with `psi4 <http://www.psicode.org/>`_. The plugin has been tested on late March 2017 psi4 pre-1.1.
.. comment  `tag <https://github.com/psi4/psi4/releases/tag/0.5>`_ (released February 17, 2016).

Note that as of late June 2016, DMRG keywords in `psi4 <http://www.psicode.org/>`_ have been realigned to those of the chemps2 executable. A `translation table <https://github.com/psi4/psi4/issues/150#issuecomment-228951911>`_ is available.

To perform DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations, build `psi4 <http://www.psicode.org/>`_ or use the conda binary, and then run:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 --plugin-name dmrg
    $ cd dmrg

Now, replace the file ``plugin.cc`` with ``/sourcefolder/chemps2/integrals/psi4plugins/dmrg.cc`` (file can be named ``plugin.cc`` or ``dmrg.cc``, but if the latter, change ``CMakeLists.txt`` to match). To compile the plugin, the CMakeLists.txt should be adjusted. Change the line

.. code-block:: bash

    find_package(psi4 1.0 REQUIRED)

to

.. code-block:: bash

    find_package(psi4 1.0 REQUIRED COMPONENTS chemps2)

Now call ``psi4 --plugin-compile`` and execute the result. Additional variables can be passed to the ``cmake`` command (including `-DCheMPS2_DIR`), but none should be necessary. Avoid building against a Psi4 with PCMSolver enabled, as this will cause trouble with capturing stdout. To compile the plugin, run:

.. code-block:: bash

    $ make

An example input file to perform a DMRG-CI calculation with the ``dmrg`` plugin:

.. literalinclude:: H2O.dmrgci.in

Note that the option ``dmrg_max_iter`` has been set to ``1``, so that only one active space calculation is performed. This file (``H2O.dmrgci.in``) should be placed in the folder ``/mypsi4plugins/dmrg``. The DMRG-CI calculation can then be started with:

.. code-block:: bash

    $ cd /mypsi4plugins/dmrg
    $ psi4 H2O.dmrgci.in H2O.dmrgci.out

An example input file to perform a DMRG-SCF calculation with the ``dmrg`` plugin:

.. literalinclude:: O2.dmrgscf.in

This file (``O2.dmrgscf.in``) should be placed in the folder ``/mypsi4plugins/dmrg``. The DMRG-SCF calculation can then be started with:

.. code-block:: bash

    $ cd /mypsi4plugins/dmrg
    $ psi4 O2.dmrgscf.in O2.dmrgscf.out

An example input file to perform a DMRG-CASPT2 calculation with the ``dmrg`` plugin:

.. literalinclude:: N2.caspt2.in

This file (``N2.caspt2.in``) should be placed in the folder ``/mypsi4plugins/dmrg``. The DMRG-CASPT2 calculation can then be started with:

.. code-block:: bash

    $ cd /mypsi4plugins/dmrg
    $ psi4 N2.caspt2.in N2.caspt2.out

Since April 2015, CheMPS2 is also an integral part of `psi4 <http://www.psicode.org/>`_. Please consult `psi4 <http://www.psicode.org/>`_'s documentation on how to run DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations with `psi4 <http://www.psicode.org/>`_.


pyscf
-----

`pyscf <http://sunqm.github.io/pyscf/>`_ is a new quantum chemistry package, in which all layers are written or interfaced in python. In the future, the package will be able to perform DMRG-CI and DMRG-SCF calculations using PyCheMPS2: `chemps2.py <https://github.com/sunqm/pyscf/blob/master/future/dmrgscf/chemps2.py>`_.

Examples of how to extract MO integrals from `pyscf <http://sunqm.github.io/pyscf/>`_ to perform DMRG-CI calculations with PyCheMPS2 can be found in:

#. ``/sourcefolder/chemps2/integrals/pyscf/example.py``
#. ``/sourcefolder/chemps2/integrals/pyscf/example2.py``
#. ``/sourcefolder/chemps2/integrals/pyscf/example3.py``
#. ``/sourcefolder/chemps2/integrals/pyscf/dmrgci.py``
#. ``/sourcefolder/chemps2/integrals/pyscf/call_chemps2.py``

Please remember to append the correct pyscf and PyCheMPS2 directories to ``sys.path`` at the top of these files.


