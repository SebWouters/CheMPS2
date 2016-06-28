.. index:: Matrix elements
.. index:: psi4
.. index:: pyscf

Interfaces to psi4 and pyscf
============================

CheMPS2 is currently interfaced with two ab initio quantum chemistry packages:

#. The ``dmrg`` plugin to `psi4 <http://www.psicode.org/>`_ allows to perform DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations with CheMPS2, by using `psi4 <http://www.psicode.org/>`_ input files. Since April 2015, CheMPS2 is also an integral part of psi4.

#. There are also DMRG-CI and DMRG-SCF interfaces between `pyscf <http://sunqm.github.io/pyscf/>`_ and CheMPS2.


psi4 ``dmrg`` plugin
--------------------

DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations can be performed directly with `psi4 <http://www.psicode.org/>`_. The plugin has been tested on `psi4-0.5 <https://github.com/psi4/psi4/releases/tag/0.5>`_ (released February 17, 2016).

Note that as of late June 2016, DMRG keywords in Psi4 have been realigned to those of the chemps2 executable. A `translation table <https://github.com/psi4/psi4/issues/150#issuecomment-228951911>`_ is available.

To perform DMRG-CI, DMRG-SCF, and DMRG-CASPT2 calculations, build `psi4 <http://www.psicode.org/>`_ or use the conda binary, and then run:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 --new-plugin dmrg
    $ cd dmrg

Now, replace the file ``dmrg.cc`` with ``/sourcefolder/chemps2/integrals/psi4plugins/dmrg.cc``. To compile the plugin, the Makefile should be adjusted. Change the line

.. code-block:: bash

    $(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS)

to

.. code-block:: bash

    $(CXX) $(LDFLAGS) -o $@ $^ $(CXXDEFS) -lchemps2

Remember to add the library and include paths to the Makefile as well, if ``libchemps2`` is not installed in a standard location. For debian/sid, the HDF5 headers are located in the folder ``/usr/include/hdf5/serial``. It might be necessary to add it to the ``INCLUDES`` variable in the Makefile.
To compile the plugin, run:

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


