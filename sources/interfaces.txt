.. index:: Matrix elements
.. index:: psi4
.. index:: pyscf

Interfaces to psi4 and pyscf
============================

CheMPS2 is currently interfaced with two ab initio quantum chemistry packages:

#. The ``dmrgci`` and ``dmrgscf`` plugins to `psi4 <http://www.psicode.org/>`_ allow to perform DMRG-CI and DMRG-SCF calculations with CheMPS2, by using `psi4 <http://www.psicode.org/>`_ input files. Since April 2015, CheMPS2 is also an integral part of psi4.

#. There are also DMRG-CI and DMRG-SCF interfaces between `pyscf <http://sunqm.github.io/pyscf/>`_ and CheMPS2.


psi4 ``dmrgci`` plugin
----------------------

DMRG-CI calculations can be performed directly with `psi4 <http://www.psicode.org/>`_. The plugin has been tested on commit ``14c78eabdca86f8e094576890518d93d300d2500`` (February 27, 2015) on https://github.com/psi4/psi4public, and should work on later versions as well.

To perform DMRG-CI calculations, build `psi4 <http://www.psicode.org/>`_ with the plugin option, and then run:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 --new-plugin dmrgci
    $ cd dmrgci

Now, replace the file ``dmrgci.cc`` with ``/sourcefolder/chemps2/integrals/psi4plugins/dmrgci.cc``. To compile the plugin, the Makefile should be adjusted. Change the line

.. code-block:: bash

    PSIPLUGIN = -L$(OBJDIR)/lib -lplugin

to

.. code-block:: bash

    PSIPLUGIN = -L$(OBJDIR)/lib -lplugin -lchemps2

Remember to add the library and include paths to the Makefile as well, if ``libchemps2`` is not installed in a standard location.
To compile the plugin, run:

.. code-block:: bash

    $ make

An example input file to use the ``dmrgci`` plugin:

.. literalinclude:: H2O.dmrgci.in

This file (``H2O.dmrgci.in``) should be placed in the folder ``/mypsi4plugins``. The DMRG-CI calculation can then be started with:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 H2O.dmrgci.in H2O.dmrgci.out
    
Since April 2015, CheMPS2 is also an integral part of `psi4 <http://www.psicode.org/>`_. Please consult `psi4 <http://www.psicode.org/>`_'s documentation on how to run DMRG-CI calculations with `psi4 <http://www.psicode.org/>`_.


psi4 ``dmrgscf`` plugin
-----------------------

DMRG-SCF calculations can be performed directly with `psi4 <http://www.psicode.org/>`_. The plugin has been tested on commit ``14c78eabdca86f8e094576890518d93d300d2500`` (February 27, 2015) on https://github.com/psi4/psi4public, and should work on later versions as well.

To perform DMRG-SCF calculations, build `psi4 <http://www.psicode.org/>`_ with the plugin option, and then run:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 --new-plugin dmrgscf
    $ cd dmrgscf

Now, replace the file ``dmrgscf.cc`` with ``/sourcefolder/chemps2/integrals/psi4plugins/dmrgscf.cc``. To compile the plugin, the Makefile should be adjusted. Change the line

.. code-block:: bash

    PSIPLUGIN = -L$(OBJDIR)/lib -lplugin

to

.. code-block:: bash

    PSIPLUGIN = -L$(OBJDIR)/lib -lplugin -lchemps2

Remember to add the library and include paths to the Makefile as well, if ``libchemps2`` is not installed in a standard location.
To compile the plugin, run:

.. code-block:: bash

    $ make

An example input file to use the ``dmrgscf`` plugin:

.. literalinclude:: O2.dmrgscf.in

This file (``O2.dmrgscf.in``) should be placed in the folder ``/mypsi4plugins``. The DMRG-CI calculation can then be started with:

.. code-block:: bash

    $ cd /mypsi4plugins
    $ psi4 O2.dmrgscf.in O2.dmrgscf.out

Since April 2015, CheMPS2 is also an integral part of `psi4 <http://www.psicode.org/>`_. Please consult `psi4 <http://www.psicode.org/>`_'s documentation on how to run DMRG-SCF calculations with `psi4 <http://www.psicode.org/>`_.


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


