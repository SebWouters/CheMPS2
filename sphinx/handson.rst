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

.. index:: Hands-on

DMRG workshop (12-jul-2016): hands-on session
=============================================

Introduction
------------

The geometry of tetracene was optimized at the restricted B3LYP/6-31G* level of theory, and can be found in the file `tetracene.fcidump.in <https://github.com/sebwouters/chemps2/raw/master/sphinx/tetracene.fcidump.in>`_:

.. literalinclude:: tetracene.fcidump.in

The goal of this afternoon is to calculate the vertical singlet-triplet gap with DMRG(18, 18)-CASPT2/6-31G*.

`chemps2 <https://github.com/sebwouters/chemps2>`_ is a C++ library for spin-adapted DMRG calculations which can be incorporated in quantum chemistry packages. This has been done for `psi4 <http://www.psicode.org/>`_. Alternatively, the same functionality can be used with the binary, when the required matrix elements have been generated in ``FCIDUMP`` format. We will follow the second route this afternoon. The advantage of the latter route is that you are not tied to `psi4 <http://www.psicode.org/>`_ to obtain matrix elements. In the future you can use `molcas <http://www.molcas.org/>`_, `molpro <https://www.molpro.net/>`_, `dalton <http://www.daltonprogram.org/>`_... The disadvantage is that a full-rank ``FCIDUMP`` file is required, and that less virtual (secondary) orbitals can be used than with density-fitted DMRG-SCF and DMRG-CASPT2.

Please read an entire section before starting the instructions. Then you will have all the useful information you need!

UGent HPC
---------

Follow the `instructions <http://hpc.ugent.be/userwiki/index.php/User:VscConnect>`_ to log in to the UGent HPC.

Submit an interactive job in the XYZ queue:

.. code-block:: bash

    $ qsub -I -W x=FLAGS:ADVRES:dmrg.198 -l walltime=06:00:00 -l nodes=1:ppn=8

Once you are on the node, change directory to, for example, the following folder:

.. code-block:: bash

    $ cd $VSC_SCRATCH_NODE
    $ mkdir dmrg_workshop
    $ cd dmrg_workshop/

.. note::

    Please keep in mind that you will need about 20 GB of disk for the ``FCIDUMP`` file and 1 GB of disk for the `chemps2 <https://github.com/sebwouters/chemps2>`_ checkpoints!

``FCIDUMP`` and ``MOLDEN``
--------------------------

We will first use a plugin to `psi4 <http://www.psicode.org/>`_ to generate the RHF matrix elements in ``FCIDUMP`` format, as well as the corresponding ``MOLDEN`` file. As said before, any other program which is able to generate these two types of files can be used as well. Load the `psi4 <http://www.psicode.org/>`_ module and generate a new plugin called ``fcidump``:

.. code-block:: bash

    $ module load PSI4/1.0-intel-2016a-mt-Python-2.7.11
    $ psi4 --new-plugin fcidump

Overwrite the dummy file ``fcidump.cc`` and compile:

.. code-block:: bash

    $ cd fcidump/
    $ rm fcidump.cc
    $ wget 'https://github.com/sebwouters/chemps2/raw/master/integrals/psi4plugins/fcidump.cc'
    $ make
    $ cd ../

The required ``FCIDUMP`` file and the corresponding ``MOLDEN`` file can now be generated with `psi4 <http://www.psicode.org/>`_:

.. code-block:: bash

    $ wget 'https://github.com/sebwouters/chemps2/raw/master/sphinx/tetracene.fcidump.in'
    $ OMP_NUM_THREADS=8 psi4 -n 8 tetracene.fcidump.in &
    $ tail -n 3000 -f tetracene.fcidump.out
    $ ls -alh TETRACENE.FCIDUMP
    $ ls -alh tetracene.molden

.. note::

    The specified symmetry group in ``tetracene.fcidump.in`` was ``csz``, a subgroup of ``d2h``. In the ``csz`` symmetry group, the 18 active space :math:`\pi`-orbitals can be localized to the carbon atoms. This is not the case for the ``d2h`` symmetry group.

.. note::

    While you are waiting for the ``FCIDUMP`` file of size 20 GB, you can already proceed with the next section.

.. note::

    You can also use the precreated files from the folder ``/apps/gent/tutorials/DMRG/`` instead:
    
    .. code-block:: bash
    
        $ ls -al /apps/gent/tutorials/DMRG/tetracene.fcidump.in
        $ ls -al /apps/gent/tutorials/DMRG/tetracene.fcidump.out
        $ ls -al /apps/gent/tutorials/DMRG/TETRACENE.FCIDUMP
        $ ls -al /apps/gent/tutorials/DMRG/tetracene.molden

Basis choice
------------

Now that you have the required matrix elements in ``FCIDUMP`` format and the corresponding ``MOLDEN`` file, we can perform calculations with `chemps2 <https://github.com/sebwouters/chemps2>`_ v1.7.2. This module should have been loaded together with the `psi4 <http://www.psicode.org/>`_ module. If this was not the case, you can load it with:

.. code-block:: bash

    $ module load CheMPS2/1.7.2-intel-2016a

Study the options of the binary:

.. code-block:: bash

    $ chemps2 --version
    $ chemps2 --help

Perform each calculation in a separate folder. This way checkpoint files will not get mixed up. Create a folder ``ci_input_orbs/`` and in that folder an input file ``ci_input_orbs.in`` for `chemps2 <https://github.com/sebwouters/chemps2>`_ with the following options:

- Target the singlet ground state
- Use an (18, 18) active space
- Switch off the CASPT2 calculation
- Overwrite the tmp folder with the existing path ``/local/NUMBER.master15.delcatty.gent.vsc/``, where ``NUMBER`` is the job number which you see with

.. code-block:: bash

    $ qstat -n

- Perform one DMRG-SCF iteration, which corresponds to DMRG-CI
- The active space orbitals should be the RHF molecular orbitals (i.e. the input orbitals)
- Use the convergence scheme

 +-------------------+------------------+-----------------+------------------------+-----------------+
 | :math:`D_{SU(2)}` | :math:`E_{conv}` | :math:`N_{max}` | :math:`\gamma_{noise}` | :math:`r_{tol}` |
 +===================+==================+=================+========================+=================+
 | 200               | 1e-6             | 10              | 0.05                   | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 400               | 1e-6             | 10              | 0.05                   | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 600               | 1e-6             | 10              | 0.05                   | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 600               | 1e-8             | 3               | 0.0                    | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 400               | 1e-8             | 3               | 0.0                    | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 200               | 1e-8             | 3               | 0.0                    | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+

- Set the option ``SCF_MOLDEN`` to the corresponding molden file

When you have created the input file, you can double check with the :ref:`solution <first_ptr_solution>`.

Run the calculation:

.. code-block:: bash

    $ cd ci_input_orbs/
    $ OMP_NUM_THREADS=4 chemps2 --file=ci_input_orbs.in &> ci_input_orbs.out &
    $ cd ../

.. note::

    You have now only used 4 of the 8 available cores. Proceed with the instructions below while waiting for the calculation to finish.

Create a folder ``ci_local_orbs/`` and in that folder an input file ``ci_local_orbs.in`` for `chemps2 <https://github.com/sebwouters/chemps2>`_, which is identical to ``ci_input_orbs.in``, except for the active space orbitals. These should now be localized orbitals. When you have created the input file, you can double check with the :ref:`solution <second_ptr_solution>`.

Run the calculation:

.. code-block:: bash

    $ cd ci_local_orbs/
    $ OMP_NUM_THREADS=4 chemps2 --file=ci_local_orbs.in &> ci_local_orbs.out &
    $ cd ../
    $ tail -n 300 ci_input_orbs/ci_input_orbs.out
    $ tail -n 300 ci_local_orbs/ci_local_orbs.out

When the calculations are finished, take a look at the files

 - ``ci_input_orbs/tetracene.molden.rotated``
 - ``ci_local_orbs/tetracene.molden.rotated``

with your favourite visualization software. Do the first 18 ``App`` or ``A"`` orbitals have the desired shape? How are they ordered? Once you have formulated your own answer, you can double check with the :ref:`solution <third_ptr_solution>`.

Compare the energies of the last three sweep instructions as a function of :math:`D_{SU(2)}` for both calculations. Thereto you can grep for:

.. code-block:: bash

    $ grep "Minimum energy encountered during the last sweep" ci_input_orbs/ci_input_orbs.out
    $ grep "Minimum energy encountered during the last sweep" ci_local_orbs/ci_local_orbs.out

What do you observe? Can you explain it? Once you have formulated your own answer, you can double check with the :ref:`solution <fourth_ptr_solution>`.

DMRG-SCF
--------

Use localized orbitals for the active space from now on. Perform the DMRG-SCF orbital optimization for the singlet and the triplet. Also put DIIS on when the update norm is smaller than 1e-2, switch ``PRINT_CORR`` to ``TRUE``, and remove the ``SCF_MOLDEN`` line. Use the following convergence scheme:

 +-------------------+------------------+-----------------+------------------------+-----------------+
 | :math:`D_{SU(2)}` | :math:`E_{conv}` | :math:`N_{max}` | :math:`\gamma_{noise}` | :math:`r_{tol}` |
 +===================+==================+=================+========================+=================+
 | 250               | 1e-6             | 8               | 0.05                   | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 500               | 1e-8             | 8               | 0.05                   | 1e-5            |
 +-------------------+------------------+-----------------+------------------------+-----------------+
 | 750               | 1e-10            | 8               | 0.0                    | 1e-8            |
 +-------------------+------------------+-----------------+------------------------+-----------------+

Why is the reduced virtual dimension not lowered at the end of the DMRG calculation? Why is the last :math:`r_{tol}` smaller? When you have created the input files, you can double check with the solution for the :ref:`singlet <fifth_ptr_solution>` and the :ref:`triplet <sixth_ptr_solution>`.

Run the calculation in separate folders:

.. code-block:: bash

    $ cd scf_singlet/
    $ OMP_NUM_THREADS=4 chemps2 --file=scf_singlet.in &> scf_singlet.out &
    $ cd ../scf_triplet/
    $ OMP_NUM_THREADS=4 chemps2 --file=scf_triplet.in &> scf_triplet.out &
    $ cd ../
    $ tail -n 300 scf_singlet/scf_singlet.out
    $ tail -n 300 scf_triplet/scf_triplet.out

What is the DMRG-SCF singlet-triplet gap you obtain? Double check with the :ref:`solution <seventh_ptr_solution>`.

Do you see polyradical character in the natural orbital occupation numbers for the singlet and/or triplet? How can you observe this in the correlation functions? Tip: It might be interesting to read 

J. Hachmann, J. J. Dorando, Michael Avilés and Garnet Kin-Lic Chan, *Journal of Chemical Physics* **127**, 134309 (2007): `doi link <http://dx.doi.org/10.1063/1.2768362>`_ or `arXiv <https://arxiv.org/abs/0707.3120>`_

.. note::

    If you are in a hurry or immediately want to start with the DMRG-CASPT2 calculations, you can also use the precreated checkpoints from the folder ``/apps/gent/tutorials/DMRG/``:
    
    .. code-block:: bash
    
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_CASSCF.h5.singlet scf_singlet/.
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_CASSCF.h5.triplet scf_triplet/.

DMRG-CASPT2
-----------

.. note::

    DMRG-CASPT2 checkpoints can be used when you kill a DMRG-CASPT2 calculation before it is finished, or to redo the DMRG-CASPT2 calculation with another IPEA or IMAG shift. In case you would like to use checkpoints for the DMRG-CASPT2 calculations, it is important that for subsequent runs **exactly** the same orbitals are used. Therefore, start from the converged DMRG-SCF checkpoint ``CheMPS2_CASSCF.h5`` and do the following things:
    
     - Put ``SCF_DIIS_THR`` to ``0.0``
     - Delete any checkpoints named ``CheMPS2_DIIS.h5``
     - Switch ``SCF_ACTIVE_SPACE`` to ``I``
    
    This ensures that for the subsequent DMRG-CASPT2 runs, **exactly** the orbitals from ``CheMPS2_CASSCF.h5`` are used.

How large is the singlet-triplet gap with DMRG-CASPT2 when an IPEA shift of 0.0 and an IMAG shift of 0.0 are used? Is it best to use ``A`` or ``P`` for the option ``CASPT2_ORBS``, and why? In your input file, also switch on the DMRG-CASPT2 checkpoint, because later we will redo the calculation with an IPEA shift of 0.25. Use the same convergence scheme as for the DMRG-SCF calculations.

.. note::

    Sometimes a larger virtual dimension can be required for DMRG-CASPT2 as compared to DMRG-SCF, because the excited wavefunctions
    
    .. math::
    
        \left| sz, \alpha, \beta \right\rangle = \left[ \alpha \left( \hat{E}_{sz} + \hat{E}_{zs} \right) + \beta \right] \left| \Psi_0 \right\rangle
    
    are a linear combination over three matrix product states: :math:`\left| \Psi_0 \right\rangle`, :math:`\hat{E}_{sz} \left| \Psi_0 \right\rangle`, and :math:`\hat{E}_{zs} \left| \Psi_0 \right\rangle`. In practice, you should therefore check how the DMRG-CASPT2 second order energy in `chemps2 <https://github.com/sebwouters/chemps2>`_ varies with :math:`D_{SU(2)}`!

When you have created the input files, you can double check with the solution for the :ref:`singlet <eigth_ptr_solution>` and the :ref:`triplet <nineth_ptr_solution>`.

Run the calculations, but please remember to copy over the converged DMRG-SCF orbitals:

.. code-block:: bash

    $ cd pt2_singlet/
    $ cp ../scf_singlet/CheMPS2_CASSCF.h5 .
    $ OMP_NUM_THREADS=4 chemps2 --file=pt2_singlet.in &> pt2_singlet.out &
    $ cd ../pt2_triplet/
    $ cp ../scf_triplet/CheMPS2_CASSCF.h5 .
    $ OMP_NUM_THREADS=4 chemps2 --file=pt2_triplet.in &> pt2_triplet.out &
    $ cd ../
    $ tail -n 300 pt2_singlet/pt2_singlet.out
    $ tail -n 300 pt2_triplet/pt2_triplet.out

How large is the singlet-triplet gap with DMRG-CASPT2 when an IPEA shift of 0.0 and an IMAG shift of 0.0 are used?

And with an IPEA shift of 0.25 and an IMAG shift of 0.0?

You can double check with the :ref:`solution <tenth_ptr_solution>`.

.. note::

    You will see
    
    .. code-block:: bash
    
        CheMPS2::DMRG::Symm4RDM( X , Y ) : Elapsed wall time = Z seconds.
    
    appear in the output, with X and Y integers, and Z a floating point number. An estimate for the total wall time for the contraction of the 4-RDM with the Fock matrix is :math:`\frac{18 (18 + 1)}{2} Z` seconds.
    
    **So this last exercise is homework!**
    
    Compile `chemps2 <https://github.com/sebwouters/chemps2>`_ on your institution's HPC (or ask your admin or Sebastian to), and submit a non-interactive job for the DMRG-CASPT2 calculations.
    
    Yes, I have tricked you into using `chemps2 <https://github.com/sebwouters/chemps2>`_ in the future!

.. note::

    You can also find the precreated DMRG-CASPT2 checkpoints and the corresponding output in the folder ``/apps/gent/tutorials/DMRG/``:
    
    .. code-block:: bash
    
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_f4rdm.h5.singlet pt2_singlet/.
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_MPS0.h5.singlet pt2_singlet/.
        $ less /apps/gent/tutorials/DMRG/pt2_singlet.out.0.0
        $ less /apps/gent/tutorials/DMRG/pt2_singlet.out.0.25
        
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_f4rdm.h5.triplet pt2_triplet/.
        $ cp /apps/gent/tutorials/DMRG/CheMPS2_MPS0.h5.triplet pt2_triplet/.
        $ less /apps/gent/tutorials/DMRG/pt2_triplet.out.0.0
        $ less /apps/gent/tutorials/DMRG/pt2_triplet.out.0.25

To study an example DMRG-CASPT2 output during this workshop, perform a small active space calculation. For example:

.. code-block:: bash

    $ mkdir small_caspt2/
    $ cd small_caspt2/
    $ wget 'https://github.com/sebwouters/chemps2/raw/master/tests/matrixelements/N2.CCPVDZ.FCIDUMP'
    $ wget 'https://github.com/sebwouters/chemps2/raw/master/tests/test14.input'
    $ sed -i "s/\/path\/to/./" test14.input
    $ sed -i "s/\/tmp/\/local\/NUMBER.master15.delcatty.gent.vsc\//" test14.input
    $ cat test14.input
    $ chemps2 --file=test14.input &> test14.output &
    $ tail -n 3000 -f test14.output

Do you know the difference between the diagonal, non-variational, and variational second order perturbation energies? How is the reference weight calculated and what does it mean?

.. code-block:: bash

    $ grep "E2" test14.output
    $ grep "Reference weight" test14.output

Solutions
---------

.. _first_ptr_solution:

ci_input_orbs.in
~~~~~~~~~~~~~~~~

.. code-block:: bash

    FCIDUMP = /path/to/TETRACENE.FCIDUMP

    GROUP        = 3
    MULTIPLICITY = 1
    NELECTRONS   = 120
    IRREP        = 0
    EXCITATION   = 0

    SWEEP_STATES       = 200,  400,  600,  600,  400,  200
    SWEEP_ENERGY_CONV  = 1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8
    SWEEP_MAX_SWEEPS   = 10,   10,   10,   3,    3,    3
    SWEEP_NOISE_PREFAC = 0.05, 0.05, 0.05, 0.0,  0.0,  0.0
    SWEEP_DVDSON_RTOL  = 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5

    NOCC = 51,  0
    NACT = 0,   18
    NVIR = 171, 54

    SCF_STATE_AVG    = FALSE
    SCF_DIIS_THR     = 0.0
    SCF_GRAD_THR     = 1e-6
    SCF_MAX_ITER     = 1
    SCF_ACTIVE_SPACE = I
    SCF_MOLDEN       = /path/to/tetracene.molden

    CASPT2_CALC    = FALSE
    CASPT2_ORBS    = A
    CASPT2_IPEA    = 0.0
    CASPT2_IMAG    = 0.0
    CASPT2_CHECKPT = FALSE
    CASPT2_CUMUL   = FALSE

    PRINT_CORR = TRUE

    TMP_FOLDER = /local/NUMBER.master15.delcatty.gent.vsc/

.. _second_ptr_solution:

ci_local_orbs.in
~~~~~~~~~~~~~~~~

Difference with :ref:`input orbitals <first_ptr_solution>`:

.. code-block:: bash

    SCF_ACTIVE_SPACE = L

.. _third_ptr_solution:

tetracene.molden.rotated for the localized active space orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: handson_orbitals.png

For ``ci_local_orbs/tetracene.molden.rotated``, the active space orbitals are localized on the carbon atoms, and are ordered according to the one-dimensional topology of the molecule.

.. _fourth_ptr_solution:

Molecular vs. localized orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: handson_comparison.png

.. _fifth_ptr_solution:

scf_singlet.in
~~~~~~~~~~~~~~

.. code-block:: bash

    FCIDUMP = /path/to/TETRACENE.FCIDUMP

    GROUP        = 3
    MULTIPLICITY = 1
    NELECTRONS   = 120
    IRREP        = 0
    EXCITATION   = 0

    SWEEP_STATES       = 250,  500,  750
    SWEEP_ENERGY_CONV  = 1e-6, 1e-8, 1e-10
    SWEEP_MAX_SWEEPS   = 8,    8,    8
    SWEEP_NOISE_PREFAC = 0.05, 0.05, 0.0
    SWEEP_DVDSON_RTOL  = 1e-5, 1e-5, 1e-8

    NOCC = 51,  0
    NACT = 0,   18
    NVIR = 171, 54

    SCF_STATE_AVG    = FALSE
    SCF_DIIS_THR     = 1e-2
    SCF_GRAD_THR     = 1e-6
    SCF_MAX_ITER     = 100
    SCF_ACTIVE_SPACE = L

    CASPT2_CALC    = FALSE
    CASPT2_ORBS    = A
    CASPT2_IPEA    = 0.0
    CASPT2_IMAG    = 0.0
    CASPT2_CHECKPT = FALSE
    CASPT2_CUMUL   = FALSE

    PRINT_CORR = TRUE

    TMP_FOLDER = /local/NUMBER.master15.delcatty.gent.vsc/

.. _sixth_ptr_solution:

scf_triplet.in
~~~~~~~~~~~~~~

Difference with :ref:`singlet <fifth_ptr_solution>`:

.. code-block:: bash

    MULTIPLICITY = 3

.. _seventh_ptr_solution:

DMRG-SCF singlet-triplet gap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both DMRG-SCF calculations are converged with 8 macro-iterations. The gap is

.. math::

    \Delta E = E_{triplet} - E_{singlet} = -688.803387 - (-688.867150) ~ E_{h} = 63.763 ~ mE_{h} = 40.012 ~ kcal/mol

.. _eigth_ptr_solution:

pt2_singlet.in
~~~~~~~~~~~~~~

.. code-block:: bash

    FCIDUMP = /path/to/TETRACENE.FCIDUMP

    GROUP        = 3
    MULTIPLICITY = 1
    NELECTRONS   = 120
    IRREP        = 0
    EXCITATION   = 0

    SWEEP_STATES       = 250,  500,  750
    SWEEP_ENERGY_CONV  = 1e-6, 1e-8, 1e-10
    SWEEP_MAX_SWEEPS   = 8,    8,    8
    SWEEP_NOISE_PREFAC = 0.05, 0.05, 0.0
    SWEEP_DVDSON_RTOL  = 1e-5, 1e-5, 1e-8

    NOCC = 51,  0
    NACT = 0,   18
    NVIR = 171, 54

    SCF_STATE_AVG    = FALSE
    SCF_DIIS_THR     = 0.0
    SCF_GRAD_THR     = 1e-6
    SCF_MAX_ITER     = 100
    SCF_ACTIVE_SPACE = I

    CASPT2_CALC    = TRUE
    CASPT2_ORBS    = A
    CASPT2_IPEA    = 0.0
    CASPT2_IMAG    = 0.0
    CASPT2_CHECKPT = TRUE
    CASPT2_CUMUL   = FALSE

    PRINT_CORR = TRUE

    TMP_FOLDER = /local/NUMBER.master15.delcatty.gent.vsc/

.. _nineth_ptr_solution:

pt2_triplet.in
~~~~~~~~~~~~~~

Difference with :ref:`singlet <eigth_ptr_solution>`:

.. code-block:: bash

    MULTIPLICITY = 3

.. _tenth_ptr_solution:

DMRG-CASPT2 singlet-triplet gaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For IPEA = 0.0 a.u. the gap is

.. math::

    \Delta E = E_{triplet} - E_{singlet} = -690.945663 - (-691.000608) ~ E_{h} = 54.945 ~ mE_{h} = 34.479 ~ kcal/mol

and for IPEA = 0.25 a.u. the gap is

.. math::

    \Delta E = E_{triplet} - E_{singlet} = -690.923604 - (-690.987560) ~ E_{h} = 63.956 ~ mE_{h} = 40.133 ~ kcal/mol


