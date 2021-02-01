[![Build Status](https://travis-ci.org/SebWouters/CheMPS2.svg?branch=master)](https://travis-ci.org/SebWouters/CheMPS2)

![CheMPS2 Logo](https://raw.githubusercontent.com/sebwouters/CheMPS2/master/CheMPS2/CheMPS2logo.png?raw=true)

CheMPS2: a spin-adapted implementation of DMRG for ab initio quantum chemistry
==============================================================================

Copyright (C) 2013-2021 Sebastian Wouters <sebastianwouters@gmail.com>

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
library is illustrated in the
[c++](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-libchemps2) and
[python](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-pychemps2)
tests.

The CheMPS2 binary allows to perform DMRG-SCF and DMRG-CASPT2 calculations
based on input files. Molecular orbital matrix elements should then be
provided in FCIDUMP format. Usage of the binary is illustrated in the
[binary](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-the-chemps2-binary)
example.


Publications
------------

To acknowledge CheMPS2, please cite the publications in
[sphinx/publications.rst](sphinx/publications.rst). The
corresponding compiled html documentation can be consulted
[here](http://sebwouters.github.io/CheMPS2/publications.html).


Documentation
-------------

Information on CheMPS2 can be found in the following documents:

1. [publications](http://sebwouters.github.io/CheMPS2/publications.html)
2. [user manual](http://sebwouters.github.io/CheMPS2/index.html)
3. [workshop video](https://www.youtube.com/watch?v=U96atV5Akx4)
4. [doxygen html output](http://sebwouters.github.io/CheMPS2/doxygen/index.html)

The [user manual](http://sebwouters.github.io/CheMPS2/index.html) contains
elaborate information on

* the installation of CheMPS2
* the DMRG, DMRG-SCF, and DMRG-CASPT2 algorithms
* the symmetries which are exploited in CheMPS2
* how to generate matrix elements with plugins to [psi4](http://www.psicode.org)
* how to perform DMRG, DMRG-SCF, and DMRG-CASPT2 calculations
* the interfaces to [psi4](http://www.psicode.org) and [pyscf](https://github.com/sunqm/pyscf)


Installation
------------

Information on how to build and install CheMPS2 from source can
be found at [sphinx/sourcecode.rst](sphinx/sourcecode.rst). The
corresponding compiled html documentation can be consulted
[here](http://sebwouters.github.io/CheMPS2/sourcecode.html).


List of files
-------------

The files in the CheMPS2 library, as well as the files required
to perform test runs, can be found in [FILES.md](FILES.md).


