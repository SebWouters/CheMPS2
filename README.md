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
library is illustrated in the
[c++](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-libchemps2) and
[python](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-pychemps2)
tests.

The CheMPS2 binary allows to perform DMRG-SCF and DMRG-CASPT2 calculations
based on input files. Molecular orbital matrix elements should then be
provided in FCIDUMP format. Usage of the binary is illustrated in the
[binary](http://sebwouters.github.io/CheMPS2/sourcecode.html#test-the-chemps2-binary)
example.


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

4. S. Wouters, V. Van Speybroeck and D. Van Neck,
   DMRG-CASPT2 study of the longitudinal static second hyperpolarizability of all-trans polyenes,
   *Preprint* (2016),
   [arXiv:1605.05526](http://arxiv.org/abs/1605.05526)

        @article{CheMPS2cite4,
            author = {Sebastian Wouters and Veronique {Van Speybroeck} and Dimitri {Van Neck}},
            title = {DMRG-CASPT2 study of the longitudinal static second hyperpolarizability of all-trans polyenes},
            journal = {arXiv:1605.05526},
            year = {2016},
            url = {http://arxiv.org/abs/1605.05526}
        }


Installation
------------

Information on how to build and install CheMPS2 from source can be found at
[sphinx/sourcecode.rst](sphinx/sourcecode.rst), or alternatively in compiled
form at the
[user manual](http://sebwouters.github.io/CheMPS2/sourcecode.html).


List of files in CheMPS2
------------------------

The files in the CheMPS2 library, as well as the files required to perform
test runs, can be found in [FILES.md](FILES.md).


