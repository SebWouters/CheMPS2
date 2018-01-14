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

=====================
CheMPS2 documentation
=====================

CheMPS2 is a scientific library which contains a spin-adapted implementation of the density matrix renormalization group (DMRG) for ab initio quantum chemistry. This wavefunction method allows to obtain numerical accuracy in active spaces beyond the capabilities of full configuration interaction (FCI), and allows to extract the 2-, 3-, and 4-particle reduced density matrices (2-, 3- and 4-RDM) of the active space.

For general active spaces up to 40 electrons in 40 orbitals can be handled with DMRG, and for one-dimensional active spaces up to 100 electrons in 100 orbitals. The 2-RDM of these active spaces can also be easily extracted, while the 3- and 4-RDM are limited to about 28 orbitals.

When the active space size becomes prohibitively expensive for FCI, DMRG can be used to replace the FCI solver in the complete active space self consistent field (CASSCF) method and the corresponding complete active space second order perturbation theory (CASPT2). The corresponding methods are called DMRG-SCF and DMRG-CASPT2, respectively. For DMRG-SCF the active space 2-RDM is required, and for DMRG-CASPT2 the active space 4-RDM.

CheMPS2 is designed for high-performance computers, with a hybrid parallelization for mixed distributed and shared memory architectures, realized with MPI and OpenMP.

CheMPS2 is distributed under the GNU General Public License version 2, and can be downloaded from https://github.com/sebwouters/chemps2.

Contents
========
.. toctree::
   :maxdepth: 2
   :numbered:

   sourcecode
   publications
   method
   symmetry
   matrixelements
   inoutput
   resources
   dmrgscf
   dmrgscfcalcs
   caspt2
   interfaces
   handson
   
For more elaborate information on CheMPS2, please consult the :ref:`publications <label-publications>`, `workshop video <https://www.youtube.com/watch?v=U96atV5Akx4>`_, and the `doxygen html output <http://sebwouters.github.io/CheMPS2/doxygen/index.html>`_.

