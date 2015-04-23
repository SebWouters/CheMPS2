=====================
CheMPS2 documentation
=====================

CheMPS2 is a scientific library which contains a spin-adapted implementation of the density matrix renormalization group (DMRG) for ab initio quantum chemistry. This method allows to obtain numerical accuracy in active spaces beyond the capabilities of full configuration interaction (FCI):

* up to 40 electrons in 40 orbitals for general active spaces
* up to 100 electrons in 100 orbitals for one-dimensional active spaces, such as the :math:`\pi`-system of all-trans polyenes

In addition, DMRG allows to obtain the 2-RDM of the active space efficiently. The method is therefore ideal to replace the FCI solver in the complete active space configuration interaction (CASCI) and complete active space self consistent field (CASSCF) methods, when the active space sizes become prohibitively expensive for FCI. The corresponding methods are called DMRG-CI and DMRG-SCF, respectively. Because DMRG can handle the abovementioned active space sizes, it allows to obtain FCI energies for small systems such as dimers, while for larger systems it is ideal to treat the static/strong correlation in a large active space.

The design philosophy for CheMPS2 is to be a lightweight, efficient, and stable library. For an input Hamiltonian and targeted symmetry sector, the library performs successive DMRG sweeps according to a user-defined convergence scheme. As output, the library returns the minimal encountered energy as well as the 2-RDM of the active space. With the latter, various molecular properties can be calculated, as well as the gradient and Hessian for orbital rotations or nuclear displacements. In addition, several correlation functions can be obtained to investigate the electronic structure in the active space.

In the future, the parallelization of CheMPS2 for shared memory architectures will be extended to a hybrid scheme with both shared (OpenMP) and distributed (MPI) memory parallelization. In addition, the calculation of the 2-RDM of the active space will be extended to the 3-RDM.

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
   interfaces
   
   
For more elaborate information on CheMPS2, please consult the :ref:`publications <label-publications>` and the `doxygen html output <http://sebwouters.github.io/CheMPS2/doxygen/index.html>`_.

