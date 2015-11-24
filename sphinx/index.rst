=====================
CheMPS2 documentation
=====================

CheMPS2 is a scientific library which contains a spin-adapted implementation of the density matrix renormalization group (DMRG) for ab initio quantum chemistry. This method allows to obtain numerical accuracy in active spaces beyond the capabilities of full configuration interaction (FCI):

* up to 40 electrons in 40 orbitals for general active spaces
* up to 100 electrons in 100 orbitals for one-dimensional active spaces, such as the :math:`\pi`-system of all-trans polyenes

In addition, DMRG allows to obtain the 2-RDM of the active space efficiently. The method is therefore ideal to replace the FCI solver in the complete active space self consistent field (CASSCF) method, when the active space size becomes prohibitively expensive for FCI. The corresponding method is called DMRG-SCF. Because DMRG can handle the abovementioned active space sizes, it allows to obtain FCI energies for small systems such as dimers, while for larger systems it is ideal to treat the static/strong correlation in a large active space.

CheMPS2 is designed to be a high-performance library. For an input Hamiltonian and targeted symmetry sector, the library performs successive DMRG sweeps according to a user-defined convergence scheme. As output, the library returns the minimal encountered energy as well as the 2-RDM, 3-RDM, and various orbital correlation functions in the active space. With the 2-RDM, various molecular properties can be calculated, as well as the gradient and Hessian for orbital rotations or nuclear displacements.

CheMPS2 is parallelized for shared memory architectures with OpenMP and for distributed memory architectures with MPI. A hybrid combination of both parallelization strategies is supported.

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

