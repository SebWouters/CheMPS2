=====================
CheMPS2 documentation
=====================

.. image:: CheMPS2logo.png

CheMPS2 is a scientific library in the field of ab initio quantum chemistry. It contains a spin-adapted implementation of the density matrix renormalization group (DMRG) for ab initio quantum chemistry. This method allows to obtain numerical accuracy in active spaces beyond the capabilities of full configuration interaction (FCI).

For an input Hamiltonian and targeted symmetry sector, the library performs successive DMRG sweeps according to a user-defined convergence scheme. As output, the library returns the minimal encountered energy as well as the 2-RDM of the active space. With the latter, various molecular properties can be calculated, as well as the gradient and Hessian for orbital rotations or nuclear displacements. In addition, several correlation functions can be obtained to investigate the electronic structure in the active space.

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
   dmrgscf
   dmrgscfcalcs
   interfaces
   
   
For more elaborate information on CheMPS2, please consult the :ref:`publications <label-publications>` and the `doxygen html output <http://sebwouters.github.io/CheMPS2/doxygen/index.html>`_.

