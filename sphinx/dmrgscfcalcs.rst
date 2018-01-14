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

DMRG-SCF calculations
=====================

``CheMPS2::DMRGSCFoptions``
---------------------------

An instance of the ``CheMPS2::DMRGSCFoptions`` class is required to perform DMRG-SCF calculations. It allows to overwrite the default algorithmic choices.

.. code-block:: c++

    CheMPS2::DMRGSCFoptions::DMRGSCFoptions()
    void CheMPS2::DMRGSCFoptions::setDoDIIS( const bool DoDIIS_in=false )
    void CheMPS2::DMRGSCFoptions::setDIISGradientBranch( const double DIISGradientBranch_in=1e-2 )
    void CheMPS2::DMRGSCFoptions::setNumDIISVecs( const int NumDIISVecs_in=7 )
    void CheMPS2::DMRGSCFoptions::setStoreDIIS( const bool StoreDIIS_in=true )
    void CheMPS2::DMRGSCFoptions::setDIISStorageName( const string DIISStorageName_in="CheMPS2_DIIS.h5" )
    void CheMPS2::DMRGSCFoptions::setMaxIterations( const int MaxIterations_in=100 )
    void CheMPS2::DMRGSCFoptions::setGradientThreshold( const double GradientThreshold_in=1e-6 )
    void CheMPS2::DMRGSCFoptions::setStoreUnitary( const bool StoreUnitary_in=true )
    void CheMPS2::DMRGSCFoptions::setUnitaryStorageName( const string UnitaryStorageName_in="CheMPS2_CASSCF.h5" )
    void CheMPS2::DMRGSCFoptions::setStateAveraging( const bool StateAveraging_in=true )
    void CheMPS2::DMRGSCFoptions::setWhichActiveSpace( const int WhichActiveSpace_in=0 )
    void CheMPS2::DMRGSCFoptions::setDumpCorrelations( const bool DumpCorrelations_in=false )
    void CheMPS2::DMRGSCFoptions::setStartLocRandom( const bool StartLocRandom_in=false )

* The variable ``DoDIIS_in`` allows to switch on the DIIS acceleration for DMRG-SCF calculations.
* If ``DoDIIS_in==true``, the DIIS acceleration starts when the update norm :math:`\|\vec{x}\|_2` is smaller than ``DIISGradientBranch_in``.
* The variable ``NumDIISVecs_in`` allows to set the size of the DIIS list.
* The variable ``StoreDIIS_in`` allows to switch off storing the DIIS checkpoint file to disk.
* The variable ``DIISStorageName_in`` contains the filename of the DIIS checkpoint.
* The variable ``MaxIterations_in`` contains the maximum number of DMRG-SCF iterations to be performed.
* The variable ``GradientThreshold_in`` specifies the convergence threshold for the norm of the orbital rotation gradient :math:`\|\vec{g}\|_2`.
* The variable ``StoreUnitary_in`` allows to switch off storing the orbital rotation checkpoint file to disk.
* The variable ``UnitaryStorageName_in`` contains the filename of the orbital rotation checkpoint.
* The variable ``StateAveraging_in`` allows to switch off state-averaged DMRG-SCF calculations (to perform state-specific DMRG-SCF calculations).
* The variable ``WhichActiveSpace_in`` allows to switch between active space orbital choices (please remember the remark in the section :ref:`chemps2_orbitalchoiceordering`): 
    * input orbitals (``0``) : No orbital rotations in addition to the occupied-active, active-virtual, and occupied-virtual rotations are performed.
    * natural orbitals (``1``) : After each DMRG calculation, and before the calculation of the gradient and Hessian for occupied-active, active-virtual, and occupied-virtual rotations, the active space orbitals are rotated to natural orbitals. This options allows to perform an initial approximate large active space DMRG calculation, from which a smaller active space can subsequently be selected based on the natural orbital occupation numbers.
    * localized (Edmiston-Ruedenberg) and ordered (Fiedler vector of the exchange matrix) orbitals (``2``) : Per irrep, the orbitals are localized with an augmented Hessian Newton-Raphson optimization of the Edmiston-Ruedenberg cost function. The orbitals are ordered per irrep based on the Fiedler vector of the exchange matrix. Orbitals belonging to different irreps are not mingled in the process.
    * ordered orbitals only (based on the Fiedler vector of the exchange matrix) (``3``) : The orbitals are ordered based on the Fiedler vector of the exchange matrix. Orbitals belonging to different irreps are mingled in the process.
* The variable ``DumpCorrelations_in`` allows to switch on printing the correlation functions defined in the section :ref:`chemps2_dmrg_object` after each DMRG calculation during the DMRG-SCF iterations.
* The variable ``StartLocRandom_in`` allows to start the localization of the orbitals (if ``WhichActiveSpace_in==2``) from a random orbital rotation. To study the :math:`\pi`-orbitals of polyenes, for example, the matrix elements should be generated in the :math:`\mathsf{C_s}` subgroup of the molecule's point group in order to be able to localize them. In order to break the symmetry of the R(O)HF orbitals during localization, it is important to start from a random orbital rotation.


.. _label-casscf-calculations-api:

``CheMPS2::CASSCF``
-------------------

Augmented Hessian Newton-Raphson DMRG-SCF calculations can be performed with the ``CheMPS2::CASSCF`` class:

.. code-block:: c++

    CheMPS2::CASSCF::CASSCF( CheMPS2::Hamiltonian * ham_in, int * docc, int * socc, int * nocc, int * ndmrg, int * nvirt, const string tmp_folder )
    double CheMPS2::CASSCF::solve( const int Nelectrons, const int TwoS, const int Irrep, ConvergenceScheme * OptScheme, const int rootNum, DMRGSCFoptions * scf_options )

The variables ``docc``, ``socc``, ``nocc``, ``ndmrg``, ``nvirt`` are arrays which have as length the number of irreps of the point group of the Hamiltonian (see the section :ref:`chemps2_psi4irrepconventions`). The first two arrays contain the number of doubly and singly occupied orbitals per irrep of the initial R(O)HF calculation. The last three arrays define for each of the irreps the number of occupied, active, and virtual orbitals, respectively.

The augmented Hessian Newton-Raphson DMRG-SCF calculation is started by calling the function ``solve``. The variables ``Nelectrons``, ``TwoS``, and ``Irrep`` define the active space symmetry sector. Note that ``Nelectrons`` is the total number of electrons, i.e. also the electrons in the doubly occupied core orbitals. ``TwoS`` is twice the targeted spin (multiplicity minus one). The numbering convention for the irreps can be found in the section :ref:`chemps2_psi4irrepconventions`. The variable ``rootNum`` defines how many states should be calculated during each DMRG calculation: ``rootNum==1`` means ground state only, ``rootNum==2`` means ground state and first excited state, etc. The DMRG instructions are passed in the ``CheMPS2::ConvergenceScheme`` object and the DMRG-SCF algorithmic choices in the ``CheMPS2::DMRGSCFoptions`` object. After completion, the function ``solve`` returns the DMRG-SCF energy.

For DMRG-SCF calculations, the number of reduced virtual basis states should not be decreased in the ``CheMPS2::ConvergenceScheme`` object. It is however advised to perform a few sweeps without noise and with very small residual norm tolerance (1e-8 to 1e-10) at the largest value of :math:`D_{\mathsf{SU(2)}}`. When you want to extrapolate the energy in the converged active space, it is better to create an orbital rotation checkpoint, and restart the DMRG-SCF calculation for one iteration (``MaxIterations_in==1``) with a different ``CheMPS2::ConvergenceScheme``.
   

