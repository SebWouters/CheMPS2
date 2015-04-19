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
    * input orbitals (``0``)
    * natural orbitals (``1``)
    * localized (Edmiston-Ruedenberg) and ordered (Fiedler vector of the exchange matrix) orbitals (``2``)
* The variable ``DumpCorrelations_in`` allows to switch on printing the correlation functions defined in the section :ref:`chemps2_dmrg_object` after each DMRG calculation during the DMRG-SCF iterations.
* The variable ``StartLocRandom_in`` allows to start the localization of the orbitals (if ``WhichActiveSpace_in==2``) from a random orbital rotation. To study the :math:`\pi`-orbitals of polyenes, for example, the matrix elements should be generated in the :math:`\mathsf{Cs}` subgroup of the molecule's point group in order to be able to localize them. In order to break the symmetry of the R(O)HF orbitals during localization, it is important to start from a random orbital rotation.


``CheMPS2::CASSCF``
-------------------

Augmented Hessian Newton-Raphson DMRG-SCF calculations can be performed with the ``CheMPS2::CASSCF`` class:

.. code-block:: c++

    CheMPS2::CASSCF::CASSCF( CheMPS2::Hamiltonian * HamIn, int * DOCCin, int * SOCCin )
    void CheMPS2::CASSCF::setupStart( int * NoccIn, int * NDMRGIn, int * NvirtIn )
    double CheMPS2::CASSCF::doCASSCFnewtonraphson( const int Nelectrons, const int TwoS, const int Irrep, CheMPS2::ConvergenceScheme * OptScheme, const int rootNum, CheMPS2::DMRGSCFoptions * theDMRGSCFoptions )

The variables ``DOCCin``, ``SOCCin``, ``NoccIn``, ``NDMRGIn``, ``NvirtIn`` are arrays which have as length the number of irreps of the point group of the Hamiltonian (see the section :ref:`chemps2_psi4irrepconventions`). The first two arrays contain the number of doubly and singly occupied orbitals per irrep of the initial R(O)HF calculation. The last three arrays define for each of the irreps the number of occupied, active, and virtual orbitals, respectively.

The three functions above should be called in the order they are displayed. The augmented Hessian Newton-Raphson DMRG-SCF calculation is started by calling the function ``doCASSCFnewtonraphson``. The variables ``Nelectrons``, ``TwoS``, and ``Irrep`` define the active space symmetry sector. Note that ``TwoS`` is twice the targeted spin (multiplicity minus one). The numbering convention for the irreps can be found in the section :ref:`chemps2_psi4irrepconventions`. The variable ``rootNum`` defines how many states should be calculated during each DMRG calculation: ``rootNum==1`` means ground state only, ``rootNum==2`` means ground state and first excited state, etc. The DMRG instructions are passed in the ``CheMPS2::ConvergenceScheme`` object and the DMRG-SCF algorithmic choices in the ``CheMPS2::DMRGSCFoptions`` object. After completion, the function ``doCASSCFnewtonraphson`` returns the DMRG-SCF energy.

For DMRG-SCF calculations, the number of reduced virtual basis states should not be descreased in the ``CheMPS2::ConvergenceScheme`` object. It is however advised to perform a few sweeps without noise at the largest value of :math:`D_{\mathsf{SU(2)}}`. When you want to extrapolate the energy in the converged active space, it is better to create an orbital rotation checkpoint, and restart the DMRG-SCF calculation for one iteration (``MaxIterations_in==1``).
   

