SET (BUILD_SHARED_LIBS ON)

include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/)

link_directories ($ENV{LD_LIBRARY_PATH})

add_library (CheMPS2 CASSCF.cpp CASSCFdebug.cpp CASSCFhamiltonianrotation.cpp CASSCFnewtonraphson.cpp ConvergenceScheme.cpp DMRG.cpp DMRGmpsio.cpp DMRGoperators.cpp DMRGtechnics.cpp FourIndex.cpp Hamiltonian.cpp Heff.cpp HeffDiagonal.cpp HeffDiagrams1.cpp HeffDiagrams2.cpp HeffDiagrams3.cpp HeffDiagrams4.cpp HeffDiagrams5.cpp Irreps.cpp PrintLicense.cpp Problem.cpp Sobject.cpp SyBookkeeper.cpp TensorA.cpp TensorB.cpp TensorC.cpp TensorD.cpp TensorDiag.cpp TensorF0Cbase.cpp TensorF0.cpp TensorF1.cpp TensorF1Dbase.cpp TensorL.cpp TensorO.cpp TensorQ.cpp TensorS0Abase.cpp TensorS0.cpp TensorS1Bbase.cpp TensorS1.cpp TensorSwap.cpp TensorT.cpp TensorX.cpp TwoDM.cpp TwoIndex.cpp)

if (INTELCOMPILER)
   target_link_libraries (CheMPS2 mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread hdf5 gsl gslcblas)
else (INTELCOMPILER)
   target_link_libraries (CheMPS2 blas lapack hdf5 gsl gslcblas)
endif (INTELCOMPILER)


