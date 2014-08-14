SET (BUILD_SHARED_LIBS ON)

link_directories ($ENV{LD_LIBRARY_PATH})

if (MKL)
    set( ENV{BLA_VENDOR} "Intel10_64lp" )
endif (MKL)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CheMPS2_SOURCE_DIR}/CMake/")

find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)

include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/ ${HDF5_INCLUDE_DIRS})

set (CHEMPS2LIB_SOURCE_FILES "CASSCF.cpp" "CASSCFdebug.cpp" "CASSCFnewtonraphson.cpp" "ConvergenceScheme.cpp" "Correlations.cpp" "DIIS.cpp" "DMRG.cpp" "DMRGmpsio.cpp" "DMRGoperators.cpp" "DMRGSCFindices.cpp" "DMRGSCFunitary.cpp" "DMRGSCFVmatRotations.cpp" "DMRGtechnics.cpp" "FourIndex.cpp" "Hamiltonian.cpp" "Heff.cpp" "HeffDiagonal.cpp" "HeffDiagrams1.cpp" "HeffDiagrams2.cpp" "HeffDiagrams3.cpp" "HeffDiagrams4.cpp" "HeffDiagrams5.cpp" "Irreps.cpp" "PrintLicense.cpp" "Problem.cpp" "Sobject.cpp" "SyBookkeeper.cpp" "TensorA.cpp" "TensorB.cpp" "TensorC.cpp" "TensorD.cpp" "TensorDiag.cpp" "TensorF0Cbase.cpp" "TensorF0.cpp" "TensorF1.cpp" "TensorF1Dbase.cpp" "TensorGYZ.cpp" "TensorK.cpp" "TensorL.cpp" "TensorM.cpp" "TensorO.cpp" "TensorQ.cpp" "TensorS0Abase.cpp" "TensorS0.cpp" "TensorS1Bbase.cpp" "TensorS1.cpp" "TensorSwap.cpp" "TensorT.cpp" "TensorX.cpp" "TwoDM.cpp" "TwoDMstorage.cpp" "TwoIndex.cpp")

add_library (CheMPS2 ${CHEMPS2LIB_SOURCE_FILES})

target_link_libraries (CheMPS2 ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES})

