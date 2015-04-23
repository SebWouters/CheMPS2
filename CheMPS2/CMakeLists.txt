include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/ ${HDF5_INCLUDE_DIRS})

set (CHEMPS2LIB_SOURCE_FILES "CASSCF.cpp" "CASSCFdebug.cpp" "CASSCFnewtonraphson.cpp" "ConvergenceScheme.cpp" "Correlations.cpp" "Davidson.cpp" "DIIS.cpp" "DMRG.cpp" "DMRGmpsio.cpp" "DMRGoperators.cpp" "DMRGSCFindices.cpp" "DMRGSCFintegrals.cpp" "DMRGSCFmatrix.cpp" "DMRGSCFoptions.cpp" "DMRGSCFunitary.cpp" "DMRGSCFVmatRotations.cpp" "DMRGSCFwtilde.cpp" "DMRGtechnics.cpp" "EdmistonRuedenberg.cpp" "FCI.cpp" "FourIndex.cpp" "Hamiltonian.cpp" "Heff.cpp" "HeffDiagonal.cpp" "HeffDiagrams1.cpp" "HeffDiagrams2.cpp" "HeffDiagrams3.cpp" "HeffDiagrams4.cpp" "HeffDiagrams5.cpp" "Initialize.cpp" "Irreps.cpp" "PrintLicense.cpp" "Problem.cpp" "Sobject.cpp" "SyBookkeeper.cpp" "TensorA.cpp" "TensorB.cpp" "TensorC.cpp" "TensorD.cpp" "TensorDiag.cpp" "TensorF0Cbase.cpp" "TensorF0.cpp" "TensorF1.cpp" "TensorF1Dbase.cpp" "TensorGYZ.cpp" "TensorK.cpp" "TensorL.cpp" "TensorM.cpp" "TensorO.cpp" "TensorQ.cpp" "TensorS0Abase.cpp" "TensorS0.cpp" "TensorS1Bbase.cpp" "TensorS1.cpp" "TensorSwap.cpp" "TensorT.cpp" "TensorX.cpp" "TwoDM.cpp" "TwoDMstorage.cpp" "TwoIndex.cpp")

add_library (chemps2-base OBJECT ${CHEMPS2LIB_SOURCE_FILES})

if (NOT STATIC_ONLY)
    set_target_properties (chemps2-base PROPERTIES POSITION_INDEPENDENT_CODE 1)
    add_library (chemps2-shared SHARED $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-shared ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES})
    set_target_properties (chemps2-shared PROPERTIES SOVERSION ${CheMPS2_LIB_SOVERSION} CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME "chemps2")
endif (NOT STATIC_ONLY)

add_library (chemps2-static STATIC $<TARGET_OBJECTS:chemps2-base>)
target_link_libraries (chemps2-static ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES})
set_target_properties (chemps2-static PROPERTIES CLEAN_DIRECT_OUTPUT 1 OUTPUT_NAME "chemps2")

if (NOT STATIC_ONLY)
    install (TARGETS chemps2-shared DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif (NOT STATIC_ONLY)
install (TARGETS chemps2-static DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
install (DIRECTORY include/ DESTINATION ${CheMPS2_INC_RELPATH} FILES_MATCHING PATTERN "*.h")

