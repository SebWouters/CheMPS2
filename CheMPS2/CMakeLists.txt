cmake_policy (SET CMP0022 NEW)
cmake_policy (SET CMP0028 NEW)

set (CHEMPS2LIB_SOURCE_FILES "CASPT2.cpp"
                             "CASSCF.cpp"
                             "CASSCFdebug.cpp"
                             "CASSCFnewtonraphson.cpp"
                             "CASSCFpt2.cpp"
                             "ConjugateGradient.cpp"
                             "ConvergenceScheme.cpp"
                             "Correlations.cpp"
                             "Cumulant.cpp"
                             "Davidson.cpp"
                             "DIIS.cpp"
                             "DMRG.cpp"
                             "DMRGfock.cpp"
                             "DMRGmpsio.cpp"
                             "DMRGoperators.cpp"
                             "DMRGoperators3RDM.cpp"
                             "DMRGSCFindices.cpp"
                             "DMRGSCFintegrals.cpp"
                             "DMRGSCFmatrix.cpp"
                             "DMRGSCFoptions.cpp"
                             "DMRGSCFrotations.cpp"
                             "DMRGSCFunitary.cpp"
                             "DMRGSCFwtilde.cpp"
                             "DMRGtechnics.cpp"
                             "EdmistonRuedenberg.cpp"
                             "Excitation.cpp"
                             "FCI.cpp"
                             "FourIndex.cpp"
                             "Hamiltonian.cpp"
                             "Heff.cpp"
                             "HeffDiagonal.cpp"
                             "HeffDiagrams1.cpp"
                             "HeffDiagrams2.cpp"
                             "HeffDiagrams3.cpp"
                             "HeffDiagrams4.cpp"
                             "HeffDiagrams5.cpp"
                             "Initialize.cpp"
                             "Irreps.cpp"
                             "Molden.cpp"
                             "PrintLicense.cpp"
                             "Problem.cpp"
                             "Sobject.cpp"
                             "SyBookkeeper.cpp"
                             "Tensor3RDM.cpp"
                             "TensorF0.cpp"
                             "TensorF1.cpp"
                             "TensorGYZ.cpp"
                             "TensorKM.cpp"
                             "TensorL.cpp"
                             "TensorO.cpp"
                             "TensorOperator.cpp"
                             "TensorQ.cpp"
                             "TensorS0.cpp"
                             "TensorS1.cpp"
                             "TensorT.cpp"
                             "TensorX.cpp"
                             "ThreeDM.cpp"
                             "TwoDM.cpp"
                             "TwoIndex.cpp"
                             "Wigner.cpp")

add_library (chemps2-base OBJECT ${CHEMPS2LIB_SOURCE_FILES})
target_include_directories (chemps2-base PRIVATE ${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2
                                                 $<TARGET_PROPERTY:tgt::hdf5,INTERFACE_INCLUDE_DIRECTORIES>)

if (BUILD_FPIC OR NOT STATIC_ONLY)
    set_target_properties (chemps2-base PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()

if (NOT STATIC_ONLY)
    add_library           (chemps2-shared SHARED  $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-shared PRIVATE ${LIBC_INTERJECT})
    target_link_libraries (chemps2-shared PUBLIC  tgt::lapack
                                                  tgt::hdf5)
    set_target_properties (chemps2-shared PROPERTIES SOVERSION ${CheMPS2_LIB_SOVERSION}
                                                     MACOSX_RPATH ON
                                                     OUTPUT_NAME "chemps2"
                                                     EXPORT_NAME "chemps2")
endif()

if (NOT SHARED_ONLY)
    add_library           (chemps2-static STATIC  $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-static PRIVATE ${LIBC_INTERJECT})
    target_link_libraries (chemps2-static PUBLIC  tgt::lapack
                                                  tgt::hdf5)
    set_target_properties (chemps2-static PROPERTIES OUTPUT_NAME "chemps2"
                                                     EXPORT_NAME "chemps2")
endif()

if (STATIC_ONLY)
    add_library (chemps2-lib ALIAS chemps2-static)
else ()
    add_library (chemps2-lib ALIAS chemps2-shared)
endif()

add_executable             (chemps2-bin executable.cpp)
target_link_libraries      (chemps2-bin chemps2-lib ${LIBC_INTERJECT})
set_target_properties      (chemps2-bin PROPERTIES OUTPUT_NAME "chemps2")

# <<<  Install  >>>

if (NOT STATIC_ONLY)
    install (TARGETS chemps2-shared chemps2-bin
                     EXPORT "${PROJECT_NAME}Targets-shared"
                     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if (NOT SHARED_ONLY)
    install (TARGETS chemps2-static chemps2-bin
                     EXPORT "${PROJECT_NAME}Targets-static"
                     ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                     RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install (DIRECTORY include/chemps2/
                   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/chemps2
                   FILES_MATCHING PATTERN "*.h")

# <<<  Export interface  >>>

if (NOT STATIC_ONLY)
    if (APPLE)
        set_target_properties(chemps2-shared PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif()
    target_compile_definitions (chemps2-shared INTERFACE USING_${PROJECT_NAME})
    target_include_directories (chemps2-shared INTERFACE $<BUILD_INTERFACE:${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2>
                                                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                                                         $<TARGET_PROPERTY:tgt::hdf5,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

if (NOT SHARED_ONLY)
    target_compile_definitions (chemps2-static INTERFACE USING_${PROJECT_NAME})
    target_include_directories (chemps2-static INTERFACE $<BUILD_INTERFACE:${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2>
                                                         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                                                         $<TARGET_PROPERTY:tgt::hdf5,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

# <<<  Export Config  >>>

    # explicit "share" not "DATADIR" for CMake search path
set (CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")
if (NOT STATIC_ONLY)
    install (EXPORT "${PROJECT_NAME}Targets-shared"
                     NAMESPACE "${PROJECT_NAME}::"
                     DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

if (NOT SHARED_ONLY)
    install (EXPORT "${PROJECT_NAME}Targets-static"
                     NAMESPACE "${PROJECT_NAME}::"
                     DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

