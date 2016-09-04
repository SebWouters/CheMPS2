cmake_policy(SET CMP0022 NEW)
set(PN ${PROJECT_NAME})

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
target_include_directories(chemps2-base PRIVATE ${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2 ${HDF5_INCLUDE_DIRS})

if (BUILD_FPIC OR NOT STATIC_ONLY)
    set_target_properties (chemps2-base PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()

if (NOT STATIC_ONLY)
    add_library (chemps2-shared SHARED $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-shared PRIVATE ${LIBC_INTERJECT})
    set_target_properties (chemps2-shared PROPERTIES SOVERSION ${CheMPS2_LIB_SOVERSION}
                                                     OUTPUT_NAME "chemps2"
                                                     EXPORT_NAME "chemps2")
endif ()

if (NOT SHARED_ONLY)
    add_library (chemps2-static STATIC $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-static PRIVATE ${LIBC_INTERJECT})
    set_target_properties (chemps2-static PROPERTIES OUTPUT_NAME "chemps2"
                                                     EXPORT_NAME "chemps2")
endif ()

if (STATIC_ONLY)
    add_library (chemps2 ALIAS chemps2-static)
else ()
    add_library (chemps2 ALIAS chemps2-shared)
endif ()

add_executable (chemps2bin executable.cpp)
target_link_libraries (chemps2bin chemps2 ${LIBC_INTERJECT})
set_target_properties (chemps2bin PROPERTIES OUTPUT_NAME "chemps2")
target_include_directories(chemps2bin PRIVATE ${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2 ${HDF5_INCLUDE_DIRS})

# <<<  Install  >>>

if (NOT STATIC_ONLY)
    install (TARGETS chemps2-shared chemps2bin
             EXPORT "${PN}Targets-shared"
             LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (NOT SHARED_ONLY)
    install (TARGETS chemps2-static chemps2bin
             EXPORT "${PN}Targets-static"
             ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
             RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

install (DIRECTORY include/chemps2/
         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/chemps2
         FILES_MATCHING PATTERN "*.h")

# <<<  Export Interface  >>>

if (NOT STATIC_ONLY)
    target_compile_definitions(chemps2-shared INTERFACE USING_${PN})
    target_link_libraries     (chemps2-shared INTERFACE ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})
    target_include_directories(chemps2-shared INTERFACE
                               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                               $<INSTALL_INTERFACE:${HDF5_INCLUDE_DIRS}>)
endif ()

if (NOT SHARED_ONLY)
    target_compile_definitions(chemps2-static INTERFACE USING_${PN})
    target_link_libraries     (chemps2-static INTERFACE ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})
    target_include_directories(chemps2-static INTERFACE
                               $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                               $<INSTALL_INTERFACE:${HDF5_INCLUDE_DIRS}>)
endif ()

# <<<  Export Config  >>>

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/${PN}")
if (NOT STATIC_ONLY)
    install(EXPORT "${PN}Targets-shared"
            NAMESPACE "${PN}::"
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()

if (NOT SHARED_ONLY)
    install(EXPORT "${PN}Targets-static"
            NAMESPACE "${PN}::"
            DESTINATION ${CMAKECONFIG_INSTALL_DIR})
endif()
