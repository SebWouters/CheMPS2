set(PN ${PROJECT_NAME})

include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2/ ${HDF5_INCLUDE_DIRS})

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

if (BUILD_FPIC OR NOT STATIC_ONLY)
    set_target_properties (chemps2-base PROPERTIES POSITION_INDEPENDENT_CODE 1)
endif()

if (NOT STATIC_ONLY)
    add_library (chemps2-shared SHARED $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-shared ${LIBC_INTERJECT} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})
    set_target_properties (chemps2-shared PROPERTIES SOVERSION ${CheMPS2_LIB_SOVERSION} OUTPUT_NAME "chemps2")
endif (NOT STATIC_ONLY)

if (NOT SHARED_ONLY)
    add_library (chemps2-static STATIC $<TARGET_OBJECTS:chemps2-base>)
    target_link_libraries (chemps2-static ${LIBC_INTERJECT} ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES})
    set_target_properties (chemps2-static PROPERTIES OUTPUT_NAME "chemps2")
endif (NOT SHARED_ONLY)

if (NOT STATIC_ONLY)
    install (TARGETS chemps2-shared
             EXPORT "${PN}Targets"
             LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif (NOT STATIC_ONLY)
if (NOT SHARED_ONLY)
    install (TARGETS chemps2-static
             EXPORT "${PN}Targets"
             ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
endif (NOT SHARED_ONLY)
install (DIRECTORY include/chemps2/ DESTINATION include/chemps2 FILES_MATCHING PATTERN "*.h")

if (NOT STATIC_ONLY)
    add_library(chemps2 ALIAS chemps2-shared)
else ()
    add_library(chemps2 ALIAS chemps2-static)
endif ()

add_executable (chemps2bin executable.cpp)
target_link_libraries (chemps2bin chemps2)
target_link_libraries (chemps2bin ${LIBC_INTERJECT})
target_link_libraries (chemps2bin ${LAPACK_LIBRARIES})
target_link_libraries (chemps2bin ${HDF5_LIBRARIES})
set_target_properties (chemps2bin PROPERTIES OUTPUT_NAME "chemps2")

install (TARGETS chemps2bin
         EXPORT "${PN}Targets"
         RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})

# <<<  Export Interface  >>>

#target_compile_definitions(chemps2-static INTERFACE USING_${PN})
#target_include_directories(efp INTERFACE
#                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)


# <<<  Export Config  >>>

#include(CMakePackageConfigHelpers)

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/${PN}")
install(EXPORT "${PN}Targets"
        NAMESPACE "${PN}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

