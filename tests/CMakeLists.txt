include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2/ ${HDF5_INCLUDE_DIRS})

link_directories (${CMAKE_BINARY_DIR}/CheMPS2)

file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/tests)

if (WITH_MPI)
    set (TESTLIST "test1" "test2" "test3" "test4" "test5" "test7" "test9" "test10" "test11" "test12")
else (WITH_MPI)
    set (TESTLIST "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9" "test10" "test11" "test12" "test13" "test14")
endif (WITH_MPI)

foreach (ITEM ${TESTLIST})
    configure_file (${CMAKE_SOURCE_DIR}/tests/${ITEM}.cpp.in ${CMAKE_BINARY_DIR}/tests/tests/${ITEM}.cpp)
    add_executable (${ITEM} ${CMAKE_BINARY_DIR}/tests/tests/${ITEM}.cpp)
    if (STATIC_ONLY)
        add_dependencies (${ITEM} chemps2-static)
        target_link_libraries (${ITEM} chemps2 ${LAPACK_LIBRARIES} ${HDF5_LIBRARIES} ${GSL_LIBRARIES})
    else (STATIC_ONLY)
        add_dependencies (${ITEM} chemps2-shared)
        target_link_libraries (${ITEM} chemps2)
    endif (STATIC_ONLY)
    add_test(${ITEM} ${ITEM})
endforeach (ITEM)

