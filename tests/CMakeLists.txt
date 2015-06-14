include_directories (${CheMPS2_SOURCE_DIR}/CheMPS2/include/chemps2/
                     ${HDF5_INCLUDE_DIRS}
                     ${HDF5_INCLUDE_DIRS}/hdf5/serial/)
# Last line to accommodate for the libhdf5-dev changes https://packages.debian.org/sid/amd64/libhdf5-dev/filelist until CMake catches up

link_directories (${CMAKE_BINARY_DIR}/CheMPS2)

file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/tests)

set (TESTLIST "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9" "test10")

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

