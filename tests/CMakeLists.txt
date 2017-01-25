file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/tests)

set (TESTLIST "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9" "test10" "test11" "test12" "test13" "test14")

foreach (ITEM ${TESTLIST})
    configure_file (${CMAKE_SOURCE_DIR}/tests/${ITEM}.cpp.in ${CMAKE_BINARY_DIR}/tests/tests/${ITEM}.cpp)
    add_executable (${ITEM} ${CMAKE_BINARY_DIR}/tests/tests/${ITEM}.cpp)
    target_link_libraries (${ITEM} chemps2-lib)
    add_test (${ITEM} ${ITEM})
endforeach()

