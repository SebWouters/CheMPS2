set (PyCheMPS2_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/PyCheMPS2")

configure_file (${CMAKE_SOURCE_DIR}/PyCheMPS2/setup.py.in ${CMAKE_BINARY_DIR}/setup.py)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN "*.pyx")
file (COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILES_MATCHING PATTERN "*.pxd")

file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/PyTests)
set (PYCHEMPS2LIST "ReadinHamiltonianPsi4" "test1" "test2" "test3" "test4" "test5" "test6" "test7" "test8" "test9")

foreach (ITEM ${PYCHEMPS2LIST})
    configure_file (${CMAKE_SOURCE_DIR}/PyCheMPS2/${ITEM}.py.in ${CMAKE_BINARY_DIR}/PyTests/${ITEM}.py)
endforeach (ITEM)

