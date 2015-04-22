find_package(Sphinx)
if (NOT SPHINX_FOUND)
    message (FATAL_ERROR "Sphinx is needed to build the user manual. Please install it correctly.")
endif (NOT SPHINX_FOUND)
set(SPHINX_DIR "${CheMPS2_BINARY_DIR}/sphinx")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in" "${SPHINX_DIR}/conf.py")
add_custom_target(sphinx ${SPHINX_EXECUTABLE} -q -b html
                                              -c "${SPHINX_DIR}"
                                              -d "${SPHINX_DIR}/.doctree"
                                              "${CMAKE_CURRENT_SOURCE_DIR}"
                                              "${SPHINX_DIR}/html")

