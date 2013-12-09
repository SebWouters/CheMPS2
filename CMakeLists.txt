cmake_minimum_required (VERSION 2.6)
project (CheMPS2)

option (INTELCOMPILER "Compile with the intel compiler" ON)
option (BUILD_DOCUMENTATION "Use Doxygen to create a PDF manual" ON)
set (CMAKE_VERBOSE_MAKEFILE OFF)

if (INTELCOMPILER)
   set (CMAKE_CXX_COMPILER icpc)
   set (CMAKE_CXX_FLAGS "-g -Wall -openmp -xHost -ipo")
else (INTELCOMPILER)
   set (CMAKE_CXX_COMPILER g++)
   set (CMAKE_CXX_FLAGS "-g -Wall -fopenmp -O3")
endif (INTELCOMPILER)

add_subdirectory (CheMPS2)
add_subdirectory (tests)

if (BUILD_DOCUMENTATION)
   find_package (Doxygen)
   if (NOT DOXYGEN_FOUND)
      message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
   endif (NOT DOXYGEN_FOUND)
   configure_file (${CheMPS2_SOURCE_DIR}/CheMPS2/Doxyfile.in ${CheMPS2_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
   add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE} ${CheMPS2_BINARY_DIR}/Doxyfile SOURCES ${CheMPS2_BINARY_DIR}/Doxyfile)
endif (BUILD_DOCUMENTATION)


