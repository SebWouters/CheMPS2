project (CheMPS2)
set (CheMPS2_MMP_VERSION   "1.8.10")
set (CheMPS2_VERSION       "${CheMPS2_MMP_VERSION}-3 (2021-02-01)")
#set (CheMPS2_VERSION       "${CheMPS2_MMP_VERSION}-15 (2017-05-22)")
set (CheMPS2_AUTHORS       "Sebastian Wouters")
set (CheMPS2_DESCRIPTION   "A spin-adapted implementation of DMRG for ab initio quantum chemistry")
set (CheMPS2_EMAIL         "sebastianwouters@gmail.com")
set (CheMPS2_URL           "https://github.com/SebWouters/CheMPS2")
set (CheMPS2_LICENSE       "GNU General Public License, version 2")
set (CheMPS2_LIB_SOVERSION 3)

cmake_minimum_required (VERSION 3.0.2)

# <<<  CMake includes  >>>

include (GNUInstallDirs)
include (CheckCXXCompilerFlag)
include (CMakePackageConfigHelpers)

# <<<  Default CMake options  >>>

option (MKL                  "Compile using the MKL"                   OFF)
option (BUILD_DOXYGEN        "Use Doxygen to create a HTML/PDF manual" OFF)
option (BUILD_SPHINX         "Build the user manual with Sphinx"       OFF)
option (STATIC_ONLY          "Compile only the static library"         OFF)
option (SHARED_ONLY          "Compile only the shared library"         OFF)
option (ENABLE_TESTS         "Compile the tests"                       ON)
option (ENABLE_XHOST         "Enable processor-specific optimizations" ON)
option (ENABLE_GENERIC       "Enable mostly static linking in shared library" OFF)
option (ENABLE_OPENMP        "Enable OpenMP parallelization"           ON)
option (WITH_MPI             "Build the library with MPI"              OFF)
option (BUILD_FPIC           "Static library in STATIC_ONLY will be compiled with position independent code" OFF)

set (CMAKE_VERBOSE_MAKEFILE  OFF)
set (CMAKE_MODULE_PATH       ${CMAKE_MODULE_PATH} "${CheMPS2_SOURCE_DIR}/CMake/")

if (STATIC_ONLY AND SHARED_ONLY)
    message (FATAL_ERROR "The options STATIC_ONLY=ON and SHARED_ONLY=ON are conflicting." )
endif()

# <<<  Build with HDF5 1.10 API: CheMPS2 does not support HDF5 1.12 API  >>>
add_definitions(-DH5_USE_110_API)

# <<<  Enable MPI  >>>

if (WITH_MPI)
    add_definitions (-DCHEMPS2_MPI_COMPILATION)
endif ()

# <<< Pass CheMPS2 version  >>>

add_definitions(-DCHEMPS2_VERSION="${CheMPS2_VERSION}")

# <<<  Check build type  >>>

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
message("-- The CMake build type is ${CMAKE_BUILD_TYPE}")

# <<<  Enable OpenMP  >>>

if (ENABLE_OPENMP)
    find_package (OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# <<<  Enable host specific optimizations  >>>

if (ENABLE_XHOST)
    check_cxx_compiler_flag (-xHost HAS_XHOST)
    check_cxx_compiler_flag (-march=native HAS_MARCH_NATIVE)
    if (HAS_XHOST)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")
    elseif (HAS_MARCH_NATIVE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
endif()

# <<<  Enable inter-parts / link-time optimization  >>>

check_cxx_compiler_flag (-ipo HAS_IPO)
if (HAS_IPO)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo")
endif()

check_cxx_compiler_flag (-flto HAS_FLTO)
if (HAS_FLTO)
    set (CMAKE_CXX_FLAGS "-flto ${CMAKE_CXX_FLAGS}")
endif()

# <<<  Enable generic  >>>

if (ENABLE_GENERIC)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")
    check_cxx_compiler_flag (-static-intel HAS_INTEL_COMPILERS)
    if (HAS_INTEL_COMPILERS)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-intel -wd10237")
    endif()
endif()

# <<<  Find LAPACK and BLAS  >>>

link_directories ($ENV{LD_LIBRARY_PATH})
if (MKL)
    set (ENV{BLA_VENDOR} "Intel10_64lp")
endif()

find_package (TargetLAPACK REQUIRED)

get_property(_lapack_libraries TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES)
foreach (_l IN LISTS _lapack_libraries)
    if (${_l} MATCHES "mkl")
        set (MKL "ON")
    endif()
endforeach()

if (MKL)
    add_definitions (-DCHEMPS2_MKL)
endif ()

# <<<  Find HDF5  >>>

find_package (TargetHDF5 REQUIRED)

# <<<  Add source files, tests, sphinx, and doxygen  >>>

enable_testing ()
add_subdirectory (CheMPS2)
if (ENABLE_TESTS)
    add_subdirectory (tests)
endif()

if (BUILD_SPHINX)
    add_subdirectory (sphinx)
endif()

if (BUILD_DOXYGEN)
    find_package (Doxygen)
    if (NOT DOXYGEN_FOUND)
        message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
    endif()
    configure_file (${CheMPS2_SOURCE_DIR}/CheMPS2/Doxyfile.in ${CheMPS2_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
    add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE} ${CheMPS2_BINARY_DIR}/Doxyfile SOURCES ${CheMPS2_BINARY_DIR}/Doxyfile)
endif()

# <<<  Export config  >>>

# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set (CMAKECONFIG_INSTALL_DIR "share/cmake/${PROJECT_NAME}")
configure_package_config_file (${CheMPS2_SOURCE_DIR}/CMake/${PROJECT_NAME}Config.cmake.in
                                    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                                    INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
write_basic_package_version_file (${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                      VERSION ${CheMPS2_MMP_VERSION}
                                      COMPATIBILITY AnyNewerVersion)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
               ${CheMPS2_SOURCE_DIR}/CMake/FindTargetHDF5.cmake
               DESTINATION ${CMAKECONFIG_INSTALL_DIR})

