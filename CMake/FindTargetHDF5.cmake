# FindTargetHDF5.cmake
# --------------------
#
# HDF5 cmake module to wrap FindHDF5.cmake in a target.
#
# This module sets the following variables in your project: ::
#
#   TargetHDF5_FOUND - true if HDF5 and all required components found on the system
#   TargetHDF5_VERSION - HDF5 version in format Major.Minor.Release
#   TargetHDF5_MESSAGE - status message with HDF5 library path list and version
#
# Note that components are passed along to find_package(HDF5 (untested) but not checked in the direct TargetHDF5Config
# Note that version checking/attaching not working yet
#
# This module *unsets* the following conventional HDF5 variables so as
#   to force using the target: ::
#
#   HDF5_FOUND
#   HDF5_VERSION
#   HDF5_INCLUDE_DIRS
#   HDF5_LIBRARIES
#
# Exported targets::
#
# If HDF5 is found, this module defines the following :prop_tgt:`IMPORTED`
# target. ::
#
#   tgt::hdf5 - the HDF5 libraries with headers attached.
#
# Suggested usage::
#
#   find_package(TargetHDF5)
#   find_package(TargetHDF5 1.8.16 REQUIRED)
#
#
# The following variables can be set to guide the search for this package::
#
#   TargetHDF5_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package

set(PN TargetHDF5)

# 1st precedence - libraries passed in through -DHDF5_LIBRARIES
if (HDF5_LIBRARIES AND HDF5_INCLUDE_DIRS)
    if (HDF5_VERSION)
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "HDF5 detection suppressed.")
        endif()

        add_library (tgt::hdf5 INTERFACE IMPORTED)
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
        set (${PN}_VERSION ${HDF5_VERSION})
    else()
        message (FATAL_ERROR "Humor the build system - pass in the version, too (for example, -DHDF5_VERSION=1.8.17).")
    endif()
else()
    # 2nd precedence - target already prepared and findable in TargetHDF5Config.cmake
    find_package (TargetHDF5 QUIET CONFIG)
    if ((TARGET tgt::hdf5) AND (${PN}_VERSION))
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "TargetHDF5Config detected.")
        endif()
    else()
        # 3rd precedence - usual variables from FindHDF5.cmake
        find_package (HDF5 QUIET COMPONENTS ${HDF5_FIND_COMPONENTS})
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "HDF5 detected.")
        endif()
    
        add_library (tgt::hdf5 INTERFACE IMPORTED)
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES ${HDF5_LIBRARIES})
        set_property (TARGET tgt::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
        set (${PN}_VERSION ${HDF5_VERSION})

        unset (HDF5_FOUND)
        unset (HDF5_VERSION)
        unset (HDF5_LIBRARIES)
        unset (HDF5_INCLUDE_DIRS)
    endif()
endif()    

get_property(_ill TARGET tgt::hdf5 PROPERTY INTERFACE_LINK_LIBRARIES)
get_property(_iid TARGET tgt::hdf5 PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
set(${PN}_MESSAGE "Found HDF5: ${_ill} (found version ${${PN}_VERSION})")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(${PN}
                                  REQUIRED_VARS ${PN}_MESSAGE
                                  VERSION_VAR ${PN}_VERSION)
