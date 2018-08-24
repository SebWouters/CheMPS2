# FindTargetLAPACK.cmake
# ----------------------
#
# LAPACK cmake module to wrap FindLAPACK.cmake in a target.
#
# This module sets the following variables in your project: ::
#
#   TargetLAPACK_FOUND - true if BLAS/LAPACK found on the system
#   TargetLAPACK_MESSAGE - status message with BLAS/LAPACK library path list
#
# This module *unsets* the following conventional LAPACK variables so as
#   to force using the target: ::
#
#   LAPACK_FOUND
#   LAPACK_LIBRARIES
#
# In order of decreasing precedence, this module returns in a target ``tgt::lapack``
#  (1) the libraries passed through CMake variable LAPACK_LIBRARIES,
#  (2) the libraries defined in a detectable TargetLAPACKConfig.cmake file
#      (skip via DISABLE_FIND_PACKAGE_TargetLAPACK), or
#  (3) the libraries detected by the usual FindLAPACK.cmake module.
#

set(PN TargetLAPACK)

# 1st precedence - libraries passed in through -DLAPACK_LIBRARIES
if (LAPACK_LIBRARIES)
    if (NOT ${PN}_FIND_QUIETLY)
        message (STATUS "LAPACK detection suppressed.")
    endif()

    add_library (tgt::lapack INTERFACE IMPORTED)
    set_property (TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
else()
    # 2nd precedence - target already prepared and findable in TargetLAPACKConfig.cmake
    if (NOT "${DISABLE_FIND_PACKAGE_${PN}}")
        find_package (TargetLAPACK QUIET CONFIG)
    endif()
    if (TARGET tgt::lapack)
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "TargetLAPACKConfig detected.")
        endif()
    else()
        # 3rd precedence - usual variables from FindLAPACK.cmake
        find_package (LAPACK QUIET MODULE)
        if (NOT ${PN}_FIND_QUIETLY)
            message (STATUS "LAPACK detected.")
        endif()
    
        add_library (tgt::lapack INTERFACE IMPORTED)
        set_property (TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})

        unset (LAPACK_FOUND)
        unset (LAPACK_LIBRARIES)
    endif()
endif()    

get_property (_ill TARGET tgt::lapack PROPERTY INTERFACE_LINK_LIBRARIES)
set (${PN}_MESSAGE "Found LAPACK: ${_ill}")
if ((TARGET tgt::blas) AND (TARGET tgt::lapk))
    get_property (_illb TARGET tgt::blas PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property (_illl TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
    set (${PN}_MESSAGE "Found LAPACK: ${_illl};${_illb}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (${PN} DEFAULT_MSG ${PN}_MESSAGE)
