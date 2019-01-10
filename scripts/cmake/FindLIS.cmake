# - Try to find LIS
# Once done, this will define
#
#  LIS_FOUND
#  LIS_INCLUDE_DIRS
#  LIS_LIBRARIES

find_path( LIS_INCLUDE_DIR
    NAMES lis.h
    PATHS ${OGS_PRECOMPILED_LIBS_DIR_FOUND})

if ( UNIX )
    # Tell if the unix system is on 64-bit base
    if(CMAKE_SIZEOF_VOID_P MATCHES "8")
        find_library(LIS_LIBRARIES
            NAMES lis-64
            PATHS ${OGS_PRECOMPILED_LIBS_DIR_FOUND}
        )
    else ()
        find_library(LIS_LIBRARIES
            NAMES lis-32
            PATHS ${OGS_PRECOMPILED_LIBS_DIR_FOUND}
        )
    endif ()
else ()
    find_library(LIS_LIBRARIES
        NAMES lisomp
        PATHS ${OGS_PRECOMPILED_LIBS_DIR_FOUND}
    )
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIS
    REQUIRED_VARS LIS_INCLUDE_DIR LIS_LIBRARIES
    FOUND_VAR LIS_FOUND
)
