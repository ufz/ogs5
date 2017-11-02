# - Try to find freesteam
# Once done, this will define
#
#  FREESTEAM_FOUND
#  FREESTEAM_INCLUDE_DIRS
#  FREESTEAM_LIBRARIES
set(DEFAULT_FREESTEAM_INCLUDE_DIR /usr/local/)
find_path( FREESTEAM_INCLUDE_DIR
	NAMES steam.h
	PATHS ${DEFAULT_FREESTEAM_INCLUDE_DIR}/include/freesteam/)

if ( UNIX OR CYGWIN OR MINGW)
	find_library(FREESTEAM_LIBRARIES
			NAMES libfreesteam.so
		PATHS ${DEFAULT_FREESTEAM_INCLUDE_DIR}/lib)
else ()
	message(FATAL_ERROR  "freesteam is only available on UNIX/CYGWIN/MINGW")
endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FREESTEAM
	REQUIRED_VARS FREESTEAM_INCLUDE_DIR FREESTEAM_LIBRARIES
	FOUND_VAR FREESTEAM_FOUND
)
