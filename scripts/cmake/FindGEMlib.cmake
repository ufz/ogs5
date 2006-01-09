# - Try to find lib for the Gibbs Energy Minimization (GEM)
# Once done, this will define
#
#  GEMlib_FOUND
#  GEMlib_INCLUDE_DIRS
#  GEMlib_LIBRARIES

if (NOT GEMlib_FOUND)

	include(LibFindMacros)

	find_path( GEMlib_INCLUDE_DIR
		NAMES node.h
		PATHS ${CMAKE_SOURCE_DIR}/LIB/ )

	if ( UNIX )
		find_library(GEMlib_LIBRARIES
			NAMES libgemipm2k
			PATHS ${CMAKE_SOURCE_DIR}/LIB/ )
	else ()
		find_library(GEMlib_LIBRARIES
			NAMES GEMS3_rl
			PATHS ${CMAKE_SOURCE_DIR}/LIB/ )
		find_library(GEMlib_LIBRARIES
			NAMES qd
			PATHS ${CMAKE_SOURCE_DIR}/LIB/ )
	endif ()

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
	if (NOT GEMlib_LIBRARIES STREQUAL "GEMlib_LIBRARIES-NOTFOUND" AND NOT GEMlib_INCLUDE_DIR STREQUAL "GEMlib_INCLUDE_DIR-NOTFOUND")
		set(GEMlib_PROCESS_INCLUDES GEMlib_INCLUDE_DIR)
		set(GEMlib_PROCESS_LIBS GEMlib_LIBRARIES)
		libfind_process(GEMlib)
	else ()
		message (STATUS "Warning: GEMlib not found!")
	endif ()

endif ()
