# Find Intel Math Kernel Library
#
# This sets the following variables:
#
#   MKL_INCLUDES   - Directory containing MKL include files
#   MKL_LIBRARIES  - Path to MKL libraries
#
# Users may wish to set:
#  MKL_USE_STATIC  - If ON search for static libraries
#  MKL_DIR         - Location to start searching for MKL libraries
#

set(MKL_DIR "${MKL_DIR}" CACHE PATH "MKL root diretory")

set(CMAKE_FIND_LIBRARY_SUFFIXES_ORG ${CMAKE_FIND_LIBRARY_SUFFIXES})
if(MKL_USE_STATIC)
	set(CMAKE_FIND_LIBRARY_SUFFIXES .a) # use static libs
endif()

find_path(MKL_INCLUDES NAMES mkl.h
	HINTS ${MKL_DIR} PATH_SUFFIXES include
)

# find architecture dependent libraries
if (${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	#64 bit
	set(MKL_PATH_SUFFIXES "lib" "lib/intel64" "lib/em64t" "64")
	find_library(MKL_INTEL_LIBRARY mkl_intel_lp64 HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	find_library(MKL_CORE_LIBRARY mkl_core HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
else()
	#32 bit
	set(MKL_PATH_SUFFIXES "lib/32" "32")
	find_library(MKL_INTEL_LIBRARY mkl_intel HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	find_library(MKL_CORE_LIBRARY mkl_core HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	#find_library(MKL_CORE mkl_solver HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
endif()

list(APPEND MKL_LIBRARIES ${MKL_INTEL_LIBRARY} ${MKL_CORE_LIBRARY})

# find threading libraries
if (PARALLEL_USE_OPENMP)
	if (UNIX AND NOT APPLE)
		if(CMAKE_C_COMPILER EQUAL "icc")
			set(MKL_USE_INTEL_THREAD ON)
		else()
			set(MKL_USE_GNU_THREAD ON)
		endif()
	else()
		set(MKL_USE_INTEL_THREAD ON)
	endif()
else()
	set(MKL_USE_SEQ_THREAD ON)
endif()

if(MKL_USE_INTEL_THREAD)
	find_library(MKL_INTEL_THREAD_LIBRARY mkl_intel_thread HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	find_library(MKL_OMP5_LIBRARY iomp5 HINTS "${MKL_DIR}/../lib" PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	find_library(PTHREAD_LIBRARY pthread)
	list(APPEND MKL_LIBRARIES ${MKL_INTEL_THREAD_LIBRARY} ${MKL_OMP5_LIBRARY} ${PTHREAD_LIBRARY})
elseif(MKL_USE_GNU_THREAD)
	set(MKL_USE_INTEL_THREAD OFF)
	find_library(MKL_GNU_THREAD_LIBRARY mkl_gnu_thread HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	list(APPEND MKL_LIBRARIES ${MKL_GNU_THREAD_LIBRARY})
elseif(MKL_USE_SEQ_THREAD)
	find_library(MKL_SEQUENTIAL_LIBRARY mkl_sequential HINTS ${MKL_DIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
	list(APPEND MKL_LIBRARIES ${MKL_SEQUENTIAL_LIBRARY})
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_ORG})

message (STATUS "MKL libraries found: ${MKL_LIBRARIES}")

#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDES MKL_LIBRARIES)
if(MKL_FOUND)
	mark_as_advanced(MKL_DIR MKL_INCLUDES MKL_LIBRARIES)
endif()
