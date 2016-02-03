#### OGS-5 Project ####

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.8.3)

# Project name
project(OGS)

###############
### Options ###
###############
# Configurations
set(OGS_CONFIG "FEM" CACHE STRING "OpenGeoSys configuration")
set_property(CACHE OGS_CONFIG PROPERTY STRINGS
	FEM # Default
	SP
	GEMS
	JFNK
	PETSC
	PETSC_GEMS
	IPQC
	PQC
	BRNS
	LIS
	CAP
	MKL
	MPI
)

# Options
option(BLUE_G "Blue/G optimization" OFF)
option(OGS_USE_CVODE "Use sundials cvode for TNEQ/TES" OFF)
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, defaults to native.")
option(OGS_PYTHON "Enable python scripting interface" OFF)
option(OGS_BUILD_UTILITIES "Build additional utility programs?" OFF)
option(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
option(OGS_ONE_BIG_GTEST "Should all gtests be collected to one ctest" ON)
option(OGS_BUILD_INFO "Should build information be generated" ON)
option (OGS_BUILD_TESTS "Enables building of tests." OFF)
option (OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)
option (OGS_DELETE_EDGES_AFTER_INIT "Delete mesh edges after initialization if possible" OFF)
option (OGS_OUTPUT_TIME_ALLPCS "Output clock time for all processes" OFF)

include(scripts/cmake/CMakeSetup.cmake)
include(scripts/cmake/CompilerSetup.cmake)
include (scripts/cmake/Find.cmake)
include(scripts/cmake/SubmoduleSetup.cmake)
include(scripts/cmake/GeneralProjectSetup.cmake)

if(GCC)
	if(GPROF_PATH)
		option(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
		if(OGS_PROFILE AND DOT_TOOL_PATH)
			option(OGS_OUTPUT_PROFILE "Enables pdf output of profiling results." OFF)
		endif()
	endif()

	OPTION_REQUIRES(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OGS_BUILD_TESTS)
endif()

mark_as_advanced(FORCE OGS_PYTHON OGS_CMAKE_DEBUG OGS_ONE_BIG_GTEST OGS_BUILD_INFO CMAKE_CMD_ARGS OGS_DELETE_EDGES_AFTER_INIT)

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
set (OGS_VERSION_MAJOR 5)
set (OGS_VERSION_MINOR 7)
set (OGS_VERSION_PATCH 0)
set (OGS_RELEASE_PERSONS "WH/WW/LB")
#set (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
set (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}(${OGS_RELEASE_PERSONS})")
# should we use ISO dates: 2010-12-31 ?
set (OGS_DATE "07.07.2015")


###################
### Fem configs ###
###################
message (STATUS	"OGS_CONFIG: ${OGS_CONFIG}" )
if(OGS_CONFIG STREQUAL FEM)
	set(SOLVER_PKG_NAME "Default")

elseif(OGS_CONFIG STREQUAL SP)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	set(SOLVER_PKG_NAME "Alternative linear solver package (NEW_EQS)")

elseif(OGS_CONFIG STREQUAL JFNK)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	add_definitions(-DJFNK_H2M)

elseif(OGS_CONFIG STREQUAL PETSC)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/scripts/cmake/cmake/findPETSC)

	find_package(PETSc REQUIRED QUIET)
	if(PETSC_VERSION VERSION_GREATER 3.3)
		add_definitions(-DUSEPETSC34 -DUSE_PETSC)
		message(STATUS "found version greater 3.3, version is ${PETSC_VERSION}")
	else()
		add_definitions(-DUSE_PETSC)
	endif()
	include_directories( ${PETSC_INCLUDES} )
	set(SOLVER_PKG_NAME "PETSc linear solver package")

elseif(OGS_CONFIG STREQUAL PETSC_GEMS)
	include_directories( ${Boost_INCLUDE_DIR} )
	link_directories( ${Boost_LIBRARY_DIR} )
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/scripts/cmake/cmake/findPETSC)

	find_package(PETSc REQUIRED QUIET)
	if(PETSC_VERSION VERSION_GREATER 3.3)
		add_definitions(-DUSEPETSC34 -DUSE_PETSC -DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN)
		message(STATUS "found version greater 3.3, version is ${PETSC_VERSION}")
	else()
		add_definitions(-DUSE_PETSC -DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN)
	endif()
	include_directories( ${PETSC_INCLUDES} )

	if (MSVC)
		set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	endif ()

elseif(OGS_CONFIG STREQUAL MPI)
	set(PARALLEL_USE_MPI ON CACHE INTERNAL "Use Message Passing Interface")
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	set(SOLVER_PKG_NAME "Parallel linear solver package (MPI)")

elseif(OGS_CONFIG STREQUAL MPI_KRC)
	add_definitions(-DUSE_MPI_KRC)
	set(PARALLEL_USE_MPI ON CACHE INTERNAL "Use Message Passing Interface")
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")

elseif(OGS_CONFIG STREQUAL GEMS)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	add_definitions(-DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN )
	if (MSVC)
		set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	endif ()

	include_directories( ${Boost_INCLUDE_DIR} )
	link_directories( ${Boost_LIBRARY_DIR} )

	if(PARALLEL_USE_MPI)
		add_definitions(-DUSE_MPI_GEMS)
	endif()

elseif(OGS_CONFIG STREQUAL IPQC)
	add_definitions( -DOGS_FEM_IPQC )
	if(PARALLEL_USE_MPI)
		set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	endif()

elseif(OGS_CONFIG STREQUAL PQC)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	add_definitions(-DLIBPHREEQC -D_POSIX_SOURCE -DNO_ERROR_CONTROL)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_definitions ( "-DMDL_DEBUG" )
	endif() # Debug

elseif(OGS_CONFIG STREQUAL BRNS)
	if (UNIX)
		find_package(LAPACK)
	endif()
	set(OGS_USE_NEW_EQS ON BOOL)
	add_definitions(-DBRNS -DNO_ERROR_CONTROL -D_POSIX_SOURCE)
	if(PARALLEL_USE_MPI)
		add_definitions(-DUSE_MPI_BRNS)
	endif()

elseif(OGS_CONFIG STREQUAL LIS)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	add_definitions(-DIPMGEMPLUGIN -DLIS -o3)
	set(PARALLEL_USE_OPENMP ON CACHE INTERNAL "Use Open Multi-Processing")
	set(LIS ON CACHE INTERNAL "Use Library of Iterative Solvers")
	set(SOLVER_PKG_NAME "LIS linear solver package")

elseif(OGS_CONFIG STREQUAL MKL)
	set(OGS_USE_NEW_EQS ON CACHE INTERNAL "")
	add_definitions(-DNO_ERROR_CONTROL -DMKL -DLIS -o3 -DIPMGEMPLUGIN)
	set(PARALLEL_USE_OPENMP ON CACHE INTERNAL "Use Open Multi-Processing")
	set(LIS ON INTERNAL "Use Library of Iterative Solvers")
	set(MKL ON INTERNAL "Use Math Kernel Library")
	mark_as_advanced(LIS MKL)
	set(SOLVER_PKG_NAME "MKL direct linear solver package")

elseif(OGS_CONFIG STREQUAL CAP)
	set( PRJ_EXT "-CAP" )
	add_definitions(-DOGS_FEM_CAP)
	if (UNIX) # Only supported on Linux
		include_directories(/usr/lib64)
		link_directories(/usr/lib64)
	endif ()

endif()

### End fem configs ###

if(OGS_USE_NEW_EQS)
	message(STATUS "Use New EQS class")
	add_definitions(-DNEW_EQS)
endif()

if(PARALLEL_USE_MPI)
	message(STATUS "Enabling MPI")
	find_package(MPI REQUIRED)
	add_definitions(-DUSE_MPI)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

if(PARALLEL_USE_OPENMP)
	find_package( OpenMP REQUIRED )
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
endif()

if(OGS_DELETE_EDGES_AFTER_INIT)
	add_definitions(-DOGS_DELETE_EDGES_AFTER_INIT)
endif()

if(OGS_OUTPUT_TIME_ALLPCS)
	add_definitions(-DOUTPUT_TIME_FOR_ALLPCS)
endif()

# Add subdirectories with the projects
add_subdirectory( ThirdParty )
include_directories (SYSTEM ${PROJECT_SOURCE_DIR}/ThirdParty)
if(OGS_USE_CVODE)
	add_definitions(-DOGS_USE_CVODE)
	include_directories (SYSTEM
		${PROJECT_SOURCE_DIR}/ThirdParty/cvode/include
		${PROJECT_BINARY_DIR}/ThirdParty/cvode/include
	)
endif()

add_subdirectory( Base )
add_subdirectory( MathLib )
add_subdirectory( GEO )
add_subdirectory( MSH )
add_subdirectory( MSHGEOTOOLS )
add_subdirectory( FEM )
add_subdirectory( GCC )
if (OGS_CONFIG STREQUAL CAP)
	add_subdirectory( CMP )
endif ()
add_subdirectory( FileIO )
add_subdirectory( OGS )


# Create Configure.h
configure_file (Base/BuildInfo.cpp.in ${PROJECT_BINARY_DIR}/Base/BuildInfo.cpp)

if (EXAMPLEDATA_DIR_FOUND)
	install (DIRECTORY ${EXAMPLEDATA_DIR_FOUND} DESTINATION . PATTERN .svn EXCLUDE)
endif ()

include(scripts/cmake/test/Test.cmake)
if (OGS_BUILD_TESTS)
	add_subdirectory(tests)
endif ()

include( scripts/cmake/packaging/Pack.cmake)

if (OGS_BUILD_UTILITIES)
	add_subdirectory (UTL/MSHGEOTOOLS/)
	add_subdirectory (UTL/GIS2FEM/)
endif ()

## Documentation ##
if(DOXYGEN_FOUND)
	option(DOCS_GENERATE_DIAGRAMS "Use the DOT tool to generate class diagrams." OFF)
	option(DOCS_GENERATE_CALL_GRAPHS "Generate call dependency graphs." OFF)
	option(DOCS_GENERATE_COLLABORATION_GRAPHS "Generate collaboration graphs." OFF)
	if(DOT_TOOL_PATH AND DOCS_GENERATE_DIAGRAMS)
		set(DOT_FOUND "YES" CACHE INTERNAL "")
	endif() # DOT_TOOL_PATH AND DOCS_GENERATE_DIAGRAMS
	if(DOCS_GENERATE_CALL_GRAPHS)
		set(DOCS_GENERATE_CALL_GRAPHS_STRING "YES" CACHE INTERNAL "")
	endif() # DOCS_GENERATE_CALL_GRAPHS
	if(DOCS_GENERATE_COLLABORATION_GRAPHS)
		set(DOCS_GENERATE_COLLABORATION_GRAPHS_STRING "YES" CACHE INTERNAL "")
	endif() # DOCS_GENERATE_COLLABORATION_GRAPHS
	get_filename_component(DOT_TOOL_PATH_ONLY ${DOT_TOOL_PATH} PATH)
	configure_file(scripts/doxygen/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)
	add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating source code documentation with Doxygen." VERBATIM)
endif() # DOXYGEN_FOUND
