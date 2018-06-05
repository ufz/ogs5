#### OGS-5 Project ####

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.8.3)

if(CMAKE_VERSION VERSION_GREATER 3.1)
	cmake_policy(SET CMP0054 NEW)
endif()
cmake_policy(SET CMP0011 OLD)

# Project name
project(OGS-5)

###############
### Options ###
###############

# Chemical solver
set(OGS_CHEMSOLVER "DEFAULT" CACHE STRING "OpenGeoSys geochemical solver")
set_property(CACHE OGS_CHEMSOLVER PROPERTY STRINGS
	DEFAULT # Default
	NONE
	GEMS
	IPQC
	PQC
	BRNS
	CAP
)

# Linear solver
set(OGS_LSOLVER "DEFAULT" CACHE STRING "OpenGeoSys linear solver type")
set_property(CACHE OGS_LSOLVER PROPERTY STRINGS
	DEFAULT # Default
	RF      # RF
	SP      # New sparse solver
	PETSC   # PETSc
)

# Other options
option(BLUE_G "Blue/G optimization" OFF)
option(PARALLEL_USE_MPI "Use MPI parallization" OFF)
option(OGS_USE_JFNK "Use Jacobain free method for solving H2M" OFF)
option(OGS_USE_LIS "Use Lis solver" OFF)
option(OGS_USE_MKL "Use PARDISO in MKL" OFF)
option(OGS_USE_CVODE "Use sundials cvode for TNEQ/TES" OFF)
set(OGS_CPU_ARCHITECTURE "native" CACHE STRING "Processor architecture, defaults to native.")
option(OGS_PYTHON "Enable python scripting interface" OFF)
option(OGS_BUILD_UTILITIES "Build additional utility programs?" OFF)
option(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
option(OGS_ONE_BIG_GTEST "Should all gtests be collected to one ctest" ON)
option(OGS_BUILD_INFO "Should build information be generated" ON)
option(OGS_BUILD_TESTS "Enables building of tests." OFF)
option(OGS_NO_EXTERNAL_LIBS "Builds OGS without any external dependencies." OFF)
option(OGS_DELETE_EDGES_AFTER_INIT "Delete mesh edges after initialization if possible" OFF)
option(OGS_OUTPUT_TIME_ALLPCS "Output clock time for all processes" OFF)

##########################
### Configuration type ###
##########################
set(OGS_CONFIG "FEM" CACHE STRING "OpenGeoSys configuration")
set_property(CACHE OGS_CONFIG PROPERTY STRINGS
	FEM  # RF
	SP   # SP
	JFNK # SP + JFNK
	GEMS # SP + GEMS
	IPQC # RF + IPQC
	PQC  # SP + PQC
	BRNS # SP + BRNS
	CAP  # RF + CAP
	LIS  # SP + LIS
	MKL  # SP + LIS + MKL
	MPI  # SP(MPI)
	PETSC # PETSC(MPI)
	PETSC_GEMS # PETSC(MPI) + GEMS
)

# default setting for each configuration
set(OGS_CONFIG_DEFAULT_LSOLVER "RF" CACHE INTERNAL "OGS default linear solver type for the given configuration")
set(OGS_CONFIG_DEFAULT_CHEMSOLVER "NONE" CACHE INTERNAL "OGS default geochemical solver library")

if (OGS_CONFIG STREQUAL SP)
	set(OGS_CONFIG_DEFAULT_LSOLVER　"SP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL JFNK)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL GEMS)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "GEMS" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL IPQC)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "IPQC" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL PQC)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "PQC" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL BRNS)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "BRNS" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL CAP)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "CAP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL LIS)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL MKL)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL MPI)
	set(OGS_CONFIG_DEFAULT_LSOLVER "SP" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL PETSC)
	set(OGS_CONFIG_DEFAULT_LSOLVER "PETSC" CACHE INTERNAL "" FORCE)
elseif (OGS_CONFIG STREQUAL PETSC_GEMS)
	set(OGS_CONFIG_DEFAULT_LSOLVER "PETSC" CACHE INTERNAL "" FORCE)
	set(OGS_CONFIG_DEFAULT_CHEMSOLVER "GEMS" CACHE INTERNAL "" FORCE)
endif()

# overwrite options if necessary
if (OGS_LSOLVER STREQUAL DEFAULT OR NOT OGS_CONFIG EQUAL OGS_CONFIG_old)
	set(OGS_LSOLVER ${OGS_CONFIG_DEFAULT_LSOLVER} CACHE STRING "" FORCE)
endif()
if (OGS_CHEMSOLVER STREQUAL DEFAULT OR NOT OGS_CONFIG EQUAL OGS_CONFIG_old)
	set(OGS_CHEMSOLVER ${OGS_CONFIG_DEFAULT_CHEMSOLVER} CACHE STRING "" FORCE)
endif()
set(OGS_CONFIG_old ${OGS_CONFIG} CACHE INTERNAL "Copy of OGS_CONFIG")

if (OGS_CONFIG STREQUAL JFNK)
	set(OGS_USE_JFNK ON)
elseif (OGS_CONFIG STREQUAL LIS)
	set(OGS_USE_LIS ON)
elseif (OGS_CONFIG STREQUAL MKL)
	set(OGS_USE_LIS ON)
	set(OGS_USE_MKL ON)
elseif (OGS_CONFIG STREQUAL MPI)
	set(PARALLEL_USE_MPI ON)
elseif (OGS_CONFIG STREQUAL PETSC)
	set(PARALLEL_USE_MPI ON)
elseif (OGS_CONFIG STREQUAL PETSC_GEMS)
	set(PARALLEL_USE_MPI ON)
endif()

#########################
### Finding libraries ###
#########################

if (OGS_CHEMSOLVER STREQUAL GEMS)
	# need to set this before Find.cmake
	set(OGS_USE_BOOST ON CACHE INTERNAL "Use Boost library")
endif()

include(scripts/cmake/CMakeSetup.cmake)
include(scripts/cmake/CompilerSetup.cmake)
include(scripts/cmake/Find.cmake)
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
### OGS version information is now obtained from Git tags. Adjust authors & date.
###########################################################################
set (OGS_RELEASE_PERSONS "LB/NW/TN/WW")
set (OGS_DATE "09.08.2016")

###################
### Fem configs ###
###################
message (STATUS	"OGS_CONFIG: ${OGS_CONFIG}" )
message (STATUS	"OGS_CHEMSOLVER: ${OGS_CHEMSOLVER}" )
message (STATUS	"OGS_LSOLVER: ${OGS_LSOLVER}" )

# chemical solvers
if(OGS_CHEMSOLVER STREQUAL GEMS)
	add_definitions(-DNO_ERROR_CONTROL -DGEM_REACT -DGEMlib -DIPMGEMPLUGIN )
	if (MSVC)
		set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	endif ()

	include_directories( ${Boost_INCLUDE_DIR} )
	link_directories( ${Boost_LIBRARY_DIR} )

elseif(OGS_CHEMSOLVER STREQUAL IPQC)
	add_definitions( -DOGS_FEM_IPQC )

elseif(OGS_CHEMSOLVER STREQUAL PQC)
	add_definitions(-DLIBPHREEQC -D_POSIX_SOURCE -DNO_ERROR_CONTROL)
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_definitions ( "-DMDL_DEBUG" )
	endif() # Debug

elseif(OGS_CHEMSOLVER STREQUAL BRNS)
	if (UNIX)
		find_package(LAPACK)
	endif()
	add_definitions(-DBRNS -DNO_ERROR_CONTROL -D_POSIX_SOURCE)

elseif(OGS_CHEMSOLVER STREQUAL CAP)
	set( PRJ_EXT "-CAP" )
	add_definitions(-DOGS_FEM_CAP)
	if (UNIX) # Only supported on Linux
		include_directories(/usr/lib64)
		link_directories(/usr/lib64)
	endif ()

endif()

# linear solvers
if(OGS_USE_JFNK)
	set(OGS_LSOLVER "SP")
	add_definitions(-DJFNK_H2M)
endif()

if(OGS_USE_LIS)
	set(OGS_LSOLVER "SP")
	add_definitions(-DIPMGEMPLUGIN -DLIS -o3)
	set(PARALLEL_USE_OPENMP ON CACHE INTERNAL "Use Open Multi-Processing")
	set(LIS ON CACHE INTERNAL "Use Library of Iterative Solvers")
endif()

if(OGS_USE_MKL)
	set(OGS_LSOLVER "SP")
	add_definitions(-DNO_ERROR_CONTROL -DMKL -o3 -DIPMGEMPLUGIN)
	set(PARALLEL_USE_OPENMP ON CACHE INTERNAL "Use Open Multi-Processing")
	set(MKL ON INTERNAL "Use Math Kernel Library")
endif()

if(OGS_LSOLVER STREQUAL RF)
	set(SOLVER_PKG_NAME "Default")

elseif(OGS_LSOLVER STREQUAL SP)
	add_definitions(-DNEW_EQS)
	if(PARALLEL_USE_MPI)
		add_definitions(-DUSE_MPI)
	endif()
	set(SOLVER_PKG_NAME "Alternative linear solver package (NEW_EQS)")

elseif(OGS_LSOLVER STREQUAL PETSC)
	message(STATUS "Configuring for PETSc")

	option(FORCE_PETSC_EXECUTABLE_RUNS "Force CMake to accept a given PETSc configuration" ON)

	##Force CMake to accept a given PETSc configuration in case the failure of MPI tests
	##This may cause the compilation broken.
	if(FORCE_PETSC_EXECUTABLE_RUNS)
		set(PETSC_EXECUTABLE_RUNS YES)
	endif()

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/scripts/cmake/cmake/findPETSC)

	find_package(PETSc REQUIRED QUIET)
	message(STATUS "PETSc version is ${PETSC_VERSION}")
	include_directories( ${PETSC_INCLUDES} )
	add_definitions(-DUSE_PETSC)
	add_definitions(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)
	set(SOLVER_PKG_NAME "PETSc linear solver package")
endif()

### End fem configs ###

if(PARALLEL_USE_MPI)
	message(STATUS "Enabling MPI")
	find_package(MPI REQUIRED)
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

include(scripts/cmake/packaging/Pack.cmake)

if (OGS_CONFIG STREQUAL FEM OR OGS_CONFIG STREQUAL SP)
    set(OGS_BUILD_UTILITIES ON)
endif ()

if (OGS_BUILD_UTILITIES)
	add_subdirectory (UTL/MSHGEOTOOLS/)
	add_subdirectory (UTL/mHM2OGS/)
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
