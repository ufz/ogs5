if(DEFINED EXAMPLEDATA_DIR)
	find_path (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
else()
	find_path (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
endif()

######################
### Find tools     ###
######################

find_package (PythonInterp)
find_package(Subversion)
find_package(Git REQUIRED)

find_program(BASH_TOOL_PATH bash DOC "The bash executable")
find_program(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")
find_package(Doxygen)
find_program(GPROF_PATH gprof DOC "GNU profiler gprof")
find_package(cppcheck)

# Find Exuberant ctags or BBEdit for code completion
find_program(CTAGS_TOOL_PATH ctags DOC "Exuberant ctags")
if(CTAGS_TOOL_PATH)
	add_custom_target(ctags
		ctags -R --fields=+iamS -f ${CMAKE_SOURCES_DIR}/../tags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Creating tags..." VERBATIM
	)
endif()

find_program(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
find_program(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
find_program(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Hide these variables for the CMake user
mark_as_advanced(DOT_TOOL_PATH GPROF_PATH CTAGS_TOOL_PATH BBEDIT_TOOL_PATH
	UNZIP_TOOL_PATH
)
########################
### Find other stuff ###
########################

# Check if on Jenkins
if(NOT $ENV{JENKINS_URL} STREQUAL "")
	set(JENKINS_URL $ENV{JENKINS_URL})
	set(JENKINS_JOB_NAME $ENV{JOB_NAME})
endif()


######################
### Find libraries ###
######################
if(OGS_NO_EXTERNAL_LIBS)
	return()
endif()

find_path(OGS_LIBS_DIR_FOUND .dummy PATHS
    ${CMAKE_SOURCE_DIR}/../Libs
    /opt/ogs/ogs5-libs
    NO_DEFAULT_PATH)
if(OGS_LIBS_DIR_FOUND)
	# Find precompiled libraries (for BRNS GEMS LIS)
	find_path (OGS_PRECOMPILED_LIBS_DIR_FOUND BrnsDll.lib ${OGS_LIBS_DIR_FOUND}/precompiled)
	if (OGS_PRECOMPILED_LIBS_DIR_FOUND)
		include_directories (${OGS_LIBS_DIR_FOUND}/precompiled)
		link_directories (${OGS_LIBS_DIR_FOUND}/precompiled)
	else ()
		if (WIN32)
			if (OGS_CHEMSOLVER STREQUAL BRNS OR OGS_CHEMSOLVER STREQUAL GEMS)
				message (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
			endif ()
		else ()
			if (LIS OR MKL)
				message (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
			endif ()
		endif ()
	endif ()
endif()


## pthread ##
set ( CMAKE_THREAD_PREFER_PTHREAD ON CACHE BOOL "" )
find_package( Threads )
if ( CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREADS)
	set (HAVE_PTHREADS TRUE CACHE BOOL "Is PThreads found.")
	message (STATUS "pthread library found." )
endif ()
if(HAVE_PTHREADS)
	add_definitions(-DHAVE_PTHREADS)
endif()
mark_as_advanced(CMAKE_THREAD_PREFER_PTHREAD)

## boost (see FindBoost.cmake for more options) ##
if(OGS_USE_BOOST)
	##kg44 this configuration works for boost and petsc on a cray
	option(Boost_USE_STATIC_LIBS "" ON)
	option(Boost_USE_MULTITHREADED "" ON)
	option(Boost_USE_STATIC_RUNTIME "" ON)
	mark_as_advanced(Boost_USE_STATIC_LIBS Boost_USE_MULTITHREADED Boost_USE_STATIC_RUNTIME)

	find_package(Boost 1.50.0 COMPONENTS filesystem system thread regex REQUIRED)
endif()

# Find Math Kernel Library (MKL)
if(OGS_USE_MKL)
	if(OGS_LIBS_DIR_FOUND AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND NOT DEFINED MKL_DIR)
		message(STATUS "Using MKL from OGS Libs-folder. This implies PARALLEL_USE_OPENMP=ON.")
		if(NOT IS_DIRECTORY ${OGS_LIBS_DIR_FOUND}/MKL/include)
			message(STATUS "Extracting MKL include files...")
			execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${OGS_LIBS_DIR_FOUND}/MKL/mkl-include.tgz
				WORKING_DIRECTORY ${OGS_LIBS_DIR_FOUND}/MKL)
		endif()
		if(NOT IS_DIRECTORY ${OGS_LIBS_DIR_FOUND}/MKL/64)
			message(STATUS "Extracting MKL libraries...")
			execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf ${OGS_LIBS_DIR_FOUND}/MKL/mkl-64.tgz
				WORKING_DIRECTORY ${OGS_LIBS_DIR_FOUND}/MKL)
		endif()
		set(MKL_DIR ${OGS_LIBS_DIR_FOUND}/MKL)
		set(PARALLEL_USE_OPENMP ON CACHE BOOL "" FORCE)
	endif()
	find_package( MKL REQUIRED )
	include_directories (${MKL_INCLUDE_DIR})
endif()

if(OGS_USE_LIS)
	find_package( LIS REQUIRED )
endif()
