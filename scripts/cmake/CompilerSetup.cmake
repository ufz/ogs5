include(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
include(SetDefaultBuildType)
include(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Release)
include(MSVCMultipleProcessCompile) # /MP Switch for VS
set(CMAKE_OSX_ARCHITECTURES "x86_64")

if(OGS_CPU_ARCHITECTURE STREQUAL "generic")
	set(CPU_FLAGS "-mtune=generic")
else()
	set(CPU_FLAGS "-march=${OGS_CPU_ARCHITECTURE}")
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
	set(COMPILER_IS_CLANG TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(COMPILER_IS_GCC TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	set(COMPILER_IS_INTEL TRUE)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	set(COMPILER_IS_MSVC TRUE)
endif() # CMAKE_CXX_COMPILER_ID

if (WIN32)
	if (COMPILER_IS_MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS)
		# Sets warning level 3 and ignores some warnings
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")
		set(GCC OFF)

		DisableCompilerFlag(DEBUG /RTC1)
	endif ()
endif ()

if(COMPILER_IS_GCC)
	set(GCC ON)
	if( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		if (OGS_CONFIG STREQUAL PETSC_GEMS)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
		endif()
	endif()

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPU_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")

	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	string(REPLACE "\n" "" GCC_VERSION ${GCC_VERSION})
	message(STATUS "GCC_VERSION: ${GCC_VERSION}")
	if (NOT (GCC_VERSION VERSION_LESS 4.8) )
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") # suppress warnings in Eigen
	endif()

	add_definitions(-DGCC)

	if (OGS_PROFILE)
		if( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
			message(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
		endif()
		set(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
	endif ()
endif()

if(COMPILER_IS_CLANG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPU_FLAGS}")
endif()

if(BLUE_G)
	add_definitions(-O3 -qstrict -qarch=qp -qtune=qp)
endif()

if(OGS_COVERAGE)
  include(CodeCoverage)
endif()
