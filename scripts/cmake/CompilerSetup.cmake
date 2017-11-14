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
	if (MSVC)
		# enable parallel compilation
		# specify Exception Handling Model in msvc
		# disable unknown pragma warning (4068)
		# disable unsafe function warning (4996)
		# disable decorated name length exceeded, name was truncated (4503)
		# disable conversion from 'size_t' to 'type', possible loss of data (4267)
		# disable qualifier applied to function type has no meaning; ignored (4180)
		# disable C++ exception specification ignored except to indicate a function is not __declspec(nothrow) (4290)
		# disable conversion from 'type1' to 'type2', possible loss of data (4244)
		# disable forcing value to bool 'true' or 'false' (performance warning) (4800)
		# define miniupnp static library
		add_compile_options(/MP /EHsc /wd4068 /wd4996 /wd4503 /wd4267 /wd4180 /wd4290 /wd4244 /wd4800)

		add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS)
		# Sets warning level 3 and ignores some warnings
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
		set(GCC OFF)

		DisableCompilerFlag(DEBUG /RTC1)
	endif ()
endif ()

if(COMPILER_IS_GCC)
	set(GCC ON)
	if( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		if (OGS_CHEMSOLVER STREQUAL GEMS)
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
