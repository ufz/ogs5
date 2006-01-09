include(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
include(SetDefaultBuildType)
include(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Release)
include(MSVCMultipleProcessCompile) # /MP Switch for VS

if (WIN32)
	## For Visual Studio compiler
	if (MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS)
		# Sets warning level 3 and ignores some warnings
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")
		set(GCC OFF)

		DisableCompilerFlag(DEBUG /RTC1)

    # Set $PATH to Visual Studio bin directory. Needed for finding dumpbin.exe
    if (MSVC80)
      set(ENV{PATH} "$ENV{PATH};$ENV{VS80COMNTOOLS}..\\..\\VC\\bin")
    endif ()
    if (MSVC90)
      set(ENV{PATH} "$ENV{PATH};$ENV{VS90COMNTOOLS}..\\..\\VC\\bin")
    endif ()
    if (MSVC10)
      set(ENV{PATH} "$ENV{PATH};$ENV{VS100COMNTOOLS}..\\..\\VC\\bin")
    endif ()

	else ()
#FOR CYGWIN.  25.02.2010. WW
		message (STATUS "Might be GCC under cygwin.")
		set(GCC ON)
#		message (FATAL_ERROR "Aborting: On Windows only the Visual Studio compiler is supported!")
	endif ()
endif ()

### For GNU C/CXX. WW
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	set(GCC ON)
	if( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		if (OGS_FEM_PETSC_GEMS)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
		endif()
	endif()
	# -g
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wall -Wextra -fno-nonansi-builtins")

	execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	string(REPLACE "\n" "" GCC_VERSION ${GCC_VERSION})
	message(STATUS "GCC_VERSION: ${GCC_VERSION}")
	if (NOT (GCC_VERSION VERSION_LESS 4.8) )
	  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") # suppress warnings in Eigen
	endif()

	# would be cool: -Woverloaded-virtual, would be overkill: -Weffc++
	add_definitions(-DGCC)

	if (OGS_PROFILE)
		if( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
			message(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
		endif()
		set(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
	endif ()
endif() # CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC

if(BLUE_G)
	add_definitions(-O3 -qstrict -qarch=qp -qtune=qp)
endif()

if(OGS_COVERAGE)
  include(CodeCoverage)
endif()
