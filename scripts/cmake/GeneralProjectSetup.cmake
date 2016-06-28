# Set build directories
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
set( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
if (MSVC)
	set(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/release/ogs)
else ()
	set(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/ogs)
endif ()

# Collect build information such as revision/commit and timestamp
if (OGS_BUILD_INFO)
	if(GIT_FOUND)
		include(GetGitRevisionDescription)
		GET_GIT_HEAD_REVISION(GIT_REFSPEC GIT_SHA1)
		string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1_SHORT)
	endif()

	execute_process(
		COMMAND ${DATE_TOOL_PATH} "+%Y-%m-%d" # %H:%M:%S"
		OUTPUT_VARIABLE BUILD_TIMESTAMP
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

endif () # OGS_BUILD_INFO

# This is for Configure.h which is generated later
include_directories( ${PROJECT_BINARY_DIR}/Base )

# Check for number of processors
if (NOT DEFINED PROCESSOR_COUNT)
	include(ProcessorCount)
	ProcessorCount(PROCESSOR_COUNT)
endif()
if(PROCESSOR_COUNT EQUAL 0)
	message(WARNING "Processor count could not be detected. Setting to one processor.")
	set(PROCESSOR_COUNT 1)
else()
	message(STATUS "Number of processors: ${PROCESSOR_COUNT}")
endif()

# Release and git info
include(GetGitRevisionDescription)
GET_GIT_HEAD_REVISION(GIT_REFSPEC GIT_SHA1)
string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHA1_SHORT)

if($ENV{CI})
	set(OGS_VERSION 5.x.x) # Dummy version for CI-environment (Travis) or subproject
elseif(IS_SUBPROJECT)
	set(OGS_VERSION x.x.x)
else()
	GIT_GET_TAG(GIT_DESCRIBE)
	if(GIT_DESCRIBE)
		string(REGEX MATCH ^[0-9|\\.]+ GIT_TAG ${GIT_DESCRIBE})
		set(OGS_VERSION ${GIT_TAG})
        set(OGS_IS_RELEASE FALSE)
		if(GIT_DESCRIBE MATCHES ".*-.*-.*")
			# Commit is not a tag
            set(OGS_VERSION ${GIT_DESCRIBE})
		else()
			set(OGS_IS_RELEASE TRUE)
		endif()
        message(STATUS "OGS version: ${OGS_VERSION}")
	else()
		message(WARNING "Git repository contains no tags! Please run: git fetch --tags")
	endif()
endif()
