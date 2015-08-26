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
		# Get git commit
		execute_process(
			COMMAND ${GIT_EXECUTABLE} "log" "--name-status" "HEAD^..HEAD"
			COMMAND ${GREP_TOOL_PATH} "-m" "1" "commit"
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_COMMIT_INFO
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
	endif() # GIT_FOUND

	find_path(HIDDEN_SVN_DIR entries ${CMAKE_SOURCE_DIR}/.svn)
	if(Subversion_FOUND AND HIDDEN_SVN_DIR)
		Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
		set(SVN_REVISION ${Project_WC_REVISION})
	endif() # Subversion_FOUND AND HIDDEN_SVN_DIR
	unset(HIDDEN_SVN_DIR)

	execute_process(
		COMMAND ${DATE_TOOL_PATH} "+%Y-%m-%d" # %H:%M:%S"
		OUTPUT_VARIABLE BUILD_TIMESTAMP
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

endif () # OGS_BUILD_INFO

# This is for Configure.h which is generated later
include_directories( ${PROJECT_BINARY_DIR}/Base )

# Check for number of processors
include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(PROCESSOR_COUNT EQUAL 0)
	message(WARNING "Processor count could not be detected. Setting to one processor.")
	set(PROCESSOR_COUNT 1)
else()
	message(STATUS "Number of processors: ${PROCESSOR_COUNT}")
endif()
