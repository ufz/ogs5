#
# Benchmark
# -------
#
# Creates a benchmark run. Order of arguments can be arbitrary.
#
# Benchmark(
#   AUTHOR <Author short code>
#   PATH <benchmark path incl. name>  # relative to the Benchmarks folder
#   CONFIG <ogs configuration> # defaults to FEM
#   REQUIRED_CMAKE_OPTIONS <0..N CMake options> # optional
#   OUTPUT_FILES <0..N names of result files to compare>
#   NUM_PROCESSORS <number of processors to run this test on> # Is passed to mpirun
#   TIMEOUT <Time in minutes> # Build is aborted if timeout is exceeded
#   RUNTIME <Expected runtime in seconds>
# )
#
# Example:
#
# Benchmark(AUTHOR YS
#   PATH RWPT/Forchheimer/forchheimer_rwpt
#   CONFIG FEM
#   OUTPUT_FILES
#     forchheimer_rwpt_GROUNDWATER_FLOW60.vtu
#     forchheimer_rwpt_domain_ele_GROUNDWATER_FLOW.tec
#   RUNTIME 5
# )
#

function(Benchmark)

	# parse arguments
	set(options NONE)
	set(oneValueArgs AUTHOR PATH NUM_PROCESSORS TIMEOUT RUNTIME CONFIG)
	set(multiValueArgs REQUIRED_CMAKE_OPTIONS OUTPUT_FILES)
	cmake_parse_arguments(Benchmark "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

	# set defaults
	if(NOT Benchmark_CONFIG)
		set(Benchmark_CONFIG FEM)
	endif()
	if(NOT Benchmark_NUM_PROCESSORS)
		set(Benchmark_NUM_PROCESSORS 1)
	endif()
	if(NOT Benchmark_TIMEOUT)
		set(Benchmark_TIMEOUT ${BENCHMARK_TIMEOUT})
	endif()
	if(NOT Benchmark_RUNTIME)
		set(Benchmark_RUNTIME 1)
	endif()

	# Check required CMake configuration
	if(NOT Benchmark_CONFIG STREQUAL OGS_CONFIG)
		return()
	endif()
	foreach(REQUIRED_CMAKE_OPTION ${Benchmark_REQUIRED_CMAKE_OPTIONS})
		if(NOT ${REQUIRED_CMAKE_OPTION})
			return()
			# message("Disabling benchmark ${Benchmark_NAME} because ${REQUIRED_CMAKE_OPTION} = ${${REQUIRED_CMAKE_OPTION}}")
		endif()
	endforeach()

	# set a label for benchmark categorization
	if(Benchmark_RUNTIME GREATER 600)    # > 10 minutes
		set(TEST_LABEL exceeding)
	elseif(Benchmark_RUNTIME GREATER 60) # 1 - 10 minutes
		set(TEST_LABEL long)
	elseif(Benchmark_RUNTIME GREATER 10) # 10 - 60 seconds
		set(TEST_LABEL normal)
	else()
		set(TEST_LABEL short)            # < 10 seconds
	endif()

	string(REGEX MATCH "[^/]+$" benchmarkStrippedName ${Benchmark_PATH})
	string(LENGTH ${Benchmark_PATH} benchmarkNameLength)
	string(LENGTH ${benchmarkStrippedName} benchmarkStrippedNameLength)
	math(EXPR substringLength ${benchmarkNameLength}-${benchmarkStrippedNameLength})
	string(SUBSTRING ${Benchmark_PATH} 0 ${substringLength} benchmarkDir)
	string(REGEX REPLACE "/$" "" benchmarkDir ${benchmarkDir}) # Remove trailing slash

	if(BENCHMARK_REF_DIR_FOUND)
		# Lists must be passed as a space-separated string, separate_arguments() restores the list.
		string (REPLACE ";" " " Benchmark_OUTPUT_FILES "${Benchmark_OUTPUT_FILES}")
	else()
		set(Benchmark_OUTPUT_FILES "")
	endif()

	# Adds a benchmark run. This calls AddTest.cmake to execute several steps.
	add_test(
		${Benchmark_AUTHOR}-${Benchmark_PATH}
		${CMAKE_COMMAND}
		-DBENCHMARK_DIR=${BENCHMARK_DIR_FOUND}
		-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
		-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
		-DEXECUTABLE_OUTPUT_PATH=${EXECUTABLE_OUTPUT_PATH}
		-DbenchmarkStrippedName=${benchmarkStrippedName}
		-DbenchmarkDir=${benchmarkDir}
		-DOGS_PROFILE=${OGS_PROFILE}
		-DOGS_OUTPUT_PROFILE=${OGS_OUTPUT_PROFILE}
		-DGPROF_PATH=${GPROF_PATH}
		-DDOT_TOOL_PATH=${DOT_TOOL_PATH}
		-DNUMDIFF_TOOL_PATH=${NUMDIFF_TOOL_PATH}
		-DNUM_PROCESSES=${Benchmark_NUM_PROCESSORS}
		-DOUTPUT_FILES=${Benchmark_OUTPUT_FILES}
		-DBENCHMARK_REF_DIR=${BENCHMARK_REF_DIR_FOUND}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/RunBenchmark.cmake
	)
	set_tests_properties(${Benchmark_AUTHOR}-${Benchmark_PATH} PROPERTIES
		PROCESSORS ${Benchmark_NUM_PROCESSORS}
		TIMEOUT ${Benchmark_TIMEOUT}
		COST ${Benchmark_RUNTIME}
		LABELS "${TEST_LABEL};all;${Benchmark_AUTHOR}"
	)

endfunction()
