# This script is called from AddBenchmark in Macros.cmake
# It deletes the benchmark output files and then runs the benchmark.
if (WIN32)

	separate_arguments(FILES_TO_DELETE_VARS WINDOWS_COMMAND ${FILES_TO_DELETE})

	execute_process (
		COMMAND del /S /Q ${FILES_TO_DELETE_VARS}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks)

	execute_process (
		COMMAND ${EXECUTABLE_OUTPUT_PATH}/Release/ogs ${benchmarkStrippedName}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
		TIMEOUT ${BENCHMARK_TIMEOUT}
		RESULT_VARIABLE EXIT_CODE)

else ()

	separate_arguments(FILES_TO_DELETE_VARS UNIX_COMMAND ${FILES_TO_DELETE})

	if(OGS_FEM_CONFIG STREQUAL "OGS_FEM_PETSC" OR OGS_FEM_CONFIG STREQUAL "OGS_FEM_PETSC_GEMS" OR OGS_FEM_CONFIG STREQUAL "OGS_FEM_MPI")
		set(MPI_RUN_COMMAND "mpirun" "-np" "${NUM_PROCESSES}")
	else()
		set(MPI_RUN_COMMAND "")
	endif()

	execute_process (
		COMMAND rm -f ${FILES_TO_DELETE_VARS}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks)
	if(OGS_PROFILE AND NOT (OGS_FEM_CONFIG STREQUAL "OGS_FEM_PETSC" OR OGS_FEM_CONFIG STREQUAL "OGS_FEM_PETSC_GEMS" OR OGS_FEM_CONFIG STREQUAL "OGS_FEM_MPI"))
		message(STATUS "Profiling benchmark")
		if(OGS_OUTPUT_PROFILE)
			message(STATUS "Executing gprof2dot.py")
			execute_process (
				COMMAND ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				TIMEOUT ${BENCHMARK_TIMEOUT}
				RESULT_VARIABLE EXIT_CODE)

			# Run gprof2dot.py
			execute_process (
				COMMAND ${GPROF_PATH} ${EXECUTABLE_OUTPUT_PATH}/ogs
				COMMAND ${PROJECT_SOURCE_DIR}/scripts/gprof2dot.py -s -n 5.0 -e 1.0
				COMMAND ${DOT_TOOL_PATH} -Tpng -o ${PROJECT_SOURCE_DIR}/../benchmarks/results/${benchmarkStrippedName}.png
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir})

		else()
			execute_process (
				COMMAND ${GPROF_PATH} ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				TIMEOUT ${BENCHMARK_TIMEOUT}
				RESULT_VARIABLE EXIT_CODE)
		endif()
	else()
		execute_process (
			COMMAND ${MPI_RUN_COMMAND} ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
			TIMEOUT ${BENCHMARK_TIMEOUT}
			RESULT_VARIABLE EXIT_CODE)
	endif()

endif ()

if(EXIT_CODE GREATER 0)
	message(FATAL_ERROR "Benchmark exited with code: ${EXIT_CODE}")
endif()
