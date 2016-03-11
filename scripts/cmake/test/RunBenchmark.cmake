set(BENCHMARK_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks/${benchmarkDir})
foreach(OUTPUT_FILE ${OUTPUT_FILES})
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${BENCHMARK_OUTPUT_DIRECTORY}/${OUTPUT_FILE})
endforeach()
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${BENCHMARK_OUTPUT_DIRECTORY})

if (WIN32)

	execute_process (
		COMMAND ${EXECUTABLE_OUTPUT_PATH}/Release/ogs --output-directory ${BENCHMARK_OUTPUT_DIRECTORY} ${benchmarkStrippedName}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
		TIMEOUT ${BENCHMARK_TIMEOUT}
		RESULT_VARIABLE EXIT_CODE)

else ()

	if(PARALLEL_USE_MPI)
		set(MPI_RUN_COMMAND "mpirun" "-np" "${NUM_PROCESSES}")
	else()
		set(MPI_RUN_COMMAND "")
	endif()

	if(OGS_PROFILE AND NOT PARALLEL_USE_MPI)
		message(STATUS "Profiling benchmark")
		if(OGS_OUTPUT_PROFILE)
			message(STATUS "Executing gprof2dot.py")
			execute_process (
				COMMAND ${EXECUTABLE_OUTPUT_PATH}/ogs --output-directory ${BENCHMARK_OUTPUT_DIRECTORY} ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				RESULT_VARIABLE EXIT_CODE)

			# Run gprof2dot.py
			execute_process (
				COMMAND ${GPROF_PATH} ${EXECUTABLE_OUTPUT_PATH}/ogs --output-directory ${BENCHMARK_OUTPUT_DIRECTORY}
				COMMAND ${PROJECT_SOURCE_DIR}/scripts/gprof2dot.py -s -n 5.0 -e 1.0
				COMMAND ${DOT_TOOL_PATH} -Tpng -o ${PROJECT_SOURCE_DIR}/../benchmarks/results/${benchmarkStrippedName}.png
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir})

		else()
			execute_process (
				COMMAND ${GPROF_PATH} ${EXECUTABLE_OUTPUT_PATH}/ogs --output-directory ${BENCHMARK_OUTPUT_DIRECTORY} ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				RESULT_VARIABLE EXIT_CODE)
		endif()
	else()
		execute_process (
			COMMAND ${MPI_RUN_COMMAND} ${EXECUTABLE_OUTPUT_PATH}/ogs --output-directory ${BENCHMARK_OUTPUT_DIRECTORY} ${benchmarkStrippedName}
			OUTPUT_FILE ${BENCHMARK_OUTPUT_DIRECTORY}/${benchmarkStrippedName}.log
			ERROR_FILE ${BENCHMARK_OUTPUT_DIRECTORY}/${benchmarkStrippedName}.error
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
			RESULT_VARIABLE EXIT_CODE)
	endif()

endif ()

if(EXIT_CODE GREATER 0)
	message(FATAL_ERROR "Benchmark exited with code: ${EXIT_CODE}")
endif()

# file comparison with numdiff or cmake
if(NUMDIFF_TOOL_PATH)
	set(TESTER_COMMAND ${NUMDIFF_TOOL_PATH})
	set(TESTER_ARGS --absolute-tolerance=1e-5 --relative-tolerance=1e-4)
	message("Using numdiff")
else()
	set(TESTER_COMMAND ${CMAKE_COMMAND} -E compare_files)
	set(TESTER_ARGS "")
	message("Using cmake diff")
endif()
set(TESTER_COMMAND ${TESTER_COMMAND} )

set(SCRIPT_EXIT_CODE 0)
set(NUMDIFF_OUTPUT_FILE ${BENCHMARK_OUTPUT_DIRECTORY}/${benchmarkStrippedName}.numdiff)
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NUMDIFF_OUTPUT_FILE})
separate_arguments(OUTPUT_FILES) # reconstruct list
foreach(OUTPUT_FILE ${OUTPUT_FILES})
	execute_process(COMMAND ${TESTER_COMMAND} ${TESTER_ARGS}
		${BENCHMARK_REF_DIR}/${benchmarkDir}/${OUTPUT_FILE}
		${BENCHMARK_OUTPUT_DIRECTORY}/${OUTPUT_FILE}
		RESULT_VARIABLE EXIT_CODE
		OUTPUT_QUIET
	)

	if(EXIT_CODE GREATER 0)
		if(NUMDIFF_TOOL_PATH)
			execute_process(COMMAND ${TESTER_COMMAND} ${TESTER_ARGS} -E -S
				${BENCHMARK_REF_DIR}/${benchmarkDir}/${OUTPUT_FILE}
				${BENCHMARK_OUTPUT_DIRECTORY}/${OUTPUT_FILE}
				ERROR_VARIABLE NUMDIFF_ERROR
				OUTPUT_VARIABLE NUMDIFF_OUT
			)
			file(APPEND ${NUMDIFF_OUTPUT_FILE} "### Numdiff error ###\n${NUMDIFF_ERROR}\n\n")
			file(APPEND ${NUMDIFF_OUTPUT_FILE} "### Numdiff output ###\n${NUMDIFF_OUT}\n\n")
		endif()
		set(SCRIPT_EXIT_CODE 1)
		message(STATUS "Benchmark file compare of ${OUTPUT_FILE} failed.")
	endif()
endforeach()

if(SCRIPT_EXIT_CODE GREATER 0)
	if(NUMDIFF_TOOL_PATH)
		message(FATAL_ERROR "Benchmark ${benchmarkDir}/${benchmarkStrippedName} failed.\n See ${NUMDIFF_OUTPUT_FILE} for details.")
	else()
		message(FATAL_ERROR "Benchmark ${benchmarkDir}/${benchmarkStrippedName} failed.")
	endif()
endif()
