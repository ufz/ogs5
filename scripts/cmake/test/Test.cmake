if(DEFINED BENCHMARK_DIR)
	find_path (BENCHMARK_DIR_FOUND .benchmarks ${BENCHMARK_DIR})
else()
    find_path (BENCHMARK_DIR_FOUND .benchmarks
        ${PROJECT_SOURCE_DIR}/../benchmarks
        ${PROJECT_SOURCE_DIR}/benchmarks
    )
endif()

if(DEFINED BENCHMARK_REF_DIR)
	find_path (BENCHMARK_REF_DIR_FOUND .benchmarks_ref ${BENCHMARK_REF_DIR})
else()
    find_path (BENCHMARK_REF_DIR_FOUND .benchmarks_ref
        ${PROJECT_SOURCE_DIR}/../benchmarks_ref
        ${PROJECT_SOURCE_DIR}/benchmarks_ref
    )
endif()

if(DEFINED TESTDATA_DIR)
	find_path(TESTDATA_DIR_FOUND testdata.dummy ${TESTDATA_DIR})
else()
	find_path(TESTDATA_DIR_FOUND testdata.dummy ${PROJECT_SOURCE_DIR}/../testdata)
endif()

find_program(NUMDIFF_TOOL_PATH numdiff)

enable_testing()

if(NOT BENCHMARK_DIR_FOUND)
	return()
endif()

if(DEFINED ENV{CTEST_NUM_THREADS})
    set(NUM_CTEST_PROCESSORS $ENV{CTEST_NUM_THREADS})
elseif(DEFINED ENV{NUM_THREADS})
    set(NUM_CTEST_PROCESSORS $ENV{NUM_THREADS})
else()
    set(NUM_CTEST_PROCESSORS 3)
endif()

# Print benchmark config summary
if(BENCHMARK_REF_DIR_FOUND)
	if(NUMDIFF_TOOL_PATH)
		message(STATUS "Running benchmarks with numdiff file comparison.")
	else()
		message(STATUS "Running benchmarks with strict file comparison.")
	endif()
else()
	message(STATUS "Benchmark reference files not found (Specify non-default location with BENCHMARK_REF_DIR). No file comparison!")
endif()

set(BENCHMARK_TIMEOUT 3600 CACHE INTERNAL "") # in s, 60 minutes timeout on normal benchmarks
set(EXCEEDING_BENCHMARK_TIMEOUT 86400 CACHE INTERNAL "") # 1 day timeout on exceeding benchmarks
file(GLOB BENCHMARK_CONFIGS ${PROJECT_SOURCE_DIR}/scripts/cmake/benchmarks/*.cmake)
foreach(BENCHMARK_CONFIG ${BENCHMARK_CONFIGS})
	include(${BENCHMARK_CONFIG})
endforeach()

# Test targets
if(CMAKE_CONFIGURATION_TYPES)
	set(CONFIG_PARAMETER --build-config "$<CONFIGURATION>")
endif()

foreach(LABEL all short normal long exceeding)
	if(LABEL STREQUAL "all")
		set(TARGET_NAME benchmarks)
	else()
		set(TARGET_NAME benchmarks-${LABEL})
	endif()
	add_custom_target(${TARGET_NAME}
        COMMAND ${CMAKE_CTEST_COMMAND} -T Test
            --label-regex ${LABEL} --force-new-ctest-process --output-on-failure
            ${CONFIG_PARAMETER} --parallel ${NUM_CTEST_PROCESSORS}
        DEPENDS ogs
        USES_TERMINAL
	)
endforeach()
add_custom_target(benchmarks-short-normal
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
        --label-regex "short\\|normal" --force-new-ctest-process --output-on-failure
        ${CONFIG_PARAMETER} --parallel ${NUM_CTEST_PROCESSORS}
    DEPENDS ogs
    USES_TERMINAL
)
add_custom_target(benchmarks-short-normal-long
    COMMAND ${CMAKE_CTEST_COMMAND} -T Test
        --label-regex "short\\|normal\\|long" --force-new-ctest-process
        --output-on-failure
        ${CONFIG_PARAMETER} --parallel ${NUM_CTEST_PROCESSORS}
    DEPENDS ogs
    USES_TERMINAL
)

set_directory_properties(PROPERTIES
	ADDITIONAL_MAKE_CLEAN_FILES ${PROJECT_BINARY_DIR}/benchmarks
)
