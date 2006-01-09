# check 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	set( HAVE_64_BIT 0 )
	set( BITS 32 )
else()
	set( HAVE_64_BIT 1 )
	add_definitions(-DHAVE_64_BIT)
	set( BITS 64)
endif()

# Visual Studio detection
if(${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
	OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")

	set(VS32 TRUE)
	set(VS64 FALSE)
	message(STATUS "Generator: Visual Studio 32 Bit")

endif()

if(${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")

	set(VS32 FALSE)
	set(VS64 TRUE)
	message(STATUS "Generator: Visual Studio 64 Bit")

endif()

# Convert environment variables
if(NOT $ENV{LIBRARIES_DIR} STREQUAL "")
	string(REGEX REPLACE "\\\\" "/" LIBRARIES_DIR $ENV{LIBRARIES_DIR})
endif()

macro(COPY_FILE_IF_CHANGED in_file out_file target)
	if(${in_file} IS_NEWER_THAN ${out_file})
		#message("Copying file: ${in_file} to: ${out_file}")
		add_custom_command (
				TARGET     ${target}
				POST_BUILD
				COMMAND    ${CMAKE_COMMAND}
				ARGS       -E copy ${in_file} ${out_file}
		)
	endif()
endmacro()

macro(COPY_FILE_INTO_DIRECTORY_IF_CHANGED in_file out_dir target)
		get_filename_component(file_name ${in_file} NAME)
		COPY_FILE_IF_CHANGED(${in_file} ${out_dir}/${file_name}
${target})
endmacro()

#Copies all the files from in_file_list into the out_dir.
# sub-trees are ignored (files are stored in same out_dir)
macro(COPY_FILES_INTO_DIRECTORY_IF_CHANGED in_file_list out_dir target)
	foreach(in_file ${in_file_list})
				COPY_FILE_INTO_DIRECTORY_IF_CHANGED(${in_file}
${out_dir} ${target})
		endforeach()
endmacro()

macro(COPY_FILE_INTO_EXECUTABLE_DIRECTORY in_file target)
	if(WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}/Debug
			${target}
		)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			"${CMAKE_CURRENT_SOURCE_DIR}/${in_file}"
			${EXECUTABLE_OUTPUT_PATH}/Release
			${target}
		)
	else()
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}
			${target}
		)
	endif()
endmacro()

# Adds a benchmark run.
# authorName Your short name
# benchmarkName Relative path in benchmarks directory
# ogsConfiguration E.g. "OGS_FEM"
# numProcesses Number of processes for mpirun, Please set to 1 for non-MPI benchmarks
# Additional arguments add output files to compare
function (ADD_BENCHMARK authorName benchmarkName ogsConfiguration numProcesses)

  set(CONFIG_MATCH FALSE)
  if("${ogsConfiguration}" STREQUAL "OGS_FEM" AND OGS_FEM )
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_SP" AND OGS_FEM_SP)
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_GEMS" AND OGS_FEM_GEMS)
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_BRNS" AND OGS_FEM_BRNS)
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_PQC" AND OGS_FEM_PQC)
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_IPQC" AND OGS_FEM_IPQC)
	set(CONFIG_MATCH TRUE)
  endif()

  if("${ogsConfiguration}" STREQUAL "OGS_FEM_CAP" AND OGS_FEM_CAP)
    set(CONFIG_MATCH TRUE)
  endif()

  if (UNIX) # Only supported on Linux
	if("${ogsConfiguration}" STREQUAL "OGS_FEM_LIS" AND OGS_FEM_LIS)
	  set(CONFIG_MATCH TRUE)
	endif()
	if("${ogsConfiguration}" STREQUAL "OGS_FEM_MKL" AND OGS_FEM_MKL)
	  set(CONFIG_MATCH TRUE)
	endif()
	if("${ogsConfiguration}" STREQUAL "OGS_FEM_MPI" AND OGS_FEM_MPI)
	  set(CONFIG_MATCH TRUE)
	endif()
	if("${ogsConfiguration}" STREQUAL "OGS_FEM_PETSC" AND OGS_FEM_PETSC)
	  set(CONFIG_MATCH TRUE)
	endif()
	if("${ogsConfiguration}" STREQUAL "OGS_FEM_PETSC_GEMS" AND OGS_FEM_PETSC_GEMS)
	  set(CONFIG_MATCH TRUE)
	endif()
  endif ()

  if (CONFIG_MATCH)
	if(WIN32)
	  set(ogsExe ${EXECUTABLE_OUTPUT_PATH}/Release/ogs)
	else()
	  set(ogsExe ${EXECUTABLE_OUTPUT_PATH}/ogs)
	endif()

	# Set timeout
	string(REGEX MATCH "EXCEEDING" benchmarkExceeding ${benchmarkName})
	if(benchmarkExceeding)
		set(THIS_BENCHMARK_TIMEOUT ${EXCEEDING_BENCHMARK_TIMEOUT})
	else()
		set(THIS_BENCHMARK_TIMEOUT ${BENCHMARK_TIMEOUT})
	endif()

	string (REGEX MATCH "[^/]+$" benchmarkStrippedName ${benchmarkName})
	string (LENGTH ${benchmarkName} benchmarkNameLength)
	string (LENGTH ${benchmarkStrippedName} benchmarkStrippedNameLength)
	math (EXPR substringLength ${benchmarkNameLength}-${benchmarkStrippedNameLength})
	string (SUBSTRING ${benchmarkName} 0 ${substringLength} benchmarkDir)
	string (REPLACE "/" "_" benchmarkNameUnderscore ${benchmarkName})
	string (REPLACE "_LONG_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})
	string (REPLACE "_EXCEEDING_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})

	# Delete output files on testing
	foreach (entry ${ARGN})
	set (FILES_TO_DELETE "${FILES_TO_DELETE} ${entry}")
	endforeach ()
	if(OGS_PROFILE)
		set(FILES_TO_DELETE "${FILES_TO_DELETE} gmon.out")
	endif()

	# Adds a benchmark run. This calls AddTest.cmake to execute several steps.
	add_test(
		${authorName}_BENCHMARK_${benchmarkName}
		${CMAKE_COMMAND}
		-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
		-DEXECUTABLE_OUTPUT_PATH=${EXECUTABLE_OUTPUT_PATH}
		-DbenchmarkStrippedName=${benchmarkStrippedName}
		-DbenchmarkDir=${benchmarkDir}
		-DFILES_TO_DELETE=${FILES_TO_DELETE}
		-DOGS_PROFILE=${OGS_PROFILE}
		-DOGS_OUTPUT_PROFILE=${OGS_OUTPUT_PROFILE}
		-DGPROF_PATH=${GPROF_PATH}
		-DDOT_TOOL_PATH=${DOT_TOOL_PATH}
		-DBENCHMARK_TIMEOUT=${THIS_BENCHMARK_TIMEOUT}
		-DOGS_FEM_CONFIG=${ogsConfiguration}
		-DNUM_PROCESSES=${numProcesses}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/AddBenchmark.cmake
	)

	# compare file differences with python script
	if(PYTHONINTERP_FOUND)
		file (REMOVE ${PROJECT_SOURCE_DIR}/../benchmarks/results/temp/temp_${benchmarkNameUnderscore}.txt)
		foreach (entry ${ARGN})
			file (APPEND ${PROJECT_SOURCE_DIR}/../benchmarks/results/temp/temp_${benchmarkNameUnderscore}.txt "${entry}\n")
		endforeach ()
		add_test(
			${authorName}_FILECOMPARE_${benchmarkName}
			${CMAKE_COMMAND} -E chdir ${PROJECT_SOURCE_DIR}/../benchmarks/results
			${PYTHON_EXECUTABLE}
			${PROJECT_SOURCE_DIR}/scripts/compare.py
			temp/temp_${benchmarkNameUnderscore}.txt
			../../benchmarks_ref/
			${authorName}_${benchmarkNameUnderscore}.html
			../
		)
	endif()

  # copy benchmark output files to reference directory
  if(COPY_BENCHMARKS_TO_REF)
	foreach (entry ${ARGN})
	  configure_file( ${PROJECT_SOURCE_DIR}/../benchmarks/${entry} ${PROJECT_SOURCE_DIR}/../benchmarks_ref/${entry} COPYONLY)
	endforeach ()
  endif()

  endif()

endfunction ()

# Checks if a valid ogs configuration is given
function(CHECK_CONFIG)

	set(configs
		"${OGS_FEM}"
		"${OGS_FEM_SP}"
		"${OGS_FEM_MPI}"
		"${OGS_FEM_GEMS}"
		"${OGS_FEM_BRNS}"
		"${OGS_FEM_MKL}"
		"${OGS_FEM_PQC}"
		"${OGS_FEM_IPQC}"
		"${OGS_FEM_LIS}"
		"${OGS_FEM_CHEMAPP}"
		"${OGS_FEM_PETSC}"
		"${OGS_FEM_PETSC_GEMS}"
        "${OGS_FEM_CAP}")
	set(counter 0)

	foreach(config ${configs})
		if(config)
			math(EXPR counter "${counter} + 1")
		endif()
#		message(STATUS "config test ${config} found total ${counter}")
	endforeach()
	if(counter EQUAL 0)
		message(STATUS "No configuration specified. Assuming default configuration...")
		set(OGS_FEM ON)
                set(counter 1)
	endif()

	if(counter GREATER 1)
		message(FATAL_ERROR "Error: More than one OGS configuration given (${counter}). Please use only one of the following configurations:
			OGS_FEM (Default FEM configuration)
			OGS_FEM_SP
			OGS_FEM_MPI
			OGS_FEM_GEMS
			OGS_FEM_BRNS
			OGS_FEM_MKL
			OGS_FEM_PQC
			OGS_FEM_IPQC
			OGS_FEM_LIS
			OGS_FEM_CHEMAPP
			OGS_FEM_PETSC
			OGS_FEM_PETSC_GEMS
			OGS_FEM_CAP")
	endif()

endfunction()

# Creates one ctest for each googletest found in source files passed as arguments
# number two onwards. Argument one specifies the testrunner executable.
macro(ADD_GOOGLE_TESTS executable)
	foreach ( source ${ARGN} )
		file(READ "${source}" contents)
		string(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
		foreach(hit ${found_tests})
			string(REGEX REPLACE ".*\\(([A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+)\\).*" "\\1.\\2" test_name ${hit})
			add_test(${test_name} ${executable}  --gtest_output=xml --gtest_filter=${test_name} ${MI3CTestingDir})
			# message ("Adding test: ${test_name}")
		endforeach()
	endforeach()
endmacro()

# copies the model files to the build dir and adds them as targets so that
# the build files are re-built when the source model files change
macro ( UPDATE_MODEL_FILES dirOUT fileLIST )
	get_filename_component( _tdir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
	#message (STATUS "Copying files to ${dirOUT} from ${fileLIST}.\n")
	foreach ( _file1 ${${fileLIST}} )
		set( _file ${CMAKE_CURRENT_SOURCE_DIR}/${_file1} )
		get_filename_component( _fdest ${_file} NAME )
		set( dest ${dirOUT}/${_fdest} )
		#message( STATUS "Copying ${_file} to ${dest} \n" )

		#message (STATUS "Adding targets ${_tdir}.${_fdest}\n").\n")
		add_custom_target( ${_tdir}.${_fdest}
			${CMAKE_COMMAND} -E copy_if_different
			${_file} ${dest})

		add_dependencies( testrunner ${_tdir}.${_fdest} )
	endforeach()
endmacro ()
