# check 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	set( HAVE_64_BIT 0 )
	set( BITS 32 )
else()
	set( HAVE_64_BIT 1 )
	add_definitions(-DHAVE_64_BIT)
	set( BITS 64)
endif()

# Convert environment variables
if(NOT $ENV{LIBRARIES_DIR} STREQUAL "")
	string(REGEX REPLACE "\\\\" "/" LIBRARIES_DIR $ENV{LIBRARIES_DIR})
endif()

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
