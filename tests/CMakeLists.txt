if( MSVC ) # VS2012 doesn't support correctly the tuples yet
	add_definitions( /D _VARIADIC_MAX=10 )
endif()

add_subdirectory( gtest )

# Add tests which do not need testdata here:
set ( SOURCES
	testrunner.cpp
	testBase.cpp
)

# Add tests here if they need testdata
if(TESTDATA_DIR_FOUND)
	set(SOURCES ${SOURCES}
		GEO/TestPointInTriangle.cpp
		GEO/TestPointInPolygon.cpp
	)
endif()

if (NOT MSVC)
	set ( SOURCES ${SOURCES}
		testOgs.cxx
		testBmSkel.cxx)
endif ()

set ( SOURCES ${SOURCES}
	Matrix/testMatrix.cpp
    )

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/tests/gtest
)

add_executable (testrunner ${SOURCES})
target_link_libraries(testrunner
	Base
	FEM
	FileIO
	GEO
	GTest
	MathLib
	MSH
)

if (MSVC)
	set(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/testrunner)
	set(QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/QtTestrunner)
else ()
	set(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/testrunner)
	if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
		target_link_libraries( testrunner pthread )
	endif() # "${CMAKE_SYSTEM_NAME}" MATCHES "Linux"
	set (QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/QtTestrunner)
endif() # MSVC

# Add targets for running tests (with coverage reports)
if(OGS_COVERAGE)
	# Use cobertura output on Jenkins
	if(JENKINS_URL)
		SETUP_TARGET_FOR_COVERAGE_COBERTURA(gtests_cobertura testrunner "gtests_coverage_results")
		if(BENCHMARK_DIR_FOUND)
			SETUP_TARGET_FOR_COVERAGE_COBERTURA(benchmark_coverage ctest "benchmark_coverage_results" "-E;'Tests|FILE|EXCEED';-j;${PROCESSOR_COUNT}")
		endif() # BENCHMARK_DIR_FOUND
	# Use lcov output else
	else()
		SETUP_TARGET_FOR_COVERAGE(gtests testrunner "gtests_coverage_results")
		if(BENCHMARK_DIR_FOUND)
			SETUP_TARGET_FOR_COVERAGE(benchmark_coverage ctest "benchmark_coverage_results" "-E;'Tests|FILE|EXCEED';-j;${PROCESSOR_COUNT}")
		endif() # BENCHMARK_DIR_FOUND
	endif() # JENKINS_URL

else() # without coverage reports
	add_custom_target(gtests testrunner	DEPENDS testrunner)
	set_target_properties(testrunner gtests PROPERTIES FOLDER Testing)
endif() # OGS_COVERAGE

# Either collect all gtests in one ctest run of testrunner executable
# or run each gtest individually
if (OGS_ONE_BIG_GTEST)
	add_test (GTests ${testrunnerExe}  --gtest_output=xml:./)
else ()
	ADD_GOOGLE_TESTS ( ${testrunnerExe} ${SOURCES})
endif ()

# add new test models here not all of these models will be needed, depending
# on compile flags but it is safe to copy them all across
if (NOT MSVC)
	add_subdirectory( data/bmskel )
endif ()
