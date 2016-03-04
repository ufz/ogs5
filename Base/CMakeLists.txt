set( HEADERS
	binarySearch.h
	break.h
	BuildInfo.h
	DateTools.h
	display.h
	FileFinder.h
	FileTools.h
	makros.h
	memory.h
	MemWatch.h
	printList.h
	quicksort.h
	StringTools.h
	swap.h
	timer.h
	uniqueListInsert.h
	wait.h
)

set( SOURCES
	binarySearch.cpp
	break.cpp
	${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp
	DateTools.cpp
	display.cpp
	FileTools.cpp
	makros.cpp
	memory.cpp
	MemWatch.cpp
	StringTools.cpp
	timer.cpp
)

# Create the library
add_library( Base STATIC ${HEADERS} ${SOURCES} )

set_target_properties(Base PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
)

if(OGS_BUILD_CORRECT_SCIENTIFIC)
	# Create utility executable
	add_executable( correctScientificNotation CorrectScientificNotationMain.cpp )
endif()
