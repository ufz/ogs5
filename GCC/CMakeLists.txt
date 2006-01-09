set( HEADERS
	Density.h
	DiffAL.h
	DiffOH.h
	eqbrm.h
	HKF.h
	IAPWS-IF97.h
	IO.h
	logk.h
	NR.h
	PITZdata.h
	species.h
	VLE.h
	WATEQdata.h
)

set ( SOURCES
	Density.cpp
	DiffOH.cpp
	HKF.cpp
	IAPWS-IF97.cpp
	IO.cpp
	logk.cpp
	NR.cpp
	PITZdata.cpp
	VLE.cpp
)

# Create the library
add_library( GCC STATIC ${HEADERS} ${SOURCES} )
target_link_libraries( GCC BrentLib )
include_directories( ${CMAKE_SOURCE_DIR}/GCC )

if (WIN32)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN")
endif ()

if (UNIX)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIPMGEMPLUGIN -D__unix")
endif ()
