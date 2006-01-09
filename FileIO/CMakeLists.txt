set( HEADERS
	FEMIO/BoundaryConditionIO.h
	FEMIO/GeoIO.h
	FEMIO/ProcessIO.h
	MathIO/CRSIO.h
	MeshIO/LegacyVtkInterface.h
	MeshIO/OGSMeshIO.h
	OGSIOVer4.h
	readNonBlankLineFromInputStream.h
	StationIO.h
	Writer.h
	XmlIO/RapidXMLInterface.h
)

set( SOURCES
	FEMIO/BoundaryConditionIO.cpp
	FEMIO/GeoIO.cpp
	FEMIO/ProcessIO.cpp
	MeshIO/LegacyVtkInterface.cpp
	MeshIO/OGSMeshIO.cpp
	OGSIOVer4.cpp
	readNonBlankLineFromInputStream.cpp
	StationIO.cpp
	Writer.cpp
	XmlIO/RapidXMLInterface.cpp
)

if(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	include_directories(${Boost_INCLUDE_DIR})
endif()

# Create the library
add_library( FileIO STATIC ${HEADERS} ${SOURCES} )

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/GeoIO
	${CMAKE_SOURCE_DIR}/FileIO/MathIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/GEM
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
)

target_link_libraries( FileIO
	Base
	GEO
	FEM
)

if (Shapelib_FOUND)
	target_link_libraries( FileIO ${Shapelib_LIBRARIES} )
endif ()
