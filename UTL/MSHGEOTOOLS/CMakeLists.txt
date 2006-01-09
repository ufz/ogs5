
include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
	${CMAKE_SOURCE_DIR}/UTL/MSHGEOTOOLS
)

# Create executables
add_executable( ExtractMeshNodes mainExtractMeshNodes.cpp )
add_executable( ExtractMeshNodeIDs mainExtractMeshNodeIDs.cpp )
add_executable( ModifyMeshProperties
	mainModifyMeshProperties.cpp
	ModifyMeshProperties.h
	ModifyMeshProperties.cpp )
add_executable( filterMeshNodes filterMeshNodes.cpp )
add_executable( convertGLIVerticalSurfaceToPolygon mainConvertGLIVerticalSurfaceToPolygon.cpp )

set_target_properties(ExtractMeshNodeIDs ExtractMeshNodes ModifyMeshProperties filterMeshNodes convertGLIVerticalSurfaceToPolygon
  PROPERTIES FOLDER Utilities)

target_link_libraries( ExtractMeshNodeIDs
	FEM
	FileIO
	GEO
	MSH
	MSHGEOTOOLS
)

target_link_libraries( ExtractMeshNodes
	FEM
	FileIO
	GEO
	MSH
	MSHGEOTOOLS
)


target_link_libraries( convertGLIVerticalSurfaceToPolygon
	FEM
	FileIO
	GEO
	MSH
	MSHGEOTOOLS
)


target_link_libraries( ModifyMeshProperties
	FEM
	FileIO
	GEO
	MathLib
	MSH
	MSHGEOTOOLS
)

target_link_libraries( filterMeshNodes
	Base
	MSH
)

add_executable( testMeshSearchAlgorithms testMeshSearchAlgorithms.cpp )
target_link_libraries( testMeshSearchAlgorithms
	Base
	FileIO
	GEO
	MSH
	MSHGEOTOOLS
)
