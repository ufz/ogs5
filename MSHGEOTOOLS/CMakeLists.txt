include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
)

# create library
add_library( MSHGEOTOOLS STATIC ExtractMeshNodes.h ExtractMeshNodes.cpp )

target_link_libraries (
	MSHGEOTOOLS
	${CMAKE_THREAD_LIBS_INIT}
	GEO
	MSH
)

