##-------------------------------------------------------
##
##    Tool to convert GIS shape data into FEM data
##                By Wenqing Wang
##
##                     03.2010
##-------------------------------------------------------

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

set(GIS2FEM_SOURCES main.cpp misc.h misc.cpp )

add_executable(fem_toolkits ${GIS2FEM_SOURCES})

target_link_libraries( fem_toolkits
	FileIO
	MSH
	MSHGEOTOOLS
)
