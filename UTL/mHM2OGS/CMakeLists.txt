##-------------------------------------------------------
##
##    Tool to process mHM data for OGS
##        (Based on a function in GIS2FEM (03.2010))
##                By Wenqing Wang
##
##                     05.2018
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

set(mHM2OGS_SOURCES main.cpp misc.h misc.cpp )

add_executable(mHM2OGS ${mHM2OGS_SOURCES})

target_link_libraries( mHM2OGS
	FileIO
	MSH
	MSHGEOTOOLS
)
