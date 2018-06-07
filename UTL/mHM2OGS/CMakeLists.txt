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

FILE(GLOB mHM2OGS_H RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR}   *.h)
FILE(GLOB mHM2OGS_CPP RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_executable(mHM2OGS ${mHM2OGS_H} ${mHM2OGS_CPP})

target_link_libraries( mHM2OGS
	FileIO
	MSH
	MSHGEOTOOLS
)

if (WIN32)
	if (MSVC)
		target_link_libraries( mHM2OGS
			winmm.lib # for timer
		)
	endif()
endif()
