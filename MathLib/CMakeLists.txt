set( HEADERS
	AnalyticalGeometry.h
	EarClippingTriangulation.h
	InterpolationAlgorithms/CubicSpline.h
	InterpolationAlgorithms/InverseDistanceInterpolation.h
	InterpolationAlgorithms/LinearIntervalInterpolation.h
	InterpolationAlgorithms/PiecewiseLinearInterpolation.h
	LinAlg/GaussAlgorithm.h
	LinAlg/TriangularSolve.h
	LinkedTriangle.h
	MathTools.h
	Matrix.h
	max.h
	Vector3.h
)

set( SOURCES
	AnalyticalGeometry.cpp
	EarClippingTriangulation.cpp
	InterpolationAlgorithms/CubicSpline.cpp
	InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp
	LinkedTriangle.cpp
	MathTools.cpp
)

if(OGS_LSOLVER STREQUAL PETSC)
	set( SOURCES ${SOURCES} PETSC/PETScLinearSolver.h PETSC/PETScLinearSolver.cpp)
endif()

add_library( MathLib STATIC ${HEADERS} ${SOURCES} )

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
	${CMAKE_SOURCE_DIR}/MathLib/PETSC
)

target_link_libraries (
	MathLib
	Base
)

