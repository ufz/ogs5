
include_directories(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
)

# Create executables
add_executable( testPointInTriangle testPointInTriangle.cpp )

target_link_libraries( testPointInTriangle
	FEM
        MSH
       	GEO
        MathLib
       	FileIO
        Base
	${QT_LIBRARIES}
)
