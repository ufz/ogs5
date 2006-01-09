set( HEADERS
	AxisAlignedBoundingBox.h
	BruteForceClosestPair.h
	ClosestPair.h
	Color.h
	geo_dom.h
	geo_lib.h
	geo_lin.h
	geo_mathlib.h
	geo_ply.h
	geo_pnt.h
	geo_sfc.h
	geo_vol.h
	GeoObject.h
	GEOObjects.h
	GeoType.h
	Grid.h
	Point.h
	PointVec.h
	PointWithID.h
	Polygon.h
	Polyline.h
	PolylineVec.h
	PolylineWithSegmentMarker.h
	ProjectData.h
	PropertyBounds.h
	QuadTree.h
	Raster.h
	SensorData.h
	SimplePolygonTree.h
	Station.h
	Surface.h
	SurfaceGrid.h
	SurfaceVec.h
	TemplatePoint.h
	TemplateVec.h
	Triangle.h
)

set( SOURCES
	AxisAlignedBoundingBox.cpp
	BruteForceClosestPair.cpp
	Color.cpp
	geo_dom.cpp
	geo_lib.cpp
	geo_lin.cpp
	geo_mathlib.cpp
	geo_ply.cpp
	geo_pnt.cpp
	geo_sfc.cpp
	geo_vol.cpp
	GEOObjects.cpp
	GeoType.cpp
	Point.cpp
	PointVec.cpp
	Polygon.cpp
	Polyline.cpp
	PolylineWithSegmentMarker.cpp
	ProjectData.cpp
	Raster.cpp
	SensorData.cpp
	SimplePolygonTree.cpp
	Station.cpp
	Surface.cpp
	SurfaceGrid.cpp
	Triangle.cpp
)

# Create the library
add_library( GEO STATIC ${SOURCES} )


include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
	${CMAKE_SOURCE_DIR}/MSH
)


target_link_libraries (
	GEO
	Base
	FEM
	MathLib
)

