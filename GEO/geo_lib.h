/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Control Level
   Task:
   Programing:
   07/2003 TK/OK Implementation
   11/2003 WW Geo_Read class
   11/2003 OK -> GeoLib class
   09/2005 CC GeoLib2
**************************************************************************/

#ifndef geolib_INC
#define geolib_INC

enum GEO_TYPE {GS_POINT,GS_POLYLINE,GS_SURFACE,GS_VOLUME};

// C++ STL
#include <string>

class GeoLib
{
public:
	GeoLib() {}
	~GeoLib() {}
private:
	friend class CGLPolyline;
	friend class Surface;
};

extern void GEOLIB_Read_GeoLib (const std::string &file_name_path_base);
extern void GEOLIB_Clear_GeoLib_Data ();
//OK
#define GLI_FILE_EXTENSION ".gli"

extern void GEO_Delete_DoublePoints(); // TODO: Implement new function from TK
extern void GEO_Serialize_Point_Numbers();
extern double GEO_Get_Min_PolySeg_length();
extern double GEO_Get_Max_PolySeg_length();
extern void GEO_Set_Poly_Seg_Length(double min_seg_length, double max_seg_length);
extern void GEO_Write_GMSH_Input_File(const char* file_name);
extern void GEOWrite(std::string); //OK4.1
extern void GEO_Get_Min_Max_Distance_of_polyline_neighbor_points();
extern void GEO_Copy_Min_Distance_Of_Neighbor_Points_To_Mesh_Density();
extern void GEO_Get_Min_Distance_of_neighbor_points();
extern void GEO_Mesh_Density_BiSect();
#endif
