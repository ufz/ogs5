/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef rf_sfc_INC

#define rf_sfc_INC

#include "geo_ply.h"

//-------------------------------------------------------------------------
class CTriangle
{
public:
	long number;
	double x[3];
	double y[3];
	double z[3];
	long msh_numbers[3];
};
//-------------------------------------------------------------------------
class CTIN
{
public:
	CTIN() {}
	~CTIN();
	std::string name;
	std::vector<CTriangle*> Triangles;
};

//-------------------------------------------------------------------------
class Surface
{
public:
	Surface();
	~Surface();
	// ID
	long id; // CC
	std::string name;
	// Properties
	int type;
	std::string type_name;
	std::string data_name;
	// int data_type;
	double epsilon;
	double mesh_density;
	int mat_group; // MMP
	std::string mat_group_name;
	double Radius; // Radius of cylinder. WW
	// display
	int m_color[3];
	int display_mode_2d;
	int display_mode_3d;
	int display_mode_bc;
	bool highlighted; // CC
	// topology
	bool order;
	bool createtins;
	double center_point[3];
	// TIN
	CTIN* TIN;
	// point vector
	std::vector<CGLPoint*> polygon_point_vector;
	// polylines
	// list<CGLPolyline*> polyline_of_surface_list;
	std::vector<CGLPolyline*> polyline_of_surface_vector;
	std::vector<int> polyline_of_surface_orient;
	std::vector<double*> nodes_coor_vector;
	// MSH
	int meshing_allowed; // TK
	//----------------------------------------------------------------
	// Method
	// I/O
	void output(FILE* geo_file, int& p_index, int& l_index, int& pl_index, int& s_index);
	void Write(const std::string&);
	std::ios::pos_type Read(std::ifstream*);
	// Topology
	void PolylineOrientation(); // CC
	void ReArrangePolylineList();
	void PolygonPointVector(); // OK/CC
	// point
	void CalcCenterPoint(void);
	// display
	void AssignColor(); // CC
	// TIN
	void CreateTIN(void);
	void ReadTIN(const std::string&); // CC
	void WriteTIN(const std::string&); // CC
	void WriteTINTecplot(const std::string&); // CC
	// Tecplot
	void WriteTecplot(std::fstream*);
	bool PointInSurface(CGLPoint*); // OK
	// material
	long profile_code; // YD

private:
	//
	friend class CGLLine; // WW
};
// vector
extern std::vector<Surface*> surface_vector; // CC
extern std::vector<Surface*> GetSurfaceVector(void); // CC
extern void GEOCreateSurfacePointVector(void); // CC
// Access
extern Surface* GEOGetSFCByName(const std::string&);
// I/O
extern void GEOReadSurfaces(const std::string& file_name_path_base);
extern void GEOWriteSurfaces(const std::string&); // C
// Remove
extern void GEORemoveAllSurfaces(); // CC
extern void GEORemoveSurface(long); // CC
extern void GEORemoveSFC(Surface* m_sfc);
// Topology
extern void GEOSurfaceTopology(void);
extern void GEOUnselectSFC(); // OK
// TIN
#define TIN_FILE_EXTENSION ".tin"

extern int sfc_ID_max;
// MSH
void MSHUnselectSFC(); // OK
extern int SFCGetMaxMATGroupNumber(); // OK
extern void SFCAssignMATColors(); // OK
extern void SFCAssignColors(); // OK

#endif
