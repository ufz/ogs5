/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Volume
   Task:
   Programing:
   07/2003 OK/TK Implementation
   09/2005 CC GeoLib2
**************************************************************************/
#ifndef gs_vol_INC
#define gs_vol_INC
#include "geo_sfc.h"

class CGLVolume
{
	//----------------------------------------------------------------------
	// Data
private:
public:
	// ID
	std::string name;
	int type;
	std::string type_name;
	// GEO
	int display_mode;
	int data_type;
	Surface* m_sfc;
	std::vector<Surface*> surface_vector; // CC
	// list<string> surface_name_list;// todo CC
	int layer; // OK
	// MAT
	int mat_group;
	std::string mat_group_name;
	bool selected;
	// kg44 needed for domain decomposition with PETSC
	bool for_ic;
	//----------------------------------------------------------------------
	// Methods
private:
public:
	// Construction
	CGLVolume(void);
	~CGLVolume(void);
	// Access
	long GEOInsertVolume2List(CGLVolume*);
	// By PCH
	void AddSurface(Surface* P_Surface);
	// I/O
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::string); // CC
	// GEO
	bool PointInVolume(CGLPoint*, int);
	// kg44 needed for domain decomposition with PETSC
	void SetConditionTypeIC(const bool value) { for_ic = value; };
	bool GetConditionTypeIC() const { return for_ic; };
	//----------------------------------------------------------------------
};

extern std::vector<CGLVolume*> volume_vector; // CC
extern std::vector<CGLVolume*> GEOGetVolumes(void);
// I/O
extern void GEOReadVolumes(std::string file_name_path_base);
extern void GEOWriteVolumes(std::string); // CC
void GEOWriteVOL(FILE*);
// Access
extern CGLVolume* GEOGetVOL(std::string); // CC//OK
extern CGLVolume* GetVolume(long nsel); // CC
// Remove
extern void GEORemoveAllVolumes();
extern void GEORemoveVOL(CGLVolume*); // OK
#endif
