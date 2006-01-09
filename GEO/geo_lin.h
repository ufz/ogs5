/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef gl_lin_INC
#define gl_lin_INC

#include "geo_pnt.h"

/*---------------------------------------------------------------*/
class CGLLine
{
private:
	bool marked; // For mesh
	friend class Surface; //WW
	// kg44 needed for domain decomposition with PETSC
	bool for_ic;

public:
	CGLLine(void);
	~CGLLine(void);
	CGLPoint* m_point1;
	CGLPoint* m_point2;
	std::string name;
	double value;
	long mesh_index, gli_line_id; // mesh_index: unique index for mesh
	long point1, point2;
	int orientation;
	double epsilon;
	long* msh_nodes;
	long no_msh_nodes;
	int mat_group;
	int display_mode;
	//MSH
	std::vector<double*> nodes_coor_vector;
	// kg44 needed for domain decomposition with PETSC
	void SetConditionTypeIC(const bool value) { for_ic=value;} ;
	bool GetConditionTypeIC() const {return for_ic;} ;
	//Method
	CGLLine* GEOGetLine(long);
	CGLLine* CheckLineOutPut();
	CGLLine* Exists();
	// void SetRFIPointsClose();//MSH + geomathlib
	//void CreateMSHLines(void);//MSH
};

extern std::vector<CGLLine*> GEOLIB_GetGLILines_Vector(void);
extern std::vector<CGLLine*> gli_lines_vector;
extern std::vector<CGLLine*> gli_file_lines_vector;
extern void Clear_LineVector();

#endif
