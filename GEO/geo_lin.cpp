/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib-Method: CGLLine
   Task:
   Programing:
   01/2004 OK Implementation
   09/2005 CC GeoLib2
**************************************************************************/

#include "files0.h"
#include "geo_lin.h"
#include "geo_mathlib.h" //CC
#include "geo_pnt.h"

using namespace std;

vector<CGLLine*> gli_lines_vector;
vector<CGLLine*> gli_file_lines_vector;

/**************************************************************************
   GeoLib-Method: CGLLine
   Task:
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CGLLine::CGLLine(void)
{
	mesh_index = -1;
	gli_line_id = -1; // mesh_index: unique index for mesh
	marked = false;
	m_point1 = NULL;
	m_point2 = NULL;
	point1 = -1;
	point2 = -1;
	orientation = -1;
	epsilon = 0.1;
	no_msh_nodes = 0;
	msh_nodes = NULL;
	mat_group = 0;
}

// deconstructor
CGLLine::~CGLLine(void)
{
	int i = 0;
	for (i = 0; i < (int)nodes_coor_vector.size(); i++)
		delete nodes_coor_vector[i];
	nodes_coor_vector.clear();
}

/*************************************************************************
   GeoLib- Funktion: GEOLIB_GetGLILines_Vector

   Aufgabe: Pointer f�r externen Vektorenzugriff

   Programmaenderungen:
   08/2003     TK        Erste Version
 **************************************************************************/
vector<CGLLine*> GEOLIB_GetGLILines_Vector(void)
{
	return gli_lines_vector;
}
/*************************************************************************
   GeoLib- Funktion: Clear_LineVector

   Aufgabe: Leeren des GLI Line-Vektors

   Programmaenderungen:
   08/2003     TK        Erste Version
 **************************************************************************/
void Clear_LineVector()
{
	gli_lines_vector.clear();
}

/**************************************************************************
   GeoLib-Method: GEOGetLine
   Task: select volume instance from list by name
   Programing:
   07/2003 OK Implementation
**************************************************************************/
CGLLine* CGLLine::GEOGetLine(long number)
{
	return gli_lines_vector[number];
}

CGLLine* CGLLine::Exists()
{
	long i = 0;
	long node1, node2, node3, node4;
	long line_vector_length = (long)gli_lines_vector.size();

	for (i = 0; i < line_vector_length; i++)
	{
		node1 = gli_lines_vector[i]->point1;
		node2 = gli_lines_vector[i]->point2;
		node3 = point1;
		node4 = point2;
		if ((node1 == node3) && (node2 == node4))
		{
			gli_line_id = gli_lines_vector[i]->gli_line_id;
			m_point1 = gli_lines_vector[i]->m_point1;
			m_point2 = gli_lines_vector[i]->m_point2;
			point1 = gli_lines_vector[i]->point1;
			point2 = gli_lines_vector[i]->point2;
			orientation = 1;
			return gli_lines_vector[i];
		}
		// if ( ((node1==point1)&&(node2==point2))||((node2==point1)&&(node1==point2)) )
		// return  gli_lines_vector[i];
	}
	return NULL;
}

/**************************************************************************
   GeoLib-Method:
   Task: Check if the line is output in order to avoid the multi-output
   Programing:
   03/2004 WW Implementation
**************************************************************************/
CGLLine* CGLLine::CheckLineOutPut()
{
	const double DistTol = 1.0e-10;
	CGLLine* CGLn = NULL;

	int Size = (long)gli_lines_vector.size();
	for (int j = 0; j < Size; j++)
	{
		CGLn = gli_lines_vector[j];
		if (CGLn->marked) // Already being output

			if (((m_point1->PointDis(CGLn->m_point1) < DistTol) && (m_point2->PointDis(CGLn->m_point2) < DistTol))
			    || ((m_point1->PointDis(CGLn->m_point2) < DistTol) && (m_point2->PointDis(CGLn->m_point1) < DistTol)))
				return CGLn;
	}
	return NULL;
}
