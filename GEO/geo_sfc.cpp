/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Surface
   Task:
   Programing:
   03/2003 WW Implementation
   11/2003 WW Pointers to the polylines instead of names of polylines are set
           as the member of surface class
   11/2003 OK Get function
   12/2003 CC Implementation/modification
   03/2004 OK TINs
   09/2005 CC GeoLib2
**************************************************************************/
// MFC
#include <cstdlib>
// GEOLib
#include "files0.h"
#include "geo_lib.h"
#include "geo_mathlib.h"
#include "geo_ply.h"
#include "geo_pnt.h"
#include "geo_sfc.h"

// File path. 11.08.2011 WW
#include "makros.h"

using namespace std;

// GSP
/*----------------------------------------------------------------------*/
// vector
std::vector<Surface*> surface_vector; // CC
int sfc_ID_max = 0; // OK
/*constructor*/
Surface::Surface(void) : Radius(0.0)
{
	order = false;
	type = -1;
	m_color[0] = 0;
	m_color[1] = 0;
	m_color[2] = 0;
	mat_group = -1;
	TIN = NULL;
	createtins = false;
	epsilon = 1e-5; // OK
	highlighted = false; // CC
	meshing_allowed = 0; // OK
	data_name = "DOMAIN";
	type_name = "POLYLINE";
	mesh_density = 100;
	id = sfc_ID_max;
	sfc_ID_max++;
}

/* Destructor*/
Surface::~Surface()
{
	int i;

	while (polyline_of_surface_vector.size()) // CC
		polyline_of_surface_vector.pop_back();

	for (i = 0; i < (int)nodes_coor_vector.size(); i++)
		delete nodes_coor_vector[i];
	nodes_coor_vector.clear();
}
CTIN::~CTIN()
{
	int i;
	int no_triangles = (int)Triangles.size();
	// name.empty();
	name.clear();
	for (i = 0; i < no_triangles; i++)
		delete Triangles[i];
	while (Triangles.size())
		Triangles.pop_back();
	// this == NULL;
}
/**************************************************************************
   GeoLib-Method: Surface::output(FILE* geo_file, const int index)
   Task: Output the information of surface to a geo syntax file
   Programing:
   11/2003 WW Implement
   03/2004 WW Avoid multiplie defined entities
   08/2004 CC see version: Geosys3909_ac10ww.zip
   07/2005 CC .gli to .geo for multi polylines of one surface
**************************************************************************/
void Surface::output(FILE* geo_file, int& p_index, int& l_index, int& pl_index, int& s_index)
{
	int orient = 1;
	int i, j, Size0, Size, count, ortL0 = 0, ortL;
	double dist = 0.0;
	const double DistTol = 1.0e-10;
	// bool firstPline = true;
	std::vector<int> lineIndex;

	CGLPolyline* p_pline = NULL;

	CGLPoint* CGPnt = NULL;
	CGLPoint* CGPnt0 = NULL;
	CGLLine* CGLn0 = NULL;
	// WW CGLLine *CGLn1 = NULL;
	CGLLine* CGLn = NULL;

	// First entity: points
	// list<CGLPolyline*>::const_iterator p  = polyline_of_surface_list.begin();//CC
	std::vector<CGLPolyline*>::iterator p = polyline_of_surface_vector.begin();

	FilePrintString(geo_file, "// Points");
	LineFeed(geo_file);
	while (p != polyline_of_surface_vector.end()) // CC
	{
		p_pline = *p;
		Size = (int)p_pline->point_vector.size();
		for (i = 0; i < Size; i++)
		{
			CGPnt = p_pline->point_vector[i];

			// Unique the indeces
			if (CGPnt->GetIndex() < 0)
			{
				// Compare with other points
				Size0 = (int)gli_points_vector.size();
				for (j = 0; j < Size0; j++)
				{
					CGPnt0 = gli_points_vector[j];
					dist = CGPnt->PointDis(CGPnt0);
					if (dist < DistTol && CGPnt0->GetIndex() >= 0)
					{
						CGPnt->SetIndex(CGPnt0->GetIndex());
						break;
					}
				}
				// If the mesh index of the point is still not prescribed
				if (CGPnt->GetIndex() < 0)
				{
					p_index++;
					CGPnt->SetIndex(p_index);
					FilePrintString(geo_file, "Point(");
					FilePrintInt(geo_file, p_index);
					FilePrintString(geo_file, ")= {");
					FilePrintDouble(geo_file, CGPnt->x);
					FilePrintString(geo_file, ", ");
					FilePrintDouble(geo_file, CGPnt->y);
					FilePrintString(geo_file, ", ");
					FilePrintDouble(geo_file, CGPnt->z);
					FilePrintString(geo_file, ", ");
					FilePrintDouble(geo_file, CGPnt->mesh_density);
					// FilePrintString(geo_file, "DensityScale");
					FilePrintString(geo_file, "};");
					LineFeed(geo_file);
				}
			}
		}
		p++;
	}
	// Second entity: lines
	FilePrintString(geo_file, "// Lines");
	LineFeed(geo_file);
	p = polyline_of_surface_vector.begin(); // CC
	while (p != polyline_of_surface_vector.end()) // CC
	{
		p_pline = *p;
		std::vector<CGLLine*>::iterator pl // CC
		    = p_pline->line_vector.begin();
		while (pl != p_pline->line_vector.end()) // CC
		{
			CGLn = *pl;
			// OK test 11.09.2004
			if (!CGLn->m_point1 && !CGLn->m_point2)
			{
				std::cout << "Error in Surface::output: no line points"
				          << "\n";
				FilePrintString(geo_file, "Error in Surface::output: no line points");
				return;
			}
			if (CGLn->mesh_index < 0)
			{
				Size0 = (int)gli_lines_vector.size();
				for (j = 0; j < Size0; j++)
				{
					CGLn0 = gli_lines_vector[j];
					if (CGLn0->mesh_index >= 0)
						if (((CGLn->m_point1->PointDis(CGLn0->m_point1) < DistTol)
						     && (CGLn->m_point2->PointDis(CGLn0->m_point2) < DistTol))
						    || ((CGLn->m_point1->PointDis(CGLn0->m_point2) < DistTol)
						        && (CGLn->m_point2->PointDis(CGLn0->m_point1) < DistTol)))
						{
							CGLn->mesh_index = CGLn0->mesh_index;
							break;
						}
				}

				// If the mesh index of the line is still not prescribed
				if (CGLn->mesh_index < 0)
				{
					l_index++;
					CGLn->mesh_index = l_index;
					FilePrintString(geo_file, "Line(");
					FilePrintLong(geo_file, l_index);
					FilePrintString(geo_file, ")= {");
					FilePrintInt(geo_file, CGLn->m_point1->GetIndex());
					FilePrintString(geo_file, ", ");
					FilePrintInt(geo_file, CGLn->m_point2->GetIndex());
					FilePrintString(geo_file, "};");
					LineFeed(geo_file);
					CGLn->marked = true;
				}
			}
			pl++;
		}
		p++;
	}

	// Third entity: curves
	FilePrintString(geo_file, "// Curve");
	LineFeed(geo_file);
	p = polyline_of_surface_vector.begin(); // CC
	int LocalPLyIndex = 0;
	while (p != polyline_of_surface_vector.end()) // CC
	{
		pl_index++;
		lineIndex.push_back(pl_index);

		p_pline = *p;
		FilePrintString(geo_file, "Line Loop(");
		FilePrintInt(geo_file, pl_index);
		FilePrintString(geo_file, ")={");

		orient = polyline_of_surface_orient[LocalPLyIndex];

		std::vector<CGLLine*>::iterator pl = p_pline->line_vector.begin(); // WW NULL;
		if (orient > 0)
		{
			pl = p_pline->line_vector.begin();
			count = 0;
			CGLn0 = NULL;
			// WW CGLn1 = NULL;
			while (pl != p_pline->line_vector.end())
			{
				CGLn = *pl;

				// Determine the local orentation of this line
				ortL = 1; //
				CGLn = CGLn->CheckLineOutPut();
				if (count == 0)
				{
					CGPnt = p_pline->point_vector[0];
					if (CGLn->m_point2->PointDis(CGPnt) < DistTol)
						ortL = -1;
				}
				else
				{
					if (ortL0 > 0)
					{
						if (CGLn0->m_point2->PointDis(CGLn->m_point1) < DistTol)
							ortL = 1;
						else
							ortL = -1;
					}
					else
					{
						if (CGLn0->m_point1->PointDis(CGLn->m_point1) < DistTol)
							ortL = 1;
						else
							ortL = -1;
					}
				}
				ortL0 = ortL;
				CGLn0 = CGLn;

				FilePrintInt(geo_file, ortL * CGLn->mesh_index);
				++pl;
				if (pl != p_pline->line_vector.end())
					FilePrintString(geo_file, ", ");
				count++;
			}
		}
		else
		{
			pl = p_pline->line_vector.end(); // CC
			count = (int)p_pline->line_vector.size(); // CC
			Size0 = (int)p_pline->point_vector.size();
			CGLn0 = NULL;
			// WW      CGLn1 = NULL;
			while (pl != p_pline->line_vector.begin()) // CC
			{
				--pl;
				count--;
				// Determine the local orentation of this line
				ortL = 1; //
				CGLn = *pl;

				CGLn = CGLn->CheckLineOutPut();
				int p_pline_line_list_size = (int)p_pline->line_vector.size(); // CC
				if (count == p_pline_line_list_size - 1)
					CGPnt = p_pline->point_vector[Size0 - 1];
				// if(CGLn->m_point1->PointDis(CGPnt)<DistTol)
				// ortL = -1; //CC

				else
				{
					if (ortL0 > 0)
					{
						if (CGLn0->m_point2->PointDis(CGLn->m_point1) < DistTol)
							ortL = 1;
						// else ortL = -1;//CC
					}
					else if (CGLn0->m_point1->PointDis(CGLn->m_point1) < DistTol)
						ortL = 1;
					// else ortL = -1;//CC
				}
				ortL0 = ortL;
				CGLn0 = CGLn;

				FilePrintInt(geo_file, ortL * orient * CGLn->mesh_index);
				if (pl != p_pline->line_vector.begin())
					FilePrintString(geo_file, ", "); // CC
			}
		}

		FilePrintString(geo_file, "};");
		LineFeed(geo_file);

		LocalPLyIndex++;
		++p;
	}

	// Write surface
	FilePrintString(geo_file, "// Surface");
	LineFeed(geo_file);
	FilePrintString(geo_file, "Plane Surface(");
	FilePrintInt(geo_file, s_index);
	FilePrintString(geo_file, ")={");
	int lineIndex_size = (int)lineIndex.size();
	for (i = 0; i < lineIndex_size; i++)
	{
		FilePrintInt(geo_file, lineIndex[i]);
		int lineIndex_size = (int)lineIndex.size();
		if (i != lineIndex_size - 1)
			FilePrintString(geo_file, ", ");
	}
	FilePrintString(geo_file, "};");
	LineFeed(geo_file);
	FilePrintString(geo_file, "// Physical Surface");
	LineFeed(geo_file);
	FilePrintString(geo_file, "Physical Surface(");
	FilePrintInt(geo_file, s_index);
	FilePrintString(geo_file, ")={");
	FilePrintInt(geo_file, s_index);
	FilePrintString(geo_file, "};");
	LineFeed(geo_file);
}

// OKmsh
/**************************************************************************
   GeoLib-Method: GEOSurfaceGLI2GEO
   Task: convert PLI to GEO data
   Programing:
   01/2003 WW
   09/2003 OK for GeoLib
   11/2003 OK/MB mesh_surface
**************************************************************************/
void GEOSurfaceGLI2GEO(FILE* geo_file)
{
	// OK3908
	int mesh_surface(0);
	CGLPolyline* polyline;
	/* Write the forth entity, surface */
	std::vector<CGLPolyline*>::iterator p = polyline_vector.begin(); // CC
	while (p != polyline_vector.end())
	{
		polyline = *p;
		if (polyline->getType() == 2)
		{
			mesh_surface++;
			FilePrintString(geo_file, "// Surface");
			LineFeed(geo_file);
			FilePrintString(geo_file, "Plane Surface(");
			FilePrintInt(geo_file, mesh_surface);
			FilePrintString(geo_file, ")={");
			FilePrintInt(geo_file, mesh_surface);
			FilePrintString(geo_file, "};");
			LineFeed(geo_file);

			FilePrintString(geo_file, "// Physical Surface");
			LineFeed(geo_file);
			FilePrintString(geo_file, "Physical Surface(");
			FilePrintInt(geo_file, mesh_surface);
			FilePrintString(geo_file, ")={");
			FilePrintInt(geo_file, mesh_surface);
			FilePrintString(geo_file, "};");
			LineFeed(geo_file);
		}
		++p;
	}
}
/**************************************************************************
   GeoLib-Method: Surface::Write
   Task: Write surface data
   Programing:
   11/2003 OK/WW Implementation
   03/2004 OK for TINs
   07/2005 CC id
   09/2005 CC move from Geo_lib to geo_sfc
   11/05 CC Write function
**************************************************************************/
void Surface::Write(const std::string& path_name)
{
	const char* char_string;
	//-----------------------------------------------------------------------
	// GLI File
	FILE* gli_file = NULL;
	std::string gli_file_name;
	gli_file_name = path_name + ".gli";
	const char* gli_file_name_char = 0;
	gli_file_name_char = gli_file_name.data();
	gli_file = fopen(gli_file_name_char, "a");
	if (gli_file == NULL)
		return;
	else
	{
		fprintf(gli_file, "%s\n", "#SURFACE");
		//-----------------------------------------------------------------
		fprintf(gli_file, " %s\n", "$ID");
		fprintf(gli_file, "  %ld\n", id);
		fprintf(gli_file, " %s\n", "$NAME");
		char_string = name.data();
		fprintf(gli_file, "  %s\n", char_string);
		//-----------------------------------------------------------------
		// OK fprintf(gli_file," %s\n","$TYPE");
		// OK fprintf(gli_file,"  %i\n",m_surface->type);
		//-----------------------------------------------------------------
		fprintf(gli_file, " %s\n", "$MAT_GROUP");
		fprintf(gli_file, "  %i\n", mat_group);
		//-----------------------------------------------------------------
		int polyline_of_surface_list_length = (int)polyline_of_surface_vector.size(); // CC
		if (polyline_of_surface_list_length > 0)
		{
			fprintf(gli_file, " %s\n", "$POLYLINES");
			CGLPolyline* m_polyline = NULL;
			//   list<CGLPolyline*>::iterator p = m_surface->polyline_of_surface_list.begin();//CC
			std::vector<CGLPolyline*>::iterator p = polyline_of_surface_vector.begin();
			while (p != polyline_of_surface_vector.end()) // CC
			{
				m_polyline = *p;
				fprintf(gli_file, "  %s\n", m_polyline->name.c_str());
				++p;
			}
		}
		//-----------------------------------------------------------------
		if (TIN)
		{
			int surface_TIN_length = (int)TIN->Triangles.size();
			if (surface_TIN_length > 0)
			{
				fprintf(gli_file, " %s\n", "$TIN");
				fprintf(gli_file, "  %s.tin\n", TIN->name.c_str());
			}
		}
		//-----------------------------------------------------------------
	} // dat_in
	fclose(gli_file);
}
/**************************************************************************
   GeoLib-Method: Surface::ReArragePolylineList()
   Task: To arrange the polylines in a clock/anti-clock-wise way
   Programing:
   12/2003 WW Implementation
   01/2004 CC Modification
   10/2005 CC Modification
**************************************************************************/
void Surface::ReArrangePolylineList()
{
	bool Firstply = true;
	int top = 0;
	CGLPolyline* p_pline = NULL;
	CGLPolyline* pline_b = NULL;
	CGLLine* line_b_0 = NULL;
	CGLLine* line_b_1 = NULL;

	CGLLine* line_end = NULL;
	CGLLine* line_begin = NULL;
	std::vector<CGLPolyline*> buff;
	int ite = 0;
	int ite2 = 0;
	// list<CGLPolyline*>::const_iterator p = polyline_of_surface_list.begin();//CC
	std::vector<CGLPolyline*>::iterator p = polyline_of_surface_vector.begin();

	// The first polyline
	p_pline = *p;
	if (p_pline->data_type == 1) // OK not a toplogical polyline
		return;
	buff.push_back(p_pline);
	// polyline_of_surface_vector.remove(p_pline);//CC
	polyline_of_surface_vector.erase(polyline_of_surface_vector.begin() + ite);

	while (polyline_of_surface_vector.size())
	{
		p = polyline_of_surface_vector.begin();

		while (p != polyline_of_surface_vector.end())
		{
			top = (int)buff.size() - 1;
			p_pline = *p;
			// OK
			if (p_pline->line_vector.size() == 0) // CC
				return;

			line_begin = *p_pline->line_vector.begin(); // CC
			// line_end = *(--p_pline->line_vector.end()); //CC
			long sizea = (long)p_pline->line_vector.size();
			line_end = p_pline->line_vector[sizea - 1];

			pline_b = buff[top];
			line_b_0 = *pline_b->line_vector.begin();
			//  line_b_1 = *(--pline_b->line_vector.end());
			long sizeb = (long)pline_b->line_vector.size();
			line_b_1 = pline_b->line_vector[sizeb - 1];

			if (Firstply == true)
			{
				if ((abs(line_begin->point1) == abs(line_b_1->point2))
				    || (abs(line_end->point2) == abs(line_b_1->point2)))

				{
					buff.push_back(p_pline);
					// Strike out the compared polylines
					//    polyline_of_surface_vector.remove(p_pline);
					polyline_of_surface_vector.erase(polyline_of_surface_vector.begin() + ite2);

					ite2 = 0;
					Firstply = false;
					break;
				}
			} // end of if Firstply
			else if ((abs(line_begin->point1) == abs(line_b_1->point2))
			         || (abs(line_end->point2) == abs(line_b_1->point2))
			         || (abs(line_begin->point1) == abs(line_b_0->point1))
			         || (abs(line_end->point2) == abs(line_b_0->point1)))
			{
				buff.push_back(p_pline);
				// Strike out the compared polylines
				//              polyline_of_surface_vector.remove(p_pline);
				polyline_of_surface_vector.erase(polyline_of_surface_vector.begin() + ite2);
				ite2 = 0;
				break;
			}

			++ite2;
			++p;
		} // end of while
	} // end of while
	//

	////////////////////////////////////////////////////////////////////////////
	long buff_size = (long)buff.size();
	for (int i = 0; i < buff_size; i++)
	{
		pline_b = buff[i];
		polyline_of_surface_vector.push_back(pline_b);
	}

	// Release the buff
	while (buff.size())
		buff.pop_back();
}

/**************************************************************************
   GeoLib-Method: Surface::PolylineOrientation(void)
   Task: check the polyline orientation of surface
   Programing:
   12/2003 CC Implementation
   01/2004 CC Modification
**************************************************************************/
void Surface::PolylineOrientation()
{
	int orient = 1;
	bool firstPline = true;
	CGLPolyline* p_pline = NULL;
	CGLLine* m_line = NULL;

	CGLLine* line_end = NULL;
	// WW CGLLine *line_begin=NULL;

	// list<CGLPolyline*>::const_iterator p = polyline_of_surface_list.begin();
	std::vector<CGLPolyline*>::iterator p = polyline_of_surface_vector.begin();

	p_pline = *p;
	if (p_pline->data_type == 1) // OK not a toplogical polyline
		return;

	while (p != polyline_of_surface_vector.end())
	{
		p_pline = *p;

		std::vector<CGLLine*>::iterator pl = p_pline->line_vector.begin(); // CC

		// Determine the orientation
		if (!firstPline)
		{
			m_line = *pl;
			if (m_line->point1 == line_end->point2 || m_line->point1 == line_end->point1

			    )
				orient = 1;

			else
				orient = -1;
		}

		polyline_of_surface_orient.push_back(orient);
		// if it is the first polyline
		// WW   line_begin = *pl;
		if (orient == 1)
			while (pl != p_pline->line_vector.end()) // CC
			{
				m_line = *pl;

				++pl;
			}
		else
		{
			std::vector<CGLLine*>::iterator pl2 = p_pline->line_vector.end(); // CC
			--pl2;
			while (pl2 != p_pline->line_vector.begin()) // CC
			{
				m_line = *pl2;

				--pl2;
			}
			m_line = *pl2;
		}

		line_end = m_line;
		firstPline = false;

		++p;
	} // end of while
}
/**************************************************************************
   GeoLib-Method: GEORemoveAllSurfaces()
   Task:
   Programing:
   12/2003 CC Implementation
   01/2005 OK List destructor
   CCToDo Surface destructor
   08/2005 CC
   03/2006 CC
   05/2006 TK
**************************************************************************/
void GEORemoveAllSurfaces()
{
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		delete surface_vector[i];
		surface_vector[i] = NULL;
	}
	surface_vector.clear();
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   04/2005 CC   Implementation
   07/2005 CC Modification remove polyline
   03/2060 CC destructor
**************************************************************************/
void GEORemoveSurface(long nSel)
{
	Surface* m_sfc = NULL;
	m_sfc = surface_vector[nSel];
	delete m_sfc;
	surface_vector.erase(surface_vector.begin() + nSel);
}

/**************************************************************************
   GeoLib-Method:  GEOSurfaceTopology(long)
   Task:
   Programing:
   12/2003 CC Implementation
   03/2004 OK only if polyline_of_surface_list exists
   08/2004 CC Modification
   09/2004 OK for what surfaces (data_type)
   01/2006 OK Compact
   10/2010 TF some improvements in code structure
**************************************************************************/
void GEOSurfaceTopology(void)
{
	for (size_t i = 0; i < surface_vector.size(); i++)
	{
		Surface* sfc = surface_vector[i];
		if (sfc->polyline_of_surface_vector.size() == 0)
			continue;
		if (sfc->order == false)
		{
			for (size_t j = 0; j < sfc->polyline_of_surface_vector.size(); j++)
			{
				CGLPolyline* ply = sfc->polyline_of_surface_vector[j];
				if (ply->getDataType() == 1) // OK not a toplogical polyline
					continue;
				// compute lines and put into polyline line, here all orientation is 1
				ply->ComputeLines();
			}
			//..................................................................
			// Reorder the polyline list
			sfc->ReArrangePolylineList();
			//..................................................................
			// Check start and end point of polyline, change the orientation value
			sfc->PolylineOrientation();
			// surface polyline --point rearrangement
			sfc->order = false;
		}
		//----------------------------------------------------------------------
	}
}

/**************************************************************************
   GeoLib-Method: GetPointsIn
   Task:
   Programing:
   01/2004 OK Implementation
   08/2005 CC Modification
**************************************************************************/
/*long* Surface::GetPointsIn(long* number_of_nodes)
   {
   long *nodes = NULL;
   long i;
   double *xp=NULL,*yp=NULL,*zp=NULL;
   long anz_relevant = 0;
   CGLPoint* m_pnt = NULL;

   if(!polygon_point_vector.empty()) {
    xp = (double*) Malloc(((long)polygon_point_vector.size())*sizeof(double));
    yp = (double*) Malloc(((long)polygon_point_vector.size())*sizeof(double));
    zp = (double*) Malloc(((long)polygon_point_vector.size())*sizeof(double));
    long polygon_point_vector_length = (long)polygon_point_vector.size();
    for(i=0;i<polygon_point_vector_length;i++) {
      xp[i] = polygon_point_vector[i]->x;
      yp[i] = polygon_point_vector[i]->y;
      zp[i] = polygon_point_vector[i]->z;
    }

    for(i=0;i<NodeListSize();i++) {
      if (GetNode(i)==NULL) continue;
      m_pnt->x = GetNodeX(i);
      m_pnt->y = GetNodeY(i);
      m_pnt->z = GetNodeZ(i);
      if(m_pnt->IsInsidePolygonPlain(
                                   xp,yp,zp,\
                                   (long)polygon_point_vector.size())) {
        anz_relevant++;

        nodes = (long *) Realloc(nodes,sizeof(long)*anz_relevant);

        nodes[anz_relevant-1] = i;
      }
    }

   }
   // Destructions
   // nodes extern
   xp = (double*) Free(xp);
   yp = (double*) Free(yp);
   zp = (double*) Free(zp);
   //
   *number_of_nodes = anz_relevant;
   return nodes;
   }*/
/**************************************************************************
   GeoLib-Method: CalcCenterPoint
   Task:
   Programing:
   01/2004 OK Implementation
   last modification:
**************************************************************************/
void Surface::CalcCenterPoint(void)
{
	long i;
	long polygon_point_vector_length = (long)polygon_point_vector.size();
	center_point[0] = 0.0;
	center_point[1] = 0.0;
	center_point[2] = 0.0;
	for (i = 0; i < polygon_point_vector_length; i++)
	{
		center_point[0] += polygon_point_vector[i]->x;
		center_point[1] += polygon_point_vector[i]->y;
		center_point[2] += polygon_point_vector[i]->z;
	}
	center_point[0] /= (double)polygon_point_vector_length;
	center_point[1] /= (double)polygon_point_vector_length;
	center_point[2] /= (double)polygon_point_vector_length;
}
/**************************************************************************
   GeoSys-GUI Function
   Programing:
   01/2004 OK Implementation
   09/2005 CC Modification
**************************************************************************/
/*bool Surface::IsPointInSurface(CGLPoint *m_point)
   {
   bool ok = false;
   m_point = m_point;
   return ok;
   }*/
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   01/2005 OK TIN based on layer polylines
   last modification:
**************************************************************************/
void Surface::CreateTIN(void)
{
	long i;
	int j;
	CGLPoint m_point;
	CTriangle* m_triangle;
	//----------------------------------------------------------------------
	TIN = new CTIN();
	TIN->name = name;
	//----------------------------------------------------------------------
	// OK41
	if (type == 1) // TIN based on 2 layer polylines
	{ //  list<CGLPolyline*>::const_iterator  p = polyline_of_surface_list.begin();
		std::vector<CGLPolyline*>::iterator p = polyline_of_surface_vector.begin();
		CGLPolyline* m_polyline1 = NULL;
		CGLPolyline* m_polyline2 = NULL;
		while (p != polyline_of_surface_vector.end())
		{
			m_polyline1 = *p;
			++p;
			m_polyline2 = *p;
			break;
		}
		long no_points = (long)m_polyline1->point_vector.size();
		double xp[3], yp[3], zp[3];
		for (i = 0; i < no_points - 1; i++)
		{
			// first triangle of quad
			xp[0] = m_polyline1->point_vector[i]->x;
			yp[0] = m_polyline1->point_vector[i]->y;
			zp[0] = m_polyline1->point_vector[i]->z;
			xp[1] = m_polyline1->point_vector[i + 1]->x;
			yp[1] = m_polyline1->point_vector[i + 1]->y;
			zp[1] = m_polyline1->point_vector[i + 1]->z;
			xp[2] = m_polyline2->point_vector[i]->x;
			yp[2] = m_polyline2->point_vector[i]->y;
			zp[2] = m_polyline2->point_vector[i]->z;
			m_triangle = new CTriangle;
			m_triangle->number = (long)TIN->Triangles.size();
			for (j = 0; j < 3; j++)
			{
				m_triangle->x[j] = xp[j];
				m_triangle->y[j] = yp[j];
				m_triangle->z[j] = zp[j];
			}
			TIN->Triangles.push_back(m_triangle);
			// second triangle of quad
			xp[0] = m_polyline2->point_vector[i]->x;
			yp[0] = m_polyline2->point_vector[i]->y;
			zp[0] = m_polyline2->point_vector[i]->z;
			xp[1] = m_polyline2->point_vector[i + 1]->x;
			yp[1] = m_polyline2->point_vector[i + 1]->y;
			zp[1] = m_polyline2->point_vector[i + 1]->z;
			xp[2] = m_polyline1->point_vector[i + 1]->x;
			yp[2] = m_polyline1->point_vector[i + 1]->y;
			zp[2] = m_polyline1->point_vector[i + 1]->z;
			m_triangle = new CTriangle;
			m_triangle->number = (long)TIN->Triangles.size();
			for (j = 0; j < 3; j++)
			{
				m_triangle->x[j] = xp[j];
				m_triangle->y[j] = yp[j];
				m_triangle->z[j] = zp[j];
			}
			TIN->Triangles.push_back(m_triangle);
		} // no_points
	}

	//----------------------------------------------------------------------
	// Memory
	if (TIN->Triangles.size() == 0)
		delete TIN;
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   last modification:08/2005 CC
**************************************************************************/
void GEOWriteSurfaceTINs(string file_path)
{
	vector<Surface*>::iterator ps = surface_vector.begin();
	Surface* m_surface = NULL;
	while (ps != surface_vector.end())
	{
		m_surface = *ps;
		if (!m_surface->TIN)
		{
			++ps;
			continue;
		}
		else
		{
			m_surface->WriteTIN(file_path);
			++ps;
		}
	}
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   last modification:
**************************************************************************/
void GEOWriteSurfaceTINsTecplot(string file_path)
{
	vector<Surface*>::iterator ps = surface_vector.begin();
	Surface* m_surface = NULL;
	while (ps != surface_vector.end())
	{
		m_surface = *ps;
		if (!m_surface->TIN)
		{
			++ps;
			continue;
		}
		else
		{
			m_surface->WriteTINTecplot(file_path);
			++ps;
		}
	}
}

/*
   void Surface::operator = (const Surface m_surface)
   {
   name = m_surface.name;
   }
 */

/**************************************************************************
   GEOLib-Method:
   Task: Surface read function
   Programing:
   03/2004 OK Implementation
   05/2004 CC Modification
**************************************************************************/
void GEOReadSurfaces(const std::string& file_name_path_base)
{
	Surface* m_surface = NULL;
	char line[MAX_ZEILEN];
	string sub_line;
	string line_string;
	string gli_file_name;
	ios::pos_type position;
	//========================================================================
	// File handling
	gli_file_name = file_name_path_base + ".gli";
	ifstream gli_file(gli_file_name.data(), ios::in);
	if (!gli_file.good())
		return;
	gli_file.seekg(0L, ios::beg); // rewind?
	//========================================================================
	// Keyword loop
	while (!gli_file.eof())
	{
		gli_file.getline(line, MAX_ZEILEN);
		line_string = line;
		if (line_string.find("#STOP") != string::npos) // 11.08.2011. WW
			break;

		//----------------------------------------------------------------------
		if (line_string.find("#SURFACE") != string::npos) // keyword found
		{
			m_surface = new Surface();
			m_surface->AssignColor();
			// OK not here m_surface->TIN = new CTIN;
			position = m_surface->Read(&gli_file);
			surface_vector.push_back(m_surface);
			gli_file.seekg(position, ios::beg);
		} // keyword found
	} // eof
	gli_file.close();
}

/**************************************************************************
   GeoLib-Method: Read
   Task: Read volume data from file
   Programing:
   03/2004 OK Implementation
   05/2005 OK EPSILON
   09/2005 CC file_path_base
**************************************************************************/
std::ios::pos_type Surface::Read(std::ifstream* gli_file)
{
	char line[MAX_ZEILEN];
	std::string line_string;

	bool new_keyword = false;
	std::ios::pos_type position;
	CGLPolyline* m_polyline = NULL;
	CGLPoint* m_point = NULL;
	bool ok_true = true;
	type_name = "POLYLINE"; // CC8888
	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		position = gli_file->tellg();
		gli_file->getline(line, MAX_ZEILEN);
		line_string = line;
		if (line_string.find("#") != string::npos)
		{
			new_keyword = true;
			break;
		}
		//.....................................................................
		if (line_string.find("$ID") != string::npos) // subkeyword found CC
		{
			gli_file->getline(line, MAX_ZEILEN);
			line_string = line;
			remove_white_space(&line_string);
			id = strtol(line_string.data(), NULL, 0);
			continue;
		}
		//....................................................................
		if (line_string.find("$NAME") != string::npos) // subkeyword found
		{
			gli_file->getline(line, MAX_ZEILEN);
			line_string = line;
			remove_white_space(&line_string);
			name = line_string.substr(0);
			continue;
		} // subkeyword found
		//....................................................................
		if (line_string.find("$EPSILON") != string::npos) // subkeyword found
		{
			gli_file->getline(line, MAX_ZEILEN);
			line_string = line;
			remove_white_space(&line_string);
			epsilon = strtod(line_string.data(), NULL);
			continue;
		} // subkeyword found
		//....................................................................
		if (line_string.find("$TYPE") != string::npos) // subkeyword found
		{
			gli_file->getline(line, MAX_ZEILEN);
			line_string = line;
			remove_white_space(&line_string);
			type = strtol(line_string.data(), NULL, 0);
			if (type == 100)
			{
				int p_index;
				polygon_point_vector.clear();
				*gli_file >> p_index; // Center of top surface of cylinder
				m_point = GEOGetPointById(p_index); // CC
				polygon_point_vector.push_back(m_point);
				*gli_file >> p_index; // Center of bottom surface of cylinder
				m_point = GEOGetPointById(p_index); // CC
				polygon_point_vector.push_back(m_point);
				*gli_file >> Radius >> epsilon;
				gli_file->ignore(MAX_ZEILEN, '\n');
			}
		} // subkeyword found
		//....................................................................
		if (line_string.find("$TIN") != string::npos) // subkeyword found
		{
			(*gli_file) >> line_string;
			remove_white_space(&line_string);
			TIN = new CTIN();
			// TIN->name = name;
			ReadTIN(line_string);
			std::string cut_string;
			TIN->name = get_sub_string2(line_string, ".", &cut_string);
			type = 1; // OK41
			type_name = "TIN"; // CC8888
			data_name = line_string; // CC8888
		} // subkeyword found
		//....................................................................
		if (line_string.find("$MAT_GROUP") != string::npos) // subkeyword found
		{
			gli_file->getline(line, MAX_ZEILEN);
			line_string = line;
			mat_group = strtol(line_string.data(), NULL, 0);
		} // subkeyword found
		//....................................................................
		if (line_string.find("$POLYLINES") != string::npos) // subkeyword found
			// int pos1 = 0, pos2 = 0;
			while (ok_true)
			{
				position = gli_file->tellg();
				gli_file->getline(line, MAX_ZEILEN);
				line_string = line;
				type = 0;
				remove_white_space(&line_string);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					gli_file->seekg(position, ios::beg);
					break;
				}
				m_polyline = GEOGetPLYByName(line_string); // CC
				if (m_polyline)
					polyline_of_surface_vector.push_back(m_polyline);
				// type = m_polyline->data_type; //MB Please set surface type with surface type, not polyline type
				else
					cout << "Warning in Surface::Read, polyline not found: " << line_string << "\n";
			}
	}
	return position;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   09/2005 CC file_path_base
   10/2005 OK Path 2nd
**************************************************************************/
void Surface::ReadTIN(const std::string& tin_file_name)
{
	string cut_string;
	string delimiter_type(" ");
	char line[MAX_ZEILEN];
	string line_string;
	string sub_string;
	string buffer;
	CTriangle* m_triangle = NULL;
	//----------------------------------------------------------------------
	// File handling
	string tin_file_name_path;
	/*  //11.08.2011 WW
	   //  tin_file_name_path = FileName;  //WW/JOD // LB: commented out to compile GEO without dependency on FEM
	   basic_string <char>::size_type indexCh1a;
	   indexCh1a = tin_file_name_path.find_last_of("\\");

	   if( indexCh1a < tin_file_name_path.size()) { // JOD 4.7.10,  \ exists in path, DEBUG case
	   string stra = tin_file_name_path.substr (0, indexCh1a);
	   tin_file_name_path.clear();
	   tin_file_name_path = stra+"\\"+tin_file_name;

	   }
	   else     //   \ does not exist, RELEASE case
	   tin_file_name_path = tin_file_name;
	 */

	tin_file_name_path = FilePath + tin_file_name; // 11.08.2011. WW
	ifstream tin_file(tin_file_name_path.data(), ios::in);
	if (!tin_file.good())
		return;
	tin_file.seekg(0L, ios::beg);
	//========================================================================
	while (!tin_file.eof())
	{
		tin_file.getline(line, MAX_ZEILEN);
		buffer = line;
		if (buffer.size() == 0)
			continue;
		//
		m_triangle = new CTriangle;
		// triangle number
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->number = strtol(sub_string.data(), NULL, 0);
		// point 1
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->x[0] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->y[0] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->z[0] = strtod(sub_string.data(), NULL);
		// point 2
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->x[1] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->y[1] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->z[1] = strtod(sub_string.data(), NULL);
		// point 3
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->x[2] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->y[2] = strtod(sub_string.data(), NULL);
		//
		buffer = cut_string;
		sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
		m_triangle->z[2] = strtod(sub_string.data(), NULL);
		//
		TIN->Triangles.push_back(m_triangle);
	}
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   01/2005 OK File handling
   08/2005 CC string file_path
**************************************************************************/
void Surface::WriteTIN(const std::string& file_path)
{
	string delimiter(" ");
	long i;
	//----------------------------------------------------------------------
	if (TIN && ((long)TIN->Triangles.size() == 0))
		return;
	//----------------------------------------------------------------------
	// File handling
	string tin_path;
	tin_path = file_path;
	string tin_file_name = TIN->name + TIN_FILE_EXTENSION;
	string tin_path_base_type = tin_path + tin_file_name;
	fstream tin_file(tin_path_base_type.data(), ios::trunc | ios::out);
	tin_file.setf(ios::scientific, ios::floatfield);
	tin_file.precision(12);
	//----------------------------------------------------------------------
	if (!tin_file.good())
		return;
	tin_file.seekg(0L, ios::beg);
	for (i = 0; i < (long)TIN->Triangles.size(); i++)
		tin_file << TIN->Triangles[i]->number << delimiter << TIN->Triangles[i]->x[0] << delimiter
		         << TIN->Triangles[i]->y[0] << delimiter << TIN->Triangles[i]->z[0] << delimiter
		         << TIN->Triangles[i]->x[1] << delimiter << TIN->Triangles[i]->y[1] << delimiter
		         << TIN->Triangles[i]->z[1] << delimiter << TIN->Triangles[i]->x[2] << delimiter
		         << TIN->Triangles[i]->y[2] << delimiter << TIN->Triangles[i]->z[2] << "\n";
	// tin_file.close();
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   last modification:08/2005 CC
**************************************************************************/
void GEOCreateLayerSurfaceTINs(int nb_prism_layers, double thickness_prism_layers)
{
	long i, j;
	long no_TIN_Triangles;
	vector<Surface*>::iterator ps = surface_vector.begin();
	Surface* m_surface = NULL;
	Surface* m_surface_new = NULL;
	char layer_number[10];
	CTriangle* m_triangle = NULL;
	vector<Surface*> new_surfaces;

	// Go through the existing surfaces (GEO_type)
	while (ps != surface_vector.end())
	{
		m_surface = *ps;
		if (!m_surface->TIN)
		{
			++ps;
			continue;
		}
		no_TIN_Triangles = (long)m_surface->TIN->Triangles.size();
		if (m_surface->type == 0)
			// Go through the layers
			for (i = 0; i < nb_prism_layers; i++)
			{
				m_surface_new = new Surface;
				sprintf(layer_number, "%ld", i + 1);
				m_surface_new->name = m_surface->name + "_L" + layer_number;
				m_surface_new->type = m_surface->type;
				m_surface_new->TIN = new CTIN;
				m_surface_new->TIN->name = m_surface_new->name;
				for (j = 0; j < no_TIN_Triangles; j++)
				{
					m_triangle = new CTriangle;
					m_triangle->number = m_surface->TIN->Triangles[j]->number;
					m_triangle->x[0] = m_surface->TIN->Triangles[j]->x[0];
					m_triangle->x[1] = m_surface->TIN->Triangles[j]->x[1];
					m_triangle->x[2] = m_surface->TIN->Triangles[j]->x[2];
					m_triangle->y[0] = m_surface->TIN->Triangles[j]->y[0];
					m_triangle->y[1] = m_surface->TIN->Triangles[j]->y[1];
					m_triangle->y[2] = m_surface->TIN->Triangles[j]->y[2];
					m_triangle->z[0] = m_surface->TIN->Triangles[j]->z[0] - (i + 1) * thickness_prism_layers;
					m_triangle->z[1] = m_surface->TIN->Triangles[j]->z[1] - (i + 1) * thickness_prism_layers;
					m_triangle->z[2] = m_surface->TIN->Triangles[j]->z[2] - (i + 1) * thickness_prism_layers;
					m_surface_new->TIN->Triangles.push_back(m_triangle);
				}
				new_surfaces.push_back(m_surface_new);
			}
		++ps;
	}
	//
	int new_surfaces_vector_length = (int)new_surfaces.size();
	for (i = 0; i < new_surfaces_vector_length; i++)
	{
		m_surface_new = new_surfaces[i];
		surface_vector.push_back(m_surface_new);
	}
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   01/2005 OK File handling
   last modification:
   09/2005 CC file_path
**************************************************************************/
void Surface::WriteTINTecplot(const std::string& file_path)
{
	string delimiter(", ");
	long i;
	if (TIN && (TIN->Triangles.size() > 0))
	{
		//--------------------------------------------------------------------
		// File handling
		string tin_path;
		//  CGSProject* m_gsp = GSPGetMember("gli");
		// if(m_gsp)
		tin_path = file_path; // CC
		string tin_file_name = TIN->name + TEC_FILE_EXTENSIONS;
		string tin_path_base_type = tin_path + tin_file_name;
		fstream tin_file(tin_path_base_type.data(), ios::trunc | ios::out);
		tin_file.setf(ios::scientific, ios::floatfield);
		tin_file.precision(12);
		//--------------------------------------------------------------------
		if (!tin_file.good())
			return;
		tin_file.seekg(0L, ios::beg);
		tin_file << "VARIABLES = X,Y,Z"
		         << "\n";
		long no_triangles = (long)TIN->Triangles.size();
		long no_nodes = 3 * no_triangles;
		tin_file << "ZONE T = " << TIN->name << delimiter << "N = " << no_nodes << delimiter << "E = " << no_triangles
		         << delimiter << "F = FEPOINT" << delimiter << "ET = TRIANGLE"
		         << "\n";
		for (i = 0; i < no_triangles; i++)
		{
			tin_file << TIN->Triangles[i]->x[0] << " " << TIN->Triangles[i]->y[0] << " " << TIN->Triangles[i]->z[0]
			         << "\n";
			tin_file << TIN->Triangles[i]->x[1] << " " << TIN->Triangles[i]->y[1] << " " << TIN->Triangles[i]->z[1]
			         << "\n";
			tin_file << TIN->Triangles[i]->x[2] << " " << TIN->Triangles[i]->y[2] << " " << TIN->Triangles[i]->z[2]
			         << "\n";
		}
		for (i = 0; i < no_triangles; i++)
			tin_file << 3 * i + 1 << " " << 3 * i + 2 << " " << 3 * i + 3 << "\n";
	}
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   04/2004 OK Implementation
   01/2005 OK TIN case
   last modification:
**************************************************************************/
/*vector<long> Surface::GetMSHNodesClose()
   {
   long i,j;
   long no_nodes = 0;
   long *nodes_array = NULL;
   long no_points;
   CGLPoint m_node;
   CGLPolyline* m_polyline = NULL;
   CGLPolyline* m_polyline1 = NULL;
   CGLPolyline* m_polyline2 = NULL;
   vector<long>nodes_vector;
   double xp[3],yp[3],zp[3];
   CTriangle *m_triangle = NULL;
   long no_triangles;

   list<CGLPolyline*>::const_iterator p = polyline_of_surface_list.begin();
   nodes_vector.clear();

   switch(data_type) {
    //--------------------------------------------------------------------
    case 0: // surface polygon
      nodes_array = GetPointsIn(this,&no_nodes);//change this into m_sfc todo
      for(i=0;i<no_nodes;i++)
        nodes_vector.push_back(nodes_array[i]);
      break;
    //--------------------------------------------------------------------
    case 1: // TIN
   //OK41
      no_triangles = (long)TIN->Triangles.size();
      for(i=0;i<no_triangles;i++){
        m_triangle = TIN->Triangles[i];
        xp[0] = m_triangle->x[0];
        yp[0] = m_triangle->y[0];
        zp[0] = m_triangle->z[0];
        xp[1] = m_triangle->x[1];
        yp[1] = m_triangle->y[1];
        zp[1] = m_triangle->z[1];
        xp[2] = m_triangle->x[2];
        yp[2] = m_triangle->y[2];
        zp[2] = m_triangle->z[2];
        for(j=0;j<NodeListSize();j++){
          m_node.x = GetNodeX(j);
          m_node.y = GetNodeY(j);
          m_node.z = GetNodeZ(j);
          if(m_node.IsInsideTriangle(xp,yp,zp)){
            nodes_vector.push_back(j);
          }
        }
      }
      break;
    //--------------------------------------------------------------------
    case 2: // 2 vertical polylines //OK
      // .................................................................
      // nodes close to first polyline
      p = polyline_of_surface_list.begin();
      while(p!=polyline_of_surface_list.end()) {
        m_polyline = *p;
        nodes_array = MSHGetNodesClose(&no_nodes,m_polyline);//CC
        break;
      }
      // .....................................................................
      // using triangles
      p = polyline_of_surface_list.begin();
      while(p!=polyline_of_surface_list.end()) {
        m_polyline1 = *p;
   ++p;
        m_polyline2 = *p;
        break;
      }
      no_points = (long)m_polyline1->point_vector.size();

      for(j=0;j<no_nodes;j++) {
        m_node.x = GetNodeX(nodes_array[j]);
        m_node.y = GetNodeY(nodes_array[j]);
        m_node.z = GetNodeZ(nodes_array[j]);
        for(i=0;i<no_points-1;i++) {
          // first triangle of quad
          xp[0] = m_polyline1->point_vector[i]->x;
          yp[0] = m_polyline1->point_vector[i]->y;
          zp[0] = m_polyline1->point_vector[i]->z;
          xp[1] = m_polyline1->point_vector[i+1]->x;
          yp[1] = m_polyline1->point_vector[i+1]->y;
          zp[1] = m_polyline1->point_vector[i+1]->z;
          xp[2] = m_polyline2->point_vector[i]->x;
          yp[2] = m_polyline2->point_vector[i]->y;
          zp[2] = m_polyline2->point_vector[i]->z;
          if(m_node.IsInsideTriangle(xp,yp,zp)) {
            nodes_vector.push_back(nodes_array[j]);
          }
          // second triangle of quad
          xp[0] = m_polyline2->point_vector[i]->x;
          yp[0] = m_polyline2->point_vector[i]->y;
          zp[0] = m_polyline2->point_vector[i]->z;
          xp[1] = m_polyline2->point_vector[i+1]->x;
          yp[1] = m_polyline2->point_vector[i+1]->y;
          zp[1] = m_polyline2->point_vector[i+1]->z;
          xp[2] = m_polyline1->point_vector[i+1]->x;
          yp[2] = m_polyline1->point_vector[i+1]->y;
          zp[2] = m_polyline1->point_vector[i+1]->z;
          if(m_node.IsInsideTriangle(xp,yp,zp)) {
            nodes_vector.push_back(nodes_array[j]);
          }
        } // no_points
      } // no_nodes

      break;
   }
   return nodes_vector;
   }*/

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   05/2004 CC Implementation Assign randam color to surface
   last modification:
**************************************************************************/
void Surface::AssignColor()
{
	long r_number = rand();
	m_color[0] = (long)(r_number / (double)RAND_MAX * 256);
	r_number = rand();
	m_color[1] = (long)(r_number / (double)RAND_MAX * 256);
	r_number = rand();
	m_color[2] = (long)(r_number / (double)RAND_MAX * 256);
}

/**************************************************************************
   GeoLib-Method: GetSurfaceVector
   Task:
   Programing:
   05/2004 TK Implementation
   08/2005 CC vector
**************************************************************************/
vector<Surface*> GetSurfaceVector(void)
{
	return surface_vector;
}
/**************************************************************************
   FEMLib-Method:
   Task: write function based on CBoundaryCondition::WriteTecplot
   Programing:
   01/2005 OK Implementation
   last modification:
**************************************************************************/
void Surface::WriteTecplot(fstream* tec_file)
{
	long i;
	CGLPolyline* m_polyline1 = NULL;
	CGLPolyline* m_polyline2 = NULL;
	// list<CGLPolyline*>::const_iterator p;
	vector<CGLPolyline*>::iterator p;
	long no_points = 0;
	vector<CTriangle*> triangle_vector;
	//----------------------------------------------------------------------
	// Write header
	*tec_file << "VARIABLES = X,Y,Z"
	          << "\n";

	switch (type)
	{
		case 22:
			p = polyline_of_surface_vector.begin();
			while (p != polyline_of_surface_vector.end())
			{
				m_polyline1 = *p;
				++p;
				m_polyline2 = *p;
				break;
			}
			no_points = (long)m_polyline1->point_vector.size();
			break;
	}
	long no_nodes = 2 * no_points;
	// long no_elements = triangle_vector.size();
	long no_elements = 2 * (no_points - 1);
	// Write
	*tec_file << "ZONE T = " << name << ", "
	          << "N = " << no_nodes << ", "
	          << "E = " << no_elements << ", "
	          << "F = FEPOINT"
	          << ", "
	          << "ET = TRIANGLE"
	          << "\n";
	//----------------------------------------------------------------------
	// Write data
	if (m_polyline1)
		for (i = 0; i < no_points; i++)
			*tec_file << m_polyline1->point_vector[i]->x << " " << m_polyline1->point_vector[i]->y << " "
			          << m_polyline1->point_vector[i]->z << " "
			          << "\n";

	if (m_polyline2)
		for (i = 0; i < no_points; i++)
			*tec_file << m_polyline2->point_vector[i]->x << " " << m_polyline2->point_vector[i]->y << " "
			          << m_polyline2->point_vector[i]->z << " "
			          << "\n";

	for (i = 0; i < no_points - 1; i++)
		*tec_file << i + 1 << " " << i + 1 + 1 << " " << no_points + i + 1 << "\n";
	for (i = 0; i < no_points - 1; i++)
		*tec_file << no_points + i + 1 << " " << no_points + i + 1 + 1 << " " << i + 1 + 1 << "\n";
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/
Surface* GEOGetSFCByName(const std::string& name)
{
	Surface* m_sfc = NULL;
	vector<Surface*>::iterator p = surface_vector.begin(); // CC
	while (p != surface_vector.end())
	{
		m_sfc = *p;
		if (m_sfc->name.compare(name) == 0)
			return m_sfc;
		++p;
	}
	return NULL;
}

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   08/2005 OK Implementation
   08/2005 CC Modification
   ToDo:
   polygon_point_vector.resize()
**************************************************************************/
void Surface::PolygonPointVector()
{
	int LocalPlyIndex = 0;
	if (type == 100) // Cylinder surface WW
		return;
	polygon_point_vector.clear();
	long i = 0;
	CGLLine* gl_line = NULL;
	CGLPoint* gl_point = NULL;
	CGLPolyline* gl_polyline = NULL;
	// CC 02/2008-----------------------------------------------begin
	vector<CGLPolyline*>::iterator p1 = polyline_of_surface_vector.begin();
	if (polyline_of_surface_vector.size() == 1)
	{
		// Closed polyline
		gl_polyline = *p1;
		vector<CGLPoint*>::iterator pt = gl_polyline->point_vector.begin();
		while (pt != gl_polyline->point_vector.end())
		{
			gl_point = *pt;

			polygon_point_vector.push_back(gl_point);
			++pt;
		}
		vector<CGLPoint*>::iterator pn = gl_polyline->point_vector.begin();
		gl_point = *pn;
		polygon_point_vector.push_back(gl_point);
	}
	else
	{
		// CC-------------------------------------------------------end
		vector<CGLPolyline*>::iterator p1 = polyline_of_surface_vector.begin();
		//----------------------------------------------------------------------
		while (p1 != polyline_of_surface_vector.end())
		{
			gl_polyline = *p1;
			//....................................................................
			if (gl_polyline->line_vector.size() > 0)
			{
				vector<CGLLine*>::iterator pb = gl_polyline->line_vector.begin();
				vector<CGLLine*>::iterator pe = gl_polyline->line_vector.end();
				if (gl_polyline->line_vector.size() == 0)
					return;
				if (polyline_of_surface_orient[LocalPlyIndex] == 1)
					while (pb != gl_polyline->line_vector.end())
					{
						gl_line = *pb;
						polygon_point_vector.push_back(GEOGetPointById(gl_line->point2));

						++i;
						++pb;
					}
				else
					while (pe != gl_polyline->line_vector.begin())
					{
						gl_line = *(--pe);
						polygon_point_vector.push_back(GEOGetPointById(gl_line->point1));
						++i;
					}
			}
			//....................................................................
			else // point_vector exists
			{
				vector<CGLPoint*>::iterator pv = gl_polyline->point_vector.begin();
				while (pv != gl_polyline->point_vector.end())
				{
					gl_point = *pv;
					polygon_point_vector.push_back(gl_point);
					++pv;
					++i;
				}
			}
			++p1;
			LocalPlyIndex++;
		}
		// CC--------------------------------------------------------begin
	}
	// CC--------------------------------------------------------end
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   09/2005 CC Modification
**************************************************************************/
void GEOCreateSurfacePointVector(void)
{
	vector<Surface*>::iterator ps = surface_vector.begin();
	Surface* m_surface = NULL;
	while (ps != surface_vector.end())
	{
		m_surface = *ps;
		m_surface->PolygonPointVector(); // CC
		++ps;
	}
}

/**************************************************************************
   GeoSys-GUI Function
   Programing:
   01/2004 OK Implementation
**************************************************************************/
// REMOVE CANDIDATE
// bool Surface::PointInSurface(CGLPoint *m_point)
//{
//  bool ok = false;
//  m_point = m_point;
//  return ok;
//}

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
void GEOUnselectSFC()
{
	Surface* m_sfc = NULL;
	vector<Surface*>::const_iterator p_sfc;
	p_sfc = surface_vector.begin();
	while (p_sfc != surface_vector.end())
	{
		m_sfc = *p_sfc;
		m_sfc->highlighted = false;
		++p_sfc;
	}
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   11/2005 CC Modification
**************************************************************************/
void GEOWriteSurfaces(const std::string& path_name)
{
	Surface* m_sfc = NULL;
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		m_sfc = surface_vector[i];
		m_sfc->Write(path_name);
	}
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
void GEORemoveSFC(Surface* m_sfc)
{
	Surface* m_sfc_this = NULL;
	// WW vector<Surface*>::const_iterator p_sfc = surface_vector.begin();
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		m_sfc_this = surface_vector[i];
		if (m_sfc_this->name.compare(m_sfc->name) == 0)
		{
			delete m_sfc_this;
			surface_vector.erase(surface_vector.begin() + i);
			// i--;
			return;
		}
	}
}
/* List example
   vector<Surface*>::const_iterator p_sfc = surface_vector.begin();
   while(p_sfc!=surface_vector.end()){
    p_sfc = surface_vector.begin();
    while(p_sfc!=surface_vector.end()) {
      m_sfc = *p_sfc;
      if(m_sfc->name.find(sfc_lay_name)!=string::npos){
        delete m_sfc;
        surface_vector.remove(*p_sfc);
        break;
      }
   ++p_sfc;
    }
   }
 */

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
void MSHUnselectSFC()
{
	Surface* m_sfc = NULL;
	vector<Surface*>::const_iterator p_sfc;
	p_sfc = surface_vector.begin();
	while (p_sfc != surface_vector.end())
	{
		m_sfc = *p_sfc;
		m_sfc->meshing_allowed = 0;
		++p_sfc;
	}
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
int SFCGetMaxMATGroupNumber()
{
	int sfc_max_mat_group_number = -1;
	Surface* m_sfc = NULL;
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		m_sfc = surface_vector[i];
		if (m_sfc->mat_group > sfc_max_mat_group_number)
			sfc_max_mat_group_number = m_sfc->mat_group;
	}
	return sfc_max_mat_group_number;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   12/2005 OK Color of surface MAT groups
**************************************************************************/
void SFCAssignMATColors()
{
	long r_number;
	vector<long*> aux_vector;
	long* m_color;
	//----------------------------------------------------------------------
	int sfc_max_mat_group_number = SFCGetMaxMATGroupNumber();
	for (int i = 0; i < (int)sfc_max_mat_group_number + 1; i++)
	{
		m_color = new long[3];
		r_number = rand();
		m_color[0] = (long)(r_number / (double)RAND_MAX * 256);
		r_number = rand();
		m_color[1] = (long)(r_number / (double)RAND_MAX * 256);
		r_number = rand();
		m_color[2] = (long)(r_number / (double)RAND_MAX * 256);
		aux_vector.push_back(m_color);
	}
	//----------------------------------------------------------------------
	Surface* m_sfc = NULL;
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		m_sfc = surface_vector[i];
		if (m_sfc->mat_group > 0)
		{
			m_color = aux_vector[m_sfc->mat_group];
			m_sfc->m_color[0] = m_color[0];
			m_sfc->m_color[1] = m_color[1];
			m_sfc->m_color[2] = m_color[2];
		}
	}
	//----------------------------------------------------------------------
	for (int i = 0; i < (int)aux_vector.size(); i++)
	{
		m_color = aux_vector[i];
		delete m_color;
	}
	aux_vector.clear();
	//----------------------------------------------------------------------
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   12/2005 OK Color of surface MAT groups
**************************************************************************/
void SFCAssignColors()
{
	Surface* m_sfc = NULL;
	for (int i = 0; i < (int)surface_vector.size(); i++)
	{
		m_sfc = surface_vector[i];
		m_sfc->AssignColor();
	}
}
