/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// GeoLib.cpp : Definiert den Einsprungpunkt f�r die Konsolenanwendung.
//
/*-------------------------------------------------------------------------
   GeoLib
   Programming:
   09/2003 OK/CC encapsulate Read function
   last modified:
   09/2005 CC GeoLib2
   -------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
// C++
#include <iostream>
#include <string>
using namespace std;
// GL
#include "geo_lib.h"
#include "geo_ply.h"
#include "geo_pnt.h"
#include "geo_sfc.h"
#include "geo_vol.h"

extern void remove_white_space(string* buffer);

enum GEO_TYPES {GEO_POINT = 'P', \
	        GEO_LINE = 'L',
	        GEO_POLYLINE = 'O',
	        GEO_SURFACE = 'S',
	        GEO_VOLUME = 'V',
	        GEO_DOMAIN = 'D'};

/**************************************************************************
   GeoLib- Funktion: GEOLIB_Read_GeoLib
   Aufgabe: Ansteuerung der GeoLib-Lesefunktionen
   Programmaenderungen:
   09/2003 TK Erste Version
   11/2003 OK Read polylines
   11/2003 WW Read surfaces
   09/2004 TK Path info // please check
   10/2004 OK path_name_slash
   09/2005 CC delete line lesen function
**************************************************************************/
void GEOLIB_Read_GeoLib(const std::string &file_name_path_base)
{
	// Points
	GEORemoveAllPoints();
	GEOReadPoints (file_name_path_base);
	// Polylines
	GEORemoveAllPolylines();
	GEOReadPolylines(file_name_path_base);
	// Surfaces
	GEORemoveAllSurfaces();
	GEOReadSurfaces(file_name_path_base);
	// Volumes
	GEORemoveAllVolumes();
	GEOReadVolumes(file_name_path_base);
	// Determine dependencies between GEO objects
	GEOSurfaceTopology();
	GEOCreateSurfacePointVector(); //OK
}

/*************************************************************************
   GeoLib- Funktion: GEOLIB_Clear_GeoLib_Data

   Aufgabe: L�schen aller GeoLib-Datenstrukturen

   Programmaenderungen:
   09/2003     TK        Erste Version
   01/2005 OK GEORemovePolylines, GEORemoveAllSurfaces, GEORemoveVolumes
 **************************************************************************/
void GEOLIB_Clear_GeoLib_Data ()
{
	GEORemoveAllPoints();
	Clear_LineVector();
	GEORemoveAllPolylines(); //OK41
	GEORemoveAllSurfaces(); //OK41 CC change
	GEORemoveAllVolumes(); //OK41
}
/**************************************************************************/
/* GEOLIB - Funktion: GEO_Delete_DoublePoints

   Aufgabe:

   Ergebnis:
   - void -
   Programmaenderungen:
   02/2004    TK      Erste Version
 */
/**************************************************************************/
void GEO_Delete_DoublePoints()
{
	int j = 0, k = 0;
	long pointsvectorsize;
	//WW long pointsvectorsize, linesvectorsize;
	long number_of_polylinepoints;
	long check_point;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	pointsvectorsize = (long) gli_points_vector.size();

	vector<CGLLine*> gli_lines_vector;
	gli_lines_vector = GEOLIB_GetGLILines_Vector();
	//WW linesvectorsize = (long) gli_lines_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	size_t polylinesvectorsize (polyline_vector.size());
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	/*Deleting Parts of Polyline*/
	for (size_t l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();

		number_of_polylinepoints = (long) gl_polyline->point_vector.size();
		for (j = 0; j < number_of_polylinepoints; j++)
		{
			//check_point = gl_polyline->point_vector[j]->old_id;
			//gl_polyline->point_vector[j]->id = gli_points_vector[check_point]->new_id;

			check_point = gl_polyline->point_vector[j]->id;
			check_point = gli_points_vector[check_point]->old_id;
			check_point = gli_points_vector[check_point]->new_id;
			gl_polyline->point_vector[j]->id = check_point;

			//check_point = gli_points_vector[check_point]->old_id;
		}

		for (j = 0; j < number_of_polylinepoints - 1; j++)
		{
			for (k = 0; k < number_of_polylinepoints; k++)
				check_point = gl_polyline->point_vector[k]->id;

			if (gl_polyline->point_vector[j]->id == gl_polyline->point_vector[j
			                                                                  + 1]->id)
			{
				gl_polyline->point_vector.erase(
				        gl_polyline->point_vector.begin() + j + 1);
				j--;
				number_of_polylinepoints
				        = (long) gl_polyline->point_vector.size();
			}
		}
		++p;
	}

	/*Deleting Points*/
	for (j = 0; j < pointsvectorsize; j++)
	{
		check_point = gli_points_vector[j]->new_id;
		check_point = gli_points_vector[j]->old_id;
		check_point = gli_points_vector[j]->first_identical_id;
		check_point = gli_points_vector[j]->id;

		if (gli_points_vector[j]->old_id
		    != gli_points_vector[j]->first_identical_id)
		{
			gli_points_vector.erase(gli_points_vector.begin() + j);
			pointsvectorsize = (long) gli_points_vector.size();
			j--;
		}
		else
			gli_points_vector[j]->id = gli_points_vector[j]->new_id;

	}

	GEOLIB_SetGLIPoints_Vector(gli_points_vector);
	gli_points_vector = GetPointsVector();
	pointsvectorsize = (long) gli_points_vector.size();
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Serialize_Point_Numbers

   Aufgabe: Renumbering of GLI File. Point numbers start with 0.
   Change also the pointnumbers of lines and polylines.

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Serialize_Point_Numbers()
{
	int j = 0, k = 0, l = 0;
	long pointsvectorsize, linesvectorsize, polylinesvectorsize;
	long number_of_polylinepoints;
	long check_point;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	pointsvectorsize = (long)gli_points_vector.size();

	vector<CGLLine*> gli_lines_vector;
	gli_lines_vector = GEOLIB_GetGLILines_Vector();
	linesvectorsize = (long)gli_lines_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	for (j = 0; j < pointsvectorsize; j++)
	{
		check_point = gli_points_vector[j]->new_id = j;
		check_point = gli_points_vector[j]->old_id = gli_points_vector[j]->id;
		check_point = gli_points_vector[j]->id = gli_points_vector[j]->new_id;
	}

	for (k = 0; k < linesvectorsize; k++)
	{
		check_point = gli_lines_vector[k]->point1;
		for (j = 0; j < pointsvectorsize; j++)
			if (check_point == gli_points_vector[j]->old_id)
			{
				gli_lines_vector[k]->point1 = gli_points_vector[j]->new_id;
				break;
			}
		check_point = gli_lines_vector[k]->point2;
		for (j = 0; j < pointsvectorsize; j++)
			if (check_point == gli_points_vector[j]->old_id)
			{
				gli_lines_vector[k]->point2 = gli_points_vector[j]->new_id;
				break;
			}
	}

	for (l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();

		number_of_polylinepoints = (long)gl_polyline->point_vector.size();
		for (j = 0; j < number_of_polylinepoints; j++)
		{
			check_point = gl_polyline->point_vector[j]->old_id;
			//check_point = gl_polyline->point_vector[j]->new_id;
			//check_point = gl_polyline->point_vector[j]->id;
			for (k = 0; k < pointsvectorsize; k++)
				if (check_point == gli_points_vector[k]->old_id)
				{
					gl_polyline->point_vector[j]->id =
					        gli_points_vector[k]->new_id;
					gl_polyline->point_vector[j]->old_id =
					        gli_points_vector[k]->new_id;
					gl_polyline->point_vector[j]->new_id =
					        gli_points_vector[k]->new_id;
					break;
				}
		}
		++p;
	}
	GEOLIB_SetGLIPoints_Vector(gli_points_vector);
	gli_points_vector = GetPointsVector();
	pointsvectorsize = (long)gli_points_vector.size();
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Get_Min_PolySeg_length

   Aufgabe: Min Max Polyline Segment Length for each point

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Get_Min_Max_Distance_of_polyline_neighbor_points()
{
	int i = 0, j = 0, l = 0;
	//WW long pointsvectorsize
	long polylinesvectorsize;
	long number_of_polylinepoints;
	//WW long check_point;
	double x1 = 0.0,y1 = 0.0,z1 = 0.0,seg_length = 0.0;
	double x2 = 0.0,y2 = 0.0,z2 = 0.0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	//WW pointsvectorsize =(long)gli_points_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	for (i = 0; i < (long)gli_points_vector.size(); i++)
	{
		gli_points_vector[i]->min_seg_length = 0.0;
		gli_points_vector[i]->max_seg_length = 0.0;
	}

	for (i = 0; i < (long)gli_points_vector.size(); i++)
	{
		p = polyline_vector.begin();

		for (l = 0; l < polylinesvectorsize; l++)
		{
			gl_polyline = *p;
			Name = gl_polyline->getName();

			number_of_polylinepoints = (long)gl_polyline->point_vector.size();
			for (j = 0; j < number_of_polylinepoints; j++)
				if (gli_points_vector[i]->id ==  gl_polyline->point_vector[j]->id)
				{
					if (j < number_of_polylinepoints - 1)
					{
						//WW check_point = gl_polyline->point_vector[j]->id;
						x1 = gl_polyline->point_vector[j]->x;
						y1 = gl_polyline->point_vector[j]->y;
						z1 = gl_polyline->point_vector[j]->z;
						//WW check_point = gl_polyline->point_vector[j+1]->id;
						x2 = gl_polyline->point_vector[j + 1]->x;
						y2 = gl_polyline->point_vector[j + 1]->y;
						z2 = gl_polyline->point_vector[j + 1]->z;
						seg_length = EuklVek3dDistCoor ( x1,
						                                 y1,
						                                 z1,
						                                 x2,
						                                 y2,
						                                 z2 );
						if (gli_points_vector[i]->min_seg_length == 0.0)
							gli_points_vector[i]->min_seg_length =
							        seg_length;
						if (gli_points_vector[i]->max_seg_length == 0.0)
							gli_points_vector[i]->max_seg_length =
							        seg_length;

						if (gli_points_vector[i]->min_seg_length >
						    seg_length)
							gli_points_vector[i]->min_seg_length =
							        seg_length;
						if (gli_points_vector[i]->max_seg_length <
						    seg_length)
							gli_points_vector[i]->max_seg_length =
							        seg_length;
					}
					if (j > 0)
					{
						//WW check_point = gl_polyline->point_vector[j]->id;
						x1 = gl_polyline->point_vector[j]->x;
						y1 = gl_polyline->point_vector[j]->y;
						z1 = gl_polyline->point_vector[j]->z;
						//WW check_point = gl_polyline->point_vector[j-1]->id;
						x2 = gl_polyline->point_vector[j - 1]->x;
						y2 = gl_polyline->point_vector[j - 1]->y;
						z2 = gl_polyline->point_vector[j - 1]->z;
						seg_length = EuklVek3dDistCoor ( x1,
						                                 y1,
						                                 z1,
						                                 x2,
						                                 y2,
						                                 z2 );
						if (gli_points_vector[i]->min_seg_length == 0.0)
							gli_points_vector[i]->min_seg_length =
							        seg_length;
						if (gli_points_vector[i]->max_seg_length == 0.0)
							gli_points_vector[i]->max_seg_length =
							        seg_length;

						if (gli_points_vector[i]->min_seg_length >
						    seg_length)
							gli_points_vector[i]->min_seg_length =
							        seg_length;
						if (gli_points_vector[i]->max_seg_length <
						    seg_length)
							gli_points_vector[i]->max_seg_length =
							        seg_length;
					}
				}

			++p;
		}
	}
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Get_Min_Distance_of_neighbor_points

   Aufgabe: Min Distance for all points

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Get_Min_Distance_of_neighbor_points()
{
	int i = 0, j = 0;
	//WW long check_point;
	double x1 = 0.0,y1 = 0.0,z1 = 0.0,seg_length = 0.0;
	double x2 = 0.0,y2 = 0.0,z2 = 0.0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();

	for (i = 0; i < (long)gli_points_vector.size(); i++)
		gli_points_vector[i]->min_seg_length = 0.0;


	for (i = 0; i < (long)gli_points_vector.size(); i++)
	{
		for (j = 0; j < (long)gli_points_vector.size(); j++)

			if (gli_points_vector[i]->id !=  gli_points_vector[j]->id)
			{
				//WW check_point = gli_points_vector[i]->id;
				x1 = gli_points_vector[i]->x;
				y1 = gli_points_vector[i]->y;
				z1 = gli_points_vector[i]->z;
				//WW check_point = gli_points_vector[j]->id;
				x2 = gli_points_vector[j]->x;
				y2 = gli_points_vector[j]->y;
				z2 = gli_points_vector[j]->z;
				seg_length = EuklVek3dDistCoor ( x1, y1, z1, x2, y2, z2 );
				if (gli_points_vector[i]->min_seg_length == 0.0)
					gli_points_vector[i]->min_seg_length = seg_length;
				if (gli_points_vector[i]->min_seg_length > seg_length)
					gli_points_vector[i]->min_seg_length = seg_length;
			}


	}
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Get_Min_PolySeg_length

   Aufgabe: Min Max Polyline Segment Length for each point

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Copy_Min_Distance_Of_Neighbor_Points_To_Mesh_Density()
{
	int i = 0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();

	/*for (i=0;i<(long)gli_points_vector.size();i++)
	   {
	    gli_points_vector[i]->mesh_density = gli_points_vector[i]->min_seg_length;
	   }*/

	for (i = 0; i < (long)gli_points_vector.size(); i++)
		if (gli_points_vector[i]->mesh_density > gli_points_vector[i]->min_seg_length)
			gli_points_vector[i]->mesh_density = gli_points_vector[i]->min_seg_length;

}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Mesh_Density_BiSect

   Aufgabe: halving of the mesh density (maybe in meshlib)

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Mesh_Density_BiSect()
{
	int i = 0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();

	for (i = 0; i < (long)gli_points_vector.size(); i++)
		gli_points_vector[i]->mesh_density = 0.5 * gli_points_vector[i]->mesh_density;
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Get_Min_PolySeg_length

   Aufgabe: Analyse Min Max Segment Length

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
double GEO_Get_Min_PolySeg_length()
{
	int j = 0, l = 0, hit = 0;
	//WW long pointsvectorsize,
	long polylinesvectorsize;
	long number_of_polylinepoints;
	//WW long check_point;
	double x1 = 0.0,y1 = 0.0,z1 = 0.0,seg_length = 0.0;
	double x2 = 0.0,y2 = 0.0,z2 = 0.0,seg_length_saved = 0.0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	//WW pointsvectorsize =(long)gli_points_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	for (l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();

		number_of_polylinepoints = (long)gl_polyline->point_vector.size();
		for (j = 0; j < number_of_polylinepoints - 1; j++)
		{
			//WW check_point = gl_polyline->point_vector[j]->id;
			x1 = gl_polyline->point_vector[j]->x;
			y1 = gl_polyline->point_vector[j]->y;
			z1 = gl_polyline->point_vector[j]->z;

			//WW check_point = gl_polyline->point_vector[j+1]->id;
			x2 = gl_polyline->point_vector[j + 1]->x;
			y2 = gl_polyline->point_vector[j + 1]->y;
			z2 = gl_polyline->point_vector[j + 1]->z;

			seg_length = EuklVek3dDistCoor ( x1, y1, z1, x2, y2, z2 );

			if (j == 0 && hit == 0)
			{
				seg_length_saved = seg_length;
				hit++;
			}
			else if (seg_length < seg_length_saved)
				seg_length_saved = seg_length;
		}
		++p;
	}
	return seg_length_saved;
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Get_Min_PolySeg_length

   Aufgabe: Analyse Min Max Segment Length

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
double GEO_Get_Max_PolySeg_length()
{
	int j = 0, l = 0, hit = 0;
	//WW long pointsvectorsize,
	long polylinesvectorsize;
	long number_of_polylinepoints;
	//WW long check_point;
	double x1 = 0.0,y1 = 0.0,z1 = 0.0,seg_length = 0.0;
	double x2 = 0.0,y2 = 0.0,z2 = 0.0,seg_length_saved = 0.0;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	//WW pointsvectorsize =(long)gli_points_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	for (l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();

		number_of_polylinepoints = (long)gl_polyline->point_vector.size();
		for (j = 0; j < number_of_polylinepoints - 1; j++)
		{
			//WW check_point = gl_polyline->point_vector[j]->id;
			x1 = gl_polyline->point_vector[j]->x;
			y1 = gl_polyline->point_vector[j]->y;
			z1 = gl_polyline->point_vector[j]->z;

			//WW check_point = gl_polyline->point_vector[j+1]->id;
			x2 = gl_polyline->point_vector[j + 1]->x;
			y2 = gl_polyline->point_vector[j + 1]->y;
			z2 = gl_polyline->point_vector[j + 1]->z;

			seg_length = EuklVek3dDistCoor ( x1, y1, z1, x2, y2, z2 );

			if (j == 0 && hit == 0)
			{
				seg_length_saved = seg_length;
				hit++;
			}
			else if (seg_length > seg_length_saved)
				seg_length_saved = seg_length;
		}
		++p;
	}
	return seg_length_saved;
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Polylines_per_Point

   Aufgabe:     Number of polylines which share the point

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Polylines_per_Point()
{
	int j = 0, l = 0;
	int nb_of_ply;
	//WW long pointsvectorsize,
	long polylinesvectorsize;
	long number_of_polylinepoints;
	long check_point;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	//WW pointsvectorsize =(long)gli_points_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	for (l = 0; l < (long)gli_points_vector.size(); l++)
		gli_points_vector[l]->nb_of_ply = 0;


	for (l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();
		number_of_polylinepoints = (long)gl_polyline->point_vector.size();
		for (j = 0; j < number_of_polylinepoints; j++)
		{
			check_point = gl_polyline->point_vector[j]->id;
			nb_of_ply = gli_points_vector[check_point]->nb_of_ply;
			nb_of_ply = gli_points_vector[check_point]->nb_of_ply = nb_of_ply + 1;

			if(gl_polyline->point_vector[j]->id == gl_polyline->point_vector[0]->id &&
			   j != 0)
				gli_points_vector[check_point]->nb_of_ply = nb_of_ply - 1;

			nb_of_ply = gli_points_vector[check_point]->nb_of_ply;
		}
		++p;
	}
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Set_Poly_Seg_Length

   Aufgabe: Setzt die L�nge der Polyliniensegmente

   Ergebnis:
   - void -
   Programmaenderungen:
   08/2005    TK      Erste Version
 */
/**************************************************************************/
void GEO_Set_Poly_Seg_Length(double min_seg_length, double max_seg_length)
{
	int j = 0, l = 0, k = 0, i = 0;
	int j_old = 0;
	//long pointsvectorsize,
	long polylinesvectorsize;
	//WW long number_of_polylinepoints;
	long check_point;
	int nb_of_ply, nb_of_ply_1st, nb_of_ply_2nd;
	double x1 = 0.0,y1 = 0.0,z1 = 0.0,seg_length = 0.0,seg_length_old = 0.0;
	double x2 = 0.0,y2 = 0.0,z2 = 0.0;
	double e_x,e_y,e_z,pmax_x,pmax_y,pmax_z; //WW ,pmin_x,pmin_y,pmin_z;
	long check_point_id_1,check_point_id_2,check_point_id_3,check_point_id_4;
	CGLPoint* m_point = NULL;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	//WW pointsvectorsize =(long)gli_points_vector.size();

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector;
	polyline_vector = GetPolylineVector();
	polylinesvectorsize = (long)polyline_vector.size();
	vector<CGLPolyline*>::const_iterator p = polyline_vector.begin();

	string Name;

	/*Number of polylines which share the point*/
	GEO_Polylines_per_Point();

	/*Adding and deleting*/
	for (l = 0; l < polylinesvectorsize; l++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();

		for (j = 0; j < (long)gl_polyline->point_vector.size() - 1; j++)
		{
			/*Calculation Segment length*/
			check_point = gl_polyline->point_vector[j]->id;
			x1 = gl_polyline->point_vector[j]->x;
			y1 = gl_polyline->point_vector[j]->y;
			z1 = gl_polyline->point_vector[j]->z;
			check_point = gl_polyline->point_vector[j + 1]->id;
			x2 = gl_polyline->point_vector[j + 1]->x;
			y2 = gl_polyline->point_vector[j + 1]->y;
			z2 = gl_polyline->point_vector[j + 1]->z;
			seg_length = EuklVek3dDistCoor ( x1, y1, z1, x2, y2, z2 );

			if (seg_length_old == seg_length && j_old == j)
				break;
			j_old = j;
			seg_length_old = seg_length;
			/*Calculation Einheitsvektor*/
			e_x = (x2 - x1) / seg_length;
			e_y = (y2 - y1) / seg_length;
			e_z = (z2 - z1) / seg_length;
			/*Calculation Point with max distance*/
			pmax_x = x1 + (max_seg_length * e_x);
			pmax_y = y1 + (max_seg_length * e_y);
			pmax_z = z1 + (max_seg_length * e_z);
			/*Calculation Point with min distance*/
			//pmin_x = x1+(min_seg_length* e_x);
			//pmin_y = y1+(min_seg_length* e_y);
			//pmin_z = z1+(min_seg_length* e_z);

			if (seg_length < min_seg_length) /*Delete */
			{
				if (gl_polyline->point_vector[j + 1]->id !=
				    gl_polyline->point_vector[0]->id)
				{
					check_point = gl_polyline->point_vector[j + 1]->id;
					nb_of_ply = gli_points_vector[check_point]->nb_of_ply;
					gl_polyline->point_vector.erase(
					        gl_polyline->point_vector.begin() + j + 1);
					j--;
					//WW number_of_polylinepoints = (long)gl_polyline->point_vector.size();
				}
				else
				{
					check_point = gl_polyline->point_vector[j]->id;
					nb_of_ply = gli_points_vector[check_point]->nb_of_ply;
					gl_polyline->point_vector.erase(
					        gl_polyline->point_vector.begin() + j);
					j--;
					//WW number_of_polylinepoints = (long)gl_polyline->point_vector.size();
				}

				if (nb_of_ply > 1)
				{
					vector<CGLPolyline*>::const_iterator pp =
					        polyline_vector.begin();

					for (i = 0; i < (long)polyline_vector.size(); i++)
					{
						gl_polyline = *pp;
						for (k = 0;
						     k < (long)gl_polyline->point_vector.size();
						     k++)
							if (gl_polyline->point_vector[k]->id ==
							    check_point)
							{
								gl_polyline->point_vector.erase(
								        gl_polyline->point_vector.
								        begin() +
								        k);
								k--;
							}
						pp++;
					}
				}
				gl_polyline = *p;
			}

			if (seg_length > max_seg_length) /*Add*/
			{
				check_point = gl_polyline->point_vector[j]->id;
				nb_of_ply_1st = gli_points_vector[check_point]->nb_of_ply;
				check_point = gl_polyline->point_vector[j + 1]->id;
				nb_of_ply_2nd = gli_points_vector[check_point]->nb_of_ply;

				check_point_id_1 = gl_polyline->point_vector[j]->id;
				check_point_id_2 = gl_polyline->point_vector[j + 1]->id;

				m_point = new CGLPoint();
				m_point->id = (long)gli_points_vector.size();
				m_point->x = pmax_x;
				m_point->y = pmax_y;
				m_point->z = pmax_z;
				m_point->nb_of_ply = 1;
				gli_points_vector.push_back(m_point);

				gl_polyline->point_vector.insert(
				        gl_polyline->point_vector.begin() + j + 1,
				        gli_points_vector.end() - 1,gli_points_vector.end());

				//WW number_of_polylinepoints = (long)gl_polyline->point_vector.size();

				if (nb_of_ply_1st > 1 && nb_of_ply_2nd > 1)
				{
					vector<CGLPolyline*>::const_iterator pp =
					        polyline_vector.begin();

					for (i = 0; i < (long)polyline_vector.size(); i++)
					{
						gl_polyline = *pp;
						for (k = 0;
						     k < (long)gl_polyline->point_vector.size() - 1;
						     k++)
						{
							check_point_id_3 =
							        gl_polyline->point_vector[k]->id;
							check_point_id_4 =
							        gl_polyline->point_vector[k +
							                                  1]->id;

							if ((check_point_id_3 ==
							     check_point_id_1 &&
							     check_point_id_4 ==
							     check_point_id_2) ||
							    (check_point_id_3 ==
							     check_point_id_2 &&
							     check_point_id_4 == check_point_id_1))
							{
								gl_polyline->point_vector.insert(
								        gl_polyline->point_vector.
								        begin() +
								        k + 1,
								        gli_points_vector.end() - 1,
								        gli_points_vector.end());

								check_point =
								        (long)gli_points_vector.
								        size() - 1;
								check_point =
								        gli_points_vector[
								                check_point]->id;
								check_point =
								        (long)gli_points_vector.
								        size() - 1;
								nb_of_ply =
								        gli_points_vector[
								                check_point]->
								        nb_of_ply;
								check_point =
								        (long)gli_points_vector.
								        size() - 1;
								gli_points_vector[check_point]->
								nb_of_ply = nb_of_ply++;

								k--;
							}
						}
						pp++;
					}
				}
				gl_polyline = *p;
				j--;
				if((pmax_x == x2 && pmax_y == y2 &&
				    pmax_z == z2) || (pmax_x == x1 && pmax_y == y1 && pmax_z == z1))
					j++;
			}
		}
		++p;
	}

	GEOLIB_SetGLIPoints_Vector(gli_points_vector);

	/*Deleting the points in the point vector*/
	GEO_Polylines_per_Point();
	for (i = 0; i < (long)gli_points_vector.size(); i++)
		if (gli_points_vector[i]->nb_of_ply < 1)
		{
			gli_points_vector.erase(gli_points_vector.begin() + i);
			i--;
		}

	GEOLIB_SetGLIPoints_Vector(gli_points_vector);
}

/**************************************************************************/
/* GEOLIB - Funktion: GEO_Write_GMSH_Input_File()

   Aufgabe:

   Ergebnis:
   - void -
   Programmaenderungen:
   02/2004    TK      Erste Version
 */
/**************************************************************************/
void GEO_Write_GMSH_Input_File(const char* file_name)
{
	FILE* geo_file = NULL;
	geo_file = fopen(file_name, "w+t");

	/*GMSH-Properties:*/
	//fprintf(geo_file,"%s\n","Mesh.Algorithm = 3;");
	fprintf(geo_file,"%s\n","Mesh.Smoothing = 8;");

	/*Points:*/
	long pointsvectorsize;
	int i = 0, k = 0;
	long point1, point2;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector();
	pointsvectorsize = (long)gli_points_vector.size();
	fprintf(geo_file,"%s\n","//Points");
	for (i = 0; i < pointsvectorsize; i++)
	{
		fprintf(geo_file,"%s","Point(");
		fprintf(geo_file,"%i",i);
		fprintf(geo_file,"%s",") = {");
		fprintf(geo_file,"%g",gli_points_vector[i]->x);
		fprintf(geo_file,"%s",", ");
		fprintf(geo_file,"%g",gli_points_vector[i]->y);
		fprintf(geo_file,"%s",", ");
		fprintf(geo_file,"%g",gli_points_vector[i]->z);
		fprintf(geo_file,"%s",", ");
		fprintf(geo_file,"%g",gli_points_vector[i]->mesh_density);
		fprintf(geo_file,"%s\n","};");
	}
	/*Lines:*/
	long linesvectorsize;
	int j = 0;
	vector<CGLLine*> gli_lines_vector;
	gli_lines_vector = GEOLIB_GetGLILines_Vector();
	linesvectorsize = (long)gli_lines_vector.size();
	fprintf(geo_file,"%s\n","//Lines");

	/*Polylines:*/
	long check_point;
	long polylinesvectorsize, number_of_polylinepoints, number_of_polylinelines,
	     surfacepolyline_vectorlength;
	int line_id_counter = 0;
	int line_id_counter_begin = 0;
	string Name;

	CGLPolyline* gl_polyline = NULL;
	vector<CGLPolyline*> polyline_vector; //CC
	polyline_vector = GetPolylineVector(); //CC
	polylinesvectorsize = (long)polyline_vector.size(); //CC
	vector<CGLPolyline*>::iterator p = polyline_vector.begin(); //CC
	fprintf(geo_file,"%s\n","//Polylines");

	for (i = 0; i < polylinesvectorsize; i++)
	{
		gl_polyline = *p;
		Name = gl_polyline->getName();
		number_of_polylinepoints = (long)gl_polyline->point_vector.size();
		number_of_polylinelines = (long)((gl_polyline->getLineVector()).size());
		/*if $Points*/
		if (number_of_polylinepoints > 0)
		{
			for (j = 0; j < number_of_polylinepoints - 1; j++)
			{
				line_id_counter++;
				if (j == 0)
					line_id_counter_begin = line_id_counter;
				point1 = gl_polyline->point_vector[j]->id;
				point2 = gl_polyline->point_vector[j + 1]->id;
				fprintf(geo_file,"%s","Line(");
				fprintf(geo_file,"%i",line_id_counter);
				fprintf(geo_file,"%s",") = {");
				fprintf(geo_file,"%ld",point1);
				fprintf(geo_file,"%s",", ");
				fprintf(geo_file,"%ld",point2);
				fprintf(geo_file,"%s\n","};");
			}

			fprintf(geo_file,"%s","Line Loop(");
			fprintf(geo_file,"%i",i);
			fprintf(geo_file,"%s",") = {");
			line_id_counter = line_id_counter_begin;
			gl_polyline->setID(i);

			for (j = 0; j < number_of_polylinepoints - 1; j++)
			{
				fprintf(geo_file,"%i",line_id_counter);
				line_id_counter++;
				if (j < number_of_polylinepoints - 2)
					fprintf(geo_file,"%s",", ");
			}
			fprintf(geo_file,"%s\n","};");
		}
		/*if $Lines*/ //Todo: benchmark test
		if (number_of_polylinelines > 0)
		{
			for (j = 0; j < linesvectorsize; j++)
			{
				fprintf(geo_file,"%s","Line(");
				fprintf(geo_file,"%ld",gli_lines_vector[i]->gli_line_id);
				fprintf(geo_file,"%s",") = {");
				fprintf(geo_file,"%ld",gli_lines_vector[i]->point1);
				fprintf(geo_file,"%s",", ");
				fprintf(geo_file,"%ld",gli_lines_vector[i]->point2);
				fprintf(geo_file,"%s\n","};");
			}
			for (j = 0; j < number_of_polylinelines; j++)
			{
				check_point = (gl_polyline->getLineVector())[j]->gli_line_id;
				fprintf(geo_file,"%ld",check_point);
				if (j < number_of_polylinelines - 1)
					fprintf(geo_file,"%s",", ");
			}
			fprintf(geo_file,"%s\n","};");
		}
		++p;
	}
	line_id_counter = 0;

	/*Surfaces:*/
	long surfacesvectorsize;
	Surface* gl_surface = NULL;
	vector<Surface*> surface_vector;
	CGLPolyline* sfc_polyline = NULL;
	surface_vector = GetSurfaceVector(); //CC
	surfacesvectorsize = (long)surface_vector.size(); //CC
	vector<Surface*>::iterator ps = surface_vector.begin();

	fprintf(geo_file,"%s\n","//Surfaces");
	for (i = 0; i < surfacesvectorsize; i++)
	{
		gl_surface = *ps;
		fprintf(geo_file,"%s","Plane Surface(");
		fprintf(geo_file,"%i",i);
		fprintf(geo_file,"%s",") = {");

		//list<CGLPolyline*>::const_iterator pp = gl_surface->polyline_of_surface_list.begin();//CC
		vector<CGLPolyline*>::iterator pp = gl_surface->polyline_of_surface_vector.begin();
		surfacepolyline_vectorlength = (long)gl_surface->polyline_of_surface_vector.size(); //CC
		vector<CGLPolyline*>::iterator p2 = polyline_vector.begin(); //CC
		sfc_polyline = *pp;
		for (j = 0; j < surfacepolyline_vectorlength; j++)
		{
			for (k = 0; k < polylinesvectorsize; k++)
			{
				gl_polyline = *p2;
				if (sfc_polyline->getName() == gl_polyline->getName())
					fprintf(geo_file,"%ld",
					        static_cast<long>(gl_polyline->getID()));
				++p2;
			}

			if (j < surfacepolyline_vectorlength - 1)
				fprintf(geo_file,"%s",", ");
			++pp;
		}
		fprintf(geo_file,"%s\n","};");
		++ps;
	}

	fclose(geo_file);
}

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   03/2004 TK Implementation
   Task: filenamebase = Pfadname ohne extension
   01/2005 OK Volumes
   09/2005 CC surface and volume write function
   11/2005 CC Write
**************************************************************************/
void GEOWrite(string file_name)
{
	char file_name_char[1024];
	strcpy(file_name_char,file_name.c_str());
	// Points
	GEOWritePoints(file_name_char); //CC
	// Polylines
	GEOWritePolylines(file_name_char); //CC
	// Surfaces
	GEOWriteSurfaces(file_name);
	// Volumes
	GEOWriteVolumes(file_name);
	// Domain
	// STOP
	FILE* gli_file = NULL;
	string gli_file_name;
	gli_file_name = file_name + ".gli";
	const char* gli_file_name_char = 0;
	gli_file_name_char = gli_file_name.data();
	gli_file = fopen(gli_file_name_char, "a");
	fprintf(gli_file,"%s","#STOP");
	fclose(gli_file);
}
