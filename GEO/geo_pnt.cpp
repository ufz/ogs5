/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Point
   Task:
   Programing:
   07/2003 OK/TK Implementation
   12/2003 CC Modification
   08/2005 CC Modification
   09/2005 CC GeoLib2
**************************************************************************/

#include <cstdlib>
// GeoLib
#include "geo_pnt.h"
#include "mathlib.h"

using namespace std;

/*----------------------------------------------------------------------*/
std::vector<CGLPoint*> gli_points_vector;
std::vector<CGLPoint*> pnt_properties_vector;
double pnt_x_min;
double pnt_x_max;
double pnt_y_min;
double pnt_y_max;
double pnt_z_min;
double pnt_z_max;
long pnt_id_last = 0;
/*----------------------------------------------------------------------*/
// constructor
/**************************************************************************
   GEOLib-Method:
   Programing:
   10/2005 OK PNT name
**************************************************************************/
CGLPoint::CGLPoint(void) : x(data[0]), y(data[1]), z(data[2]), _propert(0.0)
{
	index_msh = -1;
	circle_pix = 3; // CC
	mesh_density = 100.;
	value = 0.0;
	data[0] = 0.0;
	data[1] = 0.0;
	data[2] = 0.0;
	//  highlighted = false;
	//  m_color[0] = 0;
	//  m_color[1] = 0;
	//  m_color[2] = 0;
	x_pix = 0;
	y_pix = 0;
	mat = -1; // CC9999 // material properties
	// PNT name //OK
	//......................................................................
	name = "POINT";
	std::ostringstream os; // TF 01/2010
	if (gli_points_vector.size() > 0)
		os << (*(--gli_points_vector.end()))->id + 1; // TF 01/2010
	else
		os << 0; // TF 01/2010
	name += os.str(); // TF 01/2010
	//......................................................................
}

CGLPoint::CGLPoint(double x1, double y1, double z1) : x(data[0]), y(data[1]), z(data[2])
{
	data[0] = x1;
	data[1] = y1;
	data[2] = z1;
}

CGLPoint::CGLPoint(const double* coordinates) : x(data[0]), y(data[1]), z(data[2])
{
	for (size_t k(0); k < 3; k++)
		data[k] = coordinates[k];
}

// deconstructor
CGLPoint::~CGLPoint(void)
{
}

/**************************************************************************
   GeoLib-Method:
   Task: Check if the point exists in .gli file.
   Description: If there is the same point in the file, the function return
   the index of the same point. If not, the function returns -100.
   Programing:
   03/2005 PCH Implementation
   05/2005 PCH Modified to make this function independent with GUI
   08/2005 CC Modification Move from Geo_ply to Geo_Pnt
**************************************************************************/
int CGLPoint::IsPointExist()
{
	// The tolerance for comparing two doubles
	double tolerance = 1e-6; // This is very subjective in digital numbers that represent real numbers.
	// So, what ever points in what scale, always use OpenGL coordinates here.
	int YesThereIs = -100; // Set this to mean there is no duplicate in the file

	int sizeOfPointList = (int)gli_points_vector.size();
	for (int i = 0; i < sizeOfPointList; ++i)
	{
		double Same = sqrt((gli_points_vector[i]->x - x) * (gli_points_vector[i]->x - x)
		                   + (gli_points_vector[i]->y - y) * (gli_points_vector[i]->y - y)
		                   + (gli_points_vector[i]->z - z) * (gli_points_vector[i]->z - z));
		// If the tolerance is bigger than same,
		if (Same < tolerance)
			YesThereIs = gli_points_vector[i]->id;
	}

	return YesThereIs;
}
/**************************************************************************
   GeoLib-Method: PointDis
   Task:two points distance calculation
   Programing:
   11/2003 CC Implementation
   04/2006 CC z = 0
   11/2006 TK Sorry but we can't delete one dimension here
**************************************************************************/
double CGLPoint::PointDis(CGLPoint* m_p2)
{
	return sqrt((x - m_p2->x) * (x - m_p2->x) + (y - m_p2->y) * (y - m_p2->y) + (z - m_p2->z) * (z - m_p2->z));
	// return sqrt((x-m_p2->x)*(x-m_p2->x)+(y-m_p2->y)*(y-m_p2->y)); //CC
}
/**************************************************************************
   GeoLib-Method: Get
   Task: select point instance from list by  id
   Programing:
   07/2005 CC Implementation
**************************************************************************/
CGLPoint* GEOGetPointById(long number)
{
	if (number < 0)
		return NULL;
	else
	{
		std::vector<CGLPoint*>::iterator p = gli_points_vector.begin();
		CGLPoint* m_point = NULL;
		while (p != gli_points_vector.end())
		{
			m_point = *p;
			if (m_point->id == number)
				return m_point;
			++p;
		}
	}
	return NULL;
}
/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   07/2005 CC Implementation
   10/2005 OK Test
**************************************************************************/
//// 06/2010 TF REMOVE CANDIDATE
// CGLPoint* GEOGetPointByName(const std::string &point_name)
//{
//	//CGLPoint *m_pnt = NULL;
//	std::vector<CGLPoint*>::iterator it (gli_points_vector.begin());
//	for (; it != gli_points_vector.end(); it++)
//		if ((*it)->name.compare (point_name) == 0) return *it;
//	std::cout << "Warning: point " << point_name << " not found" << "\n";
//	return NULL;
//}

/**************************************************************************
   GeoLib- Funktion: GetPointsVector

   Aufgabe: Pointer fuer externen Vektorenzugriff

   Programmaenderungen:
   08/2003     TK        Erste Version
**************************************************************************/
std::vector<CGLPoint*> GetPointsVector(void)
{
	return gli_points_vector;
}

/**************************************************************************
   GeoLib- Funktion: GEOLIB_SetGLIPoints_Vector

   Aufgabe: Pointer fuer externen Vektorenzugriff

   Programmaenderungen:
   03/2004     TK        Erste Version
**************************************************************************/
std::vector<CGLPoint*> GEOLIB_SetGLIPoints_Vector(std::vector<CGLPoint*> gl_point)
{
	gli_points_vector = gl_point;
	return gli_points_vector;
}
/**************************************************************************
   GeoLib- Funktion: GEOReadPoints
   Aufgabe: Lesen der GLI Points und schreiben in einen Vector
   08/2005 CC Implementation
**************************************************************************/
void GEOReadPoints(const std::string& file_name_path_base)
{
	string gli_file_name;
	char line[MAX_ZEILEN];
	string line_string;
	gli_file_name = file_name_path_base + ".gli";
	int ok = 1;
	std::stringstream in; // OK
	ios::pos_type position;
	CGLPoint* m_gli_points = NULL;
	ifstream dat_in(gli_file_name.data(), ios::in);
	if (!dat_in.good())
		return;
	dat_in.seekg(0L, ios::beg); // rewind?
	//----------------------------------------------------------------------
	while (!dat_in.eof())
	{
		dat_in.getline(line, MAX_ZEILEN);
		line_string = line;

		if (line_string.find("#STOP") != string::npos) // 11.08.2011. WW
			break;
		if (line_string.find("#POINTS") != string::npos)
			while (ok)
			{
				m_gli_points = new CGLPoint();
				position = m_gli_points->Read(&dat_in, ok);
				if (ok && m_gli_points->id >= 0) // CC8888
					gli_points_vector.push_back(m_gli_points);
				else // CC8888
					delete m_gli_points; // CC8888
				dat_in.seekg(position, ios::beg);
			}
	}
	//----------------------------------------------------------------------
	dat_in.close(); // OK41
}

/**************************************************************************
   GeoLib-Method: Read
   Task:
   Programing:
   08/2005 CC Implementation
   01/2010 TF few changess
**************************************************************************/
ios::pos_type CGLPoint::Read(ifstream* gli_file, int& ok)
{
	std::stringstream in;
	char line[MAX_ZEILEN];
	ios::pos_type position;
	position = gli_file->tellg();
	gli_file->getline(line, MAX_ZEILEN);
	std::string line_string(line);
	if (line_string.find("#") != std::string::npos)
	{
		ok = 0;
		return position;
	}
	in.str(line);
	ok = sscanf(line, "%ld", &id);
	if (ok == 1)
	{
		in >> id >> x >> y >> z;
		size_t pos1 = 0;
		if (line_string.find("$MD") != std::string::npos)
		{
			pos1 = line_string.find_first_of("M");
			in.str(line_string.substr(pos1 + 2, std::string::npos));
			in >> mesh_density;
		}
		if (line_string.find("$ID") != std::string::npos) // OK
		{
			pos1 = line_string.find_first_of("I");
			in.str(line_string.substr(pos1 + 2, std::string::npos));
			in >> name;
		}
		in.clear();
	}
	return gli_file->tellg();
}
/*************************************************************************
   GeoLib- Funktion: Clear_PointVector

   Aufgabe: Leeren des GLI Points-Vektors

   Programmaenderungen:
   08/2003     TK        Erste Version
   03/2006 CC Modification
   05/2006 TK BUGFIX
 **************************************************************************/
void GEORemoveAllPoints()
{
	// CGLPoint * m_pnt = NULL;
	for (int i = 0; i < (int)gli_points_vector.size(); i++)
	{
		// m_pnt = gli_points_vector[0]; TK What's that?
		// delete m_pnt;
		delete gli_points_vector[i];
		gli_points_vector[i] = NULL;
	}
	gli_points_vector.clear();
}
/**************************************************************************
   GeoLib-Method: GEOPointRemove
   Task: remove point from the list
   Programing:
   10/2003 CC Implementation
   03/2006 CC Modification destructor
**************************************************************************/
void GEORemovePoint(long nSel)
{
	CGLPoint* m_pnt = NULL;
	m_pnt = gli_points_vector[nSel];
	delete m_pnt;
	gli_points_vector.erase(gli_points_vector.begin() + nSel);
}
/**************************************************************************/
/* GEOLIB - Funktion: GEO_Search_DoublePoints

   Aufgabe:

   Ergebnis:
   - void -
   Programmaenderungen:
   02/2004    TK      Erste Version
 */
/**************************************************************************/
void GEO_Search_DoublePoints(double tolerance)
{
	int j = 0, k = 0;
	long pointsvectorsize;
	long first_id, double_id;
	double x2check = 0.0, y2check = 0.0, z2check = 0.0;
	double x2check0 = 0.0, y2check0 = 0.0, z2check0 = 0.0;
	long hits = 0, allhits = 0;
	long numbering = 0;
	double point_distance;
	string coordinate;
	string new_coordinate;
	vector<CGLPoint*> gli_points_vector;
	gli_points_vector = GetPointsVector(); // CC
	pointsvectorsize = (long)gli_points_vector.size();
	if (tolerance < MKleinsteZahlen)
		tolerance = 0.0000000000000001;
	for (j = 0; j < pointsvectorsize; j++)
	{
		double_id = gli_points_vector[j]->id;
		first_id = double_id;
		x2check = gli_points_vector[j]->x;
		y2check = gli_points_vector[j]->y;
		z2check = gli_points_vector[j]->z;
		hits = 0;
		for (k = 0; k < j; k++)
		{
			x2check0 = gli_points_vector[k]->x;
			y2check0 = gli_points_vector[k]->y;
			z2check0 = gli_points_vector[k]->z;

			point_distance = EuklVek3dDistCoor(x2check, y2check, z2check, x2check0, y2check0, z2check0);

			if (point_distance < tolerance)
			{
				first_id = gli_points_vector[k]->first_identical_id;
				hits++;
				allhits++;
			}
			if (hits > 0 && fabs(x2check0 - x2check) < MKleinsteZahlen && fabs(y2check0 - y2check) < MKleinsteZahlen
			    && fabs(z2check0 - z2check) < MKleinsteZahlen)
				if (hits > 1)
				{
					hits++;
					allhits++;
				}
		}

		gli_points_vector[j]->first_identical_id = first_id;
		gli_points_vector[j]->number_of_doubled_points = hits;
		if (hits > 0)
			gli_points_vector[j]->new_id = first_id;
		else
		{
			gli_points_vector[j]->new_id = numbering;
			numbering++;
		}
		gli_points_vector[j]->old_id = gli_points_vector[j]->id;
	}

	GEOLIB_SetGLIPoints_Vector(gli_points_vector);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK/JG Implementation
   03/2005 OK Extensions
**************************************************************************/
void GEOReadPointProperties(const std::string& file_name_base)
{
	CGLPoint* m_pnt = NULL;
	char line[MAX_ZEILEN];
	std::string sub_string;
	std::string line_string;
	//  ios::pos_type position;
	std::vector<std::string> pnt_properties_name_vector;
	string delimiter_type(";");
	int pos1, pos2;
	char pnt_name[80];

	cout << "Read PNT properties from " << file_name_base << "_springs"
	     << "\n";
	//========================================================================
	// File handling
	string csv_file_name = file_name_base + CSV_FILE_EXTENSIONS; // OK4105
	ifstream csv_file(csv_file_name.data(), ios::in);
	if (!csv_file.good())
		return;
	csv_file.seekg(0L, ios::beg);
	//========================================================================
	// Keyword loop
	csv_file.getline(line, MAX_ZEILEN);
	line_string = line;
	//----------------------------------------------------------------------
	// evaluate header
	//----------------------------------------------------------------------
	while (!csv_file.eof())
	{
		csv_file.getline(line, MAX_ZEILEN);
		line_string = line;
		if (line_string.empty())
			return;
		sprintf(pnt_name, "PNT_PROPERTY_%i", (int)gli_points_vector.size());
		m_pnt = new CGLPoint();
		m_pnt->name = (string)pnt_name;
		gli_points_vector.push_back(m_pnt);
		m_pnt->id = (long)gli_points_vector.size() - 1;
		// lese zeile No	Name	X	Y	Z	m3/s ...
		// No
		pos1 = 0;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->name = sub_string;
		// X
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->x = strtod(sub_string.data(), NULL);
		// Y
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->y = strtod(sub_string.data(), NULL);
		// Z
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->z = strtod(sub_string.data(), NULL);
		// LENGTH
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->length = strtod(sub_string.data(), NULL);
		// ST
		pos1 = pos2 + 1;
		sub_string = get_sub_string(line_string, delimiter_type, pos1, &pos2);
		sub_string = line_string.substr(pos1, pos2 - pos1);
		m_pnt->value = strtod(sub_string.data(), NULL);
		m_pnt->value /= m_pnt->length;
		//
		pnt_properties_vector.push_back(m_pnt);
	} // eof
}
/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   07/2005 OK Implementation
**************************************************************************/
void GEOCalcPointMinMaxCoordinates()
{
	pnt_x_min = 1.e+19;
	pnt_x_max = -1.e+19;
	pnt_y_min = 1.e+19;
	pnt_y_max = -1.e+19;
	pnt_z_min = 1.e+19;
	pnt_z_max = -1.e+19;
	double value;
	//----------------------------------------------------------------------
	size_t size(gli_points_vector.size());
	for (size_t i = 0; i < size; i++)
	{
		value = gli_points_vector[i]->x;
		if (value < pnt_x_min)
			pnt_x_min = value;
		if (value > pnt_x_max)
			pnt_x_max = value;
		value = gli_points_vector[i]->y;
		if (value < pnt_y_min)
			pnt_y_min = value;
		if (value > pnt_y_max)
			pnt_y_max = value;
		value = gli_points_vector[i]->z;
		if (value < pnt_z_min)
			pnt_z_min = value;
		if (value > pnt_z_max)
			pnt_z_max = value;
	}
}

/**************************************************************************
   GeoLib - Funktion: IsInsidePolygonPlain
   Programmaenderungen:
   10/1998     AH         Erste Version
   08/2005 CC Modification
**************************************************************************/
bool CGLPoint::IsInsidePolygonPlain(double* xp, double* yp, double* zp, long np)
{
	long k;
	double dp[3], dn[3];
	double phi, cosphi, norm_dn, norm_dp;
	double sumphi = 0.0;
	double pi = 4.0 * atan(1.);
	double vorz = 0.;
	double nxp[3], pxp[3];
	double d_eps = 10. * MKleinsteZahlen;
	double a_eps = 1000. * MKleinsteZahlen;

	if (np <= 2)
		return true;

	dn[0] = xp[0] - x;
	dn[1] = yp[0] - y;
	dn[2] = zp[0] - z;
	dp[0] = xp[1] - x;
	dp[1] = yp[1] - y;
	dp[2] = zp[1] - z;
	M3KreuzProdukt(dn, dp, pxp);
	for (k = 1; k < np; k++)
	{
		dp[0] = xp[k] - x;
		dp[1] = yp[k] - y;
		dp[2] = zp[k] - z;
		norm_dn = MBtrgVec(dn, 3);
		norm_dp = MBtrgVec(dp, 3);
		if (norm_dn <= d_eps || norm_dp <= d_eps)
			return true; /* auf Rand */
		cosphi = MSkalarprodukt(dn, dp, 3) / norm_dn / norm_dp;
		phi = acos(cosphi);
		M3KreuzProdukt(dn, dp, nxp);
		if (MSkalarprodukt(nxp, pxp, 3) >= 0.0)
			vorz = 1.0;
		else
			vorz = -1.0;
		sumphi += vorz * phi;
		dn[0] = dp[0];
		dn[1] = dp[1];
		dn[2] = dp[2];
	}

	if (fabs(sumphi - 2. * pi) <= a_eps)
		return true; /* Im Gebiet */
	else if (fabs(sumphi - pi) <= a_eps)
		return true; /* auf Rand */
	else
		return false;
}
/**************************************************************************
   MSHLib-Method: NodeInTriangle
   Task:
   Programing:
   11/2003 OK Implementation
   08/2005 CC Modification CGLPoint*
**************************************************************************/
bool CGLPoint::IsInsideTriangle(double* xv, double* yv, double* zv)
{
	double vec1[3], vec2[3], vec3[3];
	double area, area1, area2, area3;
	bool ok = false;
	double eps;
	double dist;
	// triangle area
	vec1[0] = xv[1] - xv[0];
	vec1[1] = yv[1] - yv[0];
	vec1[2] = zv[1] - zv[0];
	vec2[0] = xv[0] - xv[2];
	vec2[1] = yv[0] - yv[2];
	vec2[2] = zv[0] - zv[2];
	M3KreuzProdukt(vec1, vec2, vec3);
	area = fabs(0.5 * MBtrgVec(vec3, 3));
	eps = area * 1e-2;
	// partial area1
	vec1[0] = xv[0] - x;
	vec1[1] = yv[0] - y;
	vec1[2] = zv[0] - z;
	vec2[0] = xv[1] - x;
	vec2[1] = yv[1] - y;
	vec2[2] = zv[1] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area1 = fabs(0.5 * MBtrgVec(vec3, 3));
	// partial area2
	vec1[0] = xv[1] - x;
	vec1[1] = yv[1] - y;
	vec1[2] = zv[1] - z;
	vec2[0] = xv[2] - x;
	vec2[1] = yv[2] - y;
	vec2[2] = zv[2] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area2 = fabs(0.5 * MBtrgVec(vec3, 3));
	// partial area
	vec1[0] = xv[2] - x;
	vec1[1] = yv[2] - y;
	vec1[2] = zv[2] - z;
	vec2[0] = xv[0] - x;
	vec2[1] = yv[0] - y;
	vec2[2] = zv[0] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area3 = fabs(0.5 * MBtrgVec(vec3, 3));
	// check
	dist = fabs(area - area1 - area2 - area3);
	if (dist < eps)
		ok = true;
	return ok;
}
/**************************************************************************
   MSHLib-Method: IsInsideRectangle
   Task:
   Programing:
   11/2003 OK Implementation
   08/2005 CC Modification CGLPoint*
**************************************************************************/
bool CGLPoint::IsInsideRectangle(double* xv, double* yv, double* zv)
{
	double vec1[3], vec2[3], vec3[3];
	double area, area1, area2, area3, area4;
	bool ok = false;
	double eps;
	double dist;
	// rectangle area
	vec1[0] = xv[1] - xv[0];
	vec1[1] = yv[1] - yv[0];
	vec1[2] = zv[1] - zv[0];
	vec2[0] = xv[0] - xv[2];
	vec2[1] = yv[0] - yv[2];
	vec2[2] = zv[0] - zv[2];
	M3KreuzProdukt(vec1, vec2, vec3);
	area = fabs(0.5 * MBtrgVec(vec3, 3));
	vec1[0] = xv[2] - xv[0];
	vec1[1] = yv[2] - yv[0];
	vec1[2] = zv[2] - zv[0];
	vec2[0] = xv[0] - xv[3];
	vec2[1] = yv[0] - yv[3];
	vec2[2] = zv[0] - zv[3];
	M3KreuzProdukt(vec1, vec2, vec3);
	area += fabs(0.5 * MBtrgVec(vec3, 3));
	eps = area * 1e-2;
	// partial area1
	vec1[0] = xv[0] - x;
	vec1[1] = yv[0] - y;
	vec1[2] = zv[0] - z;
	vec2[0] = xv[1] - x;
	vec2[1] = yv[1] - y;
	vec2[2] = zv[1] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area1 = fabs(0.5 * MBtrgVec(vec3, 3));
	// partial area2
	vec1[0] = xv[1] - x;
	vec1[1] = yv[1] - y;
	vec1[2] = zv[1] - z;
	vec2[0] = xv[2] - x;
	vec2[1] = yv[2] - y;
	vec2[2] = zv[2] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area2 = fabs(0.5 * MBtrgVec(vec3, 3));
	// partial area3
	vec1[0] = xv[2] - x;
	vec1[1] = yv[2] - y;
	vec1[2] = zv[2] - z;
	vec2[0] = xv[3] - x;
	vec2[1] = yv[3] - y;
	vec2[2] = zv[3] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area3 = fabs(0.5 * MBtrgVec(vec3, 3));
	// partial area4
	vec1[0] = xv[3] - x;
	vec1[1] = yv[3] - y;
	vec1[2] = zv[3] - z;
	vec2[0] = xv[0] - x;
	vec2[1] = yv[0] - y;
	vec2[2] = zv[0] - z;
	M3KreuzProdukt(vec1, vec2, vec3);
	area4 = fabs(0.5 * MBtrgVec(vec3, 3));
	// check
	dist = fabs(area - area1 - area2 - area3 - area4);
	if (dist < eps)
		ok = true;
	return ok;
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   11/2003 OK Implementation
   08/2005 CC Modification CGLPoint* - Move from GeoLib to MSHLib
**************************************************************************/
bool CGLPoint::IsInsidePrism(double* xv, double* yv, double* zv)
{
	bool ok = false;
	double t_vol1, t_vol2;
	double p_vol1, p_vol2, p_vol3;
	double sum, prism;

	double xt[4], yt[4], zt[4];
	xt[0] = xv[0];
	xt[1] = xv[1];
	xt[2] = xv[2];
	xt[3] = x;
	yt[0] = yv[0];
	yt[1] = yv[1];
	yt[2] = yv[2];
	yt[3] = y;
	zt[0] = zv[0];
	zt[1] = zv[1];
	zt[2] = zv[2];
	zt[3] = z;
	t_vol1 = CalcTetraederVolume(xt, yt, zt);
	xt[0] = xv[3];
	xt[1] = xv[4];
	xt[2] = xv[5];
	xt[3] = x;
	yt[0] = yv[3];
	yt[1] = yv[4];
	yt[2] = yv[5];
	yt[3] = y;
	zt[0] = zv[3];
	zt[1] = zv[4];
	zt[2] = zv[5];
	zt[3] = z;
	t_vol2 = CalcTetraederVolume(xt, yt, zt);

	double xp[5], yp[5], zp[5];
	xp[0] = xv[0];
	xp[1] = xv[1];
	xp[2] = xv[4];
	xp[3] = xv[3];
	xp[4] = x;
	yp[0] = yv[0];
	yp[1] = yv[1];
	yp[2] = yv[4];
	yp[3] = yv[3];
	yp[4] = y;
	zp[0] = zv[0];
	zp[1] = zv[1];
	zp[2] = zv[4];
	zp[3] = zv[3];
	zp[4] = z;
	p_vol1 = CalcPyramidVolume(xp, yp, zp);
	xp[0] = xv[1];
	xp[1] = xv[2];
	xp[2] = xv[5];
	xp[3] = xv[4];
	xp[4] = x;
	yp[0] = yv[1];
	yp[1] = yv[2];
	yp[2] = yv[5];
	yp[3] = yv[4];
	yp[4] = y;
	zp[0] = zv[1];
	zp[1] = zv[2];
	zp[2] = zv[5];
	zp[3] = zv[4];
	zp[4] = z;
	p_vol2 = CalcPyramidVolume(xp, yp, zp);
	xp[0] = xv[2];
	xp[1] = xv[0];
	xp[2] = xv[3];
	xp[3] = xv[5];
	xp[4] = x;
	yp[0] = yv[2];
	yp[1] = yv[0];
	yp[2] = yv[3];
	yp[3] = yv[5];
	yp[4] = y;
	zp[0] = zv[2];
	zp[1] = zv[0];
	zp[2] = zv[3];
	zp[3] = zv[5];
	zp[4] = z;
	p_vol3 = CalcPyramidVolume(xp, yp, zp);

	sum = t_vol1 + t_vol2 + p_vol1 + p_vol2 + p_vol3;
	prism = CalcPrismVolume(xv, yv, zv);
	double eps = pow(prism, 0.33) * 1e-3;
	if (fabs(sum - prism) < eps)
		ok = true;

	return ok;
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   11/2003 OK Implementation
   08/2005 CC Modification Move
**************************************************************************/
bool CGLPoint::IsInTriangleXYProjection(double* xv, double* yv)
{
	double area, area1, area2, area3;
	bool ok = false;
	double eps;
	double dist;
	double mat3x3[9];
	// triangle area
	mat3x3[0] = 1.0;
	mat3x3[1] = 1.0;
	mat3x3[2] = 1.0;
	mat3x3[3] = xv[0];
	mat3x3[4] = xv[1];
	mat3x3[5] = xv[2];
	mat3x3[6] = yv[0];
	mat3x3[7] = yv[1];
	mat3x3[8] = yv[2];
	area = fabs(0.5 * M3Determinante(mat3x3));
	eps = area * 1e-6;
	// partial area1
	mat3x3[0] = 1.0;
	mat3x3[1] = 1.0;
	mat3x3[2] = 1.0;
	mat3x3[3] = x;
	mat3x3[4] = xv[1];
	mat3x3[5] = xv[2];
	mat3x3[6] = y;
	mat3x3[7] = yv[1];
	mat3x3[8] = yv[2];
	area1 = fabs(0.5 * M3Determinante(mat3x3));
	// partial area2
	mat3x3[0] = 1.0;
	mat3x3[1] = 1.0;
	mat3x3[2] = 1.0;
	mat3x3[3] = xv[0];
	mat3x3[4] = x;
	mat3x3[5] = xv[2];
	mat3x3[6] = yv[0];
	mat3x3[7] = y;
	mat3x3[8] = yv[2];
	area2 = fabs(0.5 * M3Determinante(mat3x3));
	// partial area3
	mat3x3[0] = 1.0;
	mat3x3[1] = 1.0;
	mat3x3[2] = 1.0;
	mat3x3[3] = xv[0];
	mat3x3[4] = xv[1];
	mat3x3[5] = x;
	mat3x3[6] = yv[0];
	mat3x3[7] = yv[1];
	mat3x3[8] = y;
	area3 = fabs(0.5 * M3Determinante(mat3x3));
	// check
	dist = fabs(area - area1 - area2 - area3);
	if (dist < eps)
		ok = true;
	return ok;
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
double CGLPoint::PointDisXY(CGLPoint* m_pnt)
{
	return sqrt((x - m_pnt->x) * (x - m_pnt->x) + (y - m_pnt->y) * (y - m_pnt->y));
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
CGLPoint* CGLPoint::Exist()
{
	double tolerance = 1e-3;
	return isPointInPointVector(gli_points_vector, this, tolerance);
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
int GEOMaxPointID()
{
	int max_pnt_id = -1;
	CGLPoint* m_pnt = NULL;
	for (int i = 0; i < (long)gli_points_vector.size(); i++)
	{
		m_pnt = gli_points_vector[i];
		if (m_pnt->id > max_pnt_id)
			max_pnt_id = m_pnt->id;
	}
	return max_pnt_id;
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   11/2005 CC Implementation
**************************************************************************/
void GEOWritePoints(char* file_name)
{
	FILE* f = NULL;
	string gli_file_name;
	long gli_points_vector_size;

	long i;
	const char* gli_file_name_char = 0;
	// File handling
	gli_file_name = (string)file_name + ".gli";
	gli_file_name_char = gli_file_name.data();
	f = fopen(gli_file_name_char, "w+");
	if (f)
	{
		fprintf(f, "#POINTS\n");
		gli_points_vector_size = (long)gli_points_vector.size();
		for (i = 0; i < gli_points_vector_size; i++)
			if (gli_points_vector[i]->id > -1) // OK
				fprintf(f, " %ld %#16.11g %#16.11g %#16.11g $MD %g $ID %s\n", gli_points_vector[i]->id,
				        gli_points_vector[i]->x, gli_points_vector[i]->y, gli_points_vector[i]->z,
				        gli_points_vector[i]->mesh_density, gli_points_vector[i]->name.c_str());
		fclose(f);
	}
	else
	{
	}
}
/**************************************************************************
   GeoLib-Method:
   Task: Calculates angle sum in a tringle.
   Only Point Input => point method.
   Programing:
   09/2005 TK Implementation
**************************************************************************/
double AngleSumPointInsideTriangle(double* point, double* tri_p1, double* tri_p2, double* tri_p3, double tolerance)
{
	double sum = 0.0;
	// OK char sum_check[56];
	double pi = 3.14159265359;
	double a_quantum, b_quantum, c_quantum;
	double cos_alpha[3], alpha_rad[3], alpha_deg[3];
	double a[3], b[3], c[3];

	/* a-Vektor-Coordinates (a = P0P1)*/
	a[0] = tri_p1[0] - point[0];
	a[1] = tri_p1[1] - point[1];
	a[2] = tri_p1[2] - point[2];
	a_quantum = sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));
	/* b-Vektor-Coordinates (a = P0P2)*/
	b[0] = tri_p2[0] - point[0];
	b[1] = tri_p2[1] - point[1];
	b[2] = tri_p2[2] - point[2];
	b_quantum = sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]));
	/* c-Vektor-Coordinates (a = P0P3)*/
	c[0] = tri_p3[0] - point[0];
	c[1] = tri_p3[1] - point[1];
	c[2] = tri_p3[2] - point[2];
	c_quantum = sqrt((c[0] * c[0]) + (c[1] * c[1]) + (c[2] * c[2]));

	/*Angle: at P0 -P1P2*/
	cos_alpha[0] = ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])) / (a_quantum * b_quantum);
	if (cos_alpha[0] > 1.0)
		cos_alpha[0] = 1.0;
	if (cos_alpha[0] < -1.0)
		cos_alpha[0] = -1.0;
	alpha_rad[0] = acos(cos_alpha[0]);
	alpha_deg[0] = (180 / pi) * acos(cos_alpha[0]);
	/*Angle: at P0 -P2P3*/
	cos_alpha[1] = ((b[0] * c[0]) + (b[1] * c[1]) + (b[2] * c[2])) / (b_quantum * c_quantum);
	if (cos_alpha[1] > 1.0)
		cos_alpha[1] = 1.0;
	if (cos_alpha[1] < -1.0)
		cos_alpha[1] = -1.0;
	alpha_rad[1] = acos(cos_alpha[1]);
	alpha_deg[1] = (180 / pi) * acos(cos_alpha[1]);
	/*Angle: at P0 -P3P1*/
	cos_alpha[2] = ((c[0] * a[0]) + (c[1] * a[1]) + (c[2] * a[2])) / (c_quantum * a_quantum);
	if (cos_alpha[2] > 1.0)
		cos_alpha[2] = 1.0;
	if (cos_alpha[2] < -1.0)
		cos_alpha[2] = -1.0;
	alpha_rad[2] = acos(cos_alpha[2]);
	alpha_deg[2] = (180 / pi) * acos(cos_alpha[2]);

	sum = (cos_alpha[0] + cos_alpha[1] + cos_alpha[2]);
	sum = (alpha_rad[0] + alpha_rad[1] + alpha_rad[2]);
	sum = (alpha_deg[0] + alpha_deg[1] + alpha_deg[2]);

	if (cos_alpha[0] == 1.0 || cos_alpha[1] == 1.0 || cos_alpha[2] == 1.0 || cos_alpha[0] == -1.0
	    || cos_alpha[0] == -1.0
	    || cos_alpha[0] == -1.0)
	{
		if (cos_alpha[0] == 1.0 || cos_alpha[1] == 1.0 || cos_alpha[2] == 1.0)
			sum = 180;
		else
			sum = 360;
	}

	if ((point[0] == tri_p1[0] && point[1] == tri_p1[1] && point[2] == tri_p1[2])
	    || (point[0] == tri_p2[0] && point[1] == tri_p2[1] && point[2] == tri_p2[2])
	    || (point[0] == tri_p3[0] && point[1] == tri_p3[1] && point[2] == tri_p3[2]))
		sum = 360;

	if (a_quantum < tolerance || b_quantum < tolerance || c_quantum < tolerance)
		sum = 360;

	// if (alpha_deg[0]<1 || alpha_deg[1]<1 || alpha_deg[2]<1)sum=360;
	return sum; // 359.99999999997630
}

double sqrDist(const CGLPoint* p0, const CGLPoint* p1)
{
	double v[3] = {p1->getX() - p0->getX(), p1->getY() - p0->getY(), p1->getZ() - p0->getZ()};
	return MathLib::scpr(v, v, 3);
}

CGLPoint* isPointInPointVector(const std::vector<CGLPoint*>& vec, const CGLPoint* const p, double tol)
{
	size_t size(vec.size());
	for (size_t i = 0; i < size; ++i)
		if (sqrDist(p, vec[i]) < tol)
			return vec[i];
	return NULL;
}

// LB TODO: trunk transition hack
string get_sub_string(string buffer, string delimiter, int pos1, int* pos2)
{
	int pos = 0;
	string empty_string("");
	// string sub_string_this;
	*pos2 = (int)buffer.find(delimiter, pos1);
	if (*pos2 < 0)
		return empty_string;
	while (*pos2 <= pos1)
	{
		pos1++;
		*pos2 = (int)buffer.find(delimiter, pos1);
		if (*pos2 < 0)
		{
			*pos2 = (int)buffer.size();
			break;
		}
		if (pos1 >= (int)buffer.size())
			break;
	}
	string sub_string_this = buffer.substr(pos1, *pos2);
	while (pos >= 0)
	{
		pos = (int)sub_string_this.find_first_of(" ");
		if (pos < 0)
			break;
		sub_string_this.erase(pos, 1);
	}
	return sub_string_this;
}
