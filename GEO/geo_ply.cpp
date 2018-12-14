/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <stdlib.h>
/* Objects */
#include "files0.h"
#include "geo_lib.h"
#include "geo_mathlib.h"
#include "geo_ply.h"
#include "geo_pnt.h"
#include "geo_sfc.h"

// MathLib
#include "InterpolationAlgorithms/CubicSpline.h"

using namespace std;

// GSP
std::vector<CGLPolyline*> polyline_vector;
int ply_max = -1;

std::vector<CColumn*> column_vector;
std::vector<CSoilProfile*> profile_vector;  // YD

/*----------------------------------------------------------------------*/
// constructor
CGLPolyline::CGLPolyline(void) : _set_eps(false), epsilon(0.01)
{
    name = "POLYLINE";
    //	closed = false;
    type = 2;
    computeline = false;
    mat_group = -1;
    data_type = 0;
    //	minDis = 0;
    ply_max++;  // OK
    id = ply_max;
    //	mesh_density = 100.0;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
CGLPolyline::CGLPolyline(const std::string& ply_name)
    : name(ply_name), _set_eps(false), epsilon(0.01)
{
    //	closed = false;
    type = 2;
    computeline = false;
    mat_group = -1;
    data_type = 0;
    //	minDis = 0;
    ply_max++;  // OK
    id = ply_max;
}

// deconstructor
CGLPolyline::~CGLPolyline(void)
{
    sbuffer.clear();
    ibuffer.clear();
    OrderedPoint.clear();
    msh_nodes_vector.clear();
}

const std::string& CGLPolyline::getName() const
{
    return name;
}

void CGLPolyline::setName(const std::string& nname)
{
    name = nname;
}

size_t CGLPolyline::getID() const
{
    return id;
}

void CGLPolyline::setID(size_t nid)
{
    id = nid;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/
CGLPolyline* GEOGetPLYByName(const std::string& name)
{
    for (std::vector<CGLPolyline*>::iterator it = polyline_vector.begin();
         it != polyline_vector.end();
         ++it)
    {
        if ((*it)->getName() == name)
            return *it;
    }
    return NULL;
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/
CGLPolyline* GEOGetPLYById(size_t number)
{
    std::vector<CGLPolyline*>::iterator p = polyline_vector.begin();  // CC
    while (p != polyline_vector.end())
    {
        if ((*p)->getID() == number)
            return *p;
        ++p;
    }
    return NULL;
}
/**************************************************************************
   GeoLib-Method: GEORemovePolylines
   Task:
   Programing:
   11/2003 CC Implementation
   01/2005 OK List destructor
   CCToDo Polyline destructor
   08/2005 CC Modification
   02/2006 CC destructor
   05/2006 TK BUGFIX
**************************************************************************/
void GEORemoveAllPolylines()
{
    int i;
    // CGLPolyline * m_ply = NULL;
    for (i = 0; i < (int)polyline_vector.size(); i++)
    {
        // m_ply = polyline_vector[0]; //TK: What's that Cui?
        // delete m_ply;
        delete polyline_vector[i];
        polyline_vector[i] = NULL;
    }
    polyline_vector.clear();  // CC

    for (i = 0; i < (int)gli_lines_vector.size(); i++)
    {
        delete gli_lines_vector[i];
        gli_lines_vector[i] = NULL;
    }
    gli_lines_vector.clear();
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2006 CC destructor
**************************************************************************/
void GEORemovePolyline(long nSel)
{
    CGLPolyline* m_ply = NULL;
    m_ply = polyline_vector[nSel];
    delete m_ply;
    polyline_vector.erase(polyline_vector.begin() + nSel);
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   07/2003 OK Implementation
   08/2005 CC Modification
**************************************************************************/
std::vector<CGLPolyline*> GetPolylineVector(void)
{
    return polyline_vector;
}
/**************************************************************************
   GeoLib-Method: GEOPolylineGLI2GEO
   Task: convert PLI to GEO data
   Programing:
   11/2002 WW/OK   Implementation
   12/2002 OK Meshing algorithm output
   12/2002 OK use specific polylines for meshing (type=1)
   01/2002 WW Extend to multi-polyline
   09/2003 OK for GeoLib
   12/2003 OK construct lines only once
   02/2004 WW Surface-wise output
   11/2005 TK FOR-LOOP because of GEOLIB2
**************************************************************************/
void GEOPolylineGLI2GEO(FILE* geo_file)
{
    int i;
    // vector<CGLLine*> gli_lines_vector;
    std::vector<CGLPoint*> gli_Points_vector;
    long gli_lines_vector_size;
    gli_Points_vector = GetPointsVector();
    /*---------------------------------------------------------------*/
    /*--------------------   Write geo file  ------------------------*/
    /*--------- Mesh density scaling, smoothing method --------------*/
    /*---------------------------------------------------------------*/
    /*---------------------------------------------------------------*/
    /* Meshing algorithm */
    FilePrintString(geo_file, "Mesh.Algorithm = 3;");
    LineFeed(geo_file);
    FilePrintString(geo_file, "Mesh.Smoothing = 4;");
    LineFeed(geo_file);
    //----------------------------------------------------------------------
    // Compute lines
    gli_lines_vector_size = (long)gli_lines_vector.size();
    if (gli_lines_vector_size == 0)
    {
        CGLPolyline* p_pline = NULL;
        std::vector<CGLPolyline*>::iterator p = polyline_vector.begin();  // CC
        while (p != polyline_vector.end())
        {
            p_pline = *p;
            if (p_pline->getLineVector().size() == 0)  // CC
                p_pline->ComputeLines();
            p++;
        }
    }
    //----------------------------------------------------------------------
    // Write
    int p_counter = 0;
    int l_counter = 0;
    int s_counter = 0;
    int pl_counter = 0;
    for (i = 0; i < (int)surface_vector.size(); i++)
    {
        if (!surface_vector[i]->meshing_allowed)  // OK/TK
            continue;
        s_counter++;
        surface_vector[i]->output(
            geo_file, p_counter, l_counter, pl_counter, s_counter);
    }
    //----------------------------------------------------------------------
    // Null the index
    long gli_point_vector_size = (long)gli_Points_vector.size();
    for (i = 0; i < gli_point_vector_size; i++)
        gli_Points_vector[i]->SetIndex(-1);
    for (i = 0; i < gli_lines_vector_size; i++)
        gli_lines_vector[i]->mesh_index = -1;
    //----------------------------------------------------------------------
}
/**************************************************************************
   GeoLib-Method: GEOReadPolylineNew
   Task: Read polyline data from file
   Programing:
   07/2003 OK Implementation
   12/2003 no line output
   06(2005 OK point_vector
   11/2005 CC Write function
**************************************************************************/
void CGLPolyline::Write(char* file_name)
{
    FILE* f = NULL;
    const char* filename = 0;
    std::string gli_file_name;
    long i;
    // sprintf(gli_file_name,"%s.%s",file_name,"gli");
    gli_file_name = (string)file_name + ".gli";
    filename = gli_file_name.data();
    f = fopen(filename, "a");
    fprintf(f, "#POLYLINE\n");
    fprintf(f, " $ID\n");                          // CC
    fprintf(f, "  %ld\n", static_cast<long>(id));  // CC
    fprintf(f, " $NAME\n");
    if (data_type == 1)              // CC8888
        fprintf(f, "  POLYLINE\n");  // CC8888
    else
        // CC8888
        fprintf(f, "  %s\n", name.c_str());
    fprintf(f, " $TYPE\n");
    fprintf(f, "  %d\n", type);
    fprintf(f, " $EPSILON\n");
    fprintf(f, "  %g\n", epsilon);
    fprintf(f, " $MAT_GROUP\n");
    fprintf(f, "  %d\n", mat_group);
    if (data_type == 0)
    {
        fprintf(f, " $POINTS\n");
        for (i = 0; i < (long)point_vector.size(); i++)
            fprintf(f, " %ld\n", point_vector[i]->id);
    }
    else if (data_type == 1)
    {
        fprintf(f, " $POINT_VECTOR\n");
        std::string ply_file_name("POLYLINE");
        ply_file_name += PLY_FILE_EXTENSION;

        fprintf(f, "  %s\n", ply_file_name.data());  // TK
    }
    fclose(f);
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   11/2003 OK/CC Implementation
**************************************************************************/
bool CGLPolyline::PointExists(CGLPoint* point, CGLPoint* point1)
{
    long i;
    double dist;
    long dist2;
    long ply_points_vector_length = (long)point_vector.size();
    for (i = 0; i < ply_points_vector_length; i++)
    {
        dist2 = ((point_vector[i]->x_pix - point->x_pix) *
                     (point_vector[i]->x_pix - point->x_pix) +
                 (point_vector[i]->y_pix - point->y_pix) *
                     (point_vector[i]->y_pix - point->y_pix));
        dist = sqrt((double)dist2);
        if (dist < point->circle_pix)
        {
            point1->x_pix = point_vector[i]->x_pix;
            point1->y_pix = point_vector[i]->y_pix;
            point1->x = point_vector[i]->x;
            point1->y = point_vector[i]->y;
            point1->z = point_vector[i]->z;
            return true;
        }
    }
    return false;
}
/**************************************************************************
   GeoLib-Method: CenterPoint
   Task:
   Programing:
   12/2003 OK Implementation
**************************************************************************/
CGLPoint* CGLPolyline::CenterPoint(void)
{
    CGLPoint* m_point = NULL;
    int polyline_point_vector_length = (long)point_vector.size();
    if (polyline_point_vector_length == 0)
        return NULL;
    int center_point = (int)polyline_point_vector_length / 2;
    m_point = point_vector[center_point];
    return m_point;
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   01/2005 OK File handling
   08/2005 CC
**************************************************************************/
void CGLPolyline::WritePointVector(const std::string& base)
{
    (void)base;

    string ply_path;
    string ply_path_base_type;
    long i;
    //----------------------------------------------------------------------
    long no_points = (long)point_vector.size();
    if (no_points > 0)
    {
        //----------------------------------------------------------------------
        // File handling
        std::string ply_file_name = "POLYLINE";
        ply_file_name += PLY_FILE_EXTENSION;  // CC
        ply_path_base_type = ply_path + ply_file_name;
        std::fstream ply_file(ply_path_base_type.data(),
                              std::ios::trunc | std::ios::out);
        ply_file.setf(std::ios::scientific, std::ios::floatfield);
        ply_file.precision(12);
        //--------------------------------------------------------------------
        if (!ply_file.good())
            return;
        ply_file.seekg(0L, std::ios::beg);
        for (i = 0; i < no_points; i++)
            ply_file << point_vector[i]->x << " " << point_vector[i]->y << " "
                     << point_vector[i]->z << "\n";
    }
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   03/2004 OK Implementation
   08/2005 CC file_path
   10/2005 OK Path
**************************************************************************/
void CGLPolyline::ReadPointVector(const std::string& base)
{
    std::string cut_string;
    std::string delimiter_type(" ");
    char line[MAX_ZEILEN];
    std::string line_string;
    std::string sub_string;
    std::string buffer;
    CGLPoint* m_point = NULL;
    //----------------------------------------------------------------------
    // File handling
    std::string ply_file_name = base;
    std::ifstream ply_file(base.c_str());  //,std::ios::in);
    if (!ply_file.is_open())
    {
        std::cout << "*** Warning in CGLPolyline::ReadPointVector: File "
                  << ply_file_name << " not found"
                  << "\n";
        return;
    }
    ply_file.seekg(0L, std::ios::beg);
    //----------------------------------------------------------------------
    while (!ply_file.eof())
    {
        ply_file.getline(line, MAX_ZEILEN);
        buffer = line;
        if (buffer.size() == 0)
            continue;
        m_point = new CGLPoint;
        //
        sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
        m_point->x = strtod(sub_string.data(), NULL);
        buffer = cut_string;
        //
        sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
        m_point->y = strtod(sub_string.data(), NULL);
        buffer = cut_string;
        //
        sub_string = get_sub_string2(buffer, delimiter_type, &cut_string);
        m_point->z = strtod(sub_string.data(), NULL);
        buffer = cut_string;
        //
        point_vector.push_back(m_point);
    }
}

/**************************************************************************
   GeoLib-Method: ComputeLines
   Task:
   Programing:
   11/2003 WW Implementation
   01/2004 CC Modification -- remove orientation value
   01/2005 CC Modification -- remove polyline->closed part
**************************************************************************/
void CGLPolyline::ComputeLines()
{
    // 1 creating all lines from this polyline
    // 2 put lines into polyline line
    if (!computeline)
    {
        CGLPoint* m_point1 = NULL;
        CGLPoint* m_point2 = NULL;
        CGLLine* m_line_exist = NULL;
        CGLLine* m_line2 = NULL;

        size_t number_of_polyline_points = point_vector.size();
        for (size_t i = 1; i < number_of_polyline_points; i++)
        {
            m_point1 = point_vector[i - 1];
            m_point2 = point_vector[i];
            CGLLine* m_line = new CGLLine;
            m_line->point1 = m_point1->id;
            m_line->point2 = m_point2->id;
            m_line->m_point1 = m_point1;
            m_line->m_point2 = m_point2;

            m_line_exist = m_line->Exists();
            if (m_line_exist)
                line_vector.push_back(m_line_exist);  // CC
            else
            {
                m_line2 = m_line;
                gli_lines_vector.push_back(m_line2);
                // m_line2->line_index =(long)gli_lines_vector.size();
                m_line2->gli_line_id = (long)gli_lines_vector.size();
                m_line2->orientation = 1;
                line_vector.push_back(m_line2);  // CC
            }
        }
        computeline = true;
    }
}

/**************************************************************************
   GEOLib-Method:
   Task: polyline read function
   Programing:
   03/2004 CC Implementation
   05/2004 CC Modification
   04/2005 CC Modification calculate the minimal distance between points
reference for mesh density of line element calculation 07/2005 CC read ID of
polyline 08/2005 CC parameter
**************************************************************************/
void GEOReadPolylines(const std::string& file_name_path_base)
{
    CGLPolyline* m_polyline = NULL;
    char line[MAX_ZEILEN];
    std::string sub_line;
    std::string line_string;
    std::string gli_file_name;
    std::string path_name;
    std::ios::pos_type position;
    //========================================================================
    // File handling
    gli_file_name = file_name_path_base + ".gli";
    //  pos = (int)file_name_path_base.rfind('\\'); //CC remove
    // path_name = file_name_path_base.substr(0,(pos+1));//CC remove
    std::ifstream gli_file(gli_file_name.data(), std::ios::in);
    if (!gli_file.good())
    {
        std::cerr << "stream error GEOReadPolylines "
                  << "\n";
        return;
    }
    gli_file.seekg(0L, std::ios::beg);  // rewind?
    //========================================================================
    // Keyword loop
    while (!gli_file.eof())
    {
        gli_file.getline(line, MAX_ZEILEN);
        line_string = line;

        if (line_string.find("#STOP") != string::npos)  // 11.08.2011. WW
            break;

        //----------------------------------------------------------------------
        if (line_string.find("#POLYLINE") != string::npos)  // keyword found
        {
            m_polyline = new CGLPolyline();
            position = m_polyline->Read(gli_file);  // CC8888, TF
            polyline_vector.push_back(m_polyline);
            gli_file.seekg(position, std::ios::beg);
            // OK->CC encapsulate function
            //..................................................................
            //			m_polyline->CalcMinimumPointDistance();
            //..................................................................
        }              // keyword found
    }                  // eof
    gli_file.close();  // OK41
}

/**************************************************************************
   GeoLib-Method: Read
   Task: Read polyline data from file
   Programing:
   03/2004 CC Implementation
   09/2004 OK file path for PLY files
   07/2005 CC PLY id
   08/2005 CC parameter
   09/2005 CC itoa - convert integer to string
**************************************************************************/
ios::pos_type CGLPolyline::Read(std::ifstream& gli_file)  // CC8888
{
    char line[MAX_ZEILEN];
    std::string line_string;

    bool new_keyword = false;
    std::ios::pos_type position;
    CGLPoint* m_point = NULL;

    // Schleife ueber alle Phasen bzw. Komponenten
    while (!new_keyword)
    {
        position = gli_file.tellg();
        gli_file.getline(line, MAX_ZEILEN);
        line_string = line;
        if (line_string.find("#") != string::npos)
        {
            new_keyword = true;
            break;
        }
        if (line_string.find("$ID") != string::npos)  // subkeyword found CC
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            id = strtol(line_string.data(), NULL, 0);

            continue;
        }
        //....................................................................
        if (line_string.find("$NAME") != string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            name = line_string.substr(0);
            continue;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$TYPE") != string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            type = strtol(line_string.data(), NULL, 0);
            continue;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$EPSILON") != string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            epsilon = strtod(line_string.data(), NULL);
            _set_eps = true;
            continue;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$MAT_GROUP") != string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            mat_group = strtol(line_string.data(), NULL, 0);
            continue;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$POINTS") != string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            //			ply_data = "POINTS";//CC8888
            long lvalue = strtol(line_string.data(), NULL, 0);
            int i = 1;
            while (i == 1)
            {
                m_point = GEOGetPointById(lvalue);  // CC wrong get point by id
                if (m_point)
                    AddPoint(m_point);
                // CC---------------------------
                else
                {
                    std::cout << "Error: point " << lvalue << " not found"
                              << "\n";
                    //--------------------------------------------------
                    std::string m_strname = name + ": Point not found: point ";
                    char m_strnameid[10];
                    sprintf(m_strnameid, "%li", lvalue);  // OK
                    // itoa((int)lvalue,m_strnameid,10); //CC 09/05 convert
                    // integer to string ultoa convert unsigned long to string
                    std::string error_str = m_strname + m_strnameid;
                    return position;
                }
                gli_file.getline(line, MAX_ZEILEN);
                line_string = line;
                if ((line_string.find("#") != string::npos)

                    || (line_string.find("$") != string::npos))  // OK bugfix
                {
                    i = 0;
                    new_keyword = true;
                }
                else
                    lvalue = strtol(line_string.data(), NULL, 0);
            }  // end of while
            data_type = 0;
        }  // subkeyword found
        //....................................................................
        if (line_string.find("$POINT_VECTOR") !=
            string::npos)  // subkeyword found
        {
            gli_file.getline(line, MAX_ZEILEN);
            line_string = line;
            remove_white_space(&line_string);
            ReadPointVector(line_string.substr(0));  // CC
            data_type = 1;                           // OK41
            //			ply_data = line_string;
        }  // subkeyword found
        //========================================================================
    }
    return position;
}

/**************************************************************************
   GeoLib-Method:addpoint
   Task:
   Programing:
   03/2004 CC Implementation
**************************************************************************/
void CGLPolyline::AddPoint(CGLPoint* m_point)
{
    // point_list.push_back(m_point);//CC remove
    point_vector.push_back(m_point);
}
/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   05/2004 CC Implementation Assign randam color to polyline
   last modification:
**************************************************************************/
// void CGLPolyline::AssignColor() {
//	long r_number = rand();
//	m_color[0] = (long) (r_number / (double) RAND_MAX * 256);
//	r_number = rand();
//	m_color[1] = (long) (r_number / (double) RAND_MAX * 256);
//	r_number = rand();
//	m_color[2] = (long) (r_number / (double) RAND_MAX * 256);
//}

/**************************************************************************
   FEMLib-Method: InterpolationAlongPolyline
   Task: Prescribe the boundary values to all points of a polyline by the means
   of spline
   Programing:
   02/2004 WW Implementation
   last modification:
**************************************************************************/
void InterpolationAlongPolyline(CGLPolyline* plyL,
                                std::vector<double>& bcNodalValue)
{
    // Obtain fem node for groupvector
    const size_t SizeCGLPoint(plyL->point_vector.size());
    double xp(0.0), yp(0.0), zp(0.0);

    // Prepare spline data
    double sl = 0.0;
    std::vector<double> ss0;
    std::vector<double> bVal;
    for (size_t i = 0; i < SizeCGLPoint - 1; i++)
    {
        CGLPoint* CGPa = plyL->point_vector[i];
        CGLPoint* CGPb = plyL->point_vector[i + 1];
        xp = CGPb->x - CGPa->x;
        yp = CGPb->y - CGPa->y;
        zp = CGPb->z - CGPa->z;
        // NW: why checking zero value? Zero could also be used for BC/ST.
        if (fabs(plyL->point_vector[i]->getPropert()) > MKleinsteZahlen)
        {
            bVal.push_back(plyL->point_vector[i]->getPropert());
            ss0.push_back(sl);
        }
        sl += sqrt(xp * xp + yp * yp + zp * zp);
    }
    // Last point
    if (fabs(plyL->point_vector[SizeCGLPoint - 1]->getPropert()) >
        MKleinsteZahlen)
    {
        ss0.push_back(sl);
        bVal.push_back(plyL->point_vector[SizeCGLPoint - 1]->getPropert());
    }

    if (ss0.size() == 0)
    {
        std::cout << "Error in CGLPolyline::InterpolationAlongPolyline: no PNT "
                     "data found for Spline interpolation"
                  << "\n";
        return;
    }

    // Spline interpolation
    std::cout << "MathLib::CubicSpline input data: "
              << "\n";
    for (size_t k(0); k < ss0.size(); k++)
        std::cout << "\t" << ss0[k] << " " << bVal[k] << "\n";
    std::cout << "end MathLib::CubicSpline input data"
              << "\n";

    MathLib::CubicSpline* csp = new MathLib::CubicSpline(ss0, bVal);

    // Interpolate
    size_t number_of_nodes(bcNodalValue.size()), i;
    for (i = 0; i < number_of_nodes; i++)
        bcNodalValue[plyL->getOrderedPoints()[i]] =
            csp->interpolation(plyL->getSBuffer()[i]);
    std::cout << "MathLib::CubicSpline results: "
              << "\n";
    for (size_t k(0); k < plyL->getSBuffer().size(); k++)
        std::cout << "\t" << plyL->getOrderedPoints()[k] << " "
                  << plyL->getSBuffer()[k] << " "
                  << bcNodalValue[plyL->getOrderedPoints()[k]]
                  << " (bcNodalValue[" << plyL->getOrderedPoints()[k] << "])"
                  << "\n";
    std::cout << "end MathLib::CubicSpline results"
              << "\n";

    // Release the memory
    delete csp;
    csp = NULL;
    ss0.clear();
    bVal.clear();
    //  OrderedNode.clear();
    //  plyL->sbuffer.clear();
    //  plyL->ibuffer.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   08/2005 CC file_path
   last modification:
**************************************************************************/
void CGLPolyline::WriteTecplot(const std::string& file_path)
{
    long i;
    //----------------------------------------------------------------------
    // File handling
    std::string tec_path;
    /* CGSProject* m_gsp = GSPGetMember("gli");
       if(m_gsp)
       tec_path = m_gsp->path; */

    std::string tec_file_name = file_path + name + ".tec";
    // string tec_file_name = tec_path + name + ".tec";
    std::fstream tec_file(tec_file_name.data(), ios::trunc | ios::out);
    tec_file.setf(std::ios::scientific, std::ios::floatfield);
    tec_file.precision(12);
    // Write header
    tec_file << "VARIABLES = X,Y,Z"
             << "\n";
    long no_nodes = (long)point_vector.size();
    long no_elements = no_nodes - 1;
    tec_file << "ZONE T = " << name << ", "
             << "N = " << no_nodes << ", "
             << "E = " << no_elements << ", "
             << "F = FEPOINT"
             << ", "
             << "ET = TRIANGLE"
             << "\n";
    // Write data
    for (i = 0; i < no_nodes; i++)
        tec_file << point_vector[i]->x << " " << point_vector[i]->y << " "
                 << point_vector[i]->z << " "
                 << "\n";
    for (i = 0; i < no_elements; i++)
        tec_file << i + 1 << " " << i + 1 + 1 << " " << i + 1 << "\n";
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
**************************************************************************/
void CGLPolyline::SortPointVectorByDistance()
{
    long no_points = (long)point_vector.size();
    if (no_points == 0)
        return;
    long i;
    double pt1[3], pt2[3];
    long* nodes_sorted = NULL;
    long* nodes_unsorted = new long[no_points];
    for (i = 0; i < no_points; i++)
        nodes_unsorted[i] = i;
    pt1[0] = point_vector[0]->x;
    pt1[1] = point_vector[0]->y;
    pt1[2] = point_vector[0]->z;
    double* node_distances = new double[no_points];
    for (i = 1; i < no_points; i++)
    {
        pt2[0] = point_vector[i]->x;
        pt2[1] = point_vector[i]->y;
        pt2[2] = point_vector[i]->z;
        node_distances[i] = MCalcDistancePointToPoint(pt1, pt2);
    }
    nodes_sorted = TOLSortNodes1(nodes_unsorted, node_distances, no_points);
    // Reorder point vector
    std::vector<CGLPoint*> aux_point_vector;
    CGLPoint* m_pnt;
    for (i = 0; i < no_points; i++)
    {
        m_pnt = point_vector[nodes_sorted[i]];
        aux_point_vector.push_back(m_pnt);
    }
    point_vector.clear();
    point_vector = aux_point_vector;
    // Release memory
    delete[] node_distances;
    delete[] nodes_unsorted;
}
/**************************************************************************
   FEMLib-Method: GetOrderedNodeByDistance_Polyline(const long *nodes,
   const CGLPolyline *plyL, vector<double>& Distance,
   vector<int>& OrderedNode)
   Task: Prescibe the boundary values to all points of a polyline by the means
   of spline

   Programing:
   02/2004 WW Implementation
   last modification:
   23/2004 WW
   08/2005 WW Set as a member of polyline
   04/2006 WW Fix a big for polyines have more than two points
**************************************************************************/
void CGLPolyline::GetPointOrderByDistance()
{
    //	std::cout << "CGLPolyline::GetPointOrderByDistance() polyline " <<
    //getName() << "\n"; 	for (size_t k(0); k<sbuffer.size(); k++) { 		std::cout <<
    //"\tsbuffer[" << k << "]: " << sbuffer[k] << ", ibuffer[" << k << "]: " <<
    //ibuffer[k] <<"\n";
    //	}

    long number_of_nodes, i, l;
    int j, k;

    double sl = 0.0, s0 = 0.0;
    double xp, yp, zp;

    number_of_nodes = (long)ibuffer.size();
    OrderedPoint.resize(number_of_nodes);
    // Obtain fem node for groupvector
    CGLPoint *CGPa = NULL, *CGPb = NULL;
    const int SizeCGLPoint = (int)point_vector.size();

    // Reorder the node
    for (j = 0; j < number_of_nodes; j++)
        OrderedPoint[j] = j;

    // Reorder the nodes finded along polyline
    /*
       for(i=0; i<SizeCGLPoint-1; i++)
       {
       // Reorder the node
       for(j=0; j<number_of_nodes; j++)
       {
       if(ibuffer[j]!=i) continue;
       for(k=j; k<number_of_nodes; k++)
       {
       if(ibuffer[k]!=i) continue;
       s0 = sbuffer[j];
       sl = sbuffer[k];
       if(sl<s0)
       {
       l = OrderedPoint[j];
       sbuffer[k] = s0;
       sbuffer[j] = sl;
       OrderedPoint[j] = OrderedPoint[k];
       OrderedPoint[k] = l;
       }
       }
       }
       }
     */

    // Reorder the nodes within a sector
    for (j = 0; j < number_of_nodes; j++)
    {
        i = ibuffer[j];
        for (k = j; k < number_of_nodes; k++)
        {
            if (ibuffer[k] != i)
                continue;
            s0 = sbuffer[j];
            sl = sbuffer[k];
            if (sl < s0)
            {
                l = OrderedPoint[j];
                sbuffer[k] = s0;
                sbuffer[j] = sl;
                OrderedPoint[j] = OrderedPoint[k];
                OrderedPoint[k] = l;
            }
        }
    }

    sl = 0.0;
    for (i = 0; i < SizeCGLPoint - 1; i++)
    {
        for (j = 0; j < number_of_nodes; j++)
            if (ibuffer[j] == i)
                sbuffer[j] +=
                    sl;  // distance to the first point of the polyline
        CGPa = point_vector[i];
        CGPb = point_vector[i + 1];
        xp = CGPb->x - CGPa->x;
        yp = CGPb->y - CGPa->y;
        zp = CGPb->z - CGPa->z;
        sl += sqrt(xp * xp + yp * yp + zp * zp);
    }
    //
    // Reorder all nodes found
    for (j = 0; j < number_of_nodes; j++)
        for (k = j; k < number_of_nodes; k++)
        {
            if (k == j)
                continue;
            s0 = sbuffer[j];
            sl = sbuffer[k];
            if (sl < s0)
            {
                l = OrderedPoint[j];
                sbuffer[k] = s0;
                sbuffer[j] = sl;
                OrderedPoint[j] = OrderedPoint[k];
                OrderedPoint[k] = l;
            }
        }
    //	std::cout << "CGLPolyline::GetPointOrderByDistance() polyline " <<
    //getName() << "\n"; 	for (size_t k(0); k<OrderedPoint.size(); k++) {
    //		std::cout << "OrderedPoint[" << k << "]: " << OrderedPoint[k] << ",
    //sbuffer[" << k << "]: " << sbuffer[k] <<
    //", ibuffer[" << k << "]: " << ibuffer[k] <<"\n";
    //	}
}

/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
// void GEOUnselectPLY() {
//	CGLPolyline* m_ply = NULL;
//	vector<CGLPolyline*>::const_iterator p_ply;
//	p_ply = polyline_vector.begin();
//	while (p_ply != polyline_vector.end()) {
//		m_ply = *p_ply;
//		m_ply->highlighted = false;
//		++p_ply;
//	}
//}

/**************************************************************************
   GEOLib-Method
   Programing:
   10/2005 OK Implementation
**************************************************************************/
void CGLPolyline::SetPointOrderByDistance(CGLPoint* m_pnt)
{
    int i;
    // OK  double eps = 1e-3;
    //----------------------------------------------------------------------
    CGLPoint* m_pnt_i = NULL;
    std::vector<CGLPoint*> point_vector_aux;
    for (i = 0; i < (int)point_vector.size(); i++)
    {
        m_pnt_i = point_vector[i];
        // OK if(m_pnt->PointDisXY(m_pnt_i)<epsilon){
        if (m_pnt->PointDis(m_pnt_i) < epsilon)
            point_vector_aux.push_back(m_pnt_i);
    }
    if ((int)point_vector_aux.size() == 0)
        return;
    for (i = 0; i < (int)point_vector.size(); i++)
    {
        m_pnt_i = point_vector[i];
        if (m_pnt_i == point_vector_aux[0])
            continue;
        point_vector_aux.push_back(m_pnt_i);
    }
    point_vector.clear();
    point_vector = point_vector_aux;
    SortPointVectorByDistance();
    /*
       //----------------------------------------------------------------------
       sbuffer.clear();
       ibuffer.clear();
       CGLPoint* m_pnt_0 = point_vector[0];
       for(i=0;i<(int)point_vector.size();i++){
       m_pnt_i = point_vector[i];
       ibuffer.push_back(i);
       dist = m_pnt_i->PointDisXY(m_pnt_0);
       sbuffer.push_back(dist);
       }
       //----------------------------------------------------------------------
       GetPointOrderByDistance();
     */
    //----------------------------------------------------------------------
}
/**************************************************************************
   GEOLib-Method:
   Task:
   Programing:
   11/2005 CC Implementation
**************************************************************************/
void GEOWritePolylines(char* file_name)
{
    CGLPolyline* m_ply = NULL;
    for (int i = 0; i < (int)polyline_vector.size(); i++)
    {
        m_ply = polyline_vector[i];
        m_ply->Write(file_name);
    }
}

/**************************************************************************
   GeoLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
**************************************************************************/
void GEORemovePLY(CGLPolyline* m_ply)
{
    CGLPolyline* m_ply_this = NULL;
    // WW std::vector<CGLPolyline*>::const_iterator p_ply =
    // polyline_vector.begin();
    for (int i = 0; i < (int)polyline_vector.size(); i++)
    {
        m_ply_this = polyline_vector[i];
        if (m_ply_this->getName().compare(m_ply->getName()) == 0)
        {
            delete m_ply_this;
            polyline_vector.erase(polyline_vector.begin() + i);
            // i--;
            return;
        }
    }
}

/**************************************************************************
   GeoLib-Method:
   01/2006 OK Implementation based on CCs version
**************************************************************************/
// void CGLPolyline::CalcMinimumPointDistance() {
//	CGLPoint* start_point = NULL;
//	CGLPoint* end_point = NULL;
//	double m_dXMin = 1.e+19;
//	double min;
//	for (int i = 0; i < (int) point_vector.size() - 1; i++) {
//		start_point = point_vector[i];
//		end_point = point_vector[i + 1];
//		//min =
// sqrt((start_point->x-end_point->x)*(start_point->x-end_point->x)+(start_point->y-end_point->y)*(start_point->y-end_point->y));
//		min = start_point->PointDis(end_point);
//		if (min < m_dXMin)
//			m_dXMin = min;
//	}
////	minDis = m_dXMin;
//}

/**************************************************************************
   GeoSys-GUI-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
CColumn::~CColumn()
{
    for (size_t i = 0; i < line_vector.size(); i++)
        delete line_vector[i];
    line_vector.clear();
}

void CColumn::deleteLines()
{
    for (size_t j = 0; j < line_vector.size(); j++)
        delete line_vector[j];
    line_vector.clear();
}

/**************************************************************************
   GeoSys-GUI-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
void COLDelete()
{
    CColumn* m_col = NULL;
    for (int i = 0; i < (int)column_vector.size(); i++)
    {
        m_col = column_vector[i];
        delete m_col;
    }
    column_vector.clear();
}
/**************************************************************************
   GeoSys-GUI-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
void COLDeleteLines()
{
    for (int i = 0; i < (int)column_vector.size(); i++)
        column_vector[i]->deleteLines();
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
CColumn* COLGet(int col_id)
{
    CColumn* m_col = NULL;
    for (size_t i = 0; i < column_vector.size(); i++)
    {
        m_col = column_vector[i];
        if (m_col->getMatGroup() == col_id)
            return m_col;
    }
    return NULL;
}

/**************************************************************************
   GEOLib-Method:
   Programing:
   12/2005 OK Implementation
**************************************************************************/
CColumn* COLGet(const std::string& col_name)
{
    for (size_t i = 0; i < column_vector.size(); i++)
        if (column_vector[i]->getName().compare(col_name) == 0)
            return column_vector[i];
    return NULL;
}
/**************************************************************************
   GEOLib-Method:
   Programing:
   08/2006 YD Implementation
**************************************************************************/
CSoilProfile::CSoilProfile()  // YD
{
}
/**************************************************************************
   GEOLib-Method:
   Programing:
   08/2006 YD Implementation
**************************************************************************/
CSoilProfile::~CSoilProfile()  // YD
{
    for (int i = 0; i < (int)soil_type.size(); i++)
        soil_type[i] = 0;
    for (int i = 0; i < (int)soil_layer_thickness.size(); i++)
        soil_layer_thickness[i] = 0.;
    soil_type.clear();
    soil_layer_thickness.clear();
}
#ifdef RFW_FRACTURE
/**************************************************************************
   GeoLib-Method: CalcPolylineLength();

   Task: Calculate the length of a polyline, not from end-to-end, but rather
from CGLPoint-to-CGLPoint along the entire polyline. Argument: Ergebnis: double
poly_length : Length of the polyline, along the points, not from end to end
   Programing:
   05/2005 RFW Implementierung
   05/2006 RFW �nderung
**************************************************************************/
double CGLPolyline::CalcPolylineLength()
{
    double dx, dy, dz, dist = 0, poly_length = 0;

    for (int j = 0; j < (int)point_vector.size() - 1; j++)
    {
        dx = point_vector[j]->x - point_vector[j + 1]->x;
        dy = point_vector[j]->y - point_vector[j + 1]->y;
        dz = point_vector[j]->z - point_vector[j + 1]->z;

        dist = sqrt(dx * dx + dy * dy + dz * dz);
        poly_length = poly_length + dist;
    }
    return poly_length;
}
#endif
