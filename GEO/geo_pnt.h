/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Point
   Task:
   Programing:
   07/2003 OK/CC/TK/WW GEOLib1
   07/2005 CC/OK GEOLib2 Design
**************************************************************************/
#ifndef gs_pnt_INC
#define gs_pnt_INC
// C++ STL
#include "MathTools.h"
#include "geo_mathlib.h"
#include <string>
#include <vector>
//
#include <fstream>
#include <iostream>
#include <sstream>

// LB TODO: trunk transition hack
std::string get_sub_string(std::string buffer, std::string delimiter, int pos1,
                           int* pos2);

/*---------------------------------------------------------------*/
class CGLPoint
{
private:
    /** geometry */
    double data[3];

public:
    double getPropert() const { return _propert; }
    void setPropert(double propert) { _propert = propert; }
    //----------------------------------------------------------------------
    // Properties
    // ID
    std::string name;
    long id;  // CC

    /* *** get and set *** */
    /** get the x coordinate of the point */
    double getX() const { return data[0]; }
    /** get the y coordinate of the point */
    double getY() const { return data[1]; }
    /** get the z coordinate of the point */
    double getZ() const { return data[2]; }
    /** set the x coordinate of the point */
    void setX(double v) { data[0] = v; }
    /** set the y coordinate of the point */
    void setY(double v) { data[1] = v; }
    /** set the z coordinate of the point */
    void setZ(double v) { data[2] = v; }
    /** get the x,y and z coordinate of the point */
    const double* getPoint() const { return data; }
    // do not use this attributes directly - use the getter and setter
    // methods!!!
    double& x;
    double& y;
    double& z;

    // double epsilon;
    double length;                  // OK well bore depth in 2D modells
    long first_identical_id;        // TK
    long old_id, new_id;            // TK
    int nb_of_ply;                  // TK Number of Polylines using this point
    long number_of_doubled_points;  // TK
    // Meshing
    int index_msh;          // WW
    double min_seg_length;  // TK for polylines
    double max_seg_length;  // TK
    double mesh_density;
    // Properties
    int type;   // OK4801
    long node;  // OK
    int mat;    // CC9999
    double value;

    // Display
    // bool highlighted;
    int x_pix, y_pix;
    int circle_pix;
    // int display_mode;
    int m_color[3];
    // bool selected;
    int plg_hightlight_seg;
    //----------------------------------------------------------------------
    // Methods
    // Create
    CGLPoint(void);
    // 07/2010 TF
    /**
     * constructor for CGLPoint object
     * @param x the first coordinate
     * @param y the second coordinate
     * @param z the third coordinate
     */
    CGLPoint(double x, double y, double z);
    /**
     * constructor for CGLPoint object
     * @param coordinates the field contains the coordinates
     */
    CGLPoint(const double* coordinates);

    ~CGLPoint(void);

    // Access
    void SetIndex(int L_index) { index_msh = L_index; }
    int GetIndex() { return index_msh; }
    // I/O
    // void Write(char*);
    std::ios::pos_type Read(std::ifstream*, int&);  // CC
    // GEO
    int IsPointExist();            // CC
    CGLPoint* Exist();             // OK
    double PointDis(CGLPoint*);    // CC
    double PointDisXY(CGLPoint*);  // OK
    bool IsInsidePolygonPlain(double*, double*, double*,
                              long);  // CC IsPointInsideSurface()?
    bool IsInsideTriangle(double*, double*, double*);   // CC
    bool IsInTriangleXYProjection(double*, double*);    // CC
    bool IsInsideRectangle(double*, double*, double*);  // CC
    bool IsInsidePrism(double*, double*, double*);      // CC
private:
    double _propert;
};
//------------------------------------------------------------------------
// Properties
extern std::vector<CGLPoint*> GetPointsVector(void);  // CC
extern std::vector<CGLPoint*> gli_points_vector;
extern std::vector<CGLPoint*> pnt_properties_vector;  // OK

//------------------------------------------------------------------------
// Remove
extern void GEORemoveAllPoints();
extern void GEORemovePoint(long);  // CC
//........................................................................
// I/O
extern void GEOReadPoints(const std::string& file_name_path_base);  // CC
extern void GEOReadPointProperties(const std::string&);
extern void GEOWritePoints(char* file_name);  // CC
//........................................................................
// Access
// extern CGLPoint* GEOGetPointByName(const std::string &);//CC
extern CGLPoint* GEOGetPointById(long);  // CC
//........................................................................
// GEO
extern long GEOPointID();
extern void GEO_Search_DoublePoints(double);  // TK
extern std::vector<CGLPoint*> GEOLIB_SetGLIPoints_Vector(
    std::vector<CGLPoint*> gl_point);  // TK
extern double AngleSumPointInsideTriangle(double* point, double* tri_p1,
                                          double* tri_p2, double* tri_p3,
                                          double tolerance);
extern void GEOCalcPointMinMaxCoordinates();  // OK

/** tests if the point p is near by a point of the vector vec
 * \param vec the point vector
 * \param p the point
 * \param tol value for stopping criterion (||vec[i] - p||^2 < tol)
 * \returns pointer to the point
 * */
CGLPoint* isPointInPointVector(const std::vector<CGLPoint*>& vec,
                               const CGLPoint* const p, double tol = 1e-3);

//........................................................................
// variables
extern double pnt_x_min;
extern double pnt_x_max;
extern double pnt_y_min;
extern double pnt_y_max;
extern double pnt_z_min;
extern double pnt_z_max;
// Id
extern long pnt_id_last;
#endif
