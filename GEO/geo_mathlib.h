/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object:GEO Mathlib
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/
#ifndef geomathlib_INC
#define geomathlib_INC

#include <cmath>
#include <limits>

#define MKleinsteZahlen std::numeric_limits<double>::min()
#define MAX_ZEILEN 2048 // 256 //2048 OK
#define CSV_FILE_EXTENSIONS ".csv"
#define TEC_FILE_EXTENSIONS ".tec"

extern double EuklVek3dDist(double* x, double* y);
extern double EuklVek3dDistCoor(double x1, double y1, double z1, double x2, double y2, double z2);
extern double Vek3dDistCoor(double x1, double y1, double z1, double x2, double y2, double z2, int norm);
extern int M3KreuzProdukt(double* vec1, double* vec2, double* vec);
extern double MBtrgVec(double* vec, long n);
extern double MSkalarprodukt(double* vec1, double* vec2, long g);
extern double M3Determinante(double* m);
extern double M4Determinante(double* m);
extern double CalcTetraederVolume(double* x, double* y, double* z); // CC
extern double CalcPyramidVolume(double* x, double* y, double* z); // CC
extern double CalcPrismVolume(double* x, double* y, double* z); // CC
extern double MCalcProjectionOfPointOnPlane(double* pt, double* e1, double* e2, double* e3, double* proj);
extern double MCalcProjectionOfPointOnPlane(double* pt, double* e1, double* e2, double* e3, double* proj);
extern double MCalcDistancePointToPoint(double* pt1, double* pt2);
extern long* TOLSortNodes1(long*, double*, int);
extern int MPhi2D(double* vf, double r, double s);
#ifdef RFW_FRACTURE
extern bool LineSegmentIntersection(vector<double>,
                                    vector<double>,
                                    vector<double>,
                                    vector<double>,
                                    vector<double>&); // RFW 04/2005
#endif

#endif
