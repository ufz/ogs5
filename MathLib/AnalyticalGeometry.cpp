/*
 * \file AnalyticalGeometry.cpp
 *
 *  Created on: Mar 17, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cmath>
#include <cstdlib>  // for exit
#include <fstream>
#include <limits>
#include <list>

// Base
#include "quicksort.h"
#include "swap.h"

// GEO
#include "AxisAlignedBoundingBox.h"
#include "Polyline.h"
#include "Triangle.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "GaussAlgorithm.h"
#include "MathTools.h"
#include "Matrix.h"  // for transformation matrix
#include "max.h"

namespace MathLib
{
Orientation getOrientation(const double& p0_x, const double& p0_y,
                           const double& p1_x, const double& p1_y,
                           const double& p2_x, const double& p2_y)
{
    double h1((p1_x - p0_x) * (p2_y - p0_y));
    double h2((p2_x - p0_x) * (p1_y - p0_y));

    double tol(sqrt(std::numeric_limits<double>::min()));
    if (fabs(h1 - h2) <= tol * max(fabs(h1), fabs(h2)))
        return COLLINEAR;
    if (h1 - h2 > 0.0)
        return CCW;

    return CW;
}

Orientation getOrientation(const GEOLIB::Point* p0, const GEOLIB::Point* p1,
                           const GEOLIB::Point* p2)
{
    return getOrientation((*p0)[0], (*p0)[1], (*p1)[0], (*p1)[1], (*p2)[0],
                          (*p2)[1]);
}

bool lineSegmentIntersect(const GEOLIB::Point& a, const GEOLIB::Point& b,
                          const GEOLIB::Point& c, const GEOLIB::Point& d,
                          GEOLIB::Point& s)
{
    //*** in order to make the intersection test more stable
    // compute bounding box for the points
    GEOLIB::AABB aabb;
    aabb.update(a);
    aabb.update(b);
    aabb.update(c);
    aabb.update(d);
    // transforming coordinates to interval [0,1]x[0,1]x[0,1]
    double delta(aabb.getMaxPoint()[0] - aabb.getMinPoint()[0]);
    const double tmp(aabb.getMaxPoint()[1] - aabb.getMinPoint()[1]);
    if (delta < tmp)
    {
        delta = tmp;
    }
    GEOLIB::Point a_cpy((a[0] - aabb.getMinPoint()[0]) / delta,
                        (a[1] - aabb.getMinPoint()[1]) / delta, 0.0);
    GEOLIB::Point b_cpy((b[0] - aabb.getMinPoint()[0]) / delta,
                        (b[1] - aabb.getMinPoint()[1]) / delta, 0.0);
    GEOLIB::Point c_cpy((c[0] - aabb.getMinPoint()[0]) / delta,
                        (c[1] - aabb.getMinPoint()[1]) / delta, 0.0);
    GEOLIB::Point d_cpy((d[0] - aabb.getMinPoint()[0]) / delta,
                        (d[1] - aabb.getMinPoint()[1]) / delta, 0.0);

    Matrix<double> mat(2, 2);
    mat(0, 0) = b_cpy[0] - a_cpy[0];
    mat(1, 0) = b_cpy[1] - a_cpy[1];
    mat(0, 1) = c_cpy[0] - d_cpy[0];
    mat(1, 1) = c_cpy[1] - d_cpy[1];

    // check if vectors are parallel
    double eps(sqrt(std::numeric_limits<double>::min()));
    if (fabs(mat(1, 1)) < eps)
    {
        // vector (D-C) is parallel to x-axis
        if (fabs(mat(0, 1)) < eps)
            // vector (B-A) is parallel to x-axis
            return false;
    }
    else
    {
        // vector (D-C) is not parallel to x-axis
        if (fabs(mat(0, 1)) >= eps)
            // vector (B-A) is not parallel to x-axis
            // \f$(B-A)\f$ and \f$(D-C)\f$ are parallel iff there exists
            // a constant \f$c\f$ such that \f$(B-A) = c (D-C)\f$
            if (fabs(mat(0, 0) / mat(0, 1) - mat(1, 0) / mat(1, 1)) <
                eps * fabs(mat(0, 0) / mat(0, 1)))
                return false;
    }

    double* rhs(new double[2]);
    rhs[0] = c_cpy[0] - a_cpy[0];
    rhs[1] = c_cpy[1] - a_cpy[1];

    GaussAlgorithm<Matrix<double> > lu_solver(mat, mat.getNRows());
    lu_solver.execute(rhs);
    if (0 <= rhs[0] && rhs[0] <= 1.0 && 0 <= rhs[1] && rhs[1] <= 1.0)
    {
        s[0] = a[0] + rhs[0] * (b[0] - a[0]);
        s[1] = a[1] + rhs[0] * (b[1] - a[1]);
        s[2] = a[2] + rhs[0] * (b[2] - a[2]);

        // check z component
        double z0(a[2] - d[2]),
            z1(rhs[0] * (b[2] - a[2]) + rhs[1] * (d[2] - c[2]));
        delete[] rhs;
        if (std::fabs(z0 - z1) < eps)
            return true;
        else
            return false;
    }
    else
    {
        delete[] rhs;
    }
    return false;
}

bool lineSegmentsIntersect(const GEOLIB::Polyline* ply, size_t& idx0,
                           size_t& idx1, GEOLIB::Point& intersection_pnt)
{
    size_t n_segs(ply->getNumberOfPoints() - 1);
    /**
     * computing the intersections of all possible pairs of line segments of the
     * given polyline as follows: let the segment \f$s_1 = (A,B)\f$ defined by
     * \f$k\f$-th and \f$k+1\f$-st point of the polyline and segment \f$s_2 =
     * (C,B)\f$ defined by \f$j\f$-th and \f$j+1\f$-st point of the polyline,
     * \f$j>k+1\f$
     */
    for (size_t k(0); k < n_segs - 2; k++)
        for (size_t j(k + 2); j < n_segs; j++)
        {
            if (k != 0 || j < n_segs - 1)
            {
                if (lineSegmentIntersect(*(*ply)[k], *(*ply)[k + 1], *(*ply)[j],
                                         *(*ply)[j + 1], intersection_pnt))
                {
                    idx0 = k;
                    idx1 = j;
                    return true;
                }
            }
        }
    return false;
}

double calcTriangleArea(GEOLIB::Point const& a, GEOLIB::Point const& b,
                        GEOLIB::Point const& c)
{
    double const u[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    double const v[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    double w[3];
    MathLib::crossProd(u, v, w);
    return 0.5 * (sqrt(MathLib::scpr(w, w, 3)));
}

bool isPointInTriangle(const double q[3], const double a[3], const double b[3],
                       const double c[3], double eps)
{
    if (sqrt(MathLib::sqrDist(q, a)) < eps ||
        sqrt(MathLib::sqrDist(q, b)) < eps ||
        sqrt(MathLib::sqrDist(q, c)) < eps)
    {
        return true;
    }

    MathLib::Vector const v(a, b);
    MathLib::Vector const w(a, c);

    MathLib::Matrix<double> mat(2, 2);
    mat(0, 0) = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    mat(0, 1) = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
    mat(1, 0) = mat(0, 1);
    mat(1, 1) = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
    double y[2] = {
        v[0] * (q[0] - a[0]) + v[1] * (q[1] - a[1]) + v[2] * (q[2] - a[2]),
        w[0] * (q[0] - a[0]) + w[1] * (q[1] - a[1]) + w[2] * (q[2] - a[2])};

    MathLib::GaussAlgorithm<Matrix<double> > gauss(mat, mat.getNRows());
    gauss.execute(y);

    const double lower(-eps);
    const double upper(1 + eps);

    if (lower <= y[0] && y[0] <= upper && lower <= y[1] && y[1] <= upper &&
        y[0] + y[1] <= upper)
    {
        double const q_projected[3] = {a[0] + y[0] * v[0] + y[1] * w[0],
                                       a[1] + y[0] * v[1] + y[1] * w[1],
                                       a[2] + y[0] * v[2] + y[1] * w[2]};
        if (sqrt(MathLib::sqrDist(q, q_projected)) < eps)
        {
            return true;
        }
    }

    return false;
}

bool isPointInTriangle(const GEOLIB::Point* p, const GEOLIB::Point* a,
                       const GEOLIB::Point* b, const GEOLIB::Point* c,
                       double eps)
{
    return isPointInTriangle(p->getData(), a->getData(), b->getData(),
                             c->getData(), eps);
}

// NewellPlane from book Real-Time Collision detection p. 494
void getNewellPlane(const std::vector<GEOLIB::Point*>& pnts,
                    Vector& plane_normal, double& d)
{
    d = 0;
    Vector centroid;
    size_t n_pnts(pnts.size());
    for (size_t i(n_pnts - 1), j(0); j < n_pnts; i = j, j++)
    {
        plane_normal[0] +=
            ((*(pnts[i]))[1] - (*(pnts[j]))[1]) *
            ((*(pnts[i]))[2] + (*(pnts[j]))[2]);  // projection on yz
        plane_normal[1] +=
            ((*(pnts[i]))[2] - (*(pnts[j]))[2]) *
            ((*(pnts[i]))[0] + (*(pnts[j]))[0]);  // projection on xz
        plane_normal[2] +=
            ((*(pnts[i]))[0] - (*(pnts[j]))[0]) *
            ((*(pnts[i]))[1] + (*(pnts[j]))[1]);  // projection on xy

        centroid += *(pnts[j]);
    }

    plane_normal *= 1.0 / plane_normal.Length();
    d = centroid.Dot(plane_normal) / n_pnts;
}

void rotatePointsToXY(Vector& plane_normal, std::vector<GEOLIB::Point*>& pnts)
{
    double small_value(sqrt(std::numeric_limits<double>::min()));
    if (fabs(plane_normal[0]) < small_value &&
        fabs(plane_normal[1]) < small_value)
        return;

    // *** some frequently used terms ***
    // sqrt (v_1^2 + v_2^2)
    double h0(sqrt(plane_normal[0] * plane_normal[0] +
                   plane_normal[1] * plane_normal[1]));
    // 1 / sqrt (v_1^2 + v_2^2)
    double h1(1 / h0);
    // 1 / sqrt (h0 + v_3^2)
    double h2(1.0 / sqrt(h0 + plane_normal[2] * plane_normal[2]));

    Matrix<double> rot_mat(3, 3);
    // calc rotation matrix
    rot_mat(0, 0) = plane_normal[2] * plane_normal[0] * h2 * h1;
    rot_mat(0, 1) = plane_normal[2] * plane_normal[1] * h2 * h1;
    rot_mat(0, 2) = -h0 * h2;
    rot_mat(1, 0) = -plane_normal[1] * h1;
    rot_mat(1, 1) = plane_normal[0] * h1;
    rot_mat(1, 2) = 0.0;
    rot_mat(2, 0) = plane_normal[0] * h2;
    rot_mat(2, 1) = plane_normal[1] * h2;
    rot_mat(2, 2) = plane_normal[2] * h2;

    double* tmp(NULL);
    size_t n_pnts(pnts.size());
    for (size_t k(0); k < n_pnts; k++)
    {
        tmp = rot_mat * pnts[k]->getData();
        for (size_t j(0); j < 3; j++)
            (*(pnts[k]))[j] = tmp[j];
        delete[] tmp;
    }

    tmp = rot_mat * plane_normal.getData();
    for (size_t j(0); j < 3; j++)
        plane_normal[j] = tmp[j];

    delete[] tmp;
}

void rotatePointsToXZ(Vector& n, std::vector<GEOLIB::Point*>& pnts)
{
    double small_value(sqrt(std::numeric_limits<double>::min()));
    if (fabs(n[0]) < small_value && fabs(n[1]) < small_value)
        return;

    // *** some frequently used terms ***
    // n_1^2 + n_2^2
    const double h0(n[0] * n[0] + n[1] * n[1]);
    // 1 / sqrt (n_1^2 + n_2^2)
    const double h1(1.0 / sqrt(h0));
    // 1 / sqrt (n_1^2 + n_2^2 + n_3^2)
    const double h2(1.0 / sqrt(h0 + n[2] * n[2]));

    Matrix<double> rot_mat(3, 3);
    // calc rotation matrix
    rot_mat(0, 0) = n[1] * h1;
    rot_mat(0, 1) = -n[0] * h1;
    rot_mat(0, 2) = 0.0;
    rot_mat(1, 0) = n[0] * h2;
    rot_mat(1, 1) = n[1] * h2;
    rot_mat(1, 2) = n[2] * h2;
    rot_mat(2, 0) = n[0] * n[2] * h1 * h2;
    rot_mat(2, 1) = n[1] * n[2] * h1 * h2;
    rot_mat(2, 2) = -sqrt(h0) * h2;

    double* tmp(NULL);
    size_t n_pnts(pnts.size());
    for (size_t k(0); k < n_pnts; k++)
    {
        tmp = rot_mat * pnts[k]->getData();
        for (size_t j(0); j < 3; j++)
            (*(pnts[k]))[j] = tmp[j];
        delete[] tmp;
    }

    tmp = rot_mat * n.getData();
    for (size_t j(0); j < 3; j++)
        n[j] = tmp[j];

    delete[] tmp;
}

}  // end namespace MathLib
