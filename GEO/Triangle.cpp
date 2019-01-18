/*
 * \file Triangle.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Triangle.h"

// MathLib
#include "MathTools.h"
#include "AnalyticalGeometry.h"
#include "Matrix.h"
#include "Vector3.h"
#include "GaussAlgorithm.h"
#include "mathlib.h"

namespace GEOLIB
{
Triangle::Triangle(std::vector<Point*> const& pnt_vec)
    : _pnts(pnt_vec), _initialized(false), _longest_edge(0.0)
{
    _pnt_ids[0] = std::numeric_limits<size_t>::max();
    _pnt_ids[1] = std::numeric_limits<size_t>::max();
    _pnt_ids[2] = std::numeric_limits<size_t>::max();
    for (int i = 0; i < 3; ++i)
        _normal_vector[i] = std::numeric_limits<double>::max();
}

Triangle::Triangle(std::vector<Point*> const& pnt_vec, size_t pnt_a,
                   size_t pnt_b, size_t pnt_c)
    : _pnts(pnt_vec), _initialized(true), _longest_edge(0.0)
{
    _pnt_ids[0] = pnt_a;
    _pnt_ids[1] = pnt_b;
    _pnt_ids[2] = pnt_c;
    _longest_edge = MathLib::sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
    double tmp(MathLib::sqrDist(_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
    if (tmp > _longest_edge)
        _longest_edge = tmp;
    tmp = MathLib::sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
    if (tmp > _longest_edge)
        _longest_edge = tmp;
    _longest_edge = sqrt(_longest_edge);
    calculateNormal();
}

void Triangle::setTriangle(size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
    assert(pnt_a < _pnts.size() && pnt_b < _pnts.size() &&
           pnt_c < _pnts.size());
    _pnt_ids[0] = pnt_a;
    _pnt_ids[1] = pnt_b;
    _pnt_ids[2] = pnt_c;

    _longest_edge = MathLib::sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
    double tmp(MathLib::sqrDist(_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
    if (tmp > _longest_edge)
        _longest_edge = tmp;
    tmp = MathLib::sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
    if (tmp > _longest_edge)
        _longest_edge = tmp;
    _longest_edge = sqrt(_longest_edge);
}

bool Triangle::containsPoint(const double* pnt, double eps) const
{
    GEOLIB::Point const p(pnt);
    return MathLib::isPointInTriangle(
        &p, _pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]], eps);
}

bool Triangle::containsPoint2D(const double* pnt) const
{
    GEOLIB::Point const& a(*(_pnts[_pnt_ids[0]]));
    GEOLIB::Point const& b(*(_pnts[_pnt_ids[1]]));
    GEOLIB::Point const& c(*(_pnts[_pnt_ids[2]]));

    // criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
    MathLib::Matrix<double> mat(2, 2);
    mat(0, 0) = b[0] - a[0];
    mat(0, 1) = c[0] - a[0];
    mat(1, 0) = b[1] - a[1];
    mat(1, 1) = c[1] - a[1];
    double y[2] = {pnt[0] - a[0], pnt[1] - a[1]};

    MathLib::GaussAlgorithm<MathLib::Matrix<double>> gauss(mat, mat.getNRows());
    gauss.execute(y);

    const double delta(std::numeric_limits<double>::epsilon());
    const double upper(1 + delta);

    // check if u0 and u1 fulfills the condition (with some delta)
    if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper &&
        y[0] + y[1] <= upper)
        return true;
    return false;
}

void getPlaneCoefficients(Triangle const& tri, double c[3])
{
    GEOLIB::Point const& p0(*(tri.getPoint(0)));
    GEOLIB::Point const& p1(*(tri.getPoint(1)));
    GEOLIB::Point const& p2(*(tri.getPoint(2)));
    MathLib::Matrix<double> mat(3, 3);
    mat(0, 0) = p0[0];
    mat(0, 1) = p0[1];
    mat(0, 2) = 1.0;
    mat(1, 0) = p1[0];
    mat(1, 1) = p1[1];
    mat(1, 2) = 1.0;
    mat(2, 0) = p2[0];
    mat(2, 1) = p2[1];
    mat(2, 2) = 1.0;
    c[0] = p0[2];
    c[1] = p1[2];
    c[2] = p2[2];

    MathLib::GaussAlgorithm<MathLib::Matrix<double>> gauss(mat, mat.getNRows());
    gauss.execute(c);
}

void Triangle::calculateNormal()
{
    GEOLIB::Point const& a(*(_pnts[_pnt_ids[0]]));
    GEOLIB::Point const& b(*(_pnts[_pnt_ids[1]]));
    GEOLIB::Point const& c(*(_pnts[_pnt_ids[2]]));

    const double v1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    MathLib::crossProd(v1, v2, _normal_vector);

    // normalize normal vector
    NormalizeVector(_normal_vector, 3);
}

std::ostream& operator<<(std::ostream& out, GEOLIB::Triangle const& tri)
{
    out << *tri.getPoint(0) << " " << *tri.getPoint(1) << " "
        << *tri.getPoint(2);
    return out;
}

}  // end namespace GEOLIB
