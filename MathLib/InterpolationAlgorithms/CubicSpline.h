/*
 * CubicSpline.h
 *
 *  Created on: Jul 27, 2010
 *      Author: TF (moved class CubicSpline from geo_mathlib.{h.cpp})
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_

#include <cstddef>
#include <vector>

namespace MathLib
{
class CubicSpline
{
public:
    CubicSpline(const std::vector<double>& s, const std::vector<double>& val);
    ~CubicSpline();
    double interpolation(double x) const;

private:
    const size_t n;
    double* bb;
    double* cc;
    double* dd;
    std::vector<double> xx;
    std::vector<double> yy;
    void computeCoefficents();
};
}  // end namespace MathLib

#endif /* CUBICSPLINE_H_ */
