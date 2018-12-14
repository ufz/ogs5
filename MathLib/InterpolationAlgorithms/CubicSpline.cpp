/*
 * CubicSpline.cpp
 *
 *  Created on: Jul 27, 2010
 *      Author: TF (moved class CubicSpline from geo_mathlib.{h.cpp})
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "CubicSpline.h"
#include "MathTools.h"

namespace MathLib
{
CubicSpline::CubicSpline(const std::vector<double>& s,
                         const std::vector<double>& val)
    : n(s.size()), bb(new double[n]), cc(new double[n]), dd(new double[n])
{
    xx = s;
    yy = val;
    computeCoefficents();
}

CubicSpline::~CubicSpline()
{
    delete[] bb;
    delete[] cc;
    delete[] dd;
    bb = NULL;
    cc = NULL;
    dd = NULL;
}

double CubicSpline::interpolation(double x) const
{
    double val = 0.0;
    double y_max = -1.0e14;
    double y_min = 1.0e14;
    double x0 = 0.0, y0 = 0.0, x1 = 0.0, y1 = 0.0;
    bool withinR = false;
    for (size_t i = 0; i < n - 1; i++)
        if (x >= xx[i] && x < xx[i + 1])
        {
            // 07/2010 TF
            //			val = yy[i] + bb[i] * (x - xx[i]) + cc[i] * pow(x -
            //xx[i], 2.0)
            //					+ dd[i] * pow(x - xx[i], 3.0);
            double t(x - xx[i]);
            // employing Horner-Schema in order to save multiplications
            val = yy[i] + t * (bb[i] + t * (cc[i] + t * dd[i]));

            // Check the local range
            if (yy[i] > y_max)
                y_max = yy[i];
            if (yy[i + 1] > y_max)
                y_max = yy[i + 1];
            if (yy[i] < y_min)
                y_min = yy[i];
            if (yy[i + 1] < y_min)
                y_min = yy[i + 1];

            // Linear interpolation
            if (val < y_min || val > y_max)
                val = yy[i] + (yy[i + 1] - yy[i]) * t / (xx[i + 1] - xx[i]);
            withinR = true;
            break;
        }

    if (withinR)
        return val;
    //-------------------------------------
    // Extrapolate
    //-------------------------------------
    // Compute the global range
    y_max = -1.0e14;
    y_min = 1.0e14;
    for (size_t i = 0; i < n; i++)
    {
        if (yy[i] > y_max)
            y_max = yy[i];
        if (yy[i] < y_min)
            y_min = yy[i];
    }
    if ((x >= xx[n - 1]) || (x <= xx[0]))
    {
        if (x <= xx[0])
        {
            x0 = xx[0];
            y0 = yy[0];
            x1 = xx[1];
            y1 = yy[1];
        }
        else if ((x >= xx[n - 1]))
        {
            x0 = xx[n - 2];
            y0 = yy[n - 2];
            x1 = xx[n - 1];
            y1 = yy[n - 1];
        }
        val = y0 + (y1 - y0) * (x - x0) / (x1 - x0);
        if (val < y_min)
            val = y_min;
        if (val > y_max)
            val = y_max;
    }

    return val;
}

/*******************************************************************
   the coefficients b(i), c(i), and d(i) are computed
   for a cubic interpolating spline

    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3

    for  x(i) .le. x .le. x(i+1)

   input..

    n = the number of data points or knots (n.ge.2)
    x = the abscissas of the knots in strictly increasing order
    y = the ordinates of the knots

   output..

    b, c, d  = arrays of spline coefficients as defined above.

   using  p  to denote differentiation,

    y(i) = s(x(i))
    b(i) = sp(x(i))
    c(i) = spp(x(i))/2
    d(i) = sppp(x(i))/6  (derivative from the right)

   the accompanying function subprogram  seval  can be used
   to evaluate the spline.
*******************************************************************/
void CubicSpline::computeCoefficents()
{
    if (n < 2)
    {
        printf("Dimension can not be less than 3 in spline");
        abort();
    }

    if (n == 2)
    {
        bb[0] = (yy[1] - yy[0]) / (xx[1] - xx[0]);
        cc[0] = 0.0;
        dd[0] = 0.0;
        bb[1] = bb[0];
        cc[1] = 0.0;
        dd[1] = 0.0;
    }
    else
    {
        /***************************************************
           set up tridiagonal system
           b = diagonal, d = offdiagonal, c = right hand side.
         ****************************************************/
        dd[0] = xx[1] - xx[0];
        cc[1] = (yy[1] - yy[0]) / dd[0];
        for (size_t i = 1; i < n - 1; i++)
        {
            dd[i] = xx[i + 1] - xx[i];
            bb[i] = 2.0 * (dd[i - 1] + dd[i]);
            cc[i + 1] = (yy[i + 1] - yy[i]) / dd[i];
            cc[i] = cc[i + 1] - cc[i];
        }

        /******************************************************
           end conditions.  third derivatives at  x[1]  and  x[n]
           obtained from divided differences
        ******************************************************/
        bb[0] = -dd[0];
        bb[n - 1] = -dd[n - 2];
        cc[0] = 0.0;
        cc[n - 1] = 0.0;
        if (n > 3)
        {
            cc[0] = cc[2] / (xx[3] - xx[1]) - cc[1] / (xx[2] - xx[0]);
            cc[n - 1] = cc[n - 2] / (xx[n - 1] - xx[n - 3]) -
                        cc[n - 3] / (xx[n - 2] - xx[n - 4]);
            cc[0] = cc[0] * fastpow(dd[0], 2) / (xx[3] - xx[0]);
            cc[n - 1] =
                -cc[n - 1] * fastpow(dd[n - 2], 2) / (xx[n - 1] - xx[n - 4]);
        }

        // *** forward elimination
        double t;
        for (size_t i = 1; i < n; i++)
        {
            t = dd[i - 1] / bb[i - 1];
            bb[i] = bb[i] - t * dd[i - 1];
            cc[i] = cc[i] - t * cc[i - 1];
        }
        // *** back substitution
        cc[n - 1] = cc[n - 1] / bb[n - 1];
        for (size_t ib = 0; ib < n - 1; ib++)
        {
            size_t i = n - 1 - ib;
            cc[i] = (cc[i] - dd[i] * cc[i + 1]) / bb[i];
        }
        //
        //  c[i] is now the sigma[i] of the text
        //
        //  compute polynomial coefficients
        //
        bb[n - 1] = (yy[n - 1] - yy[n - 1 - 1]) / dd[n - 1 - 1] +
                    dd[n - 1 - 1] * (cc[n - 1 - 1] + 2.0 * cc[n - 1]);
        for (size_t i = 0; i < n - 1; i++)
        {
            bb[i] =
                (yy[i + 1] - yy[i]) / dd[i] - dd[i] * (cc[i + 1] + 2.0 * cc[i]);
            dd[i] = (cc[i + 1] - cc[i]) / dd[i];
            cc[i] = 3.0 * cc[i];
        }
        cc[n - 1] = 3.0 * cc[n - 1];
        dd[n - 1] = dd[n - 2];
    }
}
}  // end namespace MathLib
