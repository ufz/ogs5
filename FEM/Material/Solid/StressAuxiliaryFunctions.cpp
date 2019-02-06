/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file StressAuxiliaryFunctions.cpp
 *
 * Created on February 6, 2019, 12:08 PM
 *
 */

#include "StressAuxiliaryFunctions.h"

#include <cmath>

namespace SolidProp
{
void getDeviatoricStess(double const* const stress,
                        const int nstress_components, double* s)
{
    for (int i = 0; i < nstress_components; i++)
        s[i] = stress[i];

    const double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
    for (int i = 0; i < 3; i++)
        s[i] -= mean_stress;
}

double getStressNorm(double const* const s, const int nstress_components)
{
    double ns = 0.0;
    for (int i = 0; i < 3; i++)
        ns += s[i] * s[i];
    for (int i = 3; i < nstress_components; i++)
        ns += 2.0 * s[i] * s[i];
    return std::sqrt(ns);
}

double* getPrincipleStresses(double const* const s, const int nstress_components)
{
    static double principle_stress[3];
    const double s11 = s[0];
    const double s22 = s[1];
    const double s33 = s[2];
    const double s12 = s[3];

    const double I1 = s11 + s22 + s33;
    double I2 = s11 * s22 + s22 * s33 + s33 * s11 - s12 * s12;
    double I3 = s11 * s22 * s33 - s12 * s12 * s33;
    if (nstress_components == 6)
    {
        const double s13 = s[4];
        const double s23 = s[5];
        I2 += -s13 * s13 - s23 * s23;
        I3 += 2 * s12 * s23 * s13 - s23 * s23 * s11 - s13 * s13 * s22;
    }

    double alpha =
        std::acos(0.5 * (2 * I1 * I1 * I1 - 9. * I1 * I2 + 27.0 * I3) /
                  std::pow(I1 * I1 - 3.0 * I2, 1.5));
    const double I1_d3 = I1 / 3.0;
    const double fac = 2. * std::sqrt(I1 * I1 - 3 * I2) / 3.0;

    principle_stress[0] = I1_d3 + fac * std::cos(alpha);
    principle_stress[1] = I1_d3 + fac * std::cos(alpha - 2.0 * pi / 3.0);
    principle_stress[2] = I1_d3 + fac * std::cos(alpha - 4.0 * pi / 3.0);

    return principle_stress;
}

}  // end of namespace
