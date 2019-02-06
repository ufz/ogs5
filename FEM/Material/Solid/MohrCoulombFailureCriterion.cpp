/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file MohrCoulombFailureCriterion.cpp
 *
 * Created on February 6, 2019, 1:40 PM
 *
 */

#include "MohrCoulombFailureCriterion.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "StressAuxiliaryFunctions.h"

namespace SolidProp
{
MohrCoulombFailureCriterion::MohrCoulombFailureCriterion(const double c,
                                                         const double angle)
    : _c(c), _angle(angle * pi / 180.0)
{
}

double MohrCoulombFailureCriterion::getFailureIndex(
    double const* const s, const int nstress_components) const
{
    double const* const principle_stresses =
        getPrincipleStresses(s, nstress_components);

    double s_min = std::numeric_limits<double>::max();
    double s_max = -std::numeric_limits<double>::max();

    for (int i = 0; i < 3; i++)
    {
        s_max = std::max(s_max, principle_stresses[i]);
        s_min = std::min(s_min, principle_stresses[i]);
    }

    const double tau_max = 0.5 * std::abs(s_max - s_min);
    const double s_n = -0.5 * (s_max + s_min);

    return tau_max / (_c + s_n * std::tan(_angle));
}

}  // end of namespace
