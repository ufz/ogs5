/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file DamageZonePermeability.cpp
 *
 * Created on February 6, 2019, 3:25 PM
 *
 */

#include "DamageZonePermeability.h"

#include <cmath>

#include "Material/Solid/MohrCoulombFailureCriterion.h"
namespace PorousMediumProperty
{
void DamageZonePermeability::computeDamageZonePermeability(
    double* intrinsic_permeability,
    SolidProp::MohrCoulombFailureCriterion const& failure_criterion,
    double const* const stress, const int dim)
{
    const int n_stresses = (dim == 3) ? 6 : 4;
    const double failure_index =
        failure_criterion.getFailureIndex(stress, n_stresses);

    if (failure_index < 1.0)
        return;

    const double extra_k = _a * std::exp(_b * failure_index);
    for (int i = 0; i < dim; i++)
    {
        intrinsic_permeability[i * dim + i] += extra_k;
    }
}
}  // End of the namespace
