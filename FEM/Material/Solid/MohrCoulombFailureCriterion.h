/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file MohrCoulombFailureCriterion.h
 *
 * Created on February 6, 2019, 1:40 PM
 *
 */

#pragma once

namespace SolidProp
{
class MohrCoulombFailureCriterion
{
public:
    MohrCoulombFailureCriterion(const double c, const double angle);

    /**
     *
     * @param s                  Stress.
     * @param nstress_components Number of stress components.
     * @return
     */
    double getFailureIndex(double const* const s,
                           const int nstress_components) const;

private:
    const double _c;      /// apparent cohesion.
    const double _angle;  /// The angle of internal friction.
};

}  // end of namespace
