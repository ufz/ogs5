/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file DamageZonePermeability.h
 *
 * Created on February 6, 2019, 3:25 PM
 *
 */

#pragma once

namespace SolidProp
{
class MohrCoulombFailureCriterion;
}
namespace PorousMediumProperty
{
/**
 * Permeability model for damage zone, which is determined by the Mohr-Coulomb
 * failure index as
 *   \f$ k=k_0 + a \mbox{e}^{b\,f} \f$
 *   where \f$k_0\f$ is the intrinsic permeability, \f$a\f$ and \f$b\f$ are
 *   coefficient, and \f$ f \f$ is the failure index.
 */
class DamageZonePermeability
{
public:
    DamageZonePermeability(const double a, const double b) : _a(a), _b(b) {}
    /**
     *
     * @param intrinsic_permeability Passed in as intrinsic permeability
     *                               tensor, and passed out the modification of
     *                               intrinsic permeability with damage zone
     *                               model.
     * @param failure_criterion      Failure criterion.
     * @param stress                 Stress tensor.
     * @param dim                    Dimension of the intrinsic permeability
     *                               tensor.
     */
    void computeDamageZonePermeability(
        double* intrinsic_permeability,
        SolidProp::MohrCoulombFailureCriterion const& failure_criterion,
        double const* const stress,
        const int dim);

private:
    const double _a;  /// Coefficient a
    const double _b;  /// Coefficient a
};
}  // End of the namespace
