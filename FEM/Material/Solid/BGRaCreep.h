/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   BGRaCreep.h
 *
 * Created on January 21, 2019, 4:43 PM
 *
 */

#pragma once

#include "matrix_class.h"

namespace SolidProp
{
class CSolidProperties;

/**
 * \brief A class for computing the BGRa creep model, which reads
 * \f[
 * \dot {\mathbf{\epsilon}}^{cr}=\sqrt{\frac{3}{2}}A \mathrm{e}^{-\frac{Q}{RT}}
 * \left(\frac{\sigma_{eff}}{\sigma_0}\right)^n\frac{\mathbf{s}}{||\mathbf{s}||}
 * \f]
 * where \f$\sigma_{eff}=\sqrt{\frac{3}{2}}||\mathbf{s}||\f$, \f$A, \sigma_0, n,
 * Q\f$ are parameter, and \f$R\f$ is the gas constant.
 *
 */
class BGRaCreep
{
public:
    BGRaCreep(const double A, const double n, const double sigma_f,
              const double Q, const double tolerance, const int max_iterations)
        : _a(A),
          _n(n),
          _sigma_f(sigma_f),
          _q(Q),
          _tolerance(tolerance),
          _max_iterations(max_iterations),
          _jacobian(Math_Group::Matrix(6, 6))
    {
    }

    ~BGRaCreep() {}
    /**
     *
     * @param dt        Time increment.
     * @param T         Temperature.
     * @param solid_properties The solid properties, which contains G and mu.
     * @param De        The elastic tensor.
     * @param Dec       The elastic creep tensor to be calculated. Assuming that
     *                  it is initialized as De.
     * @param stress    Stress of the previous time step.
     * @param dstress   Stress increment. It passes out the integrated stress.
     * @param update_s  An indicator to indicate whether the stress integration
     *                  is for the updating of stress.
     */
    void integrateStress(const double dt, const double T,
                         const CSolidProperties& solid_properties,
                         const Math_Group::Matrix& De, Math_Group::Matrix& Dec,
                         double const* const stress, double* dstress,
                         const int update_s);

private:
    const double _a;        /// A parameter determined by experiment.
    const double _n;        /// Creep rate exponent n.
    const double _sigma_f;  /// A stress scaling factor.
    const double _q;        /// Activation energy

    /// Tolerance for the convergence of nonlinear iterations.
    const double _tolerance;
    const int _max_iterations;  /// Maximum nonlinear iterations.

    // Not for multi-threads
    Math_Group::Matrix _jacobian;
};
}  // end of namspace
