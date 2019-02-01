/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   BGRaCreep.cpp
 *
 * Created on January 21, 2019, 4:43 PM
 *
 */

#include "BGRaCreep.h"

#include <cmath>
#include <limits>

#include "PhysicalConstant.h"
#include "LinAlg/GaussAlgorithm.h"
#include "rf_msp_new.h"
#include "display.h"

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

double getCreepConstantCoefficient(const double A, const double n,
                                   const double sigma0)
{
    return A * std::pow(1.5, 0.5 * (1 + n)) / std::pow(sigma0, n);
}

void BGRaCreep::integrateStress(const double dt, const double T,
                                const CSolidProperties& solid_properties,
                                const Math_Group::Matrix& De,
                                Math_Group::Matrix& Dec,
                                double const* const stress, double* dstress,
                                const int update_s)
{
    const int nstress_components = Dec.Rows();
    const int dim = (nstress_components > 4) ? 3 : 2;

    // Assign the try stress to dstress
    double try_stress[6];
    for (int i = 0; i < nstress_components; i++)
    {
        try_stress[i] = stress[i] + dstress[i];
    }

    // Deviatoric stress of the try stress, s_{n+1}=s_n+ds
    double s[6];
    getDeviatoricStess(try_stress, nstress_components, s);
    const double norm_s = getStressNorm(s, nstress_components);

    // In case |s_{try}| is zero and _n < 3 (rare case).
    if (norm_s < std::numeric_limits<double>::epsilon() *
                     solid_properties.getYoungsModulus())
    {
        if (update_s > 0)
        {
            // Assign the try stress to dstress as the converged stress.
            for (int i = 0; i < nstress_components; i++)
            {
                dstress[i] = try_stress[i];
            }
            return;
        }

        solid_properties.ElasticConsitutive(dim, &Dec);

        return;
    }

    const double b = dt * getCreepConstantCoefficient(_a, _n, _sigma_f) *
                     std::exp(-_q / (PhysicalConstant::IdealGasConstant * T));
    const double G = solid_properties.getShearModulus();

    double solution[6];
    for (int i = 0; i < nstress_components; i++)
    {
        solution[i] = try_stress[i];
    }

    _jacobian.LimitSize(nstress_components, nstress_components);
    MathLib::GaussAlgorithm<Math_Group::Matrix> linear_solver(_jacobian,
                                                              _jacobian.Rows());

    double r[6];
    int it = 0;
    double norm_r = 0.0;
    for (it = 0; it < _max_iterations; it++)
    {
        getDeviatoricStess(solution, nstress_components, s);
        const double norm_s_n1 = getStressNorm(s, nstress_components);

        double const G2b = 2.0 * b * G;
        double const B = (_n - 1) * G2b * std::pow(norm_s_n1, _n - 3);
        // B * s \otimes s
        for (int i = 0; i < nstress_components; i++)
        {
            for (int j = 0; j < nstress_components; j++)
            {
                const double fac = (i >= 3 || j >= 3) ? 2.0 : 1.0;
                _jacobian(i, j) = fac * B * s[i] * s[j];
            }
        }

        double const A = G2b * std::pow(norm_s_n1, _n - 1);
        for (int i = 0; i < nstress_components; i++)
        {
            _jacobian(i, i) += 1.0 + A;
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                _jacobian(i, j) -= A / 3.0;
            }
        }

        // r
        norm_r = 0.0;
        for (int i = 0; i < nstress_components; i++)
        {
            r[i] = solution[i] - try_stress[i] + A * s[i];

            const double fac = (i >= 3) ? 2.0 : 1.0;
            norm_r += fac * r[i] * r[i];
        }

        if (std::sqrt(norm_r) < _tolerance)
            break;

        linear_solver.execute(r);
        // Solution of the linear solver is in r.

        // Update solutions
        for (int i = 0; i < nstress_components; i++)
        {
            solution[i] -= r[i];
        }
    }

    if (it == _max_iterations - 1)
    {
        Display::ScreenMessage(
            "The local Newton method did not converge within the given "
            "number of iterations. Iteration: %d, residual norm: "
            "%g",
            it, norm_r);
		abort();
    }

    for (int i = 0; i < nstress_components; i++)
    {
        dstress[i] = solution[i];
    }
    if (update_s > 0)
        return;

    // Compute Dec = J^{-1}_{r} De
    Dec.LimitSize(nstress_components, nstress_components);
    for (int j = 0; j < nstress_components; j++)
    {
        for (int i = 0; i < nstress_components; i++)
        {
            r[i] = De(i, j);
        }
        // the i_th column of the invJac matrix
        linear_solver.executeWithExistedElimination(r);
        for (int i = 0; i < nstress_components; i++)
            Dec(i, j) = r[i];
    }
}

}  // namespace SolidProp
