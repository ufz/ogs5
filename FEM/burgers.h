/**
 * \copyright
 * Copyright (c) 2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BURGERS_H
#define BURGERS_H

#include <cfloat>
#include "Eigen/Dense"
// FEM-Makros
#include "makros.h"
#include "tools.h"
#include "invariants.h"

namespace Burgers
{
typedef Eigen::Matrix<double, 6, 1> KVec;
typedef Eigen::Matrix<double, 6, 6> KMat;

class SolidBurgers
{
public:
    explicit SolidBurgers(const Math_Group::Matrix& data);
    ~SolidBurgers() {}
    // basic material parameters
    double GK0, GM0, KM0, etaK0, etaM0, mK, mvK, mvM, l0, T_ref, m_GM, m_KM, B,
        Q;
    // solution dependent values
    double GM, KM, GK, etaK, etaM;

    void UpdateBurgersProperties(double s_eff, const double Delta_T);
    void CalResidualBurgers(const double dt, const KVec& strain_curr,
                            const KVec& stress_curr, KVec& strain_Kel_curr,
                            const KVec& strain_Kel_t, KVec& strain_Max_curr,
                            const KVec& strain_Max_t,
                            Eigen::Matrix<double, 18, 1>& res);
    void CalJacobianBurgers(const double dt, Eigen::Matrix<double, 18, 18>& Jac,
                            const double s_eff, const KVec& sig_i,
                            const KVec& eps_K_i);
    void CaldGdEBurgers(Eigen::Matrix<double, 18, 6>& dGdE);
};
}  // namespace Burgers
#endif  // BURGERS_H
