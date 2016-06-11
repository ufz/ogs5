#ifndef MINKLEY_H
#define MINKLEY_H

#include <cfloat>
// FEM-Makros
#include "makros.h"
#include "tools.h"
#include "rf_msp_new.h"
#include "Eigen/Dense"
#include "invariants.h"

namespace Math_Group
{class Matrix;
}

namespace SolidMath
{class Invariants;
}

namespace Minkley{
using Math_Group::Matrix;
using SolidMath::Invariants;
class SolidMinkley
{
public:
    SolidMinkley(const Matrix* data);
    ~SolidMinkley();
    //basic material parameters
    double GK0, GM0, KM0, etaK0, etaM0, mvM;
    //solution dependent values
    double etaM;
    //Parameters for Minkley model
    double nvM, coh0, hard, phi, psi, thetaT, eta_reg;
    //solution dependent
    double coh;
    double l0; // temperature parameter for Maxwell viscosity
    double T0; // reference temperature for Maxwell viscosity

    Invariants* smath;

    void UpdateMinkleyProperties(double s_eff, const double eps_p_eff, double Temperature);
    double YieldMohrCoulomb(const Eigen::Matrix<double,6,1> &sig);
    void CalViscoelasticResidual(const double dt, const Eigen::Matrix<double,6,1> &dstrain_curr, const double e_curr, const double e_p_curr,
                                 const Eigen::Matrix<double,6,1> &stress_curr, const Eigen::Matrix<double,6,1> &dstrain_Kel_curr,
                                 const Eigen::Matrix<double,6,1> &dstrain_Kel_t, const Eigen::Matrix<double,6,1> &dstrain_Max_curr,
                                 const Eigen::Matrix<double,6,1> &dstrain_Max_t, const Eigen::Matrix<double,6,1> &dstrain_p_curr, Eigen::Matrix<double,18,1> &res);
    void CalViscoelasticJacobian(const double dt, const Eigen::Matrix<double,6,1> &stress_curr, const double sig_eff,Eigen::Matrix<double,18,18> &Jac);
    void CaldGdE(Eigen::Matrix<double,18,6> &dGdE);
    void CalViscoplasticResidual(const double dt, const Eigen::Matrix<double,6,1> &dstrain_curr, const double e_curr,
                                 const Eigen::Matrix<double,6,1> &stress_curr, const Eigen::Matrix<double,6,1> &dstrain_Kel_curr,
                                 const Eigen::Matrix<double,6,1> &dstrain_Kel_t, const Eigen::Matrix<double,6,1> &dstrain_Max_curr,
                                 const Eigen::Matrix<double,6,1> &dstrain_Max_t, const Eigen::Matrix<double,6,1> &dstrain_pl_curr,
                                 const Eigen::Matrix<double,6,1> &dstrain_pl_t, const double e_pl_vol_curr, const double e_pl_vol_t,
                                 const double e_pl_eff_curr, const double e_pl_eff_t, const double lam_curr, Eigen::Matrix<double,27,1> &res);
    void CalViscoplasticJacobian(const double dt, const Eigen::Matrix<double,6,1> &stress_curr, const double sig_eff,
                                               const double lam_curr, Eigen::Matrix<double,27,27> &Jac);
    void CalEPdGdE(Eigen::Matrix<double,27,6> &dGdE);
    void NumericalJacobian(const double dt, const Eigen::Matrix<double,6,1> &dstrain_curr, const double e_curr,
                                         const Eigen::Matrix<double,6,1> &stress_curr, const Eigen::Matrix<double,6,1> &dstrain_Kel_curr,
                                         const Eigen::Matrix<double,6,1> &dstrain_Kel_t, const Eigen::Matrix<double,6,1> &dstrain_Max_curr,
                                         const Eigen::Matrix<double,6,1> &dstrain_Max_t, const Eigen::Matrix<double,6,1> &dstrain_pl_curr,
                                         const Eigen::Matrix<double,6,1> &dstrain_pl_t, const double e_pl_vol_curr, const double e_pl_vol_t,
                                         const double e_pl_eff_curr, const double e_pl_eff_t, const double lam_curr, Eigen::Matrix<double,27,27> &Jac_num);

private:
    double A(const double theta, const double alpha);
    double B(const double theta, const double alpha);
    Eigen::Matrix<double,6,1> DetaM_Dsigma(double sig_eff, const Eigen::Matrix<double,6,1> &sigd_i);
    double DG_DI1(const double alpha);
    double DG_DJ2(const double theta, const double J2, const double alpha);
    double DG_Dtheta(const double theta, const double J2, const double alpha);
    double Dtheta_DJ2(const double theta, const double J2);
    double Dtheta_DJ3(const double theta, const double J3);
    Eigen::Matrix<double,6,6> s_odot_s(const Eigen::Matrix<double,6,1> &vec);
    double DDG_DDJ2(const double theta, const double J2, const double alpha);
    double DDG_DJ2_Dtheta(const double theta, const double J2, const double alpha);
    double DDG_DDtheta(const double theta, const double J2, const double alpha);
    double DDtheta_DDJ2(const double theta, const double J2);
    double DDtheta_DJ2_Dtheta(const double theta, const double J2);
    double DDtheta_DJ3_Dtheta(const double theta, const double J3);
    double DDtheta_DDJ3(const double theta, const double J3);
};
}
#endif // MINKLEY_H
