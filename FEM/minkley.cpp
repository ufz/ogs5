/**
 * \copyright
 * Copyright (c) 2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "minkley.h"
#include "PhysicalConstant.h"

namespace Minkley
{
// Template to implement Signum Function
template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

SolidMinkley::SolidMinkley(const Math_Group::Matrix& data)
{
    GK0 = data(0);    // Kelvin shear modulus
    etaK0 = data(1);  // Kelvin viscosity
    GM0 = data(2);    // Maxwell shear modulus
    KM0 = data(3);    // Maxwell bulk modulus
    etaM0 = data(4);  // Maxwell viscosity
    if (std::abs(data(5)) < 1.e-3)
    {
        mvM = std::log(
            1. +
            std::sqrt(
                2.));  // m -- effective stress factor in viscosity relationship
        nvM = 0.;  // n -- effective stress exponent in viscosity relationship
        std::cout
            << "WARNING. Viscosity parameter m was chosen lower than 1e-3. "
               "Setting model parameters to achieve constant viscosity.";
    }
    else
    {
        mvM =
            data(5);  // m -- effective stress factor in viscosity relationship
        nvM = data(
            6);  // n -- effective stress exponent in viscosity relationship
    }
    coh0 = data(7);                 // initial cohesion
    hard = data(8);                 // hardening/softening modulus
    phi = data(9) * PI / 180.;      // friction angle
    psi = data(10) * PI / 180.;     // dilatancy angle
    thetaT = data(11) * PI / 180.;  // transition angle
    eta_reg = data(12);             // viscosity for viscoplastic regularisation
    m_GM = data(13);   // slope of elesticity temperature dependence
    m_KM = data(14);   // slope of elesticity temperature dependence
    T_ref = data(15);  // reference temperature dependency parameter for "
    Bt = data(16);     // constant factor for Arrhenius term
    Q = data(17);      // activation energy in Arrhenius term
    hard2 = data(18);  // second order hardening term
    hard4 = data(19);  // second order hardening term

    etaM = etaM0;
    coh = coh0;
    GM = GM0;
    KM = KM0;

    if (!T_Process)
    {
        m_GM = 0.;
        m_KM = 0.;
        Bt = 1.;
        Q = 0.;  // for cutting off Arrhenius term
        T_ref = PhysicalConstant::CelsiusZeroInKelvin;
    }
}

/**************************************************************************
   FEMLib-Method: Minkley::UpdateMinkleyProperties()
   Task: Updates BURGERS material parameters in LUBBY2 fashion
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void SolidMinkley::UpdateMinkleyProperties(double s_eff, const double eps_p_eff,
                                           double Temperature)
{
    const double dT(Temperature - T_ref);
    GM = GM0 + m_GM * dT;
    KM = KM0 + m_KM * dT;

    s_eff *= GM;

    if (s_eff > DBL_EPSILON)
        etaM =
            etaM0 /
            std::sinh(mvM * std::pow(s_eff, nvM));  // viscosity function update
    else
        etaM = etaM0;
    etaM *=
        Bt *
        std::exp(Q * (-dT) /
                 (PhysicalConstant::IdealGasConstant * Temperature * T_ref));

    coh =
        coh0 *
        (1. +
         eps_p_eff *
             (hard +
              eps_p_eff *
                  (hard2 +
                   eps_p_eff * eps_p_eff *
                       hard4)));  // fourth order isotropic hardening/softening
                                  //	if (etaM / etaM0 < 1.e-2)
    //		std::cout << "WARNING: Maxwell viscosity sank to 100th of original
    //value." << std::endl;
}

/**************************************************************************
   FEMLib-Method: Minkley::A()
   Task: Expression A in Sloan's yield function
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::A(const double theta, const double alpha)
{
    double A;
    A = std::cos(thetaT) / 3. *
        (3. + std::tan(thetaT) * std::tan(3. * thetaT) +
         1. / std::sqrt(3.) * sgn(theta) *
             (std::tan(3. * thetaT) - 3. * std::tan(thetaT)) * std::sin(alpha));
    return A;
}

/**************************************************************************
   FEMLib-Method: Minkley::B()
   Task: Expression B in Sloan's yield function
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::B(const double theta, const double alpha)
{
    double B;
    B = 1. / (3. * std::cos(3. * thetaT)) *
        (sgn(theta) * std::sin(thetaT) +
         std::sin(alpha) * std::cos(thetaT) / std::sqrt(3.));
    return B;
}

/**************************************************************************
   FEMLib-Method: Minkley::YieldMohrCoulomb()
   Task: Yield function Mohr Coulomb with corner smoothing (Sloan et al.)
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::YieldMohrCoulomb(const KVec& sig)
{
    double F(0.);
    const KVec sigd(SolidMath::P_dev * sig);
    const double theta(SolidMath::CalLodeAngle(sigd));

    if (std::abs(theta) < thetaT)
    {
        F = SolidMath::CalI1(sig) / 3. * std::sin(phi) +
            std::sqrt(SolidMath::CalJ2(sigd)) *
                (std::cos(theta) -
                 1. / std::sqrt(3.) * std::sin(phi) * std::sin(theta)) -
            coh * std::cos(phi);
    }
    else
        F = SolidMath::CalI1(sig) / 3. * std::sin(phi) +
            std::sqrt(SolidMath::CalJ2(sigd)) *
                (A(theta, phi) - B(theta, phi) * std::sin(3. * theta)) -
            coh * std::cos(phi);

    return F;
}

/**************************************************************************
   FEMLib-Method: Minkley::DetaM_Dsigma()
   Task: Derivative of Maxwell viscosity with respect to normalised deviatoric
stress Programing: 06/2015 TN Implementation
**************************************************************************/
KVec SolidMinkley::DetaM_Dsigma(double sig_eff, const KVec& sigd_i)
{
    KVec res;
    if (sig_eff < DBL_EPSILON || std::abs(nvM) < DBL_EPSILON)
        return sigd_i * 0.;
    else
    {
        res = 3. / 2. * sigd_i * GM;  // sig_eff in denominator lumped into pow
                                      // function (-2 instead of -1)
        res *= -etaM * mvM * nvM *
               std::pow(sig_eff, nvM - 2.)  // etaM = etaM0/sinh(..)
               / std::tanh(mvM * std::pow(sig_eff, nvM));
        return res;
    }
}

/**************************************************************************
   FEMLib-Method: Minkley::CalViscoelasticResidual()
   Task: Calculates the 12x1 residual vector. Implementation fully implicit
only. Programing: 06/2015 TN Implementation
**************************************************************************/
void SolidMinkley::CalViscoelasticResidual(
    const double dt, const KVec& dstrain_curr, const double e_curr,
    const double e_p_curr, const KVec& stress_curr,
    const KVec& dstrain_Kel_curr, const KVec& dstrain_Kel_t,
    const KVec& dstrain_Max_curr, const KVec& dstrain_Max_t,
    const KVec& dstrain_p_curr, Eigen::Matrix<double, 18, 1>& res)
{
    const KVec dstress_curr(GM * SolidMath::P_dev * stress_curr);

    // calculate stress residual
    res.block<6, 1>(0, 0) =
        stress_curr - (2. * (dstrain_curr - dstrain_Kel_curr -
                             dstrain_Max_curr - dstrain_p_curr) +
                       KM / GM * (e_curr - e_p_curr) * SolidMath::ivec);
    // calculate Kelvin strain residual
    res.block<6, 1>(6, 0) =
        (dstrain_Kel_curr - dstrain_Kel_t) / dt -
        (dstress_curr - 2. * GK0 * dstrain_Kel_curr) / (2. * etaK0);
    // calculate Kelvin strain residual
    res.block<6, 1>(12, 0) =
        (dstrain_Max_curr - dstrain_Max_t) / dt - dstress_curr / (2. * etaM);
}

/**************************************************************************
   FEMLib-Method: Minkley::CalViscoelasticJacobian()
   Task: Calculates the 12x2 Jacobian matrix. Implementation fully implicit
only. Programing: 06/2015 TN Implementation
**************************************************************************/
void SolidMinkley::CalViscoelasticJacobian(const double dt,
                                           const KVec& stress_curr,
                                           const double sig_eff,
                                           Eigen::Matrix<double, 18, 18>& Jac)
{
    // 6x6 submatrices of the Jacobian
    const KVec sigd_curr(GM * SolidMath::P_dev * stress_curr);
    const KVec dmu_vM = DetaM_Dsigma(sig_eff * GM, sigd_curr);

    // Check Dimension of Jacobian
    // assert(Jac.cols() == 18 && Jac.rows() == 18);
    Jac.setZero(18, 18);

    // build G_11
    Jac.block<6, 6>(0, 0) = SolidMath::ident;

    // build G_12
    Jac.block<6, 6>(0, 6) = 2. * SolidMath::ident;

    // build G_13
    Jac.block<6, 6>(0, 12) = 2. * SolidMath::ident;

    // build G_21
    Jac.block<6, 6>(6, 0) = -GM / (2. * etaK0) * SolidMath::P_dev;

    // build G_22
    Jac.block<6, 6>(6, 6) = (1. / dt + GK0 / etaK0) * SolidMath::ident;

    // nothing to do for G_23

    // build G_31
    Jac.block<6, 6>(12, 0) =
        -1. / (2. * etaM) *
        (GM * SolidMath::P_dev - sigd_curr * dmu_vM.transpose() / etaM);

    // nothing to do for G_32

    // build G_33
    Jac.block<6, 6>(12, 12) = 1. / dt * SolidMath::ident;
}

/**************************************************************************
   FEMLib-Method: Burgers::CaldGdE()
   Task: Calculates the 12x6 derivative of the residuals with respect to total
strain. Implementation fully implicit only. Programing: 06/2015 TN
Implementation
**************************************************************************/
void SolidMinkley::CaldGdE(Eigen::Matrix<double, 18, 6>& dGdE)
{
    // Check Dimension of dGdE
    // assert(dGdE.cols() == 6 && dGdE.rows() == 18);
    dGdE.setZero(18, 6);
    dGdE.block<6, 6>(0, 0) =
        -2. * SolidMath::P_dev - 3. * KM / GM * SolidMath::P_sph;
}

/**************************************************************************
   FEMLib-Method: Minkley::DG_DI1()
   Task: \partial G / \partial I_1
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DG_DI1(const double alpha)
{
    return std::sin(alpha) / 3.;
}

/**************************************************************************
   FEMLib-Method: Minkley::DG_DJ2()
   Task: \partial G / \partial J2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DG_DJ2(const double theta, const double J2,
                            const double alpha)
{
    if (std::abs(theta) < thetaT)
        return (std::cos(theta) -
                std::sin(alpha) * std::sin(theta) / std::sqrt(3.)) /
               (2. * std::sqrt(J2));
    else
        return (A(theta, alpha) - B(theta, alpha) * std::sin(3. * theta)) /
               (2. * std::sqrt(J2));
}

/**************************************************************************
   FEMLib-Method: Minkley::DG_Dtheta()
   Task: \partial G / \partial theta
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DG_Dtheta(const double theta, const double J2,
                               const double alpha)
{
    if (std::abs(theta) < thetaT)
        return -std::sqrt(J2) *
               (std::sin(theta) +
                std::sin(alpha) * std::cos(theta) / std::sqrt(3.));
    else
        return -std::sqrt(J2) * 3. * B(theta, alpha) * std::cos(3. * theta);
}

/**************************************************************************
   FEMLib-Method: Minkley::Dtheta_DJ2()
   Task: \partial theta / \partial J_2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::Dtheta_DJ2(const double theta, const double J2)
{
    return -std::tan(3. * theta) / (2. * J2);
}

/**************************************************************************
   FEMLib-Method: Minkley::Dtheta_DJ3()
   Task: \partial theta / \partial J_3
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::Dtheta_DJ3(const double theta, const double J3)
{
    return std::tan(3. * theta) / (3. * J3);
}

/**************************************************************************
   FEMLib-Method: Minkley::CalViscoplasticResidual()
   Task: Calculates the 21x1 residual vector. Implementation fully implicit
only. Programing: 06/2015 TN Implementation
**************************************************************************/
void SolidMinkley::CalViscoplasticResidual(
    const double dt, const KVec& dstrain_curr, const double e_curr,
    const KVec& stress_curr, const KVec& dstrain_Kel_curr,
    const KVec& dstrain_Kel_t, const KVec& dstrain_Max_curr,
    const KVec& dstrain_Max_t, const KVec& dstrain_pl_curr,
    const KVec& dstrain_pl_t, const double e_pl_vol_curr,
    const double e_pl_vol_t, const double e_pl_eff_curr,
    const double e_pl_eff_t, const double lam_curr,
    Eigen::Matrix<double, 27, 1>& res)
{
    const KVec sigd_curr(GM * SolidMath::P_dev * stress_curr);
    const double J_2(SolidMath::CalJ2(sigd_curr)),
        J_3(SolidMath::CalJ3(sigd_curr)),
        theta(SolidMath::CalLodeAngle(sigd_curr));
    const KVec dev_sigd_curr_inv(SolidMath::P_dev *
                                 SolidMath::InvertVector(sigd_curr));
    const double vol_flow(3. * DG_DI1(psi));

    KVec dev_flow;
    if (std::abs(J_3) > 0.)
        dev_flow = (DG_DJ2(theta, J_2, psi) +
                    DG_Dtheta(theta, J_2, psi) * Dtheta_DJ2(theta, J_2)) *
                       sigd_curr +
                   (DG_Dtheta(theta, J_2, psi) * Dtheta_DJ3(theta, J_3) * J_3) *
                       dev_sigd_curr_inv;
    else
        dev_flow.setZero(6);

    // calculate stress residual
    res.block<6, 1>(0, 0) =
        stress_curr - (2. * (dstrain_curr - dstrain_Kel_curr -
                             dstrain_Max_curr - dstrain_pl_curr) +
                       KM / GM * (e_curr - e_pl_vol_curr) * SolidMath::ivec);

    // calculate deviatoric Kelvin strain residual
    res.block<6, 1>(6, 0) =
        (dstrain_Kel_curr - dstrain_Kel_t) / dt -
        (sigd_curr - 2. * GK0 * dstrain_Kel_curr) / (2. * etaK0);

    // calculate deviatoric Maxwell strain residual
    res.block<6, 1>(12, 0) =
        (dstrain_Max_curr - dstrain_Max_t) / dt - sigd_curr / (2. * etaM);

    // calculate deviatoric plastic strain residual
    res.block<6, 1>(18, 0) =
        (dstrain_pl_curr - dstrain_pl_t) / dt - lam_curr * dev_flow;

    // calculate volumetric plastic strain residual
    res.block<1, 1>(24, 0)(0) =
        (e_pl_vol_curr - e_pl_vol_t) / dt - lam_curr * vol_flow;

    // calculate effective plastic strain residual
    res.block<1, 1>(25, 0)(0) =
        (e_pl_eff_curr - e_pl_eff_t) / dt -
        std::sqrt(2. / 3. * lam_curr * lam_curr *
                  (double)(dev_flow.transpose() * dev_flow));

    // yield function with viscoplastic regularisation
    res.block<1, 1>(26, 0)(0) =
        YieldMohrCoulomb(stress_curr * GM) / GM - lam_curr * eta_reg;
}

/**************************************************************************
   FEMLib-Method: Minkley::s_odot_s()
   Task: vec \odot vec (all in Kelvin mapping)
   Programing:
   06/2015 TN Implementation
**************************************************************************/
/*Kelvin mapping of odot changed according to symbolic conversion.
Note: the factors of 2 and sqrt(2) come from the fact that the
incoming quantities are transformed to actual tensor coordinates
and the resulting 4th order tensor is then remapped with the
Kelvin scheme*/
KMat SolidMinkley::s_odot_s(const KVec& vec)
{
    KMat odot;

    odot(0, 0) = -vec(0) * vec(0);
    odot(0, 1) = odot(1, 0) = -vec(3) * vec(3) / 2.;
    odot(0, 2) = odot(2, 0) = -vec(5) * vec(5) / 2.;
    odot(0, 3) = odot(3, 0) = -vec(0) * vec(3);
    odot(0, 4) = odot(4, 0) = -vec(3) * vec(5) / std::sqrt(2.);
    odot(0, 5) = odot(5, 0) = -vec(0) * vec(5);

    odot(1, 1) = -vec(1) * vec(1);
    odot(1, 2) = odot(2, 1) = -vec(4) * vec(4) / 2.;
    odot(1, 3) = odot(3, 1) = -vec(3) * vec(1);
    odot(1, 4) = odot(4, 1) = -vec(1) * vec(4);
    odot(1, 5) = odot(5, 1) = -vec(3) * vec(4) / std::sqrt(2.);

    odot(2, 2) = -vec(2) * vec(2);
    odot(2, 3) = odot(3, 2) = -vec(5) * vec(4) / std::sqrt(2.);
    odot(2, 4) = odot(4, 2) = -vec(4) * vec(2);
    odot(2, 5) = odot(5, 2) = -vec(5) * vec(2);

    odot(3, 3) = -vec(0) * vec(1) - vec(3) * vec(3) / 2.;
    odot(3, 4) = odot(4, 3) =
        -vec(3) * vec(4) / 2. - vec(5) * vec(1) / std::sqrt(2.);
    odot(3, 5) = odot(5, 3) =
        -vec(0) * vec(4) / std::sqrt(2.) - vec(3) * vec(5) / 2.;

    odot(4, 4) = -vec(1) * vec(2) - vec(4) * vec(4) / 2.;
    odot(4, 5) = odot(5, 4) =
        -vec(3) * vec(2) / std::sqrt(2.) - vec(5) * vec(4) / 2.;

    odot(5, 5) = -vec(0) * vec(2) - vec(5) * vec(5) / 2.;
    return odot;
}

/**************************************************************************
   FEMLib-Method: Minkley::DDG_DDJ2()
   Task: \partial^2 G / \partial J_2^2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDG_DDJ2(const double theta, const double J2,
                              const double alpha)
{
    if (std::abs(theta) < thetaT)
        return (std::cos(theta) -
                std::sin(alpha) * std::sin(theta) / std::sqrt(3.)) /
               (-4. * std::pow(J2, 1.5));
    else
        return (A(theta, alpha) - B(theta, alpha) * std::sin(3. * theta)) /
               (-4. * std::pow(J2, 1.5));
}

/**************************************************************************
   FEMLib-Method: Minkley::DDG_DJ2_Dtheta()
   Task: \partial^2 G / (\partial J_2 \partial theta)
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDG_DJ2_Dtheta(const double theta, const double J2,
                                    const double alpha)
{
    if (std::abs(theta) < thetaT)
        return (std::sin(theta) +
                std::sin(alpha) * std::cos(theta) / std::sqrt(3.)) /
               (-2. * std::sqrt(J2));
    else
        return 3. * B(theta, alpha) * std::cos(3. * theta) /
               (-2. * std::sqrt(J2));
}

/**************************************************************************
   FEMLib-Method: Minkley::DDG_DDtheta()
   Task: \partial^2 G / \partial theta^2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDG_DDtheta(const double theta, const double J2,
                                 const double alpha)
{
    if (std::abs(theta) < thetaT)
        return (std::cos(theta) -
                std::sin(alpha) * std::sin(theta) / std::sqrt(3.)) *
               (-std::sqrt(J2));
    else
        return 9. * B(theta, alpha) * std::sin(3. * theta) * std::sqrt(J2);
}

/**************************************************************************
   FEMLib-Method: Minkley::DDtheta_DDJ2()
   Task: \partial^2 theta / \partial J_2^2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDtheta_DDJ2(const double theta, const double J2)
{
    return std::tan(3. * theta) / (2. * J2 * J2);
}

/**************************************************************************
   FEMLib-Method: Minkley::DDtheta_DJ2_Dtheta()
   Task: \partial^2 theta / (\partial J_2 \partial theta)
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDtheta_DJ2_Dtheta(const double theta, const double J2)
{
    return -3. / (2. * J2 * std::pow(std::cos(3. * theta), 2.));
}

/**************************************************************************
   FEMLib-Method: Minkley::DDtheta_DJ3_Dtheta()
   Task: \partial^2 theta / (\partial J_3 \partial theta)
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDtheta_DJ3_Dtheta(const double theta, const double J3)
{
    return 1. / (J3 * std::pow(std::cos(3. * theta), 2.));
}

/**************************************************************************
   FEMLib-Method: Minkley::DDtheta_DDJ3()
   Task: \partial^2 theta / \partial J_3^2
   Programing:
   06/2015 TN Implementation
**************************************************************************/
double SolidMinkley::DDtheta_DDJ3(const double theta, const double J3)
{
    return -std::tan(3. * theta) / (3. * J3 * J3);
}

/**************************************************************************
   FEMLib-Method: Minkley::CalViscoplasticJacobian()
   Task: Calculates the 27x27 Jacobian matrix. Implementation fully implicit
only. Programing: 06/2015 TN Implementation
**************************************************************************/
void SolidMinkley::CalViscoplasticJacobian(const double dt,
                                           const KVec& stress_curr,
                                           const double sig_eff,
                                           const double lam_curr,
                                           const double e_eff_i,
                                           Eigen::Matrix<double, 27, 27>& Jac)
{
    // submatrices of the Jacobian
    const KVec sigd_curr(GM * SolidMath::P_dev * stress_curr);
    const double J_2(SolidMath::CalJ2(sigd_curr)),
        J_3(SolidMath::CalJ3(sigd_curr)),
        theta(SolidMath::CalLodeAngle(sigd_curr));
    const KVec sigd_curr_inv(SolidMath::InvertVector(sigd_curr));
    const KVec dev_sigd_curr_inv(SolidMath::P_dev * sigd_curr_inv);
    const double vol_flow(3. * DG_DI1(psi));
    const KVec dmu_vM = DetaM_Dsigma(sig_eff * GM, sigd_curr);
    KVec dev_flow;
    KMat Ddev_flowDsigma;
    const double DthetaDJ2(Dtheta_DJ2(theta, J_2));
    const double DthetaDJ3(Dtheta_DJ3(theta, J_3));

    // Check Dimension of Jacobian
    // assert(Jac.cols() == 27 && Jac.rows() == 27);
    Jac.setZero(27, 27);

    if (std::abs(J_3) < 0.)
    {
        dev_flow.Zero(6, 1);
        Ddev_flowDsigma.Zero(6, 6);
    }
    else
    {
        const double DGDtheta(DG_Dtheta(theta, J_2, psi));
        dev_flow =
            (DG_DJ2(theta, J_2, psi) + DGDtheta * DthetaDJ2) * sigd_curr +
            (DGDtheta * DthetaDJ3 * J_3) * dev_sigd_curr_inv;
        Ddev_flowDsigma =
            ((DG_DJ2(theta, J_2, psi) + DGDtheta * DthetaDJ2) *
                 SolidMath::P_dev +
             (DGDtheta * DthetaDJ3 * J_3) * SolidMath::P_dev *
                 s_odot_s(sigd_curr_inv) * SolidMath::P_dev +
             (DDG_DDJ2(theta, J_2, psi) +
              DDG_DJ2_Dtheta(theta, J_2, psi) * DthetaDJ2 +
              DGDtheta * DDtheta_DDJ2(theta, J_2)) *
                 sigd_curr * sigd_curr.transpose() +
             (DDG_DJ2_Dtheta(theta, J_2, psi) +
              DDG_DDtheta(theta, J_2, psi) * DthetaDJ2 +
              DGDtheta * DDtheta_DJ2_Dtheta(theta, J_2)) *
                 (DthetaDJ2 * sigd_curr * sigd_curr.transpose() +
                  DthetaDJ3 * J_3 * sigd_curr * dev_sigd_curr_inv.transpose()) +
             DDG_DJ2_Dtheta(theta, J_2, psi) * DthetaDJ3 * J_3 *
                 dev_sigd_curr_inv * sigd_curr.transpose() +
             DGDtheta * J_3 * (DthetaDJ3 + DDtheta_DDJ3(theta, J_3) * J_3) *
                 dev_sigd_curr_inv * dev_sigd_curr_inv.transpose() +
             J_3 *
                 (DDG_DDtheta(theta, J_2, psi) * DthetaDJ3 +
                  DGDtheta * DDtheta_DJ3_Dtheta(theta, J_3)) *
                 (DthetaDJ2 * dev_sigd_curr_inv * sigd_curr.transpose() +
                  DthetaDJ3 * J_3 * dev_sigd_curr_inv *
                      dev_sigd_curr_inv.transpose())) *
            GM;
    }

    const double eff_flow =
        std::sqrt(2. / 3. * lam_curr * lam_curr *
                  (double)(dev_flow.transpose() * dev_flow));

    // build G_11
    // G_66 = SolidMath::ident/dt +
    Jac.block<6, 6>(0, 0) = SolidMath::ident;

    // build G_12, G_13 and G_14
    Jac.block<6, 6>(0, 6) = 2. * SolidMath::ident;
    Jac.block<6, 6>(0, 12) = 2. * SolidMath::ident;
    Jac.block<6, 6>(0, 18) = 2. * SolidMath::ident;

    // build G_15
    Jac.block<6, 1>(0, 24) = KM / GM * SolidMath::ivec;

    // G_16 and G_17 remain zeros and are not set separately

    // build G_21
    Jac.block<6, 6>(6, 0) = -GM / (2. * etaK0) * SolidMath::P_dev;

    // build G_22
    Jac.block<6, 6>(6, 6) = (1. / dt + GK0 / etaK0) * SolidMath::ident;

    // G_23 through G_27 are zero

    // build G_31
    Jac.block<6, 6>(12, 0) =
        -1. / (2. * etaM) *
        (GM * SolidMath::P_dev - sigd_curr * dmu_vM.transpose() / etaM);

    // G_32 is 0

    // G_33
    Jac.block<6, 6>(12, 12) = SolidMath::ident / dt;

    // G_34 through G_37 are zero

    // G_41
    Jac.block<6, 6>(18, 0) = -lam_curr * Ddev_flowDsigma;

    // build G_42 and G_43 are zero

    // build G_44
    Jac.block<6, 6>(18, 18) = SolidMath::ident / dt;

    // G_45 and G_46 are zero

    // build G_47
    Jac.block<6, 1>(18, 26) = -dev_flow;

    // G_51 to G_54 are zero

    // build G_55
    Jac.block<1, 1>(24, 24)(0) = 1. / dt;

    // G_56 is zero

    // G_57
    Jac.block<1, 1>(24, 26)(0) = -vol_flow;

    // Conditional build of G_61 und G_67
    if (std::abs(eff_flow) > 0.)
    {
        Jac.block<1, 6>(25, 0) = -2. * lam_curr * lam_curr / (3. * eff_flow) *
                                 dev_flow.transpose() *
                                 Ddev_flowDsigma.transpose();
        Jac.block<1, 1>(25, 26)(0) = -2. / (3. * eff_flow) *
                                     std::abs(lam_curr) *
                                     (double)(dev_flow.transpose() * dev_flow);
    }
    // G_62 to G_64 and G_65 are zero

    // build G_66
    Jac.block<1, 1>(25, 25)(0) = 1. / dt;

    // Yield surface derivatives
    if (std::abs(J_3) > 0.)
        dev_flow =
            (DG_DJ2(theta, J_2, phi) + DG_Dtheta(theta, J_2, phi) * DthetaDJ2) *
                sigd_curr +
            (DG_Dtheta(theta, J_2, phi) * DthetaDJ3 * J_3) * dev_sigd_curr_inv;

    else
        dev_flow.Zero(6, 1);

    // build G_71
    Jac.block<1, 6>(26, 0) =
        (DG_DI1(phi) * SolidMath::ivec + dev_flow).transpose();

    // G_72 - G_75 zero

    // build G_76
    Jac.block<1, 1>(26, 25)(0) =
        -coh0 *
        (hard + 2. * e_eff_i * (hard2 + 2. * e_eff_i * e_eff_i * hard4)) *
        std::cos(phi) / GM;

    // build G_77
    Jac.block<1, 1>(26, 26)(0) = -eta_reg;
}

/**************************************************************************
   FEMLib-Method: Burgers::CalEPdGdE()
   Task: Calculates the 27x6 derivative of the residuals with respect to total
strain. Implementation fully implicit only. Programing: 06/2015 TN
Implementation
**************************************************************************/
void SolidMinkley::CalEPdGdE(Eigen::Matrix<double, 27, 6>& dGdE)
{
    // Check Dimension of dGdE
    // assert(dGdE.cols() == 6 && dGdE.rows() == 27);
    dGdE.setZero(27, 6);
    dGdE.block<6, 6>(0, 0) =
        -2. * SolidMath::P_dev - 3. * KM / GM * SolidMath::P_sph;
}

}  // namespace Minkley
