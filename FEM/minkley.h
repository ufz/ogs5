/**
 * \copyright
 * Copyright (c) 2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MINKLEY_H
#define MINKLEY_H

#include <cfloat>
#include "Eigen/Dense"
// FEM-Makros
#include "makros.h"
#include "tools.h"
#include "rf_msp_new.h"
#include "invariants.h"
#include "PhysicalConstant.h"

namespace Minkley
{
typedef Eigen::Matrix<double, 6, 1> KVec;
typedef Eigen::Matrix<double, 6, 6> KMat;

class SolidMinkley
{
public:
	explicit SolidMinkley(const Math_Group::Matrix& data);
	~SolidMinkley() {}
	// basic material parameters
	double GK0, GM0, KM0, etaK0, etaM0, mvM;
	// solution dependent values
	double etaM, GM, KM;
	// Parameters for Minkley model
	double nvM, coh0, hard, phi, psi, thetaT, eta_reg, hard2, hard4;
	// solution dependent
	double coh;
	double m_GM; // slope of elesticity temperature dependence
	double m_KM; // slope of elesticity temperature dependence
	double T_ref; // reference temperature dependency parameter for "
	double Bt; // constant factor for Arrhenius term
	double Q; // activation energy in Arrhenius term

	void UpdateMinkleyProperties(double s_eff, const double eps_p_eff, double Temperature);
	double YieldMohrCoulomb(const KVec& sig);
	void CalViscoelasticResidual(const double dt, const KVec& dstrain_curr, const double e_curr, const double e_p_curr,
	                             const KVec& stress_curr, const KVec& dstrain_Kel_curr, const KVec& dstrain_Kel_t,
	                             const KVec& dstrain_Max_curr, const KVec& dstrain_Max_t, const KVec& dstrain_p_curr,
	                             Eigen::Matrix<double, 18, 1>& res);
	void CalViscoelasticJacobian(const double dt, const KVec& stress_curr, const double sig_eff,
	                             Eigen::Matrix<double, 18, 18>& Jac);
	void CaldGdE(Eigen::Matrix<double, 18, 6>& dGdE);
	void CalViscoplasticResidual(const double dt, const KVec& dstrain_curr, const double e_curr,
	                             const KVec& stress_curr, const KVec& dstrain_Kel_curr, const KVec& dstrain_Kel_t,
	                             const KVec& dstrain_Max_curr, const KVec& dstrain_Max_t, const KVec& dstrain_pl_curr,
	                             const KVec& dstrain_pl_t, const double e_pl_vol_curr, const double e_pl_vol_t,
	                             const double e_pl_eff_curr, const double e_pl_eff_t, const double lam_curr,
	                             Eigen::Matrix<double, 27, 1>& res);
	void CalViscoplasticJacobian(const double dt, const KVec& stress_curr, const double sig_eff, const double lam_curr,
	                             const double e_eff_i, Eigen::Matrix<double, 27, 27>& Jac);
	void CalEPdGdE(Eigen::Matrix<double, 27, 6>& dGdE);

private:
	double A(const double theta, const double alpha);
	double B(const double theta, const double alpha);
	KVec DetaM_Dsigma(double sig_eff, const KVec& sigd_i);
	double DG_DI1(const double alpha);
	double DG_DJ2(const double theta, const double J2, const double alpha);
	double DG_Dtheta(const double theta, const double J2, const double alpha);
	double Dtheta_DJ2(const double theta, const double J2);
	double Dtheta_DJ3(const double theta, const double J3);
	KMat s_odot_s(const KVec& vec);
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
