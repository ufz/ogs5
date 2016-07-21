#include "burgers.h"
#include "PhysicalConstant.h"
namespace Burgers
{
SolidBurgers::SolidBurgers(const Math_Group::Matrix& data)
{
	GK0 = data(0); // Kelvin shear modulus
	mK = data(1); // dependency parameter for "
	etaK0 = data(2); // Kelvin viscosity
	mvK = data(3); // dependency parameter for "
	GM0 = data(4); // Maxwell shear modulus
	KM0 = data(5); // Maxwell bulk modulus
	etaM0 = data(6); // Maxwell viscosity
	mvM = data(7); // dependency parameter for "
	m_GM = data(8); // slope of elesticity temperature dependence
	m_KM = data(9); // slope of elesticity temperature dependence
	T_ref = data(10); // reference temperature dependency parameter for "
	B = data(11); // constant factor for Arrhenius term
	Q = data(12); // activation energy in Arrhenius term

	GM = GM0;
	KM = KM0;
	GK = GK0;
	etaK = etaK0;
	etaM = etaM0;

	if (!T_Process)
	{
		m_GM = 0.;
		m_KM = 0.;
		B = 1.;
		Q = 0.; // for cutting off Arrhenius term
		T_ref = 273.15;
	}
}

/**************************************************************************
   FEMLib-Method: Burgers::UpdateBurgersProperties()
   Task: Updates BURGERS material parameters in LUBBY2 fashion
   Programing:
   07/2014 TN Implementation
**************************************************************************/
void SolidBurgers::UpdateBurgersProperties(double s_eff, const double Temperature)
{
	const double dT(Temperature - T_ref);
	GM = GM0 + m_GM * dT;
	KM = KM0 + m_KM * dT;

	s_eff *= GM;

	GK = GK0 * std::exp(mK * s_eff);
	etaK = etaK0 * std::exp(mvK * s_eff);
	etaM = etaM0 * std::exp(mvM * s_eff) * B
	       * std::exp(Q * (-dT) / (PhysicalConstant::IdealGasConstant * Temperature * T_ref));
	//	if (etaM / etaM0 < 1.e-2)
	//		std::cout << "WARNING: Maxwell viscosity sank to 100th of original value." << std::endl;
}

/**************************************************************************
   FEMLib-Method: Burgers::CalResidualBurgers()
   Task: Calculates the 12x1 residual vector. Implementation fully implicit only.
   Programing:
   06/2014 TN Implementation
**************************************************************************/
void SolidBurgers::CalResidualBurgers(const double dt, const KVec& strain_curr, const KVec& stress_curr,
                                      KVec& strain_Kel_curr, const KVec& strain_Kel_t, KVec& strain_Max_curr,
                                      const KVec& strain_Max_t, Eigen::Matrix<double, 18, 1>& res)
{
	// calculate stress residual
	res.block<6, 1>(0, 0) = stress_curr - 2. * (strain_curr - strain_Kel_curr - strain_Max_curr);

	// calculate Kelvin strain residual
	res.block<6, 1>(6, 0) = 1. / dt * (strain_Kel_curr - strain_Kel_t)
	                        - 1. / (2. * etaK) * (GM * stress_curr - 2. * GK * strain_Kel_curr);

	// calculate Maxwell strain residual
	res.block<6, 1>(12, 0) = 1. / dt * (strain_Max_curr - strain_Max_t) - GM / (2. * etaM) * stress_curr;
}

/**************************************************************************
   FEMLib-Method: Burgers::CalJacobianBurgers()
   Task: Calculates the 12x12 Jacobian. Implementation fully implicit only.
   Programing:
   06/2014 TN Implementation
**************************************************************************/
void SolidBurgers::CalJacobianBurgers(const double dt, Eigen::Matrix<double, 18, 18>& Jac, const double s_eff,
                                      const KVec& sig_i, const KVec& eps_K_i)
{
	// Check Dimension of Jacobian
	// assert(Jac.cols() == 18 && Jac.rows() && 18);
	Jac.setZero(18, 18);

	// build G_11
	Jac.block<6, 6>(0, 0) = SolidMath::ident;

	// build G_12
	Jac.block<6, 6>(0, 6) = 2. * SolidMath::ident;

	// build G_13
	Jac.block<6, 6>(0, 12) = 2. * SolidMath::ident;

	// build G_21
	Jac.block<6, 6>(6, 0) = -GM / (2. * etaK) * SolidMath::ident;

	// build G_22
	Jac.block<6, 6>(6, 6) = (1. / dt + GK / etaK) * SolidMath::ident;

	// nothing to do for G_23

	// build G_31
	Jac.block<6, 6>(12, 0) = -GM / (2. * etaM) * SolidMath::ident;

	// nothing to do for G_32

	// build G_33
	Jac.block<6, 6>(12, 12) = 1. / dt * SolidMath::ident;

	if (s_eff > 0.)
	{
		const KVec dG_K = mK * 3. * GK * GM / (2. * s_eff) * sig_i;
		const KVec dmu_vK = mvK * 3. * GM * etaK / (2. * s_eff) * sig_i;
		const KVec dmu_vM = mvM * 3. * GM * etaM / (2. * s_eff) * sig_i;
		const KVec eps_K_aid = 1. / (etaK * etaK) * (GM * sig_i - 2. * GK * eps_K_i);

		// build G_21
		Jac.block<6, 6>(6, 0) += 0.5 * eps_K_aid * dmu_vK.transpose() + 1. / etaK * eps_K_i * dG_K.transpose();

		// build G_31
		Jac.block<6, 6>(12, 0) += GM / (2. * etaM * etaM) * sig_i * dmu_vM.transpose();
	}
}

/**************************************************************************
   FEMLib-Method: Burgers::CaldGdEBurgers()
   Task: Calculates the 12x6 derivative of the residuals with respect to total strain. Implementation fully implicit
only.
   Programing:
   06/2014 TN Implementation
**************************************************************************/
void SolidBurgers::CaldGdEBurgers(Eigen::Matrix<double, 18, 6>& dGdE)
{
	// Check Dimension of dGdE
	// assert(dGdE.cols() == 6 && dGdE.rows() == 18);
	dGdE.setZero(18, 6);
	dGdE.block<6, 6>(0, 0) = -2. * SolidMath::P_dev;
}
}
