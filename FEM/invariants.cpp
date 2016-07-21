#include "invariants.h"

namespace SolidMath
{
KMat Initialise_ident()
{
	KMat I;
	I.setIdentity();
	return I;
}

KMat Initialise_P_dev()
{
	KMat P;
	P.setIdentity();
	const double third(1. / 3.);
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
		{
			P(i, j) -= third;
		}
	return P;
}

KMat Initialise_P_sph()
{
	KMat P;
	P.setZero(6, 6);
	const double third(1. / 3.);
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
		{
			P(i, j) = third;
		}
	return P;
}

KVec Initialise_ivec()
{
	KVec I;
	for (size_t j = 0; j < 3; j++)
	{
		I(j) = 1.;
		I(j + 3) = 0.;
	}
	return I;
}

/**************************************************************************
   Voigt_to_Kelvin_Strain
   Task: Maps a strain vector in Voigt notation into on in Kelvin notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
KVec Voigt_to_Kelvin_Strain(const std::vector<double>& voigt_strain)
{
	KVec kelvin_strain;
	for (size_t i = 0; i < 3; i++)
	{
		// Normal components
		kelvin_strain(i) = voigt_strain[i];
		// Shear components
		kelvin_strain(i + 3) = voigt_strain[i + 3] / sqrt(2.);
	}
	return kelvin_strain;
}

/**************************************************************************
   Voigt_to_Kelvin_Stress
   Task: Maps a stress vector in Voigt notation into on in Kelvin notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
KVec Voigt_to_Kelvin_Stress(const std::vector<double>& voigt_stress)
{
	KVec kelvin_stress;
	for (size_t i = 0; i < 3; i++)
	{
		// Normal components
		kelvin_stress(i) = voigt_stress[i];
		// Shear components
		kelvin_stress(i + 3) = voigt_stress[i + 3] * sqrt(2.);
	}
	return kelvin_stress;
}

/**************************************************************************
   Kelvin_to_Voigt_Strain()
   Task: Maps a strain vector in Kelvin notation into on in Voigt notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void Kelvin_to_Voigt_Strain(const KVec& kelvin_strain, std::vector<double>& voigt_strain)
{
	for (size_t i = 0; i < 3; i++)
	{
		// Normal components
		voigt_strain[i] = kelvin_strain(i);
		// Shear components
		voigt_strain[i + 3] = kelvin_strain(i + 3) * sqrt(2.);
	}
	return;
}

/**************************************************************************
   Kelvin_to_Voigt_Stress()
   Task: Maps a stress vector in Kelvin notation into on in Voigt notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void Kelvin_to_Voigt_Stress(const KVec& kelvin_stress, std::vector<double>& voigt_stress)
{
	for (size_t i = 0; i < 3; i++)
	{
		// Normal components
		voigt_stress[i] = kelvin_stress(i);
		// Shear components
		voigt_stress[i + 3] = kelvin_stress(i + 3) / sqrt(2.);
	}
	return;
}

// Maps a 6D Kelvin vector back into 3D Tensor coordinates
Eigen::Matrix<double, 3, 3> KelvinVectorToTensor(const KVec& vec)
{
	Eigen::Matrix<double, 3, 3> tens;
	tens(0, 0) = vec(0);
	tens(1, 1) = vec(1);
	tens(2, 2) = vec(2);
	tens(0, 1) = tens(1, 0) = vec(3) / std::sqrt(2.);
	tens(1, 2) = tens(2, 1) = vec(4) / std::sqrt(2.);
	tens(0, 2) = tens(2, 0) = vec(5) / std::sqrt(2.);
	return tens;
}

// Maps a 2nd order 3D Tensor into Kelvin representation
KVec TensorToKelvinVector(const Eigen::Matrix<double, 3, 3>& tens)
{
	KVec vec;
	vec(0) = tens(0, 0);
	vec(1) = tens(1, 1);
	vec(2) = tens(2, 2);
	vec(3) = tens(0, 1) * std::sqrt(2.);
	vec(4) = tens(1, 2) * std::sqrt(2.);
	vec(5) = tens(0, 2) * std::sqrt(2.);
	return vec;
}

// Task: calculates second deviatoric invariant. Note that this routine requires a
// Kelvin mapped DEVIATORIC vector.
// A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ2(const KVec& dev_vec)
{
	double s(0.);
	s = dev_vec.transpose() * dev_vec; // Kelvin mapping, deviator
	return 0.5 * s;
}

// Task: calculates effective stress. Note that this routine requires a
// Kelvin mapped DEVIATORIC stress vector
double CalEffectiveStress(const KVec& dev_stress)
{
	double s(0.);
	s = 3. * CalJ2(dev_stress); // Kelvin mapped deviatoric stress has to be used
	return sqrt(s);
}

// Task: calculates third deviatoric invariant. Note that this routine requires a
// Kelvin mapped DEVIATORIC vector
// A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ3(const KVec& dev_vec)
{
	Eigen::Matrix<double, 3, 3> tens;
	tens = KelvinVectorToTensor(dev_vec); // TN: Not strictly necessary. Can be written explicitly for vector
	                                      // coordinates
	return tens.determinant();
}

// Takes 2nd Order tensor in Kelvin representation and returns its
// inverse in Kelvin representation
KVec InvertVector(const KVec& vec)
{
	Eigen::Matrix<double, 3, 3> tens = KelvinVectorToTensor(vec);
	return TensorToKelvinVector(tens.inverse());
}

// calculates Lode angle. Note that this routine requires a
// Kelvin mapped DEVIATORIC vector
// A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalLodeAngle(const KVec& dev_vec)
{
	const double J2(CalJ2(dev_vec));
	double theta, thetaR;

	if (std::abs(J2) > DBL_EPSILON) // catch zero-stress state
		theta = -3. * std::sqrt(3.) * CalJ3(dev_vec) / (2. * std::pow(J2, 1.5));
	else
		theta = 0.;

	thetaR = std::min(std::max(theta, -1. + DBL_EPSILON), 1. - DBL_EPSILON); // strict limits for trigonometric
	                                                                         // functions
	return 1. / 3. * std::asin(thetaR);
}

// calculates first invariant (trace) of a vector-mapped 2nd order tensor.
double CalI1(const KVec& vec)
{
	return vec(0) + vec(1) + vec(2);
}
}
