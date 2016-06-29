
#include "invariants.h"

#include <cfloat>
// FEM-Makros
#include "makros.h"
#include "tools.h"

namespace SolidMath
{

void InitialiseProjectionTensors()
{
	//set identity matrix here
	ident.setIdentity();

	//set identity vector (Kelvin mapping of 2nd order Identity)
	for (size_t i=0; i<3; i++)
	{
		ivec(i) = 1.;
		ivec(i+3) = 0.;
	}

	const double third(1./3.);
	//deviatoric projection
	P_dev.setIdentity();
	P_sph.setZero(6,6);
	for (size_t i=0; i<3; i++)
		for (size_t j=0; j<3; j++){
			P_dev(i,j) -= third;
			P_sph(i,j) = third;
		}
}

//Maps a 6D Kelvin vector back into 3D Tensor coordinates
Eigen::Matrix<double,3,3>  KelvinVectorToTensor(const Eigen::Matrix<double,6,1> &vec)
{
	Eigen::Matrix<double,3,3> tens;
	tens(0,0) = vec(0);
	tens(1,1) = vec(1);
	tens(2,2) = vec(2);
	tens(0,1) = tens(1,0) = vec(3)/std::sqrt(2.);
	tens(1,2) = tens(2,1) = vec(4)/std::sqrt(2.);
	tens(0,2) = tens(2,0) = vec(5)/std::sqrt(2.);
	return tens;
}

//Maps a 2nd order 3D Tensor into Kelvin representation
Eigen::Matrix<double,6,1>  TensorToKelvinVector(const Eigen::Matrix<double,3,3> &tens)
{
	Eigen::Matrix<double,6,1> vec;
	vec(0) = tens(0,0);
	vec(1) = tens(1,1);
	vec(2) = tens(2,2);
	vec(3) = tens(0,1) * std::sqrt(2.);
	vec(4) = tens(1,2) * std::sqrt(2.);
	vec(5) = tens(0,2) * std::sqrt(2.);
	return vec;
}

//Task: calculates second deviatoric invariant. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector.
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ2(const Eigen::Matrix<double,6,1> &dev_vec)
{
	double s(0.);
	s = dev_vec.transpose() * dev_vec; //Kelvin mapping, deviator
	return 0.5 * s;
}

//Task: calculates effective stress. Note that this routine requires a
//Kelvin mapped DEVIATORIC stress vector
double CalEffectiveStress(const Eigen::Matrix<double,6,1> &dev_stress)
{
	double s(0.);
	s = 3.* CalJ2(dev_stress); //Kelvin mapped deviatoric stress has to be used
	return sqrt(s);
}

//Task: calculates third deviatoric invariant. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ3(const Eigen::Matrix<double,6,1> &dev_vec)
{
	Eigen::Matrix<double,3,3> tens;
	tens = KelvinVectorToTensor(dev_vec); //TN: Not strictly necessary. Can be written explicitly for vector coordinates
	return tens.determinant();
}

//Takes 2nd Order tensor in Kelvin representation and returns its
//inverse in Kelvin representation
Eigen::Matrix<double,6,1> InvertVector(const Eigen::Matrix<double,6,1> &vec)
{
	Eigen::Matrix<double,3,3> tens = KelvinVectorToTensor(vec);
	return TensorToKelvinVector(tens.inverse());
}

//calculates Lode angle. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalLodeAngle(const Eigen::Matrix<double,6,1> &dev_vec)
{
	const double J2(CalJ2(dev_vec));
	double theta,thetaR;

	if (std::abs(J2) > DBL_EPSILON) //catch zero-stress state
		theta = -3. * std::sqrt(3.) * CalJ3(dev_vec) / (2. * std::pow(J2,1.5));
	else
		theta = 0.;

	thetaR = std::min(std::max(theta,-1.+DBL_EPSILON),1.-DBL_EPSILON); //strict limits for trigonometric functions
	return 1./3. * std::asin(thetaR);
}

//calculates first invariant (trace) of a vector-mapped 2nd order tensor.
double CalI1(const Eigen::Matrix<double,6,1> &vec)
{
	return vec(0) + vec(1) + vec(2);
}

}
