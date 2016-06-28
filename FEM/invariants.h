#ifndef INVARIANTS_H
#define INVARIANTS_H

#include <cfloat>
// FEM-Makros
#include "makros.h"
#include "tools.h"

#include "Eigen/Dense"

namespace Math_Group
{class Matrix;
}

namespace SolidMath{
using Math_Group::Matrix;

class Invariants
{
public:
	Invariants()
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
    ~Invariants();
	const double CalEffectiveStress(const Eigen::Matrix<double,6,1> &dev_stress);
	const double CalJ2(const Eigen::Matrix<double,6,1> &dev_vec);
	const double CalJ3(const Eigen::Matrix<double,6,1> &dev_vec);
	const Eigen::Matrix<double,3,3>  KelvinVectorToTensor(const Eigen::Matrix<double,6,1> &vec);
	const Eigen::Matrix<double,6,1>  TensorToKelvinVector(const Eigen::Matrix<double,3,3> &tens);
	const Eigen::Matrix<double,6,1> InvertVector(const Eigen::Matrix<double,6,1> &vec);
	const double CalLodeAngle(const Eigen::Matrix<double,6,1> &dev_vec);
	const double CalI1(const Eigen::Matrix<double,6,1> &vec);
	Eigen::Matrix<double,6,6> P_dev; //deviatoric projection matrix
	Eigen::Matrix<double,6,6> P_sph; //spherical projection matrix
	Eigen::Matrix<double,6,6> ident; //identity matrix
	Eigen::Matrix<double,6,1> ivec; //Kelvin mapping of 2nd order identity
};
}
#endif // INVARIANTS_H

