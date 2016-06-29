#ifndef INVARIANTS_H
#define INVARIANTS_H

#include <Eigen/Dense>

namespace SolidMath
{

static Eigen::Matrix<double,6,6> P_dev; //deviatoric projection matrix
static Eigen::Matrix<double,6,6> P_sph; //spherical projection matrix
static Eigen::Matrix<double,6,6> ident; //identity matrix
static Eigen::Matrix<double,6,1> ivec; //Kelvin mapping of 2nd order identity

void InitialiseProjectionTensors();

//Maps a 6D Kelvin vector back into 3D Tensor coordinates
Eigen::Matrix<double,3,3>  KelvinVectorToTensor(const Eigen::Matrix<double,6,1> &vec);

//Maps a 2nd order 3D Tensor into Kelvin representation
Eigen::Matrix<double,6,1>  TensorToKelvinVector(const Eigen::Matrix<double,3,3> &tens);

//Task: calculates second deviatoric invariant. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector.
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ2(const Eigen::Matrix<double,6,1> &dev_vec);

//Task: calculates effective stress. Note that this routine requires a
//Kelvin mapped DEVIATORIC stress vector
double CalEffectiveStress(const Eigen::Matrix<double,6,1> &dev_stress);

//Task: calculates third deviatoric invariant. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalJ3(const Eigen::Matrix<double,6,1> &dev_vec);

//Takes 2nd Order tensor in Kelvin representation and returns its
//inverse in Kelvin representation
Eigen::Matrix<double,6,1> InvertVector(const Eigen::Matrix<double,6,1> &vec);

//calculates Lode angle. Note that this routine requires a
//Kelvin mapped DEVIATORIC vector
//A deviatoric mapping was not done here in order to avoid unnecessary calculations
double CalLodeAngle(const Eigen::Matrix<double,6,1> &dev_vec);

//calculates first invariant (trace) of a vector-mapped 2nd order tensor.
double CalI1(const Eigen::Matrix<double,6,1> &vec);

}
#endif // INVARIANTS_H


