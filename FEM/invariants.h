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
    Invariants();
    ~Invariants();
    double CalEffectiveStress(const Eigen::Matrix<double,6,1> &dev_stress);
    double CalJ2(const Eigen::Matrix<double,6,1> &dev_vec);
    double CalJ3(const Eigen::Matrix<double,6,1> &dev_vec);
    Eigen::Matrix<double,3,3>  KelvinVectorToTensor(const Eigen::Matrix<double,6,1> &vec);
    Eigen::Matrix<double,6,1>  TensorToKelvinVector(const Eigen::Matrix<double,3,3> &tens);
    Eigen::Matrix<double,6,1> InvertVector(const Eigen::Matrix<double,6,1> &vec);
    double CalLodeAngle(const Eigen::Matrix<double,6,1> &dev_vec);
    double CalI1(const Eigen::Matrix<double,6,1> &vec);
    Eigen::Matrix<double,6,6> P_dev; //deviatoric projection matrix
    Eigen::Matrix<double,6,6> P_sph; //spherical projection matrix
    Eigen::Matrix<double,6,6> ident; //identity matrix
    Eigen::Matrix<double,6,1> ivec; //Kelvin mapping of 2nd order identity
};
}
#endif // INVARIANTS_H

