#ifndef BURGERS_H
#define BURGERS_H

#include <cfloat>
// FEM-Makros
#include "makros.h"
#include "tools.h"
#include "invariants.h"
#include "Eigen/Dense"

namespace Math_Group
{class Matrix;
}

namespace SolidMath
{class Invariants;
}

namespace Burgers{
using Math_Group::Matrix;
using SolidMath::Invariants;
typedef Eigen::Matrix<double,6,1> KVec;

class SolidBurgers
{
public:
    SolidBurgers(const Matrix* data);
	~SolidBurgers()
	{
		delete smath;
		smath = NULL;
	}
    //basic material parameters
    double GK0, GM0, KM0, etaK0, etaM0, mK, mvK, mvM, l0, T_ref, m_GM, m_KM, B, Q;
    //solution dependent values
    double GM, KM, GK, etaK, etaM;

    void UpdateBurgersProperties(const double s_eff, const double Delta_T);
	void CalResidualBurgers(const double dt, const KVec &strain_curr,
							const KVec &stress_curr, KVec &strain_Kel_curr,
							const KVec &strain_Kel_t, KVec &strain_Max_curr,
							const KVec &strain_Max_t, Eigen::Matrix<double,18,1> &res);
	void CalJacobianBurgers(const double dt, Eigen::Matrix<double,18,18> &Jac, const double s_eff,
							const KVec &sig_i, const KVec &eps_K_i);
    void CaldGdEBurgers(Eigen::Matrix<double,18,6> &dGdE);

    Invariants* smath;

};
}
#endif // BURGERS_H
