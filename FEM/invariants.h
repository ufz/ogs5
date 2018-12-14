/**
 * \copyright
 * Copyright (c) 2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef INVARIANTS_H
#define INVARIANTS_H

#include <cfloat>
#include "Eigen/Dense"
// FEM-Makros
#include "makros.h"
#include "tools.h"

namespace SolidMath
{
typedef Eigen::Matrix<double, 6, 1> KVec;
typedef Eigen::Matrix<double, 6, 6> KMat;

KMat Initialise_ident();
KMat Initialise_P_dev();
KMat Initialise_P_sph();
KVec Initialise_ivec();

static const KMat ident(Initialise_ident());  // identity matrix
static const KMat P_dev(Initialise_P_dev());  // deviatoric projection matrix
static const KMat P_sph(Initialise_P_sph());  // spherical projection matrix
static const KVec ivec(
    Initialise_ivec());  // Kelvin mapping of 2nd order identity

// Kelvin/Voigt mapping routines for 6D vectors
KVec Voigt_to_Kelvin_Stress(const std::vector<double>& voigt_stress);
KVec Voigt_to_Kelvin_Strain(const std::vector<double>& voigt_strain);
void Kelvin_to_Voigt_Stress(const KVec& kelvin_stress,
                            std::vector<double>& voigt_stress);
void Kelvin_to_Voigt_Strain(const KVec& kelvin_strain,
                            std::vector<double>& voigt_strain);

// Maps a 6D Kelvin vector back into 3D Tensor coordinates
Eigen::Matrix<double, 3, 3> KelvinVectorToTensor(const KVec& vec);

// Maps a 2nd order 3D Tensor into Kelvin representation
KVec TensorToKelvinVector(const Eigen::Matrix<double, 3, 3>& tens);

// Task: calculates second deviatoric invariant. Note that this routine requires
// a Kelvin mapped DEVIATORIC vector. A deviatoric mapping was not done here in
// order to avoid unnecessary calculations
double CalJ2(const KVec& dev_vec);

// Task: calculates effective stress. Note that this routine requires a
// Kelvin mapped DEVIATORIC stress vector
double CalEffectiveStress(const KVec& dev_stress);

// Task: calculates third deviatoric invariant. Note that this routine requires
// a Kelvin mapped DEVIATORIC vector A deviatoric mapping was not done here in
// order to avoid unnecessary calculations
double CalJ3(const KVec& dev_vec);

// Takes 2nd Order tensor in Kelvin representation and returns its
// inverse in Kelvin representation
KVec InvertVector(const KVec& vec);

// calculates Lode angle. Note that this routine requires a
// Kelvin mapped DEVIATORIC vector
// A deviatoric mapping was not done here in order to avoid unnecessary
// calculations
double CalLodeAngle(const KVec& dev_vec);

// calculates first invariant (trace) of a vector-mapped 2nd order tensor.
double CalI1(const KVec& vec);
}  // namespace SolidMath
#endif  // INVARIANTS_H
