/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   testGaussAlgorithm.cpp
 *
 * Created on January 21, 2019, 9:40 AM
 *
 */
#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "matrix_class.h"
#include "LinAlg/GaussAlgorithm.h"

using Math_Group::Matrix;
using Math_Group::Vec;

TEST(LinAlg, GaussAlgorithm)
{
    Matrix m(3, 3);
    m(0, 0) = 1.0;
    m(0, 1) = 1.0;
    m(0, 2) = 1.0;
    m(1, 0) = 0.0;
    m(1, 1) = 2.0;
    m(1, 2) = 5.0;
    m(2, 0) = 2.0;
    m(2, 1) = 5.0;
    m(2, 2) = -1.0;

    // A copy of m to test its inverse.
    Matrix m0(m);

    MathLib::GaussAlgorithm<Math_Group::Matrix> linear_solver(m, m.Rows());

    Vec b(3);
    b(0) = 6;
    b(1) = -4;
    b(2) = 27;

    linear_solver.execute(b.getEntryArray());

    ASSERT_NEAR(5.0, b(0), 1.e-10);
    ASSERT_NEAR(3.0, b(1), 1.e-10);
    ASSERT_NEAR(-2.0, b(2), 1.e-10);

    double x[3];
    x[0] = 6;
    x[1] = -4;
    x[2] = 20;

    linear_solver.executeWithExistedElimination(x);
    ASSERT_NEAR(6.0, x[0], 1.e-10);
    ASSERT_NEAR(4.0 / 3.0, x[1], 1.e-10);
    ASSERT_NEAR(-4.0 / 3.0, x[2], 1.e-10);

    // Test inverse of matrix
    double r[3];
    Matrix inverse_m(3, 3);  // m^{-1}
    for (std::size_t j = 0; j < m.Rows(); j++)
    {
        for (std::size_t i = 0; i < m.Rows(); i++)
        {
            r[i] = 0.0;
            if (i == j)
                r[i] = 1.0;
        }
        // the i_th column of the invJac matrix
        linear_solver.executeWithExistedElimination(r);
        for (std::size_t i = 0; i < m.Rows(); i++)
            inverse_m(i, j) = r[i];
    }

    // m0 * m^{-1} --> m
    m = 0.0;
    m0.multi(inverse_m, m);

    for (std::size_t i = 0; i < m.Rows(); i++)
    {
        for (std::size_t j = 0; j < m.Rows(); j++)
        {
            if (i == j)
                ASSERT_NEAR(1.0, m(i, j), 1.e-10);
            else
                ASSERT_NEAR(0.0, m(i, j), 1.e-10);
        }
    }
}
