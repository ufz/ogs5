/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*!
  \file testMatrix.cpp
  \date 2014-04
  \author Wenqing

  Test matrix classes
 */

#include <../gtest/gtest.h>

#include <string>
#include <vector>

#include "matrix_class.h"

using Math_Group::Matrix;
using Math_Group::SymMatrix;
using Math_Group::Vec;

TEST(MATRIX, test_class_Matrix)
{
    Matrix m(3, 3);

    m(0, 0) = 1.;
    m(0, 1) = 2.;
    m(0, 2) = 3.;
    //
    m(1, 0) = 4.;
    m(1, 1) = 5.;
    m(1, 2) = 6.;
    //
    m(2, 0) = 7.;
    m(2, 1) = 8.;
    m(2, 2) = 9.;

    double x[3];
    double b[3];

    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;

    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;

    m.multi(x, b);
    ASSERT_EQ(14., b[0]);
    ASSERT_EQ(32., b[1]);
    ASSERT_EQ(50, b[2]);

    // Copy constructor
    Matrix m1(m);
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
    m1.multi(x, b);
    ASSERT_EQ(14., b[0]);
    ASSERT_EQ(32., b[1]);
    ASSERT_EQ(50, b[2]);

    Matrix m2(3, 3);
    m2 = 0.;
    m1 = 1.0;
    m.multi(m1, m2);
    ASSERT_EQ(6., m2(0, 1));
    ASSERT_EQ(15., m2(1, 2));
    ASSERT_EQ(24., m2(2, 0));

    Matrix m3(3, 3);
    m2 = 1.0;
    m3 = 0.;
    m.multi(m1, m2, m3);
    ASSERT_EQ(3 * 6., m3(0, 1));
    ASSERT_EQ(3 * 15., m3(1, 2));
    ASSERT_EQ(3 * 24., m3(2, 0));

    // Transpose
    m.GetTranspose(m1);
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
    m1.multi(x, b);
    ASSERT_EQ(12., b[0]);
    ASSERT_EQ(15., b[1]);
    ASSERT_EQ(18., b[2]);

    //
    m = 1.0;
    m *= 10.;
    m /= 2.0;
    m += 5.0;
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
    m.multi(x, b);
    ASSERT_EQ(30., b[0]);
    ASSERT_EQ(30., b[1]);
    ASSERT_EQ(30., b[2]);

    m1 = m;
    m1 += m;
    m1 -= m;
    b[0] = 0.0;
    b[1] = 0.0;
    b[2] = 0.0;
    m1.multi(x, b);
    ASSERT_EQ(30., b[0]);
    ASSERT_EQ(30., b[1]);
    ASSERT_EQ(30., b[2]);
}

TEST(MATRIX, test_class_SymMatrix)
{
    SymMatrix m(3);

    m(0, 0) = 1.;
    //
    m(1, 0) = 4.;
    m(1, 1) = 5.;
    //
    m(2, 0) = 7.;
    m(2, 1) = 8.;
    m(2, 2) = 9.;

    ASSERT_EQ(4., m(0, 1));
    ASSERT_EQ(7., m(0, 2));
    ASSERT_EQ(8., m(1, 2));

    Vec x(3);
    Vec b(3);

    x = 1.0;
    b = 0.0;
    m.multi(x.getEntryArray(), b.getEntryArray());
    ASSERT_EQ(12., b[0]);
    ASSERT_EQ(17., b[1]);
    ASSERT_EQ(24., b[2]);

    // Copy constructor
    SymMatrix m1(m);
    b = 0.0;
    m1.multi(x.getEntryArray(), b.getEntryArray());
    ASSERT_EQ(12., b[0]);
    ASSERT_EQ(17., b[1]);
    ASSERT_EQ(24., b[2]);

    Matrix m2(3, 3);
    m2 = 0.;
    m1 = 1.0;
    m.multi(m1, m2);
    ASSERT_EQ(12., m2(0, 1));
    ASSERT_EQ(17., m2(1, 2));
    ASSERT_EQ(24., m2(2, 0));

    Matrix m3(3, 3);
    m2 = 1.0;
    m3 = 0.;
    m.multi(m1, m2, m3);
    ASSERT_EQ(3 * 12., m3(0, 1));
    ASSERT_EQ(3 * 17., m3(1, 2));
    ASSERT_EQ(3 * 24., m3(2, 0));
}

TEST(MATRIX, test_class_Vec)
{
    Vec a(3);

    a = 1.;
    ASSERT_EQ(1., a[2]);

    a[2] = 5;
    ASSERT_EQ(5., a[2]);

    a *= 2.;
    ASSERT_EQ(10., a[2]);

    a /= 5.;
    ASSERT_EQ(2., a[2]);

    a += 5.;
    ASSERT_EQ(7., a[2]);

    Vec b(3);

    b = a;
    ASSERT_EQ(7., b[2]);

    b += a;
    ASSERT_EQ(14., b[2]);

    b -= a;
    ASSERT_EQ(7., b[2]);
}
