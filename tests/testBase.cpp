/**
 * \file Test_Base.cpp
 * 29/4/2010 LB Initial implementation
 *
 * Tests for the Base directory
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ** INCLUDES **
#include "gtest.h"

#include "swap.h"

TEST(Base, SwapInt)
{
    int arg0 = 5;
    int arg1 = 10;
    BASELIB::swap(arg0, arg1);
    ASSERT_EQ(arg0, 10);
    ASSERT_EQ(arg1, 5);
}

TEST(Base, SwapDouble)
{
    double arg0 = 5.0;
    double arg1 = 10.0;
    BASELIB::swap(arg0, arg1);
    ASSERT_EQ(arg0, 10.0);
    ASSERT_EQ(arg1, 5.0);
}

TEST(Base, SwapString)
{
    std::string arg0 = "5";
    std::string arg1 = "10";
    BASELIB::swap(arg0, arg1);
    ASSERT_EQ(arg0, std::string("10"));
    ASSERT_EQ(arg1, std::string("5"));
}
