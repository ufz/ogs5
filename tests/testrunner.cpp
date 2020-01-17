/**
 * \file testrunner.cpp
 * 29/4/2010 LB Initial implementation
 *
 * Implementation of the googletest testrunner
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ** INCLUDES **
#include "gtest.h"

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
