/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <math.h>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <string>
#include "NR.h"
using namespace std;


NR::NR(void){}
NR::~NR(void){}

double dfridr(double func(double), const double x/*, const double h, double &err*/)
{
	// Removed
	std::cout << "Error: NR::dfridr not implemented!" << std::endl;
	return numeric_limits<double>::min();
}

double NR::dfridrX(double func(double,double), const double x, const double y/*, const double h, double &err*/)
{
	// Removed
	std::cout << "Error: NR::dfridrX not implemented!" << std::endl;
	return numeric_limits<double>::min();
}

double NR::dfridrY(double func(double,double), const double x, const double y/*, const double h, double &err*/)
{
	// Removed
	std::cout << "Error: NR::dfridrY not implemented!" << std::endl;
	return numeric_limits<double>::min();
}
