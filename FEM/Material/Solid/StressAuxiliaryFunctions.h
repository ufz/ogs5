/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file StressAuxiliaryFunctions.h
 *
 * Created on February 6, 2019, 12:08 PM
 *
 */

#pragma once

namespace SolidProp
{
void getDeviatoricStess(double const* const stress,
                        const int nstress_components, double* s);

double getStressNorm(double const* const s, const int nstress_components);

double* getPrincipleStresses(double const* const s,
                             const int nstress_components);

const double pi = 3.14159265359;

}  // end of namespace
