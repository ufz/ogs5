/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \author Wenqing Wang
 *  \file   TheoreticalPorosity.h
 *  Created on June 28, 2018, 12:17 PM
 */

#pragma once

class CRFProcess;

namespace FiniteElement
{
class CElement;
}

namespace MaterialLib
{
/// The porosity is calculated by solving
/// dn/dt=(a-n)(dp/dt/K-3*alpha*dT/dt + d(e_v)/dt)
class TheoreticalPorosity
{
public:
	TheoreticalPorosity(const CRFProcess& process_T, const CRFProcess& process_H,
                     const CRFProcess& process_M, const double n0,
                     const double n_min, const double n_max,
                     const double K, const double alpha_B, double alpha_T);

	double getPorosity(const FiniteElement::CElement& fem_assembler, const double t, const double dt);

private:
	const CRFProcess& _process_T;
	const CRFProcess& _process_H;
	const CRFProcess& _process_M;

	const double _n_min; /// Minimum porosity
	const double _n_max; /// Maximum porosity
	const double _n0; /// Initial porosity.
	const double _Ks; /// Bulk modulus of solid.
	const double _alpha_B; /// Biots constant.
	const double _alpha_T; /// Linear thermal expansion.
	const int _index_T;
	const int _index_p;
	int _index_u[3];
};
} // end of name space.
