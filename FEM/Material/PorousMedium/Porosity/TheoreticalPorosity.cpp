/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \author Wenqing Wang
 *  \file   TheoreticalPorosity.cpp
 *  Created on June 28, 2018, 12:17 PM
 */

#include "TheoreticalPorosity.h"

#include <cmath>

#include "rf_pcs.h"
#include "fem_ele.h"

namespace MaterialLib
{
TheoreticalPorosity::TheoreticalPorosity(const CRFProcess& process_T, const CRFProcess& process_H,
	                    const CRFProcess& process_M, const double n0, const double K,
	                    const double alpha_B, double alpha_T)
    : _process_T(process_T), _process_H(process_H), _process_M(process_M), _n_min(0.001), _n_max(0.999), _n0(n0),
      _Ks( K / (1.0 - alpha_B)), _alpha_B(alpha_B), _alpha_T(alpha_T),
      _index_T(_process_T.GetNodeValueIndex(_process_T.GetPrimaryVName(0))),
      _index_p(_process_H.GetNodeValueIndex(_process_H.GetPrimaryVName(0)))
{
	for (std::size_t i = 0; i < process_M.GetPrimaryVNumber(); i++)
	{
		_index_u[i] = _process_M.GetNodeValueIndex(_process_M.GetPrimaryVName(i));
	}
}

double TheoreticalPorosity::getPorosity(const FiniteElement::CElement& fem_assembler, const double t, const double dt)
{
	double dot_T
	    = (fem_assembler.interpolate(_index_T + 1, &_process_T) - fem_assembler.interpolate(_index_T, &_process_T))
	      / dt;
	double dot_p
	    = (fem_assembler.interpolate(_index_p + 1, &_process_H) - fem_assembler.interpolate(_index_p, &_process_H))
	      / dt;

	double d_e_v = 0.0;
	for (std::size_t i = 0; i < _process_M.GetPrimaryVNumber(); i++)
	{
		d_e_v += fem_assembler.grad(_index_u[i] + 1, i, &_process_M)

		           - fem_assembler.grad(_index_u[i], i, &_process_M);
	}
	if (fem_assembler.isAxisymmetry())
	{
		d_e_v += (fem_assembler.interpolate(_index_u[0] + 1, &_process_M)
		            - fem_assembler.interpolate(_index_u[0], &_process_M))
		           / fem_assembler.getRadius();
	}

	const double a = dot_p / _Ks - 3.0 * _alpha_T * dot_T + d_e_v / dt;
	const double val = std::max(_n_min, _alpha_B - (_alpha_B - _n0) * std::exp(-a * t));
	return std::min(_n_max, val);
}

} // end of name space.
