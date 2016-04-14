/*! \file ShapeFunctionPool.cpp
    \brief Compute shape functions and their gradients with respect to
     the local coodinates, and store the results

     \author Wenqing Wang
     \date Feb. 2015
 
     \copyright
      Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
             See accompanying file LICENSE.txt or
             http://www.opengeosys.org/project/license
*/

#include "ShapeFunctionPool.h"

#include <cassert>     /* assert */
#include "fem_ele.h"

namespace FiniteElement
{
ShapeFunctionPool::ShapeFunctionPool(
			const std::vector<MshElemType::type>& elem_types,
			CElement& quadrature, const int num_sample_gs_pnts)
{
	const std::size_t n_ele_types = elem_types.size();
	_shape_function.reserve(n_ele_types);
	_shape_function_size.reserve(n_ele_types);
	_shape_function_center.reserve(n_ele_types);
	_grad_shape_function.reserve(n_ele_types);
	_grad_shape_function_center.reserve(n_ele_types);

	for (std::size_t i=0; i<n_ele_types; i++)
	{
		_shape_function.push_back(NULL);
		_shape_function_size.push_back(0);
		_shape_function_center.push_back(NULL);
		_grad_shape_function.push_back(NULL);
		_grad_shape_function_center.push_back(NULL);
	}

	int num_elem_nodes[2][MshElemType::LAST];
	int dim_elem[MshElemType::LAST];

	int id = static_cast<int>(MshElemType::LINE) - 1;
	num_elem_nodes[0][id] = 2;
	num_elem_nodes[1][id] = 3;
	dim_elem[id] = 1;

	id = static_cast<int>(MshElemType::QUAD) - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 9;
	dim_elem[id] = 2;

	id = static_cast<int>(MshElemType::QUAD8) - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 8;
	dim_elem[id] = 2;

	id = static_cast<int>(MshElemType::TRIANGLE) - 1;
	num_elem_nodes[0][id] = 3;
	num_elem_nodes[1][id] = 6;
	dim_elem[id] = 2;

	id = static_cast<int>(MshElemType::HEXAHEDRON) - 1;
	num_elem_nodes[0][id] = 8;
	num_elem_nodes[1][id] = 20;
	dim_elem[id] = 3;

	id = static_cast<int>(MshElemType::TETRAHEDRON) - 1;
	num_elem_nodes[0][id] = 4;
	num_elem_nodes[1][id] = 10;
	dim_elem[id] = 3;

	id = static_cast<int>(MshElemType::PRISM) - 1;
	num_elem_nodes[0][id] = 6;
	num_elem_nodes[1][id] = 15;
	dim_elem[id] = 3;

	id = static_cast<int>(MshElemType::PYRAMID) - 1;
	num_elem_nodes[0][id] = 5;
	num_elem_nodes[1][id] = 13;
	dim_elem[id] = 3;

	//std::vector<int> elem_type_ids
	for (std::size_t i=0; i<elem_types.size(); i++)
	{
		const MshElemType::type e_type = elem_types[i];
		if (e_type == MshElemType::INVALID)
			continue;
		// Set number of integration points.
		quadrature.SetGaussPointNumber(num_sample_gs_pnts);
		quadrature.SetIntegrationPointNumber(e_type);

		const int type_id = static_cast<int>(e_type) - 1;
		int num_int_pnts = quadrature.GetNumGaussPoints();
		const int num_nodes = num_elem_nodes[quadrature.getOrder() - 1][type_id];
		_shape_function_center[type_id] = new double[num_nodes];
		const int size_shape_fct = num_nodes * num_int_pnts;
		_shape_function[type_id] = new double[size_shape_fct];
		_shape_function_size[type_id] = size_shape_fct;
		_grad_shape_function[type_id] = new double[dim_elem[type_id] * size_shape_fct];
		_grad_shape_function_center[type_id] = new double[dim_elem[type_id] * num_nodes];
	}

	computeQuadratures(elem_types, num_elem_nodes, dim_elem, quadrature, num_sample_gs_pnts);
}

ShapeFunctionPool::~ShapeFunctionPool()
{
	for (std::size_t i=0; i<_shape_function.size(); i++)
	{
		if (_shape_function[i])
			delete [] _shape_function[i];
		_shape_function[i] = NULL;
	}

	for (std::size_t i=0; i<_shape_function_center.size(); i++)
	{
		if (_shape_function_center[i])
			delete [] _shape_function_center[i];
		_shape_function_center[i] = NULL;
	}

	for (std::size_t i=0; i<_shape_function.size(); i++)
	{
		if (_grad_shape_function[i])
			delete [] _grad_shape_function[i];
		_grad_shape_function[i] = NULL;
	}

	for (std::size_t i=0; i<_grad_shape_function_center.size(); i++)
	{
		if (_grad_shape_function_center[i])
			delete [] _grad_shape_function_center[i];
		_grad_shape_function_center[i] = NULL;
	}
}

void ShapeFunctionPool::
	computeQuadratures(const std::vector<MshElemType::type>& elem_types,
	                   const int num_elem_nodes[2][MshElemType::LAST],
					   const int dim_elem[],
		               CElement& quadrature, const int num_sample_gs_pnts)
{
	const int order = quadrature.getOrder();
	for (std::size_t i=0; i<elem_types.size(); i++)
	{
		const MshElemType::type e_type = elem_types[i];
		if (e_type == MshElemType::INVALID)
			continue;

		const int type_id = static_cast<int>(e_type) - 1;
		quadrature.ConfigShapefunction(e_type);

		const int nnodes = num_elem_nodes[order-1][type_id];
		const int elem_dim = dim_elem[type_id];

		double* shape_function_center_values = _shape_function_center[type_id];
		quadrature.SetCenterGP(e_type);
		quadrature.ComputeShapefct(order, shape_function_center_values);
		double* grad_shape_function_center_values = _grad_shape_function_center[type_id];
		quadrature.computeGradShapefctLocal(order, grad_shape_function_center_values);

		double* shape_function_values = _shape_function[type_id];
		double* dshape_function_values = _grad_shape_function[type_id];
		// Set number of integration points.
		quadrature.SetGaussPointNumber(num_sample_gs_pnts);
		quadrature.SetIntegrationPointNumber(e_type);

		for (int gp = 0; gp < quadrature.GetNumGaussPoints(); gp++)
	    {
			int gp_r, gp_s, gp_t;
			quadrature.SetGaussPoint(e_type, gp, gp_r, gp_s, gp_t);
			double* shape_function_values_gs
				    = &shape_function_values[gp * nnodes];
			quadrature.ComputeShapefct(order, shape_function_values_gs);

			double* dshape_function_values_gs
				    = &dshape_function_values[gp * nnodes * elem_dim];
			quadrature.computeGradShapefctLocal(order, dshape_function_values_gs);
		}
	}
}

double* ShapeFunctionPool::
getShapeFunctionValues(const MshElemType::type elem_type) const
{
	assert(_shape_function[static_cast<int>(elem_type)-1]);
	return _shape_function[static_cast<int>(elem_type)-1];
}

unsigned ShapeFunctionPool::
getShapeFunctionArraySize(const MshElemType::type elem_type) const
{
	return _shape_function_size[static_cast<int>(elem_type)-1];
}

double* ShapeFunctionPool::
getShapeFunctionCenterValues(const MshElemType::type elem_type) const
{
	assert(_shape_function_center[static_cast<int>(elem_type)-1]);
	return _shape_function_center[static_cast<int>(elem_type)-1];
}

double* ShapeFunctionPool::
getGradShapeFunctionValues(const MshElemType::type elem_type) const
{
	assert(_grad_shape_function[static_cast<int>(elem_type)-1]);
	return _grad_shape_function[static_cast<int>(elem_type)-1];
}

double* ShapeFunctionPool::
getGradShapeFunctionCenterValues(const MshElemType::type elem_type) const
{
	assert(_grad_shape_function_center[static_cast<int>(elem_type)-1]);
	return _grad_shape_function_center[static_cast<int>(elem_type)-1];
}

} // end namespace


