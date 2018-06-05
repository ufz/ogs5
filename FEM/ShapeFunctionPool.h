/*! \file ShapeFunctionPool.h
    \brief Compute shape functions and their gradients with respect to
     the local coodinates, and store the results

     \author Wenqing Wang
     \date Feb. 2015

     \copyright
      Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
             See accompanying file LICENSE.txt or
             http://www.opengeosys.org/project/license
*/
#ifndef OGS_SHAPEFUNCTIONPOOL_H
#define OGS_SHAPEFUNCTIONPOOL_H

#include <vector>

#include "MSHEnums.h"

namespace FiniteElement
{
class CElement;

class ShapeFunctionPool
{
public:
	/*!
	    \param elem_types          All involved element types.
	    \param quadrature          Numerical integration object.
	    \param num_sample_gs_pnts  Number of sample Gauss points.
	*/
	ShapeFunctionPool(const std::vector<MshElemType::type>& elem_types, CElement& quadrature,
	                  const int num_sample_gs_pnts);
	~ShapeFunctionPool() {}

	/// Get shape function values of an element type
	const double* getShapeFunctionValues(const MshElemType::type elem_type) const;
	/// Get the size of shape function array of an element type.
	unsigned getShapeFunctionArraySize(const MshElemType::type elem_type) const;

	/// Get shape function values at the element centroid of an element type
	const double* getShapeFunctionCenterValues(const MshElemType::type elem_type) const;

	/// Get the values of the gradient of shape function of an element type
	const double* getGradShapeFunctionValues(const MshElemType::type elem_type) const;

	/// Get the gradient of shape function values at the element centroid of an element type
	const double* getGradShapeFunctionCenterValues(const MshElemType::type elem_type) const;

private:
	/// Results of shape functions of all integration points.
	std::vector< std::vector<double> > _shape_function;

	/// Results of shape functions of all integration points at element centroid.
	std::vector< std::vector<double> > _shape_function_center;

	/// Results of the gradient of shape functions with respect to
	/// local coordinates of all integration points.
	std::vector< std::vector<double> > _grad_shape_function;

	/// Results of the gradient of shape functions of all integration points at element centroid.
	std::vector< std::vector<double> > _grad_shape_function_center;

	void computeQuadratures(const std::vector<MshElemType::type>& elem_types,
	                        const int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES],
	                        const int dim_elem[],
	                        CElement& quadrature,
	                        const int num_sample_gs_pnts);
};
} // end namespace

#endif
