/*! \file ShapeFunctionPool.cpp
    \brief Compute shape functions and their gradients with respect to
     the local coodinates, and store the results

     \author Wenqing Wang
     \date Feb. 2015

     \copyright
      Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
             See accompanying file LICENSE.txt or
             http://www.opengeosys.org/project/license
*/

#include "ShapeFunctionPool.h"

#include <cassert> /* assert */
#include "fem_ele.h"

namespace FiniteElement
{
ShapeFunctionPool::ShapeFunctionPool(
    const std::vector<MshElemType::type>& elem_types, CElement& quadrature,
    const int num_sample_gs_pnts)
{
    int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES];
    int dim_elem[MshElemType::NUM_ELEM_TYPES];

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

    // std::vector<int> elem_type_ids
    const std::size_t n_ele_types = elem_types.size();
    _shape_function.resize(n_ele_types);
    _shape_function_center.resize(n_ele_types);
    _grad_shape_function.resize(n_ele_types);
    _grad_shape_function_center.resize(n_ele_types);
    for (std::size_t i = 0; i < elem_types.size(); i++)
    {
        const MshElemType::type e_type = elem_types[i];
        if (e_type == MshElemType::INVALID)
            continue;
        // Set number of integration points.
        quadrature.SetGaussPointNumber(num_sample_gs_pnts);
        quadrature.SetIntegrationPointNumber(e_type);

        const int type_id = static_cast<int>(e_type) - 1;
        int num_int_pnts = quadrature.GetNumGaussPoints();
        const int num_nodes =
            num_elem_nodes[quadrature.getOrder() - 1][type_id];

        std::vector<double> elem_shape_function_center(num_nodes);
        _shape_function_center[type_id] = elem_shape_function_center;

        const int size_shape_fct = num_nodes * num_int_pnts;
        std::vector<double> elem_shape_function(size_shape_fct);
        _shape_function[type_id] = elem_shape_function;

        std::vector<double> elem_grad_shape_function(dim_elem[type_id] *
                                                     size_shape_fct);
        _grad_shape_function[type_id] = elem_grad_shape_function;

        std::vector<double> elem_grad_shape_function_center(dim_elem[type_id] *
                                                            num_nodes);
        _grad_shape_function_center[type_id] = elem_grad_shape_function_center;
    }

    computeQuadratures(elem_types, num_elem_nodes, dim_elem, quadrature,
                       num_sample_gs_pnts);
}

void ShapeFunctionPool::computeQuadratures(
    const std::vector<MshElemType::type>& elem_types,
    const int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES],
    const int dim_elem[], CElement& quadrature, const int num_sample_gs_pnts)
{
    const int order = quadrature.getOrder();
    for (std::size_t i = 0; i < elem_types.size(); i++)
    {
        const MshElemType::type e_type = elem_types[i];
        if (e_type == MshElemType::INVALID)
            continue;

        const int type_id = static_cast<int>(e_type) - 1;
        quadrature.ConfigShapefunction(e_type);

        const int nnodes = num_elem_nodes[order - 1][type_id];
        const int elem_dim = dim_elem[type_id];

        double* shape_function_center_values =
            _shape_function_center[type_id].data();
        quadrature.SetCenterGP(e_type);
        quadrature.ComputeShapefct(order, shape_function_center_values);
        double* grad_shape_function_center_values =
            _grad_shape_function_center[type_id].data();
        quadrature.computeGradShapefctLocal(order,
                                            grad_shape_function_center_values);

        double* shape_function_values = _shape_function[type_id].data();
        double* dshape_function_values = _grad_shape_function[type_id].data();
        // Set number of integration points.
        quadrature.SetGaussPointNumber(num_sample_gs_pnts);
        quadrature.SetIntegrationPointNumber(e_type);

        for (int gp = 0; gp < quadrature.GetNumGaussPoints(); gp++)
        {
            int gp_r, gp_s, gp_t;
            quadrature.SetGaussPoint(e_type, gp, gp_r, gp_s, gp_t);
            double* shape_function_values_gs =
                &shape_function_values[gp * nnodes];
            quadrature.ComputeShapefct(order, shape_function_values_gs);

            double* dshape_function_values_gs =
                &dshape_function_values[gp * nnodes * elem_dim];
            quadrature.computeGradShapefctLocal(order,
                                                dshape_function_values_gs);
        }
    }
}

const double* ShapeFunctionPool::getShapeFunctionValues(
    const MshElemType::type elem_type) const
{
    return _shape_function[static_cast<int>(elem_type) - 1].data();
}

unsigned ShapeFunctionPool::getShapeFunctionArraySize(
    const MshElemType::type elem_type) const
{
    return _shape_function[static_cast<int>(elem_type) - 1].size();
}

const double* ShapeFunctionPool::getShapeFunctionCenterValues(
    const MshElemType::type elem_type) const
{
    return _shape_function_center[static_cast<int>(elem_type) - 1].data();
}

const double* ShapeFunctionPool::getGradShapeFunctionValues(
    const MshElemType::type elem_type) const
{
    return _grad_shape_function[static_cast<int>(elem_type) - 1].data();
}

const double* ShapeFunctionPool::getGradShapeFunctionCenterValues(
    const MshElemType::type elem_type) const
{
    return _grad_shape_function_center[static_cast<int>(elem_type) - 1].data();
}

}  // namespace FiniteElement
