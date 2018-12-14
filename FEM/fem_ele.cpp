/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   The members of class Element definitions.
   Designed and programmed by WW, 06/2004
 */

#include "fem_ele.h"

#include <cfloat>
#include <cassert>

#include "msh_elem.h"
#include "rf_pcs.h"
#include "femlib.h"
#include "mathlib.h"

#ifndef USE_PETSC
#include "par_ddc.h"
#endif

#include "ShapeFunctionPool.h"

namespace FiniteElement
{
/**************************************************************************
   FEMLib-Method:
   Task: Constructor of class CElement
   Programing:
   01/2005 WW Implementation
   01/2005 OK 1D case
   01/2006 WW Axisymmetry
   Last modified:
**************************************************************************/
CElement::CElement(int CoordFlag, const int order)
    : MeshElement(NULL),
      Order(order),
      ele_dim(1),
      _ele_global_dim(1),
      nGaussPoints(1),
      nGauss(2),
      ShapeFunction(NULL),
      ShapeFunctionHQ(NULL),
      GradShapeFunction(NULL),
      GradShapeFunctionHQ(NULL),
      _is_mixed_order(false),
      T_Flag(false),
      C_Flag(false),
      F_Flag(false),
      D_Flag(0),
      RD_Flag(false),
      extrapo_method(ExtrapolationMethod::EXTRAPO_LINEAR)
{
    for (int i = 0; i < 2; i++)
    {
        _shape_function_pool_ptr[i] = NULL;
        _shape_function_result_ptr[i] = NULL;
        _grad_shape_function_result_ptr[i] = NULL;
    }

    //
    if (CoordFlag < 0)  // Axisymmetry
    {
        CoordFlag *= -1;
        axisymmetry = true;
    }
    else
        axisymmetry = false;
    //
    dim = CoordFlag / 10;
    coordinate_system = CoordFlag;
    for (int i = 0; i < 4; i++)
        unit[i] = 0.0;

    int dim_jacobian = 1;
    int size_dshapefct = 6;
    int size_dshapefctHQ = 9;
    int max_intgration_points = 27;
    switch (dim)
    {
        case 1:  // OK
            // Memory allocated for maxium 3 nodes elements
            dim_jacobian = 1;
            size_dshapefct = 3 * max_intgration_points;
            size_dshapefctHQ = 3 * max_intgration_points;
            break;
        case 2:
            // Memory allocated for maxium 9 nodes elements
            dim_jacobian = 4;
            size_dshapefct = 18 * max_intgration_points;
            size_dshapefctHQ = 18 * max_intgration_points;
            break;
        case 3:
            // Memory allocated for maxium 20 nodes elements
            dim_jacobian = 9;
            size_dshapefct = 24 * max_intgration_points;
            size_dshapefctHQ = 60 * max_intgration_points;
            //
            break;
    }

    _Jacobian = new double[dim_jacobian];
    _determinants_all = new double[max_intgration_points];
    _inv_jacobian_all = new double[max_intgration_points * dim_jacobian];
    _dshapefct_all = new double[size_dshapefct];
    _dshapefctHQ_all = new double[size_dshapefctHQ];

    time_unit_factor = 1.0;

    if (M_Process)
        D_Flag = 4;
    if (MH_Process)
        D_Flag = 41;
    F_Flag = H_Process;
    T_Flag = T_Process;
    C_Flag = MASS_TRANSPORT_Process;
    PT_Flag = 0;  // PCH Initialize to be no RWPT.
    RD_Flag = RD_Process;
    MCF_Flag = MULTI_COMPONENTIAL_FLOW_Process;

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
    idxm = NULL;        //> global indices of local matrix rows
    idxn = NULL;        //> global indices of local matrix columns
    local_idx = NULL;   //> local index for local assemble
// local_matrix = NULL; //>  local matrix
// local_vec = NULL; //>  local vector
#endif
    element_nodes_dom = NULL;
    Index = 0;
    nnodes = nnodesHQ = nNodes = 0;
    Radius = .0;
    Xi_p = .0;
}

//  Destructor of class Element
CElement::~CElement()
{
    delete[] _Jacobian;
    delete[] _determinants_all;
    delete[] _inv_jacobian_all;
    delete[] _dshapefct_all;
    delete[] _dshapefctHQ_all;

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
    if (idxm)
        delete[] idxm;
    if (idxn)
        delete[] idxn;
    if (local_idx)
        delete[] local_idx;
// if (local_idx)
//  delete [] local_matrix;
// if (local_idx)
//  delete [] local_vec;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   05/2007 WW 1D in 2D
   Last modified:
**************************************************************************/
void CElement::ConfigElement(CElem* MElement, const bool FaceIntegration)
{
    CNode* a_node = NULL;  // 07.04.2009. WW
    MeshElement = MElement;
    Index = MeshElement->GetIndex();
    nnodes = MeshElement->nnodes;
    nnodesHQ = MeshElement->nnodesHQ;
    ele_dim = MeshElement->GetDimension();
    _ele_global_dim = ele_dim;
    bool done = false;
    if (MeshElement->quadratic)
        nNodes = nnodesHQ;
    else
        nNodes = nnodes;
    // Node indices
    for (int i = 0; i < nNodes; i++)
        nodes[i] = MeshElement->nodes_index[i];
    // Put coordinates of nodes to buffer to enhance the computation
    if (!FaceIntegration)
    {
        if (dim != ele_dim)
        {
            _ele_global_dim = dim;
            //            a_node0 = MeshElement->nodes[0];      //07.04.2007. WW
            double const* const coords_node_0(MeshElement->nodes[0]->getData());
            for (int i = 0; i < nNodes; i++)
            {
                double const* const coords_node_i(
                    MeshElement->nodes[i]->getData());
                //               a_node = MeshElement->nodes[i]; //07.04.2007.
                //               WW dy = dz = 0.; dx = a_node->X()-a_node0->X();
                //               dy = a_node->Y()-a_node0->Y();
                //               dz = a_node->Z()-a_node0->Z();
                double dx(coords_node_i[0] - coords_node_0[0]);
                double dy(coords_node_i[1] - coords_node_0[1]);
                double dz(coords_node_i[2] - coords_node_0[2]);

                X[i] = (*MeshElement->transform_tensor)(0, 0) * dx +
                       (*MeshElement->transform_tensor)(1, 0) * dy +
                       (*MeshElement->transform_tensor)(2, 0) * dz;
                Y[i] = (*MeshElement->transform_tensor)(0, 1) * dx +
                       (*MeshElement->transform_tensor)(1, 1) * dy +
                       (*MeshElement->transform_tensor)(2, 1) * dz;
                //               Z[i] =  a_node->Z();
                Z[i] = coords_node_i[2];
            }
            done = true;
        }
        else
        {
            switch (dim)
            {
                case 1:
                    if (coordinate_system % 10 == 1)
                    {
                        for (int i = 0; i < nNodes; i++)
                        {
                            // 07.04.2007. WW
                            //                        a_node =
                            //                        MeshElement->nodes[i];
                            //                        X[i] = a_node->Y();
                            //                        Y[i] = a_node->X();
                            //                        Z[i] = a_node->Z();
                            double const* const coords_node_i(
                                MeshElement->nodes[i]->getData());
                            X[i] = coords_node_i[1];
                            Y[i] = coords_node_i[0];
                            Z[i] = coords_node_i[2];
                        }
                        done = true;
                    }
                    else if (coordinate_system % 10 == 2)
                    {
                        for (int i = 0; i < nNodes; i++)
                        {
                            // 07.04.2007. WW
                            //                        a_node =
                            //                        MeshElement->nodes[i];
                            //                        X[i] = a_node->Z();
                            //                        Y[i] = a_node->Y();
                            //                        Z[i] = a_node->X();
                            double const* const coords_node_i(
                                MeshElement->nodes[i]->getData());
                            X[i] = coords_node_i[2];
                            Y[i] = coords_node_i[1];
                            Z[i] = coords_node_i[0];
                        }
                        done = true;
                    }
                    break;
                case 2:
                    if (coordinate_system % 10 == 2)
                    {
                        for (int i = 0; i < nNodes; i++)
                        {
                            double const* const coords_node_i(
                                MeshElement->nodes[i]->getData());
                            X[i] = coords_node_i[0];
                            Y[i] = coords_node_i[2];
                            Z[i] = coords_node_i[1];
                        }
                        done = true;
                    }
                    break;
            }
        }
    }
    //
    if (!done)
    {
        for (int i = 0; i < nNodes; i++)
        {
            a_node = MeshElement->nodes[i];  // 07.04.2007. WW
            double const* const coords(a_node->getData());
            X[i] = coords[0];
            Y[i] = coords[1];
            Z[i] = coords[2];
        }
    }
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
    if (!FaceIntegration)
    {
        if (MeshElement
                ->g_index)  // ghost nodes pcs->pcs_number_of_primary_nvals
        {
            act_nodes = MeshElement->g_index[0];
            act_nodes_h = MeshElement->g_index[1];

            for (int i = 0; i < act_nodes_h; i++)
            {
                local_idx[i] = MeshElement->g_index[i + 2];
            }
        }
        else
        {
            act_nodes = nnodes;
            act_nodes_h = nnodesHQ;
            for (int i = 0; i < act_nodes_h; i++)
            {
                local_idx[i] = i;
            }
        }
    }
#endif

    // WW
    const MshElemType::type elem_type = MeshElement->GetElementType();
    getShapeFunctionPtr(elem_type);
    getGradShapeFunctionPtr(elem_type);

    SetIntegrationPointNumber(elem_type);
    if (_is_mixed_order)
    {
        Order = 1;
    }
    ComputeGradShapefctInElement(FaceIntegration);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   Last modified:
**************************************************************************/
void CElement::calculateRadius(const int gp)
{
    Radius = 0.0;
    getShapefunctValues(gp, 1);
    for (int i = 0; i < nnodes; i++)
        Radius += shapefct[i] * X[i];
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2005 WW Implementation
   Last modified:
**************************************************************************/
void CElement::setOrder(const int order)
{
    Order = order;
    if (Order == 1)
        nNodes = nnodes;
    else if (Order == 2)
        nNodes = nnodesHQ;
}

void CElement::SetIntegrationPointNumber(const MshElemType::type elem_type)
{
    switch (elem_type)
    {
        case MshElemType::LINE:
            // nGauss = 2;
            nGaussPoints = nGauss;
            return;
        case MshElemType::QUAD:
        case MshElemType::QUAD8:
            nGaussPoints = nGauss * nGauss;
            return;
        case MshElemType::HEXAHEDRON:
            nGaussPoints = nGauss * nGauss * nGauss;
            return;
        case MshElemType::TRIANGLE:
            nGaussPoints = 3;  // nGauss = 3; // Fixed to 3
            return;
        case MshElemType::TETRAHEDRON:
            //	   nGaussPoints = nGauss = 15;  // Fixed to 15
            nGaussPoints = 5;  // nGauss = 5; // Fixed to 5
            return;
        case MshElemType::PRISM:
            nGaussPoints = 6;  // Fixed to 6
            // nGauss = 3;               // Fixed to 3
            return;
        case MshElemType::PYRAMID:
            if (Order == 1)
                nGaussPoints = 5;  // nGauss = 5;
            else
                nGaussPoints = 8;  // nGauss = 8;  //13;
            return;
        case MshElemType::INVALID:
            std::cerr << "[CElement::ConfigNumerics] invalid element type"
                      << std::endl;
            break;
        default:
            std::cerr << "[CElement::ConfigNumerics] unknown element type"
                      << std::endl;
            break;
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   02/2005 OK Case 1: line elements
   01/2010 NW Higher order line elements
   Last modified:
**************************************************************************/
void CElement::ConfigShapefunction(MshElemType::type elem_type)
{
    SetIntegrationPointNumber(elem_type);

    switch (elem_type)
    {
        case MshElemType::LINE:
            ele_dim = 1;
            // nGauss = 2;
            nGaussPoints = nGauss;
            ShapeFunction = ShapeFunctionLine;
            ShapeFunctionHQ = ShapeFunctionLineHQ;
            GradShapeFunction = GradShapeFunctionLine;
            GradShapeFunctionHQ = GradShapeFunctionLineHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::QUAD:
            ele_dim = 2;
            nGaussPoints = nGauss * nGauss;
            ShapeFunction = ShapeFunctionQuad;
            ShapeFunctionHQ = ShapeFunctionQuadHQ;
            GradShapeFunction = GradShapeFunctionQuad;
            GradShapeFunctionHQ = GradShapeFunctionQuadHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::QUAD8:
            ele_dim = 2;
            nGaussPoints = nGauss * nGauss;
            ShapeFunction = ShapeFunctionQuad;
            ShapeFunctionHQ = ShapeFunctionQuadHQ8;
            GradShapeFunction = GradShapeFunctionQuad;
            GradShapeFunctionHQ = GradShapeFunctionQuadHQ8;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::HEXAHEDRON:
            ele_dim = 3;
            nGaussPoints = nGauss * nGauss * nGauss;
            ShapeFunction = ShapeFunctionHex;
            ShapeFunctionHQ = ShapeFunctionHexHQ;
            GradShapeFunction = GradShapeFunctionHex;
            GradShapeFunctionHQ = GradShapeFunctionHexHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::TRIANGLE:
            ele_dim = 2;
            nGaussPoints = 3;  // nGauss = 3; // Fixed to 3
            ShapeFunction = ShapeFunctionTri;
            ShapeFunctionHQ = ShapeFunctionTriHQ;
            GradShapeFunction = GradShapeFunctionTri;
            GradShapeFunctionHQ = GradShapeFunctionTriHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::TETRAHEDRON:
            ele_dim = 3;
            //	   nGaussPoints = nGauss = 15;  // Fixed to 15
            nGaussPoints = 5;  // nGauss = 5; // Fixed to 5
            ShapeFunction = ShapeFunctionTet;
            ShapeFunctionHQ = ShapeFunctionTetHQ;
            GradShapeFunction = GradShapeFunctionTet;
            GradShapeFunctionHQ = GradShapeFunctionTetHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_LINEAR;
            return;
        case MshElemType::PRISM:
            ele_dim = 3;
            nGaussPoints = 6;  // Fixed to 6
            // nGauss = 3;               // Fixed to 3
            ShapeFunction = ShapeFunctionPri;
            ShapeFunctionHQ = ShapeFunctionPriHQ;
            GradShapeFunction = GradShapeFunctionPri;
            GradShapeFunctionHQ = GradShapeFunctionPriHQ;
            extrapo_method = ExtrapolationMethod::EXTRAPO_AVERAGE;
            return;
        case MshElemType::PYRAMID:
            ele_dim = 3;
            if (Order == 1)
                nGaussPoints = 5;  // nGauss = 5;
            else
                nGaussPoints = 8;  // nGauss = 8;  //13;
            ShapeFunction = ShapeFunctionPyra;
            ShapeFunctionHQ = ShapeFunctionPyraHQ13;
            GradShapeFunction = GradShapeFunctionPyra;
            GradShapeFunctionHQ = GradShapeFunctionPyraHQ13;
            extrapo_method = ExtrapolationMethod::EXTRAPO_AVERAGE;
            return;
        case MshElemType::INVALID:
            std::cerr << "[CElement::ConfigNumerics] invalid element type"
                      << std::endl;
            break;
        default:
            std::cerr << "[CElement::ConfigNumerics] unknown element type"
                      << std::endl;
            break;
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2004 WW Implementation
   Last modified:
**************************************************************************/
double CElement::interpolate(double const* const nodalVal,
                             const int order) const
{
    const int nn = (order == 2) ? nnodesHQ : nnodes;
    double const* const inTerpo = (order == 2) ? shapefctHQ : shapefct;
    double val = 0.0;
    for (int i = 0; i < nn; i++)
        val += nodalVal[i] * inTerpo[i];
    return val;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   Last modified:
**************************************************************************/
double CElement::interpolate(const int idx, CRFProcess* m_pcs, const int order)
{
    int i;
    const int nn = (order == 2) ? nnodesHQ : nnodes;
    const double* inTerpo = (order == 2) ? shapefctHQ : shapefct;
    double val = 0.0;

    //
    for (i = 0; i < nn; i++)
        node_val[i] = m_pcs->GetNodeValue(nodes[i], idx);
    for (int i = 0; i < nn; i++)
        val += node_val[i] * inTerpo[i];
    return val;
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   Last modified:
**************************************************************************/
// double CElement::elemnt_average (const int idx, const int order)
double CElement::elemnt_average(const int idx, CRFProcess* m_pcs,
                                const int order)
{
    int i;
    int nn = nnodes;
    double val = 0.0;
    // WW    double* inTerpo = shapefct;
    if (order == 2)
        nn = nnodes;
    // WW       inTerpo = shapefctHQ;
    //
    for (i = 0; i < nn; i++)
        node_val[i] = m_pcs->GetNodeValue(nodes[i], idx);
    return val / (double)nn;
}

/**************************************************************************
   The generalized Jacobian caculation

   Arguments:
    const double *unit:           Local coordiantes
   return
        The determinate of Jacobian
   Programmaenderungen:
   06/2006     WW
   02/2005 OK case 1, line elements
   01/2006 WW Axisymmtry
   09/2006 WW 1D element in 3D
**************************************************************************/
void CElement::computeJacobian(const int gp, const int order,
                               const bool need_inverse)
{
    const int nodes_number = (order == 2) ? nnodesHQ : nnodes;
    const double* dN = (order == 2) ? dshapefctHQ : dshapefct;

    for (size_t i = 0; i < ele_dim * ele_dim; i++)
        _Jacobian[i] = 0.0;
    //--------------------------------------------------------------------
    double DetJac = 0.0;
    switch (ele_dim)
    {
        //................................................................
        case 1:
        {
            // If Line in X or Z direction, coordinate is saved in local X
            // If Line in 3D space, a transform is applied and cast coordinate
            // in local X
            const double dx = X[1] - X[0];  //+Y[1]-Y[0];
            _Jacobian[0] = 0.5 * dx;

            if (!need_inverse)
                return;

            invJacobian = &_inv_jacobian_all[gp * 1];
            invJacobian[0] = 2.0 / dx;
            DetJac = _Jacobian[0];
            // WW
            // if(MeshElement->area>0)
            // DetJac *= MeshElement->area;// Moved to
            // CFiniteElementStd::setMaterial() WW
            // DetJac*=MeshElement->GetFluxArea();//CMCD
            if (axisymmetry)
            {
                calculateRadius(gp);
                DetJac *= Radius;  // 2.0*pai*Radius;
            }
        }
        break;
        //................................................................
        case 2:
            for (int i = 0, j = nodes_number; i < nodes_number; i++, j++)
            {
                _Jacobian[0] += X[i] * dN[i];
                _Jacobian[1] += Y[i] * dN[i];
                _Jacobian[2] += X[i] * dN[j];
                _Jacobian[3] += Y[i] * dN[j];
            }

            if (!need_inverse)
                return;

            invJacobian = &_inv_jacobian_all[gp * 4];
            DetJac = _Jacobian[0] * _Jacobian[3] - _Jacobian[1] * _Jacobian[2];
            if (fabs(DetJac) < MKleinsteZahl)
            {
                std::cout << "\n*** Jacobian: Det == 0 " << DetJac
                          << " in element " << MeshElement->GetIndex() << "\n";
                abort();
            }
            invJacobian[0] = _Jacobian[3];
            invJacobian[1] = -_Jacobian[1];
            invJacobian[2] = -_Jacobian[2];
            invJacobian[3] = _Jacobian[0];
            for (size_t i = 0; i < ele_dim * ele_dim; i++)
                invJacobian[i] /= DetJac;
            //
            // By WW
            // if(MeshElement->area>0)
            // DetJac *= MeshElement->area;// Moved to
            // CFiniteElementStd::setMaterial() WW
            // DetJac*=MeshElement->GetFluxArea();//CMCD
            if (axisymmetry)
            {
                calculateRadius(gp);
                DetJac *= Radius;  // 2.0*pai*Radius;
            }
            break;
        //................................................................
        case 3:
        {
            for (int i = 0; i < nodes_number; i++)
            {
                const int j = i + nodes_number;
                const int k = i + 2 * nodes_number;

                _Jacobian[0] += X[i] * dN[i];
                _Jacobian[1] += Y[i] * dN[i];
                _Jacobian[2] += Z[i] * dN[i];

                _Jacobian[3] += X[i] * dN[j];
                _Jacobian[4] += Y[i] * dN[j];
                _Jacobian[5] += Z[i] * dN[j];

                _Jacobian[6] += X[i] * dN[k];
                _Jacobian[7] += Y[i] * dN[k];
                _Jacobian[8] += Z[i] * dN[k];
            }

            if (!need_inverse)
                return;

            invJacobian = &_inv_jacobian_all[gp * 9];
            DetJac = _Jacobian[0] * (_Jacobian[4] * _Jacobian[8] -
                                     _Jacobian[7] * _Jacobian[5]) +
                     _Jacobian[6] * (_Jacobian[1] * _Jacobian[5] -
                                     _Jacobian[4] * _Jacobian[2]) +
                     _Jacobian[3] * (_Jacobian[2] * _Jacobian[7] -
                                     _Jacobian[8] * _Jacobian[1]);

            if (fabs(DetJac) < MKleinsteZahl)
            {
                std::cout << "\n*** Jacobian: DetJac == 0 " << DetJac
                          << " in element " << MeshElement->GetIndex() << "\n";
                abort();
            }
            invJacobian[0] =
                _Jacobian[4] * _Jacobian[8] - _Jacobian[7] * _Jacobian[5];
            invJacobian[1] =
                _Jacobian[2] * _Jacobian[7] - _Jacobian[1] * _Jacobian[8];
            invJacobian[2] =
                _Jacobian[1] * _Jacobian[5] - _Jacobian[2] * _Jacobian[4];
            //
            invJacobian[3] =
                _Jacobian[5] * _Jacobian[6] - _Jacobian[8] * _Jacobian[3];
            invJacobian[4] =
                _Jacobian[0] * _Jacobian[8] - _Jacobian[6] * _Jacobian[2];
            invJacobian[5] =
                _Jacobian[2] * _Jacobian[3] - _Jacobian[5] * _Jacobian[0];
            //
            invJacobian[6] =
                _Jacobian[3] * _Jacobian[7] - _Jacobian[6] * _Jacobian[4];
            invJacobian[7] =
                _Jacobian[1] * _Jacobian[6] - _Jacobian[7] * _Jacobian[0];
            invJacobian[8] =
                _Jacobian[0] * _Jacobian[4] - _Jacobian[3] * _Jacobian[1];
            for (size_t i = 0; i < ele_dim * ele_dim; i++)
                invJacobian[i] /= DetJac;
            break;
        }  // end case 3
    }
    //--------------------------------------------------------------------
    // Use absolute value (for grids by gmsh, whose orientation is clockwise)
    _determinants_all[gp] = fabs(DetJac);
}
/***************************************************************************
   GeoSys - Funktion: CElement::RealCoordinates

   Aufgabe:
        Mapping to real coordaintes from the local ones of quadratic traingle
   element.
   Formalparameter:
           E:
             double * x         : Array of size 3, real coordiantes
             const double *u    : Array of size 2, unit coordiantes

   Programming:
   06/2003     WW        Erste Version
   07/2005     WW        Change due to geometry element object
 **************************************************************************/
void CElement::RealCoordinates(double* realXYZ)
{
    int i;
    const double* df = (Order == 2) ? shapefctHQ : shapefct;

    for (i = 0; i < 3; i++)
        realXYZ[i] = 0.0;

    for (i = 0; i < nNodes; i++)
    {
        realXYZ[0] += df[i] * X[i];
        realXYZ[1] += df[i] * Y[i];
        realXYZ[2] += df[i] * Z[i];
    }
}
/***************************************************************************
   GeoSys - Funktion: CElement::UnitCoordinates

   Aufgabe:
        Get unit coodinates from the real ones
   element.
   Formalparameter:
           E:
             double * x         : Array of size 3, real coordiantes

   Programming:
   06/2003     WW        Erste Version
 **************************************************************************/
void CElement::UnitCoordinates(double* realXYZ)
{
    setOrder(Order);

    x1buff[0] = X[0];
    x1buff[1] = Y[0];
    x1buff[2] = Z[0];

    for (int i = 1; i < nNodes; i++)
    {
        x1buff[0] += X[i];
        x1buff[1] += Y[i];
        x1buff[2] += Z[i];
    }
    for (size_t i = 0; i < 3; i++)
        x1buff[i] /= (double)nNodes;

    for (size_t i = 0; i < ele_dim; i++)
        realXYZ[i] -= x1buff[i];

    for (size_t i = 0; i < ele_dim; i++)
    {
        unit[i] = 0.0;
        for (size_t j = 0; j < ele_dim; j++)
            unit[i] += invJacobian[j * ele_dim + i] * realXYZ[j];
    }

    for (size_t i = 0; i < ele_dim; i++)
        realXYZ[i] = unit[i];
}

void CElement::SetGaussPoint(const int gp, int& gp_r, int& gp_s, int& gp_t)
{
    SetGaussPoint(MeshElement->GetElementType(), gp, gp_r, gp_s, gp_t);
}

/***************************************************************************

   08/2005     WW        Prism element
 **************************************************************************/
void CElement::SetGaussPoint(const MshElemType::type elem_type, const int gp,
                             int& gp_r, int& gp_s, int& gp_t)
{
    switch (elem_type)
    {
        case MshElemType::LINE:  // Line
            gp_r = gp;
            unit[0] = MXPGaussPkt(nGauss, gp_r);
            return;
        case MshElemType::QUAD8:  // Quadralateral
        case MshElemType::QUAD:   // Quadralateral
            gp_r = (int)(gp / nGauss);
            gp_s = gp % nGauss;
            unit[0] = MXPGaussPkt(nGauss, gp_r);
            unit[1] = MXPGaussPkt(nGauss, gp_s);
            return;
        case MshElemType::HEXAHEDRON:  // Hexahedra
            gp_r = (int)(gp / (nGauss * nGauss));
            gp_s = (gp % (nGauss * nGauss));
            gp_t = gp_s % nGauss;
            gp_s /= nGauss;
            unit[0] = MXPGaussPkt(nGauss, gp_r);
            unit[1] = MXPGaussPkt(nGauss, gp_s);
            unit[2] = MXPGaussPkt(nGauss, gp_t);
            return;
        case MshElemType::TRIANGLE:  // Triangle
            SamplePointTriHQ(gp, unit);
            break;
        case MshElemType::TETRAHEDRON:  // Tedrahedra
            // To be flexible          SamplePointTet15(gp, unit);
            SamplePointTet5(gp, unit);
            return;
        case MshElemType::PRISM:  // Prism
            // nGaussPoints = 6;         // Fixed to 6
            // nGauss = 3;               // Fixed to 3
            gp_r = gp % 3;  // gp % nGauss
            SamplePointTriHQ(gp_r, unit);
            //
            gp_s = nGaussPoints / 3;  // nGaussPoints/nGauss
            gp_t = (int)(gp / 3);     // gp/nGauss
            unit[2] = MXPGaussPkt(gp_s, gp_t);
            return;
        case MshElemType::PYRAMID:  // Pyramid
            if (Order == 1)
                SamplePointPyramid5(gp, unit);
            else
                SamplePointPyramid8(gp,
                                    unit);  // SamplePointPyramid13(gp, unit);
            return;
        default:
            std::cerr
                << "CElement::SetGaussPoint invalid mesh element type given"
                << std::endl;
            break;
    }
}
/***************************************************************************
   GeoSys - Funktion:
           CElement:: GetGaussData(const int gp)

   Aufgabe:
          Get Gauss points and weights, compute Jacobian
   Formalparameter:
           E:
             const int gp   : Gauss point index

   Programming:
   06/2004     WW        Erste Version
   08/2005     WW        Prism element
   02/2005 OK case 1
   02/2007 WW Abstract the calcultion of Gauss point in one function
 **************************************************************************/
double CElement::GetGaussData(int gp, int& gp_r, int& gp_s, int& gp_t)
{
    switch (MeshElement->GetElementType())
    {
        case MshElemType::LINE:  // Line
            gp_r = gp;
            return _determinants_all[gp] * MXPGaussFkt(nGauss, gp_r);
        case MshElemType::QUAD8:  // Quadralateral
        case MshElemType::QUAD:   // Quadralateral
            gp_r = (int)(gp / nGauss);
            gp_s = gp % nGauss;
            return _determinants_all[gp] * MXPGaussFkt(nGauss, gp_r) *
                   MXPGaussFkt(nGauss, gp_s);
        case MshElemType::HEXAHEDRON:  // Hexahedra
            gp_r = (int)(gp / (nGauss * nGauss));
            gp_s = (gp % (nGauss * nGauss));
            gp_t = gp_s % nGauss;
            gp_s /= nGauss;
            return _determinants_all[gp] * MXPGaussFkt(nGauss, gp_r) *
                   MXPGaussFkt(nGauss, gp_s) * MXPGaussFkt(nGauss, gp_t);
        case MshElemType::TRIANGLE:  // Triangle
            SamplePointTriHQ(gp, unit);
            return _determinants_all[gp] * unit[2];  // Weights
        case MshElemType::TETRAHEDRON:               // Tedrahedra
            // To be flexible          SamplePointTet15(gp, unit);
            SamplePointTet5(gp, unit);
            return _determinants_all[gp] * unit[3];  // Weights
        case MshElemType::PRISM:                     // Prism nGauss=3
            gp_r = gp % 3;                           // gp % nGauss
            //
            gp_s = nGaussPoints / 3;  // nGaussPoints/nGauss
            gp_t = (int)(gp / 3);     // gp/nGauss
            return _determinants_all[gp] * MXPGaussFktTri(3, gp_r) *
                   MXPGaussFkt(gp_s, gp_t);
        case MshElemType::PYRAMID:  // Pyramid
            if (Order == 1)
                SamplePointPyramid5(gp, unit);
            else
                SamplePointPyramid8(gp,
                                    unit);  // SamplePointPyramid13(gp, unit);
            return _determinants_all[gp] * unit[3];  // Weights
        default:
            std::cerr
                << "CElement::GetGaussData invalid mesh element type given"
                << std::endl;
            return 0.;
    }
    return 0.;
}

/***************************************************************************
   GeoSys - Funktion: FaceIntegration(const double *NodeVal)
   Task:   Used to treat Nuemann type boundary conditions (3D)
   Augument
       double *NodeVal : input, values of boundary conditions at all face node
                         Output, integration results of all face nodes
   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CElement::FaceIntegration(double* NodeVal)
{
    setOrder(Order);

    getShapeFunctionPtr(MeshElement->GetElementType());

    for (int i = 0; i < nNodes; i++)
        dbuff[i] = 0.0;
    // Loop over Gauss points
    for (gp = 0; gp < nGaussPoints; gp++)
    {
        int gp_r, gp_s, gp_t;
        double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

        getShapefunctValues(gp, Order);
        const double* sf = (Order == 1) ? shapefct : shapefctHQ;

        if (this->axisymmetry)
        {
            calculateRadius(gp);
            fkt *= Radius;  // 2.0*pai*radius;
        }

        double val = 0.0;
        // Interpolation of value at Gauss point
        for (int i = 0; i < nNodes; i++)
            val += NodeVal[i] * sf[i];
        // Integration
        for (int i = 0; i < nNodes; i++)
            dbuff[i] += val * sf[i] * fkt;
    }
    for (int i = 0; i < nNodes; i++)
        NodeVal[i] = dbuff[i];
}

void CElement::DomainIntegration(double* NodeVal)
{
    setOrder(Order);

    getShapeFunctionPtr(MeshElement->GetElementType());

    // double det = MeshElement->GetVolume();
    for (int i = 0; i < nNodes; i++)
        dbuff[i] = 0.0;
    // Loop over Gauss points
    for (int gp = 0; gp < nGaussPoints; gp++)
    {
        //---------------------------------------------------------
        //  Get local coordinates and weights
        //---------------------------------------------------------
        int gp_r, gp_s, gp_t;
        double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
        getShapefunctValues(gp, Order);
        if (this->axisymmetry)
        {
            calculateRadius(gp);
            fkt *= Radius;  // 2.0*pai*radius;
        }
        const double* sf = (Order == 1) ? shapefct : shapefctHQ;

        double val = 0.0;
        // Interpolation of value at Gauss point
        for (int i = 0; i < nNodes; i++)
            val += NodeVal[i] * sf[i];
        // Integration
        for (int i = 0; i < nNodes; i++)
            dbuff[i] += val * sf[i] * fkt;
    }
    for (int i = 0; i < nNodes; i++)
        NodeVal[i] = dbuff[i];
}

void CElement::getShapefunctValues(const int gp, const int order) const
{
    if (order == 1)
    {
        assert(static_cast<unsigned>(gp * nnodes) <
               _shape_function_pool_ptr[0]->getShapeFunctionArraySize(
                   MeshElement->GetElementType()));

        shapefct = &_shape_function_result_ptr[0][nnodes * gp];
    }
    else if (order == 2)
    {
        assert(static_cast<unsigned>(gp * nnodesHQ) <
               _shape_function_pool_ptr[1]->getShapeFunctionArraySize(
                   MeshElement->GetElementType()));

        shapefctHQ = &_shape_function_result_ptr[1][nnodesHQ * gp];
    }
}

void CElement::ComputeShapefct(const int order, double shape_function[])
{
    if (order == 1)
        ShapeFunction(shape_function, unit);
    else if (order == 2)
        ShapeFunctionHQ(shape_function, unit);
}

void CElement::computeGradShapefctLocal(const int order,
                                        double grad_shape_fucntion[])
{
    if (order == 1)
        GradShapeFunction(grad_shape_fucntion, unit);
    else if (order == 2)
        GradShapeFunctionHQ(grad_shape_fucntion, unit);
}

void CElement::getShapeFunctionCentroid()
{
    if (_shape_function_pool_ptr[0])
        shapefct =
            (_shape_function_pool_ptr[0])
                ->getShapeFunctionCenterValues(MeshElement->GetElementType());
    if (_shape_function_pool_ptr[1])
        shapefctHQ =
            (_shape_function_pool_ptr[1])
                ->getShapeFunctionCenterValues(MeshElement->GetElementType());
}

void CElement::getGradShapeFunctionCentroid()
{
    if (_shape_function_pool_ptr[0])
        dshapefct = (_shape_function_pool_ptr[0])
                        ->getGradShapeFunctionCenterValues(
                            MeshElement->GetElementType());
    if (_shape_function_pool_ptr[1])
        dshapefctHQ = (_shape_function_pool_ptr[1])
                          ->getGradShapeFunctionCenterValues(
                              MeshElement->GetElementType());
}

void CElement::getShapeFunctionPtr(const MshElemType::type elem_type)
{
    if (_shape_function_pool_ptr[0])
        _shape_function_result_ptr[0] =
            (_shape_function_pool_ptr[0])->getShapeFunctionValues(elem_type);
    if (_shape_function_pool_ptr[1])
        _shape_function_result_ptr[1] =
            (_shape_function_pool_ptr[1])->getShapeFunctionValues(elem_type);
}

void CElement::getGradShapeFunctionPtr(const MshElemType::type elem_type)
{
    if (_shape_function_pool_ptr[0])
        _grad_shape_function_result_ptr[0] =
            (_shape_function_pool_ptr[0])
                ->getGradShapeFunctionValues(elem_type);
    if (_shape_function_pool_ptr[1])
        _grad_shape_function_result_ptr[1] =
            (_shape_function_pool_ptr[1])
                ->getGradShapeFunctionValues(elem_type);
}

void CElement::getGradShapefunctValues(const int gp, const int order) const
{
    if (order == 1)
    {
        dshapefct = &_dshapefct_all[nnodes * _ele_global_dim * gp];
    }
    else
    {
        dshapefctHQ = &_dshapefctHQ_all[nnodesHQ * _ele_global_dim * gp];
    }
}

/***************************************************************************
   GeoSys - Funktion:
           CElement::ComputeGradShapefct(const double *unit, const int order)

   Aufgabe:
         Compute values of shape function at integral point unit.
   Formalparameter:
           E:
             const double *unit    : Array of size 2, unit coordiantes
             const int order    : 1, linear
                               2, quadratic

   Programming:
   06/2004     WW        Erste Version
   10/2005     WW        2D element transform in 3D space
   06/2007     WW        1D in 2D
 **************************************************************************/
void CElement::ComputeGradShapefct(const int gp, const int order,
                                   const bool is_face_integration)
{
    setOrder(order);

    double* dshp_fct = NULL;
    const double* dshp_fct_local = (order == 1) ? dshapefct : dshapefctHQ;

    if (order == 1)
    {
        dshp_fct = &_dshapefct_all[nNodes * _ele_global_dim * gp];
    }
    else
    {
        dshp_fct = &_dshapefctHQ_all[nNodes * _ele_global_dim * gp];
    }

    int j_times_ele_dim_plus_k, j_times_nNodes_plus_i;
    double Var[3];
    for (int i = 0; i < nNodes; i++)
    {
        size_t j(0);
        for (j = 0, j_times_nNodes_plus_i = i; j < ele_dim;
             j++, j_times_nNodes_plus_i += nNodes)
        {
            Var[j] = dshp_fct_local[j_times_nNodes_plus_i];
            dshp_fct[j_times_nNodes_plus_i] = 0.0;
        }
        for (j = 0, j_times_ele_dim_plus_k = 0, j_times_nNodes_plus_i = i;
             j < ele_dim; j++, j_times_nNodes_plus_i += nNodes)
        {
            for (size_t k = 0; k < ele_dim; k++, j_times_ele_dim_plus_k++)
            {
                dshp_fct[j_times_nNodes_plus_i] +=
                    invJacobian[j_times_ele_dim_plus_k] * Var[k];
            }
        }
    }

    if (is_face_integration)
        return;

    // 1D element in 3D
    if ((dim == 3 && ele_dim == 1) || (dim == 2 && ele_dim == 1))
        for (int i = 0; i < nNodes; i++)
        {
            for (size_t j = 1; j < dim; j++)
                dshp_fct[j * nNodes + i] =
                    (*MeshElement->transform_tensor)(j)*dshp_fct[i];
            dshp_fct[i] *= (*MeshElement->transform_tensor)(0);
        }
    // 2D element in 3D
    if (dim == 3 && ele_dim == 2)
    {
        const size_t n_nodes_times_ele_dim(nNodes * ele_dim);
        for (size_t i = 0; i < n_nodes_times_ele_dim; i++)
            dbuff0[i] = dshp_fct[i];
        for (int i = 0; i < nNodes; i++)
            for (size_t j = 0; j < dim; j++)
            {
                dshp_fct[j * nNodes + i] = 0.0;
                for (size_t k = 0; k < ele_dim; k++)
                    dshp_fct[j * nNodes + i] +=
                        (*MeshElement->transform_tensor)(j, k) *
                        dbuff0[k * nNodes + i];
            }
    }
}

void CElement::getLocalGradShapefunctValues(const int gp, const int order)
{
    if (order == 1)
    {
        dshapefct = &_grad_shape_function_result_ptr[0][nnodes * ele_dim * gp];
    }
    else if (order == 2)
    {
        dshapefctHQ =
            &_grad_shape_function_result_ptr[1][nnodesHQ * ele_dim * gp];
    }
}

void CElement::ComputeGradShapefctInElement(const bool is_face_integration)
{
    for (int gp = 0; gp < nGaussPoints; gp++)
    {
        getLocalGradShapefunctValues(gp, Order);
        computeJacobian(gp, Order);
        ComputeGradShapefct(gp, Order, is_face_integration);
    }
}

/***************************************************************************
   Center of reference element
   Programming:
   09/2005     WW        Erste Version
 **************************************************************************/
void CElement::SetCenterGP(const MshElemType::type elem_type)
{
    const MshElemType::type e_type = (elem_type == MshElemType::INVALID)
                                         ? MeshElement->GetElementType()
                                         : elem_type;
    // Center of the reference element
    unit[0] = unit[1] = unit[2] = 0.0;
    if (e_type == MshElemType::TRIANGLE)
        unit[0] = unit[1] = 1.0 / 3.0;
    else if (e_type == MshElemType::TETRAHEDRON)
        unit[0] = unit[1] = unit[2] = 0.25;
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::GetLocalIndex()
           For quadralateral and hexahedra element on the assumption that
           selected Gauss points form a quadralateral or hexahedra
   Aufgabe:
           Accumulate stress at each nodes
   Formalparameter:

   Programming:
   06/2004   WW
 **************************************************************************/
int CElement::GetLocalIndex(const int gp_r, const int gp_s, int gp_t)
{
    int LoIndex = -1;
    double r, s, t;

    //---------------------------------------------------------
    // Accumulate strains
    //---------------------------------------------------------
    switch (MeshElement->GetElementType())
    {
        case MshElemType::QUAD:  // Quadralateral
            r = MXPGaussPkt(nGauss, gp_r);
            s = MXPGaussPkt(nGauss, gp_s);
            if (r > 0.0 && s > 0.0)
                LoIndex = 0;
            else if (r < 0.0 && s > 0.0)
                LoIndex = 1;
            else if (r < 0.0 && s < 0.0)
                LoIndex = 2;
            else if (r > 0.0 && s < 0.0)
                LoIndex = 3;
            else if (fabs(r) < MKleinsteZahl && s > 0.0)
                LoIndex = 4;
            else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 5;
            else if (fabs(r) < MKleinsteZahl && s < 0.0)
                LoIndex = 6;
            else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                LoIndex = 7;
            else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                LoIndex = 8;
            break;
        case MshElemType::HEXAHEDRON:  // Hexahedra
            r = MXPGaussPkt(nGauss, gp_r);
            s = MXPGaussPkt(nGauss, gp_s);
            t = MXPGaussPkt(nGauss, gp_t);

            if (t > 0.0)
            {
                if (r > 0.0 && s > 0.0)
                    LoIndex = 0;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 1;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 2;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 3;
                else if (fabs(r) < MKleinsteZahl && s > 0.0)
                    LoIndex = 8;
                else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 9;
                else if (fabs(r) < MKleinsteZahl && s < 0.0)
                    LoIndex = 10;
                else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 11;
                else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                    return -1;
            }
            else if (fabs(t) < MKleinsteZahl)
            {
                if (fabs(r) < MKleinsteZahl || fabs(s) < MKleinsteZahl)
                    return -1;
                if (r > 0.0 && s > 0.0)
                    LoIndex = 16;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 17;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 18;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 19;
            }
            if (t < 0.0)
            {
                if (r > 0.0 && s > 0.0)
                    LoIndex = 4;
                else if (r < 0.0 && s > 0.0)
                    LoIndex = 5;
                else if (r < 0.0 && s < 0.0)
                    LoIndex = 6;
                else if (r > 0.0 && s < 0.0)
                    LoIndex = 7;
                else if (fabs(r) < MKleinsteZahl && s > 0.0)
                    LoIndex = 12;
                else if (r < 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 13;
                else if (fabs(r) < MKleinsteZahl && s < 0.0)
                    LoIndex = 14;
                else if (r > 0.0 && fabs(s) < MKleinsteZahl)
                    LoIndex = 15;
                else if (fabs(r) < MKleinsteZahl && fabs(s) < MKleinsteZahl)
                    return -1;
            }
            break;
        default:
            std::cerr
                << "CElement::GetLocalIndex invalid mesh element type given"
                << std::endl;
            break;
    }
    return LoIndex;
}
/***************************************************************************
   GeoSys - Funktion:
   Programming:
   02/2007   WW
 **************************************************************************/
void CElement::SetExtropoGaussPoints(const int i)
{
    int j = 0;
    MshElemType::type ElementType = MeshElement->GetElementType();
    //
    switch (ElementType)
    {
        case MshElemType::TRIANGLE:  // Triangle
            // Compute values at verteces
            // Compute values at verteces
            switch (i)
            {
                case 0:
                    unit[0] = -0.1666666666667;
                    unit[1] = -0.1666666666667;
                    break;
                case 1:
                    unit[0] = 1.6666666666667;
                    unit[1] = -0.1666666666667;
                    break;
                case 2:
                    unit[0] = -0.1666666666667;
                    unit[1] = 1.6666666666667;
                    break;
            }
            break;
        case MshElemType::QUAD8:  // Quadralateral element
        case MshElemType::QUAD:   // Quadralateral element
            // Extropolation over nodes
            switch (i)
            {
                case 0:
                    unit[0] = Xi_p;
                    unit[1] = Xi_p;
                    break;
                case 1:
                    unit[0] = -Xi_p;
                    unit[1] = Xi_p;
                    break;
                case 2:
                    unit[0] = -Xi_p;
                    unit[1] = -Xi_p;
                    break;
                case 3:
                    unit[0] = Xi_p;
                    unit[1] = -Xi_p;
                    break;
            }
            break;
        case MshElemType::HEXAHEDRON:  // Hexahedra
            if (i < 4)
            {
                j = i;
                unit[2] = Xi_p;
            }
            else
            {
                j = i - 4;
                unit[2] = -Xi_p;
            }
            switch (j)
            {
                case 0:
                    unit[0] = Xi_p;
                    unit[1] = Xi_p;
                    break;
                case 1:
                    unit[0] = -Xi_p;
                    unit[1] = Xi_p;
                    break;
                case 2:
                    unit[0] = -Xi_p;
                    unit[1] = -Xi_p;
                    break;
                case 3:
                    unit[0] = Xi_p;
                    unit[1] = -Xi_p;
                    break;
            }
            break;
        case MshElemType::TETRAHEDRON:  // Tedrahedra
            // Compute values at verteces
            switch (i)
            {
                case 0:
                    unit[0] = -0.166666666666667;
                    unit[1] = -0.166666666666667;
                    unit[2] = -0.166666666666667;
                    break;
                case 1:
                    unit[0] = 1.5;
                    unit[1] = -0.166666666666667;
                    unit[2] = -0.166666666666667;
                    break;
                case 2:
                    unit[0] = -0.166666666666667;
                    unit[1] = 1.5;
                    unit[2] = -0.166666666666667;
                    break;
                case 3:
                    unit[0] = -0.166666666666667;
                    unit[1] = -0.166666666666667;
                    unit[2] = 1.5;
                    break;
            }
            break;
        case MshElemType::LINE:
            unit[0] = -Xi_p;
            unit[1] = Xi_p;
            break;
        case MshElemType::PYRAMID:  // WW. 09.2012. WW
            SamplePointPyramid5(i, unit);
            break;
        default:
            unit[0] = unit[1] = unit[2] = 0.;  // 07.01.2011. WW
            break;
    }
}

/***************************************************************************
   GeoSys - Funktion:
   Programming:
   05/2011   NW
 **************************************************************************/
double CElement::CalcXi_p()
{
    MshElemType::type ElementType = MeshElement->GetElementType();
    Xi_p = 0.0;
    if (ElementType == MshElemType::LINE || ElementType == MshElemType::QUAD ||
        ElementType == MshElemType::QUAD8 ||
        ElementType == MshElemType::HEXAHEDRON)
    {
        for (gp = 0; gp < nGauss; gp++)
        {
            const double r = fabs(MXPGaussPkt(nGauss, gp));
            if (r > Xi_p)
                Xi_p = r;
        }
        Xi_p = 1.0 / Xi_p;
        return Xi_p;
    }
    else
        return 0.;
}

/***************************************************************************
   GeoSys - Funktion:
   Programming:
   05/2011   NW Implementation
 **************************************************************************/
double CElement::CalcAverageGaussPointValues(double* GpValues)
{
    // average
    double avg = .0;
    for (int j = 0; j < nGauss; j++)
        avg += GpValues[j];
    avg /= nGauss;

    return avg;
}

/**************************************************************************
   ElementMatrix::AllocateMemory

   Arguments:
    const int EleIndex:  Element index,
    int type          : type
                         used to get element type.
                         type = 0, for Possion type
                         type = 1, for Possion equation with deformation
coupling type = 2, for Navier equation type = 3, for Navier equation with
pressure coupling type = 4, Monlithic scheme of u-p coupling type = 5, Mass
Transport default = 0.

   Programmaenderungen:
   01/2005     WW

**************************************************************************/
void ElementMatrix::AllocateMemory(CElem* ele, int type)
{
    int nnodes, nnodesHQ, dim, size;
    size = 0;
    // The following two lines will be updated when new FEMGEO is ready
    nnodes = ele->GetVertexNumber();
    nnodesHQ = ele->GetNodesNumber_H();
    dim = ele->GetDimension();
    switch (type)
    {
        case 0:  // H || T Process
            Mass = new Matrix(nnodes, nnodes);
            //        Laplace = new SymMatrix(nnodes);
            Laplace = new Matrix(nnodes, nnodes);
            RHS = new Vec(nnodes);
            break;
        case 1:  // HM Partioned scheme, Flow
            Mass = new Matrix(nnodes, nnodes);
            //        Laplace = new SymMatrix(nnodes);
            Laplace = new Matrix(nnodes, nnodes);
            RHS = new Vec(nnodes);
            CouplingB = new Matrix(nnodes, dim * nnodesHQ);
            break;
        case 2:  // M_Process only
            size = dim * nnodesHQ;
            Stiffness = new Matrix(size, size);
            RHS = new Vec(size);
            break;
        case 3:  // MH Partioned scheme, M_Process
            size = dim * nnodesHQ;
            Stiffness = new Matrix(size, size);
            RHS = new Vec(size);
            CouplingA = new Matrix(dim * nnodesHQ, nnodes);
            break;
        case 4:  // HM monothlic scheme
            Mass = new Matrix(nnodes, nnodes);
            //        Laplace = new SymMatrix(nnodes);
            Laplace = new Matrix(nnodes, nnodes);
            size = dim * nnodesHQ;
            Stiffness = new Matrix(size, size);
            RHS = new Vec(size + nnodes);
            CouplingA = new Matrix(dim * nnodesHQ, nnodes);
            CouplingB = new Matrix(nnodes, dim * nnodesHQ);
            break;
        case 5:  // Mass Transport process
            Mass = new Matrix(nnodes, nnodes);
            Laplace = new Matrix(nnodes, nnodes);
            Advection = new Matrix(nnodes, nnodes);
            Storage = new Matrix(nnodes, nnodes);
            Content = new Matrix(nnodes, nnodes);
            RHS = new Vec(nnodes);
            break;
    }
}

/**************************************************************************
   ElementMatrix::ElementMatrix

   Arguments:
      const int EleIndex:  Element index,
                           used to get element type.
   Programmaenderungen:
   01/2005     WW

**************************************************************************/
ElementMatrix::~ElementMatrix()
{
    if (Mass)
        delete Mass;
    if (Laplace)
        delete Laplace;
    if (Advection)
        delete Advection;
    if (Storage)
        delete Storage;
    if (Content)
        delete Content;
    if (RHS)
        delete RHS;
    if (CouplingA)
        delete CouplingA;
    if (CouplingB)
        delete CouplingB;
    Mass = NULL;
    Laplace = NULL;
    Advection = NULL;
    Storage = NULL;
    Content = NULL;
    RHS = NULL;
    CouplingA = NULL;
    CouplingB = NULL;
}

/**************************************************************************
CElement::FaceNormalFluxIntegration

Used for TOTAL_FLUX calculation

Programming:
11/2014   JOD

**************************************************************************/

void CElement::FaceNormalFluxIntegration(long /*element_index*/,
                                         double* NodeVal, double* NodeVal_adv,
                                         int* /*nodesFace*/, CElem* /*face*/,
                                         CRFProcess* m_pcs,
                                         double* normal_vector)
{
    double normal_diff_flux_interpol, normal_adv_flux_interpol;
    double dbuff_adv[10], flux[3];
    // ElementValue* gp_ele = ele_gp_value[element_index];

    setOrder(Order);

    for (int i = 0; i < nNodes; i++)
    {
        dbuff[i] = 0.0;
        dbuff_adv[i] = 0.0;
    }

    // Loop over Gauss points
    for (gp = 0; gp < nGaussPoints; gp++)
    {
        int gp_r, gp_s, gp_t;
        const double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

        getShapefunctValues(gp, Order);
        double const* sf = (Order == 1) ? shapefct : shapefctHQ;

        normal_diff_flux_interpol = 0.0;

        if ((m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW) ||
            (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW))
        {
            for (int i = 0; i < nNodes; i++)  // Darcy flux
                normal_diff_flux_interpol += NodeVal[i] * sf[i];

            for (int i = 0; i < nNodes; i++)  // Integration
                dbuff[i] += normal_diff_flux_interpol * sf[i] * fkt;
        }
        else if ((m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT) ||
                 (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT))
        {
            normal_adv_flux_interpol = 0.0;

            for (int i = 0; i < nNodes; i++)
                normal_adv_flux_interpol += NodeVal_adv[i] * sf[i];

            for (int i = 0; i < nNodes; i++)
            {  // Integration
#ifdef USE_TRANSPORT_FLUX
               // Fick or Fourier diffusion
                for (int l = 0; l < 3; l++)
                    flux[l] = gp_ele->TransportFlux(l, gp);
#endif
                normal_diff_flux_interpol = PointProduction(
                    flux, normal_vector);  //    fabs(PointProduction(flux,
                                           //    normal_vector));
                dbuff[i] += normal_diff_flux_interpol * sf[i] * fkt;
                // advection
                dbuff_adv[i] += normal_adv_flux_interpol * sf[i] * fkt;
            }
        }  // end transport
    }      // end gauss points

    for (int i = 0; i < nNodes; i++)
    {
        NodeVal[i] = dbuff[i];
        NodeVal_adv[i] = dbuff_adv[i];
    }
}

}  // end namespace FiniteElement
