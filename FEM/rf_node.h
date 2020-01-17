/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
FEMLib - Object: Node
Task: class implementation
Programing:
04/2004 OK Implementation
last modified
02/2005 MB node parameter....
**************************************************************************/
#ifndef rf_node_INC
#define rf_node_INC

// C++ STL
//#include <list>
#include <iostream>
//#include <fstream>
//#include <string>
#include <vector>

// FEM
#include "FEMEnums.h"

class CNodeValue
{
public:
    CNodeValue(void);
    ~CNodeValue(void);

    void setProcessDistributionType(FiniteElement::DistributionType distype)
    {
        _node_distype = distype;
    }
    FiniteElement::DistributionType getProcessDistributionType() const
    {
        return _node_distype;
    }
    //
    long geo_node_number;
    long msh_node_number;
    double node_value;
    double node_area;

    //    int node_distype;
    double node_parameterA;
    double node_parameterB;
    double node_parameterC;
    double node_parameterD;
    double node_parameterE;
    int CurveIndex;
    int conditional;
    // std::vector<double>history_value;
    long msh_node_number_conditional;
    // JOD   st-coupling 4.7.10
    std::vector<long> msh_node_numbers_averaging;
    // JOD
    std::vector<double> msh_node_weights_averaging;
    std::string tim_type_name;
    // WW
    void Write(std::ostream& os = std::cout) const;
    void Read(std::istream& is = std::cin);  // WW
    bool check_me;                           // OK
    bool _isConstrainedSTNode;

    std::size_t getSTVectorIndex() const { return _st_vector_index; }
    void setSTVectorIndex(int index) { _st_vector_index = index; }
    std::size_t getSTVectorGroup() const { return _st_vector_group; }
    void setSTVectorGroup(int group) { _st_vector_group = group; }

private:
    FiniteElement::DistributionType _node_distype;
    std::size_t _st_vector_index;
    std::size_t _st_vector_group;
};
#endif
