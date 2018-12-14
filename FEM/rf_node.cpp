/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
FEMLib-Object: CNodeValue
Task: Functions
Programing:
01/2005 OK Implementation
**************************************************************************/

#include "rf_node.h"

/**************************************************************************
FEMLib-Method:
Task: constructor
Programing:
04/2005 OK/WW Implementation
last modification:
**************************************************************************/
CNodeValue::CNodeValue() : _node_distype(FiniteElement::INVALID_DIS_TYPE)
{
    geo_node_number = -1;
    msh_node_number = -1;
    msh_node_number_conditional = -1;
    node_value = 0.;
    node_area = 0.;
    node_parameterA = 0.;
    node_parameterB = 0.;
    node_parameterC = 0.;
    node_parameterD = 0.;
    node_parameterE = 0.;
    CurveIndex = -1;
    conditional = -1;
    check_me = true;  // OK
    _isConstrainedSTNode = false;
}

/**************************************************************************
FEMLib-Method:
Task: destructor
Programing:
04/2005 OK/WW Implementation
last modification:
**************************************************************************/
CNodeValue::~CNodeValue()
{
    check_me = false;  // OK
}

/**************************************************************************
FEMLib-Method:
Task: destructor
Programing:
03/2006 WW Implementation
last modification:
**************************************************************************/
void CNodeValue::Write(std::ostream& os) const
{
    std::string deli = "  ";
    os << geo_node_number << deli;
    os << msh_node_number << deli;
    os << CurveIndex << deli;
    os << node_value << deli;
    /*
    // This is for river flow
    // This writing will be valid for river flow when some
    // of its parameters being moved from CSourceTerm to here
    os<< node_distype <<deli;
    os<< node_area <<deli;
    os<< node_parameterA <<deli;
    os<< node_parameterB <<deli;
    os<< node_parameterC <<deli;
    os<< node_parameterD <<deli;
    os<< node_parameterE <<deli;
    os<< conditional <<deli;
    */
    os << "\n";
}

/**************************************************************************
FEMLib-Method:
Task: destructor
Programing:
03/2006 WW Implementation
last modification:
**************************************************************************/
void CNodeValue::Read(std::istream& is)
{
    is >> geo_node_number;
    is >> msh_node_number;
    is >> CurveIndex;
    is >> node_value;
    /*
    // This is for river flow
    // This writing will be valid for river flow when some
    // of its parameters being moved from CSourceTerm to here
    is>> node_distype ;
    is>> node_area ;
    is>> node_parameterA ;
    is>> node_parameterB ;
    is>> node_parameterC ;
    is>> node_parameterD ;
    is>> node_parameterE ;
    is>> conditional ;
    */
    is >> std::ws;
}
