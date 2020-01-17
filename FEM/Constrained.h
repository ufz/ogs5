/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CONSTRAINED_H_
#define CONSTRAINED_H_

#include "FEMEnums.h"

struct Constrained
{
    double constrainedValue;
    FiniteElement::ProcessType constrainedProcessType;
    FiniteElement::PrimaryVariable constrainedPrimVar;
    ConstrainedType::type constrainedDirection;
    ConstrainedVariable::type constrainedVariable;
    bool _isCompleteConstrained;
    bool _completeConstrainedStateOff;
    std::vector<bool> _constrainedNodes;
    bool _isConstrainedVelStable;
    bool _isSeepageBC;

    Constrained()
        : constrainedValue(0.0),
          constrainedProcessType(FiniteElement::INVALID_PROCESS),
          constrainedPrimVar(FiniteElement::INVALID_PV),
          constrainedDirection(ConstrainedType::INVALID_CONSTRAINED_TYPE),
          constrainedVariable(
              ConstrainedVariable::INVALID_CONSTRAINED_VARIABLE),
          _isCompleteConstrained(false),
          _completeConstrainedStateOff(false),
          _isConstrainedVelStable(false),
          _isSeepageBC(false)
    {
    }
};

#endif
