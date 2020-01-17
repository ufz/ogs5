/*
 * DistributionInfo.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: fischeth
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "DistributionInfo.h"

DistributionInfo::DistributionInfo(FiniteElement::DistributionType dt)
    : _dis_type(dt)
{
}

DistributionInfo::~DistributionInfo() {}

void DistributionInfo::setProcessDistributionType(
    FiniteElement::DistributionType dis_type)

{
    _dis_type = dis_type;
}

FiniteElement::DistributionType DistributionInfo::getProcessDistributionType()
    const
{
    return _dis_type;
}
