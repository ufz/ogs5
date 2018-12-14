/*
 * ProcessIO.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// STL
#include <sstream>

// FileIO
#include "ProcessIO.h"
#include "readNonBlankLineFromInputStream.h"

namespace FileIO
{
bool ProcessIO::readProcessInfo(std::istream& in_str,
                                FiniteElement::ProcessType& pcs_type)
{
    std::stringstream ss_in(readNonBlankLineFromInputStream(in_str));
    std::string tmp;
    ss_in >> tmp;
    pcs_type = FiniteElement::convertProcessType(tmp);
    if (pcs_type == FiniteElement::INVALID_PROCESS)
        return false;
    else
        return true;
}
}  // namespace FileIO
