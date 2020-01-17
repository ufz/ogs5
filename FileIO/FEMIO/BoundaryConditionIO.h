/*
 * BoundaryConditionIO.h
 *
 *  Created on: Apr 19, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef BOUNDARYCONDITIONIO_H_
#define BOUNDARYCONDITIONIO_H_

// STL
#include <iostream>
#include <string>

// GEOLIB
#include "GEOObjects.h"

// FEM
#include "rf_bc_new.h"

namespace FileIO
{
class BoundaryConditionIO
{
public:
    //	static CBoundaryCondition* read (std::istream& in,
    //			GEOLIB::GEOObjects const& geo_obj,
    //			std::string const& unique_fname);

    static void write(std::ostream& out, CBoundaryCondition const& bc);

    friend class CBoundaryCondition;
};
}  // namespace FileIO

#endif /* BOUNDARYCONDITIONIO_H_ */
