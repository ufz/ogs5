/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \author Wenqing Wang
 *  \file   PorosityDependentFractureRockPermeability.cpp
 *  Created on July 2, 2018, 10:23 AM
 */

#include "PorosityDependentFractureRockPermeability.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "msh_mesh.h"

namespace MaterialLib
{
PorosityDependentFractureRockPermeability::PorosityDependentFractureRockPermeability(const std::string& file_name,
                                                                                     const MeshLib::CFEMesh& mesh)
{
	std::ifstream inf(file_name.data());
	if (!inf.good())
	{
		std::cout << "File " << file_name << " for PorosityDependentFractureRockPermeability (type 9) "
		          << " is not found";
		exit(EXIT_FAILURE);
	}

	inf >> std::ws;
	while (!inf.eof())
	{
		double val;
		inf >> val >> std::ws;
		_omega.push_back(val);
	}

	if (_omega.size() != mesh.getElementVector().size())
	{
		std::cout << "The number of data in file " << file_name << " is not identical to the number of elements";
		exit(EXIT_FAILURE);
	}
}

double PorosityDependentFractureRockPermeability::getPermeability(const std::size_t element_id, const double porosity)
{
	const double omega = _omega[element_id];		
	return std::pow(2.7182818284590451,
	                (1. - omega) * std::log(4.979e-11 * std::pow(porosity, 3.11))
	                    + omega * std::log(1.143e-11 * std::pow(porosity, 0.64)));
}

} // End of name space
