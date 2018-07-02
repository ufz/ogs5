/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \author Wenqing Wang
 *  \file   PorosityDependentFractureRockPermeability.h
 *  Created on July 2, 2018, 10:23 AM
 */

#pragma once

#include <vector>
#include <string>

namespace MeshLib
{
class CFEMesh;
}
namespace MaterialLib
{
/**
 *  A porosity dependent permeability model for fractured rock by
 *  \f[
 *        \mathrm{log} k = (1-\omage) \mathrm{log} k_i + \omage \mathrm{log} k_s
 *  \f]
 * where
 * \f[
 * k_i=4.979\times 10^{-11} n ^{3.11}
 * \f]
 * \f[
 * k_s=1.143\times 10^{-11} n ^{0.64}
 * \f]
 * and $f\omage \in [0,1]$f is a factor.
 */
class PorosityDependentFractureRockPermeability
{
public:
	PorosityDependentFractureRockPermeability(const std::string& file_name, const MeshLib::CFEMesh& mesh);

	double getPermeability(const std::size_t element_id, const double porosity);

private:
	std::vector<double> _omega; /// Element wise value.
};
} // End of name space
