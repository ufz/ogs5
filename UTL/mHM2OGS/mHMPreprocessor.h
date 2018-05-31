/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   mHMPreprocessor.h
 *  Created on May 31, 2018, 1:39 PM by WW
 */

#pragma once

#include <string>

#include "MSH/msh_mesh.h"

namespace MeshLib
{
/**
 * Process the recharge data of mHM in order to used it as the Neumman BC on the
 * top surface of 3D domain.
 */
class mHMPreprocessor : public CFEMesh
{
public:
	mHMPreprocessor(std::string* geo_name) : CFEMesh(NULL, geo_name) {}
	void MarkInterface_mHM_Hydro_3D();
	void mHM2NeumannBC(const std::string output_path);
	/// Compute int {f} a dA on top surface.
	void TopSurfaceIntegration();

private:
};
} // end of namespace MeshLib
