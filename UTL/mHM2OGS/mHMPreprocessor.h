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
#include <vector>

#include "MSH/msh_mesh.h"

namespace FiniteElement
{
class ShapeFunctionPool;
}

namespace MeshLib
{
/**
 * Process the recharge data of mHM in order to used it as the Neumman BC on the
 * top surface of 3D domain.
 */
class mHMPreprocessor : public CFEMesh
{
public:
	mHMPreprocessor(std::string* geo_name) : CFEMesh(NULL, geo_name), _fem(NULL)
	{
	}

	~mHMPreprocessor();

	/*!
	   \brief Transform the precipitation data of mHM to the Neumann BC of the
	 groundwater flow equation.
	   06/2010  WW
	 */
	void transform_mHMData(const std::string& output_path);

private:
	FiniteElement::CElement* _fem;

	/*!  \brief Read GIS shapfile that stores the precipitation data.
	 * This function reads the data of mHM, which is stored in the syntax of
	 * raster file, finds face elements on the top surface, and then performs
	 * the numerical integration on the found surface elements to convert the
	 * mHM data into nodal flux values.
	 *  \param fname The input file name.
	 *  \param ofname The output file name.
	 *  \param ratio The ration of precipitation to the infiltration.
	 * 03/2010  WW
	  */
	void transfromSingle_mHMdataToNodalFlux(std::string const& fname,
	                                        std::string const& ofname,
	                                        double ratio = 0.8);
};
}  // end of namespace MeshLib
