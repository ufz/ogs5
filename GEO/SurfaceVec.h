/*
 * \file SurfaceVec.h
 *
 *  Created on: Feb 9, 2010
 *      Author: fischeth
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef SURFACEVEC_H_
#define SURFACEVEC_H_

#include "Surface.h"
#include "TemplateVec.h"

namespace GEOLIB
{
/**
 * Class SurfaceVec encapsulate a std::vector of Surfaces
 * and a name.
 * */

typedef TemplateVec<Surface> SurfaceVec;
} // end namespace

#endif /* SURFACEVEC_H_ */
