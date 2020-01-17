/*
 * PolylineVec.h
 *
 *  Created on: Feb 9, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef POLYLINEVEC_H_
#define POLYLINEVEC_H_

#include "Polyline.h"
#include "TemplateVec.h"

namespace GEOLIB
{
/**
 * \ingroup GEOLIB
 *
 * \brief class PolylineVec encapsulate a std::vector of Polylines
 * additional one can give the vector of polylines a name
 * */

typedef TemplateVec<Polyline> PolylineVec;
}  // namespace GEOLIB

#endif /* POLYLINEVEC_H_ */
