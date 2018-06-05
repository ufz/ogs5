/*
 * GeoObject.h
 *
 *  Created on: Aug 27, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef GEOOBJECT_H_
#define GEOOBJECT_H_

namespace GEOLIB
{
/**
 * \ingroup GEOLIB
 *
 * \brief Base class for classes Point, Polyline, Surface.
 */

class GeoObject
{
public:
	GeoObject() {}
	virtual ~GeoObject() {}
};
} // end namespace GEOLIB

#endif /* GEOOBJECT_H_ */
