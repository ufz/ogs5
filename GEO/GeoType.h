/*
 * GeoType.h
 *
 *  Created on: Jun 17, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef GEOTYPE_H_
#define GEOTYPE_H_

#include <string>

namespace GEOLIB
{
/**
 * \ingroup GEOLIB
 */

enum GEOTYPE
{
	INVALID = 0,
	POINT, //!< POINT
	POLYLINE, //!< POLYLINE
	SURFACE, //!< SURFACE
	VOLUME, //!< VOLUME
	GEODOMAIN //!< GEODOMAIN
};

GEOTYPE convertGeoType(const std::string& geo_type_str);

std::string convertGeoTypeToString(GEOTYPE geo_type);
} // end namespace GEOLIB

#endif /* GEOTYPE_H_ */
