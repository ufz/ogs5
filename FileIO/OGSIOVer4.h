/*
 * OGSIOVer4.h
 *
 *  Created on: Jan 14, 2010
 *      Author: TF / KR
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef OGSIOVER4_H_
#define OGSIOVER4_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Point.h"
#include "Polyline.h"
#include "StringTools.h"
#include "Surface.h"

// forward declaration
namespace GEOLIB
{
class GEOObjects;
}

namespace FileIO
{
/** I/O - routines for the OGS-4 gli file format */

/** method reads geometric objects from file in gli format */
bool readGLIFileV4(const std::string& fname,
                   GEOLIB::GEOObjects* geo,
                   std::string& unique_name,
                   std::vector<std::string>& errors);

void writeGLIFileV4(const std::string& fname, const std::string& proj_name, const GEOLIB::GEOObjects& geo);

void writeAllDataToGLIFileV4(const std::string& fname, const GEOLIB::GEOObjects& geo);

} // end namespace

#endif /* OGSIOVER4_H_ */
