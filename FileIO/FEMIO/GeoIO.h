/*
 * GeoIO.h
 *
 *  Created on: Sep 29, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef GEOIO_H_
#define GEOIO_H_

// STL
#include <fstream>
#include <string>

// FEM
#include "GeoInfo.h"

// GEO
#include "GEOObjects.h"

namespace FileIO
{
/**
 * Small class to read/write geometric information.
 * At the moment only the read-part is implemented.
 * ToDo: implement write part.
 */
class GeoIO
{
public:
    /**
     * @brief read geometric information from different files.
     *
     * Method reads geometric information from
     * - source term files (*.st, CSourceTerm::Read()),
     * - boundary condition files (*.bc, CBoundaryCondition::Read()),
     * - initial condition files (*.ic, CInitialCondition::Read()),
     * - output files (*.out, COutput::Read())
     * .
     * To store the information an object of type GeoInfo is used (CSourceTerm,
     * CBoundaryCondition, CInitialCondition and COutput inherit from GeoInfo).
     * @param geo_info object, in which the information is stored
     * @param in_str the input stream, containing the geometric information
     * @param geo_name the name of the geometric object (needed writing)
     * @param geo_obj an instance of class GEOObjects (manager of geometric
     * information)
     * @param unique_geo_name identifier of the geometric data within the
     * geometric manager GEOObjects
     * @return false if an error occured, else true
     */
    static bool readGeoInfo(GeoInfo* geo_info,
                            std::istream& in_str,
                            std::string& geo_name,
                            const GEOLIB::GEOObjects& geo_obj,
                            const std::string& unique_geo_name);
};
}  // namespace FileIO

#endif /* GEOIO_H_ */
