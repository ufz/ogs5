/**
 * \file StationIO.h
 * 23/03/2010 KR Initial implementation
 *
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef STATIONIO_H
#define STATIONIO_H

#include "Station.h"

/**
 * \brief A number of methods for data input and output for station data.
 *
 * A number of methods for data input and output for station data.
 */
class StationIO
{
public:
    /// Imports a file with station data.
    static int readStationFile(const std::string& path,
                               std::string& name,
                               std::vector<GEOLIB::Point*>* stations,
                               GEOLIB::Station::StationType type);

    /// Writes a file that contains all stratigraphies for the boreholes in the
    /// given vector.
    static void writeStratigraphyTable(
        const std::vector<GEOLIB::Point*>* boreholes,
        const std::string& filename);

private:
};

#endif  // STATIONIO_H
