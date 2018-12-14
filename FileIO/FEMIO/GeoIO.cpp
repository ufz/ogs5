/*
 * GeoIO.cpp
 *
 *  Created on: Sep 29, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <sstream>

#include "FEMIO/GeoIO.h"

// FEM
#include "readNonBlankLineFromInputStream.h"

namespace FileIO
{
bool GeoIO::readGeoInfo(GeoInfo* geo_info, std::istream& in_str,
                        std::string& geo_name,
                        const GEOLIB::GEOObjects& geo_obj,
                        const std::string& unique_geo_name)
{
    std::stringstream strstream;
    strstream.str(readNonBlankLineFromInputStream(in_str));
    std::string geo_type_name;
    strstream >> geo_type_name;

    if (geo_type_name.find("POINT") != std::string::npos)
    {
        geo_info->setGeoType(GEOLIB::POINT);
        strstream >> geo_name;

        const GEOLIB::PointVec* pnt_vec(
            geo_obj.getPointVecObj(unique_geo_name));
        if (pnt_vec)
        {
            const GEOLIB::Point* pnt(pnt_vec->getElementByName(geo_name));
            if (pnt == NULL)
            {
                std::cerr << "ERROR in GeoIO::readGeoInfo: point name \""
                          << geo_name << "\" not found!"
                          << "\n";
                exit(1);
            }
            geo_info->setGeoObj(pnt);
            return true;
        }

        std::cerr << "Error in GeoIO::readGeoInfo: point vector not found!"
                  << "\n";
        exit(1);
    }

    else if (geo_type_name.find("POLYLINE") != std::string::npos)
    {
        geo_info->setGeoType(GEOLIB::POLYLINE);
        strstream >> geo_name;
        const GEOLIB::PolylineVec* ply_vec =
            geo_obj.getPolylineVecObj(unique_geo_name);
        if (ply_vec)
        {
            const GEOLIB::Polyline* ply(ply_vec->getElementByName(geo_name));
            if (ply == NULL)
            {
                std::cerr << "error in GeoIO::readGeoInfo: polyline name \""
                          << geo_name << "\" not found!"
                          << "\n";
                exit(1);
            }
            geo_info->setGeoObj(ply);
            return true;
        }

        std::cerr << "Error in GeoIO::readGeoInfo: polyline vector not found!"
                  << "\n";
        exit(1);
    }

    else if (geo_type_name.find("SURFACE") != std::string::npos)
    {
        geo_info->setGeoType(GEOLIB::SURFACE);
        strstream >> geo_name;
        GEOLIB::SurfaceVec const* sfc_vec(
            geo_obj.getSurfaceVecObj(unique_geo_name));
        if (sfc_vec)
        {
            const GEOLIB::Surface* sfc(sfc_vec->getElementByName(geo_name));
            if (sfc == NULL)
            {
                std::cerr << "Error in GeoIO::readGeoInfo: surface name \""
                          << geo_name << "\" not found!"
                          << "\n";
                exit(1);
            }
            geo_info->setGeoObj(sfc);
            return true;
        }

        std::cerr << "Error in GeoIO::readGeoInfo: surface vector not found!"
                  << "\n";
        exit(1);
    }

    else if (geo_type_name.find("VOLUME") != std::string::npos)
    {
        geo_info->setGeoType(GEOLIB::VOLUME);
        strstream >> geo_name;
        strstream.clear();
        return true;
    }

    else if (geo_type_name.find("DOMAIN") != std::string::npos)
    {
        geo_info->setGeoType(GEOLIB::GEODOMAIN);
        strstream >> geo_name;
        strstream.clear();
        return true;
    }

    return false;
}
}  // namespace FileIO
