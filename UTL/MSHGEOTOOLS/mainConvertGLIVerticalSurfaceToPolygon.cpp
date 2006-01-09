/*
 * mainConvertGLIVerticalSurfaceToPolygon.cpp
 *
 *  Created on: Jun 5, 2012
 *      Author: fischeth
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// stl
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

// FileIO
#include "OGSIOVer4.h"
#include "OGSMeshIO.h"

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "msh_mesh.h"

// MeshLib - folder MSHGEOTOOLS
#include "ExtractMeshNodes.h"

int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cout << "Program " << argv[0]
						<< " converts a vertical surface described by a polyline embedded in a mesh to a polygon"
						<< std::endl;
		std::cout << "For this reason the polyline is projected to the bottom surface and top surface "
						"of the mesh. The resulting two polylines are linked to a closed polyline, i.e., a polygon." << std::endl;
		std::cout << "Usage: " << argv[0]
						<< " --mesh ogs_meshfile --geometry ogs_geometry_as_gli_file [--polyline-id-lower ply_id_lower --polyline-id-upper ply_id_upper]" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp(argv[1]);
	if (tmp.find("--mesh") == std::string::npos) {
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name(tmp);
	if (tmp.find(".msh") != std::string::npos)
		file_base_name = tmp.substr(0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty()) {
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh * mesh(mesh_vec[mesh_vec.size() - 1]);
	mesh->ConstructGrid();

	// *** read geometry
	tmp = argv[3];
	if (tmp.find("--geometry") == std::string::npos) {
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo(new GEOLIB::GEOObjects);
	tmp = argv[4];
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(tmp, geo, unique_name, error_strings);

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys(geo->getPolylineVec(unique_name));
	std::cout << "fetched polylines: " << std::flush << plys->size() << std::endl;
	if (!plys) {
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	//*** extract surface out of mesh
	MeshLib::ExtractMeshNodes extract_mesh_nodes(mesh);

	// *** generate a polygon from polyline
	size_t ply_id_upper(plys->size());
	size_t ply_id_lower(0);
	if (argc == 9) {
		ply_id_lower = str2number<size_t>(argv[6]);
		ply_id_upper = std::min(str2number<size_t>(argv[8]), ply_id_upper);
	}

	for (size_t k(ply_id_lower); k < ply_id_upper; k++) {
		bool closed((*plys)[k]->isClosed());
		if (!closed) {
			std::cout << "converting polyline " << k << " to polygon (closed polyline) "
							<< std::endl;
			GEOLIB::Polygon* polygon(NULL);
			extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo, unique_name, polygon);
			std::string *polygon_name(new std::string);
			geo->getPolylineVecObj(unique_name)->getNameOfElementByID(k, *polygon_name);
			(*polygon_name) += "-Polygon";
			geo->getPolylineVecObj(unique_name)->push_back(polygon, polygon_name);
		}
	}

	std::string path;
	BaseLib::extractPath(argv[4], path);
	std::string const fname_of_new_file(path + "New.gli");
	FileIO::writeGLIFileV4(fname_of_new_file, unique_name, *geo);

	delete mesh;
	delete geo;

	return 0;
}
