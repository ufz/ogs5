/*
 * modifyMeshProperties.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 *
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ModifyMeshProperties.h"

// FEM
#include "problem.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PolylineVec.h"

// FileIO
#include "MeshIO/OGSMeshIO.h"
#include "OGSIOVer4.h"

Problem* aproblem = NULL;

char program_description[]
    = "This program reads an OpenGeoSys (OGS) mesh, a geometry (in OGS gli format) and\n\
an id of a polygon described inside of the geometry. With the help of the\n\
referenced polygon a volume is described. First the polygon is projected to\n\
the bottom and top surface of the mesh and we obtain parts of the bottom and\n\
top surfaces of the mesh. The bottom and top surfaces are linked together\n\
via a third surface. The volume is defined as the part of the mesh that is\n\
inside the three surfaces. Within this volume the material ids of the elements\n\
are modified. Programmed by Tom (thomas.fischer@ufz.de)";

int main(int argc, char* argv[])
{
	if (!(argc == 9 || argc == 11))
	{
		std::cout << "********************************************************************************" << std::endl;
		std::cout << program_description << std::endl;
		std::cout << "********************************************************************************" << std::endl;
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0]
		          << "\n\t--mesh ogs_meshfile\n\t--geometry ogs_geometry\n\t--polygon_id id\n\t--material_id id"
		          << std::endl;
		std::cout << "*** or ***" << std::endl;
		std::cout << argv[0] << "\n\t--mesh ogs_meshfile\n\t--geometry ogs_geometry\n\t--polygon_id "
		                        "id\n\t--old_material_id old_id\n\t--new_material_id new_id"
		          << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp(argv[1]);
	if (tmp.find("--mesh") == std::string::npos)
	{
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name(tmp);
	std::string path;
	BaseLib::extractPath(file_base_name, path);
	if (tmp.find(".msh") != std::string::npos)
		file_base_name = tmp.substr(0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh(mesh_vec[mesh_vec.size() - 1]);

	mesh->ConstructGrid();

	// *** read geometry
	tmp = argv[3];
	if (tmp.find("--geometry") == std::string::npos)
	{
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo(new GEOLIB::GEOObjects);
	tmp = argv[4];
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(tmp, geo, unique_name, error_strings);

	// *** get polygons
	tmp = BaseLib::getFileNameFromPath(tmp, true);
	const std::vector<GEOLIB::Polyline*>* plys(geo->getPolylineVec(tmp));
	if (!plys)
	{
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	// *** get polygon id
	tmp = argv[5];
	if (tmp.find("--polygon_id") == std::string::npos)
	{
		std::cout << "could not read polygon id" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}
	size_t polygon_id(atoi(argv[6]));

	if (plys->size() <= polygon_id)
	{
		std::cout << "polyline for id " << polygon_id << " not found" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	bool closed((*plys)[polygon_id]->isClosed());
	if (!closed)
	{
		std::cout << "polyline with id " << polygon_id << " is not closed" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	size_t old_material_id(0), new_material_id(0);
	if (argc == 9)
	{
		// *** get material id
		tmp = argv[7];
		if (tmp.find("--material") == std::string::npos)
		{
			std::cout << "could not read material id" << std::endl;
			delete mesh;
			delete geo;
			return -1;
		}
		new_material_id = atoi(argv[8]);
	}
	if (argc == 11)
	{
		// *** get old material id
		tmp = argv[7];
		if (tmp.find("--old_material_id") == std::string::npos)
		{
			std::cout << "could not read old material id" << std::endl;
			delete mesh;
			delete geo;
			return -1;
		}
		old_material_id = atoi(argv[8]);
		tmp = argv[9];
		if (tmp.find("--new_material_id") == std::string::npos)
		{
			std::cout << "could not read new material id" << std::endl;
			delete mesh;
			delete geo;
			return -1;
		}
		new_material_id = atoi(argv[10]);
	}

	GEOLIB::Polygon polygon(*((*plys)[polygon_id]));
	MeshLib::ModifyMeshProperties modify_mesh_nodes(mesh);

	if (argc == 9)
	{
		modify_mesh_nodes.setMaterial(polygon, new_material_id);
	}
	else
	{
		modify_mesh_nodes.substituteMaterialID(polygon, old_material_id, new_material_id);
	}

	std::string mesh_out_fname(path + "MeshWithMaterial.msh");
	std::cout << "write " << mesh_out_fname << " ... " << std::flush;
	FileIO::OGSMeshIO mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile(mesh_out_fname);
	std::cout << "ok" << std::endl;

	return 0;
}
