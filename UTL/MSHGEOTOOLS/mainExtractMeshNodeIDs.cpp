/*
 * mainExtractMeshNodeIDs.cpp
 *
 *  Created on: Jun 13, 2012
 *      Author: fischeth
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// Base
#include "StringTools.h"

// FEM
#include "problem.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

// MSHGEOTOOLS
#include "ExtractMeshNodes.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PolylineVec.h"

// FileIO
#include "OGSIOVer4.h"

Problem* aproblem = NULL;

// MathLib
#include "Matrix.h"
#include "Vector3.h"

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		std::cout << "program " << argv[0] << " takes a *layered* mesh file and closed polylines from a geometry and computes the surface mesh nodes within the closed polyline" << std::endl;
		std::cout << "Usage: " << std::endl << argv[0] <<
		"\n\t--mesh ogs_meshfile\n\t--geometry ogs_geometry_as_gli_file" << std::endl;
		return -1;
	}

	// *** read mesh
	std::string tmp (argv[1]);
	if (tmp.find ("--mesh") == std::string::npos)
	{
		std::cout << "could not extract mesh file name" << std::endl;
		return -1;
	}

	tmp = argv[2];
	std::string file_base_name (tmp);
	if (tmp.find (".msh") != std::string::npos)
		file_base_name = tmp.substr (0, tmp.size() - 4);

	std::vector<MeshLib::CFEMesh*> mesh_vec;
	FEMRead(file_base_name, mesh_vec);
	if (mesh_vec.empty())
	{
		std::cerr << "could not read mesh from file " << std::endl;
		return -1;
	}
	MeshLib::CFEMesh* mesh (mesh_vec[mesh_vec.size() - 1]);
	mesh->ConstructGrid();

	// extract path
	std::string path;
	BaseLib::extractPath(argv[2], path);

	// *** read geometry
	tmp = argv[3];
	if (tmp.find ("--geometry") == std::string::npos)
	{
		std::cout << "could not extract geometry file name" << std::endl;
		return -1;
	}

	GEOLIB::GEOObjects* geo (new GEOLIB::GEOObjects);
	tmp = argv[4];
	std::string unique_name;
	std::vector<std::string> error_strings;
	FileIO::readGLIFileV4(tmp, geo, unique_name, error_strings);

//	{
//		const std::vector<GEOLIB::Point*>* pnts (geo->getPointVec (tmp));
//		if (pnts) {
//			std::string fname ("MeshIDs.txt");
//			std::ofstream out (fname.c_str());
//
//			std::string fname_gli ("MeshNodesAsPnts.gli");
//			std::ofstream pnt_out (fname_gli.c_str());
//			pnt_out << "#POINTS" << std::endl;
//
//			MeshLib::ExtractMeshNodes extract_mesh_nodes (mesh);
//
//			const size_t n_pnts (pnts->size());
//			for (size_t k(0); k<n_pnts; k++) {
//				extract_mesh_nodes.writeNearestMeshNodeToPoint (out, pnt_out, *((*pnts)[k]));
//			}
//			pnt_out << "#STOP" << std::endl;
//		}
//		return 0;
//	}

	// *** get Polygon
	const std::vector<GEOLIB::Polyline*>* plys (geo->getPolylineVec (unique_name));
	if (!plys)
	{
		std::cout << "could not get vector of polylines" << std::endl;
		delete mesh;
		delete geo;
		return -1;
	}

	std::vector<size_t> mesh_ids;
//	size_t ply_id (0);
//	size_t layer(1);
//	getMeshNodesFromLayerAlongPolyline(mesh, geo, unique_name, ply_id, layer, mesh_ids);
//	writeMeshNodes(mesh, mesh_ids, "MeshIDs.txt", "MeshNodesAsPoints.gli", true);

	//*** extract surface out of mesh
	MeshLib::ExtractMeshNodes extract_mesh_nodes (mesh);

//	// *** generate a polygon from polyline
////	std::vector<GEOLIB::Polyline*> polylines;
//	const size_t n_plys (plys->size());
//	for (size_t k(0); k < n_plys; k++)
//	{
//		bool closed ((*plys)[k]->isClosed());
//		if (!closed)
//		{
//			std::cout << "converting polyline " << k << " to closed polyline" << std::endl;
//			GEOLIB::Polygon* polygon(NULL);
//			extract_mesh_nodes.getPolygonFromPolyline(*((*plys)[k]), geo, unique_name, polygon);
////			polylines.push_back (polygon);
////			geo->appendPolylineVec (polylines, unique_name);
//			std::string *polygon_name(new std::string);
//			geo->getPolylineVecObj(unique_name)->getNameOfElementByID(k, *polygon_name);
//			(*polygon_name) += "-Polygon";
//			geo->getPolylineVecObj(unique_name)->push_back(polygon, polygon_name);
////			polylines.clear();
//		}
//	}
//
//	FileIO::writeGLIFileV4 ("New.gli", unique_name, *geo);

	// *** search mesh nodes for direct assigning bc, st or ic
	const size_t n_plys (plys->size());
	for (size_t k(0); k<n_plys; k++) {
		bool closed ((*plys)[k]->isClosed());
		if (!closed) {
			std::cout << "polyline " << k << " is not closed" << std::endl;
		} else {
			std::string fname (path + "MeshIDs.txt");
			std::ofstream out (fname.c_str());
			std::string fname_gli (path + "MeshNodesAsPnts.gli");
			std::ofstream pnt_out (fname_gli.c_str());
			pnt_out << "#POINTS" << std::endl;
			GEOLIB::Polygon polygon (*((*plys)[k]));
//			extract_mesh_nodes.writeMesh2DNodeIDAndArea (out, pnt_out, polygon);
			extract_mesh_nodes.writeTopSurfaceMeshNodeIDs (out, pnt_out, polygon);
			// write all nodes - not only the surface nodes
//			extract_mesh_nodes.writeMeshNodeIDs (out, pnt_out, polygon);
			pnt_out << "#STOP" << std::endl;
			out.close();
			pnt_out.close();
		}
	}

	// *** for Model Pipiripau
//	std::vector<GEOLIB::Polygon*> holes;
//	size_t bounding_polygon_id(0);
//	while (bounding_polygon_id < n_plys && ! (*plys)[bounding_polygon_id]->isClosed()) {
//		bounding_polygon_id++;
//	}
//
//	for (size_t k(bounding_polygon_id+1); k<n_plys; k++) {
//		bool closed ((*plys)[k]->isClosed());
//		if (!closed) {
//			std::cout << "polyline " << k << " is not closed" << std::endl;
//		} else {
//			holes.push_back (new GEOLIB::Polygon(*(((*plys)[k]))));
//		}
//	}
//	extract_mesh_nodes.writeMesh2DNodeIDAndArea (out, pnt_out, GEOLIB::Polygon((*((*plys)[bounding_polygon_id]))), holes);
//	for (size_t k(0); k<holes.size(); k++) {
//		delete holes[k];
//	}

//	out << "#STOP" << std::endl;
//	out.close();
//
//	pnt_out << "#STOP" << std::endl;
//	pnt_out.close ();

	delete mesh;
	delete geo;

	return 0;
}
