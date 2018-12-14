/*
 * extractMeshNodes.cpp
 *
 *  Created on: Dec 6, 2010
 *      Author: TF
 *
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// Base
#include "StringTools.h"

// FEM
#include "problem.h"

// MSH
#include "msh_lib.h"  // for FEMRead
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

int main(int argc, char* argv[])
{
    if (argc < 5)
    {
        std::cout << "program " << argv[0]
                  << " takes a mesh file and an index file and writes the mesh "
                     "node with the ids as gli points"
                  << std::endl;
        std::cout << "Usage: " << std::endl
                  << argv[0] << "\n\t--mesh ogs_meshfile\n\t--mesh-ids mesh-ids"
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

    // extract path
    std::string path;
    BaseLib::extractPath(argv[2], path);

    // *** read mesh ids
    tmp = argv[3];
    if (tmp.find("--mesh-ids") == std::string::npos)
    {
        std::cout << "could not extract mesh node id file name" << std::endl;
        return -1;
    }

    std::cout << "read mesh ids from file " << argv[4] << " ... " << std::flush;
    std::ifstream in(argv[4]);
    std::vector<size_t> mesh_ids;
    std::string id_as_str;
    size_t tmp_id;
    while (std::getline(in, id_as_str))
    {
        std::stringstream strstr(id_as_str);
        strstr >> tmp_id;
        mesh_ids.push_back(tmp_id);
    }
    std::cout << "done, " << mesh_ids.size() << " ids read" << std::endl;

    std::string out_fname(file_base_name + "-MeshNodesAsPoints.gli");
    std::ofstream os(out_fname.c_str());
    std::cout << "writing points to " << out_fname << " ... " << std::flush;
    FileIO::OGSMeshIO mesh_io;
    mesh_io.setMesh(mesh);
    os << "#POINTS" << std::endl;
    mesh_io.writeMeshNodesAsGLIPnts(mesh_ids, os);
    os << "#STOP" << std::endl;
    std::cout << "done" << std::endl;
    return 0;
}
