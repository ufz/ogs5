/*
 * filterMeshNodes.cpp
 *
 *  Created on: Jan 3, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// BaseLib
#include "StringTools.h"

// MSH
#include "msh_lib.h" // for FEMRead
#include "msh_mesh.h"

Problem* aproblem = NULL;

// Base
#include "binarySearch.h"
#include "quicksort.h"

int main(int argc, char* argv[])
{
	if (argc < 7)
	{
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << "\n\t--nodes-with-area nodes_with_area_file" << std::endl;
		std::cout << "\t--nodes-in-polygon nodes_in_polygon" << std::endl;
		std::cout << "\t--output output_file_name" << std::endl;
		std::cout << "\t[--mesh mesh]" << std::endl;
		return -1;
	}

	std::string tmp(argv[1]);
	if (tmp.find("--nodes-with-area") == std::string::npos)
	{
		std::cout << "could not find switch nodes_with_area" << std::endl;
		return -1;
	}

	tmp = argv[3];
	if (tmp.find("--nodes-in-polygon") == std::string::npos)
	{
		std::cout << "could not find switch nodes-in-polygon" << std::endl;
		return -1;
	}

	tmp = argv[5];
	if (tmp.find("--output") == std::string::npos)
	{
		std::cout << "missing switch output" << std::endl;
		return -1;
	}

	std::string input_base_fname(BaseLib::getFileNameFromPath(argv[2], false));
	std::string path;
	BaseLib::extractPath(argv[2], path);
	tmp = argv[2];
	std::ifstream in0(tmp.c_str());
	tmp = argv[4];
	std::ifstream in1(tmp.c_str());
	if (!in0 || !in1)
		std::cout << "could not open one of the input files" << std::endl;

	std::vector<size_t> ids0, ids1;
	std::vector<double> areas;

	size_t tmp_id;
	double tmp_area_val;

	/*** load first input file */
	bool all_ok(true);
	while (!in0.eof() && all_ok)
	{
		in0 >> tmp_id;
		if (in0.fail())
			all_ok = false;
		if (all_ok)
		{
			in0 >> tmp_area_val;
			if (in0.fail())
				all_ok = false;
			else
			{
				ids0.push_back(tmp_id);
				areas.push_back(tmp_area_val);
			}
		}
	}
	in0.close();
	std::cout << "read " << ids0.size() << " ids and " << areas.size() << " areas from file " << argv[2] << std::endl;

	all_ok = true;
	/*** load second input file */
	while (!in1.eof() && all_ok)
	{
		in1 >> tmp_id;
		if (!in1.fail())
			ids1.push_back(tmp_id);
		else
			all_ok = false;
	}
	in1.close();
	std::cout << "read " << ids1.size() << " ids from file " << argv[4] << std::endl;

	std::vector<size_t> perm;
	for (size_t k(0); k < ids0.size(); k++)
		perm.push_back(k);
	Quicksort<size_t>(ids0, 0, ids0.size(), perm);

	// apply permutation to areas
	std::vector<double> sorted_areas(areas.size());
	for (size_t k(0); k < perm.size(); k++)
		sorted_areas[k] = areas[perm[k]];

	// vector of flags to mark ids within vector ids0
	std::vector<bool> not_found_ids_within_ids0(ids0.size(), true); // per default

	std::vector<size_t> not_found_ids;
	std::ofstream out(argv[6]);
	if (!out)
	{
		std::cout << "could not open file " << argv[6] << " for output" << std::endl;
		return -1;
	}
	out.precision(12);
	for (size_t k(0); k < ids1.size(); k++)
	{
		size_t idx(searchElement(ids1[k], 0, ids0.size(), ids0));
		if (idx != std::numeric_limits<size_t>::max())
		{
			out << ids0[idx] << " " << sorted_areas[idx] << std::endl;
			not_found_ids_within_ids0[idx] = false;
		}
		else
		{
			not_found_ids.push_back(ids1[k]);
		}
	}
	out.close();

	std::string fname_remaining(path + input_base_fname + "-Remaining.txt");
	std::ofstream os_remaining(fname_remaining.c_str());
	std::cout << "writing to " << fname_remaining << " ... " << std::flush;
	for (size_t k(0); k < ids0.size(); k++)
	{
		if (not_found_ids_within_ids0[k])
		{
			os_remaining << ids0[k] << " " << sorted_areas[k] << std::endl;
		}
	}
	os_remaining.close();
	std::cout << "done" << std::endl;

	if (!not_found_ids.empty())
	{
		std::sort(not_found_ids.begin(), not_found_ids.end());
		for (size_t k(0); k < not_found_ids.size(); k++)
			std::cout << not_found_ids[k] << std::endl;
		std::cout << "number of not found ids: " << not_found_ids.size() << std::endl;
	}

	if (argc > 8)
	{
		tmp = argv[7];
		if (tmp.find("--mesh") != std::string::npos)
		{
			tmp = argv[8];
			if (tmp.find(".msh") != std::string::npos)
				tmp = tmp.substr(0, tmp.size() - 4);

			std::cout << "reading mesh " << tmp << " ... " << std::flush;

			std::vector<MeshLib::CFEMesh*> mesh_vec;
			FEMRead(tmp, mesh_vec);
			if (mesh_vec.empty())
			{
				std::cerr << "could not read mesh from file " << std::endl;
				return -1;
			}
			std::cout << "done" << std::endl;
			MeshLib::CFEMesh* mesh(mesh_vec[mesh_vec.size() - 1]);

			std::ofstream out_mesh_nodes_as_pnts("MeshNodesAsPnts.gli");
			if (out_mesh_nodes_as_pnts)
			{
				out_mesh_nodes_as_pnts << "#POINTS" << std::endl;
				const size_t ids0_size(ids0.size());
				for (size_t k(0); k < ids0_size; k++)
				{
					double const* const pnt(mesh->getNodeVector()[ids0[k]]->getData());
					out_mesh_nodes_as_pnts << k << " " << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
				}
				out_mesh_nodes_as_pnts << "#STOP" << std::endl;
				out_mesh_nodes_as_pnts.close();
			}
		}
	}
}
