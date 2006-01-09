/*
 * MeshNodesAlongPolyline.cpp
 *
 *  Created on: Aug 9, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// Base
#include "quicksort.h"

// MSH
#include "MeshNodesAlongPolyline.h"
#include "msh_mesh.h"

#include <algorithm>

namespace MeshLib
{
MeshNodesAlongPolyline::MeshNodesAlongPolyline(
		GEOLIB::Polyline const* const ply, MeshLib::CFEMesh const* mesh,
		double search_radius) :
	_ply(ply), _mesh(mesh), _linear_nodes (0)
{
	std::vector<CNode*> const& mesh_nodes (mesh->getNodeVector());

	size_t n_linear_order_nodes (mesh->GetNodesNumber (false));
	size_t n_nodes (mesh->GetNodesNumber (true));

	std::vector<size_t> msh_node_higher_order_ids;
	std::vector<double> dist_of_proj_higher_order_node_from_ply_start;

	//We need exactly defined polyline for DDC. If there is not any node located within the
	// threhold the ployline, we do not forced to fill the _msh_node_ids.
	//Therefore, the repeating loop is skipped.
	// repeat until at least one relevant node was found
	// WW while (_msh_node_ids.empty())
	{
		// loop over all line segments of the polyline
		for (size_t k = 0; k < ply->getNumberOfPoints() - 1; k++)
		{
			double act_length_of_ply(ply->getLength(k));
			double seg_length (sqrt(MathLib::sqrDist(_ply->getPoint(k), _ply->getPoint(k + 1))));
			double lower_lambda (- search_radius / seg_length);
			double upper_lambda (1 + search_radius / seg_length);

			// loop over all nodes
			for (size_t j = 0; j < n_nodes; j++)
			{
				double dist, lambda;

				// is the orthogonal projection of the j-th node to the
				// line g(lambda) = _ply->getPoint(k) + lambda * (_ply->getPoint(k+1) - _ply->getPoint(k))
				// at the k-th line segment of the polyline, i.e. 0 <= lambda <= 1?
				if (MathLib::calcProjPntToLineAndDists(mesh_nodes[j]->getData(),
								(_ply->getPoint(k))->getData(), (_ply->getPoint(k + 1))->getData(),
								lambda, dist) <= search_radius) {
					if (lower_lambda <= lambda && lambda <= upper_lambda) {
						if (mesh_nodes[j]->GetIndex() < n_linear_order_nodes) {
							// check if node id is already in the vector
							if (std::find(_msh_node_ids.begin(), _msh_node_ids.end(),
											mesh_nodes[j]->GetIndex()) == _msh_node_ids.end()) {
								_msh_node_ids.push_back(mesh_nodes[j] ->GetIndex());
								_dist_of_proj_node_from_ply_start. push_back(act_length_of_ply + dist);
								_linear_nodes++;
							}
						} else {
							// check if node id is already in the vector
							if (std::find(msh_node_higher_order_ids.begin(),
											msh_node_higher_order_ids.end(), mesh_nodes[j]->GetIndex())
											== msh_node_higher_order_ids.end()) {
								msh_node_higher_order_ids.push_back(mesh_nodes[j]->GetIndex());
								dist_of_proj_higher_order_node_from_ply_start .push_back(act_length_of_ply
												+ dist);
							}
						}
					} // end if lambda
				}
			} // end node loop
		} // end line segment loop

		//We need exactly defined polyline for DDC.
		//Therefore, the following two line should be dropped.
		//if (_msh_node_ids.empty())
		//	search_radius *= 2.0;
	}

	// sort the (linear) nodes along the polyline according to their distances
	Quicksort<double> (_dist_of_proj_node_from_ply_start, 0,
					_dist_of_proj_node_from_ply_start.size(), _msh_node_ids);

#ifndef NDEBUG
//	std::cout << "[DEBUG-INFO] " << "\n";
//
//	// fetch polyline name
//	std::string ply_name;
//	if (! _mesh->getGEOObjects()->getPolylineVecObj(*(_mesh->getProjectName()))->getNameOfElement(ply, ply_name)) {
//		ply_name = "unknown-ply";
//	}
//
//	std::cout << "distances of linear nodes along polyline " << ply_name <<
//	" (epsilon radius = " << search_radius << "): " << "\n";
//	for (size_t k(0); k < _dist_of_proj_node_from_ply_start.size(); k++)
//		std::cout << "\t" << _msh_node_ids[k] << " " <<
//		_dist_of_proj_node_from_ply_start[k] << "\n";
//	std::cout << "number of linear nodes along polyline " << ply_name << ": " <<
//	_dist_of_proj_node_from_ply_start.size()
//	          << ", number of higher order nodes: " << msh_node_higher_order_ids.size() <<
//	"\n";
#endif
	// assign/append higher order nodes at the end of vector _msh_node_ids
	for (size_t k(0); k < msh_node_higher_order_ids.size(); k++)
		_msh_node_ids.push_back (msh_node_higher_order_ids[k]);
	// append distances for higher order nodes at the end of vector _dist_of_proj_node_from_ply_start
	for (size_t k(0); k < dist_of_proj_higher_order_node_from_ply_start.size(); k++)
		_dist_of_proj_node_from_ply_start.push_back(dist_of_proj_higher_order_node_from_ply_start[k]);
}

const std::vector<size_t>& MeshNodesAlongPolyline::getNodeIDs () const
{
	return _msh_node_ids;
}

const GEOLIB::Polyline* MeshNodesAlongPolyline::getPolyline () const
{
	return _ply;
}

size_t MeshNodesAlongPolyline::getNumberOfLinearNodes () const
{
	return _linear_nodes;
}

const std::vector<double>& MeshNodesAlongPolyline::getDistOfProjNodeFromPlyStart() const
{
	return _dist_of_proj_node_from_ply_start;
}
} // end namespace
