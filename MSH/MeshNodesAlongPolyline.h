/*
 * MeshNodesAlongPolyline.h
 *
 *  Created on: Aug 9, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef MESHNODESALONGPOLYLINE_H_
#define MESHNODESALONGPOLYLINE_H_

// GEOLIB
#include "Polyline.h"

namespace MeshLib
{
// forward declaration
class CFEMesh;

/**
 * This class computes the ids of the mesh nodes along a polyline.
 *
 * The mesh nodes are sorted as follow:
 * [ ... ids of sorted linear nodes ... | ... ids of unsorted higher order nodes ]
 */
class MeshNodesAlongPolyline
{
public:
	MeshNodesAlongPolyline(GEOLIB::Polyline const* const ply,
		CFEMesh const* mesh, double search_radius);
	const std::vector<size_t>& getNodeIDs () const;
	const GEOLIB::Polyline* getPolyline () const;
	size_t getNumberOfLinearNodes () const;
	std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
	const GEOLIB::Polyline* _ply;
	const CFEMesh* _mesh;
	size_t _linear_nodes;
	std::vector<size_t> _msh_node_ids;
	std::vector<double> _dist_of_proj_node_from_ply_start;
};
}

#endif /* MESHNODESALONGPOLYLINE_H_ */
