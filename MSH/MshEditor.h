/**
 * \file MshEditor.h
 * 2011/06/15 KR Initial implementation
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef MSHEDITOR_H
#define MSHEDITOR_H

#include <cstddef>
#include <vector>

namespace GEOLIB {
	class PointWithID;
}

namespace MeshLib
{
class CFEMesh;
}

class GridAdapter;

/**
 * \brief A set of tools for manipulating existing meshes
 */
class MshEditor
{
public:
	MshEditor() {}
	~MshEditor() {}

	/// Returns the area assigned to each node on a surface mesh.
	static void getNodeAreas(const MeshLib::CFEMesh* mesh, std::vector<double> &node_area_vec);

	/// Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
	static MeshLib::CFEMesh* removeMeshNodes(MeshLib::CFEMesh* mesh,
	                                         const std::vector<size_t> &nodes);

	/// Returns the surface nodes of a layered mesh.
	static std::vector<GEOLIB::PointWithID*> getSurfaceNodes(const MeshLib::CFEMesh &mesh);

	/// MW: populate sort nodes vector
	static void sortNodesLexicographically(MeshLib::CFEMesh *mesh);

	/// Returns the 2d-element mesh representing the surface of the given layered mesh.
	static MeshLib::CFEMesh* getMeshSurface(const MeshLib::CFEMesh &mesh);

private:
};

#endif //MSHEDITOR_H
