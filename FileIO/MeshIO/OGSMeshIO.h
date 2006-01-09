/*
 * OGSMeshIO.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef OGSMESHIO_H_
#define OGSMESHIO_H_

#include "Writer.h"

#include <sstream>
#include <iostream>
#include <vector>

namespace MeshLib
{
class CFEMesh;
class CElem;
}

namespace FileIO
{
class OGSMeshIO : public Writer
{
public:
	/// @brief Constructor.
	OGSMeshIO();

	/// @brief Read a OGS mesh from file.
	MeshLib::CFEMesh* loadMeshFromFile(std::string const& fileName);

	/// @brief Sets the mesh.
	void setMesh(MeshLib::CFEMesh const* mesh);

	void writeMeshNodesAsGLIPnts (std::vector<size_t> const& mesh_node_ids, std::ostream & os);

protected:
	/// @brief Write functionality.
	int write(std::ostream &out);

private:
	void writeElementsExceptLines (std::vector<MeshLib::CElem*> const& ele_vec, std::ostream &out);

	MeshLib::CFEMesh const* _mesh;
};
} // end namespace FileIO

#endif /* OGSMESHIO_H_ */
