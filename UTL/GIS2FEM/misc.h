/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MISC_FEMTOOLKITS_H
#define MISC_FEMTOOLKITS_H

#include <string>

namespace MeshLib
{
class CFEMesh;
}

void writeOGSMesh(const MeshLib::CFEMesh& mesh, const std::string file_base_name);
#endif
