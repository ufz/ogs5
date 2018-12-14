/*
 * OGSMeshIO.cpp
 *
 *  Created on: Mar 3, 2011
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GEOObjects.h"
#include "MeshIO/OGSMeshIO.h"
#include "msh_lib.h"
#include "msh_mesh.h"
#include "msh_elem.h"

#include <iomanip>
#include <ctime>

namespace FileIO
{
OGSMeshIO::OGSMeshIO() : _mesh(NULL) {}

MeshLib::CFEMesh* OGSMeshIO::loadMeshFromFile(std::string const& fname)
{
    std::cout << "FEMRead ... " << std::flush;

#ifndef NDEBUG
    clock_t start(clock());
#endif

    std::vector<MeshLib::CFEMesh*> mesh_vec;
    FEMRead(fname.substr(0, fname.length() - 4), mesh_vec);

    if (!mesh_vec.empty())
    {
#ifndef NDEBUG
        clock_t start_construct(clock());
#endif

        mesh_vec[mesh_vec.size() - 1]->ConstructGrid();

        std::cout << "Nr. Nodes: "
                  << mesh_vec[mesh_vec.size() - 1]->nod_vector.size() << "\n";
        std::cout << "Nr. Elements: "
                  << mesh_vec[mesh_vec.size() - 1]->ele_vector.size() << "\n";

#ifndef NDEBUG
        clock_t end_construct(clock());
        std::cout << "constructGrid time: "
                  << (end_construct - start_construct) /
                         (double)(CLOCKS_PER_SEC)
                  << " s"
                  << "\n";
#endif

        mesh_vec[mesh_vec.size() - 1]->FillTransformMatrix();

#ifndef NDEBUG
        clock_t end(clock());
        std::cout << "total loading time: "
                  << (end - start) / (double)(CLOCKS_PER_SEC) << " s"
                  << "\n";
#endif

        return mesh_vec[mesh_vec.size() - 1];
    }

    std::cout << "Failed to load the mesh file: " << fname << "\n";
    return NULL;
}

int OGSMeshIO::write(std::ostream& out)
{
    if (!_mesh)
    {
        std::cout << "OGSMeshIO cannot write: no mesh set!"
                  << "\n";
        return 0;
    }

    setPrecision(9);

    out << "#FEM_MSH"
        << "\n";

    out << "$PCS_TYPE"
        << "\n"
        << "  " << _mesh->pcs_name << "\n";

    out << "$NODES"
        << "\n"
        << "  ";
    const size_t n_nodes(_mesh->GetNodesNumber(false));
    out << n_nodes << "\n";
    for (size_t i(0); i < n_nodes; i++)
    {
        double const* const coords(_mesh->nod_vector[i]->getData());
        out << i << " " << coords[0] << " " << coords[1] << " " << coords[2]
            << "\n";
    }

    out << "$ELEMENTS"
        << "\n"
        << "  ";
    writeElementsExceptLines(_mesh->ele_vector, out);

    out << " $LAYER"
        << "\n";
    out << "  ";
    out << _mesh->_n_msh_layer << "\n";
    out << "#STOP"
        << "\n";
    return 1;
}

void OGSMeshIO::writeElementsExceptLines(
    std::vector<MeshLib::CElem*> const& ele_vec, std::ostream& out)
{
    const size_t ele_vector_size(ele_vec.size());
    const double epsilon(std::numeric_limits<double>::epsilon());
    std::vector<bool> non_line_element(ele_vector_size, true);
    std::vector<bool> non_null_element(ele_vector_size, true);
    size_t n_elements(0);

    for (size_t i(0); i < ele_vector_size; i++)
    {
        if ((ele_vec[i])->GetElementType() == MshElemType::LINE)
        {
            non_line_element[i] = false;
            non_null_element[i] = false;
        }
        else
        {
            if (ele_vec[i]->calcVolume() < epsilon)
            {
                non_null_element[i] = false;
            }
            else
            {
                n_elements++;
            }
        }
    }
    out << n_elements << "\n";
    for (size_t i(0), k(0); i < ele_vector_size; i++)
    {
        if (non_line_element[i] && non_null_element[i])
        {
            const size_t tmp_idx(ele_vec[i]->GetIndex());
            (const_cast<MeshLib::CElem*>(ele_vec[i]))->SetIndex(k);
            ele_vec[i]->WriteIndex(out);
            (const_cast<MeshLib::CElem*>(ele_vec[i]))->SetIndex(tmp_idx);
            k++;
        }
    }
}

void OGSMeshIO::setMesh(MeshLib::CFEMesh const* mesh)
{
    _mesh = mesh;
}

void OGSMeshIO::writeMeshNodesAsGLIPnts(
    std::vector<size_t> const& mesh_node_ids, std::ostream& os)
{
    std::vector<MeshLib::CNode*> const& nodes(_mesh->getNodeVector());
    for (size_t k(0); k < mesh_node_ids.size(); k++)
    {
        MeshLib::CNode const& node(*(nodes[mesh_node_ids[k]]));
        double const* const coords(node.getData());
        _out << k << " " << coords[0] << " " << coords[1] << " " << coords[2]
             << " $NAME " << node.GetIndex() << "\n";
    }

    os << _out.str();
}

}  // end namespace FileIO
