/**
 * \file MshEditor.cpp
 * 2011/06/15 KR Initial implementation
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MshEditor.h"
#include "PointWithID.h"
#include "msh_mesh.h"
#include "GridAdapter.h"

void MshEditor::getNodeAreas(const MeshLib::CFEMesh* mesh,
                             std::vector<double>& node_area_vec)
{
    double total_area(0);

    // for each node, a vector containing all the element idget every element
    const size_t nNodes(mesh->nod_vector.size());
    for (size_t n = 0; n < nNodes; n++)
    {
        double node_area(0);

        std::vector<size_t> connected_elements(
            mesh->nod_vector[n]->getConnectedElementIDs());

        for (size_t i = 0; i < connected_elements.size(); i++)
        {
            MeshLib::CElem* Element(mesh->ele_vector[connected_elements[i]]);

            // get nodes of this element
            std::vector<MeshLib::CNode*> ElementNodes;
            Element->GetNodes(ElementNodes);

            // get area of this Element
            // first, get coordinates for each node

            GEOLIB::Point A(ElementNodes[0]->getData());
            GEOLIB::Point B(ElementNodes[1]->getData());
            GEOLIB::Point C(ElementNodes[2]->getData());

            // distances of AB, BC, and AC
            const double a = sqrt((A[0] - B[0]) * (A[0] - B[0]) +
                                  (A[1] - B[1]) * (A[1] - B[1]) +
                                  (A[2] - B[2]) * (A[2] - B[2]));
            const double b = sqrt((C[0] - B[0]) * (C[0] - B[0]) +
                                  (C[1] - B[1]) * (C[1] - B[1]) +
                                  (C[2] - B[2]) * (C[2] - B[2]));
            const double c2 = (A[0] - C[0]) * (A[0] - C[0]) +
                              (A[1] - C[1]) * (A[1] - C[1]) +
                              (A[2] - C[2]) * (A[2] - C[2]);

            // angle AC-BC
            const double cos_gamma = (c2 - a * a - b * b) / (-2 * a * b);

            // Area of tri-element
            const double Area = 0.5 * a * b * sin(acos(cos_gamma));

            node_area +=
                Area / 3.0;  // the third part of the area of each connected
                             // element adds up to the nodal area of n
            total_area += Area / 3.0;
        }

        node_area_vec.push_back(node_area);
    }

    std::cout << "Total surface Area: " << total_area << "\n";
}

MeshLib::CFEMesh* MshEditor::removeMeshNodes(MeshLib::CFEMesh* mesh,
                                             const std::vector<size_t>& nodes)
{
    MeshLib::CFEMesh* new_mesh(new MeshLib::CFEMesh(*mesh));

    // delete nodes and their connected elements and replace them with null
    // pointers
    size_t delNodes = nodes.size();
    for (size_t i = 0; i < delNodes; i++)
    {
        MeshLib::CNode* node = new_mesh->nod_vector[nodes[i]];
        std::vector<size_t> conn_elems = node->getConnectedElementIDs();
        for (size_t j = 0; j < conn_elems.size(); j++)
        {
            delete new_mesh->ele_vector[conn_elems[j]];
            new_mesh->ele_vector[conn_elems[j]] = NULL;
        }
        delete new_mesh->nod_vector[nodes[i]];
        new_mesh->nod_vector[nodes[i]] = NULL;
    }

    // create map to adjust node indices in element vector
    size_t nNodes = new_mesh->nod_vector.size();
    std::vector<int> id_map;
    size_t count = 0;
    for (size_t i = 0; i < nNodes; i++)
    {
        if (new_mesh->nod_vector[i])
        {
            new_mesh->nod_vector[i]->SetIndex(count);
            id_map.push_back(count);
            count++;
        }
        else
            id_map.push_back(-1);
    }

    // erase null pointers from node- and element vectors
    for (std::vector<MeshLib::CElem*>::iterator it =
             new_mesh->ele_vector.begin();
         it != new_mesh->ele_vector.end();)
    {
        if (*it)
            ++it;
        else
            it = new_mesh->ele_vector.erase(it);
    }

    for (std::vector<MeshLib::CNode*>::iterator it =
             new_mesh->nod_vector.begin();
         it != new_mesh->nod_vector.end();)
    {
        if (*it)
            ++it;
        else
            it = new_mesh->nod_vector.erase(it);
    }

    // re-adjust node indices
    size_t nElems = new_mesh->ele_vector.size();
    for (size_t i = 0; i < nElems; i++)
    {
        MeshLib::CElem* elem = new_mesh->ele_vector[i];
        size_t nElemNodes = elem->GetNodesNumber(false);
        for (size_t j = 0; j < nElemNodes; j++)
            elem->SetNodeIndex(j, id_map[elem->GetNodeIndex(j)]);
    }

    return new_mesh;
}

std::vector<GEOLIB::PointWithID*> MshEditor::getSurfaceNodes(
    const MeshLib::CFEMesh& mesh)
{
    std::cout << "Extracting surface nodes..."
              << "\n";
    // Sort points lexicographically
    size_t nNodes(mesh.nod_vector.size());
    std::vector<GEOLIB::PointWithID*> nodes;
    std::vector<size_t> perm;
    for (size_t j(0); j < nNodes; j++)
    {
        nodes.push_back(
            new GEOLIB::PointWithID(mesh.nod_vector[j]->getData(), j));
        perm.push_back(j);
    }
    Quicksort<GEOLIB::PointWithID*>(nodes, 0, nodes.size(), perm);

    // Extract surface points
    double eps(std::numeric_limits<double>::epsilon());
    std::vector<GEOLIB::PointWithID*> surface_pnts;
    for (size_t k(1); k < nNodes; k++)
    {
        const GEOLIB::PointWithID& p0(*(nodes[k - 1]));
        const GEOLIB::PointWithID& p1(*(nodes[k]));
        if (fabs(p0[0] - p1[0]) > eps || fabs(p0[1] - p1[1]) > eps)
            surface_pnts.push_back(nodes[k - 1]);
    }
    // Add last point
    surface_pnts.push_back(nodes[nNodes - 1]);

    return surface_pnts;
}

void MshEditor::sortNodesLexicographically(MeshLib::CFEMesh* mesh)
{
    if (mesh->nodes_are_sorted)
        return;

    const std::vector<MeshLib::CNode*>& mesh_nodes(mesh->getNodeVector());

    const size_t nNodes(mesh_nodes.size());
    std::vector<GEOLIB::PointWithID*> nodes;
    std::vector<size_t> perm;
    for (size_t j(0); j < nNodes; j++)
    {
        nodes.push_back(new GEOLIB::PointWithID(mesh_nodes[j]->getData(),
                                                mesh_nodes[j]->GetIndex()));
        perm.push_back(j);
    }

    Quicksort<GEOLIB::PointWithID*>(nodes, 0, nodes.size(), perm);

    std::vector<size_t> test;
    for (size_t j(0); j < nNodes; j++)
    {
        mesh->sorted_nodes.push_back(nodes[j]->getID());
    }

    // get x-y-coordinate intervals
    const double eps(std::numeric_limits<double>::epsilon());
    mesh->xy_change.push_back(
        -1);  // virtual last top node for reference to first bottom node
    for (size_t k(1); k < nNodes; k++)
    {
        const GEOLIB::PointWithID& p0(*(nodes[k - 1]));
        const GEOLIB::PointWithID& p1(*(nodes[k]));
        if (fabs(p0[0] - p1[0]) > eps || fabs(p0[1] - p1[1]) > eps)
            mesh->xy_change.push_back(k - 1);
    }
    // Add last point
    mesh->xy_change.push_back(nNodes - 1);

    mesh->nodes_are_sorted = true;

    const std::size_t size(nodes.size());
    for (std::size_t j(0); j < size; ++j)
        delete nodes[j];
}

MeshLib::CFEMesh* MshEditor::getMeshSurface(const MeshLib::CFEMesh& mesh)
{
    std::cout << "Extracting mesh surface..."
              << "\n";
    GridAdapter surface;
    const std::vector<GEOLIB::PointWithID*> sfc_points =
        MshEditor::getSurfaceNodes(mesh);
    const size_t nSurfacePoints(sfc_points.size());

    std::vector<GridAdapter::Element*>* elements =
        new std::vector<GridAdapter::Element*>;

    const size_t nElements = mesh.ele_vector.size();
    for (size_t j = 0; j < nElements; j++)
    {
        MeshLib::CElem* elem(mesh.ele_vector[j]);
        std::vector<size_t> elem_nodes;
        bool is_surface(true);
        for (size_t i = 0; i < 4; i++)
        {
            size_t node_index = elem->GetNodeIndex(i);
            bool node_found(false), one_node(true);
            for (size_t k = 0; k < nSurfacePoints; k++)
            {
                if (sfc_points[k]->getID() == node_index)
                {
                    node_found = true;
                    elem_nodes.push_back(k);
                    break;
                }
            }
            if (!node_found)
            {
                if (one_node == true)
                    one_node = false;
                else
                {
                    is_surface = false;
                    break;
                }
            }
        }
        if (is_surface)
        {
            GridAdapter::Element* element = new GridAdapter::Element;
            element->material = 0;
            element->type = MshElemType::TRIANGLE;
            element->nodes = elem_nodes;
            elements->push_back(element);
        }
    }

    std::vector<GEOLIB::Point*>* nodes =
        new std::vector<GEOLIB::Point*>(nSurfacePoints);
    for (size_t j = 0; j < nSurfacePoints; j++)
        //(*nodes)[sfc_points[j]->getID()]=sfc_points[j];
        (*nodes)[j] = sfc_points[j];

    surface.setNodeVector(nodes);
    surface.setElements(elements);

    MeshLib::CFEMesh* sfc_mesh = new MeshLib::CFEMesh(*surface.getCFEMesh());
    return sfc_mesh;
}
