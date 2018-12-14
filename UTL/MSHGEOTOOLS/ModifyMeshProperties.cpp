/*
 * ModifyMeshProperties.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: fischeth
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ModifyMeshProperties.h"

// BaseLib
#include "binarySearch.h"
#include "quicksort.h"

// FEM
#include "matrix_class.h"

// GEO
#include "Point.h"
#include "Polygon.h"

// MSH
#include "msh_elem.h"
#include "msh_mesh.h"
#include "msh_node.h"

// MSHGEOTOOLS
#include "ExtractMeshNodes.h"

// MathLib
#include "AnalyticalGeometry.h"

// STL
#include <fstream>

namespace MeshLib
{
ModifyMeshProperties::ModifyMeshProperties(CFEMesh* msh) : _mesh(msh) {}

ModifyMeshProperties::~ModifyMeshProperties() {}

void ModifyMeshProperties::setMaterial(const GEOLIB::Polygon& polygon,
                                       size_t mat_id)
{
    // get all nodes of mesh
    const std::vector<MeshLib::CNode*>& msh_nodes(_mesh->getNodeVector());

    // *** rotate polygon to xy_plane
    // 1 copy all points
    std::vector<GEOLIB::Point*> polygon_points;
    for (size_t k(0); k < polygon.getNumberOfPoints(); k++)
        polygon_points.push_back(new GEOLIB::Point(*(polygon[k])));
    // 2 rotate points
    MathLib::Vector plane_normal_polygon(0.0, 0.0, 0.0);
    double d_polygon(0.0);
    MathLib::getNewellPlane(polygon_points, plane_normal_polygon, d_polygon);

    //	std::cout << "plane normal: " << plane_normal_polygon << std::endl;
    MathLib::Vector tmp_plane_normal_polygon(
        plane_normal_polygon);  // NW need to keep plane_normal_polygon for
                                // later
    // use
    MathLib::rotatePointsToXY(tmp_plane_normal_polygon, polygon_points);

    // 3 create new polygon
    GEOLIB::Polyline rot_polyline(polygon_points);
    for (size_t k(0); k < polygon.getNumberOfPoints(); k++)
        rot_polyline.addPoint(k);
    rot_polyline.addPoint(0);
    GEOLIB::Polygon rot_polygon(rot_polyline);

    //	std::cout << "Polygon: " << std::endl;
    //	for (size_t k(0); k<polygon.getNumberOfPoints(); k++) {
    //		std::cout << k << ": " << *(polygon[k]) << std::endl;
    //	}
    //	std::cout << std::endl;
    //	std::cout << "rotiertes Polygon: " << std::endl;
    //	for (size_t k(0); k<rot_polygon.getNumberOfPoints(); k++) {
    //		std::cout << k << ": " << *(rot_polygon[k]) << std::endl;
    //	}
    //	std::cout << std::endl;

    // *** rotate mesh nodes to xy-plane
    // 1 copy all mesh nodes to GEOLIB::Points
    std::vector<GEOLIB::Point*> mesh_nodes_as_points;
    for (size_t j(0); j < msh_nodes.size(); j++)
        mesh_nodes_as_points.push_back(
            new GEOLIB::Point(msh_nodes[j]->getData()));
    // 2 rotate the Points
    MathLib::rotatePointsToXY(plane_normal_polygon, mesh_nodes_as_points);

    // get all elements of mesh
    const std::vector<MeshLib::CElem*>& msh_elem(_mesh->getElementVector());

    // *** perform search and modify mesh
    const size_t msh_elem_size(msh_elem.size());
    for (size_t j(0); j < msh_elem_size; j++)
    {
        // indices of nodes of the j-th element
        const Math_Group::vec<long>& nodes_indices(
            msh_elem[j]->GetNodeIndeces());
        //		size_t k;
        //		for (k = 0; k<nodes_indices.Size(); k++) {
        //			if (!
        //rot_polygon.isPntInPolygon(*(mesh_nodes_as_points[nodes_indices[k]])))
        //{ 				break;
        //			}
        //		}
        //
        //		if (k == nodes_indices.Size()) {
        //			msh_elem[j]->setPatchIndex (mat_id);
        //		}

        //		size_t cnt (0);
        //		for (k = 0; k < nodes_indices.Size(); k++)
        //			if
        //(rot_polygon.isPntInPolygon(*(mesh_nodes_as_points[nodes_indices[k]])))
        //				cnt++;
        //
        //		if (cnt >= 2)
        //			msh_elem[j]->setPatchIndex (mat_id);

        double center[3] = {0.0, 0.0, 0.0};
        for (size_t k(0); k < nodes_indices.Size(); k++)
        {
            center[0] += (*(mesh_nodes_as_points[nodes_indices[k]]))[0];
            center[1] += (*(mesh_nodes_as_points[nodes_indices[k]]))[1];
            //			center[2] +=
            //(*(mesh_nodes_as_points[nodes_indices[k]]))[2];
        }
        center[0] /= nodes_indices.Size();
        center[1] /= nodes_indices.Size();
        //		center[2] /= nodes_indices.Size();

        //		std::cout << "center of element " << j << ": " << center[0] << ",
        //" << center[1] << ", " << center[2] <<
        // std::endl;

        if (rot_polygon.isPntInPolygon(center[0], center[1], center[2]))
        {
            msh_elem[j]->setPatchIndex(mat_id);
        }
    }

    for (size_t k(0); k < polygon_points.size(); k++)
        delete polygon_points[k];
    for (size_t j(0); j < mesh_nodes_as_points.size(); j++)
        delete mesh_nodes_as_points[j];
}

void ModifyMeshProperties::substituteMaterialID(GEOLIB::Polygon const& polygon,
                                                size_t old_mat_id,
                                                size_t new_mat_id)
{
    MeshLib::ExtractMeshNodes mesh_node_extractor(_mesh);
    std::vector<size_t> mesh_node_ids;
    mesh_node_extractor.getMeshNodeIDsWithinPolygon(polygon, mesh_node_ids);
    const size_t n_mesh_node_ids(mesh_node_ids.size());
    std::vector<size_t> perm(n_mesh_node_ids);
    // init permutation for sorting
    for (size_t k(0); k < n_mesh_node_ids; k++)
        perm[k] = k;
    // sort - since we want to use binary search
    Quicksort<size_t>(mesh_node_ids, 0, n_mesh_node_ids, perm);

    //#ifndef NDEBUG
    //	std::ofstream test_out ("Points.gli");
    //	test_out << "#POINTS" << std::endl;
    //	FileIO::OGSMeshIO mesh_io;
    //	mesh_io.setMesh(_mesh);
    //	mesh_io.writeMeshNodesAsGLIPnts(mesh_node_ids, test_out);
    //	test_out << "#STOP" << std::endl;
    //#endif

    // get all nodes of the mesh
    const std::vector<MeshLib::CNode*>& mesh_nodes(_mesh->getNodeVector());
    // get all elements of the mesh
    const std::vector<MeshLib::CElem*>& mesh_elements(
        _mesh->getElementVector());

    for (size_t k(0); k < n_mesh_node_ids; k++)
    {
        std::vector<size_t> const& connected_element_ids(
            mesh_nodes[mesh_node_ids[k]]->getConnectedElementIDs());
        for (size_t j(0); j < connected_element_ids.size(); j++)
        {
            if (mesh_elements[connected_element_ids[j]]->GetPatchIndex() ==
                old_mat_id)
            {
                std::vector<size_t> connected_nodes;
                // check if all nodes of element are in the mesh_node_ids vector
                mesh_elements[connected_element_ids[j]]->getNodeIndices(
                    connected_nodes);
                bool all_found(true);
                const size_t n_connected_nodes(connected_nodes.size());
                for (size_t i(0); i < n_connected_nodes && all_found; i++)
                {
                    if (searchElement(connected_nodes[i],
                                      0,
                                      n_mesh_node_ids,
                                      mesh_node_ids) ==
                        std::numeric_limits<size_t>::max())
                    {
                        all_found = false;
                    }
                }
                if (all_found)
                {
                    mesh_elements[connected_element_ids[j]]->setPatchIndex(
                        new_mat_id);
                }
            }
        }
    }
}

}  // end namespace MeshLib
