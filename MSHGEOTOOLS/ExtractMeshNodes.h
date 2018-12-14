/*
 * ExtractMeshNodes.h
 *
 *  Created on: Dec 3, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef EXTRACTMESHNODES_H_
#define EXTRACTMESHNODES_H_

#include <iostream>

// MSH
#include "msh_lib.h"
#include "msh_mesh.h"

// GEO
#include "GEOObjects.h"
#include "PointWithID.h"
#include "Polygon.h"

namespace MeshLib
{
/**
 * This class implements an algorithm to extract mesh node ids from a given
 * *extruded* mesh.
 */
class ExtractMeshNodes
{
public:
    /**
     * constructor - takes a *extruded* mesh
     * @param msh an instance of class CFEMesh
     */
    ExtractMeshNodes(const CFEMesh* msh);
    /**
     * This method first projects a mesh node into the x-y-plane (z=0).
     * Then it checks if this mesh node is within the given polygon
     * (polygon have to be located in the x-y-plane (z=0).
     * The IDs of all (projected) mesh nodes within the polygon will
     * be written to the stream os. For visual control the mesh nodes
     * will be written (as gli points) to the stream gli_out.
     * @param os output stream for IDs
     * @param gli_out output stream for points
     * @param polygon the polygon that have to be located in the x-y-plane (z=0)
     */
    void writeMeshNodeIDs(std::ostream& os, std::ostream& gli_out,
                          const GEOLIB::Polygon& polygon);
    /**
     * This method first projects all mesh nodes into the x-y-plane (z=0).
     * Then it checks if mesh nodes are within the given polygon
     * (polygon have to be located in the x-y-plane (z=0). All the mesh
     * nodes are located within a "cylinder". The method sorts the mesh nodes
     * lexicographical (first by x, then by y and in the end by z). The id of
     * the mesh node with largest z coordinate and identical x and y coordinates
     * will be written to the stream os. For visual control the associated mesh
     * node will be written to the stream gli_out (as gli point).
     * @param os output stream for IDs
     * @param gli_out output stream for points
     * @param polygon the polygon that have to be located in the x-y-plane (z=0)
     */
    void writeTopSurfaceMeshNodeIDs(std::ostream& os, std::ostream& gli_out,
                                    const GEOLIB::Polygon& polygon);

    void writeMesh2DNodeIDAndArea(std::ostream& os, std::ostream& gli_out,
                                  const GEOLIB::Polygon& polygon);

    /**
     * Method computes the ids of mesh nodes that are inside the bounding
     * polygon and not inside a "hole". It writes these ids to the output
     * stream.
     * @param os out put stream
     * @param gli_out mesh nodes as point for visualization purposes
     * @param bounding_polygon the bounding polygon (all mesh nodes inside this
     * polygon are candidates)
     * @param holes mesh nodes insides these polygons are excluded from the
     * output
     */
    void writeMesh2DNodeIDAndArea(std::ostream& os,
                                  std::ostream& gli_out,
                                  const GEOLIB::Polygon& bounding_polygon,
                                  std::vector<GEOLIB::Polygon*> const& holes);

    void writeNearestMeshNodeToPoint(std::ostream& os, std::ostream& gli_out,
                                     GEOLIB::Point const& pnt);

    /**
     * computes the mesh nodes along a polyline belonging to the bottom surface
     * @param ply computation along the polyline ply
     * @param bottom_points the bottom mesh nodes as points
     */
    void getBottomMeshNodesAlongPolylineAsPoints(
        const GEOLIB::Polyline& ply,
        std::vector<GEOLIB::Point*>& bottom_points) const;

    /**
     * computes the mesh nodes along a polyline belonging to the top surface
     * @param polyline computation along the polyline ply
     * @param top_points the top mesh nodes as points
     */
    void getTopMeshNodesAlongPolylineAsPoints(
        const GEOLIB::Polyline& polyline,
        std::vector<GEOLIB::Point*>& top_points) const;

    /**
     * Method computes the polygon to a given polyline that is consisting of the
     * projection of this polyline to the bottom and the top surface of the mesh
     * and the links between these two polylines.
     * @param polyline the ("defining") polyline
     * @param geo_obj geometric objects manager
     * @param name the name of the group of geometric objects
     * @param polygon pointer to the resulting polygon
     *      warning: the pointer to an already existing polygon will be
     * destroyed
     */
    void getPolygonFromPolyline(const GEOLIB::Polyline& polyline,
                                GEOLIB::GEOObjects* geo_obj,
                                std::string const& name,
                                GEOLIB::Polygon*& polygon) const;

    /**
     * get a polyline that is projected to the appropriate "layer" into the mesh
     * @param ply_in polyline that can have arbitrary z coordinates
     * @param geo_obj pointer to a GEOLIB::GEOObjects object (used for returning
     * the data)
     * @param name the name of the new data
     * @param ply_out polyline with points that have the same x- and
     * y-coordinates, the z-coordinate is set from the layer
     * @param layer "layer" of the mesh
     */
    void getProjectedPolylineFromPolyline(GEOLIB::Polyline const& ply_in,
                                          GEOLIB::GEOObjects* geo_obj,
                                          std::string const& name,
                                          GEOLIB::Polyline*& ply_out,
                                          size_t layer = 0) const;

    /**
     * Method searchs the mesh nodes within a given volume described by a
     * polygon, i.e. a part of the mesh volume. The method does not take the
     * \f$z\f$ coordinates for both the polygon and the mesh nodes into
     * account.
     * @param polygon [input] the polygon that is the boundary
     * @param node_ids [output] the ids of the mesh nodes which \f$x\f$ and
     * \f$y\f$ coordinates are inside the boundary polygon
     */
    void getMeshNodeIDsWithinPolygon(GEOLIB::Polygon const& polygon,
                                     std::vector<size_t>& node_ids) const;

private:
    /**
     * This method searchs all mesh nodes with the same x and y coordinates
     * like the polyline points. It returns the mesh nodes as points and the ids
     * within objects of type PointWithID.  The mesh nodes / points are
     * sorted lexicographically.
     * @param polyline the "defining" polyline
     * @param nodes_as_points vector of GEOLIB::Point objects
     */
    void getOrthogonalProjectedMeshNodesAlongPolyline(
        GEOLIB::Polyline const& polyline,
        std::vector<GEOLIB::PointWithID>& nodes_as_points) const;
    const CFEMesh* _msh;
    /**
     * offset for gli point index
     */
    size_t _gli_pnt_offset;
};
}  // namespace MeshLib

#endif /* EXTRACTMESHNODES_H_ */
