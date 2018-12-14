/*
 * SurfaceGrid.h
 *
 *  Created on: Feb 9, 2012
 *      Author: TF
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef SURFACEGRID_H_
#define SURFACEGRID_H_

#include <vector>

// GEOLIB
#include "AxisAlignedBoundingBox.h"

namespace GEOLIB
{
// forward declaration
class Surface;
class Triangle;

class SurfaceGrid : public AABB
{
public:
    SurfaceGrid(Surface const* const sfc);
    virtual ~SurfaceGrid();

    bool isPntInSurface(const double* pnt, double eps = 0) const;

private:
#ifndef NDEBUG
#ifdef DEBUGMESHNODESEARCH
    void writeSurfaceGridData(std::ostream& os) const;
    void writeTrianglesInGridCell(std::size_t i, std::size_t j, std::size_t k,
                                  std::ostream& os) const;
#endif
#endif
    double _step_sizes[3];
    double _inverse_step_sizes[3];
    size_t _n_steps[3];
    std::vector<GEOLIB::Triangle const*>* _triangles_in_grid_box;
};

}  // end namespace GEOLIB

#endif /* SURFACEGRID_H_ */
