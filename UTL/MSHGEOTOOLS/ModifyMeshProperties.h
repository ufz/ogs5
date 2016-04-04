/*
 * ModifyMeshProperties.h
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef MODIFYMESHPROPERTIES_H_
#define MODIFYMESHPROPERTIES_H_

#include "Polygon.h"

namespace MeshLib
{
class CFEMesh;

class ModifyMeshProperties
{
public:
	ModifyMeshProperties(CFEMesh* msh);
	virtual ~ModifyMeshProperties();

	void setMaterial(const GEOLIB::Polygon& polygon, size_t mat_id);
	/**
	 * Method substitutes material ids within the given volume. The volume is described
	 * by three surfaces, the first and second surfaces are the parts of bottom and the top
	 * surface of the mesh that are bounded by the projection of the polygon to the bottom
	 * and the top surface of the mesh. These two parts are linked together via a third
	 * surface.
	 * @param polygon [input] the polygon that is projected
	 * @param old_mat_id [input] the old material id (which is changed)
	 * @param new_mat_id [input] the value of the new material id that is set
	 */
	void substituteMaterialID(GEOLIB::Polygon const& polygon, size_t old_mat_id, size_t new_mat_id);

private:
	CFEMesh* _mesh;
};
}

#endif /* MODIFYMESHPROPERTIES_H_ */
