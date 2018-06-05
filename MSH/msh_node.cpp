/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 WW/OK Encapsulated from mshlib
**************************************************************************/

#include <iomanip>
#include <string>
#include <vector>

// MSHLib
#include "msh_node.h"

//========================================================================
namespace MeshLib
{
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   08/2011 NW Implementation
**************************************************************************/
CNode::CNode(size_t Index) : CCore(Index), free_surface(-1), patch_area(-1.0), crossroad(0), eqs_index(-1)
{
	coordinate[0] = 0.0;
	coordinate[1] = 0.0;
	coordinate[2] = 0.0;
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
CNode::CNode(size_t Index, double x, double y, double z)
    : CCore(Index), free_surface(-1), patch_area(-1.0), crossroad(false), eqs_index(-1)
{
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
}

CNode::CNode(size_t Index, double const* coordinates)
    : CCore(Index), free_surface(-1), patch_area(-1.0), crossroad(false), eqs_index(-1)
{
	coordinate[0] = coordinates[0];
	coordinate[1] = coordinates[1];
	coordinate[2] = coordinates[2];
}

CNode::CNode(double x, double y, double z)
    : CCore(0), free_surface(-1), patch_area(-1.0), crossroad(false), eqs_index(-1)
{
	coordinate[0] = x;
	coordinate[1] = y;
	coordinate[2] = z;
}

CNode::CNode(double const* const coords) : CCore(0), free_surface(-1), patch_area(-1.0), crossroad(false), eqs_index(-1)
{
	coordinate[0] = coords[0];
	coordinate[1] = coords[1];
	coordinate[2] = coords[2];
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   10/2009 NW Implementation
**************************************************************************/
CNode::CNode(size_t Index, const CNode* parent)
    : CCore(Index), free_surface(-1), patch_area(-1.0), crossroad(false), eqs_index(-1)
{
	coordinate[0] = parent->coordinate[0];
	coordinate[1] = parent->coordinate[1];
	coordinate[2] = parent->coordinate[2];
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CNode::operator=(const CNode& n)
{
	boundary_type = n.boundary_type;
	index = n.index;
	mark = n.mark;
	eqs_index = n.eqs_index;
	coordinate[0] = n.coordinate[0];
	coordinate[1] = n.coordinate[1];
	coordinate[2] = n.coordinate[2];
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
bool CNode::operator==(const CNode& n)
{
	if (index == n.index)
		return true;
	else
		return false;
}
/**************************************************************************
   MSHLib-Method:
   06/2005 WW Implementation
   03/2006 OK patch_area
**************************************************************************/
void CNode::Write(std::ostream& osm) const
{
	osm.setf(std::ios::scientific, std::ios::floatfield);
	std::string deli(" ");
	std::setw(14);
	osm.precision(14);
	osm << index << deli << coordinate[0] << deli << coordinate[1] << deli << coordinate[2] << deli;

	if (patch_area > 0.0)
		osm << "$AREA" << deli << patch_area;
	osm << "\n";
}
/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   06/2005 WW Implementation
**************************************************************************/
void CNode::SetCoordinates(const double* argCoord)
{
	coordinate[0] = argCoord[0];
	coordinate[1] = argCoord[1];
	coordinate[2] = argCoord[2];
}

std::ostream& operator<<(std::ostream& os, MeshLib::CNode const& node)
{
	node.write(os);
	return os;
}

} // namespace MeshLib
