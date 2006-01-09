/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 WW/OK Encapsulation from rf_ele_msh
   last modified
**************************************************************************/
#ifndef msh_node_INC
#define msh_node_INC

#include <cassert>

#include "matrix_class.h"
#include <string>
#include <vector>

// MSHLib
#include "msh_core.h"

namespace MeshLib
{
// Class definition
class CNode : public CCore
{
public:
	int free_surface; //MB ??? mobile
	// The vector to store the representative element index.
	// This can be used to extract the norm of the plane that the element lies on.
	// Establishing this vector is done in the Fluid Momentum
	// since this is bounded by velocity.
	std::vector<long> connected_planes; // PCH

	//GUI control variables
	double patch_area; //OK4310
	bool crossroad; //KR changed to bool // PCH: Make theses privates can be done later on.
	std::vector<long> connected_faces; // BG, 09/2010, necessary for coupling to Eclipse, index of faces where the node is part of it
	std::vector<double> distance_to_connected_faces; // BG, 09/2010,  necessary for coupling to Eclipse

	/** constructor */
	CNode(size_t Index);
	CNode(double x, double y, double z);
	CNode(double const*const coords);
	CNode(size_t Index, double x, double y, double z = 0.0);
	CNode(size_t Index, double const* coordinates);
	CNode(size_t Index, const CNode* parent); //NW
	~CNode()
	{
	}

	// Operator
	void operator =(const CNode& n);
	bool operator ==(const CNode& n);

	// Get functions
	/** \brief const access operator
	 *  The access to the coordinates is like the access to a field.
	 */
	double const& operator[] (size_t idx) const
	{
		assert (idx <= 2);
		return coordinate[idx];
	}
	/** \brief access operator (see book Effektiv C++ programmieren - subsection 1.3.2 ).
	 * \sa const T& operator[] (size_t idx) const
	 */
	double& operator[] (size_t idx)
	{
		return const_cast<double&> (static_cast<const CNode&> (*this)[idx]);
	}

	/** write mesh node coordinates into stream (used from operator<<)
	 * \param os a standard output stream
	 */
	virtual void write (std::ostream &os) const
	{
		os << coordinate[0] << " " << coordinate[1] << " " << coordinate[2] << std::flush;
	}

	 double X() const { return coordinate[0]; };
	 double Y() const { return coordinate[1]; };
	 double Z() const { return coordinate[2]; };

//	 void Coordinates(double *xyz) const
//	 {
//		for (size_t i = 0; i < 3; i++)
//		   xyz[i] = coordinate[i];
//	 }

	// 04/2010 TF
	/**
	 * returns the coordinates of the mesh node
	 * @return the coordinates of this mesh node
	 */
	inline double const* getData() const
	{
		return coordinate;
	}

	// Set functions
	void SetX(double argX)
	{
		coordinate[0] = argX;
	}

	void SetY(double argY)
	{
		coordinate[1] = argY;
	}

	void SetZ(double argZ)
	{
		coordinate[2] = argZ;
	}

	void SetCoordinates(const double* argCoord);

	int GetEquationIndex() const
	{
		return eqs_index;
	}

	void SetEquationIndex(long eqIndex)
	{
		eqs_index = eqIndex;
	}

	// Output
	void Write(std::ostream& os = std::cout) const;

	std::vector<size_t> const & getConnectedElementIDs() const
	{
		return _connected_elements;
	}

	std::vector<size_t> & getConnectedElementIDs()
	{
		return _connected_elements;
	}

	std::vector<size_t> const& getConnectedNodes() const
	{
		return _connected_nodes;
	}

	std::vector<size_t> & getConnectedNodes()
	{
		return _connected_nodes;
	}

	size_t getNumConnectedNodes() const
	{
	   return _connected_nodes.size();
	}

	/*!
	 * \brief Check whether the node is non-ghost
	 * \param id_max_act_l Maximum node ID of the non-ghost node for linear element.
	 * \param id_max_l     Maximum node ID of the node of linear element.
	 * \param id_max_act_h Maximum node ID of the non-ghost node for quadratic element.
	 */
	 bool isNonGhost(const size_t id_max_act_l, const size_t id_max_l,
	              const size_t id_max_act_h)
	{
	   return (index < id_max_act_l) || (index >= id_max_l && index < id_max_act_h);
	}

private:
	double coordinate[3];
	long eqs_index; // renumber
	std::vector<size_t> _connected_nodes; //OK
	std::vector<size_t> _connected_elements;
};

std::ostream& operator<< (std::ostream &os, MeshLib::CNode const &node);

}                                                 // namespace MeshLib


#endif
