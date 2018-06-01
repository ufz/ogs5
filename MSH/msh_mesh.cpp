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
   08/2005 WW/OK Encapsulation from rf_ele_msh
   last modified
**************************************************************************/
#include <cfloat>
#include <cmath>
#include <climits>
#include <fstream>
#include <iomanip> //WW
#include <iostream>
#include <sstream>
#include <vector>

// BaseLib
#include "memory.h"
#include "Histogram.h"

#include "FileTools.h"

// GEOLib
//#include "geo_pnt.h"
//#include "geo_ply.h"

// MathLib
#include "MathTools.h"
#include "Vector3.h"
// WW #include "msh_lib.h"
#include "rf_mmp_new.h"

// MSHLib
#include "msh_mesh.h"
#ifdef BENCHMARKING
#include "benchtimer.h"
#endif

#include "rf_random_walk.h"
// For surface integration. WW. 29.02.2009

#include "mathlib.h"
// FEM
#include "fem_ele.h"
#include "files0.h"
#include "ShapeFunctionPool.h"

using FiniteElement::CElement;

// PCSLib
extern std::string GetLineFromFile1(std::ifstream*);

#define noMSH_CHECK

size_t max_dim = 0; // OK411

// class ThreadParameter
//{
// public:
//	ThreadParameter(GEOLIB::Point const* const pnt, size_t start, size_t end,
//	                std::vector<MeshLib::CNode*> const& nod_vector, size_t id) :
//		_pnt(pnt), _start(start), _end(end), _nod_vector(nod_vector), _number(
//		        start), _sqr_dist(std::numeric_limits<double>::max()), _id(id)
//	{
//	}
//
//	GEOLIB::Point const* const _pnt;
//	size_t _start;
//	size_t _end;
//	std::vector<MeshLib::CNode*> const& _nod_vector;
//	size_t _number;
//	double _sqr_dist;
//	size_t _id;
//};

// extern "C" {
// void* threadGetDist(void* ptr)
//{
//	ThreadParameter* thread_param((ThreadParameter*) (ptr));
//	size_t start(thread_param->_start);
//	size_t end(thread_param->_end);
//	std::vector<MeshLib::CNode*> const& nod_vector(thread_param->_nod_vector);
//	GEOLIB::Point const* const pnt(thread_param->_pnt);
//
//	double distmin(MathLib::sqrDist(nod_vector[start]->getData(),
//	                                pnt->getData()));
//	size_t number(start);
//	double sqr_dist(distmin);
//
//	for (size_t i = start + 1; i < end; i++)
//	{
//		sqr_dist = MathLib::sqrDist(nod_vector[i]->getData(), pnt->getData());
//		if (sqr_dist < distmin)
//		{
//			distmin = sqr_dist;
//			number = i;
//		}
//	}
//
//	thread_param->_number = number;
//	thread_param->_sqr_dist = distmin;
//
//	if (number == std::numeric_limits<size_t>::max())
//		return (void*) (-1);
//	return (void*) (number);
//}
//} // end extern "C"

//========================================================================
namespace MeshLib
{
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
**************************************************************************/
CFEMesh::CFEMesh(GEOLIB::GEOObjects* geo_obj, std::string* geo_name)
    : max_mmp_groups(0), msh_max_dim(0), _geo_obj(geo_obj), _geo_name(geo_name), _ele_type(MshElemType::INVALID),
      _n_msh_layer(0), _cross_section(false), _msh_n_lines(0), _msh_n_quads(0), _msh_n_hexs(0), _msh_n_tris(0),
      _msh_n_tets(0), _msh_n_prisms(0), _msh_n_pyras(0), _min_edge_length(1e-3), _search_length(0.0),
      NodesNumber_Linear(0), NodesNumber_Quadratic(0), useQuadratic(false), _axisymmetry(false), _mesh_grid(NULL)
{
	coordinate_system = 1;

	max_ele_dim = 0; // NW
	pcs_name = "NotSpecified"; // WW
	PT = NULL; // WW+TK
	fm_pcs = NULL; // WW
// 1.11.2007 WW
#ifdef NEW_EQS
	sparse_graph = NULL;
	sparse_graph_H = NULL;
#endif
	map_counter = 0; // 21.01.2009 WW
	has_multi_dim_ele = false; // NW
	top_surface_checked = false; // 07.06.2010.  WW

	nodes_are_sorted = false;
}

// Copy-Constructor for CFEMeshes.
// Programming: 2010/11/10 KR
CFEMesh::CFEMesh(CFEMesh const& old_mesh) : PT(NULL), _search_length(old_mesh._search_length), _mesh_grid(NULL)
{
	std::cout << "Copying mesh object ... ";

	// mesh nodes
	size_t nNodes = old_mesh.nod_vector.size();
	for (size_t i = 0; i < nNodes; i++)
		//        MeshLib::CNode* node = new MeshLib::CNode(i);
		//        double coords[3] = { old_mesh.nod_vector[i]->X(), old_mesh.nod_vector[i]->Y(),
		//        old_mesh.nod_vector[i]->Z() };
		//        node->SetCoordinates(coords);
		this->nod_vector.push_back(new MeshLib::CNode(i, old_mesh.nod_vector[i]->getData()));

	//  mesh elements
	size_t nElems = old_mesh.ele_vector.size();
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem = new MeshLib::CElem();
		elem->setElementProperties(old_mesh.ele_vector[i]->GetElementType());
		elem->SetPatchIndex(old_mesh.ele_vector[i]->GetPatchIndex());

		size_t nElemNodes = old_mesh.ele_vector[i]->nodes_index.Size();
		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, old_mesh.ele_vector[i]->GetNodeIndex(j));

		this->ele_vector.push_back(elem);
	}

	pcs_name = "NotSpecified";
	this->setNumberOfMeshLayers(old_mesh.getNumberOfMeshLayers());
	this->ConstructGrid();

	std::cout << "done."
	          << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   07/2005 WW Changes due to the geometry objects
   01/2006 YD Changes for face normal
   03/2010 changed long to size_t to avoid casts
**************************************************************************/
CFEMesh::~CFEMesh(void)
{
	// delete nodes
	size_t nNodes(nod_vector.size());
	for (size_t i = 0; i < nNodes; i++)
		delete nod_vector[i];
	nod_vector.clear();

	// Edges
	size_t nEdges(edge_vector.size());
	for (size_t i = 0; i < nEdges; i++)
		delete edge_vector[i];
	edge_vector.clear();

	// Surface faces
	size_t nFaces(face_vector.size());
	for (size_t i = 0; i < nFaces; i++)
		delete face_vector[i];
	face_vector.clear();

	// Element
	size_t nElems(ele_vector.size());
	for (size_t i = 0; i < nElems; i++)
		delete ele_vector[i];
	ele_vector.clear();

	// normal  YD
	size_t nNormals(face_normal.size());
	for (size_t i = 0; i < nNormals; i++)
		delete face_normal[i];
	face_normal.clear();

	if (PT) // WW
	{
		delete PT; // PCH
		PT = NULL;
	}

// 1.11.2007 WW
#ifdef NEW_EQS
	delete sparse_graph;
	delete sparse_graph_H;
	sparse_graph = NULL;
	sparse_graph_H = NULL;
#endif

	if (_mesh_grid)
	{
		delete _mesh_grid;
		_mesh_grid = NULL;
	}
}

void CFEMesh::setElementType(MshElemType::type type)
{
	_ele_type = type;
}

void CFEMesh::setNumberOfMeshLayers(size_t n_msh_layer)
{
	_n_msh_layer = n_msh_layer;
}

size_t CFEMesh::getNumberOfMeshLayers() const
{
	return _n_msh_layer;
}

bool CFEMesh::hasCrossSection() const
{
	return _cross_section;
}

size_t CFEMesh::getNumberOfLines() const
{
	return _msh_n_lines;
}

size_t CFEMesh::getNumberOfQuads() const
{
	return _msh_n_quads;
}

size_t CFEMesh::getNumberOfHexs() const
{
	return _msh_n_hexs;
}

size_t CFEMesh::getNumberOfTris() const
{
	return _msh_n_tris;
}

size_t CFEMesh::getNumberOfTets() const
{
	return _msh_n_tets;
}

size_t CFEMesh::getNumberOfPrisms() const
{
	return _msh_n_prisms;
}

size_t CFEMesh::getNumberOfPyramids() const
{
	return _msh_n_pyras;
}

double CFEMesh::getMinEdgeLength() const
{
	return _min_edge_length;
}

void CFEMesh::computeMinEdgeLength()
{
	const size_t s_edge_vec(edge_vector.size());
	if (s_edge_vec > 0)
	{
		_min_edge_length = edge_vector[0]->getLength();
		for (size_t k(1); k < s_edge_vec; k++)
		{
			const double kth_edge_length(edge_vector[k]->getLength());
			if (kth_edge_length < _min_edge_length)
				_min_edge_length = kth_edge_length;
		}
	}
}

void CFEMesh::setSearchLength(double len)
{
	_search_length = len;
}

double CFEMesh::getSearchLength() const
{
	return _search_length;
}

void CFEMesh::computeSearchLength(double c)
{
	const size_t n(edge_vector.size());

	if (n == 0)
	{
		std::cerr
		    << "[CFEMesh::computeSearchLength] no edges found for computing _search_length, setting _search_lenght to "
		    << 1e-3 << "\n";
		_search_length = 1e-3;
		return;
	}

	double sum(0);
	double sum_of_sqr(0);

	for (size_t k(0); k < n; k++)
	{
		const double x_k(edge_vector[k]->getLength());
		sum += x_k;
		sum_of_sqr += (x_k * x_k);
	}

	// criterion: mu - c times s, where mu is the average and s is standard deviation
	const double mu(sum / n);

	const double diff = fabs(sum_of_sqr - (sum * sum) / n);
	if (diff < std::numeric_limits<double>::epsilon())
	{
		// In case all edges have the same length
		_search_length = 0.5 * mu;
		return;
	}
	else
	{
		const double s(sqrt(1.0 / (n - 1) * diff));
		while (mu < c * s)
		{
			c *= 0.9;
		}
		_search_length = mu - c * s;
#ifndef NDEBUG
		if (c < 2)
		{
			std::cerr << "[CFEMesh::computeSearchLength] computed _search_length = " << _search_length
			          << ", the average value is: " << mu << ", standard deviation is: " << s << "\n";
		}
#endif
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   07/2005 WW Changes due to the geometry objects
   08/2005 WW/MB Keyword CrossSection
   09/2005 WW 2D-3D flag
   12/2005 OK MAT_TYPE
   03/2011 KR cleaned up code
**************************************************************************/
bool CFEMesh::Read(std::ifstream* fem_file)
{
	std::string line_string;

	bool more_mesh = false; // 12.08.2011. WW

	while (!fem_file->eof())
	{
		getline(*fem_file, line_string);

		// check keywords
		if (line_string.find("#STOP") != std::string::npos)
		{
			more_mesh = false; // 12.08.2011. WW
			break;
		}
		if (line_string.find("#FEM_MSH") != std::string::npos) // 12.08.2011. WW
		{
			more_mesh = true;
			break;
		}

		else if (line_string.find("$PCS_TYPE") != std::string::npos)
			*fem_file >> pcs_name >> std::ws; // WW
		else if (line_string.find("$GEO_NAME") != std::string::npos)
			*fem_file >> geo_name >> std::ws; // WW
		else if (line_string.find("$GEO_TYPE") != std::string::npos)
			*fem_file >> geo_type_name >> geo_name >> std::ws;
		else if (line_string.find("$AXISYMMETRY") != std::string::npos)
			_axisymmetry = true;
		else if (line_string.find("$CROSS_SECTION") != std::string::npos)
			_cross_section = true;
		else if (line_string.find("$NODES") != std::string::npos)
		{
			double x, y, z;
			size_t no_nodes, idx;
			*fem_file >> no_nodes >> std::ws;
			std::string s;
			std::ios::pos_type position = fem_file->tellg();
			for (size_t i = 0; i < no_nodes; i++)
			{
				*fem_file >> idx >> x >> y >> z;
				CNode* newNode(new CNode(idx, x, y, z));
				nod_vector.push_back(newNode);
				position = fem_file->tellg();
				*fem_file >> s;
				if (s.find("$AREA") != std::string::npos)
					*fem_file >> newNode->patch_area;
				else
					fem_file->seekg(position, std::ios::beg);
				*fem_file >> std::ws;
			}
		}
		else if (line_string.find("$ELEMENTS") != std::string::npos)
		{
			size_t no_elements;
			*fem_file >> no_elements >> std::ws;
			for (size_t i = 0; i < no_elements; i++)
			{
				CElem* newElem(new CElem(i));
				newElem->Read(*fem_file);
				setElementType(newElem->geo_type); // CC02/2006
				if (newElem->GetPatchIndex() > max_mmp_groups)
					max_mmp_groups = newElem->GetPatchIndex();
				// NW
				if (newElem->GetDimension() > this->max_ele_dim)
					this->max_ele_dim = newElem->GetDimension();
				ele_vector.push_back(newElem);
			}
		}
		else if (line_string.find("$LAYER") != std::string::npos)
			*fem_file >> _n_msh_layer >> std::ws;
	}

	return more_mesh;
}

/**************************************************************************
   FEMLib-Method: ConnectedElements2Node
   Task: Set elements connected to a node
   Programing:
   04/2007 WW Cut from Construct grid
   03/2011 KR cleaned up code
**************************************************************************/
void CFEMesh::ConnectedElements2Node(bool quadratic)
{
	size_t nNodes(nod_vector.size());
	for (size_t e = 0; e < nNodes; e++)
		nod_vector[e]->getConnectedElementIDs().clear();

	size_t nElems(ele_vector.size());
	for (size_t e = 0; e < nElems; e++)
	{
		CElem* elem = ele_vector[e];
		if (!elem->GetMark())
			continue;

		size_t nElemNodes(static_cast<size_t>(elem->GetNodesNumber(quadratic)));
		for (size_t i = 0; i < nElemNodes; i++)
		{
			MeshLib::CNode* node(nod_vector[elem->GetNodeIndex(i)]);
			node->getConnectedElementIDs().push_back(e);
		}
	}
}

/**************************************************************************
   FEMLib-Method: Construct grid
   Task: Establish topology of a grid
   Programing:
   05/2005 WW Implementation
   02/2006 YD Add 1D line neighbor element set
   01/2010 NW Changed to determine the coordinate system by domain dimensions
   03/2011 KR cleaned up code (at least a little)
**************************************************************************/
void CFEMesh::ConstructGrid()
{
	std::cout << "Executing ConstructGrid() ... " << std::flush;

	bool done;

	Math_Group::vec<CNode*> e_nodes0(20);
	Math_Group::vec<int> Edge_Orientation(15);
	Math_Group::vec<CEdge*> Edges(15);
	Math_Group::vec<CEdge*> Edges0(15);
	Math_Group::vec<CElem*> Neighbors(15);
	Math_Group::vec<CElem*> Neighbors0(15);

	Math_Group::vec<CNode*> e_edgeNodes0(3);
	Math_Group::vec<CNode*> e_edgeNodes(3);

#if !defined(USE_PETSC) // &&! defined(USE_OTHER Parallel solver lib) //WW 06.2013
	NodesNumber_Linear = nod_vector.size();
#endif

	Edge_Orientation = 1;

	// Set neighbors of node
	ConnectedElements2Node();

	// Compute neighbors and edges
	size_t e_size(ele_vector.size());

	// 2011-11-21 TF
	// initializing attributes of objects - why is this not done in the constructor?
	for (size_t e = 0; e < e_size; e++)
	{
		ele_vector[e]->InitializeMembers();
	}

	for (size_t e = 0; e < e_size; e++)
	{
		CElem* elem(ele_vector[e]);
		const Math_Group::vec<long>& node_index(elem->GetNodeIndeces());
		elem->GetNeighbors(Neighbors0);

		size_t nnodes0(elem->nnodes); // Number of nodes for linear element
		for (size_t i = 0; i < nnodes0; i++) // Nodes
			e_nodes0[i] = nod_vector[node_index[i]];

		size_t nFaces = static_cast<size_t>(elem->GetFacesNumber());
		// neighbors
		for (size_t i = 0; i < nFaces; i++) // Faces
		{
			if (Neighbors0[i])
				continue;

			done = false;
			int faceIndex_loc0[10];
			size_t nFaceNodes = static_cast<size_t>(elem->GetElementFaceNodes(i, faceIndex_loc0));
			for (size_t k = 0; k < nFaceNodes; k++) // face nodes
			{
				size_t nConnElems(e_nodes0[faceIndex_loc0[k]]->getConnectedElementIDs().size());
				for (size_t ei = 0; ei < nConnElems; ei++) // elements connected to face node
				{
					size_t ee(static_cast<size_t>(e_nodes0[faceIndex_loc0[k]]->getConnectedElementIDs()[ei]));
					if (ee == e)
						continue;
					CElem* connElem(ele_vector[ee]);
					const Math_Group::vec<long>& node_index_glb(connElem->GetNodeIndeces());
					connElem->GetNeighbors(Neighbors);
					size_t nFacesConnElem = static_cast<size_t>(connElem->GetFacesNumber());

					int faceIndex_loc[10];
					for (size_t ii = 0; ii < nFacesConnElem; ii++) // faces of elements connected to face node
					{
						size_t nFaceNodesConnElem(
						    static_cast<size_t>(connElem->GetElementFaceNodes(ii, faceIndex_loc)));
						if (nFaceNodes != nFaceNodesConnElem)
							continue;
						size_t counter(0);
						for (size_t j = 0; j < nFaceNodes; j++)
						{
							for (size_t jj = 0; jj < nFaceNodesConnElem; jj++)
								if (node_index[faceIndex_loc0[j]] == node_index_glb[faceIndex_loc[jj]])
								{
									counter++;
									break;
								}
						}
						if (counter == nFaceNodesConnElem)
						{
							Neighbors0[i] = connElem;
							Neighbors[ii] = elem;
							connElem->SetNeighbor(ii, elem);
							done = true;
							break;
						}
					}
					if (done)
						break;
				}
				if (done)
					break;
			}
		}
		elem->SetNeighbors(Neighbors0);

		if (elem->geo_type == MshElemType::LINE) // YD
		{
			size_t ii(0);
			for (size_t i = 0; i < nFaces; i++)
			{
				int faceIndex_loc0[10];
				size_t n0 = elem->GetElementFaceNodes(i, faceIndex_loc0);
				for (size_t k = 0; k < n0; k++)
				{
					size_t e_size_l = e_nodes0[faceIndex_loc0[k]]->getConnectedElementIDs().size();
					for (size_t ei = 0; ei < e_size_l; ei++)
					{
						size_t ee = e_nodes0[faceIndex_loc0[k]]->getConnectedElementIDs()[ei];
						CElem* connElem = ele_vector[ee];
						if (e_size_l == 2 && connElem->GetIndex() != elem->GetIndex())
						{
							Neighbors0[i] = connElem;
							Neighbors[ii] = connElem;
							// thisElem->SetNeighbor(ii, thisElem0);   //?? Todo YD
							ii++;
						}
					}
				}
			}
			elem->SetNeighbors(Neighbors0);
		}
		// --------------------------------

		// Edges
		size_t nedges0(elem->GetEdgesNumber());
		elem->GetEdges(Edges0);
		for (size_t i = 0; i < nedges0; i++) // edges
		{
			int edgeIndex_loc0[2];
			elem->GetLocalIndicesOfEdgeNodes(i, edgeIndex_loc0);
			// Check neighbors
			done = false;
			for (size_t k = 0; k < 2; k++) // beginning and end of edge
			{
				size_t nConnElem(e_nodes0[edgeIndex_loc0[k]]->getConnectedElementIDs().size());
				for (size_t ei = 0; ei < nConnElem; ei++) // elements connected to edge node
				{
					size_t ee(e_nodes0[edgeIndex_loc0[k]]->getConnectedElementIDs()[ei]);
					if (ee == e)
						continue;
					CElem* connElem(ele_vector[ee]);
					const Math_Group::vec<long>& node_index_glb(connElem->GetNodeIndeces());
					size_t nedges(connElem->GetEdgesNumber());
					connElem->GetEdges(Edges);
					// Edges of neighbors
					int edgeIndex_loc[2];
					for (size_t ii = 0; ii < nedges; ii++) // edges of element connected to edge node
					{
						connElem->GetLocalIndicesOfEdgeNodes(ii, edgeIndex_loc);

						if ((node_index[edgeIndex_loc0[0]] == node_index_glb[edgeIndex_loc[0]]
						     && node_index[edgeIndex_loc0[1]] == node_index_glb[edgeIndex_loc[1]])
						    || (node_index[edgeIndex_loc0[0]] == node_index_glb[edgeIndex_loc[1]]
						        && node_index[edgeIndex_loc0[1]]
						               == node_index_glb[edgeIndex_loc[0]])) // check if elements share edge

							if (Edges[ii])
							{
								Edges0[i] = Edges[ii];
								Edges[ii]->GetNodes(e_edgeNodes);
								if ((size_t)node_index[edgeIndex_loc0[0]] == e_edgeNodes[1]->GetIndex()
								    && (size_t)node_index[edgeIndex_loc0[1]]
								           == e_edgeNodes[0]->GetIndex()) // check direction of edge
									Edge_Orientation[i] = -1;
								done = true;
								break;
							}
					} //  for(ii=0; ii<nedges; ii++)
					if (done)
						break;
				} // for(ei=0; ei<e_size_l; ei++)
				if (done)
					break;
			} // for(k=0;k<2;k++)
			if (!done) // new edges and new node
			{
				Edges0[i] = new CEdge((long)edge_vector.size());
				Edges0[i]->SetOrder(false);
				e_edgeNodes0[0] = e_nodes0[edgeIndex_loc0[0]];
				e_edgeNodes0[1] = e_nodes0[edgeIndex_loc0[1]];
				e_edgeNodes0[2] = NULL;
				Edges0[i]->SetNodes(e_edgeNodes0);
				edge_vector.push_back(Edges0[i]);
			} // new edges
		} //  for(i=0; i<nedges0; i++)
		//
		// Set edges and nodes
		elem->SetOrder(false);
		elem->SetEdgesOrientation(Edge_Orientation);
		elem->SetEdges(Edges0);
		// Resize is true
		elem->SetNodes(e_nodes0, true);
	} // Over elements

	// Set faces on surfaces and others
	_msh_n_lines = 0; // Should be members of mesh
	_msh_n_quads = 0;
	_msh_n_hexs = 0;
	_msh_n_tris = 0;
	_msh_n_tets = 0;
	_msh_n_prisms = 0;
	_msh_n_pyras = 0;
	for (size_t e = 0; e < e_size; e++)
	{
		CElem* elem(ele_vector[e]);
		switch (elem->GetElementType())
		{
			case MshElemType::LINE:
				_msh_n_lines++;
				break;
			case MshElemType::QUAD:
				_msh_n_quads++;
				break;
			case MshElemType::HEXAHEDRON:
				_msh_n_hexs++;
				break;
			case MshElemType::TRIANGLE:
				_msh_n_tris++;
				break;
			case MshElemType::TETRAHEDRON:
				_msh_n_tets++;
				break;
			case MshElemType::PRISM:
				_msh_n_prisms++;
				break;
			case MshElemType::PYRAMID:
				_msh_n_pyras++;
				break;
			default:
				std::cerr << "CFEMesh::ConstructGrid MshElemType not handled"
				          << "\n";
		}
		// Compute volume meanwhile
		elem->ComputeVolume();

		if (elem->GetElementType() == MshElemType::LINE)
			continue; // line element
		//		thisElem0->GetNodeIndeces(node_index_glb0);
		//		const vec<long>& node_index_glb0 (thisElem0->GetNodeIndeces()); // compiler said: unused variable // TF
		elem->GetNeighbors(Neighbors0);
		size_t m0 = elem->GetFacesNumber();

		// Check face on surface
		for (size_t i = 0; i < m0; i++) // Faces
		{
			if (Neighbors0[i])
				continue;
			CElem* newFace = new CElem((long)face_vector.size(), elem, i);
			//          thisElem0->boundary_type='B';
			elem->no_faces_on_surface++;
			face_vector.push_back(newFace);
			Neighbors0[i] = newFace;
		}
		elem->SetNeighbors(Neighbors0);
	}
#if !defined(USE_PETSC) // &&! defined(USE_OTHER Parallel solver lib) //WW 06.2013
	NodesNumber_Quadratic = (long)nod_vector.size();
#endif
	if ((_msh_n_hexs + _msh_n_tets + _msh_n_prisms + _msh_n_pyras) > 0)
		max_ele_dim = 3;
	else if ((_msh_n_quads + _msh_n_tris) > 0)
		max_ele_dim = 2;
	else
		max_ele_dim = 1;

	// check if this mesh includes multi-dimensional elements
	if (max_ele_dim == 2 && _msh_n_lines > 0) // NW

		this->has_multi_dim_ele = true;
	else if (max_ele_dim == 3 && (_msh_n_quads + _msh_n_tris + _msh_n_lines) > 0)
		this->has_multi_dim_ele = true;

	//----------------------------------------------------------------------
	// Node information
	// 1. Default node index <---> eqs index relationship
	// 2. Coordiate system flag
	double x_sum(0.0), y_sum(0.0), z_sum(0.0);
	Eqs2Global_NodeIndex.clear();
	double xyz_max[3] = // NW
	    {-DBL_MAX, -DBL_MAX, -DBL_MAX};
	double xyz_min[3] = // NW
	    {DBL_MAX, DBL_MAX, DBL_MAX};

	for (size_t e = 0; e < nod_vector.size(); e++)
	{
#if defined(USE_PETSC) // ||defined(USE_OTHER Parallel solver lib) //WW 01.06.2012
		Eqs2Global_NodeIndex.push_back(nod_vector[e]->GetEquationIndex());
#else
		nod_vector[e]->SetEquationIndex(e);
		Eqs2Global_NodeIndex.push_back(nod_vector[e]->GetIndex());
#endif
		double const* const coords(nod_vector[e]->getData());
		x_sum += fabs(coords[0]);
		y_sum += fabs(coords[1]);
		z_sum += fabs(coords[2]);
		if (coords[0] > xyz_max[0])
			xyz_max[0] = coords[0];
		if (coords[1] > xyz_max[1])
			xyz_max[1] = coords[1];
		if (coords[2] > xyz_max[2])
			xyz_max[2] = coords[2];
		if (coords[0] < xyz_min[0])
			xyz_min[0] = coords[0];
		if (coords[1] < xyz_min[1])
			xyz_min[1] = coords[1];
		if (coords[2] < xyz_min[2])
			xyz_min[2] = coords[2];
	}
	double xyz_dim[3]; // NW
	xyz_dim[0] = xyz_max[0] - xyz_min[0];
	xyz_dim[1] = xyz_max[1] - xyz_min[1];
	xyz_dim[2] = xyz_max[2] - xyz_min[2];

	// check dimension of the domain to select appropriate coordinate system
	if (xyz_dim[0] > 0.0 && xyz_dim[1] < MKleinsteZahl && xyz_dim[2] < MKleinsteZahl) // only x-direction
		coordinate_system = 10;
	else if (xyz_dim[1] > 0.0 && xyz_dim[0] < MKleinsteZahl && xyz_dim[2] < MKleinsteZahl) // only y-direction
		coordinate_system = 11;
	else if (xyz_dim[2] > 0.0 && xyz_dim[0] < MKleinsteZahl && xyz_dim[1] < MKleinsteZahl) // only z-direction
		coordinate_system = 12;
	else if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] < MKleinsteZahl)
		coordinate_system = 21; // x & y direction
	else if (xyz_dim[0] > 0.0 && xyz_dim[2] > 0.0 && xyz_dim[1] < MKleinsteZahl)
		coordinate_system = 22; // x & z direction
	else if (xyz_dim[1] > 0.0 && xyz_dim[2] > 0.0 && xyz_dim[0] < MKleinsteZahl)
		coordinate_system = 23; // y & z direction
	else if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] > 0.0) // x, y & z direction
		coordinate_system = 32;

	// 1D in 2D
	if (_msh_n_lines > 0)
	{
		if (xyz_dim[0] > 0.0 && xyz_dim[1] > 0.0 && xyz_dim[2] < MKleinsteZahl)
			coordinate_system = 32;
		if (xyz_dim[0] > 0.0 && xyz_dim[2] > 0.0 && xyz_dim[1] < MKleinsteZahl)
			coordinate_system = 32;
	}

	max_dim = coordinate_system / 10 - 1;
	//----------------------------------------------------------------------
	// Gravity center
	for (size_t e = 0; e < e_size; e++)
		ele_vector[e]->ComputeGravityCenter(); // NW
	//----------------------------------------------------------------------

	// TEST WW
	// For sparse matrix
	ConnectedNodes(false);
	//
	e_nodes0.resize(0);
	//	node_index_glb.resize(0);
	//	node_index_glb0.resize(0);
	Edge_Orientation.resize(0);
	Edges.resize(0);
	Edges0.resize(0);
	Neighbors.resize(0);
	Neighbors0.resize(0);
	e_edgeNodes0.resize(0);
	e_edgeNodes.resize(0);
	std::cout << " done."
	          << "\n";

	// computeSearchLength();
	computeMinEdgeLength();
	setSearchLength(_min_edge_length / 2);
	constructMeshGrid();
}

void CFEMesh::constructMeshGrid()
{
	//#ifndef NDEBUG
	//	std::cout << "CFEMesh::constructMeshGrid() ... " << std::flush;
	//	clock_t start(clock());
	//#endif
	if (_mesh_grid == NULL)
		_mesh_grid = new GEOLIB::Grid<MeshLib::CNode>(this->getNodeVector(), 511);
	//#ifndef NDEBUG
	//	clock_t end(clock());
	//	std::cout << "done, took " << (end-start)/(double)(CLOCKS_PER_SEC) << " s -- " << std::flush;
	//#endif
}

/**************************************************************************
   FEMLib-Method: GenerateHighOrderNodes()
   Task:
   Programing:
   07/2007 WW Implementation
   01/2010 NW Case: a mesh with line elements
**************************************************************************/
void CFEMesh::GenerateHighOrderNodes()
{
	int j, k, ii;
	int nnodes0, nedges0, nedges;
	long e, ei, ee, e_size_l;
	int edgeIndex_loc0[2];
	bool done;
	double x0 = 0.0, y0 = 0.0, z0 = 0.0; // OK411

	// Set neighbors of node. All elements, even in deactivated subdomains, are taken into account here.
	for (e = 0; e < (long)nod_vector.size(); e++)
		nod_vector[e]->getConnectedElementIDs().clear();
	done = false;
	size_t ele_vector_size(ele_vector.size());
	for (size_t e = 0; e < ele_vector_size; e++)
	{
		CElem* thisElem0 = ele_vector[e];
		for (size_t i = 0; i < thisElem0->GetNodesNumber(false); i++)
		{
			done = false;
			long ni = thisElem0->GetNodeIndex(i);
			size_t n_connected_elements(nod_vector[ni]->getConnectedElementIDs().size());
			for (size_t j = 0; j < n_connected_elements; j++)
				if (e == nod_vector[ni]->getConnectedElementIDs()[j])
				{
					done = true;
					break;
				}
			if (!done)
				nod_vector[ni]->getConnectedElementIDs().push_back(e);
		}
	}
	//
	CNode* aNode = NULL;
	Math_Group::vec<CNode*> e_nodes0(20);
	Math_Group::vec<CNode*> e_nodes(20);
	CElem* thisElem0 = NULL;
	CElem* thisElem = NULL;
	CEdge* thisEdge0 = NULL;
	CEdge* thisEdge = NULL;
	//----------------------------------------------------------------------
	// Loop over elements (except for line elements)
	size_t e_size(ele_vector.size());
	for (size_t e = 0; e < e_size; e++)
	{
		thisElem0 = ele_vector[e];
		if (thisElem0->GetElementType() == MshElemType::LINE)
			continue; // NW

		nnodes0 = thisElem0->nnodes; // Number of nodes for linear element
		//      thisElem0->GetNodeIndeces(node_index_glb0);
		for (int i = 0; i < nnodes0; i++) // Nodes
			e_nodes0[i] = thisElem0->GetNode(i);
		// --------------------------------
		// Edges
		nedges0 = thisElem0->GetEdgesNumber();
		// Check if there is any neighbor that has new middle points
		for (int i = 0; i < nedges0; i++)
		{
			thisEdge0 = thisElem0->GetEdge(i);
			thisElem0->GetLocalIndicesOfEdgeNodes(i, edgeIndex_loc0);
			// Check neighbors
			done = false;
			for (k = 0; k < 2; k++)
			{
				e_size_l = (long)e_nodes0[edgeIndex_loc0[k]]->getConnectedElementIDs().size();
				for (ei = 0; ei < e_size_l; ei++)
				{
					size_t ee(e_nodes0[edgeIndex_loc0[k]]->getConnectedElementIDs()[ei]);
					if (ee == e)
						continue;
					thisElem = ele_vector[ee];
					nedges = thisElem->GetEdgesNumber();
					// Edges of neighbors
					for (ii = 0; ii < nedges; ii++)
					{
						thisEdge = thisElem->GetEdge(ii);
						if (*thisEdge0 == *thisEdge)
						{
							aNode = thisEdge->GetNode(2);
							if (aNode) // The middle point exist
							{
								e_nodes0[nnodes0] = aNode;
								nnodes0++;
								done = true;
								break;
							}
						}
					} //  for(ii=0; ii<nedges; ii++)
					if (done)
						break;
				} // for(ei=0; ei<e_size_l; ei++)
				if (done)
					break;
			} // for(k=0;k<2;k++)
			if (!done)
			{
				double const* const pnt0(thisEdge0->GetNode(0)->getData());
				double const* const pnt1(thisEdge0->GetNode(1)->getData());
				aNode = new CNode((long)nod_vector.size(), 0.5 * (pnt0[0] + pnt1[0]), 0.5 * (pnt0[1] + pnt1[1]),
				                  0.5 * (pnt0[2] + pnt1[2]));
				//				aNode->SetX(0.5 * (pnt0[0] + pnt1[0]));
				//				aNode->SetY(0.5 * (pnt0[1] + pnt1[1]));
				//				aNode->SetZ(0.5 * (pnt0[2] + pnt1[2]));
				e_nodes0[nnodes0] = aNode;
				thisEdge0->SetNode(2, aNode);
				nnodes0++;
				nod_vector.push_back(aNode);
			}
		} //  for(i=0; i<nedges0; i++)

		// No neighors or no neighbor has new middle point
		//
		// Quadrilateral
		if (thisElem0->GetElementType() == MshElemType::QUAD)
		{
			x0 = y0 = z0 = 0.0;
			aNode = new CNode((long)nod_vector.size());
			e_nodes0[nnodes0] = aNode;
			nnodes0 = thisElem0->nnodes;
			for (int i = 0; i < nnodes0; i++) // Nodes
			{
				double const* const pnt_i(e_nodes0[i]->getData());
				x0 += pnt_i[0];
				y0 += pnt_i[1];
				z0 += pnt_i[2];
			}
			x0 /= (double)nnodes0;
			y0 /= (double)nnodes0;
			z0 /= (double)nnodes0;
			aNode->SetX(x0);
			aNode->SetY(y0);
			aNode->SetZ(z0);
			nod_vector.push_back(aNode);
		}
		// Set edges and nodes
		thisElem0->SetOrder(true);
		// Resize is true
		thisElem0->SetNodes(e_nodes0, true);
	} // Over elements

	// Setup 1d line elements at the end
	if (_msh_n_lines > 0)
		for (size_t e = 0; e < e_size; e++)
		{
			thisElem0 = ele_vector[e];
			if (thisElem0->GetElementType() != MshElemType::LINE)
				continue;

			nnodes0 = thisElem0->nnodes;
			for (int i = 0; i < nnodes0; i++)
				e_nodes0[i] = thisElem0->GetNode(i);

			done = false;

			for (size_t i = 0; i < thisElem0->GetFacesNumber(); i++)
			{
				thisElem = thisElem0->GetNeighbor(i);
				// look for adjacent solid elements
				if (thisElem->GetElementType() == MshElemType::LINE)
					continue;

				for (j = 0; j < thisElem->nnodes; j++)
					e_nodes[j] = thisElem->GetNode(j);
				nedges = thisElem->GetEdgesNumber();
				// search a edge connecting to this line element
				for (j = 0; j < nedges; j++)
				{
					thisEdge = thisElem->GetEdge(j);
					thisElem->GetLocalIndicesOfEdgeNodes(j, edgeIndex_loc0);
					// Check neighbors
					for (k = 0; k < 2; k++)
					{
						// OK411 CNode *tmp_nod = e_nodes[edgeIndex_loc0[k]];
						e_size_l = (long)e_nodes[edgeIndex_loc0[k]]->getConnectedElementIDs().size();
						for (ei = 0; ei < e_size_l; ei++)
						{
							ee = e_nodes[edgeIndex_loc0[k]]->getConnectedElementIDs()[ei];
							if (ele_vector[ee] != thisElem0)
								continue;
							// the edge is found now
							aNode = thisEdge->GetNode(2);
							if (aNode) // The middle point exist
							{
								e_nodes0[nnodes0] = aNode;
								nnodes0++;
								done = true;
								break;
							}
							if (done)
								break;
						} // for(ei=0; ei<e_size_l; ei++)
						if (done)
							break;
					} // for(k=0;k<2;k++)
					if (done)
						break;
				} //  for(i=0; i<nedges0; i++)
				if (done)
					break;
			}
			if (!done)
			{
				aNode = new CNode((long)nod_vector.size());
				for (int i = 0; i < nnodes0; i++) // Nodes
				{
					double const* const pnt_i(e_nodes0[i]->getData());
					x0 += pnt_i[0];
					y0 += pnt_i[1];
					z0 += pnt_i[2];
				}
				x0 /= (double)nnodes0;
				y0 /= (double)nnodes0;
				z0 /= (double)nnodes0;
				aNode->SetX(x0);
				aNode->SetY(y0);
				aNode->SetZ(z0);
				e_nodes0[nnodes0] = aNode;
				nnodes0++;
				nod_vector.push_back(aNode);
			}
			thisElem0->SetOrder(true);
			thisElem0->SetNodes(e_nodes0, true);
		}
	//
	NodesNumber_Quadratic = (long)nod_vector.size();
	for (e = NodesNumber_Linear; (size_t)e < NodesNumber_Quadratic; e++)
	{
#if !defined(USE_PETSC) // && !defined(USE_OTHER Parallel solver lib)
		nod_vector[e]->SetEquationIndex(e);
#endif
		Eqs2Global_NodeIndex.push_back(nod_vector[e]->GetIndex());
	}
	for (size_t e = 0; e < e_size; e++)
	{
		thisElem0 = ele_vector[e];
		for (int i = thisElem0->nnodes; i < thisElem0->nnodesHQ; i++)
		{
			done = false;
			aNode = thisElem0->GetNode(i);
			for (k = 0; k < (int)aNode->getConnectedElementIDs().size(); k++)
				if (e == aNode->getConnectedElementIDs()[k])
				{
					done = true;
					break;
				}
			if (!done)
				aNode->getConnectedElementIDs().push_back(e);
		}
	}

	// For sparse matrix
	ConnectedNodes(true);
	ConnectedElements2Node(true);
	//
	e_nodes0.resize(0);

	// Test	WW
	/*
	   fstream n_out;
	   n_out.open("node.txt", ios::out );
	   for(e=0; e<NodesNumber_Quadratic; e++)
	   nod_vector[e]->Write(n_out);
	   n_out.close();
	 */
}

/**************************************************************************
   FEMLib-Method:
   Task:  Renumbering nodes corresponding to the activate of elements
   Programing:
   09/2005 WW Implementation
   05/2007 WW 1D in 2D
**************************************************************************/
void CFEMesh::FillTransformMatrix()
{
	CElem* elem = NULL;
	// PCH
	CRFProcess* m_pcs = PCSGet("FLUID_MOMENTUM");
	//
	if ((_msh_n_hexs + _msh_n_tets + _msh_n_prisms + _msh_n_pyras) == ele_vector.size())
		return;
	else if (coordinate_system != 32 && !this->has_multi_dim_ele)
	{
		if (m_pcs)
			; // Need to do FillTransformMatrix	// PCH
		else
			return;
	}
	bool tilted = false;
	if (coordinate_system == 32 || coordinate_system == 21 || coordinate_system == 22)
		tilted = true;
	if (!tilted)
		return;
	for (size_t i = 0; i < ele_vector.size(); i++)
	{
		elem = ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			if (coordinate_system == 21 || coordinate_system == 22)
			{
				if (elem->GetElementType() == MshElemType::LINE)
					elem->FillTransformMatrix();
			}
			else if (elem->GetElementType() == MshElemType::LINE || elem->GetElementType() == MshElemType::QUAD
			         || elem->GetElementType() == MshElemType::TRIANGLE)
				elem->FillTransformMatrix();
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:  Renumbering nodes corresponding to the activiate of elements
   Programing:
   05/2005 WW Implementation
**************************************************************************/
void CFEMesh::RenumberNodesForGlobalAssembly()
{
	CElem* elem = NULL;
	Eqs2Global_NodeIndex.clear();
	size_t nNodes(nod_vector.size());
	for (size_t l = 0; l < nNodes; l++)
		nod_vector[l]->SetEquationIndex(-1);

	size_t el_0(0);
	// Lower order
	size_t nElems(ele_vector.size());
	for (size_t l = 0; l < nElems; l++)
	{
		elem = ele_vector[l];
		if (elem->GetMark()) // Marked for use
		{
			int nVert(elem->GetVertexNumber());
			for (int i = 0; i < nVert; i++)
			{
				if (elem->nodes[i]->GetEquationIndex() < 0)
				{
					elem->nodes[i]->SetEquationIndex(el_0);
					Eqs2Global_NodeIndex.push_back(elem->nodes[i]->GetIndex());
					el_0++;
				}
				else
					continue;
			}
		}
	}
	size_t el(el_0);
	if (!getOrder())
	{
		NodesNumber_Linear = el_0;
		NodesNumber_Quadratic = el;
		return;
	}
	// High order
	for (size_t l = 0; l < nElems; l++)
	{
		elem = ele_vector[l];
		if (elem->GetMark()) // Marked for use
		{
			int nElemNodes(elem->GetNodesNumber(true));
			for (int i = elem->GetVertexNumber(); i < nElemNodes; i++)
			{
				if (elem->nodes[i]->GetEquationIndex() < 0)
				{
					elem->nodes[i]->SetEquationIndex(el);
					Eqs2Global_NodeIndex.push_back(elem->nodes[i]->GetIndex());
					el++;
				}
				else
					continue;
			}
		}
	}
	NodesNumber_Linear = el_0;
	NodesNumber_Quadratic = el;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   07/2005 WW Write by member methods of geometry objects.
   12/2005 OK MAT_TYPE
**************************************************************************/
/*
   void CFEMesh::Write(std::fstream*fem_msh_file, bool append_stop) const
   {
   long i;
   //--------------------------------------------------------------------
   //KEYWORD
   *fem_msh_file << "#FEM_MSH" << "\n";
   //--------------------------------------------------------------------
   // PCS
   *fem_msh_file << " $PCS_TYPE" << "\n";
   *fem_msh_file << "  ";
   *fem_msh_file << pcs_name << "\n";
   //--------------------------------------------------------------------
   // MAT
   if (geo_name.size() > 0)
   {
   *fem_msh_file << " $GEO_TYPE" << "\n";
   *fem_msh_file << "  ";
   //OK10_4310
   *fem_msh_file << geo_type_name << " " << geo_name << "\n";
   }
   //--------------------------------------------------------------------
   // NODES
   *fem_msh_file << " $NODES" << "\n";
   *fem_msh_file << "  ";
   //WW
   *fem_msh_file << GetNodesNumber(false) << "\n";
   for (i = 0; i < (long) nod_vector.size(); i++)
   nod_vector[i]->Write(*fem_msh_file);     //WW
   //--------------------------------------------------------------------
   // ELEMENTS
   *fem_msh_file << " $ELEMENTS" << "\n";
   *fem_msh_file << "  ";
   *fem_msh_file << (long) ele_vector.size() << "\n";
   for (i = 0; i < (long) ele_vector.size(); i++)
   {
   ele_vector[i]->SetIndex(i);              //20.01.06 WW/TK
   ele_vector[i]->WriteIndex(*fem_msh_file);//WW
   }
   //--------------------------------------------------------------------
   *fem_msh_file << " $LAYER" << "\n";
   *fem_msh_file << "  ";
   *fem_msh_file << _n_msh_layer << "\n";
   //--------------------------------------------------------------------
   if (append_stop)
   *fem_msh_file << "#STOP";
   }
 */

/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt den nahliegenden existierenden Knoten
   Programing:
   03/2010 TF implementation based on long CFEMesh::GetNODOnPNT(CGLPoint*m_pnt)
   by OK, WW
   05/2012 WW Find node in subdomains for over-lapped DDC
**************************************************************************/
long CFEMesh::GetNODOnPNT(const GEOLIB::Point* const pnt) const
{
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2012
	long node_id = -1;

	const size_t id_act_l_max = static_cast<size_t>(getNumNodesLocal());
	const size_t id_act_h_min = GetNodesNumber(false);
	const size_t id_act_h_max = getLargestActiveNodeID_Quadratic();

	double sqr_dist = 0.0;
	double distmin = getMinEdgeLength() / 10.0;
	if (distmin < 0.)
		distmin = DBL_EPSILON;

	for (size_t i = 0; i < id_act_l_max; i++)
	{
		sqr_dist = MathLib::sqrDist(nod_vector[i]->getData(), pnt->getData());
		if (sqrt(sqr_dist) < distmin)
		{
			node_id = i;
			break;
		}
	}

	if (!useQuadratic)
		return node_id;

	for (size_t i = id_act_h_min; i < id_act_h_max; i++)
	{
		sqr_dist = MathLib::sqrDist(nod_vector[i]->getData(), pnt->getData());
		if (sqrt(sqr_dist) < distmin)
		{
			node_id = i;
			break;
		}
	}
	return node_id;

#else
	MeshLib::CNode const* const node(_mesh_grid->getNearestPoint(pnt->getData()));
	return node->GetIndex();
#endif // END: if use_petsc

	//	const size_t nodes_in_usage(static_cast<size_t> (NodesInUsage()));
	//	double sqr_dist(0.0), distmin(MathLib::sqrDist (nod_vector[0]->getData(), pnt->getData()));
	//	size_t number(0);
	//	for (size_t i = 1; i < nodes_in_usage; i++)
	//	{
	//		sqr_dist = MathLib::sqrDist (nod_vector[i]->getData(), pnt->getData());
	//		if (sqr_dist < distmin)
	//		{
	//			distmin = sqr_dist;
	//			number = i;
	//		}
	//	}
	//
	//	if (number != node_idx) {
	//		double const*const data0(nod_vector[node_idx]->getData());
	//		double const*const data1(nod_vector[number]->getData());
	//		std::cout << "pnt: "<< *pnt << " mesh grid node " << nod_vector[node_idx]->GetIndex() << ": "
	//			<< data0[0] << " " << data0[1] << " " << data0[2] <<
	//			", mesh grid node (old algorithm) " << nod_vector[number]->GetIndex() << ": "
	//			<< data1[0] << " " << data1[1] << " " << data1[2] << "\n";
	//		std::cout << "bbx: " << _mesh_grid->getMinPoint() << " x " << _mesh_grid->getMaxPoint() << "\n";
	//		size_t coords[3];
	//		_mesh_grid->getGridCoords(pnt->getData(), coords);
	//		std::cout << "grid coords: " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
	//		double llf[3], urb[3];
	//		_mesh_grid->getGridCornerPoints(pnt->getData(), llf, urb);
	//		std::cout << "local bbx: " << llf[0] << " " << llf[1] << " " << llf[2] << " x " << urb[0] << " " << urb[1] <<
	//" " << urb[2] << "\n";
	//	}
	//
	//	return number;
}

/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt das nahliegende Element
   Programing:
   03/2010 TF implementation based on long CFEMesh::GetNODOnPNT(CGLPoint*m_pnt)
   by MB
**************************************************************************/
long CFEMesh::GetNearestELEOnPNT(const GEOLIB::Point* const pnt) const
{
	long nextele(-1);
	double dist(std::numeric_limits<double>::max()), dist1;

	for (size_t i = 0; i < ele_vector.size(); i++)
	{
		double const* center(ele_vector[i]->GetGravityCenter());
		dist1 = 0.0;
		for (size_t k(0); k < 3; k++)
			dist1 += (center[k] - (*pnt)[k]) * (center[k] - (*pnt)[k]);
		if (dist1 < dist)
		{
			dist = dist1;
			nextele = i;
		}
	}
	return nextele;
}

// WW. (x1-x0).(x2-x0)
inline double dotProduction(const double* x1, const double* x2, const double* x0)
{
	return (x1[0] - x0[0]) * (x2[0] - x0[0]) + (x1[1] - x0[1]) * (x2[1] - x0[1]) + (x1[2] - x0[2]) * (x2[2] - x0[2]);
}

/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt den nahliegenden existierenden Knoten
   Programing:
   03/2005 OK Implementation (based on ExecuteSourceSinkMethod11 by CT)
   07/2005 WW Node object is replaced
   10/2005 OK test
   03/2010 TF adaption to new data GEO-structures, changed the algorithm
   05/3013 WW Add restriction for the ply for the sources term
**************************************************************************/
void CFEMesh::GetNODOnPLY(const GEOLIB::Polyline* const ply,
                          std::vector<size_t>& msh_nod_vector,
                          bool automatic,
                          double eps)
{
	msh_nod_vector.clear();

	// search for nodes along polyline in previous computed polylines
	std::vector<MeshNodesAlongPolyline>::const_iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); it++)
	{
		if (it->getPolyline() == ply)
		{
			const std::vector<size_t> node_ids(it->getNodeIDs());

			size_t n_valid_nodes(0);
			if (useQuadratic)
				n_valid_nodes = node_ids.size();
			else
				n_valid_nodes = it->getNumberOfLinearNodes();

			for (size_t k(0); k < n_valid_nodes; k++)
				msh_nod_vector.push_back(node_ids[k]);
#ifndef NDEBUG
			std::string ply_name;
			if (!getGEOObjects()->getPolylineVecObj(*(getProjectName()))->getNameOfElement(ply, ply_name))
			{
				ply_name = "unknown-ply";
			}
			std::cout << "[DEBUG-INFO] access " << msh_nod_vector.size() << " buffered nodes for polyline " << ply_name
			          << "\n";
#endif
			return;
		}
	}

	// compute nodes (and supporting points) along polyline
	double search_radius(this->getMinEdgeLength()); // getSearchLength());
	if (!automatic)
		search_radius = eps;
	_mesh_nodes_along_polylines.push_back(MeshNodesAlongPolyline(ply, this, search_radius));
	const std::vector<size_t> node_ids(
	    _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size() - 1].getNodeIDs());

	size_t n_valid_nodes(0);
	if (useQuadratic)
		n_valid_nodes = node_ids.size();
	else
		n_valid_nodes = _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size() - 1].getNumberOfLinearNodes();

	for (size_t k(0); k < n_valid_nodes; k++)
		msh_nod_vector.push_back(node_ids[k]);
#ifndef NDEBUG
	std::string ply_name;
	if (!getGEOObjects()->getPolylineVecObj(*(getProjectName()))->getNameOfElement(ply, ply_name))
	{
		ply_name = "unknown-ply";
	}
	std::cout << "[DEBUG-INFO] computed " << n_valid_nodes << " nodes for polyline " << ply_name << " - "
	          << NodesInUsage() << "\n";

//	std::string fname ("MeshNodeIDsAlongPolyline"+ply_name+".txt");
//	std::ofstream os (fname.c_str());
//	for (size_t k(0); k < n_valid_nodes; k++)
//		os << node_ids[k] << "\n";
//	os.close();
#endif
}

const MeshNodesAlongPolyline& CFEMesh::GetMeshNodesAlongPolyline(const GEOLIB::Polyline* const ply)
{
	// search for nodes along polyline in previous computed polylines
	std::vector<MeshNodesAlongPolyline>::const_iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); it++)
		if (it->getPolyline() == ply)
			return *it;
	// compute nodes (and supporting points for interpolation) along polyline
	double search_radius(this->getMinEdgeLength()); // getSearchLength());
	_mesh_nodes_along_polylines.push_back(MeshNodesAlongPolyline(ply, this, search_radius));
	return _mesh_nodes_along_polylines[_mesh_nodes_along_polylines.size() - 1];
}

void CFEMesh::getPointsForInterpolationAlongPolyline(const GEOLIB::Polyline* const ply, std::vector<double>& points)
{
	// search for nodes along polyline in previous computed polylines
	std::vector<MeshNodesAlongPolyline>::const_iterator it(_mesh_nodes_along_polylines.begin());
	for (; it != _mesh_nodes_along_polylines.end(); it++)
		if (it->getPolyline() == ply)
		{
			// copy points from object into vector
			for (size_t k(0); k < it->getDistOfProjNodeFromPlyStart().size(); k++)
				points.push_back(it->getDistOfProjNodeFromPlyStart()[k]);
			return;
		}

	// compute nodes (and points according the nodes) along polyline
	double search_radius(this->getMinEdgeLength()); // getSearchLength());
	_mesh_nodes_along_polylines.push_back(MeshNodesAlongPolyline(ply, this, search_radius));
	// copy supporting points from object into vector
	for (size_t k(0); k < (_mesh_nodes_along_polylines.back()).getDistOfProjNodeFromPlyStart().size(); k++)
		points.push_back(it->getDistOfProjNodeFromPlyStart()[k]);
}

void CFEMesh::GetNODOnPLY(const GEOLIB::Polyline* const ply,
                          std::vector<long>& msh_nod_vector,
#ifdef USE_PETSC
                          const bool for_s_term,
#else
                          const bool,
#endif
                          bool automatic,
                          double search_radius)
{
	msh_nod_vector.clear(); // JOD 2014-11-10
	//----------------------------------------------------------------------
	std::vector<size_t> tmp_msh_node_vector;
	if (automatic)
	{
		search_radius = this->getMinEdgeLength() / 2;
	}
	GetNODOnPLY(ply, tmp_msh_node_vector, automatic, search_radius);

#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 08.2014
	if (for_s_term)
	{
		for (size_t k(0); k < tmp_msh_node_vector.size(); k++)
			msh_nod_vector.push_back(tmp_msh_node_vector[k]);
	}
	else
	{
		const size_t start_act_node_h = NodesNumber_Linear;
		const size_t end_act_node_h = NodesNumber_Linear
		                              + static_cast<size_t>(loc_NodesNumber_Quadratic - loc_NodesNumber_Linear);
		for (size_t k(0); k < tmp_msh_node_vector.size(); k++)
		{
			const size_t n_id = nod_vector[tmp_msh_node_vector[k]]->GetIndex();
			if (n_id < static_cast<size_t>(loc_NodesNumber_Linear)
			    || (n_id >= start_act_node_h && n_id < end_act_node_h))
				msh_nod_vector.push_back(tmp_msh_node_vector[k]);
		}
	}
#else
	for (size_t k(0); k < tmp_msh_node_vector.size(); k++)
		msh_nod_vector.push_back(tmp_msh_node_vector[k]);
#endif
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   04/2005 OK
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnSFC(Surface* m_sfc, std::vector<long>& msh_nod_vector, const bool for_s_term)
{
	msh_nod_vector.clear();
	//----------------------------------------------------------------------
	switch (m_sfc->type)
	{
		//....................................................................
		case 0: // Surface polygon
			GetNODOnSFC_PLY(m_sfc, msh_nod_vector, for_s_term);
			break;
		case 1: // TIN
			if (!m_sfc->TIN)
				return;
			GetNODOnSFC_TIN(m_sfc, msh_nod_vector);
			break;
		//....................................................................
		case 2: // 2 vertical polylines
			GetNODOnSFC_Vertical(m_sfc, msh_nod_vector);
			break;
		case 3: // projection on xy plane (all mesh points above and below the surface) //MB
			GetNODOnSFC_PLY_XY(m_sfc, msh_nod_vector);
			break;
		//....................................................................
		case 100:
			GetNodesOnCylindricalSurface(m_sfc, msh_nod_vector);
			break;
		case 4: // layer polyline, all z
			GetNODOnSFC_PLY_Z(m_sfc, msh_nod_vector);
			break;
	}
}

/**************************************************************************
   MSHLib-Method:
   Task: Get nodes on plane surface
   Programing:
   03/2010 TF
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnSFC(const GEOLIB::Surface* sfc,
                          std::vector<size_t>& msh_nod_vector,
#ifdef USE_PETSC
                          const bool for_s_term
#else
                          const bool
#endif
                          ) const
{
	msh_nod_vector.clear();

#ifdef TIME_MEASUREMENT
	clock_t begin, end;
	std::cout << "[CFEMesh::GetNODOnSFC] init SurfaceGrid ... " << std::flush;
	begin = clock();
#endif
	const_cast<GEOLIB::Surface*>(sfc)->initSurfaceGrid();
#ifdef TIME_MEASUREMENT
	end = clock();
	std::cout << "done, took " << (end - begin) / (double)(CLOCKS_PER_SEC) << " s"
	          << "\n";

	std::cout << "[CFEMesh::GetNODOnSFC] search with new algorithm ... " << std::flush;
	begin = clock();
#endif
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2012
	if (for_s_term)
	{
		const size_t nodes_in_usage = (size_t)NodesInUsage();
		for (size_t j(0); j < nodes_in_usage; j++)
		{
			if (sfc->isPntInBV((nod_vector[j])->getData(), _search_length / 2.0))
			{
				if (sfc->isPntInSfc((nod_vector[j])->getData(), _search_length / 2.0))
				{
					msh_nod_vector.push_back(nod_vector[j]->GetIndex());
				}
			}
		}
	}
	else
	{
		const size_t id_act_l_max = static_cast<size_t>(getNumNodesLocal());
		const size_t id_act_h_min = GetNodesNumber(false);
		const size_t id_act_h_max = getLargestActiveNodeID_Quadratic();

		for (size_t j = 0; j < id_act_l_max; j++)
		{
			if (sfc->isPntInBV((nod_vector[j])->getData(), _search_length / 2.0))
			{
				if (sfc->isPntInSfc((nod_vector[j])->getData(), _search_length / 2.0))
				{
					msh_nod_vector.push_back(nod_vector[j]->GetIndex());
				}
			}
		}

		if (!useQuadratic)
			return;

		for (size_t j = id_act_h_min; j < id_act_h_max; j++)
		{
			if (sfc->isPntInBV((nod_vector[j])->getData(), _search_length / 2.0))
			{
				if (sfc->isPntInSfc((nod_vector[j])->getData(), _search_length / 2.0))
				{
					msh_nod_vector.push_back(nod_vector[j]->GetIndex());
				}
			}
		}
	}

#else
	const size_t nodes_in_usage((size_t)NodesInUsage());
	for (size_t j(0); j < nodes_in_usage; j++)
	{
		if (sfc->isPntInBV((nod_vector[j])->getData(), _search_length * 0.375))
		{
			if (sfc->isPntInSfc((nod_vector[j])->getData(), _search_length * 0.375))
			{
				msh_nod_vector.push_back(nod_vector[j]->GetIndex());
			}
		}
	}
#endif

#ifdef TIME_MEASUREMENT
	end = clock();
	std::cout << "done, took " << (end - begin) / (double)(CLOCKS_PER_SEC) << " s, " << msh_nod_vector.size()
	          << "nodes found"
	          << "\n";
#endif
}

void CFEMesh::findNodesInPolygon(const double area_orig, const double tol, const size_t start_id, const size_t end_id,
                                 const CGLPolyline* ply, std::vector<long>& node_id_vector) const
{
	double x1[3];
	double x2[3];
	const size_t np = ply->point_vector.size();
	for (size_t j = start_id; j < end_id; j++)
	{
		double area_calculated = 0.0;
		for (size_t i = 0; i < np; i++)
		{
			CGLPoint* point = ply->point_vector[i];
			x1[0] = point->x;
			x1[1] = point->y;
			x1[2] = point->z;

			size_t k = i + 1;
			if (i == np - 1)
				k = 0;
			point = ply->point_vector[k];
			x2[0] = point->x;
			x2[1] = point->y;
			x2[2] = point->z;

			area_calculated += fabs(ComputeDetTri(x1, nod_vector[j]->getData(), x2));
		}

		if (fabs(area_orig - area_calculated) < tol)
			node_id_vector.push_back(nod_vector[j]->GetIndex());
	}
};

/**************************************************************************
   MSHLib-Method:
   Task: Get nodes on plane surface by comparing the area of polygon computed
   by triangles, which are formed by node and the gravity center
   with edges of polygon, respectively
   Programing:
   09/2004 WW Implementation
   04/2005 OK MSH
   07/2005 WW Node object is replaced
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnSFC_PLY(Surface const* m_sfc,
                              std::vector<long>& msh_nod_vector,
#ifdef USE_PETSC
                              const bool for_s_term
#else
                              const bool
#endif
                              ) const
{
	long i;
	int nPointsPly = 0;
	double gC[3], p1[3], p2[3];
	double Area1;
	double Tol = m_sfc->epsilon;
	CGLPolyline* m_ply = NULL;
	std::vector<CGLPolyline*>::const_iterator p_ply(m_sfc->polyline_of_surface_vector.begin()); // CC
	// Init
	msh_nod_vector.clear();

	// nodes close to first polyline
	while (p_ply != m_sfc->polyline_of_surface_vector.end())
	{
		m_ply = *p_ply;
		nPointsPly = (int)m_ply->point_vector.size();
		if (m_ply->point_vector.front() == m_ply->point_vector.back())
			nPointsPly -= 1;
		//....................................................................
		// Gravity center of this polygon
		for (i = 0; i < 3; i++)
			gC[i] = 0.0;
		for (i = 0; i < nPointsPly; i++)
		{
			gC[0] += m_ply->point_vector[i]->x;
			gC[1] += m_ply->point_vector[i]->y;
			gC[2] += m_ply->point_vector[i]->z;
		}
		for (i = 0; i < 3; i++)
			gC[i] /= (double)nPointsPly;
		//....................................................................
		// Area of this polygon by the grativity center
		Area1 = 0.0;
		for (i = 0; i < nPointsPly; i++)
		{
			p1[0] = m_ply->point_vector[i]->x;
			p1[1] = m_ply->point_vector[i]->y;
			p1[2] = m_ply->point_vector[i]->z;
			if (i < nPointsPly - 1)
			{
				p2[0] = m_ply->point_vector[i + 1]->x;
				p2[1] = m_ply->point_vector[i + 1]->y;
				p2[2] = m_ply->point_vector[i + 1]->z;
			}
			else
			{
				p2[0] = m_ply->point_vector[0]->x;
				p2[1] = m_ply->point_vector[0]->y;
				p2[2] = m_ply->point_vector[0]->z;
			}
			Area1 += fabs(ComputeDetTri(p1, gC, p2));
		}
//....................................................................
// Check nodes by comparing area
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2012
		if (for_s_term)
		{
			findNodesInPolygon(Area1, Tol, 0, NodesInUsage(), m_ply, msh_nod_vector);
		}
		else
		{
			findNodesInPolygon(Area1, Tol, 0, getNumNodesLocal(), m_ply, msh_nod_vector);

			if (useQuadratic)
			{
				findNodesInPolygon(Area1, Tol, GetNodesNumber(false), getLargestActiveNodeID_Quadratic(), m_ply,
				                   msh_nod_vector);
			}
		}
#else
		findNodesInPolygon(Area1, Tol, 0, NodesInUsage(), m_ply, msh_nod_vector);
#endif
		p_ply++;
	}
}

/**************************************************************************
   MSHLib-Method:
   Task: Get nodes on plane surface by comparing the area of polygon computed
   by triangles, which are formed by node and the gravity center
   with edges of polygon, respectively
   Programing:
   08/2005 MB based on GetNODOnSFC_PLY
   03/2009 WW Case only search the specified nodes given in msh_nod_vector
   03/2009 WW Efficiency improvement
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnSFC_PLY_XY(Surface* m_sfc, std::vector<long>& msh_nod_vector, bool givenNodesOnSurface)
{
	long i, k;
	size_t j;
	int nPointsPly = 0;
	double gC[3], p1[3], p2[3];
	double Area1, Area2;
	double Tol = m_sfc->epsilon;
	CGLPolyline* m_ply = NULL;
	std::vector<CGLPolyline*>::iterator p_ply; // CC
	// Init
	//----19.03.2009. WW
	CNode* a_node = NULL;
	std::vector<long> temp_v;
	//
	if (givenNodesOnSurface)
	{
		temp_v.resize((long)msh_nod_vector.size());
		temp_v = msh_nod_vector;
	}
	p1[2] = p2[2] = 0.;
	//----19.03.2009. WW
	//
	msh_nod_vector.clear();
	//----------------------------------------------------------------------
	// nodes close to first polyline
	// CC
	p_ply = m_sfc->polyline_of_surface_vector.begin();
	// CC
	while (p_ply != m_sfc->polyline_of_surface_vector.end())
	{
		m_ply = *p_ply;
		nPointsPly = (int)m_ply->point_vector.size();
		if (m_ply->point_vector.front() == m_ply->point_vector.back())
			nPointsPly -= 1;
		//....................................................................
		// Grativity center of this polygon
		for (i = 0; i < 3; i++)
			gC[i] = 0.0;
		for (i = 0; i < nPointsPly; i++)
		{
			gC[0] += m_ply->point_vector[i]->x;
			gC[1] += m_ply->point_vector[i]->y;
		}
		for (i = 0; i < 3; i++)
			gC[i] /= (double)nPointsPly;
		//....................................................................
		// Area of this polygon by the grativity center
		Area1 = 0.0;
		for (i = 0; i < nPointsPly; i++)
		{
			p1[0] = m_ply->point_vector[i]->x;
			p1[1] = m_ply->point_vector[i]->y;
			k = i + 1;
			if (i == nPointsPly - 1)
				k = 0;
			p2[0] = m_ply->point_vector[k]->x;
			p2[1] = m_ply->point_vector[k]->y;
			Area1 += fabs(ComputeDetTri(p1, gC, p2));
		}
		//....................................................................
		// Check nodes by comparing area
		//------- 19.03.2009. WW ----------------
		if (givenNodesOnSurface)
			for (j = 0; j < temp_v.size(); j++)
			{
				Area2 = 0.0;
				a_node = nod_vector[temp_v[j]];
				//               gC[0] = a_node->X();
				//               gC[1] = a_node->Y();
				for (i = 0; i < nPointsPly; i++)
				{
					p1[0] = m_ply->point_vector[i]->x;
					p1[1] = m_ply->point_vector[i]->y;
					k = i + 1;
					if (i == nPointsPly - 1)
						k = 0;
					p2[0] = m_ply->point_vector[k]->x;
					p2[1] = m_ply->point_vector[k]->y;
					Area2 += fabs(ComputeDetTri(p1, a_node->getData(), p2));
				}
				if (fabs(Area1 - Area2) < Tol)
					msh_nod_vector.push_back(a_node->GetIndex());
			}
		//-----------------------------------------------
		else
			for (j = 0; j < NodesInUsage(); j++)
			{
				Area2 = 0.0;
				a_node = nod_vector[j]; // 19.03.2009. WW
				//
				//               gC[0] = a_node->X();
				//               gC[1] = a_node->Y();
				for (i = 0; i < nPointsPly; i++)
				{
					p1[0] = m_ply->point_vector[i]->x;
					p1[1] = m_ply->point_vector[i]->y;
					k = i + 1;
					if (i == nPointsPly - 1)
						k = 0;
					p2[0] = m_ply->point_vector[k]->x;
					p2[1] = m_ply->point_vector[k]->y;
					Area2 += fabs(ComputeDetTri(p1, a_node->getData(), p2));
				}
				if (fabs(Area1 - Area2) < Tol)
					msh_nod_vector.push_back(a_node->GetIndex());
			}
		p_ply++;
	}
	//
	if (givenNodesOnSurface) // 19.03.2009. WW
		temp_v.clear();
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   04/2005 OK
   07/2005 WW Node object is replaced
   04/2006 TK new method
   01/2010 NW use epsilon specified in GEO as tolerance
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnSFC_TIN(Surface* m_sfc, std::vector<long>& msh_nod_vector)
{
	double dist;
	double tolerance = 0.001;
	double min_mesh_dist = 0.0;
	double tri_point1[3], tri_point2[3], tri_point3[3], checkpoint[3];
	double sfc_min[3]
	    = {m_sfc->TIN->Triangles[0]->x[0], m_sfc->TIN->Triangles[0]->y[0], m_sfc->TIN->Triangles[0]->z[0]};
	double sfc_max[3]
	    = {m_sfc->TIN->Triangles[0]->x[0], m_sfc->TIN->Triangles[0]->y[0], m_sfc->TIN->Triangles[0]->z[0]};

	CTriangle* m_triangle = NULL;
	//----------------------------------------------------------------------
	// Create Bounding BOX = MIN/MAX of X/Y/Z
	//----------------------------------------------------------------------
	// Loop over all generated triangles of surface
	for (size_t m = 0; m < m_sfc->TIN->Triangles.size(); m++)
	{
		m_triangle = m_sfc->TIN->Triangles[m];
		tri_point1[0] = m_triangle->x[0];
		tri_point1[1] = m_triangle->y[0];
		tri_point1[2] = m_triangle->z[0];
		tri_point2[0] = m_triangle->x[1];
		tri_point2[1] = m_triangle->y[1];
		tri_point2[2] = m_triangle->z[1];
		tri_point3[0] = m_triangle->x[2];
		tri_point3[1] = m_triangle->y[2];
		tri_point3[2] = m_triangle->z[2];
		if (m == 0)
		{
			sfc_min[0] = tri_point1[0];
			sfc_min[1] = tri_point1[1];
			sfc_min[2] = tri_point1[2];
			sfc_max[0] = tri_point1[0];
			sfc_max[1] = tri_point1[1];
			sfc_max[2] = tri_point1[2];
			if (tri_point1[0] < sfc_min[0])
				sfc_min[0] = tri_point1[0];
			if (tri_point2[0] < sfc_min[0])
				sfc_min[0] = tri_point2[0];
			if (tri_point3[0] < sfc_min[0])
				sfc_min[0] = tri_point3[0];
			if (tri_point1[0] > sfc_max[0])
				sfc_max[0] = tri_point1[0];
			if (tri_point2[0] > sfc_max[0])
				sfc_max[0] = tri_point2[0];
			if (tri_point3[0] > sfc_max[0])
				sfc_max[0] = tri_point3[0];
			if (tri_point1[1] < sfc_min[1])
				sfc_min[1] = tri_point1[1];
			if (tri_point2[1] < sfc_min[1])
				sfc_min[1] = tri_point2[1];
			if (tri_point3[1] < sfc_min[1])
				sfc_min[1] = tri_point3[1];
			if (tri_point1[1] > sfc_max[1])
				sfc_max[1] = tri_point1[1];
			if (tri_point2[1] > sfc_max[1])
				sfc_max[1] = tri_point2[1];
			if (tri_point3[1] > sfc_max[1])
				sfc_max[1] = tri_point3[1];
			if (tri_point1[2] < sfc_min[2])
				sfc_min[2] = tri_point1[2];
			if (tri_point2[2] < sfc_min[2])
				sfc_min[2] = tri_point2[2];
			if (tri_point3[2] < sfc_min[2])
				sfc_min[2] = tri_point3[2];
			if (tri_point1[2] > sfc_max[2])
				sfc_max[2] = tri_point1[2];
			if (tri_point2[2] > sfc_max[2])
				sfc_max[2] = tri_point2[2];
			if (tri_point3[2] > sfc_max[2])
				sfc_max[2] = tri_point3[2];
		}
		else
		{
			if (tri_point1[0] < sfc_min[0])
				sfc_min[0] = tri_point1[0];
			if (tri_point2[0] < sfc_min[0])
				sfc_min[0] = tri_point2[0];
			if (tri_point3[0] < sfc_min[0])
				sfc_min[0] = tri_point3[0];
			if (tri_point1[0] > sfc_max[0])
				sfc_max[0] = tri_point1[0];
			if (tri_point2[0] > sfc_max[0])
				sfc_max[0] = tri_point2[0];
			if (tri_point3[0] > sfc_max[0])
				sfc_max[0] = tri_point3[0];
			if (tri_point1[1] < sfc_min[1])
				sfc_min[1] = tri_point1[1];
			if (tri_point2[1] < sfc_min[1])
				sfc_min[1] = tri_point2[1];
			if (tri_point3[1] < sfc_min[1])
				sfc_min[1] = tri_point3[1];
			if (tri_point1[1] > sfc_max[1])
				sfc_max[1] = tri_point1[1];
			if (tri_point2[1] > sfc_max[1])
				sfc_max[1] = tri_point2[1];
			if (tri_point3[1] > sfc_max[1])
				sfc_max[1] = tri_point3[1];
			if (tri_point1[2] < sfc_min[2])
				sfc_min[2] = tri_point1[2];
			if (tri_point2[2] < sfc_min[2])
				sfc_min[2] = tri_point2[2];
			if (tri_point3[2] < sfc_min[2])
				sfc_min[2] = tri_point3[2];
			if (tri_point1[2] > sfc_max[2])
				sfc_max[2] = tri_point1[2];
			if (tri_point2[2] > sfc_max[2])
				sfc_max[2] = tri_point2[2];
			if (tri_point3[2] > sfc_max[2])
				sfc_max[2] = tri_point3[2];
		}
	}
	//----------------------------------------------------------------------
	// Create Local Search Vector
	// Only nodes inside searching box
	//----------------------------------------------------------------------

	CFEMesh* m_msh_aux(new CFEMesh(_geo_obj, _geo_name));

	tolerance = m_sfc->epsilon; // NW
	// NW commented out below. Minimum edge length doesn't work for some cases
	////Loop over all edges
	//     for(i=0;i<(long)edge_vector.size();i++)
	//     {
	//         if (j==0 && i==0){
	//           min_mesh_dist = edge_vector[i]->Length();
	//         }
	//         else{
	//           if (min_mesh_dist  > edge_vector[i]->Length())
	//               min_mesh_dist =  edge_vector[i]->Length();
	//         }
	//     }
	//     tolerance = min_mesh_dist;

	// TF expand bounding box with epsilon environment
	for (size_t k(0); k < 3; k++)
	{
		sfc_min[k] -= tolerance;
		sfc_max[k] += tolerance;
	}

	// Loop over all mesh nodes
	const size_t n_nodes(static_cast<size_t>(NodesInUsage()));
	for (size_t i = 0; i < n_nodes; i++)
	{
		double const* const pnt(nod_vector[i]->getData());
		if ((pnt[0] >= sfc_min[0] && pnt[0] <= sfc_max[0]) && (pnt[1] >= sfc_min[1] && pnt[1] <= sfc_max[1])
		    && (pnt[2] >= sfc_min[2] && pnt[2] <= sfc_max[2]))
			m_msh_aux->nod_vector.push_back(new CNode(i, pnt));
	}

	//----------------------------------------------------------------------
	// Search preselected Nodes within TIN Triangles
	//----------------------------------------------------------------------
	for (size_t m = 0; m < m_sfc->TIN->Triangles.size(); m++)
	{
		m_triangle = m_sfc->TIN->Triangles[m];
		tri_point1[0] = m_triangle->x[0];
		tri_point1[1] = m_triangle->y[0];
		tri_point1[2] = m_triangle->z[0];
		tri_point2[0] = m_triangle->x[1];
		tri_point2[1] = m_triangle->y[1];
		tri_point2[2] = m_triangle->z[1];
		tri_point3[0] = m_triangle->x[2];
		tri_point3[1] = m_triangle->y[2];
		tri_point3[2] = m_triangle->z[2];
		// Loop over all preselected mesh nodes
		for (size_t i = 0; i < m_msh_aux->nod_vector.size(); i++)
		{
			double const* const pnt_i(m_msh_aux->nod_vector[i]->getData());
			dist = MCalcDistancePointToPlane(pnt_i, tri_point1, tri_point2, tri_point3);
			if (dist <= tolerance && dist >= -tolerance)
				AngleSumPointInsideTriangle(checkpoint, tri_point1, tri_point2, tri_point3, min_mesh_dist);
		}
	}

	//----------------------------------------------------------------------
	// Identify the preselected nodes of the search vector and copy to msh_nod_vector
	// TODO: Works only for one mesh!!!
	//----------------------------------------------------------------------
	int index;
	// Loop over selected nodes
	for (size_t i = 0; i < m_msh_aux->nod_vector.size(); i++)
	{
		index = m_msh_aux->nod_vector[i]->GetIndex();
		double const* const pnt_index(nod_vector[index]->getData());
		double const* const pnt_i(m_msh_aux->nod_vector[i]->getData());
		if (index < (int)nod_vector.size())
			if ((m_msh_aux->nod_vector[i]->GetIndex() == nod_vector[index]->GetIndex()) && (pnt_i[0] == pnt_index[0])
			    && (pnt_i[1] == pnt_index[1])
			    && (pnt_i[2] == pnt_index[2]))
				msh_nod_vector.push_back(nod_vector[index]->GetIndex());
	}

	//----------------------------------------------------------------------
	// Delete Search Vector at the end of fem_msh_vector
	// TODO: Proper delete by MSHDelete!!!
	//----------------------------------------------------------------------
	for (size_t i = 0; i < m_msh_aux->nod_vector.size(); i++)
		delete m_msh_aux->nod_vector[i];
}

/**************************************************************************
   GeoLib-Method:
   Task: Get nodes on cylindrical surface by comparing the area of
   triangles form by nodes and two axis points
   Programing:
   10/2004 WW Implementation
   05/2005 WW Transplant to this object from GEOLIB
   04/2006 WW Case of quadratic elements
   last modification:
**************************************************************************/
void CFEMesh::GetNodesOnCylindricalSurface(Surface* m_sfc, std::vector<long>& NodesS)
{
	int l, nf;
	long i, j, m, fnode;
	const int nNodes = NodesInUsage();
	int faceIndex_loc[10];
	double p1[3], p2[3];
	double dist, R, dc1, dc2;
	CElem* elem = NULL;
	CNode* cnode = NULL;
	NodesS.clear();
	// m_sfc->epsilon = 1.0e-6;
	p1[0] = m_sfc->polygon_point_vector[0]->x;
	p1[1] = m_sfc->polygon_point_vector[0]->y;
	p1[2] = m_sfc->polygon_point_vector[0]->z;

	p2[0] = m_sfc->polygon_point_vector[1]->x;
	p2[1] = m_sfc->polygon_point_vector[1]->y;
	p2[2] = m_sfc->polygon_point_vector[1]->z;

	dist = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1])
	            + (p1[2] - p2[2]) * (p1[2] - p2[2]));

	// Check nodes by comparing area
	for (j = 0; j < nNodes; j++)
	{
		cnode = nod_vector[j];
		double const* const gc(cnode->getData());
		//         gC[0] = cnode->X();
		//         gC[1] = cnode->Y();
		//         gC[2] = cnode->Z();

		dc1 = (p2[0] - p1[0]) * (gc[0] - p1[0]) + (p2[1] - p1[1]) * (gc[1] - p1[1]) + (p2[2] - p1[2]) * (gc[2] - p1[2]);
		dc2 = (p2[0] - p1[0]) * (gc[0] - p2[0]) + (p2[1] - p1[1]) * (gc[1] - p2[1]) + (p2[2] - p1[2]) * (gc[2] - p2[2]);
		if (dc1 < 0.0)
			continue;
		if (dc2 > 0.0)
			continue;

		R = 2.0 * fabs(ComputeDetTri(p1, gc, p2)) / dist;

		if (fabs(R - m_sfc->Radius) < m_sfc->epsilon)
			NodesS.push_back(cnode->GetIndex());
	}
	bool done = false;
	int counter = 0;
	int hs;
	long NodesS_size = (long)NodesS.size();
	//
	if (useQuadratic)
		// Face elements are only in quadrilaterals or triangles
		for (i = 0; i < NodesS_size; i++)
		{
			cnode = nod_vector[NodesS[i]];
			for (j = 0; j < (long)cnode->getConnectedElementIDs().size(); j++)
			{
				elem = ele_vector[cnode->getConnectedElementIDs()[j]];
				for (size_t k = 0; k < elem->GetFacesNumber(); k++)
				{
					nf = elem->GetElementFaceNodes(k, faceIndex_loc);
					counter = 0;
					hs = (int)(nf / 2);
					for (l = 0; l < hs; l++) // loop over face vertices
					{
						fnode = elem->GetNodeIndex(faceIndex_loc[l]);
						for (m = 0; m < NodesS_size; m++)
							if (fnode == NodesS[m])
								counter++;
					}
					//
					if (counter == hs) // face is found on surface

						for (l = hs; l < nf; l++) // loop over face vertices
						{
							fnode = elem->GetNodeIndex(faceIndex_loc[l]);
							done = false;
							for (m = 0; m < (long)NodesS.size(); m++)
								if (fnode == NodesS[m])
								{
									done = true;
									break;
								}
							if (!done)
								NodesS.push_back(fnode);
						}
				}
			}
		}
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   04/2005 OK
   last modification:
   07/2010 TF small modifications concerning coding style
**************************************************************************/
void CFEMesh::GetNODOnSFC_Vertical(Surface* m_sfc, std::vector<long>& msh_nod_vector)
{
	long* nodes_array = NULL;
	long no_nodes = 0;

	// nodes close to first polyline
	std::vector<CGLPolyline*>::const_iterator p_ply(m_sfc->polyline_of_surface_vector.begin());
	while (p_ply != m_sfc->polyline_of_surface_vector.end())
	{
		// OK41 nodes_array = m_polyline->MSHGetNodesCloseXY(&no_nodes);
		// CC 10/05
		nodes_array = MSHGetNodesClose(&no_nodes, *p_ply);
		break;
	}

	// using triangles
	CGLPolyline* m_polyline1 = NULL;
	CGLPolyline* m_polyline2 = NULL;
	p_ply = m_sfc->polyline_of_surface_vector.begin();
	while (p_ply != m_sfc->polyline_of_surface_vector.end())
	{
		m_polyline1 = *p_ply;
		++p_ply;
		m_polyline2 = *p_ply;
		break;
	}
	long no_points = (long)m_polyline1->point_vector.size();

	CGLPoint m_node;
	double xp[3], yp[3], zp[3];

	long i, j;
	for (j = 0; j < no_nodes; j++)
		// OK m_node.x = GetNodeX(nodes_array[j]);
		// OK m_node.y = GetNodeY(nodes_array[j]);
		// OK m_node.z = GetNodeZ(nodes_array[j]);
		for (i = 0; i < no_points - 1; i++)
		{
			// first triangle of quad
			xp[0] = m_polyline1->point_vector[i]->x;
			yp[0] = m_polyline1->point_vector[i]->y;
			zp[0] = m_polyline1->point_vector[i]->z;
			xp[1] = m_polyline1->point_vector[i + 1]->x;
			yp[1] = m_polyline1->point_vector[i + 1]->y;
			zp[1] = m_polyline1->point_vector[i + 1]->z;
			xp[2] = m_polyline2->point_vector[i]->x;
			yp[2] = m_polyline2->point_vector[i]->y;
			zp[2] = m_polyline2->point_vector[i]->z;
			// CC 10/05
			if (m_node.IsInsideTriangle(xp, yp, zp))
				msh_nod_vector.push_back(nodes_array[j]);
			// second triangle of quad
			xp[0] = m_polyline2->point_vector[i]->x;
			yp[0] = m_polyline2->point_vector[i]->y;
			zp[0] = m_polyline2->point_vector[i]->z;
			xp[1] = m_polyline2->point_vector[i + 1]->x;
			yp[1] = m_polyline2->point_vector[i + 1]->y;
			zp[1] = m_polyline2->point_vector[i + 1]->z;
			xp[2] = m_polyline1->point_vector[i + 1]->x;
			yp[2] = m_polyline1->point_vector[i + 1]->y;
			zp[2] = m_polyline1->point_vector[i + 1]->z;
			// CC 10/05
			if (m_node.IsInsideTriangle(xp, yp, zp))
				msh_nod_vector.push_back(nodes_array[j]);
		} // no_points
	// no_nodes
}
/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt den nahliegenden existierenden Knoten
   Programing:
   03/2005 OK Implementation (based on ExecuteSourceSinkMethod11 by CT)
   last modification:
**************************************************************************/
void CFEMesh::GetNODOnPLY_XY(CGLPolyline* m_ply, std::vector<long>& msh_nod_vector)
{
	long j, k, l;
	double pt1[3], line1[3], line2[3]; // WW , pt0[3];
	double mult_eps = 1.0;
	double dist1p, dist2p, *length, laenge;
	long anz_relevant = 0;
	typedef struct
	{
		long knoten;
		double laenge;
	} INFO;
	INFO* relevant = NULL;
	int weiter;
	double w1, w2;
	long knoten_help;
	double laenge_help;
	m_ply->getSBuffer().clear();
	m_ply->getIBuffer().clear();
	msh_nod_vector.clear();
	//
	length = (double*)Malloc(sizeof(double) * (long)m_ply->point_vector.size());
	// WW pt0[0] = m_ply->point_vector[0]->x;
	// WW pt0[1] = m_ply->point_vector[0]->y;
	// WW pt0[2] = 0.0;
	/* */
	for (k = 0; k < (long)m_ply->point_vector.size() - 1; k++)
	{
		line1[0] = m_ply->point_vector[k]->x;
		line1[1] = m_ply->point_vector[k]->y;
		line1[2] = 0.0;
		line2[0] = m_ply->point_vector[k + 1]->x;
		line2[1] = m_ply->point_vector[k + 1]->y;
		line2[2] = 0.0;
		length[k] = MCalcDistancePointToPoint(line2, line1);
	}
	/* Wiederholen bis zumindest ein Knoten gefunden wurde */
	while (anz_relevant == 0)
	{
		/* Schleife ueber alle Knoten */
		for (j = 0; j < (long)nod_vector.size(); j++)
		{
			double const* const pnt_j(nod_vector[j]->getData());
			/* Schleife ueber alle Punkte des Polygonzuges */
			for (k = 0; k < (long)m_ply->point_vector.size() - 1; k++)
			{
				/* ??? */
				// TF can be moved out of the k loop in case MCalcDistancePointToLine,
				// MCalcDistancePointToPoint and MCalcDistancePointToPoint have no side effects
				pt1[0] = pnt_j[0];
				pt1[1] = pnt_j[1];
				pt1[2] = 0.0;
				line1[0] = m_ply->point_vector[k]->x;
				line1[1] = m_ply->point_vector[k]->y;
				line1[2] = 0.0;
				line2[0] = m_ply->point_vector[k + 1]->x;
				line2[1] = m_ply->point_vector[k + 1]->y;
				line2[2] = 0.0;
				/* Ist der Knoten nah am Polygonabschnitt? */
				if (MCalcDistancePointToLine(pt1, line1, line2) <= mult_eps * m_ply->epsilon)
				{
					/* Im folgenden wird mit der Projektion weitergearbeitet */
					MCalcProjectionOfPointOnLine(pt1, line1, line2, pt1);
					/* Abstand des Punktes zum ersten Punkt des Polygonabschnitts */
					dist1p = MCalcDistancePointToPoint(line1, pt1);
					/* Abstand des Punktes zum zweiten Punkt des Polygonabschnitts */
					dist2p = MCalcDistancePointToPoint(line2, pt1);
					/* Ist der Knoten innerhalb des Intervalls? */
					/* bis rf3807: if ((length[k] - dist1p - dist2p + MKleinsteZahl)/(length[k] + dist1p + dist2p +
					 * MKleinsteZahl) > -MKleinsteZahl){ */
					if ((dist1p + dist2p - length[k]) <= mult_eps * m_ply->epsilon)
					{
						// For boundara conditions. WW
						m_ply->getSBuffer().push_back(dist1p);
						m_ply->getIBuffer().push_back(k);
						anz_relevant++;
						/* Feld anpassen */
						// nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
						relevant = (INFO*)Realloc(relevant, sizeof(INFO) * anz_relevant);
						/* Ablegen von Knotennummer und Position */
						// nodes_all[anz_relevant-1] = j;
						msh_nod_vector.push_back(j);
						/* Position ermitteln */
						laenge = 0.;
						for (l = 0; l < k; l++)
							laenge += length[l];
						/* Ablegen von Knotennummer und Position */
						relevant[anz_relevant - 1].knoten = j;
						relevant[anz_relevant - 1].laenge = laenge + dist1p;
						/* Suche am Polygon abbrechen, naechster Knoten */
						k = (long)m_ply->point_vector.size();
					}
				} /* endif */
			} /* Ende Schleife ueber Polygonabschnitte */
		} /* Ende Schleife ueber Knoten */
		if (anz_relevant == 0)
			mult_eps *= 2.;
	} /* Ende Schleife Wiederholungen */
	if (mult_eps > 1.)
		std::cout << "!!! Epsilon increased in sources!"
		          << "\n";
	/* Schleife ueber alle Knoten; sortieren nach Reihenfolge auf dem Abschnitt (zyklisches Vertauschen, sehr lahm)*/
	do
	{
		weiter = 0;
		for (k = 0; k < anz_relevant - 1; k++)
		{
			w1 = relevant[k].laenge;
			w2 = relevant[k + 1].laenge;
			if (w1 > w2) /* Die Eintraege vertauschen */
			{
				knoten_help = relevant[k].knoten;
				laenge_help = relevant[k].laenge;
				relevant[k].knoten = relevant[k + 1].knoten;
				relevant[k].laenge = relevant[k + 1].laenge;
				relevant[k + 1].knoten = knoten_help;
				relevant[k + 1].laenge = laenge_help;
				weiter = 1;
			}
		}
	} while (weiter);
	relevant = (INFO*)Free(relevant);
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   09/2005 OK Implementation
   09/2005 OK Epsilon
   10/2005 OK Delete existing layer polylines
   02/2006 CC polyline id
**************************************************************************/
void CFEMesh::CreateLayerPolylines(CGLPolyline* m_ply)
{
	long i;
	CGLPolyline* m_polyline = NULL;
	char layer_number[3];
	//---------------------------------------------------------------------
	// Delete existing layer polylines
	std::string ply_lay_name = m_ply->getName() + "_L";
	for (int k = 0; k < (int)polyline_vector.size(); k++)
	{
		m_polyline = polyline_vector[k];
		if (m_polyline->getName().find(ply_lay_name) != std::string::npos)
		{
			GEORemovePLY(m_polyline);
			// GEORemovePolyline(polyline_vector.begin()+(k-l));
			k--;
		}
	}
	//---------------------------------------------------------------------
	//
	std::vector<long> ply_nod_vector;
	std::vector<long> ply_nod_vector_dummy;
	GetNODOnPLY_XY(m_ply, ply_nod_vector);
	// nodes = MSHGetNodesCloseXY(&no_nodes); //OK41
	long nodes_per_layer = (long)nod_vector.size() / (_n_msh_layer + 1);
	int ply_nod_vector_layer = (int)ply_nod_vector.size() / (_n_msh_layer + 1);
	//---------------------------------------------------------------------
	// Create layer polylines
	// polyline id CC8888---------------------------------------------------
	long size = 0;
	CGLPolyline* ms_polyline = NULL;
	long number_of_polylines = (long)polyline_vector.size();
	if (number_of_polylines == 0)
		size = 0;
	else
	{
		std::vector<CGLPolyline*>::iterator ps = polyline_vector.begin();
		while (ps != polyline_vector.end())
		{
			ms_polyline = *ps;
			++ps;
		}
		size = ms_polyline->getID() + 1;
	}

	if (ply_nod_vector_layer < 1)
		return;

	//---------------------------------------------------------------------
	// Create layer polylines
	//......................................................................
	m_polyline = new CGLPolyline;
	sprintf(layer_number, "%ld", 0L);
	// CString names =  m_ply->name + "_L" + layer_number;
	// m_polyline->name = names;

	// 10/2010 TF
	//	m_polyline->ply_name = m_ply->getName(); //CC/TK8888
	//	m_polyline->ply_name.append("_L"); //CC/TK8888
	//	m_polyline->ply_name.append(layer_number);//CC/TK8888
	std::string tmp_name(m_ply->getName() + "_L");
	tmp_name.append(layer_number);

	//	m_polyline->name = m_polyline->ply_name.data();//CC/TK8888 // TF
	m_polyline->setName(tmp_name); // TF

	m_polyline->setDataType(1);
	m_polyline->setID(size); // CC8888 / TF
	m_polyline->epsilon = m_ply->epsilon; // OK
	//	m_polyline->ply_data = m_polyline->getName () + ".ply";//CC

	for (i = 0; i < ply_nod_vector_layer; i++)
	{
		CGLPoint* point(new CGLPoint(nod_vector[ply_nod_vector[i]]->getData()));
		m_polyline->point_vector.push_back(point);
	}
	m_polyline->SetPointOrderByDistance(m_ply->point_vector[0]);
	polyline_vector.push_back(m_polyline);
	m_polyline->WritePointVector(m_polyline->getName());
	m_polyline->WriteTecplot(" "); // OK41
	//......................................................................
	for (size_t j = 1; j < (_n_msh_layer + 1); j++)
	{
		m_polyline = new CGLPolyline;
		sprintf(layer_number, "%ld", static_cast<long>(j));
		//		m_polyline->ply_name = m_ply->name.data();//CC/TK8888
		//		m_polyline->ply_name.append("_L");//CC/TK8888
		//		m_polyline->ply_name.append(layer_number);//CC/TK8888
		//		m_polyline->name = m_polyline->ply_name.data();//CC/TK8888

		std::string tmp_name(m_ply->getName() + "_L");
		tmp_name.append(layer_number);
		m_polyline->setName(tmp_name);

		m_polyline->setDataType(1);
		// OK
		m_polyline->epsilon = _min_edge_length / 2.;
		for (i = 0; i < ply_nod_vector_layer; i++)
		{
			CGLPoint* point(new CGLPoint(nod_vector[ply_nod_vector[i] + j * nodes_per_layer]->getData()));
			m_polyline->point_vector.push_back(point);
		}
		// OK    m_polyline->SortPointVectorByDistance();
		m_polyline->SetPointOrderByDistance(m_ply->point_vector[0]);
		polyline_vector.push_back(m_polyline);
		m_polyline->WritePointVector(m_polyline->getName());
		m_polyline->WriteTecplot(" "); // OK41
	}
}

/**************************************************************************
   MSHLib-Method:
   08/2006 OK Implementation
   03/2010 TF change to new data structures, changed algorithm
**************************************************************************/
void CFEMesh::GetELEOnPLY(const GEOLIB::Polyline* ply, std::vector<size_t>& ele_vector_ply, bool With1DElements)
{
	// WW EEMath_Group::vec<CEdge*> ele_edges_vector(15);
	int loc_edge_nidx[3];
	Math_Group::vec<CNode*> edge_nodes(3);

	std::vector<size_t> nodes_near_ply;

	// get mesh nodes near the polyline
	GetNODOnPLY(ply, nodes_near_ply);

	// clear the given vector
	ele_vector_ply.clear();

	// loop over all elements
	for (size_t i = 0; i < ele_vector.size(); i++)
	{
		CElem* elem = ele_vector[i]; // WW
		// WW ele_vector[i]->GetEdges (ele_edges_vector);
		size_t n_edges = elem->GetEdgesNumber();
		// Add 1D Elements for models with mixed 1D/2D Elements				BG, 11/2011
		if ((elem->GetDimension() == 1) && (With1DElements == true))
		{
			for (size_t k = 0; k < nodes_near_ply.size(); k++)
			{
				if ((static_cast<size_t>(elem->GetNodeIndex(0)) == nodes_near_ply[k])
				    || (static_cast<size_t>(elem->GetNodeIndex(1)) == nodes_near_ply[k]))
					ele_vector_ply.push_back(elem->GetIndex());
			}
		}
		// loop over all edges of the i-th element
		for (size_t j = 0; j < n_edges; j++)
		{
			// WWele_edges_vector[j]->GetNodes(edge_nodes);
			elem->GetLocalIndicesOfEdgeNodes(j, loc_edge_nidx); // WW
			size_t selected(0);
			// get all elements having an edge in common with ply
			for (size_t k = 0; k < nodes_near_ply.size(); k++)
			{
				// if (edge_nodes[0]->GetIndex() == nodes_near_ply[k])
				if (elem->GetNodeIndex(loc_edge_nidx[0]) == static_cast<long>(nodes_near_ply[k]))
					selected++;
				// if (edge_nodes[1]->GetIndex() == nodes_near_ply[k])
				if (elem->GetNodeIndex(loc_edge_nidx[1]) == static_cast<long>(nodes_near_ply[k]))
					selected++;
			}
			if (selected == 2)
				ele_vector_ply.push_back(ele_vector[i]->GetIndex());
		}
	}
}
/**************************************************************************
   MSHLib-Method:
   Task: All nodes vertical to a polyline
   02/2009 OK
**************************************************************************/
void CFEMesh::GetNODOnSFC_PLY_Z(Surface* m_sfc, std::vector<long>& msh_nod_vector)
{
	std::vector<CGLPolyline*>::iterator p_ply;
	CGLPolyline* m_ply = NULL;
	// .................................................................
	// nodes close to first polyline
	p_ply = m_sfc->polyline_of_surface_vector.begin();
	while (p_ply != m_sfc->polyline_of_surface_vector.end())
	{
		m_ply = *p_ply;
		GetNODOnPLY_XY(m_ply, msh_nod_vector);
		break;
	}
}
/**************************************************************************
   GeoSys-Method:
   Task:
   Programing:
   07/2005 OK Implementation
**************************************************************************/
std::ios::pos_type CFEMesh::GMSReadTIN(std::ifstream* tin_file)
{
	std::string line_string, s_buff;
	std::stringstream in;
	std::ios::pos_type position;
	int i;
	i = 0;
	std::string line;
	std::string sub_line;
	long no_vertexes = 0;
	CNode* m_nod = NULL;
	long no_triangles = 0;
	CElem* m_ele = NULL;
	double xyz[3];
	//========================================================================
	while (sub_line.compare("ENDT"))
	{
		getline(*tin_file, s_buff); // WW
		in.str(s_buff);
		in >> sub_line;
		//................................................................
		// TNAM "PriTIN_1gr"
		if (sub_line.find("TNAM") != std::string::npos)
			in >> sub_line;
		//................................................................
		// VERT 3173
		if (sub_line.find("VERT") != std::string::npos)
		{
			in >> no_vertexes;
			in.clear();
			for (i = 0; i < no_vertexes; i++)
			{
				m_nod = new CNode(i);
				getline(*tin_file, s_buff); // WW
				in.str(s_buff);
				in >> xyz[0] >> xyz[1] >> xyz[2];
				m_nod->SetCoordinates(xyz);
				nod_vector.push_back(m_nod);
				in.clear();
			}
		}
		//................................................................
		if (sub_line.find("TRI") != std::string::npos) // TRI 6117
		{
			// Evaluate ele_type
			in >> no_triangles;
			in.clear();
			/* OKWW
			   for(i=0;i<no_triangles;i++){
			   m_ele = new FiniteElement::CElement();
			   m_ele->type_name = "tri";
			   m_ele->ElementType = ele_type;
			   m_ele->nnodes = 3;
			   m_ele->nodes = new long[3];
			   in.str(GetLineFromFile1(tin_file)); // 3169	3168	3173
			   in >> ele_nod_number;
			   m_ele->nodes_index[0] = ele_nod_number-1;
			   in >> ele_nod_number;
			   m_ele->nodes_index[1] = ele_nod_number-1;
			   in >> ele_nod_number;
			   m_ele->nodes_index[2] = ele_nod_number-1;
			   in.clear();
			   m_ele->nodes_xyz = new double[9];
			   for(k=0;k<m_ele->nnodes;k++){
			   m_ele->nodes_xyz[k]                 = nod_vector[m_ele->nodes_index[k]]->x;
			   m_ele->nodes_xyz[k+m_ele->nnodes]   = nod_vector[m_ele->nodes_index[k]]->y;
			   m_ele->nodes_xyz[k+2*m_ele->nnodes] = nod_vector[m_ele->nodes_index[k]]->z;
			   }
			   ele_vector.push_back(m_ele);
			   }
			 */
			for (i = 0; i < no_triangles; i++)
			{
				m_ele = new CElem(i);
				m_ele->geo_type = MshElemType::TRIANGLE;
				m_ele->Read(*tin_file, 3);
				ele_vector.push_back(m_ele);
			}
		}
		//................................................................
		in.clear();
	}
	return position;
}

#ifdef ObsoleteGUI // WW 03.2012
/**************************************************************************
   GeoSys-Method:
   Task:
   Programing:
   02/2005 OK Implementation (MMP groups)
   02/2005 OK Activate from vector
   08/2005 WW Changes due to geometry objects applied
   see also: BuildActiveElementsArray
**************************************************************************/
void CFEMesh::SetActiveElements(std::vector<long>& elements_active)
{
	const size_t ele_vector_size(this->ele_vector.size());
	for (size_t i = 0; i < ele_vector_size; i++)
		ele_vector[i]->MarkingAll(false);

	const size_t elements_active_size(elements_active.size());
	for (size_t i = 0; i < elements_active_size; i++)
		ele_vector[elements_active[i]]->MarkingAll(true);
}

/**************************************************************************
   MSHLib-Method:
   Task:
   const int NLayers        :      Number of layers (start mesh)
   const int Layer          :      Layer number of the layer to be refined
   const int SUBLayers      :      Number of sublayers
   Programing:
   10/2003 WW/MB Erste Version (PrismRefine)
   08/2004 MB NElementsPerLayer = msh_no_pris / NLayers;
   09/2005 OK MSH
 ***************************************************************************/
void CFEMesh::PrismRefine(int Layer, int subdivision)
{
	const size_t nn = 6;
	int j, nes;
	size_t* element_nodes = NULL;
	// WW double nx[6], ny[6],
	double nz[6];
	// WW double dx[3], dy[3],;
	double dz[3];
	double newz;
	int row(Layer);
	int NRowsToShift;
	int NRows(_n_msh_layer + 1);
	int CountNLayers;
	CNode* m_nod = NULL;
	const int NSubLayers = subdivision + 1; // OK
	const size_t NumElement0(ele_vector.size());
	// NodeListSize() / (NLayers+1);
	size_t NNodesPerRow(nod_vector.size() / (_n_msh_layer + 1));
	// msh_no_pris / NLayers;
	//      int NElementsPerLayer = ele_vector.size() / _n_msh_layer;

	size_t nod_vector_size_add(NNodesPerRow * subdivision);
	for (size_t i = 0; i < nod_vector_size_add; i++)
	{
		m_nod = new CNode(i);
		nod_vector.push_back(m_nod);
	}
	// nod_vector.resize(nod_vector_size_new);

	// Initialisierung der Knoten flags
	const size_t nod_vector_size(nod_vector.size());
	for (size_t i = 0; i < nod_vector_size; i++)
		nod_vector[i]->SetMark(true);

	CElem* m_ele = NULL;
	CElem* m_ele_new = NULL;
	for (size_t ne = 0; ne < NumElement0; ne++)
	{
		m_ele = ele_vector[ne];
		if (m_ele->GetElementType() == MshElemType::PRISM)
		{
			// element_nodes = m_ele->nodes;
			CountNLayers = _n_msh_layer;
			for (size_t i = 0; i < nn; i++)
				// WW nx[i] = nod_vector[m_ele->nodes_index[i]]->X();
				// WW  ny[i] = nod_vector[m_ele->nodes_index[i]]->Y();
				nz[i] = nod_vector[m_ele->nodes_index[i]]->getData()[2];
			nes = 0;
			for (size_t i = 0; i < 3; i++)
				if (element_nodes[i] >= (row - 1) * NNodesPerRow && element_nodes[i] <= (row * NNodesPerRow) - 1)
					nes++;

			if (nes == 3)
			{
				for (size_t i = 0; i < 3; i++)
					// WW dx[i] = (nx[i + 3] - nx[i]) / (float) NSubLayers;
					// WW dy[i] = (ny[i + 3] - ny[i]) / (float) NSubLayers;
					dz[i] = (nz[i + 3] - nz[i]) / (float)NSubLayers;
				// Create new nodes
				// Loop over SubLayers
				for (int iii = 0; iii < NSubLayers - 1; iii++)
				{
					// neue Knoten ganz unten
					for (size_t i = 0; i < 3; i++)
						// if(NODGetFreeSurfaceFlag(element_nodes[i])==0){
						if (nod_vector[element_nodes[i]]->GetMark())
						{
							// m_nod = new CMSHNodes(); //kno = (Knoten *)CreateNodeGeometry();
							m_nod = nod_vector[(m_ele->nodes_index[i] + ((CountNLayers + 2) - row) * NNodesPerRow)];

							// TF                       xyz[0] = nod_vector[m_ele->nodes_index[i]
							// TF                          + ((CountNLayers + 1) - row) * NNodesPerRow]->X();
							// TF                       xyz[1] = nod_vector[m_ele->nodes_index[i]
							// TF                          + ((CountNLayers + 1) - row)
							// TF                          * NNodesPerRow]->Y();
							// TF                       xyz[2] = nod_vector[m_ele->nodes_index[i]
							// TF                          + ((CountNLayers + 1) - row)
							//                           * NNodesPerRow]->Z();
							// PlaceNode(kno,(element_nodes[i] + ((CountNLayers+2) - row) * NNodesPerRow));
							// TF								m_nod->SetCoordinates(xyz);
							m_nod->SetCoordinates(nod_vector[m_ele->nodes_index[i]
							                                 + ((CountNLayers + 1) - row) * NNodesPerRow]->getData());
							nod_vector[(m_ele->nodes_index[i] + ((CountNLayers + 2) - row) * NNodesPerRow)] = m_nod;
						}
					// neues Element ganz unten
					m_ele_new = new CElem();
					// m_ele_new = m_ele;
					m_ele_new->nnodes = m_ele->nnodes;
					m_ele_new->patch_index = m_ele->GetPatchIndex();
					m_ele_new->SetElementType(m_ele->GetElementType());
					// KR memory leak! knoten = new long[6];
					for (j = 0; j < 3; j++)
					{
						m_ele_new->nodes_index[j] = element_nodes[j] + ((CountNLayers + 1) - row) * NNodesPerRow;
						m_ele_new->nodes_index[j + 3] = element_nodes[j] + ((CountNLayers + 2) - row) * NNodesPerRow;
					}
					// KR fixed memory leak, knoten not needed
					// for (j = 0; j < 6; j++) {
					//	m_ele_new->nodes_index[j] = knoten[j];
					//}
					ele_vector.push_back(m_ele_new);
					/* "rowx hochziehen"   */
					/* loop ?er die betroffenen rows   */
					NRowsToShift = NRows - Layer;
					for (int i = NRowsToShift; i > 0; i--)
					{
						if (i != 1)
						{
							for (j = 0; j < 3; j++)
								// if(NODGetFreeSurfaceFlag(element_nodes[j])==0){
								if (nod_vector[element_nodes[j]]->GetMark())
								{
									m_nod = nod_vector[element_nodes[j + 3] + NNodesPerRow * (iii + i - 1)];
									m_nod->SetZ(
									    nod_vector[element_nodes[j] + NNodesPerRow * (iii + i - 1)]->getData()[2]);
									// SetNodeZ(element_nodes[j+3] + NNodesPerRow*(iii+i-1), GetNodeZ(element_nodes[j] +
									// NNodesPerRow*(iii+i-1)));
								}
						}
						else
						{
							for (j = 0; j < 3; j++)
								// if(NODGetFreeSurfaceFlag(element_nodes[j])==0)   {
								if (nod_vector[element_nodes[j]]->GetMark())
								{
									newz = nod_vector[element_nodes[j]]->getData()[2] + (dz[j] * (iii + 1));
									// newz = GetNodeZ(element_nodes[j]) + (dz[j]*(iii+1));
									nod_vector[element_nodes[j] + (i)*NNodesPerRow * (iii + 1)]->SetZ(newz);
									// SetNodeZ(element_nodes[j] + (i) * NNodesPerRow *(iii+1), newz);
								}
						}
					} /* end for Rows to shift */
					if (iii == NSubLayers - 2)
						for (j = 0; j < 3; j++)
							// NODSetFreeSurfaceFlag(element_nodes[j], 33);
							nod_vector[element_nodes[j]]->SetMark(false);

					CountNLayers++;
				} /* End Loop over SubLayers  */
			} /* End if nes==3 */
		} /* Elementtyp ==6 */
	} /* Element loop */
	_n_msh_layer += subdivision;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB Implementation based on NodeExists.
   03/2010 TF changed data type of loop variables
**************************************************************************/
bool CFEMesh::NodeExists(size_t node)
{
	size_t no_nodes(nod_vector.size());

	for (size_t i = 0; i < no_nodes; i++)
		if (node == (size_t)Eqs2Global_NodeIndex[i])
			return true;
	return false;
}
#endif // end ObsoleteGUI

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   10/2005 OK Implementation
   02/2006 WW Ordering and remove bugs
**************************************************************************/
void CFEMesh::ConnectedNodes(bool quadratic) const
{
#define noTestConnectedNodes
	bool exist = false;

	for (size_t i = 0; i < nod_vector.size(); i++)
	{
		CNode* nod = nod_vector[i];
		size_t n_connected_elements(nod->getConnectedElementIDs().size());
		for (size_t j = 0; j < n_connected_elements; j++)
		{
			CElem* ele = ele_vector[nod->getConnectedElementIDs()[j]];
			size_t n_quadratic_node(static_cast<size_t>(ele->GetNodesNumber(quadratic)));
			for (size_t l = 0; l < n_quadratic_node; l++)
			{
				exist = false;
				size_t n_connected_nodes(nod->getConnectedNodes().size());
				// WW
				for (size_t k = 0; k < n_connected_nodes; k++)
					if (nod->getConnectedNodes()[k] == static_cast<size_t>(ele->nodes_index[l]))
					{
						exist = true;
						break;
					}
				if (!exist) // WW
					nod->getConnectedNodes().push_back(ele->nodes_index[l]);
			}
		}
	}

	// Sorting. WW
	// WW
	for (size_t i = 0; i < nod_vector.size(); i++)
	{
		CNode* nod = nod_vector[i];
		size_t n_connected_nodes = nod->getConnectedNodes().size();
		for (size_t k = 0; k < n_connected_nodes; k++)
		{
			for (size_t l = k; l < n_connected_nodes; l++)
				if (nod->getConnectedNodes()[l] < nod->getConnectedNodes()[k])
				{
					const size_t n(nod->getConnectedNodes()[k]);
					nod->getConnectedNodes()[k] = nod->getConnectedNodes()[l];
					nod->getConnectedNodes()[l] = n;
				}
		}
	}
//----------------------------------------------------------------------
#ifdef TestConnectedNodes
	for (i = 0; i < (long)nod_vector.size(); i++)
	{
		nod = nod_vector[i];
		cout << (int)nod->connected_nodes.size() << ": ";
		for (m = 0; m < (int)nod->connected_nodes.size(); m++)
			cout << nod->connected_nodes[m] << " ";
		cout << "\n";
	}
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:  Renumbering nodes corresponding to the activiate of elements
   Programing:
   01/2006 YD Implementation
**************************************************************************/
void CFEMesh::FaceNormal()
{
	int idx0_face, idx1_face, idx_owner, index0 = 0, index1 = 0;

	CElem* elem(NULL);
	CElem* elem_face(NULL);

	if (!face_normal.empty())
		return; // WW
	//------------------------
	for (size_t i = 0; i < face_vector.size(); i++)
	{
		elem_face = face_vector[i];
		elem = face_vector[i]->GetOwner();

		if (elem->GetElementType() == MshElemType::LINE)
			return;

		// WW      int no_face_vertex = face_vector[i]->GetVertexNumber();
		int no_owner_vertex = face_vector[i]->GetOwner()->GetVertexNumber();
		idx0_face = face_vector[i]->GetNodeIndex(0);
		idx1_face = face_vector[i]->GetNodeIndex(1);
		double* normal = new double[3];
		for (int j = 0; j < no_owner_vertex; j++)
		{
			idx_owner = face_vector[i]->GetOwner()->GetNodeIndex(j);
			if (idx0_face == idx_owner)
				index0 = j;
		}
		for (int j = 0; j < no_owner_vertex; j++)
		{
			idx_owner = face_vector[i]->GetOwner()->GetNodeIndex(j);
			if (idx1_face == idx_owner)
				index1 = j;
		}
		if (elem->GetMark())
		{
			if ((index1 - index0) >= 1)
				elem->FaceNormal(index0, index1, normal);
			else
				elem->FaceNormal(index1, index0, normal);
		}
		face_normal.push_back(normal);
		elem_face->ComputeVolume();
	}
}

/**************************************************************************
   MSHLib-Method:
   Programing:
   03/2006 OK Implementation
**************************************************************************/
void CFEMesh::SetNODPatchAreas()
{
	long e;
	int n1 = 0, n2 = 0;
	double v1[3], v2[3], v3[3];
	double patch_area;
	double x0, y0, z0;
	CNode* m_nod = NULL;
	CElem* m_ele = NULL;
	//----------------------------------------------------------------------
	size_t nNodes(nod_vector.size());
	for (size_t i = 0; i < nNodes; i++)
	{
		m_nod = nod_vector[i]; // this node
		patch_area = 0.0;
		//....................................................................
		// triangle neighbor nodes
		for (size_t j = 0; j < m_nod->getConnectedElementIDs().size(); j++)
		{
			e = m_nod->getConnectedElementIDs()[j];
			m_ele = ele_vector[e];
			for (size_t k = 0; k < 3; k++)
				if (m_ele->GetNodeIndex(k) == static_cast<int>(i))
				{
					switch (k)
					{
						case 0:
							n1 = 2;
							n2 = 1;
							break;
						case 1:
							n1 = 0;
							n2 = 2;
							break;
						case 2:
							n1 = 1;
							n2 = 0;
							break;
					}
				}
			//..................................................................
			double const* gravity_center(m_ele->GetGravityCenter());
			double const* const pnt(m_nod->getData());
			v2[0] = gravity_center[0] - pnt[0];
			v2[1] = gravity_center[1] - pnt[1];
			v2[2] = gravity_center[2] - pnt[2];
			//..................................................................
			//				m_nod1 = nod_vector[m_ele->GetNodeIndex(n1)];
			double const* const pnt1(nod_vector[m_ele->GetNodeIndex(n1)]->getData());
			x0 = 0.5 * (pnt1[0] - pnt[0]);
			y0 = 0.5 * (pnt1[1] - pnt[1]);
			z0 = 0.5 * (pnt1[2] - pnt[2]);
			v1[0] = x0 - pnt[0];
			v1[1] = y0 - pnt[1];
			v1[2] = z0 - pnt[2];
			CrossProduction(v1, v2, v3);
			patch_area += 0.5 * MBtrgVec(v3, 3);
			//..................................................................
			//				m_nod2 = nod_vector[m_ele->GetNodeIndex(n2)];
			double const* const pnt2(nod_vector[m_ele->GetNodeIndex(n2)]->getData());
			x0 = 0.5 * (pnt2[0] - pnt[0]);
			y0 = 0.5 * (pnt2[1] - pnt[1]);
			z0 = 0.5 * (pnt2[2] - pnt[2]);
			v1[0] = x0 - pnt[0];
			v1[1] = y0 - pnt[1];
			v1[2] = z0 - pnt[2];
			CrossProduction(v1, v2, v3);
			patch_area += 0.5 * MBtrgVec(v3, 3);
			//..................................................................
		}
		m_nod->patch_area = patch_area;
	}
}

/**************************************************************************
   MSHLib-Method:
   Programing:
   03/2006 OK Implementation
**************************************************************************/
void CFEMesh::SetNetworkIntersectionNodes()
{
	long i, e;
	int j, k;
	double v3[3], nr1[3], nr2[3];
	// WW  double* gravity_center;
	CNode* m_nod = NULL;
	CElem* m_ele = NULL;
	CElem* m_ele1 = NULL;
	//----------------------------------------------------------------------
	// Is node intersection node
	for (i = 0; i < (long)nod_vector.size(); i++)
		m_nod = nod_vector[i];
	// KR m_nod->selected = false;
	// KR double eps = 1e-3;
	for (i = 0; i < (long)nod_vector.size(); i++)
	{
		m_nod = nod_vector[i]; // this node
		if ((int)m_nod->getConnectedElementIDs().size() == 0)
			continue;
		m_ele = ele_vector[m_nod->getConnectedElementIDs()[0]];
		for (k = 0; k < 3; k++)
			nr1[k] = (*m_ele->transform_tensor)(2, k);
		//....................................................................
		// Compare element normal vectors
		for (j = 1; j < (int)m_nod->getConnectedElementIDs().size(); j++)
		{
			e = m_nod->getConnectedElementIDs()[j];
			m_ele1 = ele_vector[e];
			for (k = 0; k < 3; k++)
				nr2[k] = (*m_ele1->transform_tensor)(2, k);
			CrossProduction(nr1, nr2, v3);
			/* KR
			   if (MBtrgVec(v3, 3) > eps)
			   m_nod->selected = true;
			 */
		}
	}
	// TF useless code
	//      // Count non-intersection nodes
	//      long no_non_intersection_nodes = 0;
	//      for (i = 0; i < (long) nod_vector.size(); i++)
	//      {
	//         m_nod = nod_vector[i];
	//         if (m_nod->selected)
	//            continue;
	//         no_non_intersection_nodes++;
	//      }
}

#ifdef NEW_EQS // 1.11.2007 WW
/**************************************************************************
   MSHLib-Method:
   Programing:
   11/2007 WW Implementation
   04/2011 WW CRS storage
**************************************************************************/
void CFEMesh::CreateSparseTable()
{
	Math_Group::StorageType stype;
	stype = Math_Group::JDS;
	for (int i = 0; i < (int)num_vector.size(); i++)
		if (num_vector[i]->ls_storage_method == 100)
		{
			stype = Math_Group::CRS;
			break;
		}

	// Symmetry case is skipped.
	// 1. Sparse_graph_H for high order interpolation. Up to now, deformation
	if (NodesNumber_Linear != NodesNumber_Quadratic)
		sparse_graph_H = new SparseTable(this, true, false, stype);
	// 2. M coupled with other processes with linear element
	if (sparse_graph_H)
	{
		if ((int)pcs_vector.size() > 1)
			sparse_graph = new SparseTable(this, false, false, stype);
	}
	// 3. For process with linear elements
	else
		sparse_graph = new SparseTable(this, false, false, stype);

	//  sparse_graph->Write();
	//  sparse_graph_H->Write();
	//
	// ofstream Dum("sparse.txt", ios::out);
	// sparse_graph_H->Write(Dum);
}
#endif //#ifndef NEW_EQS  // 05.03.2010 WW


#ifndef NDEBUG
GEOLIB::Grid<MeshLib::CNode> const* CFEMesh::getGrid() const
{
	return _mesh_grid;
}
#endif

/*!
   brief Find the element by a point
   YS/WW 05/2012
*/
size_t CFEMesh::FindElementByPoint(const double* xyz)
{
	double x1[3], x2[3], x3[3], x4[3], x5[3], x6[3], x7[3], x8[3];
	double a, a_sub[12];
	CNode* a_node;
	CElem* n_ele;
	double tol = 1e-9;

	const size_t ne = ele_vector.size();
	for (size_t i = 0; i < ne; i++)
	{
		n_ele = ele_vector[i];
		a = n_ele->GetVolume();

		a_node = n_ele->GetNode(0);
		x1[0] = a_node->X();
		x1[1] = a_node->Y();
		x1[2] = a_node->Z();
		a_node = n_ele->GetNode(1);
		x2[0] = a_node->X();
		x2[1] = a_node->Y();
		x2[2] = a_node->Z();

		if (n_ele->GetElementType() != MshElemType::LINE)
		{
			a_node = n_ele->GetNode(2);
			x3[0] = a_node->X();
			x3[1] = a_node->Y();
			x3[2] = a_node->Z();
		}

		if (n_ele->GetElementType() == MshElemType::QUAD || n_ele->GetElementType() == MshElemType::TETRAHEDRON)
		{
			a_node = n_ele->GetNode(3);
			x4[0] = a_node->X();
			x4[1] = a_node->Y();
			x4[2] = a_node->Z();
		}

		if (n_ele->GetElementType() == MshElemType::PYRAMID)
		{
			a_node = n_ele->GetNode(3);
			x4[0] = a_node->X();
			x4[1] = a_node->Y();
			x4[2] = a_node->Z();
			a_node = n_ele->GetNode(4);
			x5[0] = a_node->X();
			x5[1] = a_node->Y();
			x5[2] = a_node->Z();
		}

		if (n_ele->GetElementType() == MshElemType::PRISM)
		{
			a_node = n_ele->GetNode(3);
			x4[0] = a_node->X();
			x4[1] = a_node->Y();
			x4[2] = a_node->Z();
			a_node = n_ele->GetNode(4);
			x5[0] = a_node->X();
			x5[1] = a_node->Y();
			x5[2] = a_node->Z();
			a_node = n_ele->GetNode(5);
			x6[0] = a_node->X();
			x6[1] = a_node->Y();
			x6[2] = a_node->Z();
		}

		if (n_ele->GetElementType() == MshElemType::HEXAHEDRON)
		{
			a_node = n_ele->GetNode(3);
			x4[0] = a_node->X();
			x4[1] = a_node->Y();
			x4[2] = a_node->Z();
			a_node = n_ele->GetNode(4);
			x5[0] = a_node->X();
			x5[1] = a_node->Y();
			x5[2] = a_node->Z();
			a_node = n_ele->GetNode(5);
			x6[0] = a_node->X();
			x6[1] = a_node->Y();
			x6[2] = a_node->Z();
			a_node = n_ele->GetNode(6);
			x7[0] = a_node->X();
			x7[1] = a_node->Y();
			x7[2] = a_node->Z();
			a_node = n_ele->GetNode(7);
			x8[0] = a_node->X();
			x8[1] = a_node->Y();
			x8[2] = a_node->Z();
		}

		switch (n_ele->GetElementType())
		{
			case MshElemType::LINE:
				double d1, d2, d3;
				d1 = 0.;
				d2 = 0.;
				d3 = 0.;
				for (int kk = 0; kk < 3; kk++)
				{
					d1 += (x1[kk] - xyz[kk]) * (x1[kk] - xyz[kk]);
					d2 += (x2[kk] - xyz[kk]) * (x2[kk] - xyz[kk]);
					d3 += (x1[kk] - x2[kk]) * (x1[kk] - x2[kk]);
				}
				d1 = sqrt(d1);
				d2 = sqrt(d2);
				d3 = sqrt(d3);

				if (fabs((d1 + d2 - d3) / d3) < tol)
				{
					return i;
				}
				break;

			case MshElemType::TRIANGLE:
				a_sub[0] = ComputeDetTri(x1, x2, xyz);
				a_sub[1] = ComputeDetTri(x2, x3, xyz);
				a_sub[2] = ComputeDetTri(x3, x1, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] - a) / a) < tol)
				{
					return i;
				}
				break;

			case MshElemType::QUAD:
				a_sub[0] = ComputeDetTri(x1, x2, xyz);
				a_sub[1] = ComputeDetTri(x2, x3, xyz);
				a_sub[2] = ComputeDetTri(x3, x4, xyz);
				a_sub[3] = ComputeDetTri(x4, x1, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] - a) / a) < tol)
				{
					return i;
				}
				break;

			case MshElemType::TETRAHEDRON:
				a_sub[0] = ComputeDetTex(x2, x4, x3, xyz);
				a_sub[1] = ComputeDetTex(x1, x3, x4, xyz);
				a_sub[2] = ComputeDetTex(x2, x1, x4, xyz);
				a_sub[3] = ComputeDetTex(x2, x3, x1, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] - a) / a) < tol)
				{
					return i;
				}
				break;

			case MshElemType::PYRAMID:
				a_sub[0] = ComputeDetTex(x1, x2, x4, xyz);
				a_sub[1] = ComputeDetTex(x2, x3, x4, xyz);
				a_sub[2] = ComputeDetTex(x1, x5, x2, xyz);
				a_sub[3] = ComputeDetTex(x2, x5, x4, xyz);
				a_sub[4] = ComputeDetTex(x3, x5, x4, xyz);
				a_sub[5] = ComputeDetTex(x4, x5, x1, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] - a) / a) < tol)
				{
					return i;
				}
				break;

			case MshElemType::PRISM:
				a_sub[0] = ComputeDetTex(x1, x2, x3, xyz);
				a_sub[1] = ComputeDetTex(x4, x6, x5, xyz);
				a_sub[2] = ComputeDetTex(x1, x4, x2, xyz);
				a_sub[3] = ComputeDetTex(x2, x4, x5, xyz);
				a_sub[4] = ComputeDetTex(x2, x5, x3, xyz);
				a_sub[5] = ComputeDetTex(x3, x5, x6, xyz);
				a_sub[6] = ComputeDetTex(x3, x6, x1, xyz);
				a_sub[7] = ComputeDetTex(x1, x6, x4, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] + a_sub[6] + a_sub[7] - a)
				         / a) < tol)
				{
					return i;
				}
				break;

			case MshElemType::HEXAHEDRON:
				a_sub[0] = ComputeDetTex(x1, x2, x4, xyz);
				a_sub[1] = ComputeDetTex(x4, x2, x3, xyz);
				a_sub[2] = ComputeDetTex(x2, x6, x3, xyz);
				a_sub[3] = ComputeDetTex(x3, x6, x7, xyz);
				a_sub[4] = ComputeDetTex(x3, x7, x4, xyz);
				a_sub[5] = ComputeDetTex(x4, x7, x8, xyz);
				a_sub[6] = ComputeDetTex(x4, x8, x1, xyz);
				a_sub[7] = ComputeDetTex(x1, x8, x5, xyz);
				a_sub[8] = ComputeDetTex(x1, x5, x2, xyz);
				a_sub[9] = ComputeDetTex(x2, x5, x6, xyz);
				a_sub[10] = ComputeDetTex(x5, x8, x6, xyz);
				a_sub[11] = ComputeDetTex(x6, x8, x7, xyz);

				if (fabs((a_sub[0] + a_sub[1] + a_sub[2] + a_sub[3] + a_sub[4] + a_sub[5] + a_sub[6] + a_sub[7]
				          + a_sub[8] + a_sub[9] + a_sub[10] + a_sub[11] - a) / a) < tol)
				{
					return i;
				}
				break;

			default:
				// do nothing with other elements
				break;
		}
	}
	// Not find
	return -1;
}

// 09. 2012 WW
/// Free the memory occupied by edges
void CFEMesh::FreeEdgeMemory()
{
	while (edge_vector.size())
	{
		delete edge_vector[edge_vector.size() - 1];
		edge_vector[edge_vector.size() - 1] = NULL;
		edge_vector.pop_back();
	}

	const size_t ne = ele_vector.size();
	for (size_t i = 0; i < ne; i++)
	{
		ele_vector[i]->FreeEdgeMemory();
	}
}

/**************************************************************************
MSHLib-Method:
Task: For TOTAL_FLUX calculation
Programing:
11/2014 JOD
last modification:
**************************************************************************/
void CFEMesh::GetConnectedElements(std::vector<long>& nodes_on_sfc, std::vector<long>& vec_elements)
{
	std::set<long> set_nodes_on_sfc;
	long i, j, k, l;

	//----------------------------------------------------------------------
	// -->PETSC

	// init
	for (i = 0; i < (long)ele_vector.size(); i++)
	{
		ele_vector[i]->selected = 0; // TODO can use a new variable
	}
	// ---- search elements ------------------------------------------------------------------
	for (i = 0; i < (long)nodes_on_sfc.size(); i++)
	{
		k = nodes_on_sfc[i];
		for (j = 0; j < (long)nod_vector[k]->getConnectedElementIDs().size(); j++)
		{
			l = nod_vector[k]->getConnectedElementIDs()[j];
			if (ele_vector[l]->selected == 0)
				vec_elements.push_back(l);
			ele_vector[l]->selected += 1; // number of elements at surface
		}
	}
}

/*!
   Find element nodes on the top surface of a mesh domain
   07.06.2010
   By WW
 */
void CFEMesh::markTopSurfaceFaceElements3D()
{
#ifdef output_top_z
	/// For output z coordinate of all nodes on the top surface
	/// 13.08.2010. WW
	vector<bool> node_mark(NodesNumber_Linear);
	for (std::size_t i = 0; i < NodesNumber_Linear; i++)
		node_mark[i] = false;
#endif

	const double tol = sqrt(std::numeric_limits<double>::epsilon()); // 1.e-5;
	for (std::size_t i = 0; i < face_vector.size(); i++)
	{
		CElem* elem = face_vector[i];
		CElem const* const own_elem = elem->GetOwner();

		//// In element
		//		// compute center of mass
		double cent[] = {0., 0., 0.};
		for (std::size_t k = 0; k < own_elem->GetNodesNumber(false); k++)
		{
			//				node = own_elem->nodes[k];
			double const* const pnt_k(own_elem->GetNode(k)->getData());
			cent[0] += pnt_k[0];
			cent[1] += pnt_k[1];
			cent[2] += pnt_k[2];
		}
		for (int k = 0; k < 3; k++)
			cent[k] /= (double)own_elem->GetNodesNumber(false);

		//			node = elem->nodes[0];
		double const* const pnt_0(elem->GetNode(0)->getData());
		cent[0] -= pnt_0[0];
		cent[1] -= pnt_0[1];
		cent[2] -= pnt_0[2];
		NormalizeVector(cent, 3);
		elem->ComputeVolume();
		elem->FillTransformMatrix();
		/// Compute the normal to this surface element
		Math_Group::Matrix const* const transform_tensor = elem->getTransformTensor();
		double fac = cent[0] * (*transform_tensor)(0, 2) + cent[1] * (*transform_tensor)(1, 2)
		                   + cent[2] * (*transform_tensor)(2, 2);
		if (fac > 0.0)
			fac = -1.0;
		else
			fac = 1.0;
		//////

		/// If n.z>0
		if ((*transform_tensor)(2, 2) * fac > tol)
		{
			elem->SetMark(true);
			for (int k = 0; k < 3; k++)
				(*transform_tensor)(k, 2) *= fac;

#ifdef output_top_z
			for (std::size_t k = 0; k < elem->GetNodesNumber(false); k++)
				node_mark[elem->GetNodeIndex(k)] = true;
#endif
		}
		else if ((*transform_tensor)(2, 2) * fac < -tol)
		{
			elem->SetMark(false);
			for (int k = 0; k < 3; k++)
				(*transform_tensor)(k, 2) *= fac;

#ifdef output_top_z
			for (k = 0; k < elem->GetNodesNumber(quad); k++)
				node_mark[elem->GetNodeIndex(k)] = bottom;
#endif
		}
		else
			elem->SetMark(false);
	}

#ifdef output_top_z
	string ccc = FileName + "_top_head.asc";
	ofstream ofile_asci(ccc.c_str(), ios::trunc);

	for (std::size_t i = 0; i < (long)nod_vector.size(); i++)
	{
		node = nod_vector[i];
		if (!node_mark[i])
			continue;
		ofile_asci << node->GetIndex() << " " << node->Z() << "\n";
	}
	ofile_asci << "#STOP"
	           << "\n";
#endif
}

}
// namespace MeshLib
