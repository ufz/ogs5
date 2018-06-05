/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   PARLib - Object: PAR
   Task: class implementation
   Programing:
   07/2004 OK Implementation
   last modified:
**************************************************************************/
#ifndef par_dd_INC

#define par_dd_INC

// C++ STL
#include <fstream>
#include <list>
#include <string>
#include <vector>

#include "rf_pcs.h"

//---------   Declaration ---------------------WW
namespace FiniteElement
{
class CFiniteElementStd;
class CFiniteElementVec;
}
namespace process
{
class CRFProcessDeformation;
}
namespace Math_Group
{
class SparseTable;
class Linear_EQS;
}
using process::CRFProcessDeformation;
using FiniteElement::CFiniteElementVec;
using Math_Group::Linear_EQS;
using Math_Group::SparseTable;
//---------   Declaration ---------------------WW

#if defined(USE_MPI)
// WW
namespace MeshLib
{
class CFEMesh;
}
using MeshLib::CFEMesh;
#endif
void FindNodesOnInterface(CFEMesh* m_msh, bool quadr);

//-----------------------------------------------
class CPARDomain
{
private:
	std::vector<long*> element_nodes_dom; // Local DOM element nodes. WW
	long nnodes_dom;
	long nnodesHQ_dom;
	//#ifdef USE_MPI //WW
	long num_inner_nodes;
	long num_inner_nodesHQ;
	long num_boundary_nodes;
	long num_boundary_nodesHQ;
	friend void FindNodesOnInterface(CFEMesh* m_msh, bool quadr);
//#endif
#if defined(USE_MPI) // 13.12.2007 WW
	// Store global indices of all border nodes to border_nodes of the whole mesh
	// 0-->border_nodes_size, nodes for linear interpolation
	// border_nodes_size-->border_nodes_sizeH, nodes for quadratic interpolation
	std::vector<int> bnode_connected_dom; // Connected doms of border nodes
	long* t_border_nodes;
	long t_border_nodes_size;
	long t_border_nodes_sizeH;
	// For local EQS
	// For index mapping from local to global
	int dof, nq;
	long n_loc, n_bc; //, max_dimen;
	long i_start[2], i_end[2];
	long b_start[2], b_end[2], n_shift[2];
	// Data for concatenate internal entries of domain vectors
	int* receive_cnt_i;
	int* receive_disp_i;
#if defined(NEW_BREDUCE)
	// Data for concatenate border entries of domain vectors
	int* receive_cnt_b;
	int* receive_disp_b;
#endif
	//
	int* receive_cnt;
	int* receive_disp;
// friend class Math_Group::Linear_EQS;
//
#endif
	//
	long shift[5]; // WW
// Equation
#ifdef NEW_EQS // WW
	SparseTable* sparse_graph;
	SparseTable* sparse_graph_H;
	Linear_EQS* eqs; // WW
	Linear_EQS* eqsH; // WW
#endif
	friend class CRFProcess; // WW
	// WW //:: for SXC compiler
	friend class FiniteElement::CFiniteElementStd;
	// WW //:: for SXC compiler
	friend class FiniteElement::CFiniteElementVec;
	friend class process::CRFProcessDeformation; // WW //:: for SXC compiler
public:
	int ID;
	std::vector<long> elements;
	std::vector<long> nodes_inner;
	std::vector<long> nodes_halo;
	std::vector<long> nodes;
//?vector<double>matrix;
// EQS
#ifndef NEW_EQS
	LINEAR_SOLVER* eqs;
	LINEAR_SOLVER_PROPERTIES* lsp;
	char* lsp_name;
#endif
	// MSH
	CFEMesh* m_msh;
	// public:
	CPARDomain(void);
	~CPARDomain(void);
	std::ios::pos_type Read(std::ifstream*);
	void CreateNodes();
	// const long *longbuff, const bool quadr. WW
	void CreateElements(const bool quadr);
	void NodeConnectedNodes(); // WW
//
#ifdef NEW_EQS // WW
	void CreateSparseTable(); // WW
	void CreateEQS(); // WW
	void InitialEQS(CRFProcess* m_pcs); // WW
#else
	void CreateEQS(CRFProcess* m_pcs);
#endif
	void CalcElementMatrices(CRFProcess*);
	// WW   void AssembleMatrix(CRFProcess*);
	long GetDOMNode(long);
	int m_color[3]; // OK
	void WriteTecplot(std::string); // OK

	bool selected; // OK
	bool quadratic; // WW
	std::vector<long*> node_conneted_nodes; // WW
	std::vector<int> num_nodes2_node; // WW
	//
	long GetDomainNodes() const // WW
	{
		if (quadratic)
			return nnodesHQ_dom;
		else
			return nnodes_dom;
	}
	long GetDomainNodes(bool quad) const // WW
	{
		if (quad)
			return nnodesHQ_dom;
		else
			return nnodes_dom;
	}
#if defined(USE_MPI) // WW
	// long MaxDim() const {return max_dimen;}   //WW
	void ReleaseMemory();
	// WW
	void FillBorderNodeConnectDom(std::vector<int> allnodes_doms);
	long BSize() const // WW
	{
		return n_bc;
	}
	// WW
	void ConfigEQS(CNumerics* m_num, const long n, bool quad = false);
	// WW
	double Dot_Interior(const double* localr0, const double* localr1 = NULL);
	// WW
	void Global2Local(const double* global_x, double* local_x, const long n);
	// WW
	void Local2Global(const double* local_x, double* global_x, const long n);
	// WW
	void Global2Border(const double* x, double* local_x, const long n);
	// WW
	void Border2Global(const double* local_x, double* x, const long n);
	// WW
	void Local2Border(const double* local_x, double* border_x);
	// WW
	void Border2Local(const double* border_x, double* local_x);
	//
	double Dot_Border_Vec(const double* vec_x, const double* vec_y);
	//
	// WW
	void CatInnerX(double* global_x, const double* local_x, const long n);
	// WW
	void PrintEQS_CPUtime(std::ostream& os = std::cout);

#if defined(NEW_BREDUCE)
	void ReduceBorderV(double* local_x);
#endif
//
//
#endif
};

extern std::vector<CPARDomain*> dom_vector;

extern std::vector<int> node_connected_doms; // WW
extern void CountDoms2Nodes(CRFProcess* m_pcs); // WW
extern void DOMRead(std::string);
extern void DOMCreate();
//---- MPI Parallel --------------
// MH//HS
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) \
    || defined(USE_MPI_BRNS) || defined(USE_MPI_KRC)
extern char t_fname[3];
extern double time_ele_paral;

//#include <mpi.h>
// extern MPI_Comm comm_DDC;

#endif
//---- MPI Parallel --------------
#define DDC_FILE_EXTENSION ".ddc"
#endif
