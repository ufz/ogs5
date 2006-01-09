/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   Class element declaration
   class for finite element.
   Designed and programmed by WW, 06/2004
 */

#ifndef fem_INC

#define fem_INC

// C++
//#include <vector>
//#include <string>


#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
#include "prototyp.h"
#else
// MSH
#include "par_ddc.h"                              //OK //Moved from fem_ele_std.h. WW
#endif



#include "MSHEnums.h"

namespace Math_Group
{ class SymMatrix;
  class Matrix;
  typedef Matrix Vec;
}
namespace MeshLib
{ class CElem;
  class CNode;
  class CEdge;
}
class CRFProcess;

namespace  process { class CRFProcessDeformation;}

namespace FiniteElement
{
struct ExtrapolationMethod
{
	enum type
	{
		EXTRAPO_LINEAR,
		EXTRAPO_NEAREST,
		EXTRAPO_AVERAGE
	};
};

using Math_Group::SymMatrix;
using Math_Group::Matrix;
using Math_Group::Vec;

using MeshLib::CNode;
using MeshLib::CEdge;
using MeshLib::CElem;

class CElement
{
public:
	CElement (int CoordFlag, const int order = 1);
	virtual ~CElement ();
	//
	void ConfigElement(CElem* MElement, const int nquadrature_points, bool FaceIntegration = false);
	void ConfigFaceElement(CElem* MElement, bool FaceIntegration = false); // JOD 2014-11-10
	void setOrder(const int order);
	// Set Gauss point
	void SetGaussPoint(const int gp, int& gp_r, int& gp_s, int& gp_t);
	// Get Gauss integration information
	double GetGaussData(int gp, int& gp_r, int& gp_s, int& gp_t);

	// Compute values of shape function at integral point unit
	void ComputeShapefct(const int order);
	// Compute the Jacobian matrix. Return its determinate
	double computeJacobian(const int order);

	// Compute values of the derivatives of shape function at integral point
	void ComputeGradShapefct(const int order);
	// Compute the real coordinates from known unit coordinates
	void RealCoordinates(double* realXYZ);
	// Compute the unit coordinates from known unit coordinates
	void UnitCoordinates(double* realXYZ);
	// For axisymmetrical problems
	void CalculateRadius();
	//
	void setUnitCoordinates(double* u)
	{ for(int i = 0; i < 3; i++) unit[i] = u[i]; }

	// Finite element matrices and vectors
	// Compute the local finite element matrices
	void LocalAssembly(const long, const int) {}

	// Get values;
	int GetNumGaussPoints() const {return nGaussPoints; }
	int GetNumGaussSamples() const {return nGauss; }
	int Dim() const {return ele_dim; }
    double Getdshapefct(int in) {return dshapefct[in];}

	// Integrate Neumman type BC
	void FaceIntegration(double* NodeVal);

	void FaceNormalFluxIntegration(long index, double *NodeVal_adv, double *NodeVal, int* nodesFace, CElem* face, CRFProcess* m_pcs, double* normal_vector); // JOD 2014-11-10

	// Coupling
	//
	bool isTemperatureCoupling() const {return T_Flag; }
	bool isFluidPressureCoupling() const {return F_Flag; }
	int isDeformationCoupling() const {return D_Flag; }
	int isConcentrationCoupling() const {return C_Flag; }

	// Interpolate Gauss values
	double interpolate (double const * const nodalVal, const int order = 1) const;
	double interpolate (const int idx,  CRFProcess* m_pcs, const int order = 1);
	//double elemnt_average (const int idx, const int order =1);
	double elemnt_average (const int idx,  CRFProcess* m_pcs, const int order = 1);

	void SetCenterGP();
	int GetGPindex() const {return gp; }
	int GetElementIndex() const {return Index; }
	CElem* GetMeshElement() const         //OK
	{
		return MeshElement;
	}
	// For extropolating gauss value to node
	int GetLocalIndex(const int gp_r, const int gp_s, int gp_t);

	// DDC 05/2006
	void SetElementNodesDomain(long* ele_nodes)
	{element_nodes_dom = ele_nodes; }

	void SetRWPT(const int idx)           // PCH
	{
		PT_Flag = idx;
	}
protected:
	CElem* MeshElement;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	CPARDomain* m_dom;                    //OK
#endif
	long* element_nodes_dom;              //Only a pointer. For domain decomposition. WW

	friend class ::CRFProcess;
	friend class process::CRFProcessDeformation;

	// Coordinate indicator
	// 10:  X component only
	// 11: Y component only
	// 12: Z component only
	// 20:  X, Y component
	// 22:  X, Z component
	// 32:  X, Y, Z component
	int coordinate_system;
	bool axisymmetry;
	// Order of shape functions
	// Displacement, 2. Others, 1. Default, 1
	int Order;
	size_t ele_dim;                          // Dimension of element
	size_t dim;                              // Dimension of real dimension
	int nGaussPoints;                     // Number of Gauss points
	int nGauss;                           // Number of sample points for Gauss integration
	int gp;                               // Gauss point index.
	mutable double unit[4];               // Local coordintes
	double* Jacobian;                     // Jacobian matrix
	double* invJacobian;                  // Inverse of Jacobian matrix.
	double* shapefct;                     // Results of linear shape function at Gauss points
	double* shapefctHQ;                   // Results of quadratic shape function at Gauss points
	// Results of derivatives of linear shape function at Gauss points
	double* dshapefct;
	// Results of derivatives of quadratic shape function at Gauss points
	double* dshapefctHQ;
	//
	double x1buff[3],x2buff[3],x3buff[3],x4buff[3];
	// Pointer to the linear interpolation function
	VoidFuncDXCDX ShapeFunction;
	// Pointer to the quadratic interpolation function
	VoidFuncDXCDX ShapeFunctionHQ;
	// Pointer to the gradient of linear interpolation function
	VoidFuncDXCDX GradShapeFunction;
	// Pointer to the gradient of Quadratic interpolation function
	VoidFuncDXCDX GradShapeFunctionHQ;
	// Coupling
	int NodeShift[5];
	// Displacement column indeces in the node value table
	int Idx_dm0[3];
	int Idx_dm1[3];

	int idx_c0, idx_c1;

	// Coupling flag
	bool T_Flag;                          // Temperature
	bool C_Flag;                          // Concentration
	bool F_Flag;                          // Fluid
	int D_Flag;                           // Deformation
	int PT_Flag;                          // Particle Tracking Random Walk
	bool RD_Flag;                         // Dual Richards
	bool MCF_Flag;                         // Pressure Temperature coupled process
	// For extropolation
	double Xi_p;
	void SetExtropoGaussPoints(const int i); // 25.2.2007 WW
	double CalcAverageGaussPointValues(double* GpValues);
	double CalcXi_p();

	// Buffer
	int Index;
	int nNodes;
	int nnodes;
	int nnodesHQ;
	double time_unit_factor;
	double Radius;                        // For axisymmetrical problems
	long nodes[20];
	long eqs_number[20];
	double dShapefct[27];                 // Auxullary
	double X[20];
	double Y[20];
	double Z[20];
	double node_val[20];
	double dbuff[20];

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
  int act_nodes; //> activated nodes
  int act_nodes_h; //> activated nodes for high order elements
  int *idxm;  //> global indices of local matrix rows
  int *idxn;  //> global indices of local matrix columns
  int *local_idx; //> local index for local assemble
  //double *local_matrix; //>  local matrix
  //double *local_vec; //>  local vector
#endif
	ExtrapolationMethod::type extrapo_method;
	ExtrapolationMethod::type GetExtrapoMethod() {return extrapo_method; }
private:
	void ConfigNumerics(MshElemType::type elem_type, const int nquadrature_points);
};

/*------------------------------------------------------------------
   Element matrices:
   All local matrices are stored for the purpose of reducing
   compatation  time when steady state of the problem is reached
   12.01.2005. WW
   ------------------------------------------------------------------*/
class ElementMatrix
{
public:
	ElementMatrix() : Mass(NULL), Laplace(NULL),
		          Advection(NULL), Storage(NULL), Content(NULL),
		          CouplingA(NULL), CouplingB(NULL), Stiffness(NULL),
		          RHS(NULL)       //SB4200
	{
	}
	~ElementMatrix();
	// Allocate memory for strain coupling matrix
	void AllocateMemory(CElem* ele, int type = 0);
	// Set members
	void SetMass(Matrix* mass) { Mass = mass; }
	void SetMass_notsym(Matrix* mass) { Mass_notsym = mass; }
	void SetLaplace(Matrix* laplace) { Laplace = laplace; }
	void SetStiffness(Matrix* x) { Stiffness = x; }
	void SetAdvection(Matrix* x) { Advection = x; }
	void SetStorage(Matrix* x) { Storage = x; }
	void SetContent(Matrix* x) { Content = x; }
	void SetCouplingMatrixA(Matrix* cplM) {CouplingA = cplM; }
	void SetCouplingMatrixB(Matrix* cplM) {CouplingB = cplM; }
	void SetRHS(Vec* rhs) {RHS = rhs; }
	// Get members
	Matrix* GetMass() {return Mass; }
	Matrix* GetMass_notsym() {return Mass_notsym; }
	Matrix* GetLaplace() {return Laplace; }
	Matrix* GetStiffness() {return Stiffness; }
	Matrix* GetAdvection()                //SB4200
	{
		return Advection;
	}
	Matrix* GetStorage()                  //SB4200
	{
		return Storage;
	}
	Matrix* GetContent()                  //SB4200
	{
		return Content;
	}
	Matrix* GetCouplingMatrixA() {return CouplingA; }
	Matrix* GetCouplingMatrixB() {return CouplingB; }
	Vec* GetRHS() {return RHS; }
private:
	//TODO in more gernal way for the case of sym and unsym. WW      SymMatrix *Mass;
	//      SymMatrix *Laplace;
	Matrix* Mass;
	Matrix* Mass_notsym;
	Matrix* Laplace;
	Matrix* Advection;
	Matrix* Storage;
	Matrix* Content;
	Matrix* CouplingA;                    // Pressure coupling for M_Process
	Matrix* CouplingB;                    // Strain coupling gor H_Process
	Matrix* Stiffness;
	Vec* RHS;
};
}                                                 // end namespace

//=============================================
// For up coupling caculation in cel_*.cpp
// Will be removed when new FEM is ready
extern FiniteElement::CElement* elem_dm;
//=============================================
#endif
