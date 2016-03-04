/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   Task: Sparse matrix and linear equation solver
      Design and program by WW
   Programing:
   10/2007 WW/
**************************************************************************/
#ifndef eqs_class_INC
#define eqs_class_INC

// NEW_EQS To be removed
#ifdef NEW_EQS                                    //1.11.2007 WW
#include <cmath>
#include <iostream>
#include <vector>
//

#ifdef LIS
#include "lis.h"
#endif
#include "matrix_class.h"
class CNumerics;
class CRFProcess;
namespace process
{class CRFProcessDeformation;
}
namespace FiniteElement
{class CFiniteElementStd;
 class CFiniteElementVec;
}
using  ::CRFProcess;
using  process::CRFProcessDeformation;
using  FiniteElement::CFiniteElementStd;
using  FiniteElement::CFiniteElementVec;
//
#if defined(USE_MPI)
class CPARDomain;
#endif
//
namespace Math_Group
{
using namespace std;

class SparseTable;
//
class Linear_EQS
{
public:
	Linear_EQS(const SparseTable &sparse_table,
	           const long dof, bool messg = true);
#if defined(USE_MPI)
	Linear_EQS(const long size);
#endif
	~Linear_EQS();
	// Configure numerics
	void ConfigNumerics( CNumerics* m_num, const long n = 0);
	// Preconditioner;
	void Precond(double* vec_s, double* vec_r);
	void TransPrecond(double* vec_s, double* vec_r);
	//#if defined(USE_MPI)
	void Precond_Jacobi(const double* vec_s, double* vec_r);
	//#endif
	//
#ifdef JFNK_H2M
	/// GMRES. 01.09.2010. WW
	void setPCS(::CRFProcess* the_pcs) {a_pcs = the_pcs; }
	void Init_Precond_Jacobi_JFNK();
#endif
	void ComputePreconditioner();
	void ComputePreconditioner_Jacobi();
	void ComputePreconditioner_ILU() {       }
	//
	// Solver
#if defined(USE_MPI)
	int Solver(double* xg, const long n);
	int CG(double* xg, const long n);
	int BiCG(double* xg, const long n);   //02.2010. WW
	int BiCGStab(double* xg, const long n);
	int CGS(double* xg, const long n);
	double GetCPUtime() const { return cpu_time;  }
#else
#if defined(LIS) || defined(MKL)                                  //NW
	int Solver(CNumerics* num = NULL);
#else
	int Solver();
#endif
	int CG();
	int BiCG();                           //02.2010. WW
	int BiCGStab();
	int Gauss() {return -1; }
	int QMRCGStab() {return -1; }
	int CGNR() {return -1; }
	int CGS();
	int Richardson() {return -1; }
	int JOR() {return -1; }
	int SOR() {return -1; }
	int AMG1R5() {return -1; }
	int UMF() {return -1; }
	int GMRES();
#endif
	//
	void Initialize();
	void Clean();
	//
	// Access to the members
	void SetDOF(const int dof_n)          // For different processes with different DOF of OPDE. _new. 02/2010. WW
	{
		A->SetDOF(dof_n);
	}
	void SetKnownX_i(const long i, const double x_i);
	double X(const long i) const {return x[i]; }
	double RHS(const long i) const {return b[i]; }
	double NormX();
	double NormRHS() { return bNorm; }
#if defined(USE_MPI)
	int DOF() { return A->Dof(); }
	long Size() { return A->Size(); }
	void SetDomain(CPARDomain* a_dom) {dom = a_dom; }
#endif
	// Write
	void Write(std::ostream &os = std::cout);
	void WriteRHS(std::ostream &os = std::cout);
	void WriteX(std::ostream &os = std::cout);
	void Write_BIN(ostream &os);
private:                                          // Dot not remove this!
	CSparseMatrix* A;
	double* b;
	double* x;
	double* prec_M;
	//
#ifdef LIS
	// lis solver interface starts here
	LIS_MATRIX AA;
	LIS_VECTOR bb,xx;
	LIS_SOLVER solver;
#endif
#if defined(USE_MPI)
	CPARDomain* dom;
	// WW
	double* border_buffer0;
	double* border_buffer1;
	double cpu_time;
	friend class ::CPARDomain;
	//
	double dot (const double* xx,  const double* yy, const long n);
	inline void MatrixMulitVec(double* xx,  double* yy);
	inline void TransMatrixMulitVec(double* xx,  double* yy);
#endif
	//
	std::string solver_name;
	std::string precond_name;
	// Buffer
	std::vector<double*> f_buffer;
	// Controls
	int precond_type;
	int solver_type;
	bool message;
	int iter, max_iter;
	double tol, bNorm, error;
	long size_global;
	long size_A;
	// Operators
	double dot (const double* xx,  const double* yy);
	inline double Norm(const double* xx)  { return sqrt(dot(xx, xx)); }
	inline bool CheckNormRHS(const double normb_new);
#ifdef JFNK_H2M
	/// 30.06.2010. WW
	CRFProcess* a_pcs;
#endif
	/// GMRES. 30.06.2010. WW
	/// GMRES H matrix
	mutable Matrix H;
	int m_gmres;                          /// number of columns of H matrix
	void Update(double* x, int k, Matrix &h, double* s);
	void Get_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
	void Set_Plane_Rotation(double &dx, double &dy, double &cs, double &sn);
	//
	void Message();
	// Friends
	friend class ::CRFProcess;
	friend class process::CRFProcessDeformation;
	friend class FiniteElement::CFiniteElementStd;
	friend class FiniteElement::CFiniteElementVec;
	//
};
}

using Math_Group::Linear_EQS;
extern std::vector<Math_Group::Linear_EQS*> EQS_Vector;
#endif                                            // NEW_EQS
#endif
