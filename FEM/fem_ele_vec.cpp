/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   The members of class Element definitions.
   Designed and programmed by WW, 06/2004
 */

#include "fem_ele_vec.h"

#include <cfloat>
#include <algorithm>
#include "mathlib.h"
#include "matrix_class.h"
#include "matrix_routines.h"
#include "pcs_dm.h"
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
#include "tools.h" //12.2009. WW
// Equation
#if defined(NEW_EQS)
#include "equation_class.h"
using Math_Group::CSparseMatrix;
#endif

#ifndef USE_PETSC
#include "par_ddc.h"
#endif

//
#define COMP_MOL_MASS_AIR 28.96 // kg/kmol WW  28.96
#define GAS_CONSTANT 8314.41 // J/(kmol*K) WW

#if defined(USE_PETSC) // || defined(other parallel libs)//07.3013. WW
#include "PETSC/PETScLinearSolver.h"
#endif

std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;

namespace FiniteElement
{
using SolidProp::CSolidProperties;
using Math_Group::Matrix;
using Math_Group::Vec;
using ::CRFProcess;
using ::CMediumProperties;
using process::CRFProcessDeformation;
using MeshLib::CElem;
// Maximum number of nodes of linear elements
const int max_nnodes_LE = 8;
// Maximum number of nodes of 2D quadratic elements
const int max_nnodes_QE_2D = 9;
// Maximum number of nodes of 3D quadratic elements
const int max_nnodes_QE_3D = 20;

CFiniteElementVec::CFiniteElementVec(process::CRFProcessDeformation* dm_pcs,
	const int C_Sys_Flad, const int order)
    : CElement(C_Sys_Flad, order), pcs(dm_pcs), h_pcs(NULL), t_pcs(NULL), excavation(false),
	  ns((dim == 3)? 6: 4),  Flow_Type(-1),
      idx_P(-1), idx_P0(-1), idx_P1(-1), idx_P1_0(-1), idx_P2(-1),
	  idx_T0(-1), idx_T1(-1), idx_S0(-1), idx_S(-1), idx_Snw(-1), idx_pls(-1),
	  idx_p1_ini(-1), idx_p2_ini(-1),
	  PressureC(NULL), PressureC_S(NULL), PressureC_S_dp(NULL), b_rhs(NULL),
	  smat(NULL), m_mfp(NULL), m_mmp(NULL),
	  Temp(NULL), T1(NULL), Tem(273.15 + 23.0), S_Water(1.), eleV_DM(NULL),
	  _nodal_p1(NULL), _nodal_p2(NULL), _nodal_cp0(NULL), _nodal_dcp(NULL),
	  _nodal_S0(NULL), _nodal_S(NULL), AuxNodal1(NULL),
	  dynamic((dm_pcs->pcs_type_name_vector[0].find("DYNAMIC") != std::string::npos)
			? true : false),
	  Mass(NULL), Idx_Vel(NULL), dAcceleration(NULL),
	  beta2(dynamic ? dm_pcs->m_num->GetDynamicDamping_beta2() : 1.),
	  bbeta1(dynamic ? dm_pcs->m_num->GetDynamicDamping_bbeta(): 1.)
{
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	m_dom = NULL;
#endif

	for (int i = 0; i < 4; i++)
		NodeShift[i] = pcs->Shift[i];
	if (dynamic)
	{
		Mass = new Matrix(max_nnodes_QE_3D, max_nnodes_QE_3D);
		Idx_Vel = new int[3];
		dAcceleration = new Vec(max_nnodes_QE_3D * dim);

		Idx_dm0[0] = pcs->GetNodeValueIndex("ACCELERATION_X1");
		Idx_dm0[1] = pcs->GetNodeValueIndex("ACCELERATION_Y1");
		Idx_dm1[0] = Idx_dm0[0] + 1;
		Idx_dm1[1] = Idx_dm0[1] + 1;
		Idx_Vel[0] = pcs->GetNodeValueIndex("VELOCITY_DM_X");
		Idx_Vel[1] = pcs->GetNodeValueIndex("VELOCITY_DM_Y");
		//     if(problem_dimension_dm==3)
		if (dim == 3)
		{
			Idx_dm0[2] = pcs->GetNodeValueIndex("ACCELERATION_Z1");
			Idx_dm1[2] = Idx_dm1[2] + 1;
			Idx_Vel[2] = pcs->GetNodeValueIndex("VELOCITY_DM_Z");
		}
	}
	else
	{
		// Indecex in nodal value table
		Idx_dm0[0] = pcs->GetNodeValueIndex("DISPLACEMENT_X1");
		Idx_dm0[1] = pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
		Idx_dm1[0] = Idx_dm0[0] + 1;
		Idx_dm1[1] = Idx_dm0[1] + 1;

		//     if(problem_dimension_dm==3)
		if (dim == 3)
		{
			Idx_dm0[2] = pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
			Idx_dm1[2] = Idx_dm0[2] + 1;
		}
	}

	// Strain
	_nodal_strain_indices = new int[ns];
	_nodal_strain_indices[0] = pcs->GetNodeValueIndex("STRAIN_XX");
	_nodal_strain_indices[1] = pcs->GetNodeValueIndex("STRAIN_YY");
	_nodal_strain_indices[2] = pcs->GetNodeValueIndex("STRAIN_ZZ");
	_nodal_strain_indices[3] = pcs->GetNodeValueIndex("STRAIN_XY");
	// Stress
	_nodal_stress_indices = new int[ns];
	_nodal_stress_indices[0] = pcs->GetNodeValueIndex("STRESS_XX");
	_nodal_stress_indices[1] = pcs->GetNodeValueIndex("STRESS_YY");
	_nodal_stress_indices[2] = pcs->GetNodeValueIndex("STRESS_ZZ");
	_nodal_stress_indices[3] = pcs->GetNodeValueIndex("STRESS_XY");
	//
	if (dim == 3)
	{
		_nodal_strain_indices[4] = pcs->GetNodeValueIndex("STRAIN_XZ");
		_nodal_strain_indices[5] = pcs->GetNodeValueIndex("STRAIN_YZ");
		//
		_nodal_stress_indices[4] = pcs->GetNodeValueIndex("STRESS_XZ");
		_nodal_stress_indices[5] = pcs->GetNodeValueIndex("STRESS_YZ");
	}

	dstress = new double[ns];
	dstrain = new double[ns];
	stress_ne = new double[ns];
	strain_ne = new double[ns];
	stress0 = new double[ns];

	const int max_nodes =  (dim == 3) ? max_nnodes_QE_3D : max_nnodes_QE_2D;
	Sxx = new double[max_nodes];
	Syy = new double[max_nodes];
	Szz = new double[max_nodes];
	Sxy = new double[max_nodes];
	Sxz = (dim == 3) ? new double[max_nnodes_QE_3D]: NULL;
	Syz = (dim == 3) ? new double[max_nnodes_QE_3D]: NULL;
	pstr = new double[max_nodes];

	Disp =  new double[max_nodes * dim];

	B_matrix = new Matrix(ns, dim);
	B_matrix_T = new Matrix(dim, ns);
	De = new Matrix(ns, ns);
	ConsistDep = new Matrix(ns, ns);
	AuxMatrix = new Matrix(dim, dim);
	AuxMatrix2 = new Matrix(dim, ns);
	Stiffness = (pcs->Memory_Type == 0) ?
		new Matrix(max_nnodes_QE_3D * dim, max_nnodes_QE_3D * dim) : NULL;

	// For cache NW
	vec_B_matrix.resize(20);
	vec_B_matrix_T.resize(20);
	for (int i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		switch (dim)
		{
			case 2:
				vec_B_matrix[i] = new Matrix(4, 2);
				vec_B_matrix_T[i] = new Matrix(2, 4);
				break;
			case 3:
				vec_B_matrix[i] = new Matrix(6, 3);
				vec_B_matrix_T[i] = new Matrix(3, 6);
				break;
		}
	}

	RHS = (pcs->Memory_Type == 0) ? new Vec(3 * max_nnodes_QE_3D) : NULL;

	*B_matrix = 0.0;
	*B_matrix_T = 0.0;

	if (pcs->Memory_Type == 0) // Do not store local matrices
	{
		if (H_Process)
		{
			const int max_nnodes = (dim == 3) ? max_nnodes_QE_3D : max_nnodes_QE_2D;
			PressureC = new Matrix(3 * max_nnodes, max_nnodes);
		}
	}

	if (H_Process)
	{
		AuxNodal1 = (dim == 3) ? new double[max_nnodes_QE_3D * dim]
					           : new double[max_nnodes_QE_2D * dim];
	}
	// Coupling
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		//      if (pcs_vector[i]->pcs_type_name.find("FLOW") != string::npos) {
		// TF
		if (isFlowProcess(pcs_vector[i]->getProcessType()))
		{
			h_pcs = pcs_vector[i];
			// 25.04.2008, 04.09.2008  WW
			if (h_pcs->type == 1 || h_pcs->type == 41)
			{
				//				if (h_pcs->pcs_type_name.find("GROUND") != string::npos)
				// TF
				if (h_pcs->getProcessType() == GROUNDWATER_FLOW)
					Flow_Type = 10;
				else
					Flow_Type = 0;
				// 25.08.2005.  WW
			}
			else if (h_pcs->type == 14 || h_pcs->type == 22)
				Flow_Type = 1;
			else if (h_pcs->type == 1212 || h_pcs->type == 42)
			{
				Flow_Type = 2; // 25.04.2008.  WW

				// 07.2011. WW
				PressureC_S = new Matrix(60, 20);
				if (pcs->m_num->nls_method == 1 && h_pcs->type == 42) // Newton-raphson. WW
					PressureC_S_dp = new Matrix(60, 20);
			}
			// WW idx_P0 = pcs->GetNodeValueIndex("POROPRESSURE0");
			break;
			//		} else if (pcs_vector[i]->pcs_type_name.find("PS_GLOBAL") != string::npos) {
		} // TF
		else if (pcs_vector[i]->getProcessType() == PS_GLOBAL)
		{
			h_pcs = pcs_vector[i];
			if (h_pcs->type == 1313)
				Flow_Type = 3; // 05.05.2009.  PCH
			break;
		}
	}
	if (Flow_Type == 0)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		idx_p1_ini = h_pcs->GetNodeValueIndex("PRESSURE1_Ini");
		_nodal_p1 = new double[max_nnodes_LE];
		if (dynamic)
		{
			idx_P = h_pcs->GetNodeValueIndex("PRESSURE1");
			idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE_RATE1");
		}
	}
	if (Flow_Type == 10)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("HEAD") + 1;
		_nodal_p1 = new double[max_nnodes_LE];
	}
	else if (Flow_Type == 1)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		_nodal_p1 = new double[max_nnodes_LE];

		idx_P1_0 = h_pcs->GetNodeValueIndex("PRESSURE1");
		idx_p1_ini = h_pcs->GetNodeValueIndex("PRESSURE1_Ini");

		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		_nodal_S0 = new double[max_nnodes_LE];

		idx_S = h_pcs->GetNodeValueIndex("SATURATION1") + 1;
		_nodal_S = new double[max_nnodes_LE];

		_nodal_cp0 = new double[max_nnodes_LE];
		_nodal_dcp = new double[max_nnodes_LE];
	}
	else if (Flow_Type == 2)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		_nodal_p1 = new double[max_nnodes_LE];
		idx_p1_ini = h_pcs->GetNodeValueIndex("PRESSURE1_Ini");

		idx_P2 = h_pcs->GetNodeValueIndex("PRESSURE2") + 1;
		_nodal_p2 = new double[max_nnodes_LE];
		idx_p2_ini = h_pcs->GetNodeValueIndex("PRESSURE2_Ini");

		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		_nodal_S0 = new double[max_nnodes_LE];

		idx_S = h_pcs->GetNodeValueIndex("SATURATION1") + 1;
		_nodal_S = new double[max_nnodes_LE];

		_nodal_cp0 = new double[max_nnodes_LE];
		_nodal_dcp = new double[max_nnodes_LE];
	}
	else if (Flow_Type == 3)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		_nodal_p1 = new double[max_nnodes_LE];
		idx_P2 = h_pcs->GetNodeValueIndex("PRESSURE2");
		_nodal_p2 = new double[max_nnodes_LE];

		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		_nodal_S0 = new double[max_nnodes_LE];
		idx_S = idx_S0;
		_nodal_S = new double[max_nnodes_LE];

		idx_Snw = h_pcs->GetNodeValueIndex("SATURATION2") + 1;
	}

	for (size_t i = 0; i < pcs_vector.size(); i++)
		//      if (pcs_vector[i]->pcs_type_name.find("HEAT") != string::npos) {
		// TF
		if (pcs_vector[i]->getProcessType() == HEAT_TRANSPORT)
		{
			t_pcs = pcs_vector[i];
			break;
		}
	if (T_Flag)
	{
		idx_T0 = t_pcs->GetNodeValueIndex("TEMPERATURE1");
		idx_T1 = idx_T0 + 1;
		Temp = new double[max_nnodes_LE];
		T1 = new double[max_nnodes_LE];
	}

	if (enhanced_strain_dm && dim == 2)
	{
		NodesInJumpedA = new bool[max_nnodes_QE_2D];
		Ge = new Matrix(4, 2);
		Pe = new Matrix(2, 4);
		BDG = new Matrix(2, 2 * max_nnodes_QE_2D);
		PDB = new Matrix(2 * max_nnodes_QE_2D, 2);
		DtD = new Matrix(2, 2);
		PeDe = new Matrix(2, 4);
		X0 = new double[3];
		n_jump = new double[3];
		pr_stress = new double[3];
	}
	else
	{
		NodesInJumpedA = NULL;
		Ge = NULL;
		Pe = NULL;
		BDG = NULL;
		PDB = NULL;
		DtD = NULL;
		PeDe = NULL;
		X0 = NULL;
		n_jump = NULL;
		pr_stress = NULL;
	}

	//
	// Time unit factor
	time_unit_factor = pcs->time_unit_factor;

#if defined(USE_PETSC) // || defined(other parallel libs)//05~07.3013. WW
	idxm = new int[60]; //> global indices of local matrix rows
	idxn = new int[60]; //> global indices of local matrix columns
	local_idx = new int[60]; //> local index for local assemble
#endif
}

//  Constructor of class Element_DM
CFiniteElementVec::~CFiniteElementVec()
{
	delete B_matrix;
	delete B_matrix_T;
	delete[] dstress;
	delete[] dstrain;
	delete De;
	delete ConsistDep;
	delete AuxMatrix;
	delete AuxMatrix2; // NW
	delete[] Disp;
	delete[] Temp;
	delete[] T1;
	delete[] Sxx;
	delete[] Syy;
	delete[] Szz;
	delete[] Sxy;
	delete[] pstr;
	delete[] _nodal_strain_indices;
	delete[] _nodal_stress_indices;
	delete[] strain_ne;
	delete[] stress_ne;
	delete[] stress0;
	if (Sxz)
		delete[] Sxz;
	if (Syz)
		delete[] Syz;

	if (dynamic)
	{
		delete Mass;
		delete dAcceleration;
		Mass = NULL;
		dAcceleration = NULL;
	}

	if (pcs->Memory_Type == 0) // Do not store local matrices
	{
		delete Stiffness;
		delete RHS;
		Stiffness = NULL;
		RHS = NULL;
	}

	if (enhanced_strain_dm)
	{
		delete NodesInJumpedA;
		delete Ge;
		delete Pe;
		delete PeDe;
		delete BDG;
		delete PDB;
		delete DtD;

		NodesInJumpedA = NULL;
		Ge = NULL;
		Pe = NULL;
		PeDe = NULL;
		BDG = NULL;
		PDB = NULL;
		DtD = NULL;
	}

	// 11.07.2011. WW
	if (PressureC)
		delete PressureC;
	if (PressureC_S)
		delete PressureC_S;
	if (PressureC_S_dp)
		delete PressureC_S_dp;

	B_matrix = NULL;
	B_matrix_T = NULL;
	dstress = NULL;
	dstrain = NULL;
	De = NULL;
	ConsistDep = NULL;
	AuxMatrix = NULL;
	AuxMatrix2 = NULL; // NW
	Disp = NULL;
	Temp = NULL;
	T1 = NULL;
	Sxx = NULL;
	Syy = NULL;
	Szz = NULL;
	Sxy = NULL;
	Sxz = NULL;
	Syz = NULL;
	pstr = NULL;
	//  10.11.2010. WW
	if (X0)
		delete[] X0;
	if (n_jump)
		delete[] n_jump;
	if (pr_stress)
		delete[] pr_stress;
	if (Idx_Vel)
		delete[] Idx_Vel;
	X0 = n_jump = pr_stress = NULL;

	if (_nodal_p1)
		delete[] _nodal_p1;
	if (_nodal_p2)
		delete[] _nodal_p2;
	if (_nodal_cp0)
		delete[] _nodal_cp0;
	if (_nodal_dcp)
		delete[] _nodal_dcp;

	if (_nodal_S0)
		delete[] _nodal_S0;
	if (_nodal_S)
		delete[] _nodal_S;
	if (AuxNodal1)
		delete[] AuxNodal1;

	// NW
	for (int i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		delete vec_B_matrix[i];
		delete vec_B_matrix_T[i];
		vec_B_matrix[i] = NULL;
		vec_B_matrix_T[i] = NULL;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::SetMaterial(const int EleIndex)

   Aufgabe:
         Set material data for local assembly
   Formalparameter:
           E:

   Programming:
   11/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::SetMaterial()
{
	//......................................................................
	// MAT group
	int MatGroup = MeshElement->GetPatchIndex();
	//......................................................................
	// MSP
	smat = msp_vector[MatGroup];
	// WX:01.2013. time dependent E nv aniso
	if (smat->Time_Dependent_E_nv_mode == 2)
		smat->CalculateTransformMatrixFromNormalVector(ele_dim);
	smat->axisymmetry = pcs->m_msh->isAxisymmetry();
	// Single yield surface model
	if (smat->Plasticity_type == 2)
		smat->ResizeMatricesSYS(ele_dim);
	//......................................................................
	// MFP
	if (F_Flag)
	{
		m_mfp = MFPGet("LIQUID"); // YD
		if (!m_mfp)
			m_mfp = mfp_vector[0]; // OK
	}
	//......................................................................
	// MMP
	m_mmp = mmp_vector[MatGroup];
	//......................................................................
}

/**************************************************************************
   GeoSys - Function: SetMemory

   Aufgabe:
         Set memory for local matrices
   Programmaenderungen:
   01/2005   WW    Erste Version
   05/2005   WW    Dynamic analysis

**************************************************************************/
void CFiniteElementVec::SetMemory()
{
	int size = 0;
	ElementMatrix* EleMat = NULL;

	// Prepare local matrices
	if (pcs->Memory_Type == 0)
	{
		// If local matrices are not stored, resize the matrix
		size = dim * nnodesHQ;
		Stiffness->LimitSize(size, size);
		if (PressureC)
			PressureC->LimitSize(size, nnodes);
		if (PressureC_S)
			PressureC_S->LimitSize(size, nnodes);
		if (PressureC_S_dp)
			PressureC_S_dp->LimitSize(size, nnodes);
		RHS->LimitSize(size);
	}
	else
	{
		EleMat = pcs->Ele_Matrices[Index];
		Stiffness = EleMat->GetStiffness();
		RHS = EleMat->GetRHS();
		if (PressureC)
			PressureC = EleMat->GetCouplingMatrixA();
	}

	if (dynamic)
	{
		Mass->LimitSize(nnodesHQ, nnodesHQ);
		dAcceleration->LimitSize(nnodesHQ * dim);
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setB_Matrix(const int LocalIndex)

   Aufgabe:
          Form B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setB_Matrix(const int LocalIndex)
{
	switch (dim)
	{
		case 2:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_12, 0.0
			(*B_matrix)(0, 1) = 0.0;

			if (axisymmetry) // Axisymmtry
			{
				// B_21, N/r
				(*B_matrix)(1, 0) = shapefctHQ[LocalIndex] / Radius;
				// B_22, 0.0
				(*B_matrix)(1, 1) = 0.0;
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, dN/dz
				(*B_matrix)(2, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			}
			else
			{
				// B_21, 0.0
				(*B_matrix)(1, 0) = 0.0;
				// B_22, dN/dy
				(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, 0.0
				(*B_matrix)(2, 1) = 0.0;
			}
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];

			break;
		case 3:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_22, dN/dy
			(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_33, dN/dz
			(*B_matrix)(2, 2) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			//
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];
			//
			// B_51, dN/dz
			(*B_matrix)(4, 0) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_53, dN/dx
			(*B_matrix)(4, 2) = dshapefctHQ[LocalIndex];
			//
			// B_62, dN/dz
			(*B_matrix)(5, 1) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_63, dN/dy
			(*B_matrix)(5, 2) = dshapefctHQ[nnodesHQ + LocalIndex];

			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setTransB_Matrix(const int LocalIndex)

   Aufgabe:
          Form the tanspose of B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setTransB_Matrix(const int LocalIndex)
{
	setB_Matrix(LocalIndex);
	B_matrix->GetTranspose(*B_matrix_T);
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ComputeStrain(const int ip)

   Aufgabe:
          Compute strains

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::ComputeStrain(const int ip)
{
	int i, j = 0, k = 0;
	if (excavation) // WX:03.2012 if element is excavated, strain = 0
	{
		for (i = 0; i < ns; i++)
			dstrain[i] = 0.;
		return;
	}
	switch (dim)
	{
		case 2:
			for (i = 0; i < ns; i++)
				dstrain[i] = 0.0;
			if (axisymmetry)
			{
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[i] * shapefctHQ[i];
					dstrain[2] += Disp[j] * dshapefctHQ[j];
					dstrain[3] += Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}

				calculateRadius(ip);
				dstrain[1] /= Radius;
			}
			else
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[j] * dshapefctHQ[j];
					dstrain[3] += Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}
			break;
		case 3:
			for (i = 0; i < ns; i++)
				dstrain[i] = 0.0;
			for (i = 0; i < nnodesHQ; i++)
			{
				j = i + nnodesHQ;
				k = i + 2 * nnodesHQ;
				dstrain[0] += Disp[i] * dshapefctHQ[i];
				dstrain[1] += Disp[j] * dshapefctHQ[j];
				dstrain[2] += Disp[k] * dshapefctHQ[k];
				dstrain[3] += Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				dstrain[4] += Disp[i] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[i];
				dstrain[5] += Disp[j] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[j];
			}
			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: CalDensity()

   Aufgabe:
          Caculate density of porous medium
   Programming:
   05/2005     WW        Erste Version
 **************************************************************************/
double CFiniteElementVec::CalDensity()
{
	double rho;
	// OK_MFP
	//--------------------------------------------------------------------
	// MFP fluid properties
	double density_fluid = 0.0;
	double porosity = 0.0;
	// double p_g = 0.0;
	double Sw = 0.0;
	int no_phases = (int)mfp_vector.size();
	int i = 0, phase = 0;

	rho = 0.0;
	if (F_Flag)
	{
		if ((no_phases > 0) && (no_phases > phase))
			density_fluid = m_mfp->Density();

		// OK_MMP
		//--------------------------------------------------------------------
		// MMP medium properties
		porosity = m_mmp->Porosity(this);
		// Assume solid density is constant. (*smat->data_Density)(0)
		if (smat->Density() > 0.0)
		{
			Sw = 1.0; // JT, should be 1.0, unless multiphase (calculate below) (if unsaturated, fluid density would be
			          // negligible... so still works)
			if (Flow_Type > 0 && Flow_Type != 10)
			{
				Sw = 0.; // WW
				for (i = 0; i < nnodes; i++)
					Sw += shapefct[i] * _nodal_S[i];
			}
			rho = (1. - porosity) * fabs(smat->Density()) + porosity * Sw * density_fluid;

			if (Flow_Type == 2 || Flow_Type == 3)
			{
				CFluidProperties* GasProp;
				GasProp = MFPGet("GAS");
				rho += porosity * (1.0 - Sw) * GasProp->Density();
			}
		}
		else
			rho = 0.0;
	}
	else
	    // If negative value is given in the .msp file, gravity by solid is skipped
	    if (smat->Density() > 0.0)
		rho = smat->Density();
	return rho;
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: ComputeMatrix_RHS(const double fkt)
   Aufgabe:
           as the name

   Programming:
   07/2004   WW
   08/2004   OK   MFP implementation
   01/2010   NW   use chache of B,B^T matrices
 **************************************************************************/
void CFiniteElementVec::ComputeMatrix_RHS(const double fkt, const Matrix* p_D)
{
	int i, j, k, l;
	double rho, fac, fac1, fac2, dN_dx, f_buff;
	fac = fac1 = fac2 = f_buff = 0.0;
	dN_dx = 0.0;
	rho = CalDensity();
	const int nnodesHQ = this->nnodesHQ;
	const int nnodes = this->nnodes;
	const int ele_dim = this->ele_dim;
	const int ns = this->ns;

	// NW cache B, B^T
	for (i = 0; i < nnodesHQ; i++)
	{
		setTransB_Matrix(i);
		(*this->vec_B_matrix[i]) = *B_matrix;
		(*this->vec_B_matrix_T[i]) = *B_matrix_T;
	}
	Matrix* old_B_matrix = B_matrix;
	Matrix* old_B_matrix_T = B_matrix_T;

	Matrix* tmp_B_matrix = NULL;
	Matrix* tmp_B_matrix_T = NULL;
	Matrix* tmp_AuxMatrix = AuxMatrix;
	Matrix* tmp_AuxMatrix2 = AuxMatrix2;
	Matrix* tmp_Stiffness = Stiffness;

	for (i = 0; i < nnodesHQ; i++)
	{
		// NW      setTransB_Matrix(i);
		tmp_B_matrix_T = this->vec_B_matrix_T[i];
		// Local assembly of A*u=int(B^t*sigma) for Newton-Raphson method
		for (j = 0; j < ele_dim; j++)
			for (k = 0; k < ns; k++)
				(*RHS)[j * nnodesHQ + i] += (*tmp_B_matrix_T)(j, k) * (dstress[k] - stress0[k]) * fkt;
		// TEST             (*B_matrix_T)(j,k)*dstress[k]*fkt;
		if (PreLoad == 11)
			continue;
		if (excavation)
			continue; // WX:08.2011
// Local assembly of stiffness matrix, B^T C B
#ifdef JFNK_H2M
		/// If JFNK. 18.10.2010. WW
		if (pcs->m_num->nls_method == 2 && (!pcs->JFNK_precond))
			continue;
#endif
		(*tmp_AuxMatrix2) = 0.0;
		// NW
		tmp_B_matrix_T->multi(*p_D, *tmp_AuxMatrix2);
		for (j = 0; j < nnodesHQ; j++)
		{
			// NW          setB_Matrix(j);
			tmp_B_matrix = this->vec_B_matrix[j];
			// Compute stiffness matrix
			(*tmp_AuxMatrix) = 0.0;
			tmp_AuxMatrix2->multi(*tmp_B_matrix, *tmp_AuxMatrix);
			// NW          B_matrix_T->multi(*p_D, *B_matrix, *AuxMatrix);

			// Local assembly of stiffness matrix
			for (k = 0; k < ele_dim; k++)
			{
				const int kia = i + k * nnodesHQ;
				for (l = 0; l < ele_dim; l++)
					(*tmp_Stiffness)(kia, j + l * nnodesHQ) += (*tmp_AuxMatrix)(k, l) * fkt;
			}
		} // loop j
	} // loop i

	// should restore pointer NW
	B_matrix = old_B_matrix;
	B_matrix_T = old_B_matrix_T;

	//---------------------------------------------------------
	// Assemble coupling matrix
	//---------------------------------------------------------
	// LoadFactor: factor of incremental loading, prescibed in rf_pcs.cpp

	// 07.2011 WW
	if ((PressureC || PressureC_S || PressureC_S_dp) && !PreLoad)
	{
		fac = LoadFactor * fkt;

		// 07.2011. WW
		if (PressureC_S || PressureC_S_dp)
		{
			// Pressure 1
			fac2 = interpolate(_nodal_p1);
			// Saturation of phase 1
			fac1 = m_mmp->SaturationCapillaryPressureFunction(fac2);
			if (PressureC_S_dp)
				fac2 = fac1 - fac2 * m_mmp->PressureSaturationDependency(fac1, true);
			// JT: dSdP now returns actual sign (<0)
		}

		if (axisymmetry)
		{
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
					for (j = 0; j < ele_dim; j++)
					{
						dN_dx = dshapefctHQ[nnodesHQ * j + k];
						if (j == 0)
							dN_dx += shapefctHQ[k] / Radius;

						f_buff = fac * dN_dx * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
						if (PressureC_S)
							(*PressureC_S)(nnodesHQ* j + k, l) += f_buff * fac1;
						if (PressureC_S_dp)
							(*PressureC_S_dp)(nnodesHQ* j + k, l) += f_buff * fac2;
					}
			}
		}
		else
		{
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
					for (j = 0; j < ele_dim; j++)
					{
						f_buff = fac * dshapefctHQ[nnodesHQ * j + k] * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
						if (PressureC_S)
							(*PressureC_S)(nnodesHQ* j + k, l) += f_buff * fac1;
						if (PressureC_S_dp)
							(*PressureC_S_dp)(nnodesHQ* j + k, l) += f_buff * fac2;
					}
			}
		}
	}
	//---------------------------------------------------------
	// Assemble gravity force vector
	//---------------------------------------------------------
	if (rho > 0.0 && GravityForce)
	{
		// 2D, in y-direction
		// 3D, in z-direction
		i = (ele_dim - 1) * nnodesHQ;
		double timeFactor = 1.;
		if (smat->gravity_ramp)
		{
			const int curveIndex = smat->grav_curve_id;
			const int curveMethod = 0;

			int valid;
			timeFactor = GetCurveValue(curveIndex, curveMethod, aktuelle_zeit, &valid);
		}

		const double coeff = LoadFactor * rho * smat->grav_const * fkt * timeFactor;
		for (k = 0; k < nnodesHQ; k++)
		{
			(*RHS)[i + k] += coeff * shapefctHQ[k];
			//        (*RHS)(i+ka) += LoadFactor * rho * smat->grav_const * shapefctHQ[ka] * fkt;
		}
	}
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: LocalAssembly
   Aufgabe:
           Compute the local finite element matrices
   Formalparameter:
           E:
         const int update  : indicator to update stress and strain only

   Programming:
   06/2004   WW   Generalize for different element types as a member of class
   05/2005   WW   ...
   08/2010   WW   JFNK method
 **************************************************************************/
void CFiniteElementVec::LocalAssembly(const int update)
{
	int j;
	double* a_n = NULL;

	Index = MeshElement->GetIndex();
	SetMemory();
	SetMaterial();
	// 12.2009. WW
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];

#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	if (m_dom)
	{ // Moved here from GlobalAssemly. 08.2010. WW
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
		b_rhs = m_dom->eqsH->b;
#else
		b_rhs = m_dom->eqs->b;
#endif
	}
	else
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	{
#if defined(USE_PETSC) // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
		b_rhs = pcs->eqs_new->b;
#else
		b_rhs = pcs->eqs->b;
#endif
	}
	(*RHS) = 0.0;
	(*Stiffness) = 0.0;
	// 07.2011. WW
	if (PressureC)
		(*PressureC) = 0.0;
	if (PressureC_S)
		(*PressureC_S) = 0.0;
	if (PressureC_S_dp)
		(*PressureC_S_dp) = 0.0;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	if (m_dom)
	{
		// 06.2011. WW
		for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
			NodeShift[i] = m_dom->shift[i];
		for (int i = 0; i < nnodesHQ; i++)
			eqs_number[i] = element_nodes_dom[i];
	}
	else
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	{
		for (int i = 0; i < nnodesHQ; i++)
			eqs_number[i] = MeshElement->nodes[i]->GetEquationIndex();
	}
#endif
	// For strain and stress extropolation all element types
	// Number of elements associated to nodes
	for (int i = 0; i < nnodes; i++)
		dbuff[i] = (double)MeshElement->nodes[i]->getConnectedElementIDs().size();

	// Get displacement_n
	if (dynamic)
	{
		a_n = pcs->GetAuxArray();
		for (size_t i = 0; i < dim; i++)
			for (j = 0; j < nnodesHQ; j++)
			{
				// Increment of acceleration, da
				(*dAcceleration)(i* nnodesHQ + j) = pcs->GetNodeValue(nodes[j], Idx_dm0[i]);
				// Increment of displacement
				// du = v_n*dt+0.5*a_n*dt*dt+0.5*beta2*da*dt*dt
				// a_n = a_{n+1}-da
				Disp[j + i * nnodesHQ]
				    = pcs->GetNodeValue(nodes[j], Idx_Vel[i]) * dt
				      + 0.5 * dt * dt * (a_n[nodes[j] + NodeShift[i]] + beta2 * (*dAcceleration)(i* nnodesHQ + j));
			}
	}
	else
		for (size_t i = 0; i < dim; i++)
			for (j = 0; j < nnodesHQ; j++)
				// WX:03.2013 use total disp. if damage or E=f(t) is on, dstress in LocalAssembly_continumm() is also
				// changed
				if (smat->Time_Dependent_E_nv_mode > MKleinsteZahl && pcs->ExcavMaterialGroup < 0)
					Disp[j + i * nnodesHQ] = pcs->GetNodeValue(nodes[j], Idx_dm0[i])
					                         + pcs->GetNodeValue(nodes[j], Idx_dm0[i] + 1);
				else
					Disp[j + i * nnodesHQ] = pcs->GetNodeValue(nodes[j], Idx_dm0[i]);

	if (_nodal_p1)
	{
		for (int i = 0; i < nnodes; i++)
		{
			_nodal_p1[i] = h_pcs->GetNodeValue(nodes[i], idx_P1);
		}
	}
	if (_nodal_p2)
	{
		for (int i = 0; i < nnodes; i++)
		{
			_nodal_p2[i] = h_pcs->GetNodeValue(nodes[i], idx_P2);
		}
	}
	if (_nodal_S0)
	{
		for (int i = 0; i < nnodes; i++)
		{
			_nodal_S0[i] = h_pcs->GetNodeValue(nodes[i], idx_S0);
		}
	}
	if (_nodal_S)
	{
		for (int i = 0; i < nnodes; i++)
		{
			_nodal_S[i] = h_pcs->GetNodeValue(nodes[i], idx_S);
		}
	}

	// Get saturation of element nodes
	if (Flow_Type > 0 && Flow_Type != 10)
	{
		// 12.03.2008 WW
		if ((Flow_Type == 1 || Flow_Type == 2) && (smat->SwellingPressureType == 3 || smat->SwellingPressureType == 4))
		{
			double fac = 1.0;
			if (Flow_Type == 1)
				fac = -1.0;
			for (int i = 0; i < nnodes; i++)
			{
				// Pc
				_nodal_cp0[i] = fac * h_pcs->GetNodeValue(nodes[i], idx_P1 - 1);
				// dPc
				_nodal_dcp[i] = fac * (h_pcs->GetNodeValue(nodes[i], idx_P1) - h_pcs->GetNodeValue(nodes[i], idx_P1 - 1));
			}
		}
	}
	//

	// -------------------------------12.2009.  WW
	if (pcs->ite_steps == 1)
	{
		excavation = false;
		if ((smat->excavation > 0 || pcs->ExcavMaterialGroup > -1) && MeshElement->GetMark())
		{
			int valid;
			if (smat->excavation > 0)
			{
				if (GetCurveValue(smat->excavation, 0, aktuelle_zeit, &valid) < 1.0)
				{
					excavation = true;
					smat->excavated = true; // To be ... 12.2009. WW
					*(eleV_DM->Stress) = 0.;
				}
				else
					smat->excavated = false; // To be ... 12.2009. WW
			}
			// WX:07.2011
			if (static_cast<size_t>(pcs->ExcavMaterialGroup) == MeshElement->GetPatchIndex())
			{
				double const* ele_center(MeshElement->GetGravityCenter());
				if ((GetCurveValue(pcs->ExcavCurve, 0, aktuelle_zeit, &valid) + pcs->ExcavBeginCoordinate)
				        > (ele_center[pcs->ExcavDirection])
				    && (ele_center[pcs->ExcavDirection] - pcs->ExcavBeginCoordinate) > -0.001)
				{
					excavation = true;
					*(eleV_DM->Stress) = 0.;
					// MeshElement->SetExcavState(1); //WX:03.2012
				}
			}
		}
	}
	//----------------------------------------------------

	if (enhanced_strain_dm && ele_value_dm[MeshElement->GetIndex()]->Localized)
		LocalAssembly_EnhancedStrain(update);
	else
		LocalAssembly_continuum(update);

	if (update == 0)
	{
		if (dynamic)
			ComputeMass();
		GlobalAssembly();
		// Output matrices
		if (pcs->Write_Matrix)
		{
			(*pcs->matrix_file) << "### Element: " << Index << "\n";
			(*pcs->matrix_file) << "---Stiffness matrix: "
			                    << "\n";
			Stiffness->Write(*pcs->matrix_file);
			(*pcs->matrix_file) << "---RHS: "
			                    << "\n";
			RHS->Write(*pcs->matrix_file);
			if (PressureC)
			{
				(*pcs->matrix_file) << "Pressue coupling matrix: "
				                    << "\n";
				PressureC->Write(*pcs->matrix_file);
			}
			// 07.2011. WW
			if (PressureC_S)
			{
				(*pcs->matrix_file) << "Saturation depedent pressue coupling matrix: "
				                    << "\n";
				PressureC_S->Write(*pcs->matrix_file);
			}
			if (PressureC_S_dp)
			{
				(*pcs->matrix_file) << "Jacobi pressue coupling matrix: "
				                    << "\n";
				PressureC_S_dp->Write(*pcs->matrix_file);
			}
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
   12/2009   WW New excavtion approach
 **************************************************************************/
bool CFiniteElementVec::GlobalAssembly()
{
	/*// For excavation simulation. 12.2009. WW
	int valid = 0;
	if (excavation)
	{
	    excavation = true;
	    bool onExBoundary = false;

	    CNode* node;
	    CElem* elem;
	    CSolidProperties* smat_e;

	    for (int i = 0; i < nnodesHQ; i++)
	    {
	        node = MeshElement->nodes[i];
	        onExBoundary = false;
	        const size_t n_elements (node->getConnectedElementIDs().size());
	        for (size_t j = 0; j < n_elements; j++)
	        {
	            elem = pcs->m_msh->ele_vector[node->getConnectedElementIDs()[j]];
	            if (!elem->GetMark())
	                continue;

	            smat_e = msp_vector[elem->GetPatchIndex()];
	            if (smat_e->excavation > 0)
	            {
	                if (fabs(GetCurveValue(smat_e->excavation, 0,
	                                       aktuelle_zeit,
	                                       &valid) - 1.0) < DBL_MIN)
	                {
	                    onExBoundary = true;
	                    break;
	                }
	            }
	            else if(pcs->ExcavMaterialGroup > -1)
	            {
	                double const* ele_center(elem->GetGravityCenter());
	                if((GetCurveValue(pcs->ExcavCurve,0,aktuelle_zeit,
	                                  &valid) + pcs->ExcavBeginCoordinate) <
	                   (ele_center[pcs->ExcavDirection]))
	                {
	                    onExBoundary = true;
	                    break;
	                }
	                else if (elem->GetPatchIndex() != static_cast<size_t>(pcs->ExcavMaterialGroup))
	                {
	                    onExBoundary = true;
	                    break;
	                }
	            } //WX:07.2011
	            else
	            {
	                onExBoundary = true;
	                break;
	            }
	        }

	        if (!onExBoundary)
	            for (size_t j = 0; j < dim; j++)
	                (*RHS)(j * nnodesHQ + i) = 0.0;
	    }
	}*/

	GlobalAssembly_RHS();
	if (PreLoad == 11)
		return true;

	// For excavation simulation. 12.2009. WW
	if (excavation) // WX: modify
	{
		// MeshElement->MarkingAll(false);
		//*(eleV_DM->Stress) = 0.;
		//*(eleV_DM->Stress0) = 0.;
		// if (eleV_DM->Stress_j)
		//(*eleV_DM->Stress_j) = 0.0;
		return false;
	}

	GlobalAssembly_Stiffness();

	return true;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_Stiffness()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
 **************************************************************************/
#if defined(USE_PETSC) // || defined(other parallel libs)//06.2013. WW
void CFiniteElementVec::GlobalAssembly_Stiffness()
{
	const int dof = ele_dim;
	const int m_dim = nnodesHQ * dof;
	const int n_dim = m_dim;

	int non_ghost_flag = 1;
	if (act_nodes_h != nnodesHQ)
		non_ghost_flag = -1;

	petsc_group::PETScLinearSolver* eqs = pcs->eqs_new;

	// Set row and column indices.
	// For non-ghost element, all row and column indices are positive.
	// For ghost elemnent, all column indices, saved in idxn, are positive,
	// while, all row indices, saved in idxm, are negative. For the zero
	// entry of idxm, an alias of it is indroduced to represent zero for the
	// negative indexing.
	const PetscInt zero_id_alias = -eqs->Size();
	for (int i = 0; i < nnodesHQ; i++)
	{
		const int i_buff = MeshElement->nodes[i]->GetEquationIndex() * dof;
		for (int k = 0; k < dof; k++)
		{
			const int ki = k * nnodesHQ + i;
			idxm[ki] = non_ghost_flag * (i_buff + k);
			idxn[ki] = non_ghost_flag * idxm[ki]; // always positive
			if (non_ghost_flag == -1 && idxm[ki] == 0)
				idxm[ki] = zero_id_alias;
		}
	}

	// If this is a ghost element, the non-ghost entries of idxm, which holds
	// negative values, are reset to positive. The original zero entry,
	// which represents by a special nagetive value, is reset to zero again.
	if (non_ghost_flag == -1)
	{
		// Reset non-ghost entries of idxm to positive.
		for (int i = 0; i < act_nodes_h; i++)
		{
			for (int k = 0; k < dof; k++)
			{
				const int ki = k * nnodesHQ + local_idx[i];
				// Reset original zero entry to zero.
				if (idxm[ki] == zero_id_alias)
					idxm[ki] = 0;
				else
					idxm[ki] *= -1;
			}
		}
	}

	const double* local_matrix = Stiffness->getEntryArray();
	(*RHS) *= -1;
	double* local_vec = RHS->getEntryArray();

	eqs->addMatrixEntries(m_dim, idxm, n_dim, idxn, local_matrix);
	eqs->setArrayValues(1, m_dim, idxm, local_vec);
}
#else
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_Stiffness()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_Stiffness()
{
	int i, j;
	double f1, f2;
	f1 = 1.0;
	f2 = -1.0;

#if defined(NEW_EQS) && defined(JFNK_H2M)
	/// If not JFNK. 02.2011. WW
	if (pcs->m_num->nls_method == 2)
	{
		if (!pcs->JFNK_precond)
			return;

		long kk = 0;
		// Assemble Jacobi preconditioner
		for (size_t k = 0; k < ele_dim; k++)
		{
			for (size_t l = 0; l < ele_dim; l++)
				for (i = 0; i < nnodesHQ; i++)
				{
					kk = eqs_number[i] + NodeShift[k];
					for (j = 0; j < nnodesHQ; j++)
					{
						if (kk != eqs_number[j] + NodeShift[l])
							continue;
						pcs->eqs_new->prec_M[kk] += (*Stiffness)(i + k * nnodesHQ, j + l * nnodesHQ);
					}
				}
			// loop l
		} // loop k
		return;
	}
#endif
	double biot = 1.0;
	biot = smat->biot_const;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;
	if (m_dom)
		A = m_dom->eqsH->A;
	else
		A = pcs->eqs_new->A;
#endif

	if (dynamic)
	{
		f1 = 0.5 * beta2 * dt * dt;
		f2 = -0.5 * bbeta1 * dt;
		// Assemble stiffness matrix
		for (i = 0; i < nnodesHQ; i++)
		{
			for (j = 0; j < nnodesHQ; j++)
			{
				// Local assembly of stiffness matrix
				for (size_t k = 0; k < ele_dim; k++)
				{
#ifdef NEW_EQS
					(*A)(eqs_number[i] + NodeShift[k], eqs_number[j] + NodeShift[k]) += (*Mass)(i, j);
#else
					MXInc(eqs_number[i] + NodeShift[k], eqs_number[j] + NodeShift[k], (*Mass)(i, j));
#endif
				}
			} // loop j
		} // loop i
	}

	// Assemble stiffness matrix
	for (i = 0; i < nnodesHQ; i++)
	{
		for (j = 0; j < nnodesHQ; j++)
		{
			// Local assembly of stiffness matrix
			for (size_t k = 0; k < ele_dim; k++)
			{
				for (size_t l = 0; l < ele_dim; l++)
				{
#ifdef NEW_EQS
					(*A)(eqs_number[i] + NodeShift[k], eqs_number[j] + NodeShift[l])
					    += f1 * (*Stiffness)(i + k * nnodesHQ, j + l * nnodesHQ);
#else
					MXInc(eqs_number[i] + NodeShift[k],
					      eqs_number[j] + NodeShift[l],
					      f1 * (*Stiffness)(i + k * nnodesHQ, j + l * nnodesHQ));
#endif
				}
			}
		} // loop j
	} // loop i

	// TEST OUT
	// Stiffness->Write();
	if (pcs->type / 40 != 1) // Not monolithic scheme
		return;

	if (PressureC)
	{
		i = 0; // phase
		if (Flow_Type == 2) // Multi-phase-flow
			i = 1;
		GlobalAssembly_PressureCoupling(PressureC, f2 * biot, i);
	}
	// H2: p_g- S_w*p_c
	if (PressureC_S)
		GlobalAssembly_PressureCoupling(PressureC_S, -f2 * biot, 0);
	if (PressureC_S_dp)
		GlobalAssembly_PressureCoupling(PressureC_S_dp, -f2 * biot, 0);
}
#endif //#if defined(USE_PETSC) // || defined(other parallel libs)//07.2013. WW
//--------------------------------------------------------------------------
/*!
   \brief Assembe the pressure coupling matrix

    to the global stiffness matrix in the monolithic scheme

   \param pCMatrix: the matrix
   \param fct: factor
   \param phase: phasse index

    07.2011. WW
 */
#if defined(USE_PETSC) // || defined(other parallel libs)//10.3012. WW
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix*, double, const int)
{
}
#else
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix* pCMatrix, double fct, const int phase)
{
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;
	if (m_dom)
		A = m_dom->eqsH->A;
	else
		A = pcs->eqs_new->A;
#endif

	int dim_shift = dim + phase;
	// Add pressure coupling matrix to the stifness matrix
	for (int i = 0; i < nnodesHQ; i++)
	{
		for (int j = 0; j < nnodes; j++)
		{
			for (size_t k = 0; k < ele_dim; k++)
			{
#ifdef NEW_EQS
				(*A)(NodeShift[k] + eqs_number[i], NodeShift[dim_shift] + eqs_number[j])
				    += fct * (*pCMatrix)(nnodesHQ* k + i, j);
#else
				MXInc(NodeShift[k] + eqs_number[i], NodeShift[dim_shift] + eqs_number[j],
				      fct * (*pCMatrix)(nnodesHQ* k + i, j));
#endif
			}
		}
	}
}
#endif

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ComputeMass()
   Aufgabe:
           Compute the mass matrix for dynamic analyses
   Programming:
   05/2005   WW   Elastische Elemente
 **************************************************************************/
void CFiniteElementVec::ComputeMass()
{
	int i, j;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp = 0;
	gp_t = 0;
	double fkt = 0.0;

	(*Mass) = 0.0;
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		getShapefunctValues(gp, 1); // need for density calculation
		getShapefunctValues(gp, 2); // Quadratic interpolation function
		fkt *= CalDensity();

		for (i = 0; i < nnodesHQ; i++)
			for (j = 0; j < nnodesHQ; j++)
			{
				if (i > j)
					continue;
				(*Mass)(i, j) += fkt * shapefctHQ[i] * shapefctHQ[j];
			}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_RHS()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
   05/2005   WW dyn
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_RHS()
{
	bool Residual = false;
	if (Flow_Type >= 0)
	{
		if (pcs->type / 10 == 4) // Monolithic scheme
		{
			// If nonlinear deformation
			if (pcs_deformation > 100)
				Residual = true;
		}
		else // Partitioned scheme
			Residual = true;
	}
	if (dynamic)
		Residual = true;

	bool onExBoundaryState[max_nnodes_QE_3D] = {false};
	if (excavation)
	{
		int valid = 0;
		excavation = true;
		bool onExBoundary = false;
		CNode* node;
		CElem* elem;
		CSolidProperties* smat_e;

		for (int i = 0; i < nnodesHQ; i++)
		{
			node = MeshElement->nodes[i];
			onExBoundary = false;
			const std::size_t n_elements(node->getConnectedElementIDs().size());
			for (std::size_t j = 0; j < n_elements; j++)
			{
				elem = pcs->m_msh->ele_vector[node->getConnectedElementIDs()[j]];
				if (!elem->GetMark())
					continue;

				smat_e = msp_vector[elem->GetPatchIndex()];
				if (smat_e->excavation > 0)
				{
					if (fabs(GetCurveValue(smat_e->excavation, 0, aktuelle_zeit, &valid) - 1.0) < DBL_MIN)
					{
						onExBoundary = true;
						break;
					}
				}
				else if (pcs->ExcavMaterialGroup > -1)
				{
					double const* ele_center(elem->GetGravityCenter());
					if ((GetCurveValue(pcs->ExcavCurve, 0, aktuelle_zeit, &valid) + pcs->ExcavBeginCoordinate)
					    < (ele_center[pcs->ExcavDirection]))
					{
						onExBoundary = true;
						break;
					}
					else if (elem->GetPatchIndex() != static_cast<size_t>(pcs->ExcavMaterialGroup))
					{
						onExBoundary = true;
						break;
					}
				}
				else
				{
					onExBoundary = true;
					break;
				}
			}
			if (onExBoundary)
				onExBoundaryState[i] = 1;
		}
	}

	if (Residual)
	{
		const double biot = smat->biot_const;
		double nodal_pore_p[max_nnodes_LE];
		const int dim_times_nnodesHQ(dim * nnodesHQ);
		for (int i = 0; i < dim_times_nnodesHQ; i++)
			AuxNodal1[i] = 0.0;
		switch (Flow_Type)
		{
			//case 10: // Ground_flow. Will be merged to case 0
			case 0:  // Liquid flow
				// For monolithic scheme and liquid flow, the limit of positive pressure must be removed
				for (int i = 0; i < nnodes; i++)
					nodal_pore_p[i] = LoadFactor * _nodal_p1[i];
				if (excavation)
				{
					for (int i = 0; i < nnodes; i++)
					{
						if (onExBoundaryState[i] == 1) // WX:02.2013
							nodal_pore_p[i] = 0.0;
					}
				}
				if (pcs->Neglect_H_ini == 2)
				{
					for (int i = 0; i < nnodes; i++)
						nodal_pore_p[i] -= LoadFactor * h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
				}
				// If dynamic
				if (dynamic)
				{
					const double fact = bbeta1 * dt;
					double const* const	a_n = pcs->GetAuxArray();
					for (int i = 0; i < nnodes; i++)
					{
						nodal_pore_p[i] *= fact;
						nodal_pore_p[i] += dt * a_n[nodes[i] + NodeShift[problem_dimension_dm]]
						               + pcs->GetNodeValue(nodes[i], idx_P);
					}
				}
				PressureC->multi(nodal_pore_p, AuxNodal1);
				break;
			case 1: // Richards flow
				{
					if (smat->bishop_model < 0) // Without Bishop
					{
						for (int i = 0; i < nnodes; i++)
						{
							nodal_pore_p[i] = (biot < 0.0 && _nodal_p1[i] < 0.0) ?
								 0.0 : LoadFactor * S_Water * _nodal_p1[i];
						}
						if (excavation)
						{
							for (int i = 0; i < nnodes; i++)
							{
								if (onExBoundaryState[i] == 1) // WX:02.2013
									nodal_pore_p[i] = 0.0;
							}
						}
						if (pcs->Neglect_H_ini == 2)
						{
							for (int i = 0; i < nnodes; i++)
							{
								const double p0 = h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
								const double Sat0 = LoadFactor * m_mmp->SaturationCapillaryPressureFunction(-p0);
								nodal_pore_p[i] -= (p0 > 0.0) ? LoadFactor * Sat0 * p0 : 0.0;
							}
						}
						PressureC->multi(nodal_pore_p, AuxNodal1);
						break;
					}

					// Has Bishop model
					{
						for (int i = 0; i < nnodes; i++)
						{
							const double val_n = _nodal_p1[i];
							if (biot < 0.0 &&  val_n < 0.0)
								nodal_pore_p[i] = 0.0;
							else
							{
								const double S_e = m_mmp->GetEffectiveSaturationForPerm(_nodal_S[i], 0);
								nodal_pore_p[i] = LoadFactor * smat->getBishopCoefficient(S_e, val_n) * val_n;
							} // WX:12.2012 end if(biot<0.0&&val_n<0.0) else
						}

						if (excavation)
						{
							for (int i = 0; i < nnodes; i++)
							{
								if (onExBoundaryState[i] == 1)
									nodal_pore_p[i] = 0.0;
							}
						}

						if (pcs->Neglect_H_ini == 2)
						{
							for (int i = 0; i < nnodes; i++)
							{
								const double p0 = h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
								const double Sw0 = m_mmp->SaturationCapillaryPressureFunction(-p0);
								const double S_e0 =  m_mmp->GetEffectiveSaturationForPerm(Sw0, 0);
								nodal_pore_p[i] -= LoadFactor * smat->getBishopCoefficient(S_e0, p0) * p0;
							}
						}
					}
					PressureC->multi(nodal_pore_p, AuxNodal1);
				}
				break;
			case 2:
				{ // Multi-phase-flow: p_g-Sw*p_c
					// 07.2011. WW
					const int dim_times_nnodesHQ(dim * nnodesHQ);
					for (int i = 0; i < dim_times_nnodesHQ; i++)
						AuxNodal1[i] = 0.0;

					if (smat->bishop_model < 0) // No bishop model
					{
						if (pcs->Neglect_H_ini == 2)
						{
							for (int i = 0; i < nnodes; i++)
							{
								_nodal_p1[i] -= h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
								_nodal_p2[i] -= h_pcs->GetNodeValue(nodes[i], idx_p2_ini);
							}
						}

						PressureC->multi(_nodal_p2, AuxNodal1, LoadFactor);
						PressureC_S->multi(_nodal_p1, AuxNodal1, -1.0 * LoadFactor);
						break;
					}

					// Has bishop model
					for (int i = 0; i < nnodes; i++)
					{
						const double S_e = m_mmp->GetEffectiveSaturationForPerm(_nodal_S[i], 0);
						const double bishop_coef = smat->getBishopCoefficient(S_e, _nodal_p1[i]);
						const double pore_p = _nodal_p2[i] - bishop_coef * _nodal_p1[i];
						nodal_pore_p[i] = (biot < 0.0 &&  pore_p < 0.0) ? 0. : pore_p * LoadFactor;
					}

					if (excavation)
					{
						for (int i = 0; i < nnodes; i++)
						{
							if (onExBoundaryState[i] == 1) // WX:02.2013
							nodal_pore_p[i] = 0.;
						}
					}

					if (pcs->Neglect_H_ini == 2)
					{
						for (int i = 0; i < nnodes; i++)
						{
							const double p0 = h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
							const double Sw0 = m_mmp->SaturationCapillaryPressureFunction(p0);
							const double S_e0 =  m_mmp->GetEffectiveSaturationForPerm(Sw0, 0);
							const double bishop_coef = smat->getBishopCoefficient(S_e0, p0);
							const double pore_p0 = h_pcs->GetNodeValue(nodes[i], idx_p2_ini)
								   - bishop_coef * h_pcs->GetNodeValue(nodes[i], idx_p1_ini);
							nodal_pore_p[i] -= (biot < 0.0 && pore_p0 < 0.) ? 0. : bishop_coef * LoadFactor;
						}
					}

					PressureC->multi(nodal_pore_p, AuxNodal1);

					break;
				}
				case 3: // Multi-phase-flow: SwPw+SgPg	// PCH 05.05.2009
				{
					for (int i = 0; i < nnodes; i++)
					{
						double Snw = h_pcs->GetNodeValue(nodes[i], idx_Snw);
						double Sw = 1.0 - Snw;
						double Pw = _nodal_p1[i];
						double Pnw = _nodal_p2[i];
						const double val_n = Sw * Pw + Snw * Pnw;
						nodal_pore_p[i] = (biot < 0.0 && val_n < 0.0) ? 0.0 : val_n * LoadFactor;
					}
					PressureC->multi(nodal_pore_p, AuxNodal1);
					break;
				}

			} // end switch

			for (int i = 0; i < dim_times_nnodesHQ; i++)
				(*RHS)[i] -= fabs(biot) * AuxNodal1[i];
		} // End if partioned

		// If dynamic
		if (dynamic)
		{
			double const* const	a_n = pcs->GetAuxArray();
			for (std::size_t i = 0; i < dim; i++)
				for (int j = 0; j < nnodesHQ; j++)
					for (int k = 0; k < nnodesHQ; k++)
						(*RHS)[i * nnodesHQ + j]
						    += (*Mass)(j, k) * ((*dAcceleration)(i* nnodesHQ + k) + a_n[nodes[k] + NodeShift[i]]);
		}

// RHS->Write();
#if !defined(USE_PETSC) // && !defined(other parallel libs)//06.2013. WW
		for (std::size_t i = 0; i < dim; i++)
			for (int j = 0; j < nnodesHQ; j++)
				b_rhs[eqs_number[j] + NodeShift[i]] -= (*RHS)[i * nnodesHQ + j];

		if (excavation)
		{
			for (int i = 0; i < nnodesHQ; i++)
			{
				if (!onExBoundaryState[i])
				{
					for (std::size_t j = 0; j < dim; j++)
						b_rhs[eqs_number[i] + NodeShift[j]] = 0.0;
				}
			}
		}
#endif
	}
	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec:: LocalAssembly_continuum()
	   Aufgabe:
	           Compute the local finite element matrices with the framework
	        of continuum assumption
	   Formalparameter:
	           E:
	         const int update  : indicator to update stress and strain only

	   Programming:
	   02/2000   OK   Elastische Elemente
	   09/2002   WW   Local assemby of stiffness matrix of elasto-plastic
	   tangential model
	   Local assemby of residual
	   07/2003   WW   Quadratic triangle element
	   06/2004   WW   Generalize for different element types as a member of class
	   12/2005   WW   Creep
	 **************************************************************************/
	void CFiniteElementVec::LocalAssembly_continuum(const int update)
	{
		long i;

		Matrix* p_D = NULL;
		eleV_DM = ele_value_dm[MeshElement->GetIndex()];

		// ---- Gauss integral
		int gp_r = 0, gp_s = 0, gp_t;
		gp = 0;
		gp_t = 0;
		double fkt = 0.0;

		// WW double *DevStress ;
		const int PModel = smat->Plasticity_type;
		double dPhi = 0.0; // Sclar factor for the plastic strain
		//  double J2=0.0;
		double dS = 0.0;

		double ThermalExpansion = 0.0;
		double t1 = 0.0;
		bool Strain_TCS = false;

#ifdef JFNK_H2M
		bool JFNK = false;
		if (pcs->m_num->nls_method == 2) // 10.09.2010. WW
			JFNK = true;
#endif
		//
		ThermalExpansion = 0.0;
		// Thermal effect
		if (smat->Thermal_Expansion() > 0.0)
			ThermalExpansion = smat->Thermal_Expansion();

		// Get porosity model
		// ---- Material properties
		// For swelling pressure;
		double deporo = 0.0;
		// OK_MMP
		//--------------------------------------------------------------------
		// MMP medium properties
		int PoroModel = m_mmp->porosity_model;
		if (PoroModel == 4)
			// OK deporo =  PCSGetElementPorosityChangeRate(index)/(double)ele_dim;
			// MX  deporo = PCSGetELEValue(index,NULL,1.0,"n_sw_Rate")/(double)ele_dim;
			deporo = h_pcs->GetElementValue(Index, h_pcs->GetElementValueIndex("n_sw_rate")) / (double)ele_dim;
		if (T_Flag)
			for (i = 0; i < nnodes; i++)
			{
				T1[i] = t_pcs->GetNodeValue(nodes[i], idx_T1);
				Temp[i] = t_pcs->GetNodeValue(nodes[i], idx_T1) - t_pcs->GetNodeValue(nodes[i], idx_T0);
			}
		//

		if (PModel == 1 || PModel == 10 || PModel == 11) // WX:modified for DP with tension cutoff
			smat->CalulateCoefficent_DP();
		//
		if (PModel != 3 && smat->Youngs_mode != 2) // modified due to transverse isotropic elasticity: UJG 24.11.2009
		{
			if (smat->Youngs_mode < 10 || smat->Youngs_mode > 13)
			{
#ifdef RFW_FRACTURE
				smat->Calculate_Lame_Constant(GetMeshElement());
#else
			smat->Calculate_Lame_Constant();
#endif
				//
				smat->ElasticConsitutive(ele_dim, De);
			}
			else
				*De = *(smat->getD_tran()); // UJG/WW
		}
		// WX: 06.2012 E depends on stress, strain ...
		if (smat->E_Function_Model > 0)
		{
			double tmp_value = 1;
			tmp_value = smat->E_Function(ele_dim, eleV_DM, nGaussPoints);
			*De *= tmp_value;
		}

		if (PModel == 5)
			smat->CalculateCoefficent_HOEKBROWN(); // WX:02.2011
		//
		if (PoroModel == 4 || T_Flag || smat->Creep_mode > 0)
			Strain_TCS = true;
		//
		if (smat->CreepModel() == 1000) // HL_ODS
			smat->CleanTrBuffer_HL_ODS();
		// Loop over Gauss points
		for (gp = 0; gp < nGaussPoints; gp++)
		{
			//---------------------------------------------------------
			//  Get local coordinates and weights
			//  Compute Jacobian matrix and its determinate
			//---------------------------------------------------------
			fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

			//---------------------------------------------------------
			// Compute geometry
			//---------------------------------------------------------
			getGradShapefunctValues(gp, 2);
			getShapefunctValues(gp, 2);
			if (smat->Youngs_mode == 2) // WW/UJG. 22.01.2009
			{
				smat->CalcYoungs_SVV(CalcStrain_v());
				smat->ElasticConsitutive(ele_dim, De);
			}

			ComputeStrain(gp);
			if (update)
				RecordGuassStrain(gp, gp_r, gp_s, gp_t);
			if (F_Flag || T_Flag)
				getShapefunctValues(gp, 1); // Linear order interpolation function
			//---------------------------------------------------------
			// Material properties (Integration of the stress)
			//---------------------------------------------------------
			// Initial the stress vector
			if (PModel != 3)
			{
				for (i = 0; i < ns; i++)
					dstress[i] = 0.0;
				if (!excavation) // WX:07.2011 nonlinear excavation
				{
					// De->Write();
					De->multi(dstrain, dstress);
					if (smat->Time_Dependent_E_nv_mode > MKleinsteZahl && pcs->ExcavMaterialGroup < 0)
						for (i = 0; i < ns; i++)
							dstress[i] -= (*eleV_DM->Stress)(i, gp) - (*eleV_DM->Stress0)(i, gp);
				}
			}

			//---------------------------------------------------------
			// Integrate the stress by return mapping:
			//---------------------------------------------------------
			if (excavation) // if elem is excavated, only do the comp. rhs WX:08.2011
			{
				for (i = 0; i < ns; i++)
					dstress[i] += (*eleV_DM->Stress)(i, gp);
				S_Water = 1.0; // WX:02.2013
				if (Flow_Type > 0 && Flow_Type != 10)
					S_Water = interpolate(_nodal_S, 1);
			}
			else
			{
				switch (PModel)
				{
					case -1: // Pure elasticity
						// Non-linear elasticity: TE model. 10.03.2008. WW
						for (i = 0; i < ns; i++)
							dstress[i] += (*eleV_DM->Stress)(i, gp);
						break;
					case 1: // Drucker-Prager model
#ifdef JFNK_H2M
						if (smat->StressIntegrationDP(gp, eleV_DM, dstress, dPhi, update) && !JFNK)
#else
					if (smat->StressIntegrationDP(gp, eleV_DM, dstress, dPhi, update))
#endif

							// WW DevStress = smat->devS;
							smat->ConsistentTangentialDP(ConsistDep, dPhi, ele_dim);
						break;
					case 10: // Drucker-Prager model, direct integration. 02/06 WW
						if (smat->DirectStressIntegrationDP(gp, eleV_DM, dstress, update))
						{
							*ConsistDep = *De;
							smat->TangentialDP(ConsistDep);
							dPhi = 1.0;
						}
						break;
					case 11: // WX: 08.2010
					{
						double mm = 0.; // WX:09.2010. for DP with Tension.
						switch (smat->DirectStressIntegrationDPwithTension(gp, De, eleV_DM, dstress, update, mm))
						{
							case 1:
							{
								*ConsistDep = *De;
								smat->TangentialDP2(ConsistDep);
								dPhi = 1.0;
							}
							break;
							case 2:
							{
								*ConsistDep = *De;
								smat->TangentialDPwithTension(ConsistDep, mm);
								dPhi = 1.0;
							}
							break;
							case 3:
							{
								*ConsistDep = *De;
								smat->TangentialDPwithTensionCorner(ConsistDep, mm);
								dPhi = 1.0;
							}
							break;
							default:
								break;
						}
						break;
					}
					case 2: // Rotational hardening model
						// Compute stesses and plastic multi-plier
						dPhi = 0.0;
						if (smat->CalStress_and_TangentialMatrix_SYS(gp, eleV_DM, De, ConsistDep, dstress, update) > 0)
							dPhi = 1.0;
						break;
					case 3: // Generalized Cam-Clay model
						for (i = 0; i < ns; i++)
							dstress[i] = dstrain[i];
						//
						if (smat->SwellingPressureType == 3)
						{
							double suc = interpolate(_nodal_cp0);
							double dsuc = interpolate(_nodal_dcp);
							(*smat->data_Youngs)(7) = suc;
							(*smat->data_Youngs)(8) = dsuc;
							smat->CalStress_and_TangentialMatrix_CC_SubStep(gp, eleV_DM, dstress, ConsistDep, update);
							//
							// double pn = -((*eleV_DM->Stress)(0, gp)+(*eleV_DM->Stress)(1, gp)+
							//             (*eleV_DM->Stress)(2, gp))/3.0;
							// de_vsw = -smat->TEPSwellingParameter(pn,suc)*dsuc/(suc+1.0e5);
						}
						/*
						   else if (smat->SwellingPressureType==4)  //Power. 07.05.2008 WW
						   {
						   double suc = interpolate(AuxNodal1);
						   double dsuc = interpolate(AuxNodal);
						   smat->TEPSwellingParameter_kis(suc);
						   S_Water = m_mmp->SaturationCapillaryPressureFunction(suc);
						   dS = S_Water - m_mmp->SaturationCapillaryPressureFunction(suc-dsuc);
						   de_vsw = pow(S_Water, (*smat->data_Youngs)(2) )*dS;
						   }
						 */
						else
							smat->CalStress_and_TangentialMatrix_CC(gp, eleV_DM, dstress, ConsistDep, update);

						dPhi = 1.0;
						break;
					case 4: // Mohr-Coloumb	//WX:10.2010
						if (smat->DirectStressIntegrationMOHR(gp, eleV_DM, dstress, update, De, pcs->ite_steps))
						{
							*ConsistDep = *De;
							// also for tension
							smat->TangentialMohrShear(ConsistDep);
							// ConsistDep->Write();
							dPhi = 1.0;
						}
						break;
					case 44: // WX:12.2011 Mohr-Coloumb	bedding
						if (smat->StressIntegrationMOHR_Aniso(gp, eleV_DM, dstress, update, De))
						{
							*ConsistDep = *De;
							smat->TangentialMohrShear(ConsistDep); // also for tension
							// ConsistDep->Write();
							dPhi = 1.0;
						}
						break;
						/*case 5:
						   if(smat->StressIntegrationHoekBrown(gp, eleV_DM, dstress, update, De))
						   {
						   *ConsistDep = *De;
						         smat->TangentialHoekBrown(ConsistDep);		//also for tension
						         //ConsistDep->Write();
						         dPhi = 1.0;
						   }
						   break;*/
				}
			}
			// --------------------------------------------------------------------
			// Stress increment by heat, swelling, or heat
			//
			if (Strain_TCS)
			{
				if (PModel == 3)
					smat->ElasticConsitutive(ele_dim, De);
				for (i = 0; i < ns; i++)
					strain_ne[i] = 0.0;
				if (PoroModel == 4) // For swelling pressure

					for (i = 0; i < 3; i++)
						strain_ne[i] -= deporo;
				//
				if (T_Flag) // Contribution by thermal expansion
				{
					Tem = 0.0;
					t1 = 0.0;
					for (i = 0; i < nnodes; i++)
					{
						Tem += shapefct[i] * Temp[i];
						t1 += shapefct[i] * T1[i];
					}
					for (i = 0; i < 3; i++) // JT: This was commented. SHOULDN'T BE!
						strain_ne[i] -= ThermalExpansion * Tem;
				}
				// Strain increment by creep
				if (smat->Creep_mode == 1 || smat->Creep_mode == 2 || smat->Creep_mode == 3 || smat->Creep_mode == 4)
				// TN:add BGRb BGRsf,WX
				{
					for (i = 0; i < ns; i++)
						stress_ne[i] = (*eleV_DM->Stress)(i, gp);
					smat->AddStain_by_Creep(ns, stress_ne, strain_ne, t1);
				}
				if (smat->Creep_mode == 1000) // HL_ODS. Strain increment by creep
				{
					for (i = 0; i < ns; i++)
						stress_ne[i] = (*eleV_DM->Stress)(i, gp);
					smat->AddStain_by_HL_ODS(eleV_DM, stress_ne, strain_ne, t1);
				}
				// Stress deduced by thermal or swelling strain incremental:
				De->multi(strain_ne, dstress);
				for (i = 0; i < ns; i++) // JT: This was commented. It shouldn't be.
					dstrain[i] += strain_ne[i];
			}

			if (smat->Creep_mode == 1001) // BURGERS.
			{
				// get total strains at time t and current iteration i
				std::vector<double> strain_t(6), strain_curr(6);

				// get stresses as well as internal variables at time t and current iteration i (equal before local
				// Newton)
				std::vector<double> stress_curr(6), eps_K_curr(6), eps_M_curr(6);
				for (int compnt(0); compnt < ns; compnt++)
				{
					strain_t[compnt] = (*eleV_DM->Strain_t_ip)(compnt, gp);
					strain_curr[compnt] = strain_t[compnt] + dstrain[compnt];
					stress_curr[compnt] = (*eleV_DM->Stress)(compnt, gp);
					eps_K_curr[compnt] = (*eleV_DM->Strain_Kel)(compnt, gp);
					eps_M_curr[compnt] = (*eleV_DM->Strain_Max)(compnt, gp);
				}
				// Fill remaining 3D components with 0 if only 2D
				if (ns == 4)
					for (int compnt(4); compnt < 6; compnt++)
					{
						strain_t[compnt] = 0.0;
						strain_curr[compnt] = 0.0;
						stress_curr[compnt] = 0.0;
						eps_K_curr[compnt] = 0.0;
						eps_M_curr[compnt] = 0.0;
					}

				// 6x6 tangent
				Matrix ConsD(6, 6);
				// Pass as 6D vectors, i.e. set stress and strain [4] and [5] to zero for 2D and AXI as well as
				// strain[3] to zero for 2D (plane strain)
				double local_res;
				smat->LocalNewtonBurgers(dt, strain_curr, stress_curr, eps_K_curr, eps_M_curr, ConsD, t1,
										 local_res);

				// Then update (and reduce for 2D) stress increment vector and reduce (for 2D) ConsistDep, update
				// internal variables
				for (int compnt(0); compnt < ns; compnt++)
				{
					dstress[compnt] = stress_curr[compnt];
					if (update > 0)
					{
						(*eleV_DM->Strain_Kel)(compnt, gp) = eps_K_curr[compnt];
						(*eleV_DM->Strain_Max)(compnt, gp) = eps_M_curr[compnt];
						(*eleV_DM->Strain_t_ip)(compnt, gp) = strain_curr[compnt];
						(*eleV_DM->ev_loc_nr_res)(gp) = local_res;
					}
					for (int compnt2(0); compnt2 < ns; compnt2++)
						(*De)(compnt, compnt2) = ConsD(compnt, compnt2);
					(*eleV_DM->Strain)(compnt, gp) = strain_curr[compnt];
				}
			}

			if (smat->Creep_mode == 1002) // MINKLEY
			{
				// get total strains at time t and current iteration i
				std::vector<double> strain_t(6), strain_curr(6);

				std::vector<double> stress_curr(6), eps_K_curr(6), eps_M_curr(6), eps_pl_curr(6);
				for (int compnt(0); compnt < ns; compnt++)
				{
					stress_curr[compnt] = (*eleV_DM->Stress)(compnt, gp);
					eps_K_curr[compnt] = (*eleV_DM->Strain_Kel)(compnt, gp);
					eps_M_curr[compnt] = (*eleV_DM->Strain_Max)(compnt, gp);
					eps_pl_curr[compnt] = (*eleV_DM->Strain_pl)(compnt, gp);
					strain_t[compnt] = (*eleV_DM->Strain_t_ip)(compnt, gp);
					strain_curr[compnt] = strain_t[compnt] + dstrain[compnt];
				}
				// Fill remaining 3D components with 0 if only 2D
				if (ns == 4)
					for (int compnt(4); compnt < 6; compnt++)
					{
						strain_t[compnt] = 0.0;
						strain_curr[compnt] = 0.0;
						stress_curr[compnt] = 0.0;
						eps_K_curr[compnt] = 0.0;
						eps_M_curr[compnt] = 0.0;
						eps_pl_curr[compnt] = 0.0;
					}

				double e_pl_v = (*eleV_DM->e_pl)(gp);
				double e_pl_eff = (*eleV_DM->pStrain)(gp);
				double lam = (*eleV_DM->lambda_pl)(gp);//NOTE: May set starting value to zero in case of trouble with load reversals

				// 6x6 tangent
				Matrix ConsD(6, 6);

				// Pass as 6D vectors, i.e. set stress and strain [4] and [5] to zero for 2D and AXI as well as
				// strain[3] to zero for 2D (plane strain)
				double local_res(0.);
				smat->LocalNewtonMinkley(dt, strain_curr, stress_curr, eps_K_curr, eps_M_curr, eps_pl_curr, e_pl_v,
										 e_pl_eff, lam, ConsD, t1, local_res);

				// Then update (and reduce for 2D) stress increment vector and reduce (for 2D) ConsistDep, update
				// internal variables
				for (int compnt(0); compnt < ns; compnt++)
				{
					dstress[compnt] = stress_curr[compnt];
					if (update > 0)
					{
						(*eleV_DM->Strain_Kel)(compnt, gp) = eps_K_curr[compnt];
						(*eleV_DM->Strain_Max)(compnt, gp) = eps_M_curr[compnt];
						(*eleV_DM->Strain_pl)(compnt, gp) = eps_pl_curr[compnt];
						(*eleV_DM->Strain_t_ip)(compnt, gp) = strain_curr[compnt];
					}
					for (int compnt2(0); compnt2 < ns; compnt2++)
						(*De)(compnt, compnt2) = ConsD(compnt, compnt2);
					(*eleV_DM->Strain)(compnt, gp) = strain_curr[compnt];
				}

				if (update > 0)
				{
					(*eleV_DM->e_pl)(gp) = e_pl_v;
					(*eleV_DM->lambda_pl)(gp) = lam;
					(*eleV_DM->pStrain)(gp) = e_pl_eff;
					(*eleV_DM->ev_loc_nr_res)(gp) = local_res;
				}
			}


			// Fluid coupling;
			S_Water = 1.0;
			if (Flow_Type > 0 && Flow_Type != 10)
				S_Water = interpolate(_nodal_S, 1);
			// Decovalex. Swelling pressure
			if (smat->SwellingPressureType == 1)
			{
				dS = -interpolate(_nodal_S0, 1);
				dS += S_Water;
				for (i = 0; i < 3; i++)
					dstress[i] -= 2.0 * S_Water * dS * smat->Max_SwellingPressure;
			}
			else if (smat->SwellingPressureType == 2) // LBNL's model
			{
				dS = -interpolate(_nodal_S0, 1);
				dS += S_Water;
				for (i = 0; i < 3; i++)
					dstress[i] -= dS * smat->Max_SwellingPressure;
			}
			/*
			   else if(smat->SwellingPressureType==3||smat->SwellingPressureType==4) // TEP model
			   {
			   for (i = 0; i < 3; i++)
			      strain_ne[i] = -de_vsw;
			   for (i = 3; i < ns; i++)
			      strain_ne[i] = 0.;
			   smat->ElasticConsitutive(ele_dim, De);
			   De->multi(strain_ne, dstress);
			   }
			 */
			// Assemble matrices and RHS
			if (update < 1)
			{
				//---------------------------------------------------------
				// Assemble matrices and RHS
				//---------------------------------------------------------
				if (dPhi <= 0.0)
					p_D = De;
				else
					p_D = ConsistDep;
				for (i = 0; i < ns; i++)
					stress0[i] = (*eleV_DM->Stress0)(i, gp);
				ComputeMatrix_RHS(fkt, p_D);
			}
			else // Update stress

				for (i = 0; i < ns; i++)
					(*eleV_DM->Stress)(i, gp) = dstress[i];
		}
		// The mapping of Gauss point strain to element nodes
		if (update)
			ExtropolateGuassStrain();
		else if (smat->Creep_mode == 1000) // HL_ODS. Strain increment by creep
			smat->AccumulateEtr_HL_ODS(eleV_DM, nGaussPoints);
	}

	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec::RecordGuassValues()
	   Aufgabe:
	           Accumulate stress at each nodes
	   Formalparameter:
	           E:

	   Programming:
	   06/2004   WW
	 **************************************************************************/
	bool CFiniteElementVec::RecordGuassStrain(const int gp, const int gp_r, const int gp_s, int gp_t)
	{
		int LoIndex = 0;

		//---------------------------------------------------------
		// Accumulate strains
		//---------------------------------------------------------
		switch (MeshElement->GetElementType())
		{
			case MshElemType::QUAD: // Quadralateral
				LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
				Sxx[LoIndex] = dstrain[0];
				Syy[LoIndex] = dstrain[1];
				Sxy[LoIndex] = dstrain[3];
				Szz[LoIndex] = dstrain[2];
				break;
			case MshElemType::TRIANGLE: // Triangle
				Sxx[gp] = dstrain[0];
				Syy[gp] = dstrain[1];
				Szz[gp] = dstrain[2];
				Sxy[gp] = dstrain[3];
				break;
			case MshElemType::HEXAHEDRON: // Hexahedra
				LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
				if (LoIndex < 0)
					return false;
				Sxx[LoIndex] = dstrain[0];
				Syy[LoIndex] = dstrain[1];
				Szz[LoIndex] = dstrain[2];
				Sxy[LoIndex] = dstrain[3];
				Sxz[LoIndex] = dstrain[4];
				Syz[LoIndex] = dstrain[5];
				break;
			case MshElemType::TETRAHEDRON: // Tedrahedra
				Sxx[gp] = dstrain[0];
				Syy[gp] = dstrain[1];
				Szz[gp] = dstrain[2];
				Sxy[gp] = dstrain[3];
				Sxz[gp] = dstrain[4];
				Syz[gp] = dstrain[5];
				break;
			case MshElemType::PYRAMID:
				Sxx[gp] = dstrain[0];
				Syy[gp] = dstrain[1];
				Szz[gp] = dstrain[2];
				Sxy[gp] = dstrain[3];
				Sxz[gp] = dstrain[4];
				Syz[gp] = dstrain[5];
				break;
			default:
				break;
				// 3D
		}
		return false;
	}

	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec::ExtropolateGuassStrain()
	   Aufgabe:
	           Extropolate the Gauss point strains to nodes
	   Formalparameter:
	           E:

	   Programming:
	   06/2004   WW
	   02/2007   Make it work for all 2nd variables
	 **************************************************************************/
	void CFiniteElementVec::ExtropolateGuassStrain()
	{
		// WX:03.2012. if excavation dbuff changed
		if (pcs->ExcavMaterialGroup > -1)
		{
			int tmp_excavstate = -1;
			for (int i = 0; i < nnodes; i++)
			{
				for (size_t jj = 0; jj < MeshElement->nodes[i]->getConnectedElementIDs().size(); jj++)
				{
					tmp_excavstate
					    = pcs->m_msh->ele_vector[MeshElement->nodes[i]->getConnectedElementIDs()[jj]]->GetExcavState();
					if (tmp_excavstate > -1)
						dbuff[i] -= 1;
				}
				if (dbuff[i] < MKleinsteZahl) // avoid error
					dbuff[i] = 1;
			}
		}

		// l1=l2=l3=l4=0;
		MshElemType::type ElementType = MeshElement->GetElementType();
		CalcXi_p();

		//
		int i_s, i_e, ish;
		i_s = 0;
		i_e = nnodes;
		ish = 0;
		if (ElementType == MshElemType::TETRAHEDRON) // tet
		{
			i_s = 1;
			i_e = nnodes + 1;
			ish = 1;
		}
		//---------------------------------------------------------
		// Mapping Gauss point strains to nodes and update nodes
		// strains:
		//---------------------------------------------------------
		double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz;
		double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz;
		avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = 0.0;
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			avgESxx = CalcAverageGaussPointValues(Sxx);
			avgESyy = CalcAverageGaussPointValues(Syy);
			avgESzz = CalcAverageGaussPointValues(Szz);
			avgESxy = CalcAverageGaussPointValues(Sxy);
			avgESxz = CalcAverageGaussPointValues(Sxz);
			avgESyz = CalcAverageGaussPointValues(Syz);
		}

		ConfigShapefunction(ElementType);
		for (int i = 0; i < nnodes; i++)
		{
			ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = 0.0;

			// Calculate values at nodes
			if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
			{
				SetExtropoGaussPoints(i);
				ComputeShapefct(1, dbuff0); // Linear interpolation function
				//
				for (int j = i_s; j < i_e; j++)
				{
					const int k = j - ish;
					ESxx += Sxx[j] * dbuff0[k];
					ESyy += Syy[j] * dbuff0[k];
					ESxy += Sxy[j] * dbuff0[k];
					ESzz += Szz[j] * dbuff0[k];
					if (ele_dim == 3)
					{
						ESxz += Sxz[j] * dbuff0[k];
						ESyz += Syz[j] * dbuff0[k];
					}
				}
			}
			else if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
			{
				// average
				ESxx = avgESxx;
				ESyy = avgESyy;
				ESxy = avgESxy;
				ESzz = avgESzz;
				if (ele_dim == 3)
				{
					ESxz = avgESxz;
					ESyz = avgESyz;
				}
			}

			// Average value of the contribution of all neighbor elements
			ESxx /= dbuff[i];
			ESyy /= dbuff[i];
			ESxy /= dbuff[i];
			ESzz /= dbuff[i];

			ESxx += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[0]);
			ESyy += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[1]);
			ESzz += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[2]);
			ESxy += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[3]);

			pcs->SetNodeValue(nodes[i], _nodal_strain_indices[0], ESxx);
			pcs->SetNodeValue(nodes[i], _nodal_strain_indices[1], ESyy);
			pcs->SetNodeValue(nodes[i], _nodal_strain_indices[2], ESzz);
			pcs->SetNodeValue(nodes[i], _nodal_strain_indices[3], ESxy);

			if (ele_dim == 3)
			{
				ESxz /= dbuff[i];
				ESyz /= dbuff[i];
				ESxz += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[4]);
				ESyz += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[5]);
				//
				pcs->SetNodeValue(nodes[i], _nodal_strain_indices[4], ESxz);
				pcs->SetNodeValue(nodes[i], _nodal_strain_indices[5], ESyz);
			}
		}
	}

	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec::ExtropolateGuassStress()
	   Aufgabe:
	           Extropolate the Gauss point strains to nodes
	   Formalparameter:
	           E:

	   Programming:
	   06/2004   WW
	   03/2007   WW  Generize for all 2nd variables
	 **************************************************************************/
	void CFiniteElementVec::ExtropolateGuassStress()
	{
		// For strain and stress extropolation all element types
		// Number of elements associated to nodes
		nnodes = MeshElement->nnodes;
		// Node indices
		for (int i = 0; i < nnodes; i++)
		{
			nodes[i] = MeshElement->nodes[i]->GetIndex();
			dbuff[i] = (double)MeshElement->nodes[i]->getConnectedElementIDs().size();
		}
		//
		eleV_DM = ele_value_dm[MeshElement->GetIndex()];
		if (eleV_DM->pStrain) // 08.02.2008 WW
			idx_pls = pcs->GetNodeValueIndex("STRAIN_PLS");
		//
		MshElemType::type ElementType = MeshElement->GetElementType();
		SetIntegrationPointNumber(ElementType);
		for (gp = 0; gp < nGaussPoints; gp++)
		{
			int gp_r, gp_s, gp_t;
			gp_r = gp_s = gp_t = 0;
			SetGaussPoint(gp, gp_r, gp_s, gp_t);
			int i = gp;
			if (ElementType == MshElemType::QUAD || ElementType == MshElemType::HEXAHEDRON)
			{
				i = GetLocalIndex(gp_r, gp_s, gp_t);
				if (i == -1)
					continue;
			}

			Sxx[i] = (*eleV_DM->Stress)(0, gp);
			Syy[i] = (*eleV_DM->Stress)(1, gp);
			Szz[i] = (*eleV_DM->Stress)(2, gp);
			Sxy[i] = (*eleV_DM->Stress)(3, gp);
			if (eleV_DM->pStrain)
				pstr[i] = (*eleV_DM->pStrain)(gp);
			else
				pstr[i] = 0.0; // 08.02.2008 WW
			if (ele_dim == 3)
			{
				Sxz[i] = (*eleV_DM->Stress)(4, gp);
				Syz[i] = (*eleV_DM->Stress)(5, gp);
			}
		}
		//
		CalcXi_p();

		//
		int i_s, i_e, ish;
		i_s = 0;
		i_e = nnodes;
		ish = 0;
		if (ElementType == MshElemType::TETRAHEDRON) // tet
		{
			i_s = 1;
			i_e = nnodes + 1;
			ish = 1;
		}
		//---------------------------------------------------------
		// Mapping Gauss point strains to nodes and update nodes
		// strains:
		//---------------------------------------------------------
		double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz, Pls;
		double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz, avgPls;
		avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = avgPls = 0.0;
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			avgESxx = CalcAverageGaussPointValues(Sxx);
			avgESyy = CalcAverageGaussPointValues(Syy);
			avgESzz = CalcAverageGaussPointValues(Szz);
			avgESxy = CalcAverageGaussPointValues(Sxy);
			avgESxz = CalcAverageGaussPointValues(Sxz);
			avgESyz = CalcAverageGaussPointValues(Syz);
			avgPls = CalcAverageGaussPointValues(pstr);
		}

		ConfigShapefunction(ElementType);
		for (int i = 0; i < nnodes; i++)
		{
			ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = Pls = 0.0;

			// Calculate values at nodes
			if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
			{
				//
				SetExtropoGaussPoints(i);
				//
				ComputeShapefct(1, dbuff0); // Linear interpolation function
				//
				for (int j = i_s; j < i_e; j++)
				{
					int k = j - ish;
					ESxx += Sxx[j] * dbuff0[k];
					ESyy += Syy[j] * dbuff0[k];
					ESxy += Sxy[j] * dbuff0[k];
					ESzz += Szz[j] * dbuff0[k];
					Pls += pstr[j] * dbuff0[k];
					if (ele_dim == 3)
					{
						ESxz += Sxz[j] * dbuff0[k];
						ESyz += Syz[j] * dbuff0[k];
					}
				}
			}
			else if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
			{
				// average
				ESxx = avgESxx;
				ESyy = avgESyy;
				ESxy = avgESxy;
				ESzz = avgESzz;
				Pls = avgPls;
				if (ele_dim == 3)
				{
					ESxz = avgESxz;
					ESyz = avgESyz;
				}
			}

			// Average value of the contribution of ell neighbor elements
			ESxx /= dbuff[i];
			ESyy /= dbuff[i];
			ESxy /= dbuff[i];
			ESzz /= dbuff[i];
			Pls /= dbuff[i];
			//
			long node_i = nodes[i];
			ESxx += pcs->GetNodeValue(node_i, _nodal_stress_indices[0]);
			ESyy += pcs->GetNodeValue(node_i, _nodal_stress_indices[1]);
			ESzz += pcs->GetNodeValue(node_i, _nodal_stress_indices[2]);
			ESxy += pcs->GetNodeValue(node_i, _nodal_stress_indices[3]);
			if (eleV_DM->pStrain) // 08.02.2008 WW
				Pls += pcs->GetNodeValue(node_i, idx_pls);

			pcs->SetNodeValue(node_i, _nodal_stress_indices[0], ESxx);
			pcs->SetNodeValue(node_i, _nodal_stress_indices[1], ESyy);
			pcs->SetNodeValue(node_i, _nodal_stress_indices[2], ESzz);
			pcs->SetNodeValue(node_i, _nodal_stress_indices[3], ESxy);
			if (eleV_DM->pStrain) // 08.02.2008 WW
				pcs->SetNodeValue(node_i, idx_pls, fabs(Pls));

			if (ele_dim == 3)
			{
				ESxz /= dbuff[i];
				ESyz /= dbuff[i];

				ESxz += pcs->GetNodeValue(node_i, _nodal_stress_indices[4]);
				ESyz += pcs->GetNodeValue(node_i, _nodal_stress_indices[5]);

				pcs->SetNodeValue(node_i, _nodal_stress_indices[4], ESxz);
				pcs->SetNodeValue(node_i, _nodal_stress_indices[5], ESyz);
			}
		}
	}

	//==========================================================================
	// Enhanced strain element
	/***************************************************************************
	   GeoSys - Funktion:
	        CFiniteElementVec::CheckNodesInJumpedDomain(const double *tangJump)

	   Aufgabe:
	         Compute the regular enhanced strain matrix (Only 2D)
	          Ge: (See the related references)
	      (Prorior: element nodes is fetched )

	   Formalparameter:
	         E:
	   const double *tangJump   : Tangential to the jump surface

	   Programming:
	   06/2004     WW        Erste Version
	 **************************************************************************/
	void CFiniteElementVec::CheckNodesInJumpedDomain()
	{
		int i;
		double cdotpdt;
		//	static double dphi_e[3];

		// 2D
		// Get the center of the discontinuity

		X0[0] = 0.5 * ((*eleV_DM->NodesOnPath)(0, 0) + (*eleV_DM->NodesOnPath)(0, 1));
		X0[1] = 0.5 * ((*eleV_DM->NodesOnPath)(1, 0) + (*eleV_DM->NodesOnPath)(1, 1));
		X0[2] = 0.5 * ((*eleV_DM->NodesOnPath)(2, 0) + (*eleV_DM->NodesOnPath)(2, 1));

		// Determine nodes in the jumping part
		for (i = 0; i < nnodesHQ; i++)
		{
			NodesInJumpedA[i] = false;
			cdotpdt = 0.0;
			cdotpdt += n_jump[0] * (X[i] - X0[0]);
			cdotpdt += n_jump[1] * (Y[i] - X0[1]);
			if (ele_dim == 3)
				cdotpdt += n_jump[2] * (Z[i] - X0[2]);
			if (cdotpdt > 0.0)
				NodesInJumpedA[i] = true; // Nodes in $\Omega_+$
		}
	}

	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec::ComputeRESM(double * normJump)

	   Aufgabe:
	         Compute the regular enhanced strain matrix (Only 2D)
	          Ge: (See the related references)
	      (Prorior: element nodes is fetched )

	   Formalparameter:
	         E:
	   const double *tangJump   : Tangential to the jump surface

	   Programming:
	   06/2004     WW        Erste Version
	 **************************************************************************/
	void CFiniteElementVec::ComputeRESM(const double* tangJump)
	{
		static double dphi_e[3];

		for (size_t i = 0; i < ele_dim; i++)
		{
			dphi_e[i] = 0.0;
			for (int j = 0; j < nnodesHQ; j++)
				//        for(int j=0; j<nnodes; j++)

				if (NodesInJumpedA[j])
					dphi_e[i] += dshapefctHQ[i * nnodesHQ + j];
			//               dphi_e[i] += dshapefct[i*nnodes+j];
		}
		// !!! Only for 2D up to now
		tangJump = tangJump;
		// Column 1
		(*Ge)(0, 0) = n_jump[0] * dphi_e[0];
		(*Ge)(1, 0) = n_jump[1] * dphi_e[1];
		(*Ge)(2, 0) = 0.0;
		(*Ge)(3, 0) = n_jump[0] * dphi_e[1] + n_jump[1] * dphi_e[0];

		// Column 2
		(*Ge)(0, 1) = -n_jump[1] * dphi_e[0];
		(*Ge)(1, 1) = n_jump[0] * dphi_e[1];
		(*Ge)(2, 1) = 0.0;
		(*Ge)(3, 1) = -n_jump[1] * dphi_e[1] + n_jump[0] * dphi_e[0];
	}

	/***************************************************************************
	   GeoSys - Funktion:
	           CFiniteElementVec::ComputeSESM(double * normJump)

	   Aufgabe:
	         Compute the singular enhanced strain matrix (Only 2D)
	          Ge: (See the related references)
	      (Prorior: element nodes is fetched )

	   Formalparameter:
	         E:
	   const double *tangJump   : Tangential to the jump surface

	   Programming:
	   06/2004     WW        Erste Version
	 **************************************************************************/
	void CFiniteElementVec::ComputeSESM(const double* tangJump)
	{
		// !!! Only for 2D up to now
		tangJump = tangJump;
		// Column 1
		(*Pe)(0, 0) = n_jump[0] * n_jump[0];
		(*Pe)(0, 1) = n_jump[1] * n_jump[1];
		(*Pe)(0, 2) = 0.0;
		(*Pe)(0, 3) = 2.0 * n_jump[0] * n_jump[1];

		// Column 2
		(*Pe)(1, 0) = -n_jump[0] * n_jump[1];
		(*Pe)(1, 1) = n_jump[1] * n_jump[0];
		(*Pe)(1, 2) = 0.0;
		(*Pe)(1, 3) = n_jump[0] * n_jump[0] - n_jump[1] * n_jump[1];
	}

	/***************************************************************************
	   GeoSys - Funktion:
	      CFiniteElementVec::ComputePrincipleStresses(double *dEStress)
	      (2D only)
	   Aufgabe:

	   Formalparameter:
	         E:
	   const double *Stresses: Stresses

	   Return: Angle of maxium principle stress component to x direction
	   sig2<sig1

	   Programming:
	   06/2004     WW        Erste Version
	 **************************************************************************/
	double CFiniteElementVec::ComputePrincipleStresses(const double* Stresses)
	{
		double prin_ang, var;
		// Angle of the principle plane
		if (fabs(Stresses[0] - Stresses[1]) < MKleinsteZahl)
			return 0.0;
		prin_ang = atan(2.0 * Stresses[3] / (Stresses[0] - Stresses[1]));
		// Principle stress 1
		pr_stress[0] = 0.5 * (Stresses[0] + Stresses[1]) + 0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang)
		               + Stresses[3] * sin(prin_ang);
		// Principle stress 2
		pr_stress[1] = 0.5 * (Stresses[0] + Stresses[1]) - 0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang)
		               - Stresses[3] * sin(prin_ang);
		pr_stress[2] = Stresses[2];

		prin_ang *= 0.5;

		// Jump direction
		if (pr_stress[1] >= pr_stress[0])
		{
			prin_ang += 0.5 * pai;
			var = pr_stress[0];
			pr_stress[0] = pr_stress[1];
			pr_stress[1] = var;
		}
		// if(pr_stress[0]<0.0) prin_ang += pai;
		return prin_ang;
	}
	/***************************************************************************
	   GeoSys - Funktion:
	      CFiniteElementVec::ComputeJumpDirectionAngle(const double *Stresses, double *Mat)
	      (2D Drucker-Prager plasticity only.
	       cf. K. Runesson, D. Peric and S. Sture,
	      Discontinuous bifucations of elasto-plastic solutions at plane stress
	      and plane strain, Int. J. Plasticity 7(1991) 2087-2105)
	   Aufgabe:

	   Formalparameter:
	         E:
	   const double *Stresses: Stresses
	   const double *Mat     : Material parameters
	   0, dilatancy
	   1, frictional
	   2, Poission ratio
	   Return: Angle of maxium principle stress component to x direction

	   Programming:
	   06/2004     WW        Erste Version
	 **************************************************************************/
	double CFiniteElementVec::ComputeJumpDirectionAngle(const double* Mat)
	{
		double NormS, c1, c2;
		DeviatoricStress(pr_stress);
		NormS = sqrt(TensorMutiplication2(pr_stress, pr_stress, ele_dim));
		c1 = pr_stress[1] + Mat[2] * pr_stress[2] + 0.5 * (1.0 + Mat[2]) * (Mat[0] + Mat[1]) * NormS;
		c2 = pr_stress[0] + Mat[2] * pr_stress[2] + 0.5 * (1.0 + Mat[2]) * (Mat[0] + Mat[1]) * NormS;

		if (c1 >= 0.0 || c2 <= 0.0)
			NormS = -1.0;
		else
			NormS = atan(sqrt(-c2 / c1));

		return NormS; // The angle return through NormS
	}

	/**************************************************************************
	   GeoSys - Funktion:
	    void CFiniteElementVec::LocalAssembly_CheckLocalization

	   Aufgabe:
	     Trace discontinuity surface and determine the normal direction to it
	     element-wisely. (Material related)
	    (Drucker-Prager model and 2D only)

	   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
	   E :

	   Programmaenderungen:
	   06/2004   WW  Erste Version
	**************************************************************************/
	bool CFiniteElementVec::LocalAssembly_CheckLocalization(CElem * MElement)
	{
		int i, j, k;
		int MatGroup;

		double ep, p, normXi, n2;

		// For enhanced strain element
		double h_loc, detA, h_tol = 1.0e-5;
		double pr_stress_ang, loc_ang;
		static double OriJ[2], Nj[2], Aac[4], Mat[3];
		bool LOCed = false;
		//
		MeshElement = MElement;
		eleV_DM = ele_value_dm[MeshElement->GetIndex()];

		p = 0.0;
		// Get the total effective plastic strain
		ep = 0.0;
		for (i = 0; i < nGaussPoints; i++)
			ep += (*eleV_DM->pStrain)(i);
		ep /= (double)nGaussPoints;

		if (ep > 0.0) // in yield status
		{
			if (!eleV_DM->Localized)
			{
#ifdef RFW_FRACTURE
				smat->Calculate_Lame_Constant(MeshElement);
#endif
#ifndef RFW_FRACTURE
				smat->Calculate_Lame_Constant();
#endif
				smat->CalulateCoefficent_DP();

				MatGroup = MeshElement->GetPatchIndex();
				smat = msp_vector[MatGroup];

				Mat[0] = smat->Al;
				Mat[1] = smat->Xi;
				Mat[2] = smat->Poisson_Ratio();

				// Compute the average stress of this element
				for (i = 0; i < ns; i++)
					dstress[i] = 0.0;
				for (i = 0; i < nGaussPoints; i++)
					for (j = 0; j < ns; j++)
						dstress[j] += (*eleV_DM->Stress)(j, i);
				for (i = 0; i < ns; i++)
					dstress[i] /= (double)nGaussPoints;

				// Get the converged stresses of the previous step
				//--- Compute the determinate of the acoustic tensor
				pr_stress_ang = ComputePrincipleStresses(dstress);
				loc_ang = ComputeJumpDirectionAngle(Mat);

				normXi = sqrt(TensorMutiplication2(pr_stress, pr_stress, ele_dim));

				// Compute the localization condition
				h_loc = pr_stress[2] / normXi + 0.5 * (smat->Al + smat->Xi)
				        - 0.5 * sqrt(2.0 / (1.0 - smat->Poisson_Ratio())) * (smat->Al - smat->Xi);
				detA = 1.0e8;

				// Two possible jump orientation, i.e. pr_stress_ang+/-loc_ang
				OriJ[0] = pr_stress_ang - loc_ang + 0.5 * pai;
				OriJ[1] = pr_stress_ang + loc_ang + 0.5 * pai;

				// Compute the acoustic matrix
				DeviatoricStress(dstress);
				normXi = sqrt(TensorMutiplication2(dstress, dstress, ele_dim));
				if (loc_ang > 0.0)
				{
					for (i = 0; i < ns; i++)
						dstress[i] /= normXi;
					for (k = 0; k < 2; k++)
					{
						n_jump[0] = cos(OriJ[k]);
						n_jump[1] = sin(OriJ[k]);
						//
						Nj[0] = dstress[0] * n_jump[0] + dstress[3] * n_jump[1];
						Nj[1] = dstress[3] * n_jump[0] + dstress[1] * n_jump[1];

						for (i = 0; i < 4; i++)
							Aac[i] = 0.0;
						for (i = 0; i < 2; i++)
						{
							Aac[i * 2 + i] = smat->G;
							for (j = 0; j < 2; j++)
								Aac[i * 2 + j]
								    += (smat->K + smat->G / 3.0) * n_jump[i] * n_jump[j]
								       - (4.5 * smat->Al * smat->Xi * smat->K * smat->K * n_jump[i] * n_jump[j]
								          + 3.0 * smat->G * smat->K
								                * (smat->Al * Nj[i] * n_jump[j] + smat->Xi * Nj[j] * n_jump[i])
								          + 2.0 * smat->G * smat->G * Nj[i] * Nj[j])
								             / (4.5 * smat->Al * smat->Xi * smat->K + smat->G);
						}
						detA = Aac[0] * Aac[3] - Aac[1] * Aac[2];
						if (detA <= 0.0)
						{
							LOCed = true;
							break;
						}
					}
					for (i = 0; i < ns; i++)
						dstress[i] *= normXi;
				}
				//
				if (fabs(h_loc) < h_tol)
					LOCed = true;

				if (LOCed)
				{
					eleV_DM->orientation = new double[ele_dim];
					for (i = 0; i < 2; i++)
						eleV_DM->orientation[i] = OriJ[i];
					eleV_DM->Localized = true;

					// Choose one orientation. Empirical formular. 2D
					if (!Localizing)
					{
						for (i = 0; i < 3; i++)
							dstress[i] += p / 3.0;
						// WW n1 = (dstress[0]*cos(0.5*pai+OriJ[0])+dstress[1]*sin(0.5*pai+OriJ[0]))/normXi;
						n2 = (dstress[0] * cos(0.5 * pai + OriJ[1]) + dstress[1] * sin(0.5 * pai + OriJ[1])) / normXi;
						if (n2 > 0.0)
						{
							// Always use orientation[0]
							eleV_DM->orientation[0] = OriJ[1];
							eleV_DM->orientation[1] = OriJ[0];
						}
					}
				}
			}
		}
		return LOCed;
	}

	/**************************************************************************
	   GeoSys - Funktion:
	    void CFiniteElementVec::IntersectionPoint(const int O_edge, const double k,
	                                       const double *NodeA, double nodeB  )
	    2D only
	   Aufgabe:
	        Determine the second intersection point of a line and the element.
	        2D only (Geometry)
	   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
	   E :
	          const int O_edge       :   Edge will the departure point on
	   const double *NodeA    :   original node
	   double nodeB           :   Intersection point

	   Programmaenderungen:
	   06/2004   WW  Erste Version
	**************************************************************************/
	int CFiniteElementVec::IntersectionPoint(const int O_edge, const double* NodeA, double* NodeB)
	{
		int i, j, k, numf; //, nfnode;

		static double k1, k2, n1, n2, xA[3], xB[3];
		static int Face_node[8]; // Only 2D

		double Tol = 1.0e-12;
		double area0, area1;
		area0 = MeshElement->GetVolume();

		k = -1;

		eleV_DM = ele_value_dm[Index];

		numf = MeshElement->GetFacesNumber();
		for (i = 0; i < numf; i++)
		{
			k = -1;
			if (i != O_edge)
			{
				// WW nfnode =
				MeshElement->GetElementFaceNodes(i, Face_node);

				xA[0] = X[Face_node[0]];
				xA[1] = Y[Face_node[0]];
				xA[2] = Z[Face_node[0]];
				xB[0] = X[Face_node[1]];
				xB[1] = Y[Face_node[1]];
				xB[2] = Z[Face_node[1]];

				n1 = cos(eleV_DM->orientation[0]);
				n2 = sin(eleV_DM->orientation[0]);

				// parallel
				if (fabs((xB[0] - xA[0]) * n1 + (xB[1] - xA[1]) * n2) < Tol)
					continue;

				if (fabs(n2) < Tol)
				{
					NodeB[0] = NodeA[0];
					NodeB[1] = (NodeA[0] - xA[0]) * (xB[1] - xA[1]) / (xB[0] - xA[0]) + xA[1];
				}
				else if (fabs(xB[0] - xA[0]) < Tol)
				{
					NodeB[0] = xA[0];
					NodeB[1] = -n1 * (xA[0] - NodeA[0]) / n2 + NodeA[1];
				}
				else
				{
					k1 = (xB[1] - xA[1]) / (xB[0] - xA[0]);
					k2 = -n1 / n2;
					NodeB[0] = (NodeA[1] - xA[1] + k1 * xA[0] - k2 * NodeA[0]) / (k1 - k2);
					NodeB[1] = k1 * (NodeB[0] - xA[0]) + xA[1];
				}
				k = i;
			}

			// Check if this point is on an edge of this element.
			if (k >= 0) // Has intersection
			{
				area1 = 0.0;
				for (j = 0; j < numf; j++)
				{
					if (j == k)
						continue;
					// WW nfnode =
					MeshElement->GetElementFaceNodes(j, Face_node);
					xA[0] = X[Face_node[0]];
					xA[1] = Y[Face_node[0]];
					xA[2] = Z[Face_node[0]];
					xB[0] = X[Face_node[1]];
					xB[1] = Y[Face_node[1]];
					xB[2] = Z[Face_node[1]];
					area1 += ComputeDetTri(xA, xB, NodeB);
				}
				if (fabs(area0 - area1) < Tol)
					break;
			}
		}
		return k; // Local index of other intersection point
	}

	/**************************************************************************
	   GeoSys - Funktion: void CFiniteElementVec::LocalAssembly_EnhancedStrain

	   Aufgabe:
	     Local assembly within the strong discontinuity assumption
	    (Drucker-Prager model and 2D only)

	   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
	   E :
	    const int update: 1 update Gauss values
	   0 do not update

	   Ergebnis:
	   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

	   Programmaenderungen:
	   06/2004   WW  Erste Version
	**************************************************************************/
	void CFiniteElementVec::LocalAssembly_EnhancedStrain(const int update)
	{
		int gp_r(0), gp_s(0), gp_t(0);
		double fkt = 0.0, area, Jac_e = 0.0, f_j;

		// For enhanced strain element
		double f_tol = 1.0e-4; // Triangle 1.0e-8;
		double loc_dilatancy = 0.0, zeta_t0, zeta_t1;
		double sign, sj0, sj = 0.0;
		static double tt0[2], tt[2], zeta[2], D_prj[2];

		bool isLoop = true; // Used only to avoid warnings with .net

		double ThermalExpansion = 0.0, Tem = 0.0;
		gp = 0;
		BDG->LimitSize(2, 2 * nnodesHQ);
		PDB->LimitSize(2 * nnodesHQ, 2);

		if (T_Flag)
		{
			// Thermal effect
			if (smat->Thermal_Expansion() > 0.0)
				ThermalExpansion = smat->Thermal_Expansion();
			for (int i = 0; i < nnodes; i++)
				Temp[i] = t_pcs->GetNodeValue(nodes[i], idx_T1) - t_pcs->GetNodeValue(nodes[i], idx_T0);
		}

		// Elastic modulus
		smat->Calculate_Lame_Constant();

		smat->ElasticConsitutive(ele_dim, De);

		// Plasticity
		smat->CalulateCoefficent_DP();

		loc_dilatancy = smat->Al * sqrt(2.0 / (1.0 - 3.0 * smat->Al * smat->Al));
		double Hd = 0.5 * smat->Hard_Loc / (1.0 - 3.0 * smat->Al * smat->Al);

		/*
		   // Reference
		   loc_dilatancy = 1.5*sqrt(2.0)*smat->Al*sqrt(1.0-6.0*smat->Al*smat->Al);
		   double Hd = smat->Hard_Loc/(3.0-18.0*smat->Al*smat->Al);
		 */

		zeta_t0 = eleV_DM->disp_j;
		zeta_t1 = zeta_t0;

		area = MeshElement->volume;

		n_jump[0] = cos(eleV_DM->orientation[0]);
		n_jump[1] = sin(eleV_DM->orientation[0]);

		// Compute traction on the jump plane
		// Compute Pe, singular part of enhanced strain-jump matrix
		ComputeSESM();
		// PeDe=P^t *D_e
		(*PeDe) = 0.0;
		Pe->multi(*De, *PeDe);
		//

		// If this is the beginning of localization
		if (fabs(eleV_DM->tract_j) < MKleinsteZahl)
		{
			// average of stresses within an element
			for (int j = 0; j < ns; j++)
			{
				dstress[j] = 0.0;
				for (int i = 0; i < nGaussPoints; i++)
					dstress[j] += (*eleV_DM->Stress)(j, i);
			}
			for (int i = 0; i < ns; i++)
				dstress[i] /= (double)nGaussPoints;

			for (size_t i = 0; i < ele_dim; i++)
			{
				tt0[i] = 0.0;
				for (int j = 0; j < ns; j++)
					tt0[i] += (*Pe)(i, j) * dstress[j];
			}
			eleV_DM->tract_j = loc_dilatancy * tt0[0] + fabs(tt0[1]);
			/*
			   //
			   for(gp=0; gp<nGaussPoints; gp++)
			   {
			    //--------------------------------------------------------------
			    //-----------  Integrate of traction on the jump plane ---------
			    //--------------------------------------------------------------
			    for(int i=0; i<ns; i++) dstress[i] = (*eleV_DM->Stress)(i,gp);
			    fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
			    for(size_t i=0; i<ele_dim; i++)
			    {
			   tt0[i] = 0.0;
			   for(int j=0; j<ns; j++)
			   tt0[i] += (*Pe)(i,j)*dstress[j];
			   }
			   eleV_DM->tract_j = fkt*(loc_dilatancy*tt0[0]+fabs(tt0[1]));
			   }
			   eleV_DM->tract_j /= area;
			 */
		}
		//
		sj0 = eleV_DM->tract_j;
		//
		CheckNodesInJumpedDomain();
		// On discontinuity by enhanced strain
		// AuxMatrix temporarily used to store PDG
		(*AuxMatrix) = 0.0;
		// Integration of P^t*Stress^{elastic try}
		for (size_t i = 0; i < ele_dim; i++)
			tt0[i] = 0.0;
		// TEST
		for (int i = 0; i < ns; i++)
			dstress[i] = 0.0; // Test average appoach
		for (gp = 0; gp < nGaussPoints; gp++)
		{
			//--------------------------------------------------------------
			//-----------  Integrate of traction on the jump plane ---------
			//--------------------------------------------------------------
			fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
			//---------------------------------------------------------
			// Compute geometry
			//---------------------------------------------------------
			getGradShapefunctValues(gp, 2);
			ComputeStrain(gp);
			// Compute Ge, regular part of enhanced strain-jump matrix
			ComputeRESM();
			if (T_Flag) // Contribution by thermal expansion
			{
				getShapefunctValues(gp, 1); // Linear interpolation function
				Tem = 0.0;
				for (int i = 0; i < nnodes; i++)
					Tem += shapefct[i] * Temp[i];
				for (size_t i = 0; i < 3; i++)
					dstrain[i] -= ThermalExpansion * Tem;
			}
			/*
			   for(int i=0; i<ns; i++) dstress[i] = (*eleV_DM->Stress)(i,gp);
			   De->multi(dstrain, dstress);
			   // Try stress
			   Pe->multi(dstress, tt0, fkt);
			 */
			// Pe*De*Ge
			PeDe->multi(*Ge, *AuxMatrix, fkt);

			// TEST -----------Average approach -------------------------
			for (int i = 0; i < ns; i++)
				dstress[i] += (*eleV_DM->Stress)(i, gp);
			De->multi(dstrain, dstress);
			//-----------------------------------------------------------
		}

		// TEST average approach
		for (int i = 0; i < ns; i++)
			dstress[i] /= (double)nGaussPoints;
		Pe->multi(dstress, tt0, 1.0);
		for (size_t i = 0; i < ele_dim; i++)
			tt0[i] *= area;
		//-------------------------------------------------------------

		//  Local Newton iteration for discontinuity number
		while (isLoop)
		{
			zeta[1] = zeta_t1 - zeta_t0;
			zeta[0] = loc_dilatancy * fabs(zeta[1]);

			// Sign(zeta_t)
			if (fabs(zeta[1]) < MKleinsteZahl)
				sign = 1.0;
			else
				sign = zeta[1] / fabs(zeta[1]);

			//
			D_prj[0] = sign * loc_dilatancy;
			D_prj[1] = 1.0;

			sj = sj0 + Hd * zeta[1];

			//--------------------------------------------------------------
			//----------------   Local Jacobian   --------------------------
			//--------------------------------------------------------------
			for (size_t i = 0; i < ele_dim; i++)
				tt[i] = tt0[i];
			AuxMatrix->multi(zeta, tt, -1.0);

			Jac_e = 0.0;
			for (size_t i = 0; i < ele_dim; i++)
				for (size_t j = 0; j < ele_dim; j++)
					Jac_e += D_prj[i] * (*AuxMatrix)(i, j) * D_prj[j];

			Jac_e /= area;
			for (size_t i = 0; i < ele_dim; i++)
				tt[i] /= area;

			f_j = D_prj[0] * tt[0] + tt[1] - sj;

			Jac_e += Hd;
			if (fabs(f_j) < f_tol)
				break;

			zeta_t1 += f_j / Jac_e;
		} // Loop of the local Newton for enhanced parameter

		// Compute local RHS
		(*BDG) = 0.0;
		(*PDB) = 0.0;
		for (gp = 0; gp < nGaussPoints; gp++)
		{
			//--------------------------------------------------------------
			//-----------  Integrate of traction on the jump plane ---------
			//--------------------------------------------------------------
			for (int i = 0; i < ns; i++)
				dstress[i] = (*eleV_DM->Stress)(i, gp);
			fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
			//---------------------------------------------------------
			// Compute geometry
			//---------------------------------------------------------
			getGradShapefunctValues(gp, 2);
			ComputeStrain(gp);
			// Compute Ge
			ComputeRESM();

			if (T_Flag) // Contribution by thermal expansion
			{
				getShapefunctValues(gp, 1); // Linear interpolation function
				Tem = 0.0;
				for (int i = 0; i < nnodes; i++)
					Tem += shapefct[i] * Temp[i];
				for (size_t i = 0; i < 3; i++)
					dstrain[i] -= ThermalExpansion * Tem;
			}

			// Ehhanced strain:
			Ge->multi(zeta, dstrain, -1.0);

			// Compute try stress. 1. stress incremental:
			De->multi(dstrain, dstress);

			// Update stresses if convergence is reached
			if (update)
			{
				// Two Dimensional
				RecordGuassStrain(gp, gp_r, gp_s, gp_t);
				for (int i = 0; i < ns; i++)
				{
					(*eleV_DM->Stress)(i, gp) = dstress[i];
					dstrain[i] = 0.0;
				}
				Ge->multi(zeta, dstrain, 1.0);
				DeviatoricStress(dstrain);
				(*eleV_DM->pStrain)(gp) += sqrt(2.0 * TensorMutiplication2(dstrain, dstrain, 2) / 3.0);
			}
			else
			{
				// Compute stiffness matrix
				ComputeMatrix_RHS(fkt, De);
				// Stiffness contributed by enhanced strain to the stiffness matrix
				for (int i = 0; i < nnodesHQ; i++)
				{
					setTransB_Matrix(i);
					// B^T*D*G
					(*AuxMatrix) = 0.0;
					B_matrix_T->multi(*De, *Ge, *AuxMatrix);

					for (size_t k = 0; k < ele_dim; k++)
						for (size_t l = 0; l < ele_dim; l++)
							(*BDG)(k, ele_dim* i + l) += fkt * (*AuxMatrix)(k, l);
					//
					// P*D*B
					setB_Matrix(i);
					(*AuxMatrix) = 0.0;
					PeDe->multi(*B_matrix, *AuxMatrix);
					for (size_t k = 0; k < ele_dim; k++)
						for (size_t l = 0; l < ele_dim; l++)
							(*PDB)(ele_dim* i + k, l) += fkt * (*AuxMatrix)(k, l) / area;
				}
			}
		} // End of RHS assembly

		// Those contributed by enhanced strain to the stiffness matrix
		// D*D^T
		for (size_t i = 0; i < ele_dim; i++)
			for (size_t j = 0; j < ele_dim; j++)
				(*DtD)(i, j) = D_prj[i] * D_prj[j];

		//
		for (int i = 0; i < nnodesHQ; i++)
			for (int j = 0; j < nnodesHQ; j++)
			{
				// Local assembly of stiffness matrix
				for (size_t k = 0; k < ele_dim; k++)
					for (size_t l = 0; l < ele_dim; l++)
					{
						f_j = 0.0;
						for (size_t ii = 0; ii < ele_dim; ii++)
							for (size_t jj = 0; jj < ele_dim; jj++)
								f_j += (*BDG)(k, ele_dim* i + ii) * (*DtD)(ii, jj) * (*PDB)(ele_dim* j + jj, l);
						(*Stiffness)(i + k * nnodesHQ, l * nnodesHQ + j) -= f_j / Jac_e;
					}
			}

		if (update)
		{
			// Update strains.
			// The mapping of Gauss point strain to element nodes
			ExtropolateGuassStrain();
			// Update enhanced parameter
			eleV_DM->disp_j = zeta_t1;
			eleV_DM->tract_j = sj;
		}
	}


	/***************************************************************************
	   GeoSys - Funktion:
	            CFiniteElementVec:: CalcStrain_v()
	   Aufgabe:
	           Calculate effictive strain at Gauss points
	   Formalparameter:
	           E:

	   Programming:
	   01/2009   WW/UWG
	 **************************************************************************/
	double CFiniteElementVec::CalcStrain_v()
	{
		for (int j(0); j < ns; j++)
		{
			dstrain[j] = 0.0;
			for (int i = 0; i < nnodesHQ; i++)
				dstrain[j] += pcs->GetNodeValue(nodes[i], _nodal_strain_indices[j]) * shapefctHQ[i];
		}
		double val = 0;
		for (int i = 0; i < 3; i++)
			val += dstrain[i] * dstrain[i];
		for (int i = 3; i < ns; i++)
			val += 0.5 * dstrain[i] * dstrain[i];

		return sqrt(2.0 * val / 3.);
	}

	/*----------------------------------------------------------------
	   Class ElementValue_DM

	   Allocate memory for element value
	     Matrix *Mat:
	 | Index  |  Paramete |
	                          ----------------------
	 |    0   |  alpha    |
	 |    1   |  beta     |
	 |    2   |  delta    |
	 |    3   |  epsilon  |
	 |    4   |  kappa    |
	 |    5   |  gamma    |
	 |    6   |  m        |
	   ----------------------
	   -----------------------------------------------------------------*/
	ElementValue_DM::ElementValue_DM(CElem * ele, const int NGP, bool HM_Staggered)
	    : Stress(NULL), Stress_i(NULL), Stress_j(NULL), pStrain(NULL),
	      y_surface(NULL), prep0(NULL), e_i(NULL), xi(NULL), MatP(NULL),
	      Strain_Kel(NULL),	Strain_Max(NULL), Strain_pl(NULL),
	      Strain_t_ip(NULL), e_pl(NULL), ev_loc_nr_res(NULL),
	      lambda_pl(NULL), Strain(NULL), NodesOnPath(NULL),
	      orientation(NULL), scalar_aniso_comp(NULL),
	     scalar_aniso_tens(NULL)
	{
		int Plastic = 1;
		const int LengthMat = 7; // Number of material parameter of SYS model.
		int LengthBS = 4; // Number of stress/strain components
		int NGPoints = 0;
		CSolidProperties* sdp = NULL;
		int ele_dim;
		//
		MshElemType::type ele_type = ele->GetElementType();
		ele_dim = ele->GetDimension();
		sdp = msp_vector[ele->GetPatchIndex()];
		Plastic = sdp->Plastictity();

		if (ele_dim == 2)
			LengthBS = 4;
		else if (ele_dim == 3)
			LengthBS = 6;

		if (ele_type == MshElemType::TRIANGLE)
			NGPoints = 3;
		else if (ele_type == MshElemType::TETRAHEDRON)
			NGPoints = 5; // 15
		else if (ele_type == MshElemType::PRISM)
			NGPoints = 6; // 9
		else
			NGPoints = MathLib::fastpow(NGP, ele_dim);

		Stress0 = new Matrix(LengthBS, NGPoints);
		Stress_i = new Matrix(LengthBS, NGPoints);
		Stress = Stress_i;
		if (HM_Staggered)
			Stress_j = new Matrix(LengthBS, NGPoints);
		//
		if (Plastic > 0)
		{
			pStrain = new Matrix(NGPoints);
			y_surface = new Matrix(NGPoints);
			*y_surface = 0.0;
			*pStrain = 0.0;
		}
		*Stress = 0.0;

		if (Plastic == 2) // Rotational hardening model
		{
			xi = new Matrix(LengthBS - 1, NGPoints);
			MatP = new Matrix(LengthMat, NGPoints);
			*xi = 0.0;
			*MatP = 0.0;
		}
		if (Plastic == 3) // Cam-Clay
		{
			prep0 = new Matrix(NGPoints);
			e_i = new Matrix(NGPoints);
			*prep0 = 0.0;
			*e_i = 0.0;
		}
		if (sdp->CreepModel() == 1000)
		{
			xi = new Matrix(LengthBS);
			*xi = 0.0;
		}
		if (sdp->CreepModel() == 1001) // Burgers
		{
			Strain_Kel = new Matrix(6, NGPoints); // Only 3D size for now. Will work with all. Separation into special
												  // cases (LengthBS) may follow later
			*Strain_Kel = 0.0;
			Strain_Max = new Matrix(6, NGPoints); // Only 3D size for now. Will work with all. Separation into special
												  // cases (LengthBS) may follow later
			*Strain_Max = 0.0;
			Strain_t_ip = new Matrix(6, NGPoints);
			*Strain_t_ip = 0.0;
			ev_loc_nr_res = new Matrix(NGPoints);
			*ev_loc_nr_res = 0.;
			Strain = new Matrix(6, NGPoints);
			*Strain = 0.0;
		}

		if (sdp->CreepModel() == 1002) // Minkley
		{
			Strain_Kel = new Matrix(6, NGPoints); // Only 3D size for now. Will work with all. Separation into special
												  // cases (LengthBS) may follow later
			*Strain_Kel = 0.0;
			Strain_Max = new Matrix(6, NGPoints); // Only 3D size for now. Will work with all. Separation into special
												  // cases (LengthBS) may follow later
			*Strain_Max = 0.0;
			Strain_pl = new Matrix(6, NGPoints); // Only 3D size for now. Will work with all. Separation into special
												 // cases (LengthBS) may follow later
			*Strain_pl = 0.0;
			Strain_t_ip = new Matrix(6, NGPoints);
			*Strain_t_ip = 0.0;
			e_pl = new Matrix(NGPoints);
			*e_pl = 0.;
			pStrain = new Matrix(NGPoints);
			*pStrain = 0.;
			lambda_pl = new Matrix(NGPoints);
			*lambda_pl = 0.;
			ev_loc_nr_res = new Matrix(NGPoints);
			*ev_loc_nr_res = 0.;
			Strain = new Matrix(6, NGPoints);
			*Strain = 0.0;
		}
		disp_j = 0.0;
		tract_j = 0.0;
		Localized = false;
		if (sdp->Plasticity_Bedding)
		{
			scalar_aniso_comp = new Matrix(NGPoints);
			scalar_aniso_tens = new Matrix(NGPoints);
			*scalar_aniso_comp = 0.;
			*scalar_aniso_tens = 0.;
		}
	}
	// 01/2006 WW
	void ElementValue_DM::Write_BIN(std::fstream & os, const bool last_step)
	{
		if (last_step)
			Stress_i->Write_BIN(os);
		else
			Stress0->Write_BIN(os);
		Stress_i->Write_BIN(os);
		if (pStrain)
			pStrain->Write_BIN(os);
		if (y_surface)
			y_surface->Write_BIN(os);
		if (xi)
			xi->Write_BIN(os);
		if (Strain) // NB
			Strain->Write_BIN(os);
		if (Strain_Kel) // TN Burgers
			Strain_Kel->Write_BIN(os);
		if (Strain_Max) // TN Burgers
			Strain_Max->Write_BIN(os);
		if (Strain_pl) // TN
			Strain_pl->Write_BIN(os);
		if (Strain_t_ip)
			Strain_t_ip->Write_BIN(os);
		if (e_pl)
			e_pl->Write_BIN(os);
		if (lambda_pl)
			lambda_pl->Write_BIN(os);
		if (ev_loc_nr_res)
			ev_loc_nr_res->Write_BIN(os);
		if (MatP)
			MatP->Write_BIN(os);
		if (prep0)
			prep0->Write_BIN(os);
		if (e_i)
			e_i->Write_BIN(os);
		if (NodesOnPath)
			NodesOnPath->Write_BIN(os);
		if (orientation)
			os.write((char*)(orientation), sizeof(*orientation));
		os.write((char*)(&disp_j), sizeof(disp_j));
		os.write((char*)(&tract_j), sizeof(tract_j));
		os.write((char*)(&Localized), sizeof(Localized));
	}
	// 01/2006 WW
	void ElementValue_DM::Read_BIN(std::fstream & is)
	{
		Stress0->Read_BIN(is);
		Stress_i->Read_BIN(is);
		if (pStrain)
			pStrain->Read_BIN(is);
		if (y_surface)
			y_surface->Read_BIN(is);
		if (xi)
			xi->Read_BIN(is);
		if (Strain) // NB
			Strain->Read_BIN(is);
		if (Strain_Kel) // TN Burgers
			Strain_Kel->Read_BIN(is);
		if (Strain_Max) // TN Burgers
			Strain_Max->Read_BIN(is);
		if (Strain_pl) // TN
			Strain_pl->Read_BIN(is);
		if (Strain_t_ip)
			Strain_t_ip->Read_BIN(is);
		if (e_pl)
			e_pl->Read_BIN(is);
		if (lambda_pl)
			lambda_pl->Read_BIN(is);
		if (ev_loc_nr_res)
			ev_loc_nr_res->Read_BIN(is);
		if (MatP)
			MatP->Read_BIN(is);
		if (prep0)
			prep0->Read_BIN(is);
		if (e_i)
			e_i->Read_BIN(is);
		if (NodesOnPath)
			NodesOnPath->Read_BIN(is);
		if (orientation)
			is.read((char*)(orientation), sizeof(*orientation));
		is.read((char*)(&disp_j), sizeof(disp_j));
		is.read((char*)(&tract_j), sizeof(tract_j));
		is.read((char*)(&Localized), sizeof(Localized));
	}

	// 10/2011 WW
	void ElementValue_DM::ReadElementStressASCI(std::fstream & is)
	{
		size_t i, j;
		size_t ns = Stress0->Rows();
		size_t nGS = Stress0->Cols();

		for (i = 0; i < ns; i++)
		{
			is >> (*Stress0)(i, 0);
			for (j = 1; j < nGS; j++)
				(*Stress0)(i, j) = (*Stress0)(i, 0);
		}

		*Stress_i = *Stress0;
	}

	void ElementValue_DM::ResetStress(bool cpl_loop)
	{
		if (cpl_loop) // For coupling loop
		{
			(*Stress_j) = (*Stress_i);
			Stress = Stress_j;
		}
		else // Time loop
		{
			(*Stress_i) = (*Stress_j);
			Stress = Stress_i;
		}
	}

	ElementValue_DM::~ElementValue_DM()
	{
		delete Stress0;
		if (Stress_i)
			delete Stress_i;
		if (Stress_j)
			delete Stress_j;
		if (pStrain)
			delete pStrain;
		if (y_surface)
			delete y_surface;

		// Preconsolidation pressure
		if (prep0)
			delete prep0;
		if (e_i)
			delete e_i; // Void ratio
		// Variables of single yield surface model
		if (xi)
			delete xi; // Rotational hardening variables
		if (MatP)
			delete MatP; // Material parameters

		if (NodesOnPath)
			delete NodesOnPath;
		if (orientation)
			delete orientation;

		if (scalar_aniso_comp)
			delete scalar_aniso_comp; // WX:09.2011
		if (scalar_aniso_tens)
			delete scalar_aniso_tens;
		if (Strain)
			delete Strain;
		if (Strain_Kel) // TN Strain in Kelvin element
			delete Strain_Kel;
		if (Strain_Max) // TN Strain in Maxwell element
			delete Strain_Max;
		if (Strain_pl) // TN - Minkley plastic deviatoric strain
			delete Strain_pl;
		if (Strain_t_ip)
			delete Strain_t_ip;
		if (e_pl)
			delete e_pl;
		if (lambda_pl)
			delete lambda_pl;
		if (ev_loc_nr_res)
			delete ev_loc_nr_res;

		NodesOnPath = NULL;
		orientation = NULL;
		y_surface = NULL;
		Stress0 = NULL;
		Stress = NULL;
		Stress_i = NULL; // for HM coupling iteration
		Stress_j = NULL; // for HM coupling iteration
		pStrain = NULL;
		prep0 = NULL;
		e_i = NULL;
		xi = NULL;
		MatP = NULL;
		scalar_aniso_comp = NULL;
		scalar_aniso_tens = NULL;
		Strain = NULL;
		Strain_Kel = NULL;
		Strain_Max = NULL;
		Strain_pl = NULL;
		Strain_t_ip = NULL;
		e_pl = NULL;
		lambda_pl = NULL;
		ev_loc_nr_res = NULL;
	}

	double ElementValue_DM::MeanStress(const int gp)
	{
		return (*Stress)(0, gp) + (*Stress)(1, gp) + (*Stress)(2, gp);
	}

} // end namespace FiniteElement
