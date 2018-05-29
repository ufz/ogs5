/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Object: Process PCS
   Programing:
   02/2003 OK Implementation
   /2003 WW CRFProcessDeformation
   11/2003 OK re-organized
   07/2004 OK PCS2
   02/2005 WW/OK Element Assemblier and output
   12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear solver
           and parellelisation of them
   02/2008 PCH OpenMP parallelization for Lis matrix solver
**************************************************************************/
#include "FEMEnums.h"
#include "Output.h"
#include "MathTools.h"

/*--------------------- MPI Parallel  -------------------*/
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include <mpi.h>
#endif
/*--------------------- MPI Parallel  -------------------*/

/*--------------------- OpenMP Parallel ------------------*/
#if defined(LIS)
#include "lis.h"
#include <omp.h>
#endif
/*--------------------- OpenMP Parallel ------------------*/

#include "makros.h"
// C
#ifndef __APPLE__
#include <malloc.h>
#endif

// C++
#include <cfloat>
#include <iomanip> //WW
#include <iostream>
#include <algorithm>
#include <set>

#include "isnan.h"
#include "display.h"
#include "memory.h"
// GEOLib
#include "PointWithID.h"

/*------------------------------------------------------------------------*/
/* MshLib */
//#include "msh_elem.h"
//#include "msh_lib.h"
/*-----------------------------------------------------------------------*/
/* Objects */
#include "pcs_dm.h"
#include "rf_pcs.h"
#include "rf_st_new.h" // ST
//#include "rf_bc_new.h" // ST
//#include "rf_mmp_new.h" // MAT
#include "fem_ele_std.h" // ELE
#include "rf_ic_new.h" // IC
//#include "msh_lib.h" // ELE
//#include "rf_tim_new.h"
//#include "rf_out_new.h"
#include "rfmat_cp.h"
//#include "rf_mfp_new.h" // MFP
//#include "rf_num_new.h"
//#include "gs_project.h"
#include "rf_fct.h"
//#include "femlib.h"
#include "eos.h"
#include "rf_msp_new.h"
#include "rf_node.h"
#include "rf_kinreact.h"
#include "fem_ele_vec.h" //WX:08.2011

#ifdef MFC // WW
#include "rf_fluid_momentum.h"
#endif
/* Tools */
#include "mathlib.h"
//#include "files0.h"
//#include "par_ddc.h"
#include "tools.h"
//#include "rf_pcs.h"
#include "files0.h"
#ifdef GEM_REACT
// GEMS chemical solver
#include "rf_REACT_GEM.h"
REACT_GEM* m_vec_GEM;
#endif
#ifdef BRNS
// BRNS chemical solver
#include "rf_REACT_BRNS.h"
REACT_BRNS* m_vec_BRNS;
#endif

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
#include "PETSC/PETScLinearSolver.h"
// New EQS
#elif defined(NEW_EQS)
#include "equation_class.h"
#else
#include "solver.h" // ConfigRenumberProperties
#include "matrix_routines.h"
#endif
#include "problem.h"
#include "msh_faces.h"
#include "rfmat_cp.h"

#include "rf_react_int.h"
#include "VLE.h"
#include "Density.h"

// MathLib
#include "InterpolationAlgorithms/InverseDistanceInterpolation.h"
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "StringTools.h"

#include "fct_mpi.h"

using namespace std;
using namespace MeshLib;
using namespace Math_Group;

/*-------------------- ITPACKV    ---------------------------*/
extern void transM2toM6(void);
/*-------------------- ITPACKV    ---------------------------*/
/*-------------------- JAD    ---------------------------*/
extern void transM2toM5(void);
/*-------------------- JAD    ---------------------------*/
/*-----------------------------------------------------------------------*/
/* LOP */
// 16.12.2008. WW #include "rf_apl.h" // Loop...
// 16.12.2008. WW #include "loop_pcs.h"
extern VoidFuncVoid LOPCalcSecondaryVariables_USER;
//------------------------------------------------------------------------
// PCS
VoidXFuncVoidX PCSDestroyELEMatrices[PCS_NUMBER_MAX];
//------------------------------------------------------------------------
// Globals, to be checked
int pcs_no_components = 0;
bool pcs_monolithic_flow = false;
int pcs_deformation = 0;
int dm_number_of_primary_nvals = 2;
bool show_onces_adp = true;
bool show_onces_mod = true;
bool show_onces_mod_flow = true;
bool show_onces_density = true;
int memory_opt = 0;
int problem_2d_plane_dm;
int anz_nval = 0;
int anz_nval0 = 0; // WW
//
int size_eval = 0; // WW

NvalInfo* nval_data = NULL;
int anz_eval = 0;
EvalInfo* eval_data = NULL;
string project_title("New project"); // OK41

bool hasAnyProcessDeactivatedSubdomains = false; // NW
extern double gravity_constant;
//--------------------------------------------------------
// Coupling Flag. WW
bool T_Process = false; // Heat
bool H_Process = false; // Fluid
bool H2_Process = false; // Multi-phase
bool H3_Process = false; // 3-phase
bool M_Process = false; // Mechanical
bool RD_Process = false; // Richards
bool MH_Process = false; // MH monolithic scheme
bool MASS_TRANSPORT_Process = false; // Mass transport
bool FLUID_MOMENTUM_Process = false; // Momentum
bool RANDOM_WALK_Process = false; // RWPT
bool MULTI_COMPONENTIAL_FLOW_Process = false; // MCF
bool pcs_created = false;
//
int pcs_number_deformation = -1; // JT2012
int pcs_number_flow = -1; // JT2012
int pcs_number_heat = -1; // JT2012
vector<int> pcs_number_mass; // JT2012

namespace process
{
class CRFProcessDeformation;
}
using process::CRFProcessDeformation;
using MeshLib::CNode;
using MeshLib::CElem;
using FiniteElement::ElementValue;
using FiniteElement::ElementValue_DM; // WX:07.2011
using Math_Group::vec;

#define NOCHECK_EQS
#define noCHECK_ST_GROUP
#define noCHECK_BC_GROUP

extern size_t max_dim; // OK411 todo

//////////////////////////////////////////////////////////////////////////
// PCS vector
//////////////////////////////////////////////////////////////////////////
// It is better to have space between data type and data name. WW
vector<LINEAR_SOLVER*> PCS_Solver; // WW
vector<CRFProcess*> pcs_vector;
vector<double*> ele_val_vector; // PCH
// vector<string> ele_val_name_vector; // PCH
template <class T>
T* resize(T* array, size_t old_size, size_t new_size);
//////////////////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
   02/2005 WW Local elment assembly (all protected members)
   last modified:
**************************************************************************/
CRFProcess::CRFProcess(void)
    : _problem(NULL), p_var_index(NULL), num_nodes_p_var(NULL), fem(NULL), Memory_Type(0), Write_Matrix(false),
      matrix_file(NULL), WriteSourceNBC_RHS(0),
#ifdef JFNK_H2M
      JFNK_precond(false), norm_u_JFNK(NULL), array_u_JFNK(NULL), array_Fu_JFNK(NULL),
#endif
      ele_val_name_vector(std::vector<std::string>())
{
	iter_lin = 0;
	iter_lin_max = 0;
	iter_nlin = 0;
	iter_nlin_max = 0;
	iter_inner_cpl = 0;
	iter_outer_cpl = 0;
	TempArry = NULL;
	// SB:GS4  pcs_component_number=0; //SB: counter for transport components
	pcs_component_number = pcs_no_components - 1;
	//----------------------------------------------------------------------
	// NUM
	pcs_num_name[0] = NULL;
	pcs_num_name[1] = NULL;
	pcs_sol_name = NULL;
	m_num = NULL;
	cpl_type_name = "PARTITIONED"; // OK
	num_type_name = "FEM"; // OK
	rwpt_app = 0; // PCH Application types for RWPT such as Cell Dynamics, Crypto, etc.
	rwpt_count = 0; // YS 05.2013
	//
	for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
		pcs_number_mass.push_back(-1); // JT2012 (allow DOF_NUMBER_MAX potential components)
//
#if defined(USE_PETSC) // || defined(using other parallel scheme)//03.3012. WW
	eqs_new = NULL;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
	MPI_Comm_size(PETSC_COMM_WORLD, &mysize);

#elif defined(NEW_EQS) // WW 07.11.2008
	eqs_new = NULL;
	configured_in_nonlinearloop = false;
#else
	eqs = NULL; // WW
#endif
	dof = 1; // WW
	//
	// ITERATIONS AND COUPLING
	num_notsatisfied = 0;
	num_diverged = 0;
	cpl_max_relative_error = 0.0;
	cpl_max_relative_error_overall = 0.0;
	//
	// PCS designations
	isPCSDeformation = false;
	isPCSFlow = false;
	isPCSMultiFlow = false;
	isPCSHeat = false;
	isPCSMass = false;
	//
	//----------------------------------------------------------------------
	// ELE
	pcs_number_of_evals = 0;
	NumDeactivated_SubDomains = 0;
	Deactivated_SubDomain = NULL;
	//----------------------------------------------------------------------
	//
	mobile_nodes_flag = -1;
	//----------------------------------------------------------------------
	// USER
	PCSSetIC_USER = NULL;
	//----------------------------------------------------------------------
	// TIM
	tim_type = TimType::TRANSIENT;
	time_unit_factor = 1.0;
	timebuffer = 1.0e-5; // WW
	//_pcs_type_name.empty();
	adaption = false; // HS 03.2008
	this->EclipseData = NULL; // BG 09/2009, coupling to Eclipse
	this->DuMuxData = NULL; // SBG 09/2009, coupling to DuMux
	cpl_overlord = NULL;
	cpl_underling = NULL;
	pcs_is_cpl_overlord = false;
	pcs_is_cpl_underling = false;
	//----------------------------------------------------------------------
	// CPL
	for (int i = 0; i < 10; i++)
		Shift[i] = 0;
	selected = true; // OK
	// MSH OK
	m_msh = NULL;
	// Reload solutions
	reload = -1;
	nwrite_restart = 1; // kg44 write every timestep is default
	pcs_nval_data = NULL;
	pcs_eval_data = NULL;
	non_linear = false; // OK/CMCD
	cal_integration_point_value = false; // WW
	continuum = 0;
	// adaption = false; JOD removed
	compute_domain_face_normal = false; // WW
	use_velocities_for_transport = false;
	//
	additioanl2ndvar_print = -1; // WW
	flow_pcs_type = 0; // CB default: liquid flow, Sat = 1
	simulator = "GEOSYS"; // BG, 09/2009
	simulator_model_path = ""; // BG, 09/2009, folder with the Eclipse or DuMux files
	simulator_path = ""; // BG, 09/2009, Eclipse or Dumux
	simulator_well_path = ""; // KB 02/2011, Eclipse
	PrecalculatedFiles = false; // BG, 01/2011, flag for using already calculated files from Eclipse or DuMux
	SaveEclipseDataFiles = false; // Standard case: do not save files
	Phase_Transition_Model = 0; // BG, 11/2010, flag for using CO2 Phase transition (0...not used, 1...used)
	// SB redo wtp
	dissolved_co2_pcs_name
	    = ""; // SB, CB 10/2011 Name of MASS_TRANSPORT process used to store dissolved total CO2 in water
	dissolved_co2_ingas_pcs_name = "";

	//----------------------------------------------------------------------
	m_bCheck = false; // OK
	m_bCheckOBJ = false; // OK
	m_bCheckNOD = false; // OK
	m_bCheckELE = false; // OK
	m_bCheckEQS = false; // OK
	//
	write_boundary_condition = false; // 15.01.2008. WW
	OutputMassOfGasInModel = false; // 05/2012     BG
	WriteProcessed_BC = -1; // 26.08.2011. WW
	accepted = true; // 25.08.2008. WW
	accept_steps = 0; // 27.08.1008. WW
	reject_steps = 0; // 27.08.1008. WW
	ML_Cap = 0; // 23.01.2009 PCH
	PartialPS = 0; // 16.02 2009 PCH

#if defined(USE_MPI) || defined(USE_PETSC) // WW
	cpu_time_assembly = 0;
#endif
// New equation and solver WW
#ifdef NEW_EQS
	eqs_new = NULL;
	configured_in_nonlinearloop = false;
#endif
	flag_couple_GEMS = 0; // 11.2009 HS
	femFCTmode = false; // NW
	this->Gl_ML = NULL; // NW
	this->Gl_Vec = NULL; // NW
	this->Gl_Vec1 = NULL; // NW
	this->FCT_AFlux = NULL; // NW
#ifdef USE_PETSC
	this->FCT_K = NULL;
	this->FCT_d = NULL;
#endif
	ExcavMaterialGroup = -1; // 01.2010 WX
	PCS_ExcavState = -1; // WX
	Neglect_H_ini = -1; // WX
	m_conversion_rate = NULL; // WW
	isRSM = false; // WW
	eqs_x = NULL;
	_hasConstrainedBC = false;
	_hasConstrainedST = false;
	_pcs_constant_model = 0;
	_pcs_constant_value = .0;
	_pcs_constant_curve = 0;
}

void CRFProcess::setProblemObjectPointer(Problem* problem)
{
	_problem = problem;
}

Problem* CRFProcess::getProblemObjectPointer() const
{
	return _problem;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
   01/2005 WW Local assemblier as a member
   last modified:
**************************************************************************/
CRFProcess::~CRFProcess(void)
{
#ifdef USE_PETSC
	PetscPrintf(PETSC_COMM_WORLD, "\t\n>>Total Wall clock time in the assembly for %s (with PETSC):%f s\n",
	            FiniteElement::convertProcessTypeToString(this->getProcessType()).c_str(), cpu_time_assembly);

#endif
	long i;
	//----------------------------------------------------------------------
	// Finite element
	if (fem)
		delete fem; // WW
	fem = NULL;
	//----------------------------------------------------------------------
	// ELE: Element matrices
	ElementMatrix* eleMatrix = NULL;
	ElementValue* gp_ele = NULL;
	if (Ele_Matrices.size() > 0)
	{
		for (i = 0; i < (long)Ele_Matrices.size(); i++)
		{
			eleMatrix = Ele_Matrices[i];
			delete eleMatrix;
			eleMatrix = NULL;
		}
		Ele_Matrices.clear();
	}
	//----------------------------------------------------------------------
	// ELE: Element Gauss point values
	if (ele_gp_value.size() > 0)
	{
		for (i = 0; i < (long)ele_gp_value.size(); i++)
		{
			gp_ele = ele_gp_value[i];
			delete gp_ele;
			gp_ele = NULL;
		}
		ele_gp_value.clear();
	}
	//----------------------------------------------------------------------
	// OUT: Matrix output
	if (matrix_file)
	{
		matrix_file->close();
		delete matrix_file;
	}
	//----------------------------------------------------------------------
	// NOD: Release memory of node values
	for (i = 0; i < (int)nod_val_vector.size(); i++)
	{
		delete[] nod_val_vector[i]; // Add []. WW
		nod_val_vector[i] = NULL;
	}
	nod_val_vector.clear();
	//----------------------------------------------------------------------
	// ST:
	CNodeValue* m_nod_val = NULL;

	// Added &&m_nod_val for RSM model. 15.08.2011. WW
	if (!isRSM)
	{
		for (i = 0; i < (int)st_node_value.size(); i++)
		{
			m_nod_val = st_node_value[i];
			// OK delete st_node_value[i];
			// OK st_node_value[i] = NULL;
			if (m_nod_val->check_me) // OK
			{
				m_nod_val->check_me = false;
				delete m_nod_val;
				m_nod_val = NULL;
			}
		}
		st_node_value.clear();
	}
	//----------------------------------------------------------------------
	for (i = 0; i < (int)bc_node_value.size(); i++)
	{
		delete bc_node_value[i];
		bc_node_value[i] = NULL;
	}
	bc_node_value.clear();
	//----------------------------------------------------------------------
	//_pcs_type_name.clear();
	//----------------------------------------------------------------------
	// CON
	continuum_vector.clear();
	// Haibing 13112006------------------------------------------------------
	for (i = 0; i < (long)ele_val_vector.size(); i++)
		delete[] ele_val_vector[i];
	ele_val_vector.clear();
	//----------------------------------------------------------------------
	DeleteArray(Deactivated_SubDomain); // 05.09.2007 WW
	// 11.08.2010. WW
	DeleteArray(num_nodes_p_var);
	// 20.08.2010. WW
	DeleteArray(p_var_index);
#ifdef JFNK_H2M
	DeleteArray(array_u_JFNK); // 13.08.2010. WW
	DeleteArray(array_Fu_JFNK); // 31.08.2010. WW
	DeleteArray(norm_u_JFNK); // 24.11.2010. WW
#endif
	//----------------------------------------------------------------------
	if (this->m_num && this->m_num->fct_method > 0) // NW
	{
		delete this->Gl_ML;
		delete this->Gl_Vec;
		delete this->Gl_Vec1;
		delete this->FCT_AFlux;
		this->Gl_ML = NULL;
		this->Gl_Vec = NULL;
		this->Gl_Vec1 = NULL;
		this->FCT_AFlux = NULL;
#ifdef USE_PETSC
		delete this->FCT_K;
		delete this->FCT_d;
		this->FCT_K = NULL;
		this->FCT_d = NULL;
#endif
	}

	// HS, 11.2011
	if (m_conversion_rate)
		delete m_conversion_rate;

#if defined(USE_PETSC) // || defined(other parallel libs)//10.3012. WW
	PetscPrintf(PETSC_COMM_WORLD, "\n>>PETSc solver info for %s :\n",
	            FiniteElement::convertProcessTypeToString(this->getProcessType()).c_str());
	delete eqs_new;
	eqs_new = NULL;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:    Gauss point values for CFEMSH
   Programing:
   08/2005 WW Implementation
**************************************************************************/
void CRFProcess::AllocateMemGPoint()
{
	//	if (_pcs_type_name.find("FLOW") == 0)
	//		return;
	if (!ele_gp_value.empty())
		return;
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
		ele_gp_value.push_back(new ElementValue(this, m_msh->ele_vector[i]));
}

/**************************************************************************
   PCSLib-Method:
   10/2002 OK Implementation
   04/2004 WW Modification for 3D problems
   02/2005 WW New fem calculator
   04/2005 OK MSHCreateNOD2ELERelations
   07/2005 WW Geometry element objects
   02/2006 WW Removed memory leaking
   01/2006 YD MMP for each PCS
   04/2006 WW Unique linear solver for all processes if they share the same mesh
   06/2006 WW Rearrange incorporation of BC and ST. Set BC and ST for domain decomposition
**************************************************************************/
void CRFProcess::Create()
{
	if (hasAnyProcessDeactivatedSubdomains)
		CheckMarkedElement(); // 01.2014 WW

	// we need the string representation of process type at some points
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));

	// Element matrix output. WW
	if (Write_Matrix)
	{
		std::cout << "->Write Matrix" << '\n';
#if defined(USE_MPI) || defined(USE_PETSC)
		char stro[32];
		sprintf(stro, "%d", myrank);
		string m_file_name = FileName + "_" + pcs_type_name + (string)stro + "_element_matrix.txt";
#else
		std::string m_file_name = FileName + "_" + pcs_type_name + "_element_matrix.txt";
#endif
		matrix_file = new std::fstream(m_file_name.c_str(), ios::trunc | ios::out);
		if (!matrix_file->good())
			std::cout << "Warning in GlobalAssembly: Matrix files are not found"
			          << "\n";
	}
	//----------------------------------------------------------------------------
	if (m_msh) // OK->MB please shift to Config()

		//		if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0)
		if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			MSHDefineMobile(this);
	//----------------------------------------------------------------------------
	int DOF = GetPrimaryVNumber(); // OK should be PCS member variable
	//----------------------------------------------------------------------------
	// MMP - create mmp groups for each process //YD
	std::cout << "->Create MMP"
	          << "\n";
	CMediumPropertiesGroup* m_mmp_group = NULL;
	int continua = 1; // WW
	if (RD_Process)
		continua = 2;
	m_mmp_group = MMPGetGroup(pcs_type_name);
	for (int i = 0; i < continua; i++)
		if (!m_mmp_group)
		{
			m_mmp_group = new CMediumPropertiesGroup();
			m_mmp_group->pcs_type_name = pcs_type_name;
			m_mmp_group->Set(this);
			mmp_group_list.push_back(m_mmp_group);
		}
	m_mmp_group = NULL;
	//----------------------------------------------------------------------------
	// NUM_NEW
	std::cout << "->Create NUM"
	          << "\n";
	//	if (pcs_type_name.compare("RANDOM_WALK")) { // PCH RWPT does not need this.
	if (this->getProcessType() != FiniteElement::RANDOM_WALK) // PCH RWPT does not need this.
	{
		CNumerics* m_num_tmp = NULL;
		size_t no_numerics(num_vector.size());
		for (size_t i = 0; i < no_numerics; i++)
		{
			m_num_tmp = num_vector[i];

			if ((pcs_type_name.compare(m_num_tmp->pcs_type_name) == 0)
			    || (m_num_tmp->pcs_type_name.compare(pcs_primary_function_name[0]) == 0))
			{
				m_num = m_num_tmp;
				break;
			}
		}
	}
	if (!m_num)
	{
		std::cout << "Warning in CRFProcess::Create() - no NUM data"
		          << "\n";
		m_num = new CNumerics(pcs_type_name); // OK
		//		m_num = m_num_tmp;
	}
	else
	{
		if (m_num->nls_max_iterations > 1) // WW
			non_linear = true;
	}
	if (m_num->fct_method > 0) // NW
	{
// Memory_Type = 1;
#ifdef USE_PETSC
		long gl_size = m_msh->getNumNodesGlobal();
		this->FCT_K = new SparseMatrixDOK(gl_size, gl_size);
		this->FCT_d = new SparseMatrixDOK(gl_size, gl_size);
#else
		long gl_size = m_msh->GetNodesNumber(false);
#endif
		this->FCT_AFlux = new SparseMatrixDOK(gl_size, gl_size);
		this->Gl_ML = new Math_Group::Vec(gl_size);
		this->Gl_Vec = new Math_Group::Vec(gl_size);
		this->Gl_Vec1 = new Math_Group::Vec(gl_size);
	}
	//----------------------------------------------------------------------------
	// EQS - create equation system
	// WW CreateEQS();
	std::cout << "->Create EQS" << '\n';
#if !defined(USE_PETSC) // && !defined(other parallel solver lib). 04.2012 WW
#if defined(NEW_EQS)
	size_t k;
	for (k = 0; k < fem_msh_vector.size(); k++)
		if (m_msh == fem_msh_vector[k])
			break;
	// WW 02.2013. Pardiso
	int eqs_num = 3;
#ifdef USE_MPI
	eqs_num = 2;
#endif

	// if(type==4||type==41)
	//   eqs_new = EQS_Vector[2*k+1];
	if (type == 4 || (type / 10 == 4)) // 03.08.2010. WW
		eqs_new = EQS_Vector[eqs_num * k + 1];
	else
	{
// eqs_new = EQS_Vector[2*k];
#ifdef USE_MPI
		eqs_new = EQS_Vector[eqs_num * k];
#else
		if (getProcessType() == FiniteElement::MULTI_PHASE_FLOW || getProcessType() == FiniteElement::PS_GLOBAL
		    || getProcessType() == FiniteElement::TNEQ
		    || getProcessType() == FiniteElement::TES)
		{
			eqs_new = EQS_Vector[eqs_num * k + 2];
		}
		else
		{
			eqs_new = EQS_Vector[eqs_num * k];
		}
#endif
	} // WW 02.2013. Pardiso
#else
	// WW  phase=1;
	// CRFProcess *m_pcs = NULL;                      //
	// create EQS
	/// Configure EQS (old matrx) . WW 06.2011
	if (getProcessType() == FiniteElement::DEFORMATION || getProcessType() == FiniteElement::DEFORMATION_FLOW
	    || getProcessType() == FiniteElement::DEFORMATION_H2)
	{
		if (getProcessType() == FiniteElement::DEFORMATION)
			eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF, DOF * m_msh->GetNodesNumber(true));
		else if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
		{
			if (num_type_name.find("EXCAVATION") != string::npos)
				eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF - 1, DOF * m_msh->GetNodesNumber(true));
			else
				eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF,
				                            (DOF - 1) * m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false));
		}
		else if (getProcessType() == FiniteElement::DEFORMATION_H2)
			if (m_num->nls_method == 1)
				eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF,
				                            (DOF - 2) * m_msh->GetNodesNumber(true) + 2 * m_msh->GetNodesNumber(false));

		InitializeLinearSolver(eqs, m_num);
		PCS_Solver.push_back(eqs);
		size_unknowns = eqs->dim;
	}
	else
	{
		// If there is a solver existing. WW
		CRFProcess* m_pcs = NULL;
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			m_pcs = pcs_vector[i];
			if (m_pcs && m_pcs->eqs)
				//				if (m_pcs->_pcs_type_name.find("DEFORMATION") == string::npos)
				if (!isDeformationProcess(m_pcs->getProcessType()))
					break;
		}
		// If unique mesh
		if (m_pcs && m_pcs->eqs && (fem_msh_vector.size() == 1))
			eqs = m_pcs->eqs;
		else
		{
			eqs = CreateLinearSolver(m_num->ls_storage_method, m_msh->GetNodesNumber(false) * DOF);
			InitializeLinearSolver(eqs, m_num);
			PCS_Solver.push_back(eqs);
		}
		size_unknowns = eqs->dim; // WW
	}
#endif // If NEW_EQS
#endif // END: if not use PETSC
	// Set solver properties: EQS<->SOL
	// Internen Speicher allokieren
	// Speicher initialisieren

	//----------------------------------------------------------------------------
	// Time unit factor //WW
	std::cout << "->Create TIM" << '\n';
	// CTimeDiscretization* Tim = TIMGet(_pcs_type_name);
	Tim = TIMGet(pcs_type_name);
	if (!Tim)
	{
		// 21.08.2008. WW
		/* JT->WW: It doesn't seem like a good idea to give a non-existent Tim the properties of some specified [0]
		vector.
		           Why not set default values, and then let other "Tim" control the stepping?
		           In other words. If HEAT_TRANSPORT doesn't have time control,
		           we cannot assign a time control type for a FLOW process to a HEAT process, this could give incorrect
		results.
		           THE DEFAULTS ARE NOW SET UP SUCH THAT... if "Tim" doesn't exist, this process has no influence on the
		time step.
		Tim = new CTimeDiscretization(*time_vector[0], pcs_type_name);
		*/
		Tim = new CTimeDiscretization();
		Tim->pcs_type_name = pcs_type_name;
		time_vector.push_back(Tim); // 21.08.2008. WW
	}
	if (Tim->time_control_type == TimeControlType::INVALID && Tim->time_step_vector.size() > 0)
		Tim->time_control_type = TimeControlType::FIXED_STEPS;
	//
	if (Tim->time_unit.find("MINUTE") != std::string::npos)
		time_unit_factor = 60.0;
	else if (Tim->time_unit.find("HOUR") != std::string::npos)
		time_unit_factor = 3600.0;
	else if (Tim->time_unit.find("DAY") != std::string::npos)
		time_unit_factor = 86400.0;
	else if (Tim->time_unit.find("MONTH") != std::string::npos)
		time_unit_factor = 2592000.0;
	else if (Tim->time_unit.find("YEAR") != std::string::npos)
		time_unit_factor = 31536000;

	//
	if (type == 4 || type / 10 == 4)
		m_msh->SwitchOnQuadraticNodes(true);
	else
		m_msh->SwitchOnQuadraticNodes(false);

	// ELE - config and create element values
	cout << "->Config ELE values" << '\n';
	AllocateMemGPoint();

	// ELE - config element matrices
	// NOD - config and create node values
	cout << "->Config NOD values" << '\n';
	double* nod_values = NULL;
	double* ele_values = NULL; // PCH

	number_of_nvals = 2 * DOF + pcs_number_of_secondary_nvals;
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		// new time
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
		// old time //need this MB!
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
	}
	for (int i = 0; i < pcs_number_of_secondary_nvals; i++)
		// new time
		nod_val_name_vector.push_back(pcs_secondary_function_name[i]);

	long m_msh_nod_vector_size = m_msh->NodesNumber_Quadratic;
	for (long j = 0; j < number_of_nvals; j++) // Swap number_of_nvals and mesh size. WW 19.12.2012
	{
		nod_values = new double[m_msh_nod_vector_size];
		for (int i = 0; i < m_msh_nod_vector_size; i++)
			nod_values[i] = 0.0;
		nod_val_vector.push_back(nod_values);
	}
	// Create element values - PCH
	int number_of_evals = 2 * pcs_number_of_evals; // PCH, increase memory
	if (number_of_evals > 0) // WW added this "if" condition
	{
		for (int i = 0; i < pcs_number_of_evals; i++)
		{
			// new time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
			// old time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
		}
		size_t m_msh_ele_vector_size(m_msh->ele_vector.size());
		if (ele_val_vector.size() == 0)
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = new double[number_of_evals];
				size_eval += number_of_evals; // WW
				for (int i = 0; i < number_of_evals; i++)
					ele_values[i] = 0.0;
				ele_val_vector.push_back(ele_values);
			}
		else
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = ele_val_vector[j];
				ele_values = resize(ele_values, size_eval, size_eval + number_of_evals);
				size_eval += number_of_evals;
				ele_val_vector[j] = ele_values;
			}
	}
	//
	if (reload >= 2 && type != 4 && type / 10 != 4) // Modified at 03.08.2010. WW
	{
		// PCH
		cout << "Reloading the primary variables... "
		     << "\n";
		ReadSolution(); // WW
	}

	if (reload < 2) // PCH: If reload is set, no need to have ICs
	{
		// IC
		cout << "->Assign IC" << '\n';
		SetIC();
	}
	else
		// Bypassing IC
		cout << "RELOAD is set to be " << reload << ". So, bypassing IC's"
		     << "\n";

	if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos) // WW
		setIC_danymic_problems();

	// Keep all local matrices in the memory
	if (type != 55) // Not for fluid momentum. WW
	{
		if (Memory_Type != 0)
			AllocateLocalMatrixMemory();
		if (type == 4 || type / 10 == 4)
		{
			// Set initialization function
			//      CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *) (this);
			CRFProcessDeformation* dm_pcs = static_cast<CRFProcessDeformation*>(this);
			dm_pcs->Initialization();
		}
		else // Initialize FEM calculator
		{
			int Axisymm = 1; // ani-axisymmetry
			if (m_msh->isAxisymmetry())
				Axisymm = -1; // Axisymmetry is true
			fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
			fem->SetGaussPointNumber(m_num->ele_gauss_points);
		}
	}

	// Initialize the system equations
	if (PCSSetIC_USER)
		PCSSetIC_USER(pcs_type_number);

	/// Variable index for equation. 20.08.2010. WW
	if (p_var_index)
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
			p_var_index[i] = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	size_unknowns = m_msh->NodesNumber_Quadratic * pcs_number_of_primary_nvals;
#elif defined(NEW_EQS)
/// For JFNK. 01.10.2010. WW
#ifdef JFNK_H2M
	if (m_num->nls_method == 2)
	{
		size_unknowns = eqs_new->size_global;
		array_u_JFNK = new double[eqs_new->size_global];
		array_Fu_JFNK = new double[eqs_new->size_global];
	}
	else
#endif
	{
#ifdef USE_MPI
		size_unknowns = eqs_new->size_global;
#else
		size_unknowns = eqs_new->A->Dim();
#endif
	}
#endif
}

void initializeConstrainedProcesses(std::vector<CRFProcess*>& pcs_vector)
{
	// set bool for existing constrained BCs
	for (std::size_t i = 0; i < pcs_vector.size(); i++)
	{
		for (std::size_t j = 0; j < pcs_vector[i]->bc_node.size(); j++)
		{
			if (pcs_vector[i]->bc_node[j]->isConstrainedBC())
			{
				pcs_vector[i]->hasConstrainedBC(true);
				break;
			}
		}

		for (std::size_t j = 0; j < pcs_vector[i]->st_node.size(); j++)
		{
			if (pcs_vector[i]->st_node[j]->isConstrainedST())
			{
				pcs_vector[i]->hasConstrainedST(true);
				break;
			}
		}
	}

	// get the indices of velocity of flow process if contrained BC
	for (std::size_t i = 0; i < pcs_vector.size(); i++)
	{
		if (!(pcs_vector[i]->hasConstrainedBC()))
			continue;

		bool found(false);
		for (std::size_t j = 0; j < pcs_vector[i]->bc_node.size(); j++)
		{
			if (found)
				break;
			for (std::size_t k = 0; k < pcs_vector[i]->bc_node[j]->getNumberOfConstrainedBCs(); k++)
			{
				if (found)
					break;
				Constrained tmp(pcs_vector[i]->bc_node[j]->getConstrainedBC(k));
				if (tmp.constrainedVariable == ConstrainedVariable::VELOCITY)
				{
					CRFProcess* pcs = PCSGetFlow();
					pcs_vector[i]->setidxVx(pcs->GetNodeValueIndex("VELOCITY_X1", true));
					pcs_vector[i]->setidxVy(pcs->GetNodeValueIndex("VELOCITY_Y1", true));
					pcs_vector[i]->setidxVz(pcs->GetNodeValueIndex("VELOCITY_Z1", true));
					// jump out of j & k loop
					found = true;
				}
			}
		}
	}
}

void CRFProcess::SetBoundaryConditionAndSourceTerm()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));

	if (!pcs_type_name_vector.empty() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos) // WW
	{
		setBC_danymic_problems();
		setST_danymic_problems();
	}
	else
	{
		const int DOF = GetPrimaryVNumber();
		// BC - create BC groups for each process
		ScreenMessage("-> Create BC\n");
		CBoundaryConditionsGroup* m_bc_group = NULL;

		// 25.08.2011. WW
		if (WriteProcessed_BC == 2)
			Read_Processed_BC();
		else
		{
			for (int i = 0; i < DOF; i++)
			{
				// OKm_bc_group = BCGetGroup(_pcs_type_name,pcs_primary_function_name[i]);
				// OKif(!m_bc_group){
				BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
				m_bc_group = new CBoundaryConditionsGroup();
				// OK
				m_bc_group->setProcessTypeName(pcs_type_name);
				m_bc_group->setProcessPrimaryVariableName(pcs_primary_function_name[i]); // OK
				m_bc_group->Set(this, Shift[i]);

				bc_group_list.push_back(m_bc_group); // Useless, to be removed. WW
				m_bc_group = NULL;
				// OK}
			}
#ifndef USE_PETSC
			if (bc_node_value.size() < 1) // WW
				cout << "Warning: no boundary conditions specified for " << pcs_type_name << endl;
#endif
			if (WriteProcessed_BC == 1)
				Write_Processed_BC();
		}
		// ST - create ST groups for each process
		ScreenMessage("-> Create ST\n");
		CSourceTermGroup* m_st_group = NULL;

		if (WriteSourceNBC_RHS == 2) // Read from file
			ReadRHS_of_ST_NeumannBC();
		else // WW // Calculate directly
		{
			for (int i = 0; i < DOF; i++)
			{
				// OK m_st_group = m_st_group->Get(pcs_primary_function_name[i]);
				m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);

				if (!m_st_group)
				{
					m_st_group = new CSourceTermGroup();
					// OK
					m_st_group->pcs_type_name = pcs_type_name;
					// OK
					m_st_group->pcs_pv_name = pcs_primary_function_name[i];
					m_st_group->Set(this, Shift[i]);
					// Useless, to be removed. WW
					st_group_list.push_back(m_st_group);
				}
			}
			if (WriteSourceNBC_RHS == 1) // WW
				WriteRHS_of_ST_NeumannBC();
		}
		m_st_group = NULL;
	}
	// Write BC/ST nodes for vsualization.WW
	if (write_boundary_condition && WriteSourceNBC_RHS != 2)
		WriteBC();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   12/2005 WW
   03/2006 WW Write as acsi
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::WriteRHS_of_ST_NeumannBC()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}

	os << "$PCS_TYPE  "
	   << "\n";

	os << pcs_type_name << "\n";
	os << "geo_node_number  ";
	os << "msh_node_number  ";
	os << "CurveIndex ";
	os << "node_value ";
	os << "\n";
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(14);
	const size_t st_node_value_size(st_node_value.size());
	os << st_node_value_size << "\n";
	for (size_t i = 0; i < st_node_value_size; i++)
		st_node_value[i]->Write(os);
	os.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   03/2006 WW
   last modified: 04/2006
**************************************************************************/
void CRFProcess::ReadRHS_of_ST_NeumannBC()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
	std::ifstream is(m_file_name.c_str(), std::ios::in);
	if (!is.good())
	{
		cout << "File " << m_file_name << " is not found"
		     << "\n";
		abort();
	}

	std::string s_buffer;
	getline(is, s_buffer);
	getline(is, s_buffer);
	getline(is, s_buffer);
	size_t size;
	is >> size >> ws;
	st_node_value.clear();
	for (size_t i = 0; i < size; i++)
	{
		CNodeValue* cnodev = new CNodeValue();
		cnodev->Read(is);
		st_node_value.push_back(cnodev);
	}
	is.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   08/2011 WW
**************************************************************************/
void CRFProcess::Read_Processed_BC()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_eqs_BC.asc";
	std::ifstream is(m_file_name.c_str(), std::ios::in);
	if (!is.good())
	{
		cout << "File " << m_file_name << " is not found"
		     << "\n";
		abort();
	}

	std::string s_buffer;
	getline(is, s_buffer);
	getline(is, s_buffer);
	getline(is, s_buffer);
	size_t size;
	is >> size >> ws;
	bc_node_value.clear();
	for (size_t i = 0; i < size; i++)
	{
		CBoundaryConditionNode* cnodev = new CBoundaryConditionNode();
		cnodev->Read(is);
		bc_node_value.push_back(cnodev);
	}
	is.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   08/2011 WW
**************************************************************************/
void CRFProcess::Write_Processed_BC()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_eqs_BC.asc";
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}

	os << "$PCS_TYPE  "
	   << "\n";

	os << pcs_type_name << "\n";
	os << "geo_node_number  ";
	os << "msh_node_number  ";
	os << "CurveIndex ";
	os << "node_value ";
	os << "\n";
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(14);
	const size_t bc_node_value_size(bc_node_value.size());
	os << bc_node_value_size << "\n";
	for (size_t i = 0; i < bc_node_value_size; i++)
		bc_node_value[i]->Write(os);
	os.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the solution
   Programing:
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::WriteSolution()
{
	if (reload == 2 || reload <= 0)
		return;
	// kg44 write out only between nwrite_restart timesteps
	if ((aktueller_zeitschritt % nwrite_restart) > 0)
		return;

	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
#if defined(USE_PETSC) //|| defined(other parallel libs)//03.3012. WW
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::string m_file_name = FileName + "_" + pcs_type_name + "_" + pcs_primary_function_name[0] + "_primary_value_"
	                          + number2str(rank) + ".asc";

#else
	std::string m_file_name = FileName + "_" + pcs_type_name + "_" + pcs_primary_function_name[0]
	                          + "_primary_value.asc";
#endif
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}

	os.precision(15); // 15 digits accuracy seems enough? more fields are filled up with random numbers!
	os.setf(std::ios_base::scientific, std::ios_base::floatfield);

	int j;
	int* idx(new int[2 * pcs_number_of_primary_nvals]);
	for (j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			os << GetNodeValue(i, idx[j]) << "  ";
		os << "\n";
	}
	os.close();
	cout << "Write solutions for timestep " << aktueller_zeitschritt << " into file " << m_file_name << "\n";
	delete[] idx;
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the solution
   Programing:
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::ReadSolution()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
#if defined(USE_PETSC)
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::string m_file_name = FileName + "_" + pcs_type_name + "_" + pcs_primary_function_name[0] + "_primary_value_"
	                          + number2str(rank) + ".asc";

#else
	std::string m_file_name = FileName + "_" + pcs_type_name + "_" + pcs_primary_function_name[0]
	                          + "_primary_value.asc";
#endif
	std::ifstream is(m_file_name.c_str(), ios::in);
	if (!is.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}
	int j;

	int* idx(new int[2 * pcs_number_of_primary_nvals]);
	double* val(new double[2 * pcs_number_of_primary_nvals]);

	for (j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			is >> val[j];
		is >> ws;
		for (int j = 0; j < pcs_number_of_primary_nvals; j++)
		{
			SetNodeValue(i, idx[j], val[j + pcs_number_of_primary_nvals]);
			SetNodeValue(i, idx[j + pcs_number_of_primary_nvals], val[j + pcs_number_of_primary_nvals]);
		}
	}
	is.close();
	delete[] idx;
	delete[] val;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setIC_danymic_problems()
{
	const char* function_name[7];
	int i, j, nv;
	nv = 0;
	if (max_dim == 1) // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	}
	else // 3D
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	CInitialCondition* m_ic = NULL;
	long no_ics = (long)ic_vector.size();
	int nidx;
	for (i = 0; i < nv; i++)
	{
		nidx = GetNodeValueIndex(function_name[i]);
		for (j = 0; j < no_ics; j++)
		{
			m_ic = ic_vector[j];
			if (m_ic->getProcessPrimaryVariable() == FiniteElement::convertPrimaryVariable(function_name[i]))
				m_ic->Set(nidx);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setST_danymic_problems()
{
	const char* function_name[7];
	size_t nv = 0;
	if (max_dim == 1) // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	} // 3D
	else
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	// ST - create ST groups for each process
	CSourceTermGroup* m_st_group = NULL;
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	for (size_t i = 0; i < nv; i++)
	{
		m_st_group = STGetGroup(pcs_type_name, function_name[i]);
		if (!m_st_group)
		{
			m_st_group = new CSourceTermGroup();
			m_st_group->pcs_type_name = pcs_type_name;
			m_st_group->pcs_pv_name = function_name[i];
			m_st_group->Set(this, Shift[i], function_name[i]);
			st_group_list.push_back(m_st_group); // Useless, to be removed. WW
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setBC_danymic_problems()
{
	const char* function_name[7];
	size_t nv = 0;
	if (max_dim == 1) // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	} // 3D
	else
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	cout << "->Create BC" << '\n';
	CBoundaryConditionsGroup* m_bc_group = NULL;
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	for (size_t i = 0; i < nv; i++)
	{
		BCGroupDelete(pcs_type_name, function_name[i]);
		m_bc_group = new CBoundaryConditionsGroup();
		// OK
		m_bc_group->setProcessTypeName(pcs_type_name);
		// OK
		m_bc_group->setProcessPrimaryVariableName(function_name[i]);
		m_bc_group->Set(this, Shift[i], function_name[i]);
		bc_group_list.push_back(m_bc_group); // Useless, to be removed. WW
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::ConfigureCouplingForLocalAssemblier()
{
	bool Dyn = false;
	if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
		Dyn = true;
	if (fem)
		fem->ConfigureCoupling(this, Shift, Dyn);
}

/**************************************************************************
   FEMLib-Method:
   06/2003 OK Implementation
        WW 2nd version, PCS_Solver
**************************************************************************/
void PCSDestroyAllProcesses(void)
{
	CRFProcess* m_process = NULL;
	long i;
	int j;
//----------------------------------------------------------------------
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
// SOLver
#ifdef NEW_EQS // WW
#if defined(USE_MPI)
	for (j = 0; j < (int)EQS_Vector.size(); j += 2) // WW
#else
	for (j = 0; j < (int)EQS_Vector.size(); j++) // WW
#endif
	{
		if (EQS_Vector[j])
			delete EQS_Vector[j];
		EQS_Vector[j] = NULL;
#if defined(USE_MPI)
		EQS_Vector[j + 1] = NULL;
#endif
	}
#else // ifdef NEW_EQS
	// SOLDelete()
	LINEAR_SOLVER* eqs;
	for (j = 0; j < (int)PCS_Solver.size(); j++)
	{
		eqs = PCS_Solver[j];
		if (eqs->unknown_vector_indeces)
			eqs->unknown_vector_indeces = (int*)Free(eqs->unknown_vector_indeces);
		if (eqs->unknown_node_numbers)
			eqs->unknown_node_numbers = (long*)Free(eqs->unknown_node_numbers);
		if (eqs->unknown_update_methods)
			eqs->unknown_update_methods = (int*)Free(eqs->unknown_update_methods);
		eqs = DestroyLinearSolver(eqs);
	}
	PCS_Solver.clear(); // WW
#endif
//------
#endif //#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW

	//----------------------------------------------------------------------
	// PCS
	for (j = 0; j < (int)pcs_vector.size(); j++)
	{
		m_process = pcs_vector[j];
#ifdef USE_MPI // WW
		//  if(myrank==0)
		m_process->Print_CPU_time_byAssembly();
#endif
		if (m_process->pcs_nval_data)
			m_process->pcs_nval_data = (PCS_NVAL_DATA*)Free(m_process->pcs_nval_data);
		if (m_process->pcs_eval_data)
			m_process->pcs_eval_data = (PCS_EVAL_DATA*)Free(m_process->pcs_eval_data);
#ifdef PCS_NOD
		for (i = 0; i < NodeListSize(); i++)
		{
			k = GetNode(i);
			k->values[m_process->pcs_number] = (double*)Free(k->values[m_process->pcs_number]);
		}
#endif
		if (m_process->TempArry) // MX
			m_process->TempArry = (double*)Free(m_process->TempArry);
		delete (m_process);
	}
	pcs_vector.clear(); // WW
	//----------------------------------------------------------------------
	// MSH
	for (i = 0; i < (long)fem_msh_vector.size(); i++)
	{
		if (fem_msh_vector[i])
			delete fem_msh_vector[i];
		fem_msh_vector[i] = NULL;
	}
	fem_msh_vector.clear();
//----------------------------------------------------------------------

// DOM WW
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
#if defined(USE_MPI)
	// if(myrank==0)
	dom_vector[myrank]->PrintEQS_CPUtime(); // WW
#endif
	for (i = 0; i < (long)dom_vector.size(); i++)
	{
		if (dom_vector[i])
			delete dom_vector[i];
		dom_vector[i] = NULL;
	}
	dom_vector.clear();
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	//----------------------------------------------------------------------
	// ELE
	for (i = 0; i < (long)ele_val_vector.size(); i++)
		delete ele_val_vector[i];
	ele_val_vector.clear();
	//----------------------------------------------------------------------
	// IC ICDelete()
	for (i = 0; i < (long)ic_vector.size(); i++)
		delete ic_vector[i];
	ic_vector.clear();
	//----------------------------------------------------------------------
	MSPDelete(); // WW
	BCDelete(); // WW
	ICDelete(); // HS
	BCGroupDelete(); // HS
	STDelete(); // WW
	STGroupsDelete(); // HS
	GEOLIB_Clear_GeoLib_Data(); // HS
	//......................................................................
	TIMDelete(); // OK
	OUTDelete();
	NUMDelete();
	MFPDelete();
	MSPDelete();
	MMPDelete();
	MMPGroupDelete();
	MCPDelete();
	//----------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   last modified:
   10/2010 TF changed process type handling from string to enum
**************************************************************************/
bool PCSRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	char line[MAX_ZEILE];
	int indexCh1a, indexCh2a;
	std::string CommentK("//");
	std::string line_string;
	std::string pcs_file_name;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
	std::ifstream pcs_file(pcs_file_name.data(), ios::in);
	if (!pcs_file.good())
	{
		cout << "Warning: no PCS data *.pcs file is missing"
		     << "\n";
		return false;
	}

	// rewind the file
	pcs_file.clear();
	pcs_file.seekg(0, std::ios::beg);
	//========================================================================
	// Keyword loop
	std::cout << "PCSRead ... " << std::flush;
	while (!pcs_file.eof())
	{
		pcs_file.getline(line, MAX_ZEILE);
		line_string = line;
		line_string = GetUncommentedLine(line_string);
		if (line_string.find("#STOP") != string::npos)
			break;
		indexCh1a = (int)line_string.find_first_of(CommentK.c_str());
		indexCh2a = (int)line_string.find("#PROCESS");
		//----------------------------------------------------------------------
		// keyword found
		if (indexCh2a > indexCh1a && (indexCh1a == -1))
		{
			CRFProcess* m_pcs = new CRFProcess();
			m_pcs->file_name_base = file_base_name; // OK
			position = m_pcs->Read(&pcs_file);
			m_pcs->PCSReadConfigurations(); // JT

			m_pcs->pcs_number = pcs_vector.size();

			// RelocateDeformationProcess(m_pcs);
			//			if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos) { // TF
			if (isDeformationProcess(m_pcs->getProcessType()))
			{
				pcs_vector.push_back(m_pcs->CopyPCStoDM_PCS());
				pcs_vector[pcs_vector.size() - 1]->pcs_number = pcs_vector.size();
				delete m_pcs;
			}
			else
			{
				pcs_vector.push_back(m_pcs);
			}

			pcs_file.seekg(position, std::ios::beg);
		} // keyword found
	} // eof

	std::cout << "done, read " << pcs_vector.size() << " processes"
	          << "\n";

	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: Copy data to dm_pcs from PCS read function
   Programing:
   06/2007 OK/WW Implementation
   10/2010 TF many improvements
**************************************************************************/
CRFProcess* CRFProcess::CopyPCStoDM_PCS()
{
	// Numerics
	if (num_type_name.compare("STRONG_DISCONTINUITY") == 0)
		enhanced_strain_dm = 1;

	CRFProcessDeformation* dm_pcs(new CRFProcessDeformation());
	dm_pcs->setProcessType(this->getProcessType());
	dm_pcs->pcs_type_name_vector.push_back(pcs_type_name_vector[0].data());
	dm_pcs->Write_Matrix = Write_Matrix;
	dm_pcs->WriteSourceNBC_RHS = WriteSourceNBC_RHS;
	dm_pcs->num_type_name = num_type_name;
	dm_pcs->Memory_Type = Memory_Type;
	dm_pcs->NumDeactivated_SubDomains = NumDeactivated_SubDomains;
	dm_pcs->reload = reload;
	dm_pcs->nwrite_restart = nwrite_restart;
	dm_pcs->isPCSDeformation = true;
	dm_pcs->isPCSFlow = this->isPCSFlow; // JT
	dm_pcs->isPCSMultiFlow = this->isPCSMultiFlow; // JT
	// WW
	dm_pcs->write_boundary_condition = write_boundary_condition;
	if (!dm_pcs->Deactivated_SubDomain)
		dm_pcs->Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
	for (int i = 0; i < NumDeactivated_SubDomains; i++)
		dm_pcs->Deactivated_SubDomain[i] = Deactivated_SubDomain[i];
	pcs_deformation = 1;
	// WX:01.2011 for coupled excavation
	if (ExcavMaterialGroup >= 0)
	{
		dm_pcs->ExcavMaterialGroup = ExcavMaterialGroup;
		dm_pcs->ExcavDirection = ExcavDirection;
		dm_pcs->ExcavBeginCoordinate = ExcavBeginCoordinate;
		dm_pcs->ExcavCurve = ExcavCurve;
	}
	dm_pcs->Neglect_H_ini = Neglect_H_ini; // WX:08.2011
	dm_pcs->UpdateIniState = UpdateIniState; // WX:10.2011
	//
	return dynamic_cast<CRFProcess*>(dm_pcs);
}

/**************************************************************************
   FEMLib-Method:
   Task: Initial needed configurations following a read of PCS
   Programing:
   03/2012 JT
**************************************************************************/
void CRFProcess::PCSReadConfigurations()
{
	if (pcs_type_name_vector.size() > 1)
	{
		string pname = pcs_type_name_vector[0] + pcs_type_name_vector[1];
		pcs_type_name_vector.pop_back();
		if (pname.find("FLOW") != string::npos && pname.find("DEFORMATION") != string::npos)
		{
			setProcessType(FiniteElement::DEFORMATION_FLOW);
			MH_Process = true; // MH monolithic scheme
			if (pname.find("DYNAMIC") != string::npos)
				pcs_type_name_vector[0] = "DYNAMIC";
		}
	}
	else if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		// NW
		std::cout << "***Error: DEFORMATION_FLOW is vague definition."
		          << "\n";
		exit(0);
	}

	if (isFlowProcess(getProcessType()))
	{
		this->isPCSFlow = true;
		H_Process = true;
	}
	if (isMultiFlowProcess(getProcessType()))
	{
		this->isPCSMultiFlow = true;
		H2_Process = true;
	}
	if (isDeformationProcess(getProcessType()))
	{
		this->isPCSDeformation = true;
		M_Process = true;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Read deformation process
   11/2004 OK file streaming
   12/2005 OK MSH_TYPE
   01/2006 OK GEO_TYPE
**************************************************************************/
std::ios::pos_type CRFProcess::Read(std::ifstream* pcs_file)
{
	char line[MAX_ZEILE];
	string line_string;
	string CommentK("//");
	string hash("#");
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line_stream;
	std::string pcs_type_name;

	saturation_switch = false; // JOD for Richards
	//----------------------------------------------------------------------
	while (!new_keyword)
	{
		position = pcs_file->tellg();
		pcs_file->getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// WW Comment line
		if (line_string.find_first_of(CommentK.c_str()) != string::npos)
			return position;
		// SB check for comment sign ;
		line_string = GetUncommentedLine(line_string);
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position = pcs_file->tellg();
				line_string = GetLineFromFile1(pcs_file);
				if (line_string.find("#") != string::npos)
					return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line_stream.str(line_string);
				line_stream >> pcs_type_name;
				pcs_type_name_vector.push_back(pcs_type_name);
				this->setProcessType(FiniteElement::convertProcessType(pcs_type_name));
				line_stream.clear();

				if (isFlowProcess(this->getProcessType()))
				{
					H_Process = true;
					this->isPCSFlow = true; // JT2012
					pcs_number_flow = pcs_vector.size(); // JT2012
					if (this->getProcessType() == FiniteElement::PS_GLOBAL
					    || this->getProcessType() == FiniteElement::MULTI_PHASE_FLOW
					    ||
					    // this->getProcessType() == FiniteElement::TNEQ || no multiphase flow
					    pcs_type_name.find("H2") != string::npos)
					{
						this->isPCSMultiFlow = true;
					}
				}
				if (isDeformationProcess(this->getProcessType()))
				{
					M_Process = true;
					this->isPCSDeformation = true; // JT2012
					// JT: "pcs_number_deformation" is set in CRFProcessDeformation::Initialization()
				}
				if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					H_Process = true;
					MASS_TRANSPORT_Process = true;
					this->isPCSMass = true; // JT2012
					pcs_number_mass[pcs_no_components] = pcs_vector.size(); // JT2012
					pcs_no_components++;
					this->setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
				}
				if (this->getProcessType() == FiniteElement::MULTI_COMPONENTIAL_FLOW)
				{
					MULTI_COMPONENTIAL_FLOW_Process = true;
				}
				if (this->getProcessType() == FiniteElement::HEAT_TRANSPORT)
				{
					T_Process = true;
					this->isPCSHeat = true; // JT2012
					pcs_number_heat = pcs_vector.size(); // JT2012
				}
				if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
				{
					FLUID_MOMENTUM_Process = true;
				}
				if (this->getProcessType() == FiniteElement::RANDOM_WALK)
				{
					RANDOM_WALK_Process = true;
				}
			}
		//....................................................................
		// subkeyword found
		if (line_string.find("$NUM_TYPE") != string::npos)
		{
			*pcs_file >> num_type_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CPL_TYPE") != string::npos)
		{
			*pcs_file >> cpl_type_name;
			if (cpl_type_name.compare("MONOLITHIC") == 0)
			{
				pcs_monolithic_flow = true;
				pcs_deformation = 11;
			}
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIM_TYPE") != string::npos)
		{
			std::string tim_type_name;
			*pcs_file >> tim_type_name;
			tim_type = convertTimType(tim_type_name);
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$APP_TYPE") != string::npos)
		{
			*pcs_file >> rwpt_app;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// YS 05/2013
		if (line_string.find("$COUNT") != string::npos)
		{
			*pcs_file >> rwpt_count;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
		{
			*pcs_file >> primary_variable_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ELEMENT_MATRIX_OUTPUT") != string::npos)
		{
			*pcs_file >> Write_Matrix; // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// WW
		if (line_string.find("$BOUNDARY_CONDITION_OUTPUT") != string::npos)
		{
			write_boundary_condition = true;
			continue;
		}
		//....................................................................
		// BG 05/2012
		if (line_string.find("$OutputMassOfGasInModel") != string::npos)
		{
			OutputMassOfGasInModel = true;
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ST_RHS") != string::npos)
		{
			*pcs_file >> WriteSourceNBC_RHS; // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		if (line_string.find("$PROCESSED_BC") != string::npos) // 25.08.2011. WW
		{
			*pcs_file >> WriteProcessed_BC;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$MEMORY_TYPE") != string::npos)
		{
			*pcs_file >> Memory_Type; // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RELOAD") != string::npos)
		{
			*pcs_file >> reload; // WW
			if (reload == 1 || reload == 3)
				*pcs_file >> nwrite_restart; // kg44 read number of timesteps between writing restart files
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		// subkeyword found
		if (line_string.find("$DEACTIVATED_SUBDOMAIN") != string::npos)
		{
			// WW
			*pcs_file >> NumDeactivated_SubDomains >> ws;
			Deactivated_SubDomain = new int[NumDeactivated_SubDomains];
			for (int i = 0; i < NumDeactivated_SubDomains; i++)
				*pcs_file >> Deactivated_SubDomain[i] >> ws;
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MSH_TYPE") != string::npos)
		{
			*pcs_file >> msh_type_name >> ws;
			continue;
		}
		//....................................................................
		//		if (line_string.find("$GEO_TYPE") != string::npos) { //OK
		//			*pcs_file >> geo_type >> geo_type_name >> ws;
		//			continue;
		//		}
		//
		//....................................................................
		// subkeyword found
		if (line_string.find("$MEDIUM_TYPE") != string::npos)
		{
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position_subkeyword = pcs_file->tellg();
				*pcs_file >> line_string;
				if (line_string.size() == 0)
					break;
				if (line_string.find("#") != string::npos)
				{
					new_keyword = true;
					break;
				}
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.find("CONTINUUM") != string::npos)
				{
					*pcs_file >> line_string;
					// WW
					double w_m = strtod(line_string.data(), NULL);
					continuum_vector.push_back(w_m);
					// WW
					continuum_vector.push_back(1.0 - w_m);
					break; // WW
				}
				pcs_file->ignore(MAX_ZEILE, '\n');
			}
			continue;
		}
		// OK
		if (line_string.find("$SATURATION_SWITCH") != string::npos)
		{
			saturation_switch = true;
			;
			continue;
		}
		// SB4900
		if (line_string.find("$USE_VELOCITIES_FOR_TRANSPORT") != string::npos)
		{
			// Only for fluid momentum process
			if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
				use_velocities_for_transport = true;
			continue;
		}
		// Interface to Eclipse and Dumux, BG, 09/2010
		//	if(line_string.find("$SIMULATOR")!=string::npos) { //OK
		if (line_string.compare("$SIMULATOR") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator;
			continue;
		}
		if (line_string.find("$SIMULATOR_PATH") == 0) // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if (line_string.find("$SIMULATOR_MODEL_PATH") == 0)
		{
			*pcs_file >> this->simulator_model_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if (line_string.find("$USE_PRECALCULATED_FILES") == 0)
		{
			this->PrecalculatedFiles = true;
			continue;
		}
		// WTP, 04/2013, coupling to Eclipse (save data files)
		if (line_string.find("$SAVE_ECLIPSE_DATA_FILES") == 0)
		{
			this->SaveEclipseDataFiles = true;
			continue;
		}
		// KB, 02/2011, coupling to Eclipse and DuMux
		if (line_string.find("$SIMULATOR_WELL_PATH") == 0)
		{
			*pcs_file >> this->simulator_well_path;
			continue;
		}
		// BG, NB 11/2010, calculating phase transition for CO2
		if (line_string.find("$PHASE_TRANSITION") == 0)
		{
			string tempstring;
			*pcs_file >> tempstring;
			if (tempstring == "CO2_H2O_NaCl")
				this->Phase_Transition_Model = 1;
			continue;
		}
		// SB redo WTP
		if (line_string.find("$DISSOLVED_CO2_PCS_NAME") == 0)
		{ // SB, CB 10/2011
			*pcs_file >> this->dissolved_co2_pcs_name;
			continue;
		}
		if (line_string.find("$DISSOLVED_CO2_INGAS_PCS_NAME") == 0)
		{ // SB, CB 10/2011
			*pcs_file >> this->dissolved_co2_ingas_pcs_name;
			continue;
		}

		// WX:07.2011
		if (line_string.find("$TIME_CONTROLLED_EXCAVATION") == 0)
		{
			*pcs_file >> ExcavMaterialGroup >> ExcavDirection >> ExcavBeginCoordinate >> ExcavCurve;
			continue;
		}
		if (line_string.find("$NEGLECT_H_INI_EFFECT") == 0) // WX:10.2011
		{
			*pcs_file >> Neglect_H_ini;
			continue;
		}
		if (line_string.find("$UPDATE_INI_STATE") == 0) // WX:10.2011
		{
			*pcs_file >> UpdateIniState;
			continue;
		}
		if (line_string.find("$CONSTANT") == 0)
		{
			*pcs_file >> _pcs_constant_model;
			switch (_pcs_constant_model)
			{
				case 0:
					break;
				case 1:
					*pcs_file >> _pcs_constant_value;
					break;
				case 2:
					*pcs_file >> _pcs_constant_curve;
					break;
				default:
					std::cout << "-> invalid CONSTANT model " << _pcs_constant_model << std::endl;
					break;
			}
			if (_pcs_constant_model > 0)
				std::cout << "-> CONSTAT is activated." << std::endl;
			continue;
		}
		//....................................................................
	}
	//----------------------------------------------------------------------
	return position;
}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   06/2009 OK Write only if existing
**************************************************************************/
void PCSWrite(string file_base_name)
{
	if ((int)pcs_vector.size() < 1)
		return;
	//----------------------------------------------------------------------
	// File handling
	string pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
	fstream pcs_file(pcs_file_name.data(), ios::trunc | ios::out);
	pcs_file.clear();
	//----------------------------------------------------------------------
	// PCS loop
	cout << "PCSWrite"
	     << "\n";
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		pcs_file << "#PROCESS"
		         << "\n";
		m_pcs->Write(&pcs_file);
	}
	//----------------------------------------------------------------------
	pcs_file << "#STOP"
	         << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   12/2005 OK MSH_TYPE
   last modified:
**************************************************************************/
void CRFProcess::Write(std::fstream* pcs_file)
{
	*pcs_file << " $PCS_TYPE"
	          << "\n";
	*pcs_file << "  " << convertProcessTypeToString(this->getProcessType()) << "\n";

	*pcs_file << " $NUM_TYPE"
	          << "\n";
	*pcs_file << "  " << num_type_name << "\n";

	*pcs_file << " $CPL_TYPE"
	          << "\n";
	*pcs_file << "  " << cpl_type_name << "\n";

	*pcs_file << " $TIM_TYPE"
	          << "\n";
	*pcs_file << "  " << convertTimTypeToString(tim_type) << "\n";

	if (msh_type_name.size() > 0)
	{
		*pcs_file << " $MSH_TYPE"
		          << "\n";
		*pcs_file << "  " << msh_type_name << "\n";
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2004 OK Implementation
   last modified:
   10/2010 changed access to process type
**************************************************************************/
CRFProcess* PCSGet(const std::string& pcs_type_name)
{
	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(pcs_type_name));
	for (size_t i = 0; i < pcs_vector.size(); i++)
		//		m_pcs = pcs_vector[i];
		//		if(m_pcs->pcs_type_name.compare(pcs_type_name)==0) { TF
		if (pcs_vector[i]->getProcessType() == pcs_type)
			return pcs_vector[i];

	return NULL;
}

CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type)
{
	for (size_t i = 0; i < pcs_vector.size(); i++)
		if (pcs_vector[i]->getProcessType() == pcs_type)
			return pcs_vector[i];

	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2008 TK New Version with Primary Variable Comparision
   last modified:
   10/2010 TF changed access to process type
**************************************************************************/
CRFProcess* PCSGetNew(const string& pcs_type_name, const string& primary_variable_name)
{
	CRFProcess* m_pcs_return = NULL;

	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(pcs_type_name));

	int matches = 0;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* pcs = pcs_vector[i];
		//		if (pcs->pcs_type_name.compare(pcs_type_name) == 0) { TF
		if (pcs->getProcessType() == pcs_type)
		{
			for (size_t j = 0; j < pcs->GetPrimaryVNumber(); j++)
				if (primary_variable_name.compare(pcs->GetPrimaryVName(j)) == 0)
				{
					m_pcs_return = pcs;
					matches++;
					if (matches > 1)
						return NULL;
				}
		}
	}
	if (matches == 0)
		return NULL;
	else
		return m_pcs_return;
}

//////////////////////////////////////////////////////////////////////////
// Access
//////////////////////////////////////////////////////////////////////////

// OK->SB please try Get function
CRFProcess* CRFProcess::GetProcessByFunctionName(char* name)
{
	CRFProcess* m_process = NULL;
	/* Tests */
	if (!name)
		return m_process;
	int i;
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		m_process = pcs_vector[i];
		if (strcmp(StrUp(m_process->pcs_primary_function_name[0]), StrUp(name)) == 0)
			break;
	}
	return m_process;
}

// SB: new 3912
CRFProcess* CRFProcess::GetProcessByNumber(int number)
{
	CRFProcess* m_process = NULL;
	/* Tests */
	if (number < 1)
		return m_process;
	int i;
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		m_process = pcs_vector[i];
		if (m_process->pcs_number == number)
			break;
	}
	return m_process;
}

//////////////////////////////////////////////////////////////////////////
// Configuration
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2003 OK Implementation
   08/2003 WW Modified to fit monolithic scheme
   02/2005 OK Unsaturated flow (Richards model)
   02/2005 MB string
   05/2005 WW/DL Dymanic problem
   01/2006 YD Dual Richards
   OKToDo switch to char
   03/2009 PCH PS_GLOBAL
   get rid of type
**************************************************************************/
void CRFProcess::Config(void)
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	// Set mesh pointer to corresponding mesh
	m_msh = FEMGet(pcs_type_name);
	if (!m_msh)
	{
		cout << "Error in CRFProcess::Config - no MSH data"
		     << "\n";
		return;
	}

	if (hasAnyProcessDeactivatedSubdomains)
		CheckMarkedElement(); // WW

	if (continuum_vector.size() == 0) // YD
		continuum_vector.push_back(1.0);

	//	if (_pcs_type_name.compare("LIQUID_FLOW") == 0) {
	if (this->getProcessType() == FiniteElement::LIQUID_FLOW || this->getProcessType() == FiniteElement::FLUID_FLOW)
	{
		std::cout << "CRFProcess::Config LIQUID_FLOW"
		          << "\n";
		type = 1;
		ConfigLiquidFlow();
	}
	//	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0) {
	if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
	{
		type = 1;
		ConfigGroundwaterFlow();
	}
	//	if (_pcs_type_name.compare("RICHARDS_FLOW") == 0) {
	if (this->getProcessType() == FiniteElement::RICHARDS_FLOW)
	{
		if (continuum_vector.size() > 1)
		{
			RD_Process = true;
			type = 22;
		}
		else
		{
			type = 14;
		}
		ConfigUnsaturatedFlow();
	}
	//	if (_pcs_type_name.compare("OVERLAND_FLOW") == 0) {
	if (this->getProcessType() == FiniteElement::OVERLAND_FLOW)
	{
		type = 66;
		max_dim = 1;
		ConfigGroundwaterFlow();
	}
	//	if (_pcs_type_name.compare("AIR_FLOW") == 0) { //OK
	if (this->getProcessType() == FiniteElement::AIR_FLOW) // OK
	{
		type = 5;
		ConfigGasFlow();
	}
	//	if (_pcs_type_name.compare("TWO_PHASE_FLOW") == 0) {
	if (this->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
	{
		type = 12;
		ConfigMultiphaseFlow();
	}
	//	if (_pcs_type_name.compare("COMPONENTAL_FLOW") == 0) {
	//	if (COMPONENTAL_FLOW) {
	//		type = 11;
	//		ConfigNonIsothermalFlow();
	//	}
	//	if (_pcs_type_name.compare("HEAT_TRANSPORT") == 0) {
	if (this->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
		type = 3;
		ConfigHeatTransport();
	}
	//	if (_pcs_type_name.compare("MASS_TRANSPORT") == 0) {
	if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		type = 2;
		ConfigMassTransport();
	}
	//	if (_pcs_type_name.find("DEFORMATION") != string::npos)
	if (isDeformationProcess(getProcessType()))
	{
		ConfigDeformation();
	}
	//	if (_pcs_type_name.find("FLUID_MOMENTUM") != string::npos
	if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
	{
		type = 55; // WW
		ConfigFluidMomentum();
	}
	//	if (_pcs_type_name.find("RANDOM_WALK") != string::npos) {
	if (this->getProcessType() == FiniteElement::RANDOM_WALK)
	{
		type = 55; // WW
		ConfigRandomWalk();
	}
	//	if (_pcs_type_name.find("MULTI_PHASE_FLOW") != string::npos) {//24.02.2007 WW
	if (this->getProcessType() == FiniteElement::MULTI_PHASE_FLOW) // 24.02.2007 WW
	{
		type = 1212;
		ConfigMultiPhaseFlow();
	}
	//	if (_pcs_type_name.find("PS_GLOBAL") != string::npos) {//24.02.2007 WW
	if (this->getProcessType() == FiniteElement::PS_GLOBAL) // 24.02.2007 WW
	{
		type = 1313;
		ConfigPS_Global();
	}
	if (this->getProcessType() == FiniteElement::MULTI_COMPONENTIAL_FLOW) // 24.02.2007 WW
	{
		type = 1111;
		ConfigMULTI_COMPONENTIAL_FLOW();
	}
	if (this->getProcessType() == FiniteElement::TNEQ) // 24.02.2007 WW
	{
		type = 1414;
		ConfigTNEQ();
	}
	if (this->getProcessType() == FiniteElement::TES) // 24.02.2007 WW
	{
		type = 1415;
		ConfigTES();
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
   02/2005 MB Pressure version for LIQUID Flow
**************************************************************************/
void CRFProcess::ConfigLiquidFlow()
{
	// pcs_num_name[0] = "PRESSURE0";
	// pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;
	pcs_number_of_evals = 0;
	Def_Variable_LiquidFlow(); // NW

	// Output material parameters
	configMaterialParameters();
}

/**************************************************************************
   FEMLib-Method:
   03/2003 OK Implementation
        WW Splitted for processes
   02/2005 MB head version for GroundwaterFlow
   08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigGroundwaterFlow()
{
	pcs_num_name[0] = "HEAD";
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_HEAD";
	// NOD values
	pcs_number_of_primary_nvals = 1;
	pcs_primary_function_name[0] = "HEAD";
	pcs_primary_function_unit[0] = "m";
	// ELE values
	pcs_number_of_evals = 8;
	pcs_eval_name[0] = "VOLUME";
	pcs_eval_unit[0] = "m3";
	pcs_eval_name[1] = "VELOCITY1_X";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Y";
	pcs_eval_unit[2] = "m/s";
	pcs_eval_name[3] = "VELOCITY1_Z";
	pcs_eval_unit[3] = "m/s";
	pcs_eval_name[4] = "PERMEABILITY";
	pcs_eval_unit[4] = "m^2";
	pcs_eval_name[5] = "POROSITY";
	pcs_eval_unit[5] = "-";
	pcs_eval_name[6] = "PERMEABILITY_YY";
	pcs_eval_unit[6] = "m^2";
	pcs_eval_name[7] = "PERMEABILITY_ZZ";
	pcs_eval_unit[7] = "m^2";
	//----------------------------------------------------------------------
	// Secondary variables
	pcs_number_of_secondary_nvals = 5;
	pcs_secondary_function_name[0] = "FLUX";
	pcs_secondary_function_unit[0] = "m3/s";
	pcs_secondary_function_timelevel[0] = 1;
	pcs_secondary_function_name[1] = "WDEPTH";
	pcs_secondary_function_unit[1] = "m";
	pcs_secondary_function_timelevel[1] = 1;
	pcs_secondary_function_name[2] = "COUPLING"; // JOD
	pcs_secondary_function_unit[2] = "m/s";
	pcs_secondary_function_timelevel[2] = 0;
	pcs_secondary_function_name[3] = "COUPLING"; // JOD
	pcs_secondary_function_unit[3] = "m/s";
	pcs_secondary_function_timelevel[3] = 1;
	pcs_secondary_function_name[4] = "STORE"; // JOD  for subtiming 4.7.10
	pcs_secondary_function_unit[4] = "m";
	pcs_secondary_function_timelevel[4] = 1;

	pcs_number_of_secondary_nvals = 5; // WW
	// WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++; // WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++; // WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++; // WW
	//----------------------------------------------------------------------
	// WW / TF
	// Output material parameters
	configMaterialParameters();

	if (m_msh)
		m_msh->DefineMobileNodes(this);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   10/2004 OK Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigGasFlow()
{
	//----------------------------------------------------------------------
	// Primary variables - NOD values
	pcs_number_of_primary_nvals = 1;
	pcs_number_of_secondary_nvals = 0;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	//----------------------------------------------------------------------
	// Secondary variables - NOD values
	pcs_number_of_secondary_nvals = 1;
	pcs_secondary_function_name[0] = "NOD_MASS_FLUX";
	pcs_secondary_function_unit[0] = "kg/s";
	//----------------------------------------------------------------------
	// ELE values
	pcs_number_of_evals = 3;
	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
	//----------------------------------------------------------------------
	// NUM
	pcs_num_name[0] = "PRESSURE0";
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigMultiphaseFlow()
{
	switch (pcs_type_number)
	{
		case 0:
			pcs_num_name[0] = "PRESSURE0";
			pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
			break;
		case 1:
			pcs_num_name[0] = "SATURATION0";
			pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_SATURATION1";
			break;
	}
	//----------------------------------------------------------------------
	// NOD Primary variables
	pcs_number_of_primary_nvals = 1;
	switch (pcs_type_number)
	{
		case 0:
			pcs_primary_function_name[0] = "PRESSURE1";
			pcs_primary_function_unit[0] = "Pa";
			break;
		case 1:
			pcs_primary_function_name[0] = "SATURATION2";
			pcs_primary_function_unit[0] = "m3/m3";
			break;
	}
	//----------------------------------------------------------------------
	// NOD Secondary variables
	pcs_number_of_secondary_nvals = 12; // BG
	switch (pcs_type_number)
	{
		case 0:
			pcs_secondary_function_name[0] = "PRESSURE2";
			pcs_secondary_function_unit[0] = "Pa";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_name[1] = "PRESSURE2";
			pcs_secondary_function_unit[1] = "Pa";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_name[2] = "PRESSURE_CAP";
			pcs_secondary_function_unit[2] = "Pa";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_name[3] = "FLUX";
			pcs_secondary_function_unit[3] = "m3/s";
			pcs_secondary_function_timelevel[3] = 0;
			pcs_secondary_function_name[4] = "DENSITY1";
			pcs_secondary_function_unit[4] = "kg/m3";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_name[5] = "VISCOSITY1";
			pcs_secondary_function_unit[5] = "Pa s";
			pcs_secondary_function_timelevel[5] = 1;
			// BG
			pcs_secondary_function_name[6] = "VELOCITY_X1";
			pcs_secondary_function_unit[6] = "m/s";
			pcs_secondary_function_timelevel[6] = 1;
			pcs_secondary_function_name[7] = "VELOCITY_Y1";
			pcs_secondary_function_unit[7] = "m/s";
			pcs_secondary_function_timelevel[7] = 1;
			pcs_secondary_function_name[8] = "VELOCITY_Z1";
			pcs_secondary_function_unit[8] = "m/s";
			pcs_secondary_function_timelevel[8] = 1;
			pcs_secondary_function_name[9] = "VELOCITY_X2";
			pcs_secondary_function_unit[9] = "m/s";
			pcs_secondary_function_timelevel[9] = 1;
			pcs_secondary_function_name[10] = "VELOCITY_Y2";
			pcs_secondary_function_unit[10] = "m/s";
			pcs_secondary_function_timelevel[10] = 1;
			pcs_secondary_function_name[11] = "VELOCITY_Z2";
			pcs_secondary_function_unit[11] = "m/s";
			pcs_secondary_function_timelevel[11] = 1;
			break;
		case 1:
			pcs_secondary_function_name[0] = "SATURATION1";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "m3/m3";
			pcs_secondary_function_name[1] = "SATURATION1";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "m3/m3";
			pcs_secondary_function_name[2] = "PRESSURE_CAP";
			pcs_secondary_function_unit[2] = "Pa";
			pcs_secondary_function_timelevel[2] = 1;
			pcs_secondary_function_name[3] = "FLUX";
			pcs_secondary_function_unit[3] = "m3/s";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_name[4] = "DENSITY2";
			pcs_secondary_function_unit[4] = "kg/m3";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_name[5] = "VISCOSITY2";
			pcs_secondary_function_unit[5] = "Pa s";
			pcs_secondary_function_timelevel[5] = 1;
			// BG
			pcs_secondary_function_name[6] = "VELOCITY_X1";
			pcs_secondary_function_unit[6] = "m/s";
			pcs_secondary_function_timelevel[6] = 1;
			pcs_secondary_function_name[7] = "VELOCITY_Y1";
			pcs_secondary_function_unit[7] = "m/s";
			pcs_secondary_function_timelevel[7] = 1;
			pcs_secondary_function_name[8] = "VELOCITY_Z1";
			pcs_secondary_function_unit[8] = "m/s";
			pcs_secondary_function_timelevel[8] = 1;
			pcs_secondary_function_name[9] = "VELOCITY_X2";
			pcs_secondary_function_unit[9] = "m/s";
			pcs_secondary_function_timelevel[9] = 1;
			pcs_secondary_function_name[10] = "VELOCITY_Y2";
			pcs_secondary_function_unit[10] = "m/s";
			pcs_secondary_function_timelevel[10] = 1;
			pcs_secondary_function_name[11] = "VELOCITY_Z2";
			pcs_secondary_function_unit[11] = "m/s";
			pcs_secondary_function_timelevel[11] = 1;
			break;
	}
	//----------------------------------------------------------------------
	// ELE values
	pcs_number_of_evals = 7;
	switch (pcs_type_number)
	{
		case 0:
			pcs_eval_name[0] = "VELOCITY1_X";
			pcs_eval_unit[0] = "m/s";
			pcs_eval_name[1] = "VELOCITY1_Y";
			pcs_eval_unit[1] = "m/s";
			pcs_eval_name[2] = "VELOCITY1_Z";
			pcs_eval_unit[2] = "m/s";
			pcs_eval_name[3] = "POROSITY1"; // MX 03.2005
			pcs_eval_unit[3] = "-";
			pcs_eval_name[4] = "POROSITY1_IL"; // MX 03.2005
			pcs_eval_unit[4] = "-";
			pcs_eval_name[5] = "PERMEABILITY1"; // MX 03.2005
			pcs_eval_unit[5] = "-";
			pcs_eval_name[6] = "POROSITY1_SW"; // MX 03.2005
			pcs_eval_unit[6] = "-";
			break;
		case 1:
			pcs_eval_name[0] = "VELOCITY2_X";
			pcs_eval_unit[0] = "m/s";
			pcs_eval_name[1] = "VELOCITY2_Y";
			pcs_eval_unit[1] = "m/s";
			pcs_eval_name[2] = "VELOCITY2_Z";
			pcs_eval_unit[2] = "m/s";
			pcs_eval_name[3] = "POROSITY"; // MX 03.2005
			pcs_eval_unit[3] = "-";
			pcs_eval_name[4] = "POROSITY_IL"; // MX 03.2005
			pcs_eval_unit[4] = "-";
			pcs_eval_name[5] = "PERMEABILITY"; // MX 03.2005
			pcs_eval_unit[5] = "-";
			pcs_eval_name[6] = "POROSITY_SW"; // MX 03.2005
			pcs_eval_unit[6] = "-";
			break;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigNonIsothermalFlow()
{
	//----------------------------------------------------------------------
	// Primary variables
	pcs_number_of_primary_nvals = 1;
	switch (pcs_type_number)
	{
		case 0:
			pcs_num_name[0] = "PRESSURE0";
			pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
			pcs_primary_function_name[0] = "PRESSURE1";
			pcs_primary_function_unit[0] = "Pa";
			break;
		case 1:
			pcs_num_name[0] = "SATURATION0";
			pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_SATURATION1";
			pcs_primary_function_name[0] = "SATURATION2";
			pcs_primary_function_unit[0] = "m3/m3";
			break;
	}
	//----------------------------------------------------------------------
	// Secondary variables
	pcs_number_of_secondary_nvals = 6;
	switch (pcs_type_number)
	{
		case 0:
			pcs_secondary_function_name[0] = "PRESSURE2";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "Pa";
			pcs_secondary_function_name[1] = "PRESSURE2";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "Pa";
			pcs_secondary_function_name[2] = "MASS_FRACTION1";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_unit[2] = "kg/kg";
			pcs_secondary_function_name[3] = "MASS_FRACTION1";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_unit[3] = "kg/kg";
			pcs_secondary_function_name[4] = "PRESSURE_CAP";
			pcs_secondary_function_timelevel[4] = 0;
			pcs_secondary_function_unit[4] = "Pa";
			pcs_secondary_function_name[5] = "DENSITY1";
			pcs_secondary_function_timelevel[5] = 1;
			pcs_secondary_function_unit[5] = "kg/m3";
			break;
		case 1:
			pcs_secondary_function_name[0] = "SATURATION1";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "m3/m3";
			pcs_secondary_function_name[1] = "SATURATION1";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "m3/m3";
			pcs_secondary_function_name[2] = "MASS_FRACTION2";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_unit[2] = "kg/kg";
			pcs_secondary_function_name[3] = "MASS_FRACTION2";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_unit[3] = "kg/kg";
			pcs_secondary_function_name[4] = "PRESSURE_CAP";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_unit[4] = "Pa";
			pcs_secondary_function_name[5] = "DENSITY2";
			pcs_secondary_function_timelevel[5] = 1;
			pcs_secondary_function_unit[5] = "kg/m3";
			break;
	}
	// Node
	pcs_number_of_primary_nvals = 1;
	// ELE values
	pcs_number_of_evals = 14;
	pcs_eval_name[0] = "COMP_FLUX";
	pcs_eval_name[1] = "POROSITY";
	pcs_eval_name[2] = "PERMEABILITY";
	pcs_eval_name[3] = "VELOCITY1_X";
	pcs_eval_name[4] = "VELOCITY1_Y";
	pcs_eval_name[5] = "VELOCITY1_Z";
	pcs_eval_name[6] = "VELOCITY2_X";
	pcs_eval_name[7] = "VELOCITY2_Y";
	pcs_eval_name[8] = "VELOCITY2_Z";
	pcs_eval_name[9] = "POROSITY_IL";
	pcs_eval_name[10] = "VoidRatio";
	pcs_eval_name[11] = "PorosityChange";
	pcs_eval_name[12] = "n_sw_Rate";
	pcs_eval_name[13] = "POROSITY_SW";
	pcs_eval_unit[0] = "kg/s";
	pcs_eval_unit[1] = "m3/m3";
	pcs_eval_unit[2] = "m2";
	pcs_eval_unit[3] = "m/s";
	pcs_eval_unit[4] = "m/s";
	pcs_eval_unit[5] = "m/s";
	pcs_eval_unit[6] = "m/s";
	pcs_eval_unit[7] = "m/s";
	pcs_eval_unit[8] = "m/s";
	pcs_eval_unit[9] = "-";
	pcs_eval_unit[10] = "-";
	pcs_eval_unit[11] = "-";
	pcs_eval_unit[12] = "-";
	pcs_eval_unit[13] = "-";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
// void CRFProcess::ConfigNonIsothermalFlowRichards()

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2004 SB Implementation
        WW Splitted for processes
   01/2006 OK Tests
   08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigMassTransport()
{
	long comp = 1;
	/* count transport processes */
	pcs_component_number++;
	comp = pcs_component_number;
	// 1 NOD values
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 1;
	pcs_primary_function_name[0] = "";
	//  sprintf(pcs_primary_function_name[0], "%s%li","CONCENTRATION",comp);
	//----------------------------------------------------------------------
	// Tests
	// WW int size;
	// WW  size = (int)cp_vec.size();
	// int comb;                                      //OK411
	// comb = pcs_component_number;

	if ((int)cp_vec.size() < pcs_component_number + 1)
	{
		cout << "Error in CRFProcess::ConfigMassTransport - not enough MCP data" << '\n';
		return;
	}
	//----------------------------------------------------------------------
	pcs_primary_function_name[0] = cp_vec[pcs_component_number]->compname.c_str();
	// sprintf(pcs_primary_function_name[0], "%s", cp_vec[pcs_component_number]->compname.c_str());
	pcs_primary_function_unit[0] = "kg/m3"; // SB
	/* SB: Eintrag component name in Ausgabestruktur */ // SB:todo : just one phase todo:name
	/*
	      pcs_primary_function_name[0] = GetTracerCompName(0,this->pcs_component_number-1);
	      name_initial_condition_tracer_component = pcs_primary_function_name[0];
	      pcs_ic_name_mass = pcs_primary_function_name[0];
	 */
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 3; // SB3909
	pcs_secondary_function_name[0] = new char[80];
	char pcs_secondary_function_name_tmp[80];
	sprintf(pcs_secondary_function_name_tmp, "%s%li", "MASS_FLUX_", comp);
	strncpy((char*)pcs_secondary_function_name[0], pcs_secondary_function_name_tmp, 80);
	//      pcs_secondary_function_name[0] = "MASS_FLUX1";
	pcs_secondary_function_unit[0] = "kg/m3/s";
	pcs_secondary_function_timelevel[0] = 0;
	pcs_secondary_function_name[1] = new char[80];
	sprintf(pcs_secondary_function_name_tmp, "%s%li", "MASS_FLUX_", comp);
	strncpy((char*)pcs_secondary_function_name[1], pcs_secondary_function_name_tmp, 80);
	pcs_secondary_function_unit[1] = "kg/m3/s";
	pcs_secondary_function_timelevel[1] = 1;
	// KG44 added secondary function for adaptive time stepping
	string comp_name = "DELTA_" + convertPrimaryVariableToString(this->getProcessPrimaryVariable()); // JOD 2014-11-10
	pcs_secondary_function_name[2] = new char[80];
	strncpy((char*)pcs_secondary_function_name[2], comp_name.c_str(), 80);
	pcs_secondary_function_unit[2] = "kg/m3";
	pcs_secondary_function_timelevel[2] = 0;

	if (adaption)
	{
		pcs_number_of_secondary_nvals = 4;
		pcs_secondary_function_name[3] = new char[80];
		sprintf(pcs_secondary_function_name_tmp, "%s%li", "CONC_BACK_", comp);
		strncpy((char*)pcs_secondary_function_name[3], pcs_secondary_function_name_tmp, 80);
		pcs_secondary_function_unit[3] = "kg/m3";
		pcs_secondary_function_timelevel[3] = 0;
	}
	// OK  LOPCalcSecondaryVariables_USER = MTM2CalcSecondaryVariables;  //SB:todo
	// 2 ELE values
	pcs_number_of_evals = 0;
//	  pcs_eval_name[0] = "Darcy velocity";
#ifdef REACTION_ELEMENT
	pcs_number_of_evals = 1;
	pcs_eval_name[0] = pcs_primary_function_name[0];
	pcs_eval_unit[0] = "mol/kgH2O";
#endif
	// 3 ELE matrices
	// NUM
	pcs_num_name[0] = "CONCENTRATION0";
	/* SB: immer solver properties der ersten Komponente nehmen */
	// SB ??
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_CONCENTRATION1";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigHeatTransport()
{
	pcs_num_name[0] = "TEMPERATURE0";
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_TEMPERATURE1";
	// NOD
	if ((int)continuum_vector.size() == 1)
	{
		pcs_number_of_primary_nvals = 1;
		pcs_primary_function_name[0] = "TEMPERATURE1";
		pcs_primary_function_unit[0] = "K";
		pcs_number_of_secondary_nvals = 1; // JOD 2014-11-10
		pcs_secondary_function_name[0] = "DELTA_TEMPERATURE1";
		pcs_secondary_function_unit[0] = "K";
		pcs_secondary_function_timelevel[0] = 0;

#ifdef REACTION_ELEMENT
		pcs_number_of_evals = 1; // MX
		pcs_eval_name[0] = "TEMPERATURE1";
		pcs_eval_unit[0] = "K";
#endif
	}
	if ((int)continuum_vector.size() == 2)
	{
		pcs_number_of_primary_nvals = 2;
		pcs_primary_function_name[0] = "TEMPERATURE1";
		pcs_primary_function_unit[0] = "K";
		pcs_number_of_primary_nvals = 2;
		pcs_primary_function_name[1] = "TEMPERATURE2";
		pcs_primary_function_unit[1] = "K";
		pcs_number_of_secondary_nvals = 0;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 WW Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigDeformation()
{
// Generate high order nodes for all elements.
#if defined(USE_PETSC) // || defined (other parallel linear solver lib). //WW. 05.2013
	m_msh->ConfigHighOrderElements();
#else
	m_msh->GenerateHighOrderNodes(); // WW
#endif
	type = 4;

	//	if (_pcs_type_name.find("DEFORMATION") != string::npos
	//			&& _pcs_type_name.find("FLOW") != string::npos) {
	if (getProcessType() == FiniteElement::DEFORMATION_FLOW || getProcessType() == FiniteElement::DEFORMATION_H2)
	{
		type = 41;
		if (getProcessType() == FiniteElement::DEFORMATION_H2)
			type = 42;
		cpl_type_name = "MONOLITHIC";
		pcs_deformation = 11;
	}

	CNumerics* num = NULL;

	for (size_t ii = 0; ii < num_vector.size(); ii++)
	{
		num = num_vector[ii];
		if (num->pcs_type_name.find("DEFORMATION") != string::npos)
		{
			num->pcs_type_name = FiniteElement::convertProcessTypeToString(this->getProcessType());
			if (num->nls_method >= 1) // Newton-Raphson
			{
				pcs_deformation = 101;
				if (type / 10 == 4)
					pcs_deformation = 110;
			}
			break;
		}
	}

	// Prepare for restart
	// RFConfigRenumber();

	// NUM
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_DISPLACEMENT1";
	pcs_num_name[0] = "DISPLACEMENT0";
	pcs_num_name[1] = "PRESSURE0";
	if (pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
		VariableDynamics();
	else
		VariableStaticProblem();
	// OBJ names are set to PCS name
	// Geometry dimension
	problem_dimension_dm = m_msh->GetCoordinateFlag() / 10;
	problem_2d_plane_dm = 1;

	// Coupling
	int i;
	for (i = 0; i < problem_dimension_dm; i++)
		Shift[i] = i * m_msh->GetNodesNumber(true);

	/// 11-20.08.2010 WW
	long nn_H = (long)m_msh->GetNodesNumber(true);
	if (type == 4)
	{
		num_nodes_p_var = new long[problem_dimension_dm];
		p_var_index = new int[problem_dimension_dm];
		for (int i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
	}
	else if (type == 41)
	{
		num_nodes_p_var = new long[problem_dimension_dm + 1];
		p_var_index = new int[problem_dimension_dm + 1];
		for (i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
		num_nodes_p_var[problem_dimension_dm] = (long)m_msh->GetNodesNumber(false);
		Shift[problem_dimension_dm] = problem_dimension_dm * m_msh->GetNodesNumber(true);
	}
	else if (type == 42)
	{
		num_nodes_p_var = new long[problem_dimension_dm + 2];
		p_var_index = new int[problem_dimension_dm + 2];

		for (i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
		for (i = problem_dimension_dm; i < problem_dimension_dm + 2; i++)
			num_nodes_p_var[i] = (long)m_msh->GetNodesNumber(false);
		Shift[problem_dimension_dm] = problem_dimension_dm * m_msh->GetNodesNumber(true);
		Shift[problem_dimension_dm + 1] = Shift[problem_dimension_dm] + m_msh->GetNodesNumber(false);

#ifdef JFNK_H2M
		norm_u_JFNK = new double[2];
#endif
	}
}

/**************************************************************************
   FEMLib-Method: Static problems
   Task:
   Programing:
   05/2005 WW Implementation
   last modified:
**************************************************************************/
void CRFProcess::VariableStaticProblem()
{
	//----------------------------------------------------------------------
	// NOD Primary functions
	pcs_number_of_primary_nvals = 2; // OK distinguish 2/3D problems, problem_dimension_dm;
	dm_number_of_primary_nvals = 2;
	pcs_number_of_evals = 0;
	pcs_primary_function_name[0] = "DISPLACEMENT_X1";
	pcs_primary_function_name[1] = "DISPLACEMENT_Y1";
	pcs_primary_function_unit[0] = "m";
	pcs_primary_function_unit[1] = "m";
	if (max_dim == 2)
	{
		pcs_number_of_primary_nvals = 3;
		dm_number_of_primary_nvals = 3;
		pcs_primary_function_name[2] = "DISPLACEMENT_Z1";
		pcs_primary_function_unit[2] = "m";
	}
	//----------------------------------------------------------------------
	// NOD Secondary functions
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// NB new secondary variables for principle stresses; Principle stress directions might follow later
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_1"; // i.e. Sigma_1
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_2"; // i.e. Sigma_2
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_3"; // i.e. Sigma_3
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	//  pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "POROPRESSURE0";
	//  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	//  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	//  pcs_number_of_secondary_nvals++;

	if (max_dim == 2) // 3D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
	if (max_dim > 0) // >1D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}

	if (type == 41)
	{ // Monolithic scheme
		Def_Variable_LiquidFlow();

		// Output material parameters
		configMaterialParameters();
	}
	else if (type == 42) // Monolithic scheme H2M. 03.08.2010. WW
	{
		Def_Variable_MultiPhaseFlow();

		// Output material parameters
		configMaterialParameters();
	}
}

/**************************************************************************
   FEMLib-Method: Dynamic problems
   Task:
   Programing:
   05/2005 WW/LD Implementation
   last modified:
**************************************************************************/
void CRFProcess::VariableDynamics()
{
	//----------------------------------------------------------------------
	// NOD Primary functions
	pcs_number_of_primary_nvals = 2;
	dm_number_of_primary_nvals = 2;
	pcs_primary_function_name[0] = "ACCELERATION_X1";
	pcs_primary_function_name[1] = "ACCELERATION_Y1";
	pcs_primary_function_unit[0] = "m/s^2";
	pcs_primary_function_unit[1] = "m/s^2";
	if (max_dim == 2)
	{
		pcs_number_of_primary_nvals = 3;
		dm_number_of_primary_nvals = 3;
		pcs_primary_function_name[2] = "ACCELERATION_Z1";
		pcs_primary_function_unit[2] = "m/s^2";
	}
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE_RATE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa/s";
	pcs_number_of_primary_nvals++;

	//----------------------------------------------------------------------
	// NOD Secondary functions
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	//  pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "POROPRESSURE0";
	//  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	//  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	//  pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;

	// NB new secondary variables for principle stresses; Principle stress directions might follow later
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_1"; // i.e. Sigma_1
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_2"; // i.e. Sigma_2
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_3"; // i.e. Sigma_3
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;

	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_X";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_Y";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// 3D
	if (max_dim == 2) // 3D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DISPLACEMENT_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_DM_Z";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
	if (max_dim > 0) // >1D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_1_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_2_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_X"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_Y"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "NORM_STRESS_3_Z"; // i.e. Sigma_1
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 OK Implementation
   02/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigUnsaturatedFlow()
{
	if ((int)continuum_vector.size() == 1)
	{
		// 1.1 primary variables
		pcs_number_of_primary_nvals = 1;
		pcs_primary_function_name[0] = "PRESSURE1";
		pcs_primary_function_unit[0] = "Pa";
		// 1.2 secondary variables
		// OK LOPCalcSecondaryVariables_USER = MMPCalcSecondaryVariablesRichards; // p_c and S^l
		pcs_number_of_secondary_nvals = 0;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// MB
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// MB
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		if (Neglect_H_ini == 2)
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1_Ini";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++; // WX:08.2011
		}

// TEST
//#define DECOVALEX
#ifdef DECOVALEX
		// DECOVALEX Test
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_I";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
#endif
		if (adaption) // WW, JOD removed	//MW added
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STORAGE_P";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
			pcs_number_of_secondary_nvals++;
		}
		// Nodal velocity. WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		// JOD
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// JOD
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		// for calibration (unconfined groundwater flow with variable density) 5.3.07 JOD
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "HEAD";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
	else if ((int)continuum_vector.size() == 2)
	{
		dof = 2; // WW
		// 1.1 primary variables
		pcs_number_of_primary_nvals = 2; // YD
		pcs_primary_function_name[0] = "PRESSURE1";
		pcs_primary_function_unit[0] = "Pa";
		pcs_primary_function_name[1] = "PRESSURE2";
		pcs_primary_function_unit[1] = "Pa";
		// 1.2 secondary variables
		// OK LOPCalcSecondaryVariables_USER = MMPCalcSecondaryVariablesRichards; // p_c and S^l
		pcs_number_of_secondary_nvals = 0;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++; // WW
		// Nodal velocity. WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WW
		// 03.03.2008. WW
		for (size_t i = 0; i < GetPrimaryVNumber(); i++)
			Shift[i] = i * m_msh->GetNodesNumber(true);
	}

	// Output material parameters
	// WW // TF
	configMaterialParameters();

	// 2 ELE values
	pcs_number_of_evals = 8;
	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
	pcs_eval_name[3] = "POROSITY"; // MX 11.2005
	pcs_eval_unit[3] = "-";
	pcs_eval_name[4] = "POROSITY_IL"; // MX 11.2005
	pcs_eval_unit[4] = "-";
	pcs_eval_name[5] = "PERMEABILITY"; // MX 11.2005
	pcs_eval_unit[5] = "-";
	pcs_eval_name[6] = "n_sw"; // MX 11.2005
	pcs_eval_unit[6] = "-";
	pcs_eval_name[7] = "n_sw_rate"; // MX 11.2005
	pcs_eval_unit[7] = "-";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 PCH Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigFluidMomentum()
{
	// pcs_num_name[0] = "VELOCITY1_X";
	// Nothing added in terms of matrix solver.
	// Just linear solver is good enough.
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";
	// NOD values
	pcs_number_of_primary_nvals = 3;
	pcs_primary_function_name[0] = "VELOCITY1_X";
	pcs_primary_function_unit[0] = "m/s";
	pcs_primary_function_name[1] = "VELOCITY1_Y";
	pcs_primary_function_unit[1] = "m/s";
	pcs_primary_function_name[2] = "VELOCITY1_Z";
	pcs_primary_function_unit[2] = "m/s";

	// I'm adding this to initialize for Fluid Momentum process
	pcs_number_of_secondary_nvals = 0;
	pcs_number_of_evals = 3;

	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
}

/**************************************************************************/
void CRFProcess::ConfigRandomWalk()
{
	// Nothing added in terms of matrix solver.
	// Just linear solver is good enough.
	pcs_sol_name = "LINEAR_SOLVER_PROPERTIES_PRESSURE1";

	// NOD values
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;

	// 2 ELE values
	pcs_number_of_evals = 1;
	pcs_eval_name[0] = "CONCENTRATION0";
	pcs_eval_unit[0] = "kg/m3";

	for (size_t e = 0; e < m_msh->ele_vector.size(); e++)
	{
		m_msh->ele_vector[e]->AllocateMeomoryforAngle();
		m_msh->ele_vector[e]->SetAngle(0, 0.);
		m_msh->ele_vector[e]->SetAngle(1, 0.);
		m_msh->ele_vector[e]->SetAngle(2, 0.);
	}
}

////////////////////////////////////////////////////////////////////////////
//
///  Define variables of multi-phase flow model  (WW 08.2010)
//
////////////////////////////////////////////////////////////////////////////
void CRFProcess::Def_Variable_MultiPhaseFlow()
{
	// 1.1 primary variables
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE2";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	// 1.2 secondary variables
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_W";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	if (Neglect_H_ini == 2)
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1_Ini";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WX:08.2011
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE2_Ini";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WX:08.2011
	}

	// 1.3 elemental variables								// BG, 04/2012
	// pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	// CB_merge_0513 added some additional ELE values
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY";
	pcs_eval_unit[pcs_number_of_evals] = "m^2";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "POROSITY";
	pcs_eval_unit[pcs_number_of_evals] = "-";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY_YY";
	pcs_eval_unit[pcs_number_of_evals] = "m^2";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY_ZZ";
	pcs_eval_unit[pcs_number_of_evals] = "m^2";
	pcs_number_of_evals++;

	if (simulator.compare("ECLIPSE") == 0) // 02.2013. WW
	{
		// BG 01/2011, variables necessary for ECLIPSE
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DENSITY1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VISCOSITY1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DENSITY2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VISCOSITY2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;

		//// WTP: Add Temp support
		// pcs_secondary_function_name[pcs_number_of_secondary_nvals]      = "TEMPERATURE1";
		// pcs_secondary_function_unit[pcs_number_of_secondary_nvals]      = "K";
		// pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		// pcs_number_of_secondary_nvals++;
	}
}

////////////////////////////////////////////////////////////////////////////
//
///  Define variables of liquid flow model  (NW 09.2011)
//
////////////////////////////////////////////////////////////////////////////
void CRFProcess::Def_Variable_LiquidFlow()
{
	// 1.1 primary variables
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "HEAD";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// JOD 2014-11-10
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DELTA_PRESSURE1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "PA";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;

	if (Neglect_H_ini == 2)
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1_Ini";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++; // WX 08.2011
	}

	// 1.3 elemental variables
	// pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VOLUME";
	pcs_eval_unit[pcs_number_of_evals] = "m3";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "POROSITY"; // MX, test for n=n(c), 04.2005
	pcs_eval_unit[pcs_number_of_evals] = "-";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals]
	    = "PERMEABILITY"; // JT 2010 -- need this for index call of heterogeneous permeability
	pcs_eval_unit[pcs_number_of_evals] = "m2";
	pcs_number_of_evals++;

	// WTP: needed?
	// CB_merge_0513 ?? why only for eclipse??
	// if(simulator.compare("ECLIPSE") == 0) //02.2013. WW
	//    {
	// BG 01/2011, variables necessary for ECLIPSE
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DENSITY1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VISCOSITY1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	//    }
}

/**************************************************************************
   FEMLib-Method: For non-isothermal multi-phase flow
   Task:
   Programing:
   02/2007 WW Implementation
   04/2011 WW Apdate for H2M
**************************************************************************/
void CRFProcess::ConfigMultiPhaseFlow()
{
	dof = 2;
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;

	Def_Variable_MultiPhaseFlow();

	// Output material parameters
	configMaterialParameters();

	// 11.08.2010. WW
	long nn = m_msh->GetNodesNumber(false);
	//
	for (size_t i = 0; i < GetPrimaryVNumber(); i++) // 03.03.2008. WW
		Shift[i] = i * nn;

	num_nodes_p_var = new long[2];
	num_nodes_p_var[0] = num_nodes_p_var[1] = nn;
}

/**************************************************************************
   FEMLib-Method: For PS model for multiphase flow
   Task:
   Programing:
   03/2009 PCH Implementation
**************************************************************************/
void CRFProcess::ConfigPS_Global()
{
	dof = 2;
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 2;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	pcs_primary_function_name[1] = "SATURATION2";
	pcs_primary_function_unit[1] = "m3/m3";
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// CB_merge_0513 additional secondary variable
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "DENSITY2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	// CB_merge_0513 added some ELE values
	// pcs_number_of_evals = 2;
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY";
	pcs_eval_unit[pcs_number_of_evals] = "m^2";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "POROSITY";
	pcs_eval_unit[pcs_number_of_evals] = "-";
	pcs_number_of_evals++;

	//
	for (size_t i = 0; i < GetPrimaryVNumber(); i++) // 03.03.2008. WW
		Shift[i] = i * m_msh->GetNodesNumber(true);
}

/**************************************************************************
   FEMLib-Method: For Pressure-temperature-coupled flow for fluids
   Task:
   Programing:
   02/2011 AKS/NB Implementation
**************************************************************************/
void CRFProcess::ConfigMULTI_COMPONENTIAL_FLOW()
{
	const char* CmpName[6];
	CFluidProperties* m_mfp = NULL;
	m_mfp = mfp_vector[0];
	dof = 2 + m_mfp->cmpN;
	CmpName[2] = (char*)m_mfp->cmpNm1.data();
	CmpName[3] = (char*)m_mfp->cmpNm2.data();
	CmpName[4] = (char*)m_mfp->cmpNm3.data();
	CmpName[5] = (char*)m_mfp->cmpNm4.data();
	// 1.1 primary variables
	pcs_number_of_primary_nvals = dof;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_name[1] = "TEMPERATURE1";
	for (int CIndex = 2; CIndex < m_mfp->cmpN + 2; CIndex++)
		pcs_primary_function_name[CIndex] = CmpName[CIndex];

	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 0;
	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		Shift[i] = i * m_msh->GetNodesNumber(true);
}
/**************************************************************************
FEMLib-Method: For reactive thermal non-equilibrium flow
Task:
Programing:
07/2013 HS/TN Implementation
**************************************************************************/
void CRFProcess::ConfigTNEQ()
{
	dof = 4;
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 4;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	pcs_primary_function_name[1] = "TEMPERATURE1";
	pcs_primary_function_unit[1] = "K";
	pcs_primary_function_name[2] = "TEMPERATURE2";
	pcs_primary_function_unit[2] = "K";
	pcs_primary_function_name[3] = "CONCENTRATION1";
	pcs_primary_function_unit[3] = "kg/kg";
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 0;
	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// this is the nodal based reaction rate of solid density in PTC flow.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "REACT_RATE_N"; // TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// this is the nodal based solid density in PTC flow.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SOLID_DENSITY_N"; // TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// Solid heat capacity output
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SOLID_HEAT_CAP_N"; //TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "J/kgK";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	// pcs_eval_name[pcs_number_of_evals] = "POROSITY";
	// pcs_eval_unit[pcs_number_of_evals] = "-";
	// pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY";
	pcs_eval_unit[pcs_number_of_evals] = "m2";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "SOLID_DENSITY";
	pcs_eval_unit[pcs_number_of_evals] = "kg/m3";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "REACT_RATE";
	pcs_eval_unit[pcs_number_of_evals] = "kg/m3s";
	pcs_number_of_evals++;
	// pcs_eval_name[pcs_number_of_evals] = "SOLID_HEAT_CAP";
	// pcs_eval_unit[pcs_number_of_evals] = "J/kgK";
	// pcs_number_of_evals++;
	//pcs_eval_name[pcs_number_of_evals] = "T_EQUIL"; //TN
	//pcs_eval_unit[pcs_number_of_evals] = "K";
	//pcs_number_of_evals++;

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		Shift[i] = i * m_msh->GetNodesNumber(false);

	long group_count = 0;
	m_rho_s_0 = msp_vector[group_count]->Density(); // get initial solid density
	double poro = mmp_vector[group_count]->porosity;
	FiniteElement::SolidReactiveSystem react_syst = msp_vector[group_count]->getSolidReactiveSystem();
	m_conversion_rate = new conversion_rate(573.0, // T_solid, Kelvin
	                                        573.0, // T_gas, Kelvin
	                                        0.0, // p_gas
	                                        0.0, // w_water, mass fraction unitless
	                                        m_rho_s_0, // kg/m3, // rho_s_initial,
	                                        1.0 - poro, // solid volume fraction
	                                        1.0, // delta_t
	                                        react_syst);
}

/**************************************************************************
FEMLib-Method: For reactive thermal non-equilibrium flow
Task:
Programing:
07/2013 HS/TN Implementation
**************************************************************************/
void CRFProcess::ConfigTES()
{
	dof = 3;
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 3;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	pcs_primary_function_name[1] = "TEMPERATURE1";
	pcs_primary_function_unit[1] = "K";
	pcs_primary_function_name[2] = "CONCENTRATION1";
	pcs_primary_function_unit[2] = "kg/kg";
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 0;
	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// this is the nodal based reaction rate of solid density in PTC flow.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "REACT_RATE_N"; // TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// this is the nodal based solid density in PTC flow.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SOLID_DENSITY_N"; // TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "kg/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// Solid heat capacity output
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SOLID_HEAT_CAP_N"; //TN TEST
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "J/kgK";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	// pcs_eval_name[pcs_number_of_evals] = "POROSITY";
	// pcs_eval_unit[pcs_number_of_evals] = "-";
	// pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "PERMEABILITY";
	pcs_eval_unit[pcs_number_of_evals] = "m2";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "SOLID_DENSITY";
	pcs_eval_unit[pcs_number_of_evals] = "kg/m3";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "REACT_RATE";
	pcs_eval_unit[pcs_number_of_evals] = "kg/m3s";
	pcs_number_of_evals++;
	// pcs_eval_name[pcs_number_of_evals] = "SOLID_HEAT_CAP";
	// pcs_eval_unit[pcs_number_of_evals] = "J/kgK";
	// pcs_number_of_evals++;
	// pcs_eval_name[pcs_number_of_evals] = "T_EQUIL"; //TN
	// pcs_eval_unit[pcs_number_of_evals] = "K";
	// pcs_number_of_evals++;

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		Shift[i] = i * m_msh->GetNodesNumber(false);

	long group_count = 0;
	m_rho_s_0 = msp_vector[group_count]->Density(); // get initial solid density
	double poro = mmp_vector[group_count]->porosity;
	FiniteElement::SolidReactiveSystem react_syst = msp_vector[group_count]->getSolidReactiveSystem();
	m_conversion_rate = new conversion_rate(573.0, // T_solid, Kelvin
	                                        573.0, // T_gas, Kelvin
	                                        0.0, // p_gas
	                                        0.0, // w_water, mass fraction unitless
	                                        m_rho_s_0, // kg/m3, // rho_s_initial,
	                                        1.0 - poro, // solid volume fraction
	                                        1.0, // delta_t
	                                        react_syst);
}

//////////////////////////////////////////////////////////////////////////
// Configuration NOD
//////////////////////////////////////////////////////////////////////////
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 07.11.2008
/*************************************************************************
   ROCKFLOW - Function:
   Task: Config node values
   Programming: 02/2003 OK Implementation
   04/2004   WW   Modification for 3D problems
   last modified:
 **************************************************************************/
void CRFProcess::ConfigNODValues1(void)
{
	int i;
	int pcs_nval = 0;
	const int DOF = GetPrimaryVNumber();
	anz_nval0 = anz_nval;
	number_of_nvals = 2 * DOF + pcs_number_of_secondary_nvals;
	// NVAL
	pcs_nval_data = (PCS_NVAL_DATA*)Malloc(number_of_nvals * sizeof(PCS_NVAL_DATA));
	/*----------------------------------------------------------------*/
	for (i = 0; i < DOF; i++)
	{
		/* Primary variable - old time */
		// NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
		strcpy(pcs_nval_data[pcs_nval].name, pcs_primary_function_name[pcs_nval - i]);
		// Change name for the previous time level
		/*
		   char *ch = strchr(pcs_nval_data[pcs_nval].name, '1');
		   if( ch != NULL )
		   {
		   int pos = ch-pcs_nval_data[pcs_nval].name;
		   pcs_nval_data[pcs_nval].name[pos]='0';
		   }
		   else
		   strcat(pcs_nval_data[pcs_nval].name, "0");   */
		//-------------------------------------------------------------------------------
		strcpy(pcs_nval_data[pcs_nval].einheit, pcs_primary_function_unit[pcs_nval - i]);
		pcs_nval_data[pcs_nval].timelevel = 0;
		pcs_nval_data[pcs_nval].speichern = 0; // WW
		pcs_nval_data[pcs_nval].laden = 0;
		pcs_nval_data[pcs_nval].restart = 1;
		pcs_nval_data[pcs_nval].adapt_interpol = 1;
		pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
		pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
		pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
		pcs_nval++;
		/* Primary variable - new time */
		// NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
		strcpy(pcs_nval_data[pcs_nval].name, pcs_primary_function_name[pcs_nval - i - 1]);
		strcpy(pcs_nval_data[pcs_nval].einheit, pcs_primary_function_unit[pcs_nval - i - 1]);
		pcs_nval_data[pcs_nval].timelevel = 1;
		pcs_nval_data[pcs_nval].speichern = 1;

		pcs_nval_data[pcs_nval].laden = 0;
		pcs_nval_data[pcs_nval].restart = 1;
		pcs_nval_data[pcs_nval].adapt_interpol = 1;
		pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
		pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
		pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
		pcs_nval++;
	}

	/*----------------------------------------------------------------*/
	/* Secondary variables */
	for (i = 0; i < pcs_number_of_secondary_nvals; i++)
	{
		// NVAL pcs_nval_data[pcs_nval] = (PCS_NVAL_DATA *) Malloc(sizeof(PCS_NVAL_DATA));
		strcpy(pcs_nval_data[pcs_nval].name, pcs_secondary_function_name[i]);
		strcpy(pcs_nval_data[pcs_nval].einheit, pcs_secondary_function_unit[i]);
		//  pcs_nval_data[i+2]->timelevel = 1; // always at new time level
		pcs_nval_data[pcs_nval].timelevel = pcs_secondary_function_timelevel[i];
		if (pcs_nval_data[pcs_nval].timelevel == 1)
			pcs_nval_data[pcs_nval].speichern = 1;
		else
			pcs_nval_data[pcs_nval].speichern = 0;
		pcs_nval_data[pcs_nval].laden = 0;
		pcs_nval_data[pcs_nval].restart = 1;
		pcs_nval_data[pcs_nval].adapt_interpol = 1;
		pcs_nval_data[pcs_nval].vorgabe = 0.0;
#ifdef PCS_NOD
		pcs_nval_data[pcs_nval].nval_index = pcs_nval;
#else
		pcs_nval_data[pcs_nval].nval_index = anz_nval + pcs_nval;
#endif
		if (type == 4 || type == 41)
			if (dm_number_of_primary_nvals == 2 || (dm_number_of_primary_nvals == 3 && this->type == 41))
			{
				// Block:
				// STRESS_ZX, STRESS_YZ, STRAIN_ZX, STRAIN_ZY and LUMPED_STRESS
				if (i == 3 || i == 4 || i == 9 || i == 10 || i == 13)
					pcs_nval_data[pcs_nval].speichern = 0;
				if (!problem_2d_plane_dm)
					//  Block STRESS_ZZ and STRAIN_ZZ
					if (i == 5 || i == 11)
						pcs_nval_data[pcs_nval].speichern = 0;
				// if(!pcs_plasticity)
				// {
				//   if(i==12) pcs_nval_data[pcs_nval].speichern = 0;  //STRAIN_PLS
				// }
			}
		pcs_nval++;
	}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::PCSConfigNODValues
   Task: Config node values
   Programming: 02/2003 OK Implementation
   last modified:
 **************************************************************************/
void CRFProcess::ConfigNODValues2(void)
{
	int i;

	number_of_nvals = 2 * GetPrimaryVNumber() + pcs_number_of_secondary_nvals;
	for (i = 0; i < number_of_nvals; i++)
		ModelsAddNodeValInfoStructure(pcs_nval_data[i].name, pcs_nval_data[i].einheit, pcs_nval_data[i].speichern,
		                              pcs_nval_data[i].laden, pcs_nval_data[i].restart, pcs_nval_data[i].adapt_interpol,
		                              pcs_nval_data[i].vorgabe);
}
#endif //#ifndef NEW_EQS //WW. 07.11.2008
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   07/2004 OK Implementation
**************************************************************************/
int PCSGetNODValueIndex(const string& name, int timelevel)
{
	// PCS primary variables
	int pcs_vector_size = (int)pcs_vector.size();
	int i, j;
	CRFProcess* m_pcs = NULL;
	if (pcs_vector_size > 0)
		for (i = 0; i < pcs_vector_size; i++)
		{
			m_pcs = pcs_vector[i];
			for (j = 0; j < m_pcs->number_of_nvals; j++)
				if ((name.compare(m_pcs->pcs_nval_data[j].name) == 0)
				    && (m_pcs->pcs_nval_data[j].timelevel == timelevel))
					return m_pcs->pcs_nval_data[j].nval_index;
		}
	cout << "Error in PCSGetNODValueIndex: " << name << "\n";
	return -1;
}

//////////////////////////////////////////////////////////////////////////
// Configuration ELE
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Function:
   Task: Config element values
   Programming: 02/2003 OK Implementation
   last modified:
   06/2004  WW
 **************************************************************************/
void CRFProcess::ConfigELEValues1(void)
{
	int i;
	if (pcs_number_of_evals)
		pcs_eval_data = (PCS_EVAL_DATA*)Malloc(pcs_number_of_evals * sizeof(PCS_EVAL_DATA));
	for (i = 0; i < pcs_number_of_evals; i++)
	{
		// pcs_eval_data[i] = (PCS_EVAL_DATA *) Malloc(sizeof(PCS_EVAL_DATA));
		strcpy(pcs_eval_data[i].name, pcs_eval_name[i]);
		strcpy(pcs_eval_data[i].einheit, pcs_eval_unit[i]);
		pcs_eval_data[i].speichern = 1;
		pcs_eval_data[i].laden = 0;
		pcs_eval_data[i].restart = 1;
		pcs_eval_data[i].adapt_interpol = 1;
		pcs_eval_data[i].vorgabe = 0.0;
		pcs_eval_data[i].index = anz_eval + i;
		pcs_eval_data[i].eval_index = anz_eval + i; // SB
	}
}

/*************************************************************************
   ROCKFLOW - Function:
   Task: Config element values
   Programming: 04/2003 OK Implementation
   last modified:
 **************************************************************************/
void CRFProcess::ConfigELEValues2(void)
{
	int i;
	for (i = 0; i < pcs_number_of_evals; i++)
		ModelsAddElementValInfoStructure(pcs_eval_data[i].name, pcs_eval_data[i].einheit, pcs_eval_data[i].speichern,
		                                 pcs_eval_data[i].laden, pcs_eval_data[i].restart,
		                                 pcs_eval_data[i].adapt_interpol, pcs_eval_data[i].vorgabe);
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::PCSGetELEValueIndex
   Task: Provide index for element values
   Programming: 08/2003 SB Implementation
   last modified:
 **************************************************************************/
int PCSGetELEValueIndex(char* name)
{
	int i;
	CRFProcess* m_process = NULL;
	int j;
	int no_processes = (int)pcs_vector.size();
	for (j = 0; j < no_processes; j++)
	{
		m_process = pcs_vector[j];
		for (i = 0; i < m_process->pcs_number_of_evals; i++)
			if (strcmp(m_process->pcs_eval_data[i].name, name) == 0)
				return m_process->pcs_eval_data[i].eval_index;
	}
	printf("PCSGetELEValueIndex Alert\n");
	printf("%s \n", name);
	return -1;
}

//////////////////////////////////////////////////////////////////////////
// Configuration ELE matrices
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:  Activate or deactivate elements specified in .pcs file
   Programing:
   05/2005 WW Implementation
**************************************************************************/
void CRFProcess::CheckMarkedElement()
{
	size_t i, j;
	bool done;
	CElem* elem = NULL;
	CNode* node = NULL;

	size_t ele_vector_size(m_msh->ele_vector.size());

	for (size_t l = 0; l < ele_vector_size; l++)
	{
		elem = m_msh->ele_vector[l];
		done = false;
		for (i = 0; i < (size_t)NumDeactivated_SubDomains; i++)
		{
			if (elem->GetPatchIndex() == static_cast<size_t>(Deactivated_SubDomain[i]))
			{
				elem->MarkingAll(false);
				done = true;
				break;
			}
		}
		// WX:02.2013: excav with deactivated subdomain
		if (ExcavMaterialGroup == (int)elem->GetPatchIndex())
		{
			if (!elem->GetMark() || abs(elem->GetExcavState()) < MKleinsteZahl)
			{
				elem->MarkingAll(false);
				done = true;
			}
		}
		if (done)
			continue;
		else
			elem->MarkingAll(true);
	}
	size_t node_vector_size = m_msh->nod_vector.size();
	for (size_t l = 0; l < node_vector_size; l++)
		while (m_msh->nod_vector[l]->getConnectedElementIDs().size())
			m_msh->nod_vector[l]->getConnectedElementIDs().pop_back();

	for (size_t l = 0; l < ele_vector_size; l++)
	{
		elem = m_msh->ele_vector[l];
		if (!elem->GetMark())
			continue;
		for (i = 0; i < elem->GetNodesNumber(m_msh->getOrder()); i++)
		{
			done = false;
			node = elem->GetNode(i);
			for (j = 0; j < node->getConnectedElementIDs().size(); j++)
				if (l == node->getConnectedElementIDs()[j])
				{
					done = true;
					break;
				}
			if (!done)
				node->getConnectedElementIDs().push_back(l);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:  check the excavation state of each aktive element
   Programing:
   01/2011 WX Implementation
**************************************************************************/
void CRFProcess::CheckExcavedElement()
{
	int valid;
	long l;
	// bool done;
	CElem* elem = NULL;
	// CNode *node = NULL;
	for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
	{
		elem = m_msh->ele_vector[l];
		if ((int)elem->GetPatchIndex() == ExcavMaterialGroup && elem->GetExcavState() == -1) // WX:04.2012
		{
			double const* ele_center(elem->GetGravityCenter());
			if ((GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) + ExcavBeginCoordinate)
			        > (ele_center[ExcavDirection])
			    && (ele_center[ExcavDirection] - ExcavBeginCoordinate) > -0.001)
				elem->SetExcavState(1);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// PCS Execution
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   02/2003 OK Implementation
   04/2003 OK Storing time step results
   08/2004 OK PCS2
   01/2005 WW new ELE concept
   03/2005 OK MultiMSH concept
   06/2005 MB NEWTON error calculation
   05/2007 WW DOF>1, unerrelaxation for Picard, and removing the old or spurious stuff
   12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear solver
   and parellelisation of them
   07/2011 WW Add Newton-Rahson and reduce #ifdef
   3/2012  JT Clean, and correct the error obtainment for NLS and CPL
   last modified:
 **************************************************************************/
double CRFProcess::Execute()
{
	int nidx1;
	double pcs_error, nl_theta;
	long j, k, g_nnodes; // 07.01.07 WW

	pcs_error = DBL_MAX;
	g_nnodes = m_msh->GetNodesNumber(false);

	// if const, set node values directly
	if (_pcs_constant_model > 0)
	{
		double const_val = _pcs_constant_value;
		if (_pcs_constant_model == 2)
		{
			int valid;
			const_val = GetCurveValue(_pcs_constant_curve, 0, aktuelle_zeit, &valid);
		}
		std::cout << "-> set given constant value. no need to solve a linear system.\n";
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
				k = m_msh->Eqs2Global_NodeIndex[j];
				SetNodeValue(k, nidx1, const_val);
			}
		}
		return .0;
	}

#if !defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	double implicit_lim = 1.0 - DBL_EPSILON;
#ifdef NEW_EQS
	eqs_x = eqs_new->x;
#else
	eqs_x = eqs->x;
#endif
#endif

#ifdef USE_MPI // WW
	long global_eqs_dim = pcs_number_of_primary_nvals * m_msh->GetNodesNumber(false);
	CPARDomain* dom = dom_vector[myrank];
#endif

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	eqs_new->Initialize();
#elif defined(NEW_EQS) // WW
	if (!configured_in_nonlinearloop)
#if defined(USE_MPI)
	{
		// 02.2010.  WW
		dom->eqs->SetDOF(pcs_number_of_primary_nvals);
		dom->ConfigEQS(m_num, global_eqs_dim);
	}
// TEST eqs_new->Initialize();
#else
	// Also allocate temporary memory for linear solver. WW
	{
		//_new 02/2010. WW
		eqs_new->SetDOF(pcs_number_of_primary_nvals);
		eqs_new->ConfigNumerics(m_num);
	}
	eqs_new->Initialize();
#endif
#else
	SetLinearSolverType(eqs, m_num); // WW
	SetZeroLinearSolver(eqs);
#endif
	/*
	   //TEST_MPI
	   string test = "rank";
	   char stro[1028];
	   sprintf(stro, "%d",myrank);
	   string test1 = test+(string)stro+"Assemble.txt";
	   ofstream Dum(test1.c_str(), ios::out); // WW
	   dom->eqs->Write(Dum);   Dum.close();
	   MPI_Finalize();
	   exit(1);
	 */
	//----------------------------------------------------------------------
	nl_theta = 1.0 - m_num->nls_relaxation;
	if (nl_theta < DBL_EPSILON)
		nl_theta = 1.0;

	// NW. should mark active elements if any process uses deactivation
	// if(NumDeactivated_SubDomains>0)
	// TODO if it's nonlinear, CheckMarkedElement() has been already called
	if (hasAnyProcessDeactivatedSubdomains)
#ifdef NEW_EQS // WW
		if (!configured_in_nonlinearloop)
#endif
			CheckMarkedElement();
	m_msh->SwitchOnQuadraticNodes(false);

	// If not Newton-Raphson method. 20.07.2011. WW
	if (m_num->nls_method < 1)
	{
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
		InitializeRHS_with_u0();
#else
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			long const ish = ii * g_nnodes;
			for (j = 0; j < g_nnodes; j++) // WW
				eqs_x[j + ish] = GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
		}
#endif
	}

//---------------------------------------------------------------------
// Assembly
#if defined(USE_MPI) || defined(USE_PETSC) // WW
	clock_t cpu_time = 0; // WW
	cpu_time = -clock();
	if (myrank == 0)
#endif
		cout << "      Assembling equation system..."
		     << "\n";
	GlobalAssembly();
#if defined(USE_MPI) || defined(USE_PETSC) // WW
	cpu_time += clock();
	cpu_time_assembly += cpu_time;
	if (myrank == 0)
#endif
		cout << "      Calling linear solver..."
		     << "\n";
#ifdef CHECK_EQS
	std::string eqs_name = convertProcessTypeToString(this->getProcessType()) + "_EQS.txt";
	MXDumpGLS((char*)eqs_name.c_str(), 1, eqs->b, eqs->x);
#endif
//----------------------------------------------------------------------
// Execute linear solver
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	eqs_new->Solver();
	// TEST 	double x_norm = eqs_new->GetVecNormX();
	eqs_new->MappingSolution();
#elif defined(NEW_EQS) // WW
#if defined(USE_MPI)
	// 21.12.2007
	iter_lin = dom->eqs->Solver(eqs_new->x, global_eqs_dim);
#else
#if defined(LIS) || defined(MKL)
	iter_lin = eqs_new->Solver(this->m_num); // NW
#else
	iter_lin = eqs_new->Solver();
#endif
#endif
#else
	iter_lin = ExecuteLinearSolver();
#endif
	iter_lin_max = std::max(iter_lin_max, iter_lin);

	//----------------------------------------------------------------------
	// Linearized Flux corrected transport (FCT) by Kuzmin 2009
	//----------------------------------------------------------------------
	if (m_num->fct_method > 0) // NW
	{
#if defined(USE_PETSC)
		eqs_x = eqs_new->GetGlobalSolution();
		pcs_error = CalcIterationNODError(1);
#else
		pcs_error = CalcIterationNODError(m_num->getNonLinearErrorMethod(), true, false); // JT
#endif

#if defined(USE_MPI) || defined(USE_PETSC)
		if (myrank == 0)
		{
#endif
			cout << "    Relative PCS error: " << pcs_error << "\n";
			cout << "    Start FCT calculation"
			     << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
		}
#endif
		// Set u^H: use the solution as the higher-order solution
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC)
				k = m_msh->Eqs2Global_NodeIndex[j] * pcs_number_of_primary_nvals + ii;
				SetNodeValue(j, nidx1, eqs_x[k]);
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				SetNodeValue(k, nidx1, eqs_x[j + ii * g_nnodes]);
#endif
			}
		}

// Initialize the algebra system
#if defined(USE_PETSC)
		eqs_new->Initialize();
#else
#ifdef NEW_EQS // WW
		if (!configured_in_nonlinearloop)
#if defined(USE_MPI)
			dom->ConfigEQS(m_num, global_eqs_dim);
#else
			eqs_new->Initialize();
#endif
#else
		SetZeroLinearSolver(eqs);
#endif
#endif

// Set initial guess
#if !defined(USE_PETSC)
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
				k = m_msh->Eqs2Global_NodeIndex[j];
				eqs_x[j + ii * g_nnodes] = GetNodeValue(k, nidx1);
			}
		}
#endif

// Assembly
#if defined(USE_MPI) || defined(USE_PETSC)
		clock_t cpu_time = 0; // WW
		cpu_time = -clock();
#endif
		femFCTmode = true;
		GlobalAssembly();
		femFCTmode = false;
#if defined(USE_MPI) || defined(USE_PETSC)
		cpu_time += clock();
		cpu_time_assembly += cpu_time;
#endif

// Solve the algebra
#ifdef CHECK_EQS
		string eqs_name = pcs_type_name + "_EQS.txt";
		MXDumpGLS((char*)eqs_name.c_str(), 1, eqs->b, eqs->x);
#endif

#if defined(USE_PETSC)
		//		std::string eqs_output_file = FileName + number2str(aktueller_zeitschritt);
		//		eqs_new->EQSV_Viewer(eqs_output_file);
		eqs_new->Solver();
		eqs_new->MappingSolution();
#else
#ifdef NEW_EQS // WW
#if defined(USE_MPI)
		// 21.12.2007
		dom->eqs->Solver(eqs_new->x, global_eqs_dim);
#else
#if defined(LIS) || defined(MKL)
		eqs_new->Solver(this->m_num); // NW
#else
		eqs_new->Solver();
// kg44 the next lines are for debug?
//		string fname = FileName + "_equation_results.txt";
//		ofstream dum(fname.c_str(), ios::out | ios::trunc);
//		eqs_new->Write(dum);
//		exit(1);
#endif
#endif
#else // ifdef NEW_EQS
		ExecuteLinearSolver();
#endif
#endif // USE_PETSC
	}
	//----------------------------------------------------------------------
	// END OF FLUX CORRECTED TRANSPORT
	//----------------------------------------------------------------------

	// PCSDumpModelNodeValues();
	//----------------------------------------------------------------------
	// ERROR CALCULATION
	//----------------------------------------------------------------------

	//
	// Save the solution of the prevoius iteration of the nonlinear step for
	// the automatic time stepping. Modified for PETsc solver. 03.07.2012. WW
	if (Tim->GetPITimeStepCrtlType() > 0)
	{
		double* x_k = NULL;
		bool get_buffer_u_k = true;
		x_k = _problem->GetBufferArray(get_buffer_u_k);
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
#if !defined(USE_PETSC) // && !defined(other parallel libs)
			const long ish = i * g_nnodes;
#endif
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC) // || defined(other parallel libs)
				x_k[j * pcs_number_of_primary_nvals + i] = GetNodeValue(j, nidx1);
#else
				x_k[j + ish] = GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
#endif
			}
		}
	}
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	// PICARD
	//----------------------------------------------------------------------
	// Error calculation
	//----------------------------------------------------------------------
	if (m_num->nls_method_name.find("PICARD") != string::npos)
	{
		eqs_x = eqs_new->GetGlobalSolution();
		//......................................................................
		pcs_error = CalcIterationNODError(1); // OK4105//WW4117
		if (myrank == 0)
			cout << "      PCS error: " << pcs_error << "\n";

		//--------------------------------------------------------------------
		// 7 Store solution vector in model node values table
		//....................................................................
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
				k = m_msh->Eqs2Global_NodeIndex[j] * pcs_number_of_primary_nvals + i;
				SetNodeValue(j, nidx1, (1. - nl_theta) * GetNodeValue(j, nidx1) + nl_theta * eqs_x[k]);
			}
		}

	} // END PICARD
#else
	// JT: Coupling error was wrong. Now ok.
	if (iter_nlin > 0)
	{ // Just getting NL error
		pcs_error = CalcIterationNODError(m_num->getNonLinearErrorMethod(), true, false); // OK4105//WW4117//JT
	}
	else
	{ // Getting NL and CPL error
		pcs_error = CalcIterationNODError(m_num->getCouplingErrorMethod(), true, true); // JT2012
		if (m_num->getNonLinearErrorMethod()
		    != m_num->getCouplingErrorMethod()) // JT: If CPL error method is different, must call separately
			pcs_error = CalcIterationNODError(m_num->getNonLinearErrorMethod(), true,
			                                  false); // JT2012 // get the NLS error. CPL was obtained before.
	}

	//----------------------------------------------------------------------
	// PICARD
	//----------------------------------------------------------------------
	if (m_num->nls_method_name.find("PICARD") != string::npos)
	{
		if (pcs_error < 1.0) // JT: Then the solution has converged, take the final value
			nl_theta = 1.0;
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
		eqs_x = eqs_new->GetGlobalSolution();
#endif
		//
		if (nl_theta > implicit_lim) // This is most common. So go for the lesser calculations.
		{
			for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
				const long nshift = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
					k = pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[j] + ii;
					const double val_n = GetNodeValue(j, nidx1); // 03.04.2009. WW
					SetNodeValue(j, nidx1, eqs_x[k]);
					eqs_x[k] = val_n; // Used for time stepping. 03.04.2009. WW
#else
					k = m_msh->Eqs2Global_NodeIndex[j];
					const double val_n = GetNodeValue(k, nidx1); // 03.04.2009. WW
					SetNodeValue(k, nidx1, eqs_x[j + nshift]);
					eqs_x[j + nshift] = val_n; // Used for time stepping. 03.04.2009. WW
#endif
				}
			}
		}
		else
		{
			for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
				const long nshift = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
					k = pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[j] + ii;
					const double val_n = GetNodeValue(j, nidx1);
					SetNodeValue(j, nidx1, (1.0 - nl_theta) * val_n + nl_theta * eqs_x[k]);
					eqs_x[k] = val_n; // Used for time stepping. 03.04.2009. WW
#else
					k = m_msh->Eqs2Global_NodeIndex[j];
					const double val_n = GetNodeValue(k, nidx1); // 03.04.2009. WW
					SetNodeValue(k, nidx1, (1.0 - nl_theta) * val_n + nl_theta * eqs_x[j + nshift]);
					eqs_x[j + nshift] = val_n; // Used for time stepping. 03.04.2009. WW
#endif
				}
			}
		}

		// maybe works also for other processes involving velocities
		// update nod velocity if constrained BC
		// do I really need to check every single bc node, or how can I access a bc group?
		if ((this->hasConstrainedBC() || this->hasConstrainedST()) && this->accepted
		    && (this->getProcessType() == FiniteElement::RICHARDS_FLOW
		        || this->getProcessType() == FiniteElement::LIQUID_FLOW))
		{
			this->CalIntegrationPointValue();
			this->Extropolation_GaussValue();
		}
	}
//----------------------------------------------------------------------
// END OF PICARD
//----------------------------------------------------------------------
#endif

#ifdef NEW_EQS // WW
	if (!configured_in_nonlinearloop)
#if defined(USE_MPI)
		dom->eqs->Clean();
#else
		// Also allocate temporary memory for linear solver. WW
		eqs_new->Clean();
#endif
#endif

	if (BASELIB::isNAN(pcs_error))
	{
		accepted = false;
	}

	return pcs_error;
}

/*************************************************************************
   GEOSYS - Function:
   Task:
   Programming:
   08/2008 WW Implementation
   11/2008 WW Update
   last modified:
 **************************************************************************/
void CRFProcess::CopyU_n()
{
	int i, nidx1;
	long g_nnodes, j;

	double* temp_v = _problem->GetBufferArray();

	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		/// H2M with monilithic scheme. 02.2011. WW
		if (type == 42)
		{
			nidx1 = p_var_index[i] + 1;
			g_nnodes = num_nodes_p_var[i];
		}
		else
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
			g_nnodes = m_msh->GetNodesNumber(false); // DOF>1, WW
		}
#if !defined(USE_PETSC) // && !defined(other parallel libs)
		const long ish = i * g_nnodes;
#endif
		for (j = 0; j < g_nnodes; j++)
		{
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
			temp_v[j * pcs_number_of_primary_nvals + i] = GetNodeValue(j, nidx1);
#else
			temp_v[j + ish] = GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
#endif
		}
	}
}

/*************************************************************************
   ROCKFLOW - Function:
   Task: Calculate element matrices
   Programming:
   05/2003 OK Implementation
   09/2005 OK gas flow removed
   last modified:
 **************************************************************************/
void CRFProcess::CalculateElementMatrices(void)
{
	switch (this->type)
	{
		case 1: // SM
		case 2: // MTM2
		case 3: // HTM
		case 5: // Gas flow
		case 11:
		case 12: // MMP
		case 13: // MPC
		case 14: // TNEQ
		case 15: // TES
		case 66: // OF
			break;
		default:
			DisplayMsgLn("CalculateElementMatrices: no CalculateElementMatrices specified");
			abort();
	}
}

/*************************************************************************
   GeoSys-Function:
   Task: Algebraic operation for the flux corrected transport (FCT)
   Programming:
   04/2010 NW Implementation
   last modified:
   05/2013 NW Support PETSc parallelization
 **************************************************************************/
void CRFProcess::AddFCT_CorrectionVector()
{
	int idx0 = 0;
	int idx1 = idx0 + 1;
	const double theta = this->m_num->ls_theta;
	const size_t node_size = m_msh->GetNodesNumber(false);
	SparseMatrixDOK::mat_t& fct_f = this->FCT_AFlux->GetRawData();
	SparseMatrixDOK::col_t* col;
	SparseMatrixDOK::mat_t::const_iterator ii;
	SparseMatrixDOK::col_t::const_iterator jj;
	Math_Group::Vec* ML = this->Gl_ML;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL; // WW
	// if(m_dom)
	//  A = m_dom->eqs->A;
	// else
	A = this->eqs_new->A;
#endif

#ifdef USE_PETSC
	// gather K
	FCT_MPI::gatherK(FCT_MPI::ct, *FCT_K);
	// compute D
	FCT_MPI::computeD(m_msh, *FCT_K, *FCT_d);
#endif

	// List of Dirichlet nodes
	std::set<long> list_bc_nodes;
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* bc_node = bc_node_value[i];
		long nod_id = bc_node->geo_node_number;
		list_bc_nodes.insert(nod_id);
	}

	//----------------------------------------------------------------------
	// Construct global matrices: antidiffusive flux(f_ij), positivity matrix(L)
	// - f_ij = 1/dt*m_ij*(DeltaU_ij^H-DeltaU_ij^n)-theta*d_ij^H*DeltaU_ij^H-(1-theta)*d_ij^n*DeltaU_ij^n
	// - L = K + D
	// - D_ij = min(0, -K_ij, -K_ji)
	// * f_ji = -f_ij
	// * K:original coefficient matrix, D:artificial diffusion operator
	// Implementation memo:
	// - K is stored in A matrix in the element assembly.
	// - the first part of the antidiffusive flux is done in the element assembly.
	//   -> f_ij = m_ij
	//----------------------------------------------------------------------
	// f_ij*=1/dt*(DeltaU_ij^H-DeltaU_ij^n)  for i!=j
	for (size_t i = 0; i < node_size; i++)
	{
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			if (i > j)
				continue; // symmetric part, off-diagonal
			double diff_uH = this->GetNodeValue(i, idx1) - this->GetNodeValue(j, idx1);
			double diff_u0 = this->GetNodeValue(i, idx0) - this->GetNodeValue(j, idx0);
			double v = 1.0 / dt * (diff_uH - diff_u0);
			(*FCT_AFlux)(i, j) *= v; // MC is already done in local ele assembly
			(*FCT_AFlux)(j, i) *= -v; // MC is already done in local ele assembly
		}
	}

	// Complete f, L
	// Remark: Using iteration is only possible after the sparse table has been constructed.
	for (size_t i = 0; i < node_size; i++)
	{
#ifdef USE_PETSC
		const size_t i_global = FCT_GLOB_ADDRESS(i);
#endif
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			if (i > j || i == j)
				continue; // do below only for upper triangle due to symmetric

// Get artificial diffusion operator D
#ifdef USE_PETSC
			const size_t j_global = FCT_GLOB_ADDRESS(j);
			double d1 = (*FCT_d)(i_global, j_global);
#else
#if defined(NEW_EQS)
			double K_ij = (*A)(i, j);
			double K_ji = (*A)(j, i);
#else
			double K_ij = MXGet(i, j);
			double K_ji = MXGet(j, i);
#endif
			if (K_ij == 0.0 && K_ji == 0.0)
				continue;
			double d1 = GetFCTADiff(K_ij, K_ji);
#endif
			if (d1 == 0.0)
				continue;
			double d0 = d1; // TODO should use AuxMatrix at the previous time step
			// if (list_bc_nodes.find(i)!=list_bc_nodes.end() || list_bc_nodes.find(j)!=list_bc_nodes.end()) {
			//  d1 = d0 = 0.0;
			//}

			// Complete antidiffusive flux: f_ij += -theta*d_ij^H*DeltaU_ij^H - (1-theta)*d_ij^n*DeltaU_ij^n
			double diff_uH = this->GetNodeValue(i, idx1) - this->GetNodeValue(j, idx1);
			double diff_u0 = this->GetNodeValue(i, idx0) - this->GetNodeValue(j, idx0);
			double v = -(theta * d1 * diff_uH + (1.0 - theta) * d0 * diff_u0);
			(*FCT_AFlux)(i, j) += v;

			// prelimiting f
			v = (*FCT_AFlux)(i, j);
			if (this->m_num->fct_prelimiter_type == 0)
			{
				if (v * (-diff_uH) > 0.0)
					v = 0.0;
			}
			else if (this->m_num->fct_prelimiter_type == 1)
				v = MinMod(v, -d1 * diff_uH);
			else if (this->m_num->fct_prelimiter_type == 2)
				v = SuperBee(v, -d1 * diff_uH);
			(*FCT_AFlux)(i, j) = v;
#ifdef USE_PETSC
			(*FCT_AFlux)(j, i) = -v;
#else
			(*FCT_AFlux)(j, i) = v;
#endif

#ifdef USE_PETSC
			// A += theta * D
			if (i < (size_t)m_msh->getNumNodesLocal())
			{
				eqs_new->addMatrixEntry(i_global, i_global, -d1 * theta);
				eqs_new->addMatrixEntry(i_global, j_global, d1 * theta);
			}
			if (j < (size_t)m_msh->getNumNodesLocal())
			{
				eqs_new->addMatrixEntry(j_global, i_global, d1 * theta);
				eqs_new->addMatrixEntry(j_global, j_global, -d1 * theta);
			}
#else
// L = K + D
#if defined(NEW_EQS)
			(*A)(i, i) += -d1;
			(*A)(i, j) += d1;
			(*A)(j, i) += d1;
			(*A)(j, j) += -d1;
#else
			// add off-diagonal term
			MXInc(i, j, d1);
			MXInc(j, i, d1);
			// add diagonal term
			MXInc(i, i, -d1);
			MXInc(j, j, -d1);
#endif
#endif
		}
	}

	//----------------------------------------------------------------------
	// Assemble RHS: b_i += [- (1-theta) * L_ij] u_j^n
	//----------------------------------------------------------------------
	Math_Group::Vec* V1 = this->Gl_Vec1;
	Math_Group::Vec* V = this->Gl_Vec;
	(*V1) = 0.0;
	(*V) = 0.0;
#if !defined(USE_PETSC)
	double* eqs_rhs;
#ifdef NEW_EQS
	eqs_rhs = eqs_new->b;
#else
	eqs_rhs = eqs->b;
#endif
#endif
	// b = [-(1-theta) * L] u^n
	if (1.0 - theta > .0)
	{
		// u^n
		for (size_t i = 0; i < node_size; i++)
			(*V1)(i) = this->GetNodeValue(i, idx0);
		// L*u^n
		for (size_t i = 0; i < node_size; i++)
		{
#ifdef USE_PETSC
			const size_t i_global = FCT_GLOB_ADDRESS(i);
#endif
			for (size_t j = 0; j < node_size; j++)
			{
#ifdef USE_PETSC
				const size_t j_global = FCT_GLOB_ADDRESS(j);
				// b+=-(1-theta)*D*u^n
				(*V)(i) += (*FCT_d)(i_global, j_global) * (*V1)(j);
#else
#ifdef NEW_EQS
				(*V)(i) += (*A)(i, j) * (*V1)(j);
#else
				(*V)(i) += MXGet(i, j) * (*V1)(j);
#endif
#endif
			}
		}
		for (size_t i = 0; i < node_size; i++)
		{
#if defined(USE_PETSC)
			if (i < (size_t)m_msh->getNumNodesLocal())
			{
				const size_t i_global = FCT_GLOB_ADDRESS(i);
				eqs_new->add_bVectorEntry(i_global, -(1.0 - theta) * (*V)(i), ADD_VALUES);
			}
#else
			eqs_rhs[i] -= (1.0 - theta) * (*V)(i);
//(*RHS)(i+LocalShift) +=  NodalVal[i];
#endif
		}
	}

#ifndef USE_PETSC
	//----------------------------------------------------------------------
	// Assemble A matrix: 1/dt*ML + theta * L
	//----------------------------------------------------------------------
	// A matrix: theta * L
	if (theta == 0.0)
	{
#ifdef NEW_EQS
		(*A) = 0.0;
#else
		for (size_t i = 0; i < node_size; i++)
			for (size_t j = 0; j < node_size; j++)
				MXSet(i, j, 0.0);

#endif
	}
	else if (theta == 1.0)
	{
		// keep
	}
	else
	{
#ifdef NEW_EQS
		(*A) *= theta;
#else
		for (size_t i = 0; i < node_size; i++)
			for (size_t j = 0; j < node_size; j++)
				MXMul(i, j, theta);

#endif
	}
	// A matrix: += 1/dt * ML
	for (size_t i = 0; i < node_size; i++)
	{
		double v = 1.0 / dt * (*ML)(i);
#ifdef NEW_EQS
		(*A)(i, i) += v;
#else
		MXInc(i, i, v);
#endif
	}
#endif

	//----------------------------------------------------------------------
	// Assemble RHS: b += alpha * f
	//----------------------------------------------------------------------
	// Calculate R+, R-
	Math_Group::Vec* R_plus = this->Gl_Vec1;
	Math_Group::Vec* R_min = this->Gl_Vec;
	(*R_plus) = 0.0;
	(*R_min) = 0.0;
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_global = FCT_GLOB_ADDRESS(i);
		double P_plus, P_min;
		double Q_plus, Q_min;
		P_plus = P_min = 0.0;
		Q_plus = Q_min = 0.0;
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			if (i == j)
				continue;
			double f = (*jj).second; // double f = (*FCT_AFlux)(i,j);
#ifndef USE_PETSC
			if (i > j)
				f *= -1.0;
#endif
			double diff_uH = this->GetNodeValue(j, idx1) - this->GetNodeValue(i, idx1);

			P_plus += max(0.0, f);
			P_min += min(0.0, f);
			Q_plus = max(Q_plus, diff_uH);
			Q_min = min(Q_min, diff_uH);
		}
		double ml = (*ML)(i_global);

		if (P_plus == 0.0)
			(*R_plus)(i_global) = 0.0;
		else
			(*R_plus)(i_global) = min(1.0, ml * Q_plus / (dt * P_plus));
		if (P_min == 0.0)
			(*R_min)(i_global) = 0.0;
		else
			(*R_min)(i_global) = min(1.0, ml * Q_min / (dt * P_min));
	}

#ifdef USE_PETSC
	FCT_MPI::gatherR(FCT_MPI::ct, *R_plus, *R_min);
#endif

	// for Dirichlet nodes
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* bc_node = bc_node_value[i];
		long nod_id = bc_node->geo_node_number;
		const size_t i_global = FCT_GLOB_ADDRESS(nod_id);
		(*R_plus)(i_global) = 1.0;
		(*R_min)(i_global) = 1.0;
	}

	// b_i += alpha_i * f_ij
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_global = FCT_GLOB_ADDRESS(i);
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			const size_t j_global = FCT_GLOB_ADDRESS(j);
			if (i == j)
				continue;

			double f = (*jj).second; // double f = (*FCT_AFlux)(i,j);
#ifndef USE_PETSC
			if (i > j)
				f *= -1; // symmetric
#endif
			double alpha = 1.0;
			if (f > 0)
				alpha = min((*R_plus)(i_global), (*R_min)(j_global));
			else
				alpha = min((*R_plus)(j_global), (*R_min)(i_global));

			double val = .0;
			if (this->m_num->fct_const_alpha < 0.0)
				val = alpha * f;
			else
				val = this->m_num->fct_const_alpha * f;

#ifdef USE_PETSC
			if (i < (size_t)m_msh->getNumNodesLocal())
				eqs_new->add_bVectorEntry(i_global, val, ADD_VALUES);
#else
			eqs_rhs[i] += val;
#endif

			// Note: Galerkin FEM is recovered if alpha = 1 as below,
			// eqs_rhs[i] += 1.0*f;
		}
	}
}

/*************************************************************************
   GeoSys-Function:
   Task: Assemble the global system equation
   Programming:
   01/2005 WW/OK Implementation
   04/2005 OK MSH
   07/2005 WW Change due to the geometry element objects applied
   10/2005 OK DDC
   11/2005 YD time step control
   01/2006 OK/TK Tests
   12/2007 WW Spase matrix class and condensation sequential and parallel loop
   10/2010 TF changed access to process type
   06/2012 WW Node based decompostion
 **************************************************************************/
void CRFProcess::GlobalAssembly()
{
#ifdef USE_PETSC
	PetscLogDouble v1, v2;
#if (PETSC_VERSION_NUMBER > 3030)
	PetscTime(&v1);
#else
	PetscGetTime(&v1);
#endif
#endif

	// Tests
	if (!Tim)
		Tim = TIMGet(convertProcessTypeToString(this->getProcessType()));
	if (!Tim)
	{
		cout << "Error in CRFProcess::GlobalAssembly() - no TIM data"
		     << "\n";
		return;
	}

	if (!fem)
		// Which process needs this?
		// Only one instance of CFiniteElementStd is required for each process
		// Use "new" in such way will cause memory problem.
		// Please move this declaration to pcs configuration.     WW
		if (m_msh)
			fem = new CFiniteElementStd(this, m_msh->GetCoordinateFlag());

	CElem* elem = NULL;

	if (m_num->nls_method == 2) // 06.09.2010. WW
		IncorporateBoundaryConditions();
	bool Check2D3D;
	Check2D3D = false;
	if (type == 66) // Overland flow
		Check2D3D = true;
	if (this->femFCTmode) // NW
	{
		(*this->FCT_AFlux) = 0.0;
		(*this->Gl_ML) = 0.0;
		(*this->Gl_Vec) = 0.0;
		(*this->Gl_Vec1) = 0.0;
#ifdef USE_PETSC
		(*this->FCT_K) = 0.0;
		(*this->FCT_d) = .0;
#endif
	}

#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	// DDC
	if (dom_vector.size() > 0)
	{
		cout << "      Domain Decomposition" << '\n';
		CPARDomain* m_dom = NULL;
		size_t j = 0;
#if defined(USE_MPI) // WW
		j = myrank;
#else
		for (j = 0; j < dom_vector.size(); j++)
		{
#endif
		m_dom = dom_vector[j];
#ifdef NEW_EQS
		m_dom->InitialEQS(this);
#else
			SetLinearSolver(m_dom->eqs);
			SetZeroLinearSolver(m_dom->eqs);
#endif
		for (size_t i = 0; i < m_dom->elements.size(); i++)
		{
			elem = m_msh->ele_vector[m_dom->elements[i]];
			if (elem->GetMark() && elem->GetExcavState() == -1) // WX:10.2012
			{
				elem->SetOrder(false);
				// WW
				fem->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
				fem->ConfigElement(elem, Check2D3D);
				fem->m_dom = m_dom; // OK
				fem->Assembly();
			}
		}
		// m_dom->WriteMatrix();
		// MXDumpGLS("rf_pcs.txt",1,m_dom->eqs->b,m_dom->eqs->x);
		// ofstream Dum("rf_pcs.txt", ios::out); // WW
		// m_dom->eqs->Write(Dum);
		// Dum.close();
		IncorporateSourceTerms(j);
		if (m_num->nls_method != 2) // 06.09.2010. WW
			IncorporateBoundaryConditions(j);
/*
   //TEST
   string test = "rank";
   char stro[64];
   sprintf(stro, "%d",j);
   string test1 = test+(string)stro+"Assemble.txt";
   ofstream Dum(test1.c_str(), ios::out);
   m_dom->eqs->Write(Dum);
   Dum.close();
   exit(1);
 */
#ifndef USE_MPI
	}
	// Assemble global system
	DDCAssembleGlobalMatrix();
#endif
}
else
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
{ // STD
	// YDTEST. Changed to DOF 15.02.2007 WW
	for (size_t ii = 0; ii < continuum_vector.size(); ii++)
	{
		continuum = ii;
		for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
		{
			elem = m_msh->ele_vector[i];
			// Marked for use //WX: modified for coupled excavation
			if (elem->GetMark() && elem->GetExcavState() == -1)
			{
				elem->SetOrder(false);
				fem->ConfigElement(elem, Check2D3D);
				fem->Assembly();
				// NEUMANN CONTROL---------
				if (Tim->time_control_type == TimeControlType::NEUMANN)
				{
					Tim->time_step_length_neumann = MMin(Tim->time_step_length_neumann, timebuffer);
					Tim->time_step_length_neumann *= 0.5 * elem->GetVolume() * elem->GetVolume();
					if (Tim->time_step_length_neumann < MKleinsteZahl)
						Tim->time_step_length_neumann = 1.0e-5;
				}
				//------------------------------
			}
		}
	}

	if (femFCTmode) // NW
		AddFCT_CorrectionVector();

	//	          MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();
	// eqs_new->Write();
	IncorporateSourceTerms();

// MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();

#ifdef GEM_REACT
	if (getProcessType() == FiniteElement::MASS_TRANSPORT && aktueller_zeitschritt > 1)
	{ // JT->KG. New coupling system.
		if (_problem->GetCPLMaxIterations() > 1 || // Then there is a coupling on all processes
		    (this->pcs_is_cpl_overlord && this->m_num->cpl_max_iterations > 1)
		    || // This process (the overlord) controls coupling for another (the underling) process.
		    (this->pcs_is_cpl_underling
		     && this->cpl_overlord->m_num->cpl_max_iterations
		            > 1)) // This process (the underling) has a coupling with a controlling (the overlord) process
		{
			// IncorporateSourceTerms_GEMS();
		}
	}
#endif
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && !defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS                             //WW. 07.11.2008
	SetCPL(); // OK
#endif

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012.
	eqs_new->AssembleRHS_PETSc();
	eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);
#endif
	IncorporateBoundaryConditions();

// ofstream Dum("rf_pcs.txt", ios::out); // WW
// eqs_new->Write(Dum);   Dum.close();

#define nOUTPUT_EQS_BIN
#ifdef OUTPUT_EQS_BIN
	string fname = FileName + "_equation.bin";
	ofstream dum_bin(fname.c_str(), ios::out | ios::binary | ios::trunc);
	if (dum_bin.good())
	{
		eqs_new->Write_BIN(dum_bin);
		dum_bin.close();
		exit(1);
	}
#endif
//
//

//	MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); //abort();
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012.
	MPI_Barrier(MPI_COMM_WORLD);
//	eqs_new->AssembleRHS_PETSc();
// eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY );
#endif
}

#ifdef USE_PETSC
#if (PETSC_VERSION_NUMBER > 3030)
		PetscTime(&v2);
#else
		PetscGetTime(&v2);
#endif
cpu_time_assembly += v2 - v1;
#endif
}

//--------------------------------------------------------------------
/*! \brief Assmble eqiations
     for all PDEs excluding deformation;

     24.11.2010. WW
 */
void CRFProcess::GlobalAssembly_std(const bool is_mixed_order, bool Check2D3D)
{
	long i;
	CElem* elem = NULL;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark()) // Marked for use
			continue; // For OpenMP. WW

		elem->SetOrder(m_msh->getOrder());
		fem->setMixedOrderFlag(is_mixed_order);
		fem->ConfigElement(elem, Check2D3D);
		fem->Assembly();
	}
}

/*************************************************************************
   GeoSys-Function:
   Task: Integration
   Programming:
   05/2009 WW Implementation
 **************************************************************************/
void CRFProcess::Integration(vector<double>& node_velue)
{
	//----------------------------------------------------------------------
	size_t k;
	long i;
	CElem* elem = NULL;
	bool Check2D3D;
	Check2D3D = false;
	double n_val[8];

	if (type == 66) // Overland flow
		Check2D3D = true;

	vector<double> buffer((long)node_velue.size());
	for (i = 0; i < (long)buffer.size(); i++)
		buffer[i] = 0.;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark())
			continue;

		for (k = 0; k < elem->GetNodesNumber(false); k++)
			n_val[k] = node_velue[elem->GetNodeIndex(k)];

		elem->SetOrder(false);
		fem->ConfigElement(elem, Check2D3D);
		fem->FaceIntegration(n_val);

		for (k = 0; k < elem->GetNodesNumber(false); k++)
			buffer[elem->GetNodeIndex(k)] += n_val[k];
	}
	//----------------------------------------------------------------------
}

/*************************************************************************
   GeoSys-Function:
   Task: Calculate integration point velocity
   Programming:
   08/2005 WW Implementation
   last modified:
 **************************************************************************/
void CRFProcess::CalIntegrationPointValue()
{
	CElem* elem = NULL;
	continuum = 0; // 15.02.2007/
	//  cal_integration_point_value = true;
	// Currently, extropolation only valid for liquid and Richards flow.
	//	if (_pcs_type_name.find("LIQUID") != string::npos || _pcs_type_name.find(
	//			"RICHARD") != string::npos || _pcs_type_name.find(
	//			"MULTI_PHASE_FLOW") != string::npos || _pcs_type_name.find(
	//			"GROUNDWATER_FLOW") != string::npos || _pcs_type_name.find(
	//			"TWO_PHASE_FLOW") != string::npos
	//			|| _pcs_type_name.find("AIR_FLOW") != string::npos
	//			|| _pcs_type_name.find("PS_GLOBAL") != string::npos) //WW/CB

	switch (getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
		case FiniteElement::RICHARDS_FLOW:
		case FiniteElement::MULTI_PHASE_FLOW:
		case FiniteElement::GROUNDWATER_FLOW:
		case FiniteElement::TWO_PHASE_FLOW:
		case FiniteElement::DEFORMATION_H2: // 07.2011. WW
		case FiniteElement::AIR_FLOW:
		case FiniteElement::PS_GLOBAL:
		case FiniteElement::MULTI_COMPONENTIAL_FLOW: // AKS/NB
		case FiniteElement::DEFORMATION_FLOW: // NW
		case FiniteElement::TNEQ: // HS, TN
		case FiniteElement::TES: // HS, TN
		case FiniteElement::HEAT_TRANSPORT: // JOD 2014-11-10
		case FiniteElement::MASS_TRANSPORT: // JOD 2014-11-10
			cal_integration_point_value = true;
			break;
		default:
			cal_integration_point_value = false;
			return;
	}

	std::cout << "->Calculate velocity" << '\n';

	// check linear flow or not
	bool isLinearFlow = true;
	for (size_t i = 0; i < mmp_vector.size(); i++)
	{
		if (mmp_vector[i]->flowlinearity_model > 1)
		{
			isLinearFlow = false;
			break;
		}
	}
	if (isLinearFlow)
	{
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
		{
			elem = m_msh->ele_vector[i];
			if (elem->GetMark()) // Marked for use
			{
				if ((getProcessType() == FiniteElement::HEAT_TRANSPORT
				     || getProcessType() == FiniteElement::MASS_TRANSPORT) && !elem->selected)
					continue; // not selected for TOTAL_FLUX calculation JOD 2014-11-10
				fem->ConfigElement(elem);
				fem->Config(); // OK4709
				// fem->m_dom = NULL; // To be used for parallization
				if (getProcessType() == FiniteElement::MULTI_COMPONENTIAL_FLOW)
					fem->Cal_VelocityMCF();
				else
					fem->Cal_Velocity();

				// moved here from additional lower loop
				if (getProcessType() == FiniteElement::TNEQ || getProcessType() == FiniteElement::TES)
				{
					fem->CalcSolidDensityRate(); // HS, thermal storage reactions
				}
			}
		}
	}
	else
	{ // NW
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		const size_t v_itr_max(this->m_num->local_picard1_max_iterations);
		double pre_v[3] = {};
		double new_v[3] = {};
		// std::cout << "  Start local Picard iteration: tolerance = " << this->m_num->local_picard1_tolerance <<
		// std::endl;
		size_t i_itr = 0;
		double vel_error = .0;
		for (i_itr = 0; i_itr < v_itr_max; ++i_itr)
		{
			// std::cout << "  non-linear iteration: " << i_itr << "/" << v_itr_max << std::endl;
			vel_error = .0;
			for (size_t i = 0; i < mesh_ele_vector_size; i++)
			{
				elem = m_msh->ele_vector[i];
				if (elem->GetMark()) // Marked for use
				{
					ElementValue* gp_ele = ele_gp_value[i];
					gp_ele->GetEleVelocity(pre_v);

					fem->ConfigElement(elem);
					fem->Config(); // OK4709
					// fem->m_dom = NULL; // To be used for parallization

					fem->Cal_Velocity();

					gp_ele->GetEleVelocity(new_v);
					vel_error = max(vel_error, fabs(new_v[0] - pre_v[0]));
					vel_error = max(vel_error, fabs(new_v[1] - pre_v[1]));
					vel_error = max(vel_error, fabs(new_v[2] - pre_v[2]));
				}
			}
			// std::cout << "  error (max. norm): " << vel_error << std::endl;
			bool isConverged = (vel_error < this->m_num->local_picard1_tolerance);
			if (isConverged)
				break;
		}
		std::cout << "  Local Picard iteration: itr. count = " << i_itr << "/" << v_itr_max
		          << ", error(max. norm)=" << vel_error << std::endl;
	}

	// if (getProcessType() == FiniteElement::TNEQ || getProcessType() == FiniteElement::TEQ) {
	//    const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	//    for (size_t i = 0; i < mesh_ele_vector_size; i++)
	//    {
	//        elem = m_msh->ele_vector[i];
	//        if (elem->GetMark())                        // Marked for use
	//        {
	//            fem->ConfigElement(elem);
	//            fem->Config();
	//      fem->CalcSolidDensityRate(); // HS, thermal storage reactions
	//        }
	//    }
	//}

	//	if (_pcs_type_name.find("TWO_PHASE_FLOW") != string::npos) //WW/CB
	if (getProcessType() == FiniteElement::TWO_PHASE_FLOW || getProcessType() == FiniteElement::MASS_TRANSPORT
	    || getProcessType() == FiniteElement::HEAT_TRANSPORT) // WW/CB/JOD 2014-11-10
		cal_integration_point_value = false;
}

/*************************************************************************
   GeoSys-Function: CalGPVelocitiesfromFluidMomentum
   Task: Calculate gauss point velocities from fluid momentum solution
      extrapolate velocitiues from nodes to gauss points
   Programming:
   09/2009 SB Implementation

 **************************************************************************/
void CRFProcess::CalGPVelocitiesfromFluidMomentum()
{
	long i;
	MeshLib::CElem* elem = NULL;
	int i_ind[3];

	cout << "      CalGPVelocitiesfromFluidMomentum()"
	     << "\n";

	// Get fluid_momentum process
	CRFProcess* m_pcs_fm = PCSGet(FiniteElement::FLUID_MOMENTUM);

	//  check all possibilities for grid orientation (ccord_flag)
	int coordinateflag = this->m_msh->GetCoordinateFlag();
	// get index of velocity
	i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_X") + 1;
	if (coordinateflag == 11)
		// get index of velocity
		i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y") + 1;
	if (coordinateflag == 12)
		// get index of velocity
		i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;
	// get index of velocity
	i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y") + 1;
	if (coordinateflag == 22)
		// get index of velocity
		i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;
	// get index of velocity
	i_ind[2] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;

	if ((i_ind[0] < 0) || (i_ind[1] < 0) || (i_ind[2] < 0))
		cout << " Error - wrong index in Cal_GP_Velocity_FM "
		     << "\n";

	// Loop over all elements
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i]; // get element
		if (elem->GetMark()) // Marked for use
		{
			fem->ConfigElement(elem);
			fem->Cal_GP_Velocity_FM(i_ind);
		}
	} // end element loop
}

/*************************************************************************
   ROCKFLOW - Function: AllocateLocalMatrixMemory
   Task: As the function name
   Programming:
   01/2005 WW/OK Implementation
   06/2005 OK MSH
   last modified:
 **************************************************************************/
void CRFProcess::AllocateLocalMatrixMemory()
{
	long i;
	//----------------------------------------------------------------------
	int up_type = 0;
	if (!M_Process)
		up_type = 0;
	if (H_Process && M_Process)
	{
		if (type != 4 && type != 41)
			up_type = 1;
		else
		{
			if (type == 4)
				up_type = 3;
			if (type == 41)
				up_type = 4;
		}
	}
	if (!H_Process)
		up_type = 2;
	// SB for steady state element matrices in transport
	if (MASS_TRANSPORT_Process || T_Process)
		up_type = 5;
	//----------------------------------------------------------------------
	ElementMatrix* eleMatrix = NULL;
	CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		// CB THMBM
		// CB removed if condition due to conflict with deactivated subdomains
		//		if (elem->GetMark()) {     // Marked for use
		eleMatrix = new ElementMatrix();
		eleMatrix->AllocateMemory(elem, up_type);
		Ele_Matrices.push_back(eleMatrix);
		//		}
	}
}

/*************************************************************************
   FEMLib function
   Task: Assemble global system matrix
   Programming:
   05/2003 OK Implementation
   ??/???? WW Moved from AssembleSystemMatrixNew
   05/2006 WW Modified to enable dealing with the case of DOF>1
   06/2006 WW Take the advantege of sparse matrix to enhance simulation
   10/2007 WW Change for the new classes of sparse matrix and linear solver
 **************************************************************************/
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
void CRFProcess::DDCAssembleGlobalMatrix()
{
	int ii, jj, dof;
	long i, j, j0, ig, jg, ncol;
	CPARDomain* m_dom = NULL;
	long* nodes2node = NULL; // WW
	double *rhs = NULL, *rhs_dom = NULL;
	double a_ij;
	double b_i = 0.0;
	b_i = b_i; // OK411
	long no_dom_nodes;
	dof = pcs_number_of_primary_nvals; // WW
	ncol = 0; // WW
#ifndef USE_MPI
	int no_domains = (int)dom_vector.size();
	for (int k = 0; k < no_domains; k++)
	{
		m_dom = dom_vector[k];
#else
		m_dom = dom_vector[myrank];
#endif
// RHS
#if defined(NEW_EQS)
		rhs = eqs_new->b;
		if (type == 4)
			rhs_dom = m_dom->eqsH->b;
		else
			rhs_dom = m_dom->eqs->b;
#else
		rhs = eqs->b;
		rhs_dom = m_dom->eqs->b;
#endif

		no_dom_nodes = m_dom->nnodes_dom; // WW
		if (type == 4 || type == 41)
			no_dom_nodes = m_dom->nnodesHQ_dom; // WW
		if (type == 41)
			dof--;
		for (i = 0; i < no_dom_nodes; i++)
		{
			//------------------------------------------
			// Use the feature of sparse matrix of FEM
			// WW
			ig = m_dom->nodes[i];
			ncol = m_dom->num_nodes2_node[i];
			nodes2node = m_dom->node_conneted_nodes[i];
			for (j0 = 0; j0 < ncol; j0++)
			{
				j = nodes2node[j0];
				if (j >= no_dom_nodes)
					continue;
				jg = m_dom->nodes[j];
				//------------------------------------------
				// DOF loop ---------------------------WW
				for (ii = 0; ii < dof; ii++)
				{
					for (jj = 0; jj < dof; jj++)
					{
// get domain system matrix
#ifdef NEW_EQS // WW
						if (type == 4)
							a_ij = (*m_dom->eqsH->A)(i + no_dom_nodes * ii, j + no_dom_nodes * jj);
						else
							a_ij = (*m_dom->eqs->A)(i + no_dom_nodes * ii, j + no_dom_nodes * jj);
						(*eqs_new->A)(ig + Shift[ii], jg + Shift[jj]) += a_ij;
#else // ifdef  NEW_EQS
						SetLinearSolver(m_dom->eqs);
						a_ij = MXGet(i + no_dom_nodes * ii, j + no_dom_nodes * jj);
						// set global system matrix
						SetLinearSolver(eqs);
						MXInc(ig + Shift[ii], jg + Shift[jj], a_ij);
#endif
					}
				}
				// DOF loop ---------------------------WW
			}
			// set global RHS vector //OK
			for (ii = 0; ii < dof; ii++) // WW
				rhs[ig + Shift[ii]] += rhs_dom[i + no_dom_nodes * ii];
		}

		// Mono HM------------------------------------WW
		if (type != 41)
#ifndef USE_MPI
			continue;
#else
			return;
#endif
		no_dom_nodes = m_dom->nnodes_dom;
		long no_dom_nodesHQ = m_dom->nnodesHQ_dom;
		double a_ji = 0.0;
		for (i = 0; i < no_dom_nodes; i++)
		{
			ig = m_dom->nodes[i]; // WW
			ncol = m_dom->num_nodes2_node[i];
			nodes2node = m_dom->node_conneted_nodes[i];
			for (j0 = 0; j0 < ncol; j0++)
			{
				j = nodes2node[j0];
				jg = m_dom->nodes[j];
				for (ii = 0; ii < dof; ii++) // ww
				{
#if defined(NEW_EQS)
					// dom to global. WW
					a_ij = (*m_dom->eqsH->A)(i + no_dom_nodesHQ * dof, j + no_dom_nodesHQ * ii);
					a_ji = (*m_dom->eqsH->A)(j + no_dom_nodesHQ * ii, i + no_dom_nodesHQ * dof);
					(*eqs_new->A)(ig + Shift[ii], jg + Shift[problem_dimension_dm]) += a_ij;
					(*eqs_new->A)(jg + Shift[problem_dimension_dm], ig + Shift[ii]) += a_ji;
#else // if defined(NEW_EQS)
					// get domain system matrix
					SetLinearSolver(m_dom->eqs);
					a_ij = MXGet(i + no_dom_nodesHQ * dof, j + no_dom_nodesHQ * ii);
					a_ji = MXGet(j + no_dom_nodesHQ * ii, i + no_dom_nodesHQ * dof);
					// set global system matrix
					SetLinearSolver(eqs);
					MXInc(ig + Shift[ii], jg + Shift[problem_dimension_dm], a_ij);
					MXInc(jg + Shift[problem_dimension_dm], ig + Shift[ii], a_ji);
#endif
				}
			}
		}
		for (i = 0; i < no_dom_nodes; i++)
		{
			ig = m_dom->nodes[i];
			ncol = m_dom->num_nodes2_node[i];
			nodes2node = m_dom->node_conneted_nodes[i];
			for (j0 = 0; j0 < ncol; j0++)
			{
				j = nodes2node[j0];
				jg = m_dom->nodes[j];
				if (jg >= no_dom_nodes)
					continue;
// get domain system matrix
#if defined(NEW_EQS)
				// dom to global. WW
				a_ij = (*m_dom->eqsH->A)(i + no_dom_nodesHQ * dof, j + no_dom_nodesHQ * dof);
				(*eqs_new->A)(ig + Shift[problem_dimension_dm], jg + Shift[problem_dimension_dm]) += a_ij;
#else
				SetLinearSolver(m_dom->eqs);
				a_ij = MXGet(i + no_dom_nodesHQ * dof, j + no_dom_nodesHQ * dof);
				// set global system matrix
				SetLinearSolver(eqs);
				MXInc(ig + Shift[problem_dimension_dm], jg + Shift[problem_dimension_dm], a_ij);
#endif
			}
			//
			rhs[ig + Shift[problem_dimension_dm]] += rhs_dom[i + no_dom_nodesHQ * dof];
		}
// Mono HM------------------------------------WW
#ifndef USE_MPI
	}
#endif
}
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW

/*************************************************************************
   ROCKFLOW - Function:
   Task: Assemble system matrix
   Programming: 05/2003 OK Implementation
   ToDo: Prototyp function
   last modified:
 **************************************************************************/
void CRFProcess::AssembleSystemMatrixNew(void)
{
	switch (type)
	{
		case 1:
			// MakeGS_ASM_NEW(eqs->b,eqs->x,ddummy);
			// SMAssembleMatrix(eqs->b,eqs->x,ddummy,this);
			break;
		case 2: // MTM2
			break;
		case 3: // HTM
			break;
		case 5: // Gas flow
			break;
		case 11:
			break;
		case 12:
			break;
		case 13:
			break;
		case 14:
			break;
		case 66:
			// MakeGS_ASM_NEW(eqs->b,eqs->x,ddummy);
			// SMAssembleMatrix(eqs->b,eqs->x,ddummy,this);
			break;
		default:
			DisplayMsgLn("CalculateElementMatrices: no CalculateElementMatrices specified");
			abort();
	}
#ifdef PARALLEL
	DDCAssembleGlobalMatrix();
#else
	IncorporateSourceTerms();
	IncorporateBoundaryConditions();
#endif
	// SetLinearSolver(eqs);
	// MXDumpGLS("global_matrix_dd.txt",1,eqs->b,eqs->x);
}

/**************************************************************************
   FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
   Task: set PCS boundary conditions
   Programing:
   05/2006 WW Implementation
**************************************************************************/
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
void CRFProcess::SetBoundaryConditionSubDomain()
{
	int k;
	long i, j;
	CPARDomain* m_dom = NULL;
	CBoundaryConditionNode* m_bc_nv = NULL;
	CNodeValue* m_st_nv = NULL;
	//
	for (k = 0; k < (int)dom_vector.size(); k++)
	{
		m_dom = dom_vector[k];
		// BC
		for (i = 0; i < (long)bc_node_value.size(); i++)
		{
			m_bc_nv = bc_node_value[i];
			for (j = 0; j < (long)m_dom->nodes.size(); j++)
				if (m_bc_nv->geo_node_number == m_dom->nodes[j])
				{
					bc_node_value_in_dom.push_back(i);
					bc_local_index_in_dom.push_back(j);
					break;
				}
		}
		rank_bc_node_value_in_dom.push_back((long)bc_node_value_in_dom.size());
		// ST
		for (i = 0; i < (long)st_node_value.size(); i++)
		{
			m_st_nv = st_node_value[i];
			for (j = 0; j < (long)m_dom->nodes.size(); j++)
				if (m_st_nv->geo_node_number == m_dom->nodes[j])
				{
					st_node_value_in_dom.push_back(i);
					st_local_index_in_dom.push_back(j);
					break;
				}
		}
		rank_st_node_value_in_dom.push_back((long)st_node_value_in_dom.size());
	}
	long Size = (long)st_node_value.size();
	long l_index;
	for (i = 0; i < Size; i++)
	{
		l_index = st_node_value[i]->geo_node_number;
		if (node_connected_doms[l_index] > 0)
			st_node_value[i]->node_value /= (double)node_connected_doms[l_index];
	}
}

/**************************************************************************
   FEMLib-Method: CRFProcess::SetSTWaterGemSubDomain
   Task: set source/sink terms for GEMS-flow coupling
   Programing:
   05/2006 WW Implementation
   03/2010 KG44 modified to GEM
**************************************************************************/
void CRFProcess::SetSTWaterGemSubDomain(int myrank)
{
	int k;
	long i, j; // WW, dsize=0;
	CPARDomain* m_dom = NULL;
	long int m_stgem_nv = -1;
	//
	long Size = (long)Water_ST_vec.size();
	long l_index = -1;

	//	cout << "dom_vec_size: " << dom_vector.size() << "\n";
	//	for ( k=0;k< ( int ) dom_vector.size();k++ )
	//	{
	k = myrank; // do it for each domain only once!
	m_dom = dom_vector[k];
	// WW dsize=(long) m_dom->nodes.size();
	// ST
	for (i = 0; i < Size; i++)
	{
		m_stgem_nv = Water_ST_vec[i].index_node;
		for (j = 0; j < (long)m_dom->nodes.size(); j++)
			if (m_stgem_nv == m_dom->nodes[j])
			{
				// index for Water_ST_vec
				stgem_node_value_in_dom.push_back(i);
				// index for RHS
				stgem_local_index_in_dom.push_back(j);
				//	cout << "dom " << k <<  " i, j " << i << " " << j   << "\n";
			}
	}
	// only one element per domain!
	rank_stgem_node_value_in_dom.push_back((long)stgem_node_value_in_dom.size());
	//	cout << "dom " << k <<  " rank_stgem_node_value_in_dom " << (long) rank_stgem_node_value_in_dom[0]  << "\n";

	//	}

	for (i = 0; i < Size; i++)
	{
		l_index = Water_ST_vec[i].index_node;
		// cout << i << " " << node_connected_doms[l_index] << " " << "\n";
		// values for shared nodes are scaled
		Water_ST_vec[i].water_st_value /= (double)node_connected_doms[l_index];
	}
}

#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
/**************************************************************************
   FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
   Task: set PCS boundary conditions
   Programing:
   02/2004 OK Implementation
   ????    WW  and time curve
   04/2005 OK MSH
   05/2005 OK conditional BCs
   04/2006 WW Changes due to the geometry object applied
   04/2006 OK Conditions by PCS coupling OK
   05/2006 WW Re-implement
   05/2006 WW DDC
   10/2007 WW Changes for the new classes of sparse matrix and linear solver
   last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank)
{
	static long i;
	static double bc_value, fac = 1.0, time_fac = 1.0;
	long bc_msh_node = -1;
#ifndef USE_PETSC
	long bc_eqs_index;
#endif
	long shift;
	int interp_method = 0;
	int curve, valid = 0;
	int ii, idx0 = -1;
	CBoundaryConditionNode* m_bc_node; // WW
	CBoundaryCondition* m_bc; // WW
	CFunction* m_fct = NULL; // OK
	bool is_valid = false; // OK
#ifndef USE_PETSC
	bool onExBoundary = false; // WX
#endif
	bool excavated = false; // WX
	bool onDeactiveBoundary = true; // WX:09.2011
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	vector<int> bc_eqs_id;
	vector<double> bc_eqs_value;
#else
#if !defined(NEW_EQS)
	double* eqs_rhs = NULL;
#endif
	CPARDomain* m_dom = NULL;
#endif
//
#ifdef NEW_EQS
	Linear_EQS* eqs_p = NULL;
#endif
//------------------------------------------------------------WW
#ifdef JFNK_H2M
	if (m_num->nls_method == 2 && BC_JFNK.size() > 0) // 29.10.2010. WW
		return;

	/// For JFNK
	bool bc_inre_flag = false;
	double bc_init_value = 0.;
#endif

	// WW
	double Scaling = 1.0;
	if (type == 4 || type / 10 == 4)
	{
		fac = Scaling;
	}
	long begin = 0;
	long end = 0;
	long gindex = 0;
	if (rank == -1)
	{
		begin = 0;
		end = (long)bc_node_value.size();
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
// TODO
#ifdef NEW_EQS // WW
		eqs_p = eqs_new;
#else
			eqs_rhs = eqs->b;
#endif
#endif
	}
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
	else
	{
		m_dom = dom_vector[rank];
#ifdef NEW_EQS
		eqs_p = m_dom->eqs;
		if (type == 4) // WW
		{
			eqs_p = m_dom->eqsH;
		}
#else
			eqs_rhs = m_dom->eqs->b;
#endif
		if (rank == 0)
			begin = 0;
		else
			begin = rank_bc_node_value_in_dom[rank - 1];
		end = rank_bc_node_value_in_dom[rank];
	}
#endif // END: #if !defined(USE_PETSC) // && !defined(other parallel libs)

	// update nod velocity before eval of constrained BC
	if ((this->hasConstrainedBC() || this->hasConstrainedST())
	    && (this->getProcessType() == FiniteElement::RICHARDS_FLOW
	        || this->getProcessType() == FiniteElement::LIQUID_FLOW))
	{
		this->CalIntegrationPointValue();
		this->Extropolation_GaussValue();
	}

	for (i = begin; i < end; i++)
	{
		gindex = i;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
		if (rank > -1)
			gindex = bc_node_value_in_dom[i];
#endif
		m_bc_node = bc_node_value[gindex];
		m_bc = bc_node[gindex];
		//
		// WX: check if bc is aktive, when Time_Controlled_Aktive for this bc is defined
		if (m_bc->getTimeContrCurve() > 0)
			if (GetCurveValue(m_bc->getTimeContrCurve(), 0, aktuelle_zeit, &valid) < MKleinsteZahl)
				continue;

		// WX: 01.2011. for excavation bc, check if excavated and if on boundary
		if (m_bc->getExcav() > 0)
		{
			CNode* node;
// unsigned int counter;	//void warning
#ifndef USE_PETSC
			onExBoundary = false; // WX:01.2011
#endif
			excavated = false;

			node = m_msh->nod_vector[m_bc_node->geo_node_number];
			double const* node_coordinate(node->getData()); // Coordinates(node_coordinate);

			if ((node_coordinate[ExcavDirection] >= ExcavBeginCoordinate
			     && (GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) + ExcavBeginCoordinate)
			            >= node_coordinate[ExcavDirection])
			    || (node_coordinate[ExcavDirection] < ExcavBeginCoordinate
			        && (GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) + ExcavBeginCoordinate)
			               < node_coordinate[ExcavDirection]))
			{
				excavated = true;
#ifndef USE_PETSC
				for (unsigned int j = 0; j < node->getConnectedElementIDs().size(); j++)
				{
					CElem* elem = m_msh->ele_vector[node->getConnectedElementIDs()[j]];
					double const* tmp_ele_coor(elem->GetGravityCenter());
					// if(elem->GetPatchIndex()!=ExcavMaterialGroup){
					// if(elem->GetExcavState()==-1)
					if (elem->GetPatchIndex() != static_cast<size_t>(ExcavMaterialGroup))
					// to be improved for more than 1 MG excavation
					{
						onExBoundary = true;
						break;
					}
					else if (tmp_ele_coor[ExcavDirection]
					             - (GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) - ExcavBeginCoordinate)
					         > -0.001)
					{
						onExBoundary = true;
						// tmp_counter1++;
						break;
					}
				}
#endif
			}
		}

		if ((m_bc->getExcav() > 0) && !excavated) // WX:01.2011. excav bc but is not excavated jet
			continue;
		if (m_bc->getGeoTypeAsString().compare("DOMAIN") == 0) // WX
		{
			CNode* node;
			CElem* elem;
			onDeactiveBoundary = false;
			node = m_msh->nod_vector[m_bc_node->geo_node_number];
			for (size_t j = 0; j < node->getConnectedElementIDs().size(); j++)
			{
				elem = m_msh->ele_vector[node->getConnectedElementIDs()[j]];
				if (elem->GetMark())
				{
					onDeactiveBoundary = true;
					break;
				}
			}
		}
		if (!onDeactiveBoundary)
			continue;
//
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
		bc_msh_node = m_bc_node->geo_node_number;
		// Check whether the node is in this subdomain
		const long id_act_l_max = m_msh->getNumNodesLocal();
		const long id_h = static_cast<long>(m_msh->GetNodesNumber(false));
		const long id_act_h_max = static_cast<long>(m_msh->getLargestActiveNodeID_Quadratic());
		// If ghost node
		if (!((bc_msh_node < id_act_l_max) || ((bc_msh_node >= id_h) && (bc_msh_node < id_act_h_max))))
			continue;

		int dof_per_node = 0;
		if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
		{
			dof_per_node = pcs_number_of_primary_nvals;
			shift = m_bc_node->msh_node_number / m_msh->NodesNumber_Linear;
		}
		else
		{
			if (bc_msh_node < static_cast<long>(m_msh->NodesNumber_Linear))
				dof_per_node = pcs_number_of_primary_nvals;
			else
				dof_per_node = m_msh->GetCoordinateFlag() / 10;
			if (type == 4 || type / 10 == 4)
			{
				shift = m_bc_node->msh_node_number / m_msh->NodesNumber_Quadratic;
			}
			else
			{
				shift = m_bc_node->msh_node_number / m_msh->NodesNumber_Linear;
			}
		}
#else
		shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
		if (rank > -1)
		{
			bc_msh_node = bc_local_index_in_dom[i];
			int dim_space = 0;
			if (shift == 0)
				// 15.4.2008 WW if(m_msh->NodesNumber_Linear==m_msh->NodesNumber_Quadratic)
				dim_space = 0;
			else
			{
				// 02.2010. WW
				if (type == 4 || type / 10 == 4)
					dim_space = shift / m_msh->NodesNumber_Quadratic;
				else
					dim_space = shift / m_msh->NodesNumber_Linear;
			}
			shift = m_dom->shift[dim_space];
		}
		else
			bc_msh_node = m_bc_node->geo_node_number;
#endif // END: if defined(USE_PETSC) // || defined(other parallel libs
		//------------------------------------------------------------WW
		if (m_msh) // OK
			//	    if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
			//          continue;
			time_fac = 1.0;
		if (bc_msh_node >= 0)
		{
			//................................................................
			// Time dependencies - CURVE
			curve = m_bc_node->CurveIndex;
			if (curve > 0)
			{
				if (curve > 10000000) /// 16.08.2010. WW
					time_fac = GetCurveValue(curve - 10000000, interp_method, aktuelle_zeit, &valid);
				else
					time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
				if (!valid)
					continue;
			}
			else
				time_fac = 1.0;
			//................................................................
			// Time dependencies - FCT
			if (m_bc_node->fct_name.length() > 0)
			{
				m_fct = FCTGet(m_bc_node->fct_name);
				if (m_fct)
					time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
				// if(!valid) continue;
				else
					cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data"
					     << "\n";
			}
			//................................................................
			// copy values SB 09.2012
			if (m_bc_node->bc_node_copy_geom.length() > 0)
			{
				size_t geo_node_id;
				// std::string geo_file_name = (std::string)this->m_msh->_geo_name;
				GEOLIB::GEOObjects* geo_obj;
				geo_obj = this->m_msh->getGEOObjects();
				// Get geometry node
				const std::vector<GEOLIB::Point*>* pnt_vec(geo_obj->getPointVec(*this->m_msh->getProjectName()));
				// const GEOLIB::PointVec* pnt_vec(geo_obj->getPointVecObj(m_bc_node->bc_node_copy_geom));
				// const GEOLIB::Point* pnt(pnt_vec->getElementIDByName(m_bc_node->bc_node_copy_geom, msh_node_id));
				if (!((geo_obj->getPointVecObj(*this->m_msh->getProjectName()))
				          ->getElementIDByName(m_bc_node->bc_node_copy_geom_name, geo_node_id)))
				{
					cout << " Error !"
					     << "\n";
				}
				// get node Index of geometry node - only POINT so far
				long msh_node_id = m_msh->GetNODOnPNT((*pnt_vec)[geo_node_id]);
				// get value from last time step
				std::string help;
				help = m_bc_node->pcs_pv_name;
				if (this->pcs_type_name_vector[0] == "MASS_TRANSPORT")
					help = this->pcs_primary_function_name[0];
				int ind = this->GetNodeValueIndex(help) + 1; // new time level
				time_fac = this->GetNodeValue(msh_node_id, ind);
				std::cout << " Copy value " << time_fac << " from " << m_bc_node->bc_node_copy_geom << " "
				          << m_bc_node->bc_node_copy_geom_name << " to POLYLINE "
				          << "\n";
			}
			//................................................................
			// Conditions
			if (m_bc_node->conditional)
			{
				int idx_1 = -1; // 28.2.2007 WW
				for (ii = 0; ii < dof; ii++) // 28.2.2007 WW

					if (convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable())
					        .find(pcs_primary_function_name[ii]) != string::npos)
					{
						idx_1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
						break;
					}
				bc_value = time_fac * fac * GetNodeValue(m_bc_node->msh_node_number_subst, idx_1);
			}
			else
			{
				if (m_bc->getPressureAsHeadModel() == -1) // this is the default case
				{
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				}
				else // this is the PressureAsHead case
				{
					bc_value
					    = calcPressureFromHead(*m_bc, m_bc_node->msh_node_number, m_bc_node->node_value, time_fac, fac);
				}
			}
			if (m_bc->isSwitchBC())
			{
				bc_value = evaluteSwitchBC(*m_bc, *m_bc_node, time_fac, fac);
			}
			//----------------------------------------------------------------
			// MSH
			/// 16.08.2010. WW
			if (curve > 10000000 && fabs(time_fac) > DBL_EPSILON)
				bc_value = bc_value / time_fac + time_fac; /// bc_value +time_fac;

			if (m_bc->isPeriodic()) // JOD
				bc_value
				    *= sin(2 * 3.14159 * aktuelle_zeit / m_bc->getPeriodeTimeLength() + m_bc->getPeriodePhaseShift());

//----------------------------------------------------------------
#ifndef USE_PETSC
			if (rank > -1)
				bc_eqs_index = bc_msh_node;
			else
				// WW#
				bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
#endif
			//..............................................................
			// NEWTON WW   //Modified for JFNK. 09.2010. WW
			if (m_num->nls_method >= 1 // 04.08.2010. WW _name.find("NEWTON")!=string::npos
			    || type == 4
			    || type / 10 == 4)
			{ // Solution is in the manner of increment !
				idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()));
				if (m_bc_node->pcs_pv_name.find("DISPLACEMENT") != string::npos)
				{
					bc_value -= GetNodeValue(m_bc_node->geo_node_number, idx0)
					            + GetNodeValue(m_bc_node->geo_node_number, idx0 + 1);
					/// if JFNK and if the first Newton step
					/// In the successive Newton steps, node BC value taken from the previous Newton step.
					if (m_num->nls_method == 2 && ite_steps == 1)
						SetNodeValue(m_bc_node->geo_node_number, idx0, bc_value);
					if (m_bc->getNoDispIncre() > 0) // WX
						bc_value = 0;

#ifdef JFNK_H2M
					bc_inre_flag = false;
					bc_init_value = 0.;
#endif
				}
				else
				{
#ifdef JFNK_H2M
					bc_inre_flag = true;
					bc_init_value = bc_value;
#endif
					/// if JFNK and if the first Newton step. 11.11.2010. WW
					// if(m_num->nls_method==2&&ite_steps==1) /// JFNK
					/// p_{n+1} = p_b,
					//  SetNodeValue(m_bc_node->geo_node_number, idx0++, bc_value);
					/// dp = u_b-u_n
					bc_value -= GetNodeValue(m_bc_node->geo_node_number, ++idx0);
					// SetNodeValue(m_bc_node->geo_node_number, idx0, bc_value);
					// SetNodeValue(m_bc_node->geo_node_number, idx0+1, bc_value);
					// bc_value = 0.;
				}
			}
#if !defined(USE_PETSC) // && !defined(other parallel solver). //WW 04.2012. WW
			bc_eqs_index += shift;
#endif

#ifdef JFNK_H2M
			/// If JFNK method (09.2010. WW):
			if (m_num->nls_method == 2)
			{
				bc_JFNK new_bc_entry;
				new_bc_entry.var_idx = idx0 + 1;
				new_bc_entry.bc_node = m_bc_node->geo_node_number;
				new_bc_entry.bc_eqs_idx = bc_eqs_index;
				new_bc_entry.bc_value = bc_value;

				new_bc_entry.incremental = bc_inre_flag;
				new_bc_entry.bc_value0 = bc_init_value;
				BC_JFNK.push_back(new_bc_entry);
			}
			else
			{
#endif
//----------------------------------------------------------------
//----------------------------------------------------------------
/* // Make the follows as comment by WW. 04.03.2008
   //YD dual
   if(dof>1) //WW
   {
   for(ii=0;ii<dof;ii++)
   {
    if(m_bc->pcs_pv_name.find(pcs_primary_function_name[ii]) != string::npos)
    {
       //YD/WW
       //WW   bc_eqs_index += ii*(long)eqs->dim/dof;   //YD dual //DOF>1 WW
       bc_eqs_index += fem->NodeShift[ii];   //DOF>1 WW

   break;
   }
   }
   }
 */

#if !defined(USE_PETSC) // && !defined(other parallel libs)//04.3013. WW
				if (m_bc->getExcav() > 0)
				{
					if (excavated && !onExBoundary)
					{
#ifdef NEW_EQS
						eqs_p->SetKnownX_i(bc_eqs_index, bc_value); // WX:12.2012
#else
							MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
						continue; // WX:09.2011
					}
					else if (m_bc_node->pcs_pv_name.find("DISPLACEMENT") != string::npos)
						continue; //
				}
#endif // end defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW

				// save last bc_value for later usage (in constrained BC)
				m_bc_node->node_value_pre_calc = bc_value;

				if (m_bc->isConstrainedBC())
				{
					if (checkConstrainedBC(*m_bc, *m_bc_node, bc_value))
						continue;
				}
//////////////////////////////////

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
				bc_eqs_id.push_back(
				    static_cast<int>(m_msh->nod_vector[bc_msh_node]->GetEquationIndex() * dof_per_node + shift));
				bc_eqs_value.push_back(bc_value);

#elif defined(NEW_EQS) // WW
			eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
			MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
#ifdef JFNK_H2M
			}
#endif
		}
	}
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	int nbc = static_cast<int>(bc_eqs_id.size());
	if (nbc > 0)
	{
		eqs_new->setArrayValues(0, nbc, &bc_eqs_id[0], &bc_eqs_value[0], INSERT_VALUES);
		eqs_new->setArrayValues(1, nbc, &bc_eqs_id[0], &bc_eqs_value[0], INSERT_VALUES);
	}

#ifdef petsc_zero_row_test
	// We have do the following collection because MatZeroR must be called by all processes
	const int mpi_size = eqs_new->getMPI_Size();
	vector<int> r_cnt(mpi_size);
	vector<int> r_disp(mpi_size);
	vector<int> r_vec(mpi_size);
	int k;
	for (k = 0; k < mpi_size; k++)
	{
		r_cnt[k] = 1;
		r_disp[k] = k;
	}
	// Get nbc
	MPI_Allgatherv(&nbc, 1, MPI_INT, &r_vec[0], &r_cnt[0], &r_disp[0], MPI_INT, PETSC_COMM_WORLD);
	int v_disp = 0;
	for (k = 0; k < mpi_size; k++)
	{
		r_disp[k] = v_disp;
		v_disp += r_vec[k];
	}
	r_cnt.resize(v_disp);
	r_vec.resize(v_disp);
	for (k = 0; k < v_disp; k++)
	{
		r_cnt[k] = 0;
	}
	const int v_shift = r_disp[eqs_new->getMPI_Rank()];
	for (k = 0; k < nbc; k++)
	{
		r_cnt[v_shift + k] = bc_eqs_id[k];
	}
	MPI_Allreduce(&r_cnt[0], &r_vec[0], v_disp, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	eqs_new->zeroRows_in_Matrix(v_disp, &r_vec[0]);
#endif //  petsc_zero_row_test
	eqs_new->AssembleUnkowns_PETSc();
	eqs_new->AssembleRHS_PETSc();

	// TEST
	// PetscViewer viewer;
	// eqs_new->EQSV_Viewer(FileName, viewer);

	eqs_new->zeroRows_in_Matrix(nbc, &bc_eqs_id[0]);
	eqs_new->AssembleMatrixPETSc();

#endif
}

/**************************************************************************
   FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
   Task: set PCS boundary conditions for FLUID_MOMENTUM depending on axis
   Programing:
   01/2007 PCH Implementation
   last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank, const int axis)
{
	static long i;
	static double bc_value, fac = 1.0, time_fac = 1.0;
	long bc_msh_node;
	long bc_eqs_index, shift;
	int interp_method = 0;
	int curve, valid = 0;
	int idx0, idx1;
	CBoundaryConditionNode* m_bc_node; // WW
	CBoundaryCondition* m_bc; // WW
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	CPARDomain* m_dom = NULL;
#endif

	CFunction* m_fct = NULL; // OK
	bool is_valid = false; // OK
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	vector<int> bc_eqs_id;
	vector<double> bc_eqs_value;
#elif !defined(NEW_EQS)
	double* eqs_rhs = NULL;
#endif
#ifdef NEW_EQS
	Linear_EQS* eqs_p = NULL;
#endif
	//------------------------------------------------------------WW
	// WW
	double Scaling = 1.0;
	if (type == 4 || type == 41)
		fac = Scaling;

	long begin = 0;
	long end = 0;
	long gindex = 0;

	// WW CBoundaryConditionsGroup *m_bc_group = NULL;
	//  m_bc_group = BCGetGroup(this->_pcs_type_name,this->pcs_primary_function_name[axis]);
	// TF
	// WW m_bc_group =
	// BCGetGroup(convertProcessTypeToString(this->getProcessType()),this->pcs_primary_function_name[axis]);

	if (rank == -1)
	{
		begin = 0;
		end = (long)bc_node_value.size();
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
#ifdef NEW_EQS // WW
		eqs_p = eqs_new;
#else
			eqs_rhs = eqs->b;
#endif
#endif
	}
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
	else
	{
		m_dom = dom_vector[rank];
#ifdef NEW_EQS
		eqs_p = m_dom->eqs;
		if (type == 4) // WW
		{
			eqs_p = m_dom->eqsH;
		}
#else
			eqs_rhs = m_dom->eqs->b;
#endif
		if (rank == 0)
			begin = 0;
		else
			begin = rank_bc_node_value_in_dom[rank - 1];
		end = rank_bc_node_value_in_dom[rank];
	}
#endif // END: !defined(USE_PETSC) // && !defined(other parallel libs)/
	for (i = begin; i < end; i++)
	{
		gindex = i;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
		if (rank > -1)
			gindex = bc_node_value_in_dom[i];
#endif
		m_bc_node = bc_node_value[gindex];

		// PCH
		if (axis == 0 && m_bc_node->pcs_pv_name.find("VELOCITY1_X") != string::npos)
		{
			m_bc = bc_node[gindex];

#if defined(USE_PETSC) // ||defined(other parallel libs)//03~04.3012. WW
			bc_msh_node = m_bc_node->geo_node_number;

#else
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
			//
			if (rank > -1)
			{
				bc_msh_node = bc_local_index_in_dom[i];
				int dim_space = 0;
				if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
					dim_space = 0;
				else
				{
					if (shift % m_msh->NodesNumber_Quadratic == 0)
						dim_space = shift / m_msh->NodesNumber_Quadratic;
					else
						dim_space = m_msh->msh_max_dim;
				}
				shift = m_dom->shift[dim_space];
			}
			else
				bc_msh_node = m_bc_node->geo_node_number;
#endif
			//------------------------------------------------------------WW
			if (m_msh) // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
					if (!valid)
						continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data"
						     << "\n";
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value
					    = time_fac * fac
					      * GetNodeValue(m_bc_node->msh_node_number_subst,
					                     // WW  bc_value = time_fac*fac*
					                     // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					                     // // YD-----TEST---
					                     GetNodeValueIndex(pcs_primary_function_name[0]) + 1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (m_num->nls_method >= 1 //_name.find("NEWTON")!=string::npos
				    || type == 4
				    || type == 41) // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()));
					if (type == 4 || type == 41)
					{
						// 30.08.2010. WW
						if (m_bc_node->pcs_pv_name.find("PRESSURE") == string::npos)
						{
							idx1 = idx0 + 1;
							bc_value -= GetNodeValue(m_bc_node->geo_node_number, idx0)
							            + GetNodeValue(m_bc_node->geo_node_number, idx1);
						}
					}
					else
						bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == string::npos)
						continue;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif defined(NEW_EQS) // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
				MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
			}
		}
		// PCH
		else if (axis == 1 && m_bc_node->pcs_pv_name.find("VELOCITY1_Y") != string::npos)
		{
			m_bc = bc_node[gindex];
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
//
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
			if (rank > -1)
			{
				bc_msh_node = bc_local_index_in_dom[i];
				int dim_space = 0;
				if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
					dim_space = 0;
				else
				{
					if (shift % m_msh->NodesNumber_Quadratic == 0)
						dim_space = shift / m_msh->NodesNumber_Quadratic;
					else
						dim_space = m_msh->msh_max_dim;
				}
				shift = m_dom->shift[dim_space];
			}
			else
#endif
				bc_msh_node = m_bc_node->geo_node_number;
			//------------------------------------------------------------WW
			if (m_msh) // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
					if (!valid)
						continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data"
						     << "\n";
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value
					    = time_fac * fac
					      * GetNodeValue(m_bc_node->msh_node_number_subst,
					                     // WW  bc_value = time_fac*fac*
					                     // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					                     // // YD-----TEST---
					                     GetNodeValueIndex(pcs_primary_function_name[0]) + 1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (m_num->nls_method >= 1 //_name.find("NEWTON")!=string::npos
				    || type == 4
				    || type == 41) // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).c_str());
					if (type == 4 || type == 41)
					{
						idx1 = idx0 + 1;
						bc_value -= GetNodeValue(m_bc_node->geo_node_number, idx0)
						            + GetNodeValue(m_bc_node->geo_node_number, idx1);
					}
					else
						bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == string::npos)
						continue;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif defined(NEW_EQS) // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
				MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
			}
		}
		// PCH
		else if (axis == 2 && m_bc_node->pcs_pv_name.find("VELOCITY1_Z") != string::npos)
		{
			m_bc = bc_node[gindex];
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
//
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
			if (rank > -1)
			{
				bc_msh_node = bc_local_index_in_dom[i];
				int dim_space = 0;
				if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
					dim_space = 0;
				else
				{
					if (shift % m_msh->NodesNumber_Quadratic == 0)
						dim_space = shift / m_msh->NodesNumber_Quadratic;
					else
						dim_space = m_msh->msh_max_dim;
				}
				shift = m_dom->shift[dim_space];
			}
			else
#endif
				bc_msh_node = m_bc_node->geo_node_number;
			//------------------------------------------------------------WW
			if (m_msh) // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
					if (!valid)
						continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in CRFProcess::IncorporateBoundaryConditions - no FCT data"
						     << "\n";
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value
					    = time_fac * fac
					      * GetNodeValue(m_bc_node->msh_node_number_subst,
					                     // WW  bc_value = time_fac*fac*
					                     // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					                     // // YD-----TEST---
					                     GetNodeValueIndex(pcs_primary_function_name[0]) + 1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (m_num->nls_method >= 1 //_name.find("NEWTON")!=std::string::npos
				    || type == 4
				    || type == 41) // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(convertPrimaryVariableToString(m_bc->getProcessPrimaryVariable()).c_str());
					if (type == 4 || type == 41)
					{
						idx1 = idx0 + 1;
						bc_value -= GetNodeValue(m_bc_node->geo_node_number, idx0)
						            + GetNodeValue(m_bc_node->geo_node_number, idx1);
					}
					else
						bc_value = bc_value - GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(pcs_primary_function_name[continuum]) == std::string::npos)
						continue;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif defined(NEW_EQS) // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
				MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
			}
		}
	}

	//-----------------------------------------------------------------------
	/* irreg. Zeilen/Spalten regularisieren */
	/*
	   else if (GetNodeState(NodeNumber[i]) == -2 || GetNodeState(NodeNumber[i]) == -4) { // irreg.Knoten
	    if (GetRFControlGridAdapt())
	      if (AdaptGetMethodIrrNodes() == 1) {
	        MXSet(i, i, MKleinsteZahl);
	        rechts[i] = 0.0;
	      }
	    }
	   }
	 */
}

bool CRFProcess::checkConstrainedST(std::vector<CSourceTerm*>& st_vector, CSourceTerm const& st,
                                    CNodeValue const& st_node)
{
	bool return_value(false);
	for (std::size_t i = 0; i < st.getNumberOfConstrainedSTs(); i++)
	{
		bool constrained_bool(false);
		const Constrained& local_constrained(st.getConstrainedST(i));

		if (local_constrained.constrainedPrimVar == FiniteElement::PRESSURE
		    || local_constrained.constrainedPrimVar == FiniteElement::CONCENTRATION
		    || local_constrained.constrainedPrimVar == FiniteElement::TEMPERATURE)
		{
			// other process
			CRFProcess* pcs(PCSGet(local_constrained.constrainedProcessType));

			// other variable
			for (std::size_t k = 0; k < pcs->GetPrimaryVNumber(); k++)
			{
				if (pcs->GetPrimaryVName(k)
				    == FiniteElement::convertPrimaryVariableToString(local_constrained.constrainedPrimVar))
				{
					// value of PrimVar of other process at current node
					double local_value = pcs->GetNodeValue(st_node.geo_node_number, 2 * k + 1);

					// FIXME: if getPressureAsHead != -1, then the head is expected
					// however, getPressureAsHead is part of CBoundaryCondition (m_bc in IncorporateBoundaryConditions)

					if (local_constrained.constrainedDirection
					    == ConstrainedType::GREATER) // exclude greater and equal values
					{
						if (local_value
						    >= local_constrained
						           .constrainedValue) // check if calculated value (eg of other process) meets criterium
							constrained_bool = true;
					}
					else if (local_constrained.constrainedDirection
					         == ConstrainedType::SMALLER) // exclude smaller values
					{
						if (local_value < local_constrained.constrainedValue)
							constrained_bool = true;
					}
					/*else	is already checked when reading
					return false;*/
				}
			}
		}
		st_vector[st_node.getSTVectorGroup()]->setConstrainedSTNode(i, constrained_bool, st_node.getSTVectorIndex());
		if (constrained_bool)
			return_value = true;
	}

	return return_value;
}

bool CRFProcess::checkConstrainedBC(CBoundaryCondition const& bc, CBoundaryConditionNode& bc_node, double& bc_value)
{
	for (std::size_t i = 0; i < bc.getNumberOfConstrainedBCs(); i++)
	{
		const Constrained& local_constrained(bc.getConstrainedBC(i));
		if (local_constrained.constrainedVariable == ConstrainedVariable::VELOCITY)
		{
			// get velocity vector at node
			std::valarray<double> vel(0.0, 3); // ini valarray with size=3 and element values=0

			if (local_constrained._isConstrainedVelStable)
			{
				const std::size_t no_connected_nodes
				    = m_msh->nod_vector[bc_node.geo_node_number]->getNumConnectedNodes();
				for (std::size_t j = 0; j < no_connected_nodes; j++)
				{
					std::size_t const connected_node_id
					    = m_msh->nod_vector[bc_node.geo_node_number]->getConnectedNodes()[j];
					if (connected_node_id == static_cast<std::size_t>(bc_node.geo_node_number))
					{
						std::valarray<double> temp_vel(this->getNodeVelocityVector(connected_node_id));
						temp_vel *= static_cast<double>(no_connected_nodes - 1);
						vel += temp_vel;
					}
					else
						vel += this->getNodeVelocityVector(connected_node_id);
				}
			}
			else
				vel += this->getNodeVelocityVector(bc_node.geo_node_number);

			// nomalize velocity vector
			double magn_vel_v(MBtrgVec(&vel[0], 3));
			if (!(magn_vel_v > 0)) // check if velocity is not zero
			{
#ifdef DEBUG
				std::cout << "No constrained applied at node " << bc_node.msh_node_number
				          << " as magnitude of velocity " << magn_vel_v << " is not > than 0 " << std::endl;
#endif
				continue;
			}

			// NormalizeVector(&vel[0], 3); can use other way to calc normalized vector, as vector now valarray:
			vel /= magn_vel_v;

			// calculate scalar product of velocity vector and BC surface normal
			double const scalar_prod(MathLib::scpr(&(vel)[0], bc_node.GetNormalVector(), 3));

			// select case to handle BC
			if (scalar_prod < 0.05 && scalar_prod >= -1.01 // small tolerance of ~5°
			    && local_constrained.constrainedDirection == ConstrainedType::NEGATIVE) // velocity vector and bc
			                                                                            // surface normal point in
			                                                                            // opposite directions
			{
				return true; // do not apply BC (maybe later implementation: change BC to ST at this node)
			}
			else if (scalar_prod <= 1.01 && scalar_prod >= -0.05 // small tolerance of ~5°
			         && local_constrained.constrainedDirection == ConstrainedType::POSITIVE) // velocity vector points
			                                                                                 // in same direction as bc
			                                                                                 // surface normal
			{
				return true; // do not apply BC (maybe later implementation: change BC to ST at this node)
			}
		}
		if (local_constrained.constrainedPrimVar == FiniteElement::PRESSURE
		    || local_constrained.constrainedPrimVar == FiniteElement::CONCENTRATION
		    || local_constrained.constrainedPrimVar == FiniteElement::TEMPERATURE)
		{
			// other process
			CRFProcess* pcs(PCSGet(local_constrained.constrainedProcessType));

			// other variable
			for (std::size_t k = 0; k < pcs->GetPrimaryVNumber(); k++)
			{
				if (pcs->GetPrimaryVName(k)
				    == FiniteElement::convertPrimaryVariableToString(local_constrained.constrainedPrimVar))
				{
					// TODO: check case GREATER - maybe not thoroughly enough evaluated
					if (local_constrained.constrainedDirection
					    == ConstrainedType::GREATER) // exclude greater and equal values
					{
						// value of PrimVar of other process at current node
						double local_value(0.1 * pcs->GetNodeValue(bc_node.geo_node_number, 2 * k + 1)
						                   + 0.9 * pcs->GetNodeValue(bc_node.geo_node_number, 2 * k));
						if (local_value
						    >= local_constrained
						           .constrainedValue) // check if calculated value (eg of other process) meets criterium
							return true;
						if (bc.getProcessType() == local_constrained.constrainedProcessType
						    && bc_node.node_value_pre_calc
						           >= local_constrained.constrainedValue) // check if BC value meets criterium
							return true;
					}
					else if (local_constrained.constrainedDirection
					         == ConstrainedType::SMALLER) // exclude smaller values
					{
						if (bc.getProcessType()
						    == local_constrained.constrainedProcessType) // check if correct process type
						{
							if (bc_node.node_value_pre_calc
							    < local_constrained.constrainedValue) // limit BC if constrained below BC value
							{
								if (local_constrained._isSeepageBC)
								{
									// value of PrimVar of other process at current node
									double local_value(
									    0.1
									        * pcs->GetNodeValue(
									              bc_node.geo_node_number,
									              2 * k + 1) // slight averaging with last time step for stability
									    + 0.9 * pcs->GetNodeValue(bc_node.geo_node_number, 2 * k));

									if (local_value < local_constrained.constrainedValue)
										return true;
									else
										bc_value = local_constrained.constrainedValue;
								}
								else
									return true;
							}
						}
					}
				}
			}
		}
	}
	std::cout << "Non existing combination for constrained BC direction given. Using normal BC." << std::endl;
	return false;
}

/**************************************************************************
Copied & modified from
rf_kinreact.cpp (Reaction-Method:)
**************************************************************************/
std::valarray<double> CRFProcess::getNodeVelocityVector(const long node_id)
{
	CRFProcess* m_pcs = NULL;
	std::valarray<double> vel_nod(0.0, 3);

	m_pcs = PCSGetFlow();

	// Get the velocity components
	vel_nod[0] = m_pcs->GetNodeValue(node_id, this->_idxVx);
	vel_nod[1] = m_pcs->GetNodeValue(node_id, this->_idxVy);
	vel_nod[2] = m_pcs->GetNodeValue(node_id, this->_idxVz);

	// Shift entries for 2D
	const int dimensions(m_msh->GetCoordinateFlag());
	if (dimensions == 22)
	{
		vel_nod[2] = vel_nod[1];
		vel_nod[1] = 0;
	}
	else if (dimensions == 23)
	{
		vel_nod[2] = vel_nod[1];
		vel_nod[1] = vel_nod[0];
		vel_nod[0] = 0;
	}

	return vel_nod;
}

int CRFProcess::getFirstNodeBelowGWL(size_t current_node)
{
	// searches for PRESSURE1<0 and returns the node below or above node i

	int sorted_node(-1);
	for (size_t i(0); i < this->m_msh->sorted_nodes.size(); i++)
	{
		if (current_node == this->m_msh->sorted_nodes[i])
		{
			sorted_node = i;
			break;
		}
	}

	const int val_idx = this->GetNodeValueIndex("PRESSURE1", true);
	int new_GWL_st_node(-1);

	for (size_t k(1); k < this->m_msh->xy_change.size(); k++) // look for range of xy-coordinates
	{
		const long node_id = this->m_msh->xy_change[k];
		if (node_id >= sorted_node) // found range
		{
			new_GWL_st_node = this->m_msh->sorted_nodes[this->m_msh->xy_change[k - 1] + 1]; // set to bottom node first

			if (this->GetNodeValue(this->m_msh->sorted_nodes[node_id], val_idx) > 0) // top node has p>0
			{
				new_GWL_st_node = this->m_msh->sorted_nodes[node_id]; // return top node id
				break;
			}
			else
			{
				for (size_t j(node_id - 1); j > this->m_msh->xy_change[k - 1] + 1;
				     j--) // look in range from top to bottom (but skip top node)
				{
					if (this->GetNodeValue(this->m_msh->sorted_nodes[j], val_idx) > 0) // look for first node with p>0
					{
						new_GWL_st_node
						    = this->m_msh->sorted_nodes[j - 1]; // for p<0, return node id two nodes below gwl
						break;
					}
				}
			}
			break;
		}
	}

	return new_GWL_st_node;
}

// MW if set, calculate pressure from given head at bc with density, gravity constant and geodetic height for PRESSURE
// as primary variable
double CRFProcess::calcPressureFromHead(CBoundaryCondition const& bc, std::size_t node_number, double head_value,
                                        double const& time_fac, double const& fac)
{
	double local_density;
	switch (bc.getPressureAsHeadModel())
	{
		case 0:
			// use current density at node
			local_density = MFPGetNodeValue(node_number, "DENSITY", 0);
			break;
		case 1:
			// use given density
			local_density = bc.getPressureAsHeadDensity();
			break;
		default:
			std::cout
			    << "Warning! No PressureAsHeadDensity specified. Calculating density (i.e. PressureAsHeadModel 0)!"
			    << std::endl;
			local_density = MFPGetNodeValue(node_number, "DENSITY", 0);
			break;
	}

	const double local_node_elevation = this->m_msh->nod_vector[node_number]->Z();

	// pressure = gravitational_constant * density * ( head - geodetic_height )
	double pressure_value(fac * gravity_constant * local_density * (time_fac * head_value - local_node_elevation));

	return pressure_value;
}

double CRFProcess::calcHeadFromPressure(CBoundaryCondition const& bc, std::size_t node_number, double pressure_value,
                                        double const& time_fac, double const& fac)
{
	double local_density;
	switch (bc.getPressureAsHeadModel())
	{
		case 0:
			// use current density at node
			local_density = MFPGetNodeValue(node_number, "DENSITY", 0);
			break;
		case 1:
			// use given density
			local_density = bc.getPressureAsHeadDensity();
			break;
		default:
			std::cout
			    << "Warning! No PressureAsHeadDensity specified. Calculating density (i.e. PressureAsHeadModel 0)!"
			    << std::endl;
			local_density = MFPGetNodeValue(node_number, "DENSITY", 0);
			break;
	}

	const double local_node_elevation = this->m_msh->nod_vector[node_number]->Z();

	// pressure = gravitational_constant * density * ( head - geodetic_height )
	double head_value(fac * (pressure_value * time_fac / (gravity_constant * local_density) + local_node_elevation));

	return head_value;
}

double CRFProcess::evaluteSwitchBC(CBoundaryCondition const& bc, CBoundaryConditionNode const& bc_node, double time_fac,
                                   double fac)
{
	const SwitchBC& localSwitchBC(bc.getSwitchBC());
	// assume, normally to use on value
	double returnValue(time_fac * fac * localSwitchBC.switchOnValue);

	CRFProcess* pcs(PCSGet(localSwitchBC.switchProcessType));

	// other variable
	for (std::size_t i(0); i < pcs->GetPrimaryVNumber(); i++)
	{
		if (pcs->GetPrimaryVName(i) == FiniteElement::convertPrimaryVariableToString(localSwitchBC.switchPrimVar))
		{
			// value of PrimVar of other process at current node
			double localValue(0.1 * pcs->GetNodeValue(bc_node.geo_node_number, 2 * i + 1)
			                  + 0.9 * pcs->GetNodeValue(bc_node.geo_node_number, 2 * i));

			if (localValue < localSwitchBC.switchValue)
			{
				// set off value, if other value lower than switchValue
				returnValue = (time_fac * fac * localSwitchBC.switchOffValue);
			}
			break;
		}
	}
	return returnValue;
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS source terms into EQS
   Programing:
   04/2004 OK Implementation
   08/2004 WW Extension for monolithic PCS and time curve
   last modification:
   02/2005 MB River Condition and CriticalDepth
   05/2005 WW Dynamic problems
   07/2005 WW Changes due to the geometry object applied
   03/2006 WW Re-arrange
   04/2006 OK CPL
   05/2006 WW DDC
   08/2006 YD FCT use
**************************************************************************/
void CRFProcess::IncorporateSourceTerms(const int rank)
{
	double value = 0, fac = 1.0, time_fac;
	int interp_method = 0;
	int curve, valid = 0;
	long msh_node, shift;
	MshElemType::type EleType; // ii
	double q_face = 0.0;
	CElem* elem = NULL;
	CElem* face = NULL;
	ElementValue* gp_ele = NULL;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	vector<int> st_eqs_id;
	vector<double> st_eqs_value;
#else
	CPARDomain* m_dom = NULL;
	double* eqs_rhs = NULL;
	long bc_eqs_index = -1;
	int dim_space = 0; // kg44 better define here and not in a loop!
#endif
	double vel[3];
	bool is_valid; // YD
	CFunction* m_fct = NULL; // YD
	long i; //, group_vector_length;

	double Scaling = 1.0;
	if (type == 4 || type / 10 == 4)
	{
		fac = Scaling;
	}

	CNodeValue* cnodev = NULL;
	CSourceTerm* m_st = NULL;
	//
	long begin = 0;
	long end = 0;
	long gindex = 0;

	//###############################
	// NB Climate Data
	// MW Use loop for constrained ST evaluation

	for (size_t i = 0; i < st_vector.size(); i++)
	{
		// for constrainedST
		if (!(st_vector[i]->isConstrainedST()))
			continue;

		for (std::size_t j(0); j < st_vector[i]->getNumberOfConstrainedSTs(); j++)
		{
			if (!(st_vector[i]->isCompleteConstrainST(j)))
				continue;

			st_vector[i]->setCompleteConstrainedSTStateOff(false, j);
			size_t end = st_vector[i]->getNumberOfConstrainedSTNodes(j);
			for (std::size_t k(0); k < end; k++)
			{
				if (st_vector[i]->getConstrainedSTNode(j, k))
				{
					st_vector[i]->setCompleteConstrainedSTStateOff(true, j);
					break;
				}
			}
		}

		// NOTE (KR): This only works correctly if there is only ONE source term with DisType CLIMATE! TODO
		// If more are needed pls let me know
		if (st_vector[i]->getProcessDistributionType() == FiniteElement::CLIMATE)
		{
			m_st = st_vector[i];
			const MathLib::InverseDistanceInterpolation<GEOLIB::PointWithID*, GEOLIB::Station*>* distances(
			    m_st->getClimateStationDistanceInterpolation());

			// Interpolate each ST to every surface node;
			const std::vector<GEOLIB::Station*>& weather_stations(m_st->getClimateStations());
			const size_t nStations = weather_stations.size();
			vector<double> vec_etr(nStations);
			for (size_t j = 0; j < nStations; j++)
				vec_etr[j] = static_cast<GEOLIB::Station*>(weather_stations[j])
				                 ->getSensorData()
				                 ->getData(SensorDataType::RECHARGE, aktuelle_zeit, true);

			// Interpolate the point data onto the Mesh surface
			const size_t nSTNodeValues(st_node_value.size());
			const size_t nDistances(distances->getNDomainPoints());
			for (size_t j = 0; j < nSTNodeValues; j++)
			{
				// search the first node value index not set (this should be the start for the CLIMATE data
				if (st_node_value[j]->node_value == std::numeric_limits<double>::min())
				{
					// Interpolate for each surface node
					for (size_t w = 0; w < nDistances; w++)
					{
						// for each ST_node, the distance to all weather stations has been determined before
						// (InverseDistranceInterpolation)
						// Value ETR(n) at node n by inverse distance weighting of ETR(WS) for each weather station WS
						double st_value(0);
						for (size_t q = 0; q < nStations; q++)
							st_value += vec_etr[q] * distances->getDistance(w, q) / distances->getSumOfDistances(w);

						const size_t n(distances->getDomainPoint(w)->getID());
						st_node_value[j + w]->node_value = st_value * m_msh->nod_vector[n]->patch_area;
					}
					break;
				}
			}
#ifndef NDEBUG
			/* output of interpolated data */
			ofstream file;
			stringstream filename("interpolation");

			filename << FilePath << filename.str();

			if (Tim->step_current < 1000)
			{
				filename << 0;
				if (Tim->step_current < 100)
				{
					filename << 0;
					if (Tim->step_current < 10)
						filename << 0;
				}
			}

			filename << Tim->step_current << ".csv";
			cout << "      Interpolation result written in " << filename.str() << "\n";

			file.open(filename.str().c_str());
			file << "node,x,y,z,recharge,area"
			     << "\n";

			for (size_t i = 0; i < st_node_value.size(); i++)
			{
				const long z = st_node_value[i]->msh_node_number;
				const double* coords(m_msh->nod_vector[z]->getData());
				file << z << "," << coords[0] << "," << coords[1] << ","
				     << "0"
				     << "," << st_node_value[i]->node_value << "," << m_msh->nod_vector[z]->patch_area << "\n";
			}

			file.close();
#endif
		}
	}

	m_st = NULL;

	//#####################

	if (rank == -1)
	{
		begin = 0;
		end = (long)st_node_value.size();
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
#ifdef NEW_EQS // WW
		eqs_rhs = eqs_new->b; // 27.11.2007 WW
#else
			eqs_rhs = eqs->b;
#endif
#endif
	}
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
	else
	{
		m_dom = dom_vector[rank];
#ifdef NEW_EQS
		if (type == 4)
			eqs_rhs = m_dom->eqsH->b;
		else
			eqs_rhs = m_dom->eqs->b;
#else
			eqs_rhs = m_dom->eqs->b;
#endif
		if (rank == 0)
			begin = 0;
		else
			begin = rank_st_node_value_in_dom[rank - 1];
		end = rank_st_node_value_in_dom[rank];
	}
#endif // END: #if !defined(USE_PETSC) // && !defined(other parallel libs)

	for (i = begin; i < end; i++)
	{
		gindex = i;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
		if (rank > -1)
			gindex = st_node_value_in_dom[i];
#endif

		cnodev = st_node_value[gindex];

		// MW Recharge applied at GW Surface (p>0), tested only with LIQUID_FLOW or RICHARDS_FLOW and mmp
		// $PERMEABILITY_SATURATION = 10
		if (cnodev->getProcessDistributionType() == FiniteElement::RECHARGE
		    || cnodev->getProcessDistributionType() == FiniteElement::RECHARGE_DIRECT)
		{
			cnodev->msh_node_number = getFirstNodeBelowGWL(cnodev->msh_node_number);
			cnodev->geo_node_number = cnodev->msh_node_number;
		}

#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
		msh_node = cnodev->geo_node_number;
		// Check whether the node is in this subdomain
		const long id_act_l_max = m_msh->getNumNodesLocal();
		const long id_h = static_cast<long>(m_msh->GetNodesNumber(false));
		const long id_act_h_max = static_cast<long>(m_msh->getLargestActiveNodeID_Quadratic());

		// If ghost node
		if (!((msh_node < id_act_l_max) || ((msh_node >= id_h) && (msh_node < id_act_h_max))))
			continue;

		int dof_per_node = 0;
		if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
		{
			dof_per_node = pcs_number_of_primary_nvals;
			shift = cnodev->msh_node_number / m_msh->NodesNumber_Linear;
		}
		else
		{
			if (msh_node < static_cast<long>(m_msh->NodesNumber_Linear))
				dof_per_node = pcs_number_of_primary_nvals;
			else
				dof_per_node = m_msh->GetCoordinateFlag() / 10;
			if (type == 4 || type / 10 == 4)
			{
				shift = cnodev->msh_node_number / m_msh->NodesNumber_Quadratic;
			}
			else
			{
				shift = cnodev->msh_node_number / m_msh->NodesNumber_Linear;
			}
		}

#else
		shift = cnodev->msh_node_number - cnodev->geo_node_number;
		if (rank > -1)
		{
			msh_node = st_local_index_in_dom[i];
			dim_space = 0;
			if (m_msh->NodesNumber_Linear == m_msh->NodesNumber_Quadratic)
				dim_space = 0;
			else
			{
				if (shift % m_msh->NodesNumber_Quadratic == 0)
					dim_space = shift / m_msh->NodesNumber_Quadratic;
				else
					dim_space = m_msh->msh_max_dim;
			}
			shift = m_dom->shift[dim_space];
		}
		else
		{
			msh_node = cnodev->msh_node_number;
			msh_node -= shift;
		}
#endif
		value = cnodev->node_value;
		//--------------------------------------------------------------------
		// Tests
		if (msh_node < 0)
			continue;
		m_st = NULL;
		if (st_node.size() > 0 && (long)st_node.size() > i)
		{
			m_st = st_node[gindex];

			if (m_st->isConstrainedST())
			{
				bool continue_bool(false);
				continue_bool = checkConstrainedST(st_vector, *m_st, *cnodev);

				for (std::size_t temp_i(0); temp_i < m_st->getNumberOfConstrainedSTs(); temp_i++)
				{
					if (st_vector[cnodev->getSTVectorGroup()]->isCompleteConstrainST(temp_i)
					    && st_vector[cnodev->getSTVectorGroup()]->getCompleteConstrainedSTStateOff(temp_i))
					{
						continue_bool = true;
						break;
					}
				}

				if (continue_bool)
					continue;
			}

			//--------------------------------------------------------------------
			// CPL
			// if(m_st->_pcs_type_name_cond.size()>0) continue; // this is a CPL source term, JOD removed
			//--------------------------------------------------------------------
			// system dependent YD
			if (cnodev->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
			{
				long no_st_ele = (long)m_st->element_st_vector.size();
				for (long i_st = 0; i_st < no_st_ele; i_st++)
				{
					long ele_index = m_st->element_st_vector[i_st];
					elem = m_msh->ele_vector[ele_index];
					if (elem->GetMark())
					{
						fem->ConfigElement(elem);
						fem->Config();

						if (getProcessType() == FiniteElement::MULTI_COMPONENTIAL_FLOW)
							fem->Cal_VelocityMCF();
						else
							fem->Cal_Velocity();
					}
					gp_ele = ele_gp_value[ele_index];
					gp_ele->GetEleVelocity(vel);
					EleType = elem->GetElementType();
					if (EleType == MshElemType::LINE) // Line
						cnodev->node_value += vel[0];
					// Traingle & Qua
					if (EleType == MshElemType::TRIANGLE || EleType == MshElemType::QUAD)
					{
						for (size_t i_face = 0; i_face < m_msh->face_vector.size(); i_face++)
						{
							face = m_msh->face_vector[i_face];
							if ((size_t)m_st->element_st_vector[i_st] == face->GetOwner()->GetIndex())
								//
								q_face = PointProduction(vel, m_msh->face_normal[i_face]) * face->GetVolume();
							// for(i_node)
						}
						cnodev->node_value = +q_face / 2;
					}
				}
			}
			//--------------------------------------------------------------------
			// MB
			// if(m_st->conditional && !m_st->river)
			//{

			GetNODValue(value, cnodev, m_st);
		} // st_node.size()>0&&(long)st_node.size()>i
		//----------------------------------------------------------------------------------------
		//--------------------------------------------------------------------
		// Please do not move the this section
		curve = cnodev->CurveIndex;
		if (curve > 0)
		{
			// Reading Time interpolation method; BG
			if (m_st != NULL) // in some cases the m_st is not defined -> interp_method is not changed for this cases
				if (interp_method != m_st->TimeInterpolation)
					interp_method = m_st->TimeInterpolation;

			time_fac = GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
			// cout << "step: " << this->Tim->step_current << " Time: " << aktuelle_zeit << " Laenge: " <<
			// this->Tim->this_stepsize << " Beginn: " << this->Tim->time_start << " Ende " << this->Tim->time_end << "
			// Faktor: " << time_fac << "\n";
			if (!valid)
			{
				cout << "\n!!! Time dependent curve is not found. Results are not guaranteed "
				     << "\n";
				cout << " in void CRFProcess::IncorporateSourceTerms(const double Scaling)"
				     << "\n";
				time_fac = 1.0;
			}
		}
		else
			time_fac = 1.0;

		// Time dependencies - FCT    //YD
		if (m_st) // WW
		{
			// WW/YD //OK
			if (m_msh && m_msh->geo_name.find("LOCAL") != string::npos)
			{
				if (m_st->getFunctionName().length() > 0)
				{
					m_fct = FCTGet(pcs_number);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit,
						                           &is_valid,
						                           m_st->getFunctionMethod()); // fct_method. WW
					else
						cout << "Warning in CRFProcess::IncorporateSourceTerms - no FCT data"
						     << "\n";
				}
			}
			else if (m_st->getFunctionName().length() > 0)
			{
				m_fct = FCTGet(m_st->getFunctionName());
				if (m_fct)
					time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
				else
					cout << "Warning in CRFProcess::IncorporateSourceTerms - no FCT data"
					     << "\n";
			}
		}
		//----------------------------------------------------------------------------------------
		value *= time_fac * fac;
//------------------------------------------------------------------
// EQS->RHS
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW

		st_eqs_id.push_back(static_cast<int>(m_msh->nod_vector[msh_node]->GetEquationIndex() * dof_per_node + shift));
		st_eqs_value.push_back(value);

#else
		if (rank > -1)
			bc_eqs_index = msh_node + shift;
		else
			bc_eqs_index = m_msh->nod_vector[msh_node]->GetEquationIndex() + shift;
		eqs_rhs[bc_eqs_index] += value;
#endif
	}
	//====================================================================

	// if coupling to GEMS exist----------------------------------------------------
	// HS, added 11.2008
	// KG44 03/03/2010 modified to hopefully soon work with parallel solvers
	long gem_node_index = -1, glocalindex = -1;
	if (flag_couple_GEMS == 1 && aktueller_zeitschritt > 1)
	{
		begin = 0;
		if (rank == -1) // serial version and also Version for PETSC!!

			end = (long)Water_ST_vec.size();
		else // parallel version
		{
			end = 0;
			if (rank_stgem_node_value_in_dom.size() > 0)
				end = rank_stgem_node_value_in_dom[0];
		}
		// only when switch is on and not in the first time step
		// loop over the Water_ST_vec vector,
		// add the excess water to the right-hand-side of the equation
		for (i = begin; i < end; i++)
		{
			if (rank > -1) // parallel version: stgem_node_value_in_dom and stgem_local_index_in_dom contain only values
			               // for the corresponding domain  == rank
			{
				//				cout << "rank " << rank ;
				gindex = stgem_node_value_in_dom[i]; // contains indexes to water-st_vec
				//				cout << " gindex " << gindex << " i " << i << endl ;
				// contains index to node
				glocalindex = stgem_local_index_in_dom[i];
//				cout << " gem_node_index " << gem_node_index << "\n";
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW

				st_eqs_id.push_back(static_cast<int>(m_msh->nod_vector[glocalindex]->GetEquationIndex()));
				st_eqs_value.push_back(Water_ST_vec[gindex].water_st_value);
#else
				eqs_rhs[glocalindex] += Water_ST_vec[gindex].water_st_value;
#endif
			}
			else // serial version
			{
				gem_node_index = Water_ST_vec[i].index_node;
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
				st_eqs_id.push_back(static_cast<int>(m_msh->nod_vector[gem_node_index]->GetEquationIndex()));
				st_eqs_value.push_back(Water_ST_vec[i].water_st_value);

#else
				eqs_rhs[gem_node_index] += Water_ST_vec[i].water_st_value;
#endif
			}
		}
		// after finished adding to RHS, clear the vector
		Water_ST_vec.clear();
		if (rank > -1)
		{
			stgem_node_value_in_dom.clear();
			stgem_local_index_in_dom.clear();
			rank_stgem_node_value_in_dom.clear();
		}
	}
#if defined(USE_PETSC) // || defined(other parallel libs)//03~04.3012. WW
	if (st_eqs_id.size() > 0)
	{
		eqs_new->setArrayValues(1, static_cast<int>(st_eqs_id.size()), &st_eqs_id[0], &st_eqs_value[0]);
		// eqs_new->AssembleRHS_PETSc();
	}
#endif
}

#if !defined(USE_PETSC) && !defined(NEW_EQS) // || defined(other parallel libs)//03~04.3012.
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2004 OK Implementation
   11/2004 OK NUM
   07/2006 WW Parallel BiCGStab
   last modification:
**************************************************************************/
int CRFProcess::ExecuteLinearSolver(void)
{
	long iter_count = 0;
	long iter_sum = 0;
	// WW  int found = 0;
	//-----------------------------------------------------------------------
	// Set EQS
	// cout << "Before SetLinearSolver(eqs) myrank = "<< myrank<< '\n';
	SetLinearSolver(eqs);
	//--------------------------------------------------------------------
	// NUM
	// WW  found = 1;
	cg_maxiter = m_num->ls_max_iterations; // OK lsp->maxiter;
	cg_eps = m_num->ls_error_tolerance; // OK lsp->eps;
	// cg_repeat = lsp->repeat;
	vorkond = m_num->ls_precond; // OK lsp->precond;
	linear_error_type = m_num->ls_error_method; // OK lsp->criterium;

//  cout << "Before eqs->LinearSolver(eqs) myrank = "<< myrank<< '\n';

#ifdef USE_MPI
	// WW
	long dim_eqs = 0;
	if (type == 41 || type == 4) // DOF >1
	{
		dom_vector[myrank]->quadratic = true;
		if (type == 4)
			dim_eqs = pcs_number_of_primary_nvals * m_msh->GetNodesNumber(true);
		else if (type == 41)
			dim_eqs = pcs_number_of_primary_nvals * m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false);
	}
	else
	{
		dom_vector[myrank]->quadratic = false;
		dim_eqs = m_msh->GetNodesNumber(false);
	}
	iter_count = SpBICGSTAB_Parallel(dom_vector[myrank], eqs->x, dim_eqs);
#else // ifdef USE_MPI

		iter_count = eqs->LinearSolver(eqs->b, eqs->x, eqs->dim);
#endif

	eqs->master_iter = iter_count;
	if (iter_count >= cg_maxiter)
	{
		cout << "Warning in CRFProcess::ExecuteLinearSolver() - Maximum iteration number reached"
		     << "\n";
		return -1;
	}
	iter_sum += iter_count;
	//-----------------------------------------------------------------------
	// Clean results ?
	/*
	   for (i=0;i<eqs->dim;i++)
	    if (fabs(eqs->x[i])<MKleinsteZahl)
	      eqs->x[i] = 0.0;
	 */
	//-----------------------------------------------------------------------
	return iter_sum;
}
#endif
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2005 PCH Overriding
   last modification:
**************************************************************************/
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW 07.11.2008
int CRFProcess::ExecuteLinearSolver(LINEAR_SOLVER* eqs)
{
	long iter_count;
	long iter_sum = 0;
	// WW int found = 0;
	//-----------------------------------------------------------------------
	// Set EQS
	SetLinearSolver(eqs);
	//--------------------------------------------------------------------
	// NUM
	// WW found = 1;
	cg_maxiter = m_num->ls_max_iterations; // OK lsp->maxiter;
	cg_eps = m_num->ls_error_tolerance; // OK lsp->eps;
	// cg_repeat = lsp->repeat;
	vorkond = m_num->ls_precond; // OK lsp->precond;
	linear_error_type = m_num->ls_error_method; // OK lsp->criterium;

	iter_count = eqs->LinearSolver(eqs->b, eqs->x, eqs->dim);
	eqs->master_iter = iter_count;
	if (iter_count >= cg_maxiter)
	{
		cout << "Warning in CRFProcess::ExecuteLinearSolver() - Maximum iteration number reached"
		     << "\n";
		return -1;
	}
	iter_sum += iter_count;
	//-----------------------------------------------------------------------
	return iter_sum;
}
#endif
// WW
int CRFProcess::GetNODValueIndex(const string& name, int timelevel)
{
	for (int j = 0; j < number_of_nvals; j++)
		if ((name.compare(pcs_nval_data[j].name) == 0) && (pcs_nval_data[j].timelevel == timelevel))
			return pcs_nval_data[j].nval_index;
	cout << "Error in PCSGetNODValueIndex: " << name << "\n";
	return -1;
}

///////////////////////////////////////////////////////////////////////////
// Specials
///////////////////////////////////////////////////////////////////////////

/*-------------------------------------------------------------------------
   ROCKFLOW - Function: PCSRestart
   Task: Insert process to list
   Programming:
   06/2003 OK Implementation
   11/2004 OK file_name_base
   last modified:
   -------------------------------------------------------------------------*/
void PCSRestart()
{
	/*OK411
	   int j;
	   CRFProcess *m_pcs = NULL;
	   int nidx0,nidx1;
	   int i;
	   int no_processes =(int)pcs_vector.size();
	   if(no_processes==0)
	    return; //OK41
	   int ok = 0;
	   //----------------------------------------------------------------------
	   string file_name_base = pcs_vector[0]->file_name_base;
	   //OK  ok = ReadRFRRestartData(file_name_base);
	   if(ok==0){
	   cout << "RFR: no restart data" << "\n";
	   return;
	   }
	   //----------------------------------------------------------------------
	   for(i=0;i<no_processes;i++){
	   m_pcs = pcs_vector[i];
	   for(j=0;j<m_pcs->GetPrimaryVNumber();j++) {
	   // timelevel=0;
	   nidx0 = m_pcs->GetNodeValueIndex(m_pcs->GetPrimaryVName(j));
	   // timelevel= 1;
	   nidx1 = nidx0+1;
	   CopyNodeVals(nidx1,nidx0);
	   }
	   }
	 */
}

/**************************************************************************
   FEMLib-Method:
   Task: Relocate Deformation process
   Programing:
   09/2004 WW Implementation
   10/2010 TF changes due to conversion from std::string to enum for process type
**************************************************************************/
void RelocateDeformationProcess(CRFProcess* m_pcs)
{
	//   string pcs_name_dm = m_pcs->_pcs_type_name;
	FiniteElement::ProcessType pcs_name_dm(m_pcs->getProcessType());

	string num_type_name_dm;
	// Numerics
	if (m_pcs->num_type_name.compare("STRONG_DISCONTINUITY") == 0)
	{
		num_type_name_dm = m_pcs->num_type_name;
		enhanced_strain_dm = 1;
	}

	delete m_pcs;
	m_pcs = dynamic_cast<CRFProcess*>(new CRFProcessDeformation());
	m_pcs->setProcessType(pcs_name_dm);

	if (enhanced_strain_dm == 1)
		m_pcs->num_type_name = num_type_name_dm;
	pcs_deformation = 1;
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::PCSMoveNOD
   Task:
   Programming:
   08/2004 MB/OK Implementation
   last modified:
 **************************************************************************/
void CRFProcess::PCSMoveNOD(void)
{
	switch (this->type)
	{
		case 1:
			MSHMoveNODUcFlow(this);
			break;
		default:
			DisplayMsgLn("PCSMoveNOD: no valid process");
			abort();
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   10/2004 OK 2nd version
**************************************************************************/
std::string PCSProblemType()
{
	std::string pcs_problem_type;
	size_t no_processes(pcs_vector.size());

	for (size_t i = 0; i < no_processes; i++)
	{
		switch (pcs_vector[i]->getProcessType())
		{
			case FiniteElement::LIQUID_FLOW:
				pcs_problem_type = "LIQUID_FLOW";
				break;
			case FiniteElement::OVERLAND_FLOW:
				pcs_problem_type = "OVERLAND_FLOW";
				break;
			case FiniteElement::GROUNDWATER_FLOW:
				pcs_problem_type = "GROUNDWATER_FLOW";
				break;
			case FiniteElement::TWO_PHASE_FLOW:
				pcs_problem_type = "TWO_PHASE_FLOW";
				break;
			case FiniteElement::RICHARDS_FLOW: // MX test 04.2005
				pcs_problem_type = "RICHARDS_FLOW";
				break;
			case FiniteElement::DEFORMATION:
				if (pcs_problem_type.empty())
					pcs_problem_type = "DEFORMATION";
				else
					pcs_problem_type += "+DEFORMATION";
				break;
			case FiniteElement::DEFORMATION_FLOW:
				if (pcs_problem_type.empty())
					pcs_problem_type = "DEFORMATION";
				else
					pcs_problem_type += "+DEFORMATION";
				break;
			case FiniteElement::HEAT_TRANSPORT:
				if (pcs_problem_type.empty())
					pcs_problem_type = "HEAT_TRANSPORT";
				else
					pcs_problem_type += "+HEAT_TRANSPORT";
				break;
			case FiniteElement::MASS_TRANSPORT:
				if (pcs_problem_type.empty())
					pcs_problem_type = "MASS_TRANSPORT";
				else
					pcs_problem_type += "+MASS_TRANSPORT";
				break;
			case FiniteElement::FLUID_MOMENTUM:
				if (pcs_problem_type.empty())
					pcs_problem_type = "FLUID_MOMENTUM";
				else
					pcs_problem_type += "+FLUID_MOMENTUM";
				break;
			case FiniteElement::RANDOM_WALK:
				if (pcs_problem_type.empty())
					pcs_problem_type = "RANDOM_WALK";
				else
					pcs_problem_type += "+RANDOM_WALK";
				break;
			default:
				pcs_problem_type = "";
		}
	}

	//	//----------------------------------------------------------------------
	//	CRFProcess* m_pcs = NULL;
	//	// H process
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[0]) {
	//		case 'L':
	//			pcs_problem_type = "LIQUID_FLOW";
	//			break;
	//			// case 'U':
	//			//  pcs_problem_type = "UNCONFINED_FLOW";
	//			//  break;
	//		case 'O':
	//			pcs_problem_type = "OVERLAND_FLOW";
	//			break;
	//		case 'G':
	//			pcs_problem_type = "GROUNDWATER_FLOW";
	//			break;
	//		case 'T':
	//			pcs_problem_type = "TWO_PHASE_FLOW";
	//			break;
	//		case 'C':
	//			pcs_problem_type = "COMPONENTAL_FLOW";
	//			break;
	//		case 'R': //MX test 04.2005
	//			pcs_problem_type = "RICHARDS_FLOW";
	//			break;
	//		}
	//	}
	//	//----------------------------------------------------------------------
	//	// M process
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[0]) {
	//		case 'D':
	//			if (pcs_problem_type.empty())
	//				pcs_problem_type = "DEFORMATION";
	//			else
	//				pcs_problem_type += "+DEFORMATION";
	//			break;
	//		}
	//	}
	//	//----------------------------------------------------------------------
	//	// T process
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[0]) {
	//		case 'H':
	//			if (pcs_problem_type.empty())
	//				pcs_problem_type = "HEAT_TRANSPORT";
	//			else
	//				pcs_problem_type += "+HEAT_TRANSPORT";
	//			break;
	//		}
	//	}
	//	//----------------------------------------------------------------------
	//	// CB process
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[0]) {
	//		case 'M':
	//			if (pcs_problem_type.empty())
	//				pcs_problem_type = "MASS_TRANSPORT";
	//			else
	//				pcs_problem_type += "+MASS_TRANSPORT";
	//			break;
	//		}
	//	}
	//	//----------------------------------------------------------------------
	//	//----------------------------------------------------------------------
	//	// FM process
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[0]) {
	//		case 'F':
	//			if (pcs_problem_type.empty())
	//				pcs_problem_type = "FLUID_MOMENTUM";
	//			else
	//				pcs_problem_type += "+FLUID_MOMENTUM";
	//			break;
	//		}
	//	}
	//	//----------------------------------------------------------------------
	//	for (i = 0; i < no_processes; i++) {
	//		m_pcs = pcs_vector[i];
	//		switch (m_pcs->pcs_type_name[7]) { // _pcs_type_name[7] should be 'W' because 'R' is reserved for Richard
	//Flow.
	//		case 'W':
	//			if (pcs_problem_type.empty())
	//				pcs_problem_type = "RANDOM_WALK";
	//			else
	//				pcs_problem_type += "+RANDOM_WALK";
	//			break;
	//		}
	//	}

	return pcs_problem_type;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2004 OK Implementation
**************************************************************************/
// void CRFProcess::CalcELEMassFluxes(void)
//{
/*OK411
   int i;
   double e_value = -1.0;
   string e_value_name;
   double geo_factor, density;
   double velocity = 0.0;
   int e_idx;
   int phase = 0;
   long e;
   CMediumProperties* m_mmp = NULL;
   CFluidProperties* m_mfp = NULL;
   m_mfp = mfp_vector[phase]; //OK ToDo
   //======================================================================
   for(e=0;e<ElementListLength;e++){
   m_mmp = mmp_vector[ElGetElementGroupNumber(e)];
   geo_factor = m_mmp->geo_area;
   density = m_mfp->Density();
   for(i=0;i<pcs_number_of_evals;i++){
   e_value_name = pcs_eval_data[i].name;
   e_idx = PCSGetELEValueIndex(pcs_eval_data[i].name);
   if(e_value_name.find("MASS_FLUX1_X")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_X"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Y")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Y"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Z")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Z"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   }
   }
   //======================================================================
 */
//}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 OK Implementation
   02/2006 YD Dual Richards
   02/2007 WW General function for all unsaturated flow
**************************************************************************/
void CRFProcess::CalcSecondaryVariables(bool initial)
{
	//  char pcsT;
	//  pcsT = _pcs_type_name[0];
	//  if(type==1212) pcsT = 'V'; //WW
	//  switch(pcsT){
	//    case 'L':
	//      break;
	//    case 'U':
	//      break;
	//    case 'G':
	//      break;
	//    case 'T':
	//      break;
	//    case 'C':
	//      break;
	//    case 'R': // Richards flow
	//	  if(_pcs_type_name[1] == 'I')	// PCH To make a distinction with RANDOM WALK.
	//        CalcSecondaryVariablesUnsaturatedFlow(initial); // WW
	//      break;
	//    case 'D':
	//      break;
	//    case 'V':
	//      CalcSecondaryVariablesUnsaturatedFlow(initial); //WW
	//      break;
	//	case 'P':
	//      CalcSecondaryVariablesPSGLOBAL(); //WW
	//      break;
	//  }

	switch (getProcessType())
	{
		case FiniteElement::TNEQ:
			CalcSecondaryVariablesTNEQ(); // HS
			break;
		case FiniteElement::TES:
			CalcSecondaryVariablesTES(); // HS
			break;
		case FiniteElement::LIQUID_FLOW:
			if (aktueller_zeitschritt > 0)
				CalcSecondaryVariablesLiquidFlow();
			break;
		case FiniteElement::GROUNDWATER_FLOW:
			break;
		case FiniteElement::TWO_PHASE_FLOW:
			break;
		case FiniteElement::RICHARDS_FLOW: // Richards flow
		case FiniteElement::MULTI_PHASE_FLOW:
		case FiniteElement::DEFORMATION_H2: // H2M
			CalcSecondaryVariablesUnsaturatedFlow(initial);
			break;
		case FiniteElement::PS_GLOBAL:
			CalcSecondaryVariablesPSGLOBAL(); // WW
			break;
		default:
			break;
	}
}

//////////////////////////////////////////////////////////////////////////
// ReMove site
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Function: GetCompNamehelp
   Task: Namepatch, until primary function names are finally sorted out //SB:todo
   Programming:	08/2003 SB Implementation
   last modified:
   superseded by GetPFNamebyCPName() but left here, as not all files are already in the new concept
 **************************************************************************/
/* SB: namepatch
   Repariert kurzfristig die Ausgabe
   input: datafield_n[j].name
   wenn dar Name "CONCENTRATIONx" ist, wird er durch den enstprechenden Komponentennamen ersetzt, sonst bleibts */

char* GetCompNamehelp(char* inname)
{
	int comp; // WW, phase;
	char *outname, help[MAX_ZEILE];
	CRFProcess* m_pcs = NULL;
	outname = inname;
	// WW phase = 0;
	for (comp = 0; comp < GetRFProcessNumComponents(); comp++)
	{
		sprintf(help, "%s%d", "CONCENTRATION", comp + 1);
		/*  help has to be a part of inname (strstr) and also have the same length (strcmp) */
		if (strstr(inname, help) && (strcmp(inname, help) == 0))
		{
			m_pcs = m_pcs->GetProcessByFunctionName(help);
			if (m_pcs == NULL)
				break;
			//		outname = GetTracerCompName(phase,m_pcs->GetProcessComponentNumber()-1);
			//		outname = cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
			outname = (char*)cp_vec[m_pcs->GetProcessComponentNumber() - 1]->compname.data();
			return outname;
		}
	}
	return outname;
} // SB:namepatch

/*************************************************************************
   ROCKFLOW - Function: GetCPNamebyPFName(string )
   Task: Replaces CP Name by Primary function name for output input
   Input:	component property name
   Output: primary function name
   Programming:	10/2004 SB Implementation
   10/2010 TF changed access to process type
 **************************************************************************/
string GetPFNamebyCPName(string inname)
{
	int i, j; // WW, k;
	int pcs_vector_size = (int)pcs_vector.size();
	string outname;
	char help[MAX_ZEILE];
	CRFProcess* m_pcs = NULL;
	outname = "dummy";
	if (pcs_vector_size > 0)
		for (i = 0; i < pcs_vector_size; i++)
		{
			m_pcs = pcs_vector[i];
			//	if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT") == 0){ // if this is mass transport // TF
			if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				j = m_pcs->GetProcessComponentNumber();
				// WW k = cp_vec[j]->transport_phase;
				outname = cp_vec[m_pcs->GetProcessComponentNumber()]->compname;
				if (outname == inname) // right process found
				{
					sprintf(help, "%s%d", "CONCENTRATION", j);
					outname = help;
					return outname;
				}
			}
		}
	// for(i=0;comp<GetRFProcessNumComponents();i++) {
	//	sprintf(help,"%s%d","CONCENTRATION",i);
	/*  help has to be a part of inname (strstr) and also have the same length (strcmp) */
	//	if(strstr(inname, help) && (strcmp(inname,help) == 0)){
	//		m_pcs = m_pcs->GetProcessByFunctionName(help);
	//		if(m_pcs == NULL) break;
	//		outname = cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
	//		outname = (char *) cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname.data();
	//		if(outname.compare(inname) == 0)
	//			return outname;
	//	};
	// };
	// Inname is not from a mass transport process, therefore return inname
	return inname;
} // SB:namepatch

//========================================================================
// OK former model functions
int GetRFControlGridAdapt(void)
{
	// OK  return (get_rfcp_adaptive_mesh_refinement_flag(rfcp));
	if (show_onces_adp)
		cout << "GetRFControlGridAdapt - to be removed"
		     << "\n";
	show_onces_adp = false;
	return 0;
}

int GetRFControlModel(void)
{
	if (show_onces_mod)
		cout << "GetRFControlModel - to be removed"
		     << "\n";
	show_onces_mod = false;
	return -1;
}

int GetRFProcessChemicalModel(void)
{
	cout << "GetRFProcessChemicalModel - to be removed"
	     << "\n";
	return 0;
}

int GetRFProcessFlowModel(void)
{
	if (show_onces_mod_flow)
		cout << "GetRFProcessFlowModel - to be removed"
		     << "\n";
	show_onces_mod_flow = false;
	return 0;
}

int GetRFProcessHeatReactModel(void)
{
	cout << "GetRFProcessHeatReactModel - to be removed"
	     << "\n";
	return 0;
}

int GetRFProcessNumPhases(void)
{
	// DisplayMsgLn("GetRFProcessNumPhases - to be removed");
	int no_phases = (int)mfp_vector.size();
	return no_phases;
}

int GetRFProcessProcessing(char* rfpp_type)
{
	bool pcs_flow = false;
	bool pcs_deform = false;
	CRFProcess* m_pcs = NULL;
	size_t no_processes = pcs_vector.size();
	for (size_t i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos)
		if (isDeformationProcess(m_pcs->getProcessType()))
			pcs_deform = true;
		//		if (m_pcs->_pcs_type_name.find("FLOW") != string::npos)
		if (isFlowProcess(m_pcs->getProcessType()))
			pcs_flow = true;
	}

	if (strcmp(rfpp_type, "SD") == 0)
	{
		if (pcs_flow && pcs_deform)
			return 1;
	}
	else
		cout << "GetRFProcessProcessing - to be removed"
		     << "\n";
	return 0;
}

int GetRFProcessProcessingAndActivation(const char*)
{
	cout << "GetRFProcessProcessingAndActivation - to be removed"
	     << "\n";
	return 0;
}

long GetRFProcessNumComponents(void)
{
	// DisplayMsgLn("GetRFProcessNumComponents - to be removed");
	int no_components = (int)cp_vec.size();
	return no_components;
}

int GetRFControlModex(void)
{
	cout << "GetRFControlModex - to be removed"
	     << "\n";
	return 0;
}

int GetRFProcessDensityFlow(void)
{
	if (show_onces_density)
		cout << "GetRFProcessDensityFlow - to be removed"
		     << "\n";
	show_onces_density = false;
	return 0;
}

int GetRFProcessNumContinua(void)
{
	cout << "GetRFProcessNumContinua - to be removed"
	     << "\n";
	return 0;
}

int GetRFProcessNumElectricFields(void)
{
	cout << "GetRFProcessNumElectricFields - to be removed"
	     << "\n";
	return 0;
}

int GetRFProcessNumTemperatures(void)
{
	cout << "GetRFProcessNumTemperatures - to be removed"
	     << "\n";
	return -1;
}

int GetRFProcessSimulation(void)
{
	cout << "GetRFProcessSimulation - to be removed"
	     << "\n";
	return -1;
}

/**************************************************************************
   ROCKFLOW - Funktion: ModelsAddNodeValInfoStructure

   Aufgabe:
   Fuellt die Knotendaten-Infostruktur mit den zugehoerigen Modelldaten.
   Wird vom Modell der Reihe nach fuer jede Knotengroesse aufgerufen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:char *name         :Name der Knotengroesse fuer Ergebnisdatei
   E:char *einheit      :Name der phys. Einheit fuer Ergebnisdatei
   E:int speichern      :Werte sollen gespeichert werden (0/1)
   E:int laden          :Werte sollen geladen werden falls vorhanden (0/1)
   E:int restart        :Werte sollen bei Restart geladen werden (0/1)
   E:int adapt_interpol :Werte sollen beim verfeinern auf Kinder interpoliert (0/1)
   E:double vorgabe     :Vorgabe falls keine Restartdaten oder Anfangsbedingungen vorhanden sind

   Ergebnis:
   Knotenindex der gerade vergeben wurde

   Programmaenderungen:
   09/2000   CT    Erste Version

**************************************************************************/
int ModelsAddNodeValInfoStructure(char* name, char* einheit, int speichern, int laden, int restart, int adapt_interpol,
                                  double vorgabe)
{
	anz_nval++;
	nval_data = (NvalInfo*)Realloc(nval_data, anz_nval * sizeof(NvalInfo));

	nval_data[anz_nval - 1].name = NULL;
	nval_data[anz_nval - 1].einheit = NULL;

	if (name)
	{
		nval_data[anz_nval - 1].name = (char*)Malloc(((int)strlen(name) + 1) * sizeof(char));
		strcpy(nval_data[anz_nval - 1].name, name);
	}
	if (einheit)
	{
		nval_data[anz_nval - 1].einheit = (char*)Malloc(((int)strlen(einheit) + 1) * sizeof(char));
		strcpy(nval_data[anz_nval - 1].einheit, einheit);
	}

	nval_data[anz_nval - 1].speichern = speichern;
	nval_data[anz_nval - 1].laden = laden;
	nval_data[anz_nval - 1].restart = restart;
	nval_data[anz_nval - 1].adapt_interpol = adapt_interpol;
	nval_data[anz_nval - 1].vorgabe = vorgabe;

	return anz_nval - 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: ModelsAddElementValInfoStructure

   Aufgabe:
   Fuellt die Elementdaten-Infostruktur mit den zugehoerigen Modelldaten.
   Wird vom Modell der Reihe nach fuer jede Elementgroesse aufgerufen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:char *name         :Name der Elementgroesse fuer Ergebnisdatei
   E:char *einheit      :Name der phys. Einheit fuer Ergebnisdatei
   E:int speichern      :Werte sollen gespeichert werden (0/1)
   E:int laden          :Werte sollen geladen werden falls vorhanden (0/1)
   E:int restart        :Werte sollen bei Restart geladen werden (0/1)
   E:int adapt_interpol :Werte sollen beim verfeinern auf Kinder interpoliert (0/1)
   E:double vorgabe     :Vorgabe falls keine Restartdaten oder Anfangsbedingungen vorhanden sind

   Ergebnis:
   Elementindex der gerade vergeben wurde

   Programmaenderungen:
   09/2000   CT    Erste Version

**************************************************************************/
int ModelsAddElementValInfoStructure(char* name, char* einheit, int speichern, int laden, int restart,
                                     int adapt_interpol, double vorgabe)
{
	anz_eval++;
	eval_data = (EvalInfo*)Realloc(eval_data, anz_eval * sizeof(EvalInfo));

	eval_data[anz_eval - 1].name = NULL;
	eval_data[anz_eval - 1].einheit = NULL;

	if (name)
	{
		eval_data[anz_eval - 1].name = (char*)Malloc(((int)strlen(name) + 1) * sizeof(char));
		strcpy(eval_data[anz_eval - 1].name, name);
	}
	if (einheit)
	{
		eval_data[anz_eval - 1].einheit = (char*)Malloc(((int)strlen(einheit) + 1) * sizeof(char));
		strcpy(eval_data[anz_eval - 1].einheit, einheit);
	}

	eval_data[anz_eval - 1].speichern = speichern;
	eval_data[anz_eval - 1].laden = laden;
	eval_data[anz_eval - 1].restart = restart;
	eval_data[anz_eval - 1].adapt_interpol = adapt_interpol;
	eval_data[anz_eval - 1].vorgabe = vorgabe;

	return anz_eval - 1;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void PCSDelete()
{
	for (int i = 0; i < (int)pcs_vector.size(); i++)
		delete pcs_vector[i];
	pcs_vector.clear();
	pcs_no_components = 0;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modified:
**************************************************************************/
void CRFProcess::SetNodeValue(long n, int nidx, double value)
{
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::SetNodeValue() "
		     << "\n";
		abort();
	}
#endif
	// WW 11.12.2012 	nod_val_vector[n][nidx] = value;
	nod_val_vector[nidx][n] = value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
void CRFProcess::SetElementValue(long n, int nidx, double value)
{
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::SetElementValue() "
		     << "\n";
		abort();
	}
#endif
	ele_val_vector[n][nidx] = value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modified:
**************************************************************************/
double CRFProcess::GetNodeValue(size_t n, int nidx)
{
	double value;
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::GetNodeValue() "
		     << "\n";
		abort();
	}
#endif
	// WW 11.12.2012		value = nod_val_vector[n][nidx];
	value = nod_val_vector[nidx][n];
	return value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
double CRFProcess::GetElementValue(size_t n, int nidx)
{
	double value;
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in CRFProcess::GetElementValue() "
		     << "\n";
		abort();
	}
#endif
	value = ele_val_vector[n][nidx];
	return value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   01/2006 OK Test
   03/2012 JT: Checks are unnecessary. Allow reverse order.
**************************************************************************/
int CRFProcess::GetNodeValueIndex(const std::string& var_name, bool reverse_order) // JT: Allow reverse order
{
	if (!reverse_order)
	{
		for (size_t i = 0; i < nod_val_name_vector.size(); i++)
		{
			if (nod_val_name_vector[i].compare(var_name) == 0)
				return i;
		}
	}
	else
	{
		int nvals = ((int)nod_val_name_vector.size()) - 1;
		for (int j = nvals; j > -1; j--)
		{
			if (nod_val_name_vector[j].compare(var_name) == 0)
				return j;
		}
	}
	//
	return -2;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
int CRFProcess::GetElementValueIndex(const string& var_name, bool reverse_order)
{
	if (!reverse_order)
	{
		for (size_t i = 0; i < ele_val_name_vector.size(); i++)
		{
			if (ele_val_name_vector[i].compare(var_name) == 0)
				return i;
		}
	}
	else
	{
		int nvals = ((int)ele_val_name_vector.size()) - 1;
		for (int j = nvals; j > -1; j--)
		{
			if (ele_val_name_vector[j].compare(var_name) == 0)
				return j;
		}
	}
	//
	return -2;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   05/2005 OK pcs_pv_name,
   12/2005 OK RESTART
   07/2006 OK/MX MSH
**************************************************************************/
void CRFProcess::SetIC()
{
	CInitialCondition* m_ic = NULL;
	// HS, for MASS_TRANSPORT PCS,
	// it is not necessary to use PrimaryVarible as second check.
	// nidx will give the proper IC pointer.
	if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);

			// PrimaryVariable pv_i (convertPrimaryVariable(pcs_primary_function_name[i]));
			for (size_t j = 0; j < ic_vector.size(); j++)
			{
				m_ic = ic_vector[j];
				m_ic->m_msh = m_msh; // OK/MX

				if (m_ic->getProcess() == this)
				{
					m_ic->Set(nidx);
					m_ic->Set(nidx + 1);
				}
			}
		}
	else // otherwise PrimaryVariable check is still performed.

		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);
			FiniteElement::PrimaryVariable pv_i(FiniteElement::convertPrimaryVariable(pcs_primary_function_name[i]));
			for (size_t j = 0; j < ic_vector.size(); j++)
			{
				m_ic = ic_vector[j];
				m_ic->m_msh = m_msh; // OK/MX

				if (m_ic->getProcessType() != this->getProcessType())
					continue;

				m_ic->setProcess(this);
				if (m_ic->getProcessPrimaryVariable() == pv_i)
				{
					m_ic->Set(nidx);
					m_ic->Set(nidx + 1);
				} // end of if
			} // end of for j
		} // end of for i

	// end of if-else
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modified:
**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void CRFProcess::SetNODValues()
{
	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		//    SetNODValue(i,GetNODValueIndex(_pcs_type_name),eqs->x[i]);
		//    SetNODValue(i,GetNODValueIndex(_pcs_type_name)+1,eqs->x[i]);
		// WW
		SetNodeValue(m_msh->Eqs2Global_NodeIndex[i], 0, eqs->x[i]);
		// WW
		SetNodeValue(m_msh->Eqs2Global_NodeIndex[i], 1, eqs->x[i]);
	}
}
#endif

/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt den Fehler bei Iterationen
   new_iteration     : Vektor des neuen Iterationsschritts
   old_iteration_ndx : Knotenindex fuer Werte des alten Iterationsschritts
   reference_ndx     : Knotenindex fuer Werte des alten Zeitschritts (als Referenz)
   method            : Methode der Fehlerermittlung
   Programing:
   01/2005 OK NUM implementation
   05/2005 OK MSH
   08/2005 WW Re-implememtation based on NUMCalcIterationError
   01/2007 WW For DOF>1
   11/2007 WW Changes for the new classes of sparse and linear solver
   3/2012  JT Clean, add newton, add CPL vs. NLS, go to enum system
   last modification:
**************************************************************************/
#if !defined(USE_PETSC) // && !defined(other parallel libs)//02.3013. WW
double CRFProcess::CalcIterationNODError(FiniteElement::ErrorMethod method, bool nls_error, bool cpl_error)
{
	static long i, k, g_nnodes;
	static double error, error_g, val1, val2, value;
	int nidx1, ii;
	double* eqs_x = NULL; // 11.2007. WW
	int num_dof_errors = pcs_number_of_primary_nvals;
	double unknowns_norm = 0.0;
	double absolute_error[DOF_NUMBER_MAX];
#if defined(USE_PETSC) // || defined(other parallel libs)//02.2014. WW
	g_nnodes = m_msh->getNumNodesLocal();
	eqs_x = eqs_new->GetGlobalSolution();
#else
		g_nnodes = m_msh->GetNodesNumber(false);

#ifdef NEW_EQS
		eqs_x = eqs_new->x;
#else
		eqs_x = eqs->x;
#endif

#endif // if defined(USE_PETSC)

	switch (method)
	{
		//
		// --> ENORM:	|x1-x0|
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over entire solution vector (all primary variables) and checked against a single tolerance.
		//
		case FiniteElement::ENORM:
			if (m_num->nls_method > 0)
			{ // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							val1 = eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
					}
				}
			}
			else
			{ // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = GetNodeValue(i, nidx1)
						       - eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm;
			break;
		//
		// --> ERNORM:	|(x1-x0)/x0)|
		//     Norm of the solution vector delta divided by the solution vector (relative error).
		//     A single tolerance applied to all primary variables.
		//
		case FiniteElement::ERNORM:
			value = 0.0;
			if (m_num->nls_method > 0)
			{ // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//02.2018. WW
						val1 = eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
						val2 = GetNodeValue(i, nidx1);
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = eqs_x[i + ii * g_nnodes];
							val2 = GetNodeValue(k, nidx1);
#endif
						//
						unknowns_norm += val1 * val1;
						value += val2 * val2;
					}
				}
			}
			else
			{ // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = GetNodeValue(i, nidx1);
						val2 = val1 - eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = GetNodeValue(k, nidx1);
							val2 = val1 - eqs_x[i + ii * g_nnodes];
#endif
						//
						unknowns_norm += val2 * val2;
						value += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm / (sqrt(value) + DBL_EPSILON);
			break;
		//
		// --> EVNORM:	|x1-x0|
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over solution vector of each primary variable, checked againes a tolerence specific to each
		//     variable.
		//
		case FiniteElement::EVNORM:
			if (m_num->nls_method > 0)
			{ // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							val1 = eqs_x[i + ii * g_nnodes];
#endif
						error += val1 * val1;
					}
					unknowns_norm += error;
					absolute_error[ii] = sqrt(error);
				}
			}
			else
			{ // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = GetNodeValue(i, nidx1)
						       - eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
#endif
						error += val1 * val1;
					}
					unknowns_norm += error;
					absolute_error[ii] = sqrt(error);
				}
			}
			unknowns_norm = sqrt(unknowns_norm);
			break;
		//
		// --> BNORM: Get norm of solution vector, same as ENORM. RHS norm will be calculated later.
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over entire solution vector (all primary variables) and checked against a single tolerance.
		//
		case FiniteElement::BNORM:
			if (m_num->nls_method > 0)
			{ // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							val1 = eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
					}
				}
			}
			else
			{ // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
						val1 = GetNodeValue(i, nidx1)
						       - eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm;
			break;
		//
		// --> LMAX:	max(x1-x0)
		//     Local max error (across all elements) of solution vector delta (absolute error).
		//     Tolerance required for each primary variable.
		//
		case FiniteElement::LMAX:
			if (m_num->nls_method > 0)
			{ // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.3014. WW
						val1 = eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							val1 = eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
						val1 = fabs(val1);
						if (val1 > error)
							error = val1;
					}
					absolute_error[ii] = error;
				}
			}
			else
			{ // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
#if defined(USE_PETSC) // || defined(other parallel libs)//08.3014. WW
						val1 = GetNodeValue(i, nidx1)
						       - eqs_x[pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[i] + ii];
#else
							k = m_msh->Eqs2Global_NodeIndex[i];
							val1 = GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
						val1 = fabs(val1);
						if (val1 > error)
							error = val1;
					}
					absolute_error[ii] = error;
				}
			}
			unknowns_norm = sqrt(unknowns_norm);
			break;
		//
		default:
			ScreenMessage("ERROR: Invalid error method for Iteration or Coupling Node error.\n");
			return 0.0;
			//
			/*
			-----------------------------------------------------------------------------------------------
			ALTERNATIVE METHODS NOT YET IMPLEMENTED. MODIFY THEM AND ADD THEIR ENUM VALUES IF YOU WANT THEM.
			-----------------------------------------------------------------------------------------------
			// METHOD 4
			case 4:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = max_c = 0.0;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1)));
			           max_c = MMax(MMax(max_c, fabs(fabs(eqs_x[i+ii*g_nnodes]))),fabs(GetNodeValue(k, nidx1)));
			        }
			        pcs_absolute_error[ii] = error / (max_c + MKleinsteZahl);
			    }
			    break;
			//
			// METHOD 5
			case 5:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = max_c = 0.0;
			        min_c = 1.e99;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1)));
			           min_c = MMin(min_c, fabs(eqs_x[i+ii*g_nnodes]));
			           max_c = MMax(max_c, fabs(eqs_x[i+ii*g_nnodes]));
			        }
			        pcs_absolute_error[ii] = error / (max_c - min_c + MKleinsteZahl) ;
			    }
			    break;
			//
			// METHOD 6
			case 6:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = 0.0;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			        //
			        for (i = 0l; i < g_nnodes; i++) {
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -  GetNodeValue(k, nidx1))
			              / (fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1-1)) + MKleinsteZahl));
			        }
			        pcs_absolute_error[ii] = error;
			    }
			    break;
			//
			// METHOD 7
			case 7:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = change = max_c = 0.0;
			        min_c = 1.e99;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -  GetNodeValue(k, nidx1)));
			           change = MMax(change, fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k, nidx1-1)));
			        }
			        pcs_absolute_error[ii] = error / (change + MKleinsteZahl);
			    }
			    break;
			//
			*/
	}

#if defined(USE_PETSC) // || defined(other parallel libs)//08.2014. WW
	for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
		double val_temp = pcs_absolute_error[ii];
		pcs_absolute_error[ii] = 0.;
		MPI_Allreduce(&val_temp, &pcs_absolute_error[ii], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	//
	// Store the error (JT)
	// JT: now returning RELATIVE error. NECESSARY BECAUSE DOF MAY BE > 1 AND EACH DOF MAY HAVE DIFFERENT CHARACTER
	if (cpl_error)
	{ // Return coupling error
		error_g = 0.0;
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			cpl_absolute_error[ii] = absolute_error[ii];
			error = absolute_error[ii] / m_num->cpl_error_tolerance[ii];
			error_g = MMax(error_g, error); // Coupling error just stores the maximum
		}
		cpl_max_relative_error = error_g;
		cpl_num_dof_errors = num_dof_errors;
	}
	if (nls_error)
	{ // Return Non-Linear iteration error
		error_g = 0.0;
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			pcs_absolute_error[ii] = absolute_error[ii];
			pcs_relative_error[ii] = absolute_error[ii] / m_num->nls_error_tolerance[ii];
			error_g = MMax(error_g, pcs_relative_error[ii]);
		}
		nls_max_relative_error = error_g;
		pcs_num_dof_errors = num_dof_errors;
		pcs_unknowns_norm = unknowns_norm;
	}
	if (!nls_error && !cpl_error)
	{ // Then this routine called from somewhere else to get the error (i.e. time control). Store it in a temporary
		// vector for access.
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			temporary_absolute_error[ii] = absolute_error[ii];
		}
		temporary_num_dof_errors = num_dof_errors;
	}
	//
	return error_g; // Always returns the maximum relative error
}
#endif // #if !defined(USE_PETSC)  WW

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   04/2006 YD Add contiuum coupling OK???Why ere
   04/2007 WW Remove the spurious stuff
   08/2008 WW Time step size control (First)
   12/2008 WW Time step size control (Update)
   07/2011 WW Newton-Raphson method
   3/2012  JT Clean, correct error obtainment, modify Newton convergence criteria
**************************************************************************/

double CRFProcess::ExecuteNonLinear(int loop_process_number, bool print_pcs)
{
	double nonlinear_iteration_error = 0.0;
	double nl_theta, damping, norm_x0, norm_b0, norm_x, norm_b;
	double error_x1, error_x2, error_b1, error_b2 = 0, error, last_error, percent_difference;
	// double* eqs_x = NULL;     //
	bool converged, diverged;
	int ii, nidx1, num_fail = 0;
	size_t j, g_nnodes;

	string delim = " ";
	damping = 1.0;
	norm_x0 = norm_b0 = norm_x = norm_b = 0.;
	error = 1.;
	error_x2 = DBL_MAX;
	nl_theta = 1.0 - m_num->nls_relaxation; // JT
	if (nl_theta < DBL_EPSILON)
		nl_theta = 1.0;
	g_nnodes = m_msh->GetNodesNumber(false);

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	eqs_x = eqs_new->GetGlobalSolution();
#endif
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	double* eqs_b = NULL;
	int k;
#ifdef NEW_EQS
	eqs_x = eqs_new->x;
	eqs_b = eqs_new->b;
	configured_in_nonlinearloop = true;
// Also allocate temporary memory for linear solver. WW
//
#if defined(USE_MPI)
	CPARDomain* dom = dom_vector[myrank];
	dom->eqs->SetDOF(pcs_number_of_primary_nvals); //_new 02/2010 WW
	dom->ConfigEQS(m_num, pcs_number_of_primary_nvals * m_msh->GetNodesNumber(false));
#else
	eqs_new->SetDOF(pcs_number_of_primary_nvals); //_new 02/2010. WW
	eqs_new->ConfigNumerics(m_num);
#endif
//
#else // ifdef NEW_EQS
		eqs_x = eqs->x;
		eqs_b = eqs->b;
#endif
#endif
	//..................................................................
	// PI time step size control. 29.08.2008. WW
	if (Tim->GetPITimeStepCrtlType() > 0)
		CopyU_n();
	if (hasAnyProcessDeactivatedSubdomains)
		this->CheckMarkedElement(); // NW
	if (ExcavMaterialGroup > -1) // WX:07.2011 HM excavation
		this->CheckExcavedElement();
	Tim->last_dt_accepted = true; // JT2012

#if defined(USE_PETSC) || defined(USE_MPI) // || defined(other parallel libs)//01.3013. WW
	if (myrank == 0)
	{
#endif
		if (print_pcs)
		{ // JT: need check because of Regional Richards
			std::cout << "\n      ================================================"
			          << "\n";
			if (getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				std::cout << "    ->Process   " << loop_process_number << ": "
				          << convertProcessTypeToString(getProcessType()) << "\n";
				std::cout << "    ->Component " << pcs_component_number << ": " << pcs_primary_function_name[0] << "\n";
			}
			else
			{
				std::cout << "    ->Process " << loop_process_number << ": "
				          << convertProcessTypeToString(getProcessType()) << "\n";
			}
			std::cout << "      ================================================"
			          << "\n";
		}
#if defined(USE_PETSC) || defined(USE_MPI) // || defined(other parallel libs)//01.3013. WW#ifdef USE_MPI
	}
#endif

	// ------------------------------------------------------------
	// NON-LINEAR ITERATIONS (OR SINGLE ITERATION IF LINEAR)
	// ------------------------------------------------------------
	diverged = false;
	converged = false;
	accepted = true;
	last_error = 1.0;
	for (iter_nlin = 0; iter_nlin < m_num->nls_max_iterations; iter_nlin++)
	{
		cout << "    PCS non-linear iteration: " << iter_nlin << "/" << m_num->nls_max_iterations << '\n';
		nonlinear_iteration_error = Execute();
		//
		// ---------------------------------------------------
		// LINEAR SOLUTION
		// ---------------------------------------------------
		if (m_num->nls_method < 0)
		{
			PrintStandardIterationInformation(true);
			converged = true;
		}
		else
		{
			// ---------------------------------------------------
			// NON-LINEAR SOLUTION
			// ---------------------------------------------------
			//
			damping = nl_theta;
			switch (m_num->getNonLinearErrorMethod())
			{
				// For most error methods (also works for Newton)
				default:
					PrintStandardIterationInformation(true);
					//
					if (nonlinear_iteration_error <= 1.0)
					{
						converged = true;
					}
					else
					{ // Check for stagnation
						percent_difference = 100.0 * ((last_error - nonlinear_iteration_error) / last_error);
						if (iter_nlin > 0
						    && percent_difference
						           < 1.0) // less than 1% difference (or an error increase) from previous error
							num_fail++;
						else
							num_fail = 0;
						//
						if (num_fail > 1)
							diverged = true; // require 2 consecutive failures
						last_error = nonlinear_iteration_error;
					}
					break;

				// For (OGS) classic Newton error control
				case FiniteElement::BNORM:
					PrintStandardIterationInformation(false);
//
#if defined(USE_PETSC) // || defined(other parallel libs)//06.3012. WW
					norm_x = eqs_new->GetVecNormX();
					norm_b = eqs_new->GetVecNormRHS();
#else
					norm_x = pcs_unknowns_norm; // JT: this is already obtained in CalcIterationNodeError.
					norm_b = 0.0; // must calculate this
					for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
					{
						for (j = 0; j < g_nnodes; j++)
						{
							double val = eqs_b[j + ii * g_nnodes];
							norm_b += val * val;
						}
					}
					norm_b = sqrt(norm_b);
#endif
					//
					if (iter_nlin == 0)
					{
						norm_x0 = norm_x;
						norm_b0 = norm_b;
						error_x2 = error_b2 = DBL_MAX;
					}
					else
					{
						error_x1 = norm_x / norm_x0;
						error_b1 = norm_b / norm_b0;

						if (norm_x < m_num->nls_error_tolerance[0] && error_x1 > norm_x)
							error_x1 = norm_x;
						if (norm_b < m_num->nls_error_tolerance[0] && error_b1 > norm_b)
							error_b1 = norm_b;
						if (error_x1 / error_x2 > 0.1 || error_b1 / error_b2 > 0.1)
							damping *= 0.5; // take 1/2 of original theta
						//
						error = max(error_x1, error_b1);
						error_x2 = error_x1;
						error_b2 = error_b1;
						//
						// Check for divergence
						if (error > 10.0 && iter_nlin > 1)
						{
							diverged = true;
							if (Tim->GetPITimeStepCrtlType() > 0)
							{ // if PI automatic time control
								accepted = false;
								PI_TimeStepSize();
								return error;
							}
							else
							{
								cout << "Attention: Newton-Raphson step is diverged. Programme halt!"
								     << "\n";
#if defined(USE_PETSC) // || defined(other parallel libs)//08.3012. WW
								PetscFinalize();
#endif
#if defined(USE_MPI)
								MPI_Finalize();
#endif
								exit(1);
							}
						}
						//
						// Check convergence
						if (norm_x0 < m_num->nls_error_tolerance[0])
						{
							error = norm_x0;
							converged = true;
						}
						if (norm_b0 < 10 * m_num->nls_error_tolerance[0])
						{
							error = norm_b0;
							converged = true;
						}
						if (norm_b < 0.001 * norm_b0)
						{
							error = norm_b;
							converged = true;
						}
						if (error <= m_num->nls_error_tolerance[0])
						{
							converged = true;
						}
					}
// Newton information printout.
#if defined(USE_MPI) || defined(USE_PETSC)
					if (myrank == 0)
					{
#endif
						cout.width(10);
						cout.precision(3);
						cout.setf(ios::scientific);
						cout << "         NR-Error  |"
						     << "    RHS Norm|"
						     << "  Unknowns Norm|"
						     << " Damping\n";
						cout << "         " << setw(10) << error << "|  " << setw(9) << norm_b << "| ";
						cout << setw(14) << norm_x << "| " << setw(9) << damping << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
					}
#endif
					break;
			}
		}

		// CHECK FOR TIME STEP FAILURE
		// ---------------------------------------------------

		if (m_num->nls_max_iterations > 1 && !converged && iter_nlin + 1 == m_num->nls_max_iterations)
			accepted = false;

		if (!accepted || Tim->isDynamicTimeFailureSuggested(this))
		{
			accepted = false;
			Tim->last_dt_accepted = false;
			break;
		}

		// FOR NEWTON: COPY DAMPED CHANGES TO NEW TIME
		// ---------------------------------------------------
		if (m_num->nls_method > 0)
		{
			if (converged)
				damping = 1.0; // Solution has converged. Take newest values.
			//
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
				for (j = 0; j < g_nnodes; j++)
				{
					SetNodeValue(
					    j, nidx1,
					    GetNodeValue(j, nidx1)
					        + damping * eqs_x[m_msh->Eqs2Global_NodeIndex[j] * pcs_number_of_primary_nvals + ii]);
				}
#else
				const long ish = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
					k = m_msh->Eqs2Global_NodeIndex[j];
					double val = GetNodeValue(k, nidx1) + damping * eqs_x[j + ish];
					SetNodeValue(k, nidx1, val);
				}
#endif
			}
		}

		// OTHER CONFIGURATIONS AT THE END OF THIS NON-LINEAR ITERATION
		// ---------------------------------------------------
		if (mobile_nodes_flag == 1)
		{
			PCSMoveNOD();
		}
		/* JT: I don't know if this time control method is used anymore. But it relies on a single error
		       produced from CalcIterationNodeError(), but this now depends on the type of error to use.
		       Therefore, I simply provide the error of the first dof, and not depending on the error type. If
		       this time step is still used, someone will need to find another way to calculate the error it uses.
		*/
		Tim->repeat = true;
		if (converged)
		{
			Tim->repeat = false;
			Tim->nonlinear_iteration_error = pcs_absolute_error[0];
		}

		// BREAK CRITERIA
		if (converged || diverged)
		{
			break;
		}
	}
	iter_nlin_max = std::max(iter_nlin_max, iter_nlin);
	// ------------------------------------------------------------
	// NON-LINEAR ITERATIONS COMPLETE
	// ------------------------------------------------------------
	// PI time step size control. 27.08.2008. WW
	if (Tim->GetPITimeStepCrtlType() > 0)
	{
		PI_TimeStepSize(); // might also set accepted to false here.
	}
	//
	if (m_num->nls_max_iterations > 1) // only for non-linear iterations
	{
		if (diverged)
		{
			if (accepted) // only increment if not fixed by a failed time step.
				num_diverged++;
			std::cout << "\nNon-linear iteration stabilized.";
		}
		else if (!converged)
		{
			if (accepted) // only increment if not fixed by a failed time step.
				num_notsatisfied++;
			if (Tim->GetPITimeStepCrtlType() < 1) // PI has the intrinsic property of doing this. So don't print it.
				std::cout << "\nMax number of non-linear iterations reached.";
		}
	}
	//
	// Calculate secondary variables
	if (accepted)
	{
		CalcSecondaryVariables();
	}
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
// Release temporary memory of linear solver. WW
#ifdef NEW_EQS // WW
#if defined(USE_MPI)
	dom->eqs->Clean();
#else
	eqs_new->Clean(); // Release buffer momery WW
#endif
	configured_in_nonlinearloop = false;
#endif
#endif
	return nonlinear_iteration_error;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   3/2012  JT
**************************************************************************/
void CRFProcess::PrintStandardIterationInformation(bool write_std_errors)
{
	int ii;
	//
	// LINEAR SOLUTION
	if (m_num->nls_method < 0)
	{
		std::cout << "      -->LINEAR solution complete. "
		          << "\n";
		if (write_std_errors)
		{
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				std::cout << "         PCS error DOF[" << ii << "]: " << pcs_absolute_error[ii] << "\n";
			}
		}
		return;
	}
	//
	// NON-LINEAR METHODS
	if (m_num->nls_method == 0)
		std::cout << "      -->End of PICARD iteration: " << iter_nlin << "/" << m_num->nls_max_iterations << "\n";
	else
		std::cout << "      -->End of NEWTON-RAPHSON iteration: " << iter_nlin << "/" << m_num->nls_max_iterations
		          << "\n";
	//
	// Errors
	// --------------------------------------------------
	if (write_std_errors)
	{
		if (pcs_num_dof_errors == 1)
		{
			std::cout << "         PCS error: " << pcs_absolute_error[0] << "\n";
		}
		else
		{
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				std::cout << "         PCS error DOF[" << ii << "]: " << pcs_absolute_error[ii] << "\n";
			}
		}
		std::cout << "         ->Euclidian norm of unknowns: " << pcs_unknowns_norm << "\n";
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_GaussValue()
{
	int k, NS;
	long i = 0;
	int idx[3];
	// const long LowOrderNodes= m_msh->GetNodesNumber(false);
	MeshLib::CElem* elem = NULL;

	//
	NS = m_msh->GetCoordinateFlag() / 10;
	idx[0] = GetNodeValueIndex("VELOCITY_X1");
	idx[1] = GetNodeValueIndex("VELOCITY_Y1");
	idx[2] = GetNodeValueIndex("VELOCITY_Z1");
	for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		for (k = 0; k < NS; k++)
			SetNodeValue(i, idx[k], 0.0);
	if (type == 1212 || type == 1313) // Multi-phase flow
	{
		idx[0] = GetNodeValueIndex("VELOCITY_X2");
		idx[1] = GetNodeValueIndex("VELOCITY_Y2");
		idx[2] = GetNodeValueIndex("VELOCITY_Z2");
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
			for (k = 0; k < NS; k++)
				SetNodeValue(i, idx[k], 0.0);
	}
	//
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			for (k = 0; k < NS; k++)
				fem->ExtropolateGauss(*elem, this, k);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:    Calculate the material values at Gauss points and extropolate them
         to element node
   Programing:
   04/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_MatValue()
{
	//	if (_pcs_type_name.find("FLOW") == string::npos)
	if (!isFlowProcess(this->getProcessType()))
		return;
	if (additioanl2ndvar_print < 0)
		return;

	//
	int NS = m_msh->GetCoordinateFlag() / 10;
	//
	if ((additioanl2ndvar_print > 0) && (additioanl2ndvar_print < 3))
	{
		int idx[3];
		idx[0] = GetNodeValueIndex("PERMEABILITY_X1");
		idx[1] = GetNodeValueIndex("PERMEABILITY_Y1");
		if (NS > 2)
			idx[2] = GetNodeValueIndex("PERMEABILITY_Z1");
		else
			idx[2] = 0;
		for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			for (int k = 0; k < NS; k++)
				SetNodeValue(i, idx[k], 0.0);
	}
	if (additioanl2ndvar_print > 1)
	{
		int idxp = GetNodeValueIndex("POROSITY");
		for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			SetNodeValue(i, idxp, 0.0);
	}
	//
	continuum = 0;
	if (continuum_vector.size() == 2)
		continuum = 1;

	MeshLib::CElem* elem = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			fem->CalcNodeMatParatemer(*elem);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
**************************************************************************/
void PCSDelete(const std::string& m_pcs_type_name)
{
	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(m_pcs_type_name));
	CRFProcess* m_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->_pcs_type_name.compare(m_pcs_type_name) == 0) { TF
		if (m_pcs->getProcessType() == pcs_type)
		{
			delete m_pcs;
			pcs_vector.erase(pcs_vector.begin() + i);
		}
	}
}

/**************************************************************************
   GeoSys - Function: Reallocation

   Aufgabe:
        Reallocte memory by new operator
   09/2005   WW    Erste Version

**************************************************************************/
template <class T>
T* resize(T* array, size_t old_size, size_t new_size)
{
	T* temp = new T[new_size];
	for (size_t i = 0; i < old_size; i++)
		temp[i] = array[i];
	for (size_t i = old_size; i < new_size; i++)
		temp[i] = 0;
	delete[] array;
	array = temp;
	return temp;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB Implementation
**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void CRFProcess::CalcFluxesForCoupling(void)
{
	int i, j;
	double flux;
	long n_index;
	long NodeIndex_GW;
	long NodeIndex_OLF;
	int no_nodes = m_msh->getNumberOfMeshLayers() + 1;
	long no_richards_problems = (long)(m_msh->ele_vector.size() / m_msh->getNumberOfMeshLayers());
	long IndexBottomNode;
	long IndexTopNode;
	int NoOfGWNodes = 0;
	double AverageZ_GW = 0.0;
	double AverageZ_OLF = 0.0;
	double AverageH_GW = 0.0;
	double AverageH_OLF = 0.0;
	double dh;
	int idxFLUX;
	int idxHead_GW;
	int idxHead_OLF;
	MeshLib::CElem* m_ele_GW = NULL;
	MeshLib::CElem* m_ele_OLF = NULL;

	// Get processes
	CRFProcess* m_pcs_GW(PCSGet(FiniteElement::GROUNDWATER_FLOW));
	if (!m_pcs_GW) // OK
	{
		cout << "Fatal error: no GROUNDWATER_FLOW process"
		     << "\n";
		return;
	}
	CRFProcess* m_pcs_OLF(PCSGet(FiniteElement::OVERLAND_FLOW));
	if (!m_pcs_OLF) // OK
	{
		cout << "Fatal error: no OVERLAND_FLOW process"
		     << "\n";
		return;
	}

	// Get meshes
	CFEMesh* m_msh_GW = m_pcs_GW->m_msh;
	CFEMesh* m_msh_OLF = m_pcs_OLF->m_msh;

	// Get indeces
	idxHead_GW = m_pcs_GW->GetNodeValueIndex("HEAD") + 1;
	idxHead_OLF = m_pcs_OLF->GetNodeValueIndex("HEAD") + 1;
	idxFLUX = GetNodeValueIndex("FLUX") + 1;

	for (i = 0; i < no_richards_problems; i++)
	{
		IndexBottomNode = ((i + 1) * no_nodes) - 1;

		// ToDo safe somewhere else so that this has to be done only once
		//-----------------------------------------------------------------
		// Get Nearest GW and OLF Element
		GEOLIB::Point pnt(m_msh->nod_vector[IndexBottomNode]->getData());

		long EleNumber = m_msh_GW->GetNearestELEOnPNT(&pnt);

		// GW and OLF use the same Numbering !!!
		m_ele_GW = m_msh_GW->ele_vector[EleNumber];
		m_ele_OLF = m_msh_OLF->ele_vector[EleNumber];

		//-----------------------------------------------------------------
		// Get Average values for element //ToDo encapsulate //WW: CElement::elemnt_averag??e
		NoOfGWNodes = m_ele_OLF->GetNodesNumber(m_msh_GW->getOrder());
		for (j = 0; j < NoOfGWNodes; j++)
		{
			NodeIndex_GW = m_ele_GW->GetNodeIndex(j);
			NodeIndex_OLF = m_ele_OLF->GetNodeIndex(j);

			AverageZ_GW += m_pcs_GW->GetNodeValue(NodeIndex_GW, idxHead_GW);
			AverageZ_OLF += m_msh_OLF->nod_vector[NodeIndex_OLF]->getData()[2];
		}
		AverageZ_GW = AverageZ_GW / NoOfGWNodes;
		AverageZ_OLF = AverageZ_OLF / NoOfGWNodes;

		//-----------------------------------------------------------------
		// UsatZone exists -> Flux from this
		if (AverageZ_GW < AverageZ_OLF)
		{
			n_index = m_msh->Eqs2Global_NodeIndex[IndexBottomNode];
			if (m_msh->nod_vector[IndexBottomNode]->GetMark())
			{
				flux = eqs->b[IndexBottomNode];
				// FLUXES IN NEW VERSION WITH VELOCITIES !!!!!
				// WAIT FOR SEBASTIANS MASS TRANSPORT IN USAT ZONE !!!!!
				// TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
				flux = 0.00001;
				// flux = 1;
				SetNodeValue(n_index, idxFLUX, flux);
			}
		}

		//-----------------------------------------------------------------
		// No UsatZone -> Calculate Flux from leakage terms
		if (AverageZ_GW >= AverageZ_OLF)
		{
			// SetRichardsNodesToFullySaturated??
			IndexTopNode = i * no_nodes; // Top Node of Richard Column -> Flux for OLF
			// Bottom Node of Richards Column -> Flux for GW
			IndexBottomNode = ((i + 1) * no_nodes) - 1;
			//-----------------------------------------------------------------
			// Get Average values for element //ToDo encapsulate
			for (j = 0; j < NoOfGWNodes; j++)
			{
				NodeIndex_GW = m_ele_GW->GetNodeIndex(j);
				NodeIndex_OLF = m_ele_OLF->GetNodeIndex(j);
				AverageH_GW += m_pcs_GW->GetNodeValue(NodeIndex_GW, idxHead_GW);
				AverageH_OLF += m_pcs_OLF->GetNodeValue(NodeIndex_OLF, idxHead_OLF);
			}
			AverageH_GW = AverageH_GW / NoOfGWNodes;
			AverageH_OLF = AverageH_OLF / NoOfGWNodes;
			// Calculate the vertical leakage
			dh = AverageH_GW - AverageH_OLF;
			// get kf fully saturated of uppermost element ?
			// or user defined value: entry resistance / leakage factor ?
			// flux = dh * 0.001;
			flux = dh * 1.;

			// 1. Add reacharge value to GW flow -> Add to flux off IndexBottomNode
			// Achtung nur zum Testen Source für GW flow durchgehend !!!!!!
			// SetNodeValue(IndexBottomNode, idxFLUX, -flux);  //H_OLF  > H_GW -> + flux_GW
			SetNodeValue(IndexBottomNode, idxFLUX, 0.00001);

			// 2. Add reacharge value to OLF -> Add to flux off IndexTopNode
			// H_OLF  > H_GW -> - flux_OLF
			SetNodeValue(IndexTopNode, idxFLUX, flux);
			// 3. Set flag to set reacharge to Usat to zero ???
		}
	}
}
#endif
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB implementation
**************************************************************************/
void CRFProcess::CopyCouplingNODValues()
{
	// Carefull if cpl_variable = primary variable -> need extra column in NodeValueTable !
	int nidx0 = GetNodeValueIndex(m_num->cpl_variable_JOD);
	for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
		SetNodeValue(l, nidx0, GetNodeValue(l, nidx0 + 1));
	//	if (_pcs_type_name.find("RICHARDS") != string::npos) { //WW
	if (this->getProcessType() == FiniteElement::RICHARDS_FLOW) // WW
	{
		nidx0 = GetNodeValueIndex("SATURATION1");
		for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
			SetNodeValue(l, nidx0, GetNodeValue(l, nidx0 + 1));
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB implementation
   02/2006 WW Modified for the cases of high order element and saturation
   08/2008 WW Make it twofold copy: forward and backward
**************************************************************************/
void CRFProcess::CopyTimestepNODValues(bool forward)
{
	bool Quadr = false; // WW
	if (type == 4 || type == 41)
		Quadr = true;

	for (int j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		int nidx0 = GetNodeValueIndex(pcs_primary_function_name[j]);
		int nidx1 = nidx0 + 1;
		if (!forward) // 08.2008. WW
		{
			nidx0++;
			nidx1--;
		}
		for (size_t l = 0; l < m_msh->GetNodesNumber(Quadr); l++)
			SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
		// WW
		//		if (_pcs_type_name.find("RICHARDS") != string::npos || type == 1212) { //Multiphase. WW
		// Multiphase. WW
		if (this->getProcessType() == FiniteElement::RICHARDS_FLOW || type == 1212 || type == 42)
		{
			if (j == 1 && (type == 1212 || type == 42)) // Multiphase. WW
				continue;
			if (j == 0)
				nidx0 = GetNodeValueIndex("SATURATION1");
			else if (j == 1)
				nidx0 = GetNodeValueIndex("SATURATION2");
			nidx1 = nidx0 + 1;
			if (!forward) // 27.08.2008. WW
			{
				nidx0++;
				nidx1--;
			}
			//
			for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
				SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   3/2012 JT. Based on the nodal version
**************************************************************************/
void CRFProcess::CopyTimestepELEValues(bool forward)
{
	size_t j, nvals, nidx0, nidx1;
	long iel, num_ele;
	bool copy_porosity = true;
	nvals = ele_val_name_vector.size();
	if (nvals < 2)
		return; // then we don't have 2 time levels of anything
	num_ele = (long)m_msh->ele_vector.size();
//
#ifdef GEM_REACT
	// do nothing as porosity update is handled by REACT_GEMS after!! flow and transport solution
	copy_porosity = false;
#endif
	//
	for (j = 0; j < nvals - 1; j++)
	{
		if (ele_val_name_vector[j].compare(ele_val_name_vector[j + 1])
		    != 0) // If not the same, then we only have a single time slot
			continue;
		if (ele_val_name_vector[j].find("POROSITY") != string::npos && !copy_porosity) // Porosity
			continue;
		//
		nidx0 = j;
		nidx1 = j + 1;
		if (!forward)
		{
			nidx0++;
			nidx1--;
		}
		for (iel = 0; iel < num_ele; iel++)
		{
			SetElementValue(iel, nidx0, GetElementValue(iel, nidx1));
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2006 SB Implementation
   last modified:
**************************************************************************/
int PCSGetPCSIndex(const string& pcs_type_name, const string& comp_name)
{
	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(pcs_type_name));

	CRFProcess* m_pcs = NULL;
	int i, pcs_no;
	int no_processes = (int)pcs_vector.size();
	string testname;
	pcs_no = -1;
	for (i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->pcs_type_name.compare(pcs_type_name) == 0) {
		if (m_pcs->getProcessType() == pcs_type)
		{
			testname = m_pcs->pcs_primary_function_name[0];
			if (testname.compare(comp_name) == 0)
			{
				//        cout << " Found in PCSGetbyTypeandCompName for PCSType/Compname " << pcs_type_name << ", " <<
				//        comp_name;
				//        cout << " Process number " << m_pcs->pcs_number << ", compnumber " <<
				//        m_pcs->pcs_component_number << "\n";
				pcs_no = i;
				return pcs_no;
			}
		}
	}
	return pcs_no;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2005 SB Implementation
   last modified:
   10/2010 TF restructured function
**************************************************************************/
CRFProcess* PCSGet(const std::string& pcs_type_name, const std::string& comp_name)
{
	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(pcs_type_name));
	size_t no_processes(pcs_vector.size());
	for (size_t i = 0; i < no_processes; i++)
		//		if (pcs_vector[i]->pcs_type_name.compare(pcs_type_name) == 0) { // TF
		if (pcs_vector[i]->getProcessType() == pcs_type)
			if (comp_name.compare(pcs_vector[i]->pcs_primary_function_name[0]) == 0)
				return pcs_vector[i];

	return NULL;
}

CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type, const std::string& comp_name)
{
	size_t no_processes(pcs_vector.size());
	for (size_t i = 0; i < no_processes; i++)
		if (pcs_vector[i]->getProcessType() == pcs_type)
			if (comp_name.compare(pcs_vector[i]->pcs_primary_function_name[0]) == 0)
				return pcs_vector[i];

	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   12/2005 OK Implementation
**************************************************************************/
CRFProcess* PCSGet(const string& var_name, bool bdummy)
{
	string pcs_var_name;
	CRFProcess* m_pcs = NULL;
	bdummy = bdummy; // WW
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
		{
			pcs_var_name = m_pcs->pcs_primary_function_name[j];
			if (pcs_var_name.compare(var_name) == 0)
				return m_pcs;
		}
		for (size_t j = 0; j < m_pcs->GetSecondaryVNumber(); j++)
		{
			pcs_var_name = m_pcs->pcs_secondary_function_name[j];
			if (pcs_var_name.compare(var_name) == 0)
				return m_pcs;
		}
	}
	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   05/2006 CMCD Implementation
**************************************************************************/
CRFProcess* PCSGetFluxProcess()
{
	CRFProcess* m_pcs = NULL;
	bool found = false;
	const size_t no_processes(pcs_vector.size());

	for (size_t i = 0; i < no_processes; i++)
	{
		//		if (pcs_vector[i]->_pcs_type_name == "LIQUID_FLOW") { // TF
		if (pcs_vector[i]->getProcessType() == FiniteElement::LIQUID_FLOW)
		{
			m_pcs = pcs_vector[i];
			found = true;
		}
		//		if (pcs_vector[i]->_pcs_type_name == "GROUNDWATER_FLOW") {
		if (pcs_vector[i]->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		{
			m_pcs = pcs_vector[i];
			found = true;
		}
		if (pcs_vector[i]->getProcessType() == FiniteElement::RICHARDS_FLOW)
		{
			m_pcs = pcs_vector[i];
			found = true;
		}
		if (found)
			return m_pcs;
	}
	return NULL;
}

/*************************************************************************
   GeoSys-FEM Function:
   Task:
   For Richards, dual-Richards and multi-phase flow.
   Programming:
   02/2007 WW Implementation
 **************************************************************************/
void CRFProcess::CalcSecondaryVariablesUnsaturatedFlow(bool initial)
{
	int ii, jj;
	long i;
	double p_cap;
	int idxp, idxcp, idxS;
	//----------------------------------------------
	vector<string> secondSNames;
	vector<string> secondCPNames;
	secondSNames.push_back("SATURATION1");
	if (type == 1212 || type == 42) // Multiphase flow
		secondCPNames.push_back("PRESSURE1");
	else
		secondCPNames.push_back("PRESSURE_CAP1");
	if (continuum_vector.size() == 2)
	{
		secondSNames.push_back("SATURATION2");
		secondCPNames.push_back("PRESSURE_CAP2");
	}

	//
	CElem* elem = NULL;
	CFiniteElementStd* fem = GetAssembler();
	//----------------------------------------------------------------------
	for (ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		idxS = GetNodeValueIndex(secondSNames[ii].c_str()) + 1;
		if (type == 1212 || type == 42) // Multiphase flow
		{
			jj = ii;
			if (type == 42)
				jj += problem_dimension_dm;
			idxcp = GetNodeValueIndex(pcs_primary_function_name[jj]) + 1;
			idxp = GetNodeValueIndex("PRESSURE_W");
		}
		else
		{
			idxp = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			idxcp = GetNodeValueIndex(secondCPNames[ii].c_str());
		}
		//----------------------------------------------------------------------
		// Capillary pressure
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			// Copy the new saturation to the old level
			//
			if (type == 1212 || type == 42) // Multiphase flow
			{
				// p_w = p_g-p_c
				p_cap = GetNodeValue(i, idxcp + 2) - GetNodeValue(i, idxcp);
				SetNodeValue(i, idxp, p_cap);
			}
			else
			{
				p_cap = -GetNodeValue(i, idxp);
				SetNodeValue(i, idxcp, p_cap);
			}
			SetNodeValue(i, idxS, 0.0);
		}
	}
	// Cal S
	for (ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		continuum = ii;
		for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
		{
			elem = m_msh->ele_vector[i];
			if (elem->GetMark()) // Marked for use
			{
				elem->SetOrder(false);
				fem->CalcSaturation(*elem);
			}
		}
	}
	//
	if (!initial)
		return;
	//----------
	for (ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		idxS = GetNodeValueIndex(secondSNames[ii].c_str()) + 1;
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
			SetNodeValue(i, idxS - 1, GetNodeValue(i, idxS));
	}
}

/*************************************************************************
   GeoSys-FEM Function:
Task: Updating rho_s values in TNEQ

Programming:
11/2011 HS Implementation for thermal storage project
**************************************************************************/
void CRFProcess::CalcSecondaryVariablesTNEQ()
{
	CElem* elem = NULL;
	size_t i;

	// loop over all the nodes,
	// clean the nodal reaction rate values
	const size_t node_vector_size(m_msh->nod_vector.size());
	const int idx_nodal_react_rate = this->GetNodeValueIndex("REACT_RATE_N");
	const int idx_nodal_solid_density = this->GetNodeValueIndex("SOLID_DENSITY_N");
	
	const int idx_nodal_solid_cP = this->GetNodeValueIndex("SOLID_HEAT_CAP_N");

	for (size_t idx_node = 0; idx_node < node_vector_size; idx_node++)
	{
		this->SetNodeValue(idx_node, idx_nodal_react_rate, 0.0);
		this->SetNodeValue(idx_node, idx_nodal_solid_density, 0.0);
		this->SetNodeValue( idx_node, idx_nodal_solid_cP, 0.0);
	}

	// loop over all the elements and update the rho_s values.
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	for (i = 0; i < mesh_ele_vector_size; i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			fem->UpdateSolidDensity(elem->GetIndex()); // HS, thermal storage reactions
			fem->ExtrapolateGauss_ReactRate_TNEQ_TES(*elem, this); // HS added 19.02.2013
		}
	}
}

/*************************************************************************
   GeoSys-FEM Function:
Task: Updating rho_s values in TES

Programming:
11/2011 HS Implementation for thermal storage project
**************************************************************************/
void CRFProcess::CalcSecondaryVariablesTES()
{
	// TODO [CL] merge with TNEQ
	// loop over all the nodes,
	// clean the nodal reaction rate values
	const size_t node_vector_size = m_msh->nod_vector.size();
	const int idx_nodal_react_rate = GetNodeValueIndex("REACT_RATE_N");
	const int idx_nodal_solid_density = GetNodeValueIndex("SOLID_DENSITY_N");
	const int idx_nodal_solid_cP = this->GetNodeValueIndex("SOLID_HEAT_CAP_N");

	for (size_t idx_node = 0; idx_node < node_vector_size; idx_node++)
	{
		SetNodeValue(idx_node, idx_nodal_react_rate, 0.0);
		SetNodeValue(idx_node, idx_nodal_solid_density, 0.0);
		this->SetNodeValue( idx_node, idx_nodal_solid_cP, 0.0);
	}

	// loop over all the elements and update the rho_s values.
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
	{
		CElem* const elem = m_msh->ele_vector[i];
		if (elem->GetMark()) // Marked for use
		{
			fem->UpdateSolidDensity(elem->GetIndex()); // HS, thermal storage reactions
			fem->ExtrapolateGauss_ReactRate_TNEQ_TES(*elem, this); // HS added 19.02.2013
		}
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   Task: Updating secondary variables for Multiphase flow in PS_GLOBAL

   Programming:
   03/2009 PCH Implementation
 **************************************************************************/
void CRFProcess::CalcSecondaryVariablesPSGLOBAL()
{
	long i;
	int ndx_pressure1, ndx_p_cap, ndx_pressure2, ndx_s_wetting; // WW, ndx_s_nonwetting;

	// The primary variables
	ndx_pressure1 = GetNodeValueIndex("PRESSURE1");
	// WW ndx_s_nonwetting = GetNodeValueIndex("SATURATION2");

	// The secondary variables
	ndx_pressure2 = GetNodeValueIndex("PRESSURE2");
	ndx_p_cap = GetNodeValueIndex("PRESSURE_CAP");
	ndx_s_wetting = GetNodeValueIndex("SATURATION1");

	double pressure1, pressure2, p_cap, s_wetting;
	for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
	{
		pressure1 = GetNodeValue(i, ndx_pressure1 + 1); // New
		pressure2 = GetNodeValue(i, ndx_pressure1); // Old

		// Let's get capillary pressure before updating pressure2
		// by accessing the primary variable of the saturation equation
		// not the secondary variable of it.
		int ndx_sat2 = GetNodeValueIndex("SATURATION2");
		double sat2 = GetNodeValue(i, ndx_sat2 + 1);
		// Due to the iterative solution scheme in solving Snw with no
		// explicit boundary condition for non-zero flux condition,
		// Snw may become negative particularly the density difference
		// between two fluids is big. To prevent negative Snw, the
		// saturation restriction added.
		CMediumProperties* mmp = NULL;
		if (sat2 != .0 && mmp_vector.size() > 1)
		{
			double sum = 0.0;
			CNode* thisNode = m_msh->nod_vector[i];
			int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();
			// Harmonic mean
			for (int i = 0; i < NumOfNeighborElements; ++i)
			{
				// Mount neighboring elemenets and get the corresponding material group one by one.
				size_t eleIdx = thisNode->getConnectedElementIDs()[i];
				CElem* thisEle = m_msh->ele_vector[eleIdx];
				size_t matgrp = thisEle->GetPatchIndex();
				mmp = mmp_vector[matgrp];
				mmp->mode = 2;
				sum += 1.0 / sat2;
			}
			sat2 = (double)NumOfNeighborElements / sum;
		}
		else
			mmp = mmp_vector[0];
		s_wetting = 1.0 - sat2;
		// Assigning the secondary variable, Sw
		SetNodeValue(i, ndx_s_wetting, s_wetting);
		// Assigning the primary variable Snw here one more time
		// to completely bound the range of saturation
		//	SetNodeValue(i,ndx_s_nonwetting,sat2);
		//	SetNodeValue(i,ndx_s_nonwetting+1,sat2);

		// Assigning the secondary variable, Pc
		if ((mmp->capillary_pressure_model == 4 || mmp->capillary_pressure_model == 6)
		    && fabs(mmp->capillary_pressure_values[0]) < DBL_EPSILON)
		{
			p_cap = 0.;
		}
		else
		{
			if (mmp_vector.size() > 1)
				p_cap = GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0 - sat2);
			else
				p_cap = mmp->CapillaryPressureFunction(1.0 - sat2);
		}

		SetNodeValue(i, ndx_p_cap, p_cap);

		pressure2 = pressure1 + p_cap;
		// Assigning the secondary variables, Pnw
		SetNodeValue(i, ndx_pressure2, pressure2);
	}
}

/*************************************************************************
GeoSys-FEM Function:
Task: Updating secondary variables for Multiphase flow in PS_GLOBAL

Programming:
03/2009 PCH Implementation
**************************************************************************/
void CRFProcess::CalcSecondaryVariablesLiquidFlow()
{
	long i;
	int ndx_dens;
	double var[3] = {0, 0, 0};
	bool heattransport = false;

	ndx_dens = GetNodeValueIndex("DENSITY1");

	CRFProcess* m_pcs = NULL;
	if ((m_pcs = PCSGet("HEAT_TRANSPORT")))
		heattransport = true;

	CFluidProperties* m_mfp = NULL;
	m_mfp = mfp_vector[0];

	double dens;
	for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
	{
		// get pressure
		var[0] = this->GetNodeValue(i, this->GetNodeValueIndex("PRESSURE1"));
		// get temperature
		if (heattransport)
			var[1] = m_pcs->GetNodeValue(i, m_pcs->GetNodeValueIndex("TEMPERATURE1"));
		// Set salinity
		var[2] = 0.0;
		dens = m_mfp->Density(var);
		// Assigning the secondary variable
		SetNodeValue(i, ndx_dens, dens);
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate saturation on node by averaging the patches of the
         neighboring elements in three means

   0: Arithmatic mean
   1: Geomtric mean
   2: Harmonic mean

   Programing:
   03/2009 PCH Implementation
   last modification:
 *************************************************************************/
double CRFProcess::GetCapillaryPressureOnNodeByNeighobringElementPatches(int nodeIdx, int meanOption, double Sw)
{
	double p_cap = 0.0, sum = 0.0;

	CNode* thisNode = m_msh->nod_vector[nodeIdx];
	int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();

	switch (meanOption)
	{
		case 0:
			break;
		case 1:
			break;
		case 2: // Harmonic mean
			for (int i = 0; i < NumOfNeighborElements; ++i)
			{
				// Mount neighboring elemenets and get the corresponding material group one by one.
				int eleIdx = thisNode->getConnectedElementIDs()[i];
				CElem* thisEle = m_msh->ele_vector[eleIdx];
				int matgrp = thisEle->GetPatchIndex();
				CMediumProperties* mmp = mmp_vector[matgrp];
				sum += 1.0 / mmp->CapillaryPressureFunction(Sw);
			}
			p_cap = (double)NumOfNeighborElements / sum;
			break;

		default:
			cout << "Please define the option for various means!"
			     << "\n";
			cout << "The code stops at GetCapillaryPressureOnNodeByNeighobringElementPatches function!"
			     << "\n";
			abort();
			break;
	}

	return p_cap;
}

/*************************************************************************
   GeoSys-FEM Function:
   Task: Calculates saturation for richards flow,
      alternative to CRFProcess::CalcSecondaryVariablesRichards
      uses pressure of only one node
     evoked by switch case in pcs-file: SATURATION_SWITCH = true
   Programming:
   06/2007 JOD Implementation
 **************************************************************************/
void CRFProcess::CalcSaturationRichards(int timelevel, bool update)
{
	double p_cap, saturation, volume_sum;
	int idxp, idxcp, idxS, idx_tS = -1;
	const size_t number_continuum(continuum_vector.size());
	size_t i_s, i_e;

	CMediumProperties* m_mmp = NULL;
	CElem* elem = NULL;
	// WW  CFiniteElementStd* fem = GetAssembler();

	if (continuum_ic) // Create IC: for both continua
	{
		i_s = 0;
		i_e = number_continuum;
	}
	else
	{
		i_s = continuum;
		i_e = continuum + 1;
	}

	for (size_t i_pv = i_s; i_pv < i_e; i_pv++)
	{
		idxp = GetNodeValueIndex(pcs_primary_function_name[i_pv]) + timelevel;
		idxS = GetNodeValueIndex(pcs_secondary_function_name[i_pv * number_continuum]) + timelevel;
		idxcp = GetNodeValueIndex(pcs_secondary_function_name[i_pv * number_continuum + number_continuum * 2])
		        + timelevel;
		if (continuum_vector.size() > 1)
			idx_tS = GetNodeValueIndex("TOTAL_SATURATION") + timelevel;

		for (long i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			// Capillary pressure
			p_cap = -GetNodeValue(i, idxp);
			if (timelevel == 1 && update)
				SetNodeValue(i, idxcp - 1, GetNodeValue(i, idxcp));
			SetNodeValue(i, idxcp, p_cap);
			if (timelevel == 1 && update)
				SetNodeValue(i, idxS - 1, GetNodeValue(i, idxS));

			// Liquid saturation
			if (continuum_vector.size() > 1)
				SetNodeValue(i, idx_tS, 0.0);
			//

			saturation = 0., volume_sum = 0.;
			size_t elemsCnode = m_msh->nod_vector[i]->getConnectedElementIDs().size();

			for (size_t j = 0; j < elemsCnode; j++)
			{
				elem = m_msh->ele_vector[m_msh->nod_vector[i]->getConnectedElementIDs()[j]];
				m_mmp = mmp_vector[elem->GetPatchIndex()];
				volume_sum += elem->volume;
				saturation += m_mmp->SaturationCapillaryPressureFunction(p_cap) * elem->volume;
			}
			saturation /= volume_sum;
			SetNodeValue(i, idxS, saturation);
		}
	}

	if (continuum > 0)
		for (long i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			double total_S = 0;
			for (size_t j = 0; j < continuum_vector.size(); j++)
			{
				idxS = GetNodeValueIndex(pcs_secondary_function_name[j * number_continuum]) + timelevel;
				total_S += GetNodeValue(i, idxS) * continuum_vector[j];
			}
			SetNodeValue(i, idx_tS, total_S);
		}
}

/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary node values
                  of process pcs_name and for variable var_name; old and new timelevel
    01/2006   SB    Implementation
    02/2008   CB    generalization
**************************************************************************/
double PCSGetEleMeanNodeSecondary_2(long index, int pcsT, const string& var_name, int timelevel)
{
	double val = 1.0; // As this returns saturation, default is fully saturated = 1.0;
	int idx = 0, j; // OK411
	long enode;
	CRFProcess* m_pcs = NULL;
	CRFProcess* cplpcs = NULL;
	CElem* elem = NULL;

	// Get index of secondary node value
	switch (pcsT)
	{
		case 0: // Liquid_Flow
			break;
		case 1: // Groundwater Flow
			break;
		case 66: // Overland Flow
			break;
		case 5: // Air Flow
			break;
		case 11: // Componental Flow
			break;
		case 1212: // Multiphase Flow
			m_pcs = PCSGet("MULTI_PHASE_FLOW"); // SB, BG
			if (m_pcs)
			{
				idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
				cplpcs = m_pcs;
			}
			break;
		case 12: // Two_phase_Flow
			m_pcs = PCSGet(FiniteElement::TWO_PHASE_FLOW);
			if (m_pcs)
			{
				if (m_pcs->pcs_type_number == 0)
					cplpcs = pcs_vector[m_pcs->pcs_number + 1];
				else if (m_pcs->pcs_type_number == 1)
					cplpcs = pcs_vector[m_pcs->pcs_number - 1];
				idx = cplpcs->GetNodeValueIndex(var_name) + timelevel;
			}
			break;
		case 1313: // Two_phase_Flow
			m_pcs = PCSGet("PS_GLOBAL");
			if (m_pcs)
			{
				idx = m_pcs->GetNodeValueIndex(
				    var_name); // +timelevel; CB SATURATION1 in PS_GLOBAL is only for old time level
				cplpcs = m_pcs;
			}
			break;
		case 22: // Richards flow
			m_pcs = PCSGet(FiniteElement::RICHARDS_FLOW);
			if (m_pcs)
			{
				idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
				cplpcs = m_pcs;
			}
			break;
		default:
			break;
	}

	if (m_pcs)
	{
		// Get element with index index
		elem = m_pcs->m_msh->ele_vector[index];
		val = 0.0;
		for (j = 0; j < elem->GetVertexNumber(); j++) // average all adjoining nodes
		{
			enode = elem->GetNodeIndex(j);
			val += cplpcs->GetNodeValue(enode, idx);
		}
		val = val / ((double)elem->GetVertexNumber());
	}
	return val;
}

/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary node values
                  of process pcs_name and for variable var_name; old and new timelevel
    01/2006   SB    Implementation
**************************************************************************/
double PCSGetEleMeanNodeSecondary(long index, const string& pcs_name, const string& var_name, int timelevel)
{
	double val = 1.0; // As this returns saturation, default is fully saturated = 1.0;
	int idx, j;
	long enode;
	CRFProcess* m_pcs = NULL;

	// Get process by process name
	FiniteElement::ProcessType pcs_type(FiniteElement::convertProcessType(pcs_name));
	m_pcs = PCSGet(pcs_type);
	if (m_pcs)
	{
		// Get index of secondary node value
		idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
		// Get element with index index
		CElem* elem = NULL;
		elem = m_pcs->m_msh->ele_vector[index];
		val = 0.0;
		// average all adjoining nodes
		for (j = 0; j < elem->GetVertexNumber(); j++)
		{
			enode = elem->GetNodeIndex(j);
			val += m_pcs->GetNodeValue(enode, idx);
		}
		val = val / ((double)elem->GetVertexNumber());
	}
	return val;
}

/*************************************************************************
   GeoSys-FEM Function:
   01/2006 OK Implementation
 **************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void CRFProcess::SetNODFlux()
{
	long i;
	//----------------------------------------------------------------------
	int nidx;
	nidx = GetNodeValueIndex("FLUX");
	if (nidx < 0)
		return;
	double m_val;
	for (i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_val = eqs->b[i]; //? m_nod->eqs_index
		SetNodeValue(i, nidx, m_val);
	}
	//----------------------------------------------------------------------
}
#endif

/*************************************************************************
   GeoSys-FEM Function:
   01/2006 OK Implementation
 **************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void CRFProcess::AssembleParabolicEquationRHSVector()
{
	// OK  long i;
	//----------------------------------------------------------------------
	// Init
	/*	PCH & WW
	   for(i=0;i<eqs->dim;i++)
	   {
	    eqs->b[i] = 0.0;
	   }
	   //----------------------------------------------------------------------
	   CElem* m_ele = NULL;
	   for(i=0;i<(long)m_msh->ele_vector.size();i++)
	   {
	    m_ele = m_msh->ele_vector[i];
	    if(m_ele->GetMark()) // Marked for use
	   {
	   fem->ConfigElement(m_ele,false);
	   fem->AssembleParabolicEquationRHSVector();
	   //fem->AssembleParabolicEquationLHSMatrix();
	   }
	   }
	 */
	//----------------------------------------------------------------------
}
#endif //#ifndef NEW_EQS //WW. 07.11.2008
/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   02/2008 JOD removed
   03/2008 HS/KG activated for adaptive time step
   Reload primary variable
 **************************************************************************/
void CRFProcess::PrimaryVariableReload()
{
	//  char pcsT;
	//  pcsT = _pcs_type_name[0];
	//  switch(pcsT){
	//    case 'L':
	//      break;
	//    case 'U':
	//      break;
	//    case 'G':
	//      break;
	//    case 'T':
	//      break;
	//    case 'C':
	//      break;
	//    case 'M':
	//      PrimaryVariableReloadTransport();
	//      break;
	//    case 'R': // Richards flow
	//      PrimaryVariableReloadRichards();
	//      break;
	//  }

	switch (this->getProcessType())
	{
		case FiniteElement::MASS_TRANSPORT:
			PrimaryVariableReloadTransport();
			break;
		case FiniteElement::RICHARDS_FLOW:
			PrimaryVariableReloadRichards();
			break;
		default:
			break;
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   02/2008 JOD removed
   Reload primary variable of Richards Flow
 **************************************************************************/
void CRFProcess::PrimaryVariableReloadRichards()
{
	size_t i;
	int idxp, idx_storage;
	double storage_p;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
	idx_storage = GetNodeValueIndex("STORAGE_P");
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		storage_p = GetNodeValue(i, idx_storage);
		SetNodeValue(i, idxp, storage_p);
		SetNodeValue(i, idxp + 1, storage_p);
	}
	CalcSecondaryVariables(0);
	CalcSecondaryVariables(1);
}

/*************************************************************************
   GeoSys-FEM Function:
   12/2007 kg44 Implementation
   Reload primary variable for Transport
 **************************************************************************/
void CRFProcess::PrimaryVariableReloadTransport()
{
	size_t i;
	int idxp, idx_storage;
	double conc_back;
	char* mcomp_name;

	// kg44 test
	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);

	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%i", "CONC_BACK_", pcs_component_number);

	idx_storage = GetNodeValueIndex(mcomp_name);

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		conc_back = GetNodeValue(i, idx_storage);
		SetNodeValue(i, idxp, conc_back);
		SetNodeValue(i, idxp + 1, conc_back);
	}

	CalcSecondaryVariables(0);
	CalcSecondaryVariables(1);
}

/*************************************************************************
   GeoSys-FEM Function:
   12/2007 kg44 Implementation
   Reload primary variable for Transport
 **************************************************************************/
void CRFProcess::PrimaryVariableStorageTransport()
{
	size_t i;
	int idxp, idx_storage;
	double concentration;
	char* mcomp_name;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%i", "CONC_BACK_", pcs_component_number);
	//   cout << "mcomp_name"<< mcomp_name<<"\n";
	idx_storage = GetNodeValueIndex(mcomp_name);

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		concentration = GetNodeValue(i, idxp);
		SetNodeValue(i, idx_storage, concentration);
		//    SetNodeValue(i,idx_storage+1,concentration);
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   Reload primary variable of Richards Flow
 **************************************************************************/
void CRFProcess::PrimaryVariableStorageRichards()
{
	size_t i;
	int idxp, idx_storage;
	double pressure;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]) + 1;
	idx_storage = GetNodeValueIndex("STORAGE_P");
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		pressure = GetNodeValue(i, idxp);
		SetNodeValue(i, idx_storage, pressure);
		SetNodeValue(i, idx_storage + 1, pressure);
	}
}

//*************************************************************************
// GeoSys-FEM Function:
// 12/2007 kg44 Implementation
// check change of concentration and set new time step factor
//**************************************************************************/
#ifdef kg44 // WW
double CRFProcess::GetNewTimeStepSizeTransport(double mchange)
{
	size_t i;
	long mnode = -1;
	int comp;
	int idxn, idxo;
	double conc_new, conc_old, /*time_coeff,*/ max_change = 1.0e-10, tchange = 1.0;
	char* mcomp_name;

	idxo = GetNodeValueIndex(pcs_primary_function_name[0]);
	comp = pcs_component_number; // get component number
	//   cout << "comp number "<<comp<<"\n";
	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%li", "CONC_BACK_", comp);
	//   cout << "mcomp_name"<< mcomp_name<<"\n";
	idxn = GetNodeValueIndex(mcomp_name);
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		conc_old = abs(GetNodeValue(i, idxo));
		conc_new = abs(GetNodeValue(i, idxn));
		if (((conc_old) > MKleinsteZahl) && ((conc_new) > MKleinsteZahl))
		{
			max_change = MMax(max_change, conc_new / conc_old);
			mnode = i;
		}
	}
	tchange = mchange / max_change;
	if (tchange > 2.0)
		tchange = 2.0;
	cout << "Transport: max change of " << max_change << " at node " << mnode << " factor " << tchange << "\n";
	return tchange;
}
#endif

/**************************************************************************
   FEMLib-Method:
   11/2005 MB Implementation
   03/2006 OK 2nd version (binary coupling)
**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void CRFProcess::SetCPL()
{
	int i;
	double value = 0.0;
	//----------------------------------------------------------------------
	// Nothing to do
	if ((cpl_type_name.size() == 0) || (cpl_type_name.compare("PARTITIONED") == 0))
		return;
	//----------------------------------------------------------------------
	// PCS CPL
	CRFProcess* m_pcs_cpl = PCSGet(cpl_type_name);
	if (!m_pcs_cpl)
	{
		cout << "Fatal error in CRFProcess::SetCPL: no PCS data"
		     << "\n";
		return;
	}
	//----------------------------------------------------------------------
	// MSH data for PCS CPL
	CFEMesh* m_msh_cpl = m_pcs_cpl->m_msh;
	if (!m_msh_cpl)
	{
		cout << "Fatal error in CRFProcess::SetCPL: no MSH data"
		     << "\n";
		return;
	}
	//----------------------------------------------------------------------
	// GEO data for PCS CPL
	Surface* m_sfc = GEOGetSFCByName(m_msh_cpl->geo_name);
	if (!m_sfc)
	{
		cout << "Fatal error in CRFProcess::SetCPL: no GEO data"
		     << "\n";
		return;
	}
	//----------------------------------------------------------------------
	//......................................................................
	// MSH nodes of PCS CPL
	cout << "      ->CPL: " << cpl_type_name << ": ";
	vector<long> cpl_msh_nodes_vector;
	m_msh_cpl->GetNODOnSFC(m_sfc, cpl_msh_nodes_vector);
	if ((int)cpl_msh_nodes_vector.size() == 0)
		cout << "Warning in CRFProcess::SetCPL: no MSH nodes found"
		     << "\n";
	cout << "CPL nodes = " << (int)cpl_msh_nodes_vector.size() << "\n";
	//.....................................................................-
	// MSH nodes of PCS
	cout << "      ->CPL: " << convertProcessTypeToString(this->getProcessType()) << ": ";
	vector<long> msh_nodes_vector;
	m_msh->GetNODOnSFC(m_sfc, msh_nodes_vector);
	if ((int)msh_nodes_vector.size() == 0)
		cout << "Warning in CRFProcess::SetCPL: no MSH nodes found"
		     << "\n";
	cout << "CPL nodes = " << (int)msh_nodes_vector.size() << "\n";
	//----------------------------------------------------------------------
	if (m_msh_cpl->pcs_name.compare("RICHARDS_FLOW") == 0)
	{
		m_msh->SetNODPatchAreas();
		int nidx = GetNodeValueIndex("WDEPTH");
		long st_node_number;
		double st_node_value = 0.0;
		CNode* m_nod = NULL;
		for (i = 0; i < (int)msh_nodes_vector.size(); i++)
		{
			value = -2.314e-02;
			st_node_number = msh_nodes_vector[i];
			m_nod = m_msh->nod_vector[st_node_number];
			st_node_value = GetNodeValue(st_node_number, nidx);
			st_node_value /= m_nod->patch_area;
			value *= st_node_value;
			// cout << "CPL value = " << value << "\n";
			eqs->b[st_node_number] += value;
		}
	}

	//	if (_pcs_type_name.compare("RICHARDS_FLOW") == 0
	//				&& m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0) { // ToDo
	if (this->getProcessType() == FiniteElement::RICHARDS_FLOW
	    // ToDo
	    && m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0)
	{
		long msh_node_number;
		// WW long cpl_msh_nod_number;
		long cpl_msh_ele_number;
		value = 0.0;
		cout << "CPL value = " << value << "\n";
		// PCS-CON
		CRFProcess* m_pcs_cond = PCSGet(cpl_type_name);
		// int nidx = m_pcs_cond->GetNodeValueIndex(m_pcs_cond->pcs_primary_function_name[0]);
		int cpl_nidx = m_pcs_cond->GetNodeValueIndex("WDEPTH");
		//----------------------------------------------------------------------
		// ELE of PCS_CPL related to NOD of PCS
		//  CFEMesh* m_msh_this = MSHGet("RICHARDS_FLOW_LOCAL");
		CElem* m_ele_cnd = NULL;
		//----------------------------------------------------------------------
		//  CSourceTermGroup *m_st_group = NULL;
		//  CSourceTerm *m_st = NULL;
		//  m_st_group = STGetGroup(_pcs_type_name,pcs_primary_function_name[0]);
		//----------------------------------------------------------------------
		double cpl_ele_val = 0.0;
		size_t j;
		CNodeValue* cnodev = NULL;
		//  for(i=0;i<(int)m_st_group->group_vector.size();i++)
		//  ofstream st_out_file("st_out_file.txt",ios::app);
		for (i = 0; i < (int)st_node_value.size(); i++)
		{
			cnodev = st_node_value[i];
			// MSH-PCS
			// m_nod = m_msh_this->nod_vector[m_st_group->group_vector[i]->msh_node_number];
			// m_st_group->group_vector[i]->msh_node_number; //0
			msh_node_number = cnodev->msh_node_number;
			// MSH-PCS-CPL
			// WW cpl_msh_nod_number = msh_node_number;
			cpl_msh_ele_number = pcs_number; // OK:TODO
			m_ele_cnd = m_pcs_cond->m_msh->ele_vector[cpl_msh_ele_number];
			for (j = 0; j < m_ele_cnd->GetNodesNumber(false); j++)
				cpl_ele_val += m_pcs_cond->GetNodeValue(m_ele_cnd->nodes_index[j], cpl_nidx);
			cpl_ele_val /= m_ele_cnd->GetNodesNumber(false);
			// VAL-CON
			value = 2.314e-02 * cpl_ele_val * 1e-2;
			//    st_out_file << value << "\n";
			// EQS-RHS
			eqs->b[msh_node_number] += value;
		}
		//----------------------------------------------------------------------
		/*
		   CNodeValue* m_node_value = NULL;
		   m_node_value = new CNodeValue();
		   m_node_value->msh_node_number = msh_node_number;
		   m_node_value->geo_node_number = m_pnt->id;
		   m_node_value->node_value = value;
		   CSourceTermGroup *m_st_group = NULL;
		   m_st_group = STGetGroup(_pcs_type_name,pcs_primary_function_name[0]);
		   m_st_group->group_vector.push_back(m_node_value);
		   m_st_group->st_group_vector.push_back(m_st); //OK
		 */
	}

	//	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0
	//				&& m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0) { // ToDo
	if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW
	    // ToDo
	    && m_msh_cpl->pcs_name.compare("OVERLAND_FLOW") == 0)
	{
		long ie = (long)msh_nodes_vector.size() / (m_msh->getNumberOfMeshLayers() + 1);
		long of_node_number, gf_node_number;
		double of_node_value, gf_node_value;
		//  CNode* m_nod = NULL;
		int of_nidx = GetNodeValueIndex("WDEPTH");

		for (i = 0; i < ie; i++)
		{
			of_node_number = msh_nodes_vector[i]; // ToDo
			of_node_value = m_pcs_cpl->GetNodeValue(of_node_number, of_nidx);
			// m_nod = m_msh->nod_vector[gf_node_number];
			// st_node_value /= m_nod->patch_area;
			gf_node_value = of_node_value * 2e-11;
			if (gf_node_value > 1e-13)
				cout << "CPL value = " << gf_node_value << "\n";
			gf_node_number = msh_nodes_vector[i];
			eqs->b[gf_node_number] += gf_node_value;
		}
	}

	int idx = GetNodeValueIndex("FLUX") + 1;
	// for(i=0;i<(int)msh_nodes_vector.size();i++)
	for (i = 0; i < 1; i++)
		SetNodeValue(msh_nodes_vector[i], idx, value);
}
#endif

/**************************************************************************
   PCSLib-Method:
   04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateBCGroup()
{
	cout << "->Create BC" << '\n';
	CBoundaryConditionsGroup* m_bc_group = NULL;

	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
	{
		BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
		m_bc_group = new CBoundaryConditionsGroup();
		// OK
		m_bc_group->setProcessTypeName(pcs_type_name);
		// OK
		m_bc_group->setProcessPrimaryVariableName(pcs_primary_function_name[i]);
		m_bc_group->Set(this, Shift[i]);
		bc_group_list.push_back(m_bc_group);
	}
}

/**************************************************************************
   PCSLib-Method:
   04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateSTGroup()
{
	cout << "->Create ST" << '\n';
	CSourceTermGroup* m_st_group = NULL;
	// WW
	std::ifstream* iSourceNBC_RHS_file = NULL;
	std::ofstream* oSourceNBC_RHS_file = NULL;

	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));

	if (WriteSourceNBC_RHS == 1)
	{
		string m_file_name = FileName + "_" + pcs_type_name + "_source_Neumann_RHS.bin";
		oSourceNBC_RHS_file = new ofstream(m_file_name.c_str(), ios::binary);
	}
	else if (WriteSourceNBC_RHS == 2)
	{
		string m_file_name = FileName + "_" + pcs_type_name + "_source_Neumann_RHS.bin";
		iSourceNBC_RHS_file = new ifstream(m_file_name.c_str(), ios::binary);
		if (!iSourceNBC_RHS_file->good())
			cout << "_source_Neumann_RHS file is not found"
			     << "\n";
	}

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
	{
		// OK m_st_group = m_st_group->Get(pcs_primary_function_name[i]);
		m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
		if (!m_st_group)
		{
			m_st_group = new CSourceTermGroup();
			// OK
			m_st_group->pcs_type_name = pcs_type_name;
			// OK
			m_st_group->pcs_pv_name = pcs_primary_function_name[i];
			//      if(iSourceNBC_RHS_file)  // Read from data. WW
			//        m_st_group->Read(*iSourceNBC_RHS_file);
			//	  else
			m_st_group->Set(this, Shift[i]);
			st_group_list.push_back(m_st_group);
		}
	}
	if (oSourceNBC_RHS_file) // WW
		//    WriteRHS_of_ST_NeumannBC(*oSourceNBC_RHS_file);

		if (iSourceNBC_RHS_file) // WW
		{
			iSourceNBC_RHS_file->close();
			delete iSourceNBC_RHS_file;
			iSourceNBC_RHS_file = NULL;
		}
	if (oSourceNBC_RHS_file) // WW
	{
		oSourceNBC_RHS_file->close();
		delete oSourceNBC_RHS_file;
		oSourceNBC_RHS_file = NULL;
	}
}

/**************************************************************************
    PCSLib-Method:
    08/2006 OK Implementation
    04/2012 BG Extension to 2 Phases
**************************************************************************/
void CRFProcess::CalcELEFluxes(const GEOLIB::Polyline* const ply, double* result)
{
	int coordinateflag, dimension = 0, axis = 0;
	bool Edge_already_used;
	// bool Node_already_used;
	bool Use_Element;
	vector<size_t> vecConsideredEdges;
	vector<CNode*> vec_nodes_edge;
	bool Edge_on_Geo = false, Point_on_Geo;
	std::vector<size_t> nod_vector_at_geo;
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	vec<CEdge*> ele_edges_vector(15);
	double edg_normal_vector[3];
	double vn;
	double edg_length = 0.0;
	double vn_vec[3];
	double edge_vector[3];
	double f_n_sum = 0.0;
	vec<long> element_nodes(20);
	double f[3], v[3];
	double v_2[3], f_2[3];
	double f_n_sum_2 = 0.0;

	result[0] = result[1] = 0;
	FiniteElement::ProcessType pcs_type(getProcessType());

	CRFProcess* m_pcs_flow = NULL;
	//	if (_pcs_type_name.find("FLOW") != string::npos) {
	if (isFlowProcess(this->getProcessType()))
		m_pcs_flow = this;
	else
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);

	// calculates element velocity based on 1 GP
	// CalcELEVelocities();

	int v_eidx[3];
	int v_eidx_2[3];
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");

	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
	{
		v_eidx_2[0] = m_pcs_flow->GetElementValueIndex("VELOCITY2_X");
		v_eidx_2[1] = m_pcs_flow->GetElementValueIndex("VELOCITY2_Y");
		v_eidx_2[2] = m_pcs_flow->GetElementValueIndex("VELOCITY2_Z");
	}

	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			std::cout << i << " " << v_eidx[i] << "Velocity output is not specified"
			          << "\n";
			exit(0); // return 0.0;
		}
	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
		for (size_t i = 0; i < 3; i++)
			if (v_eidx_2[i] < 0)
			{
				std::cout << i << " " << v_eidx[i] << "Velocity output is not specified"
				          << "\n";
				exit(0); // return 0.0;
			}

	// determine the dimension and the orientation of the mesh
	coordinateflag = m_msh->GetCoordinateFlag();
	if (coordinateflag == 10)
	{
		dimension = 1;
		axis = 0;
	} // x only
	else if (coordinateflag == 11)
	{
		dimension = 1;
		axis = 1;
	} // y only
	else if (coordinateflag == 12)
	{
		dimension = 1;
		axis = 2;
	} // z only
	else if (coordinateflag == 21)
	{
		dimension = 2;
		axis = 1;
	} // x, y only
	else if (coordinateflag == 22)
	{
		dimension = 2;
		axis = 2;
	} // x, z only
	else if (coordinateflag == 32)
	{
		dimension = 3;
		axis = 3;
	} // x, y, z only

	// Get elements at GEO
	std::vector<size_t> ele_vector_at_geo;
	m_msh->GetELEOnPLY(ply, ele_vector_at_geo, true);
	// BG: 04/2011 nodes are needed to provide the correct edge
	m_msh->GetNODOnPLY(ply, nod_vector_at_geo);

	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		m_ele->SetNormalVector();
		m_ele->GetNodeIndeces(element_nodes);
		Use_Element = true;

		// Configure Element for interpolation of node velocities to GP velocities
		fem->ConfigElement(m_ele);
		// velocity vector
		for (size_t j = 0; j < 3; j++)
		{
			// v[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[j]);
			// Calculate Element velocity
			v[j] = fem->Get_Element_Velocity(m_ele->GetIndex(), m_pcs_flow, 0, j);
		}
		// Test mit Knotengeschwindigkeiten
		// double temp_v[3];
		// temp_v[0] = temp_v[1] = temp_v[2] = 0.0;
		// int variable_index[3];
		// variable_index[0] = m_pcs_flow->GetNodeValueIndex("VELOCITY_X1");
		// variable_index[1] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Y1");
		// variable_index[2] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Z1");
		//
		// for (size_t j = 0; j < 3; j++)
		//{
		//	for (size_t k = 0; k < m_ele->GetNodesNumber(false); k++)
		//	{
		//		temp_v[j] += m_pcs_flow->GetNodeValue(element_nodes[k], variable_index[j]);
		//	}
		//	temp_v[j] /=  m_ele->GetNodesNumber(false);
		//	v[j] = temp_v[j];
		//}

		// BG 04/2011: MassFlux Calculation is not working correctly if it is a 2D mesh in x-z-direction
		// z-velocity is stored a y-position -> this is corrected here
		if ((dimension == 2) && (axis == 2))
		{
			v[2] = v[1];
			v[1] = 0;
		}
		// edge projection // edge marked
		m_ele->GetEdges(ele_edges_vector);
		// cout << "Element: " << "\n";
		edg_length = 0;
		// loop over the edges of the element to find the edge at the polyline
		for (size_t j = 0; j < static_cast<size_t>(m_ele->GetEdgesNumber()); j++)
		{
			m_edg = ele_edges_vector[j];
			// check if edge was already used
			Edge_already_used = false;
			for (size_t k = 0; k < vecConsideredEdges.size(); k++)
				if (m_edg->GetIndex() == vecConsideredEdges[k])
					Edge_already_used = true;
			if (Edge_already_used == true)
			{
				Edge_on_Geo = false;
				Use_Element = false;
			}
			else
			{
				vec_nodes_edge.clear();
				// BG 04/2011: check if edge is completely on the polyline
				Edge_on_Geo = true;
				vec_nodes_edge.push_back(m_edg->GetNode(0));
				vec_nodes_edge.push_back(m_edg->GetNode(1));
				// loop over the nodes of the edge to check if all of them are on the polyline
				for (int k = 0; k < int(vec_nodes_edge.size()); k++)
				{
					Point_on_Geo = false;
					for (std::size_t l = 0; l < nod_vector_at_geo.size(); l++)
					{
						if (vec_nodes_edge[k]->GetIndex() == nod_vector_at_geo[l])
						{
							Point_on_Geo = true;
							l = nod_vector_at_geo.size();
						}
					}
					// cout << "     Node: " << vec_nodes_edge[k]->GetIndex() << " " << Point_on_Geo << "\n";
					if (Point_on_Geo == false)
						Edge_on_Geo = false;
				}
				if (Edge_on_Geo == true)
					j = m_ele->GetEdgesNumber();
			}
		}
		if ((m_edg->GetMark()) && (Use_Element == true) && (Edge_on_Geo == true))
		{
			vecConsideredEdges.push_back(m_edg->GetIndex()); // all edges that were already used are stored
			m_edg->SetNormalVector(m_ele->normal_vector, edg_normal_vector);
			edg_length = m_edg->getLength();
			// cout << "Element: " << m_ele->GetIndex() << " LÃ¤nge: " << edg_length << " Normalvektor: x=" <<
			// edg_normal_vector[0] << " y=" << edg_normal_vector[1] << " z=" << edg_normal_vector[2] << "\n";
			m_edg->GetEdgeVector(edge_vector);

			vn = MSkalarprodukt(v, edg_normal_vector, 3);
			for (size_t j = 0; j < 3; j++)
				vn_vec[j] = vn * edg_normal_vector[j];

			for (size_t j = 0; j < 3; j++)
				f[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea();

			f_n_sum += MBtrgVec(f, 3);

			if (pcs_type == FiniteElement::MULTI_PHASE_FLOW) // BG, 04/2012
			{
				// velocity vector
				for (size_t j = 0; j < 3; j++)
				{
					// v_2[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx_2[j]);
					// Calculate Element velocity
					v_2[j] = fem->Get_Element_Velocity(m_ele->GetIndex(), m_pcs_flow, 1, j);
				}

				// variable_index[0] = m_pcs_flow->GetNodeValueIndex("VELOCITY_X2");
				// variable_index[1] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Y2");
				// variable_index[2] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Z2");
				//
				// for (size_t j = 0; j < 3; j++)
				//{
				//	for (size_t k = 0; k < m_ele->GetNodesNumber(false); k++)
				//	{
				//		temp_v[j] += m_pcs_flow->GetNodeValue(element_nodes[k], variable_index[j]);
				//	}
				//	temp_v[j] /=  m_ele->GetNodesNumber(false);
				//	v_2[j] = temp_v[j];
				//}

				// BG 04/2011: MassFlux Calculation is not working correctly if it is a 2D mesh in x-z-direction
				// z-velocity is stored a y-position -> this is corrected here
				if ((dimension == 2) && (axis == 2))
				{
					v_2[2] = v_2[1];
					v_2[1] = 0;
				}

				vn = MSkalarprodukt(v_2, edg_normal_vector, 3);
				for (size_t j = 0; j < 3; j++)
					vn_vec[j] = vn * edg_normal_vector[j];

				for (size_t j = 0; j < 3; j++)
					f_2[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea();

				f_n_sum_2 += MBtrgVec(f_2, 3);
			}
		}
	}
	result[0] = f_n_sum;
	result[1] = f_n_sum_2;
	ele_vector_at_geo.clear();
}

/**************************************************************************
GeoSys - Function: CalcELEMassFluxes
Task: Calculate the Mass Flux for Elements at Polylines
Return: MassFlux
Programming: 05/2011 BG
Modification:
 **************************************************************************/
void CRFProcess::CalcELEMassFluxes(const GEOLIB::Polyline* const ply, std::string const& NameofPolyline, double* result)
{
	CRFProcess* m_pcs_flow = NULL;
	std::vector<size_t> ele_vector_at_geo;
	std::vector<size_t> nod_vector_at_geo;
	vector<CNode*> vec_nodes_edge;
	vec<CEdge*> ele_edges_vector(15);
	vec<long> element_nodes(20);
	vector<size_t> vecConsideredEdges;
	vector<size_t> vecConsideredNodes;
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	CompProperties* m_cp = NULL;
	CMediumProperties* m_mat_mp = NULL;
	CNode* m_node = NULL;
	bool Edge_on_Geo, Point_on_Geo;
	double edg_normal_vector[3];
	double vn;
	double edg_length = 0.0;
	double vn_vec[3];
	double edge_vector[3];
	//	double f_n_sum = 0.0; // TF 2012-08 // removed warning unused variable
	double J_adv_sum = 0.0, J_disp_sum = 0.0, J_diff_sum = 0.0, J_sum = 0.0;
	double C_ele = 0.0; // C1 = 0, C2 = 0; // TF 2012-08 // removed warning unused variable
	//	int number_nodes_at_edge = 0; // TF 2012-08 // removed warning unused variable
	double J_adv[3], J_disp[3], J_diff[3], v[3], temp_j, j_diff[3], j_disp[3], norm_v;
	double J_adv_temp, J_disp_temp, J_diff_temp, J_temp;
	int coordinateflag, dimension = 0, axis = 0;
	bool Edge_already_used, Node_already_used;
	bool Use_Element;
	double ElementConcentration[4];
	double Dm, Dm_eff, tortuosity, porosity, alpha_l, alpha_t, Disp_xx, Disp_yy, Disp_xy;
	// double Disp_zz, Disp_xz, Disp_yz;
	int group;
	//	double totalmass; // 2012-08 TF not used
	//	double totalmassflux; // 2012-08 TF / unused
	double* ConcentrationGradient(new double[3]);
	int numberPolyline = 0;

	if (this->Tim->step_current == 0)
	{
		this->PolylinesforOutput.push_back(NameofPolyline);
		numberPolyline = this->PolylinesforOutput.size() - 1;
		this->TotalMass[numberPolyline] = 0;
	}
	else
	{
		for (int i = 0; i < int(this->PolylinesforOutput.size()); i++)
			if (NameofPolyline == this->PolylinesforOutput[i])
				numberPolyline = i;
	}

	m_msh->GetELEOnPLY(ply, ele_vector_at_geo, true);
	// BG: 04/2011 nodes are needed to provide the correct edge
	m_msh->GetNODOnPLY(ply, nod_vector_at_geo);

	// get the current flow and transport process
	if (isFlowProcess(this->getProcessType()))
	{
		m_pcs_flow = this;
	}
	else
	{
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);
	}

	// get the indices of the velocity
	int v_eidx[3];
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	// check if velocity exists
	for (size_t i = 0; i < 3; i++)
	{
		if (v_eidx[i] < 0)
		{
			std::cout << i << " " << v_eidx[i] << "Velocity output is not specified"
			          << "\n";
			exit(0);
		}
	}
	// determine the dimension and the orientation of the mesh
	coordinateflag = m_msh->GetCoordinateFlag();
	if (coordinateflag == 10)
	{
		dimension = 1;
		axis = 0;
	} // x only
	else if (coordinateflag == 11)
	{
		dimension = 1;
		axis = 1;
	} // y only
	else if (coordinateflag == 12)
	{
		dimension = 1;
		axis = 2;
	} // z only
	else if (coordinateflag == 21)
	{
		dimension = 2;
		axis = 1;
	} // x, y only
	else if (coordinateflag == 22)
	{
		dimension = 2;
		axis = 2;
	} // x, z only
	else if (coordinateflag == 32)
	{
		dimension = 3;
		axis = 3;
	} // x, y, z only

	// loop over all elements at the polyline
	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		if (m_ele->GetIndex() == 4421)
			cout << i << "\n";
		m_ele->SetNormalVector();
		m_ele->GetNodeIndeces(element_nodes);
		Use_Element = true;
		// velocity vector
		for (size_t j = 0; j < 3; j++)
			v[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[j]);

		// BG 04/2011: MassFlux Calculation is not working correctly if it is a 2D mesh in x-z-direction
		// z-velocity is stored a y-position -> this is corrected here
		if ((dimension == 2) && (axis == 2))
		{
			v[2] = v[1];
			v[1] = 0;
		}
		//--------------------------------------------------------------------
		// get properties
		// get porosity of the element
		group = m_ele->GetPatchIndex();
		m_mat_mp = mmp_vector[group];
		porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1); // CB Now provides also heterogeneous porosity, model 11
		// get diffusion coefficient
		m_cp = cp_vec[this->pcs_component_number];
		Dm = m_cp->CalcDiffusionCoefficientCP(0, 0, this);
		tortuosity = m_mat_mp->TortuosityFunction(m_ele->GetIndex(), 0, 0);
		Dm_eff = Dm * tortuosity * porosity;
		// get element concentration
		C_ele = 0.0;
		for (int j = 0; j < int(m_ele->GetNodesNumber(false)); j++)
		{
			// element concentration = average over all points
			size_t nidx = GetNodeValueIndex(pcs_primary_function_name[0]) + 1;
			ElementConcentration[j] = GetNodeValue(element_nodes[j], nidx);
			// if (ElementConcentration[j] > 0)
			C_ele += ElementConcentration[j]; // BG: 04/2011 sum was not calculated correctly
		}
		C_ele /= (double)m_ele->GetNodesNumber(false);
		// calculate the norm (Betrag) of the velocity vector
		norm_v = MBtrgVec(v, 3);
		for (size_t j = 0; j < 3; j++)
		{
			J_adv[j] = 0;
			J_disp[j] = 0;
			J_diff[j] = 0;
			j_diff[j] = 0;
			j_disp[j] = 0;
		}

		//--------------------------------------------------------------------
		// Calculate Flux depending on element dimension
		if (m_ele->GetDimension() == 2)
		{
			// edge projection // edge marked
			m_ele->GetEdges(ele_edges_vector);
			// cout << "Element: " << "\n";
			edg_length = 0;
			// loop over the edges of the element to find the edge at the polyline
			for (size_t j = 0; j < static_cast<size_t>(m_ele->GetEdgesNumber()); j++)
			{
				m_edg = ele_edges_vector[j];
				// check if edge was already used
				Edge_already_used = false;
				for (int k = 0; k < int(vecConsideredEdges.size()); k++)
					if (m_edg->GetIndex() == vecConsideredEdges[k])
						Edge_already_used = true;
				if (Edge_already_used == true)
				{
					Edge_on_Geo = false;
					Use_Element = false;
				}
				else
				{
					vec_nodes_edge.clear();
					// BG 04/2011: check if edge is completely on the polyline
					Edge_on_Geo = true;
					vec_nodes_edge.push_back(m_edg->GetNode(0));
					vec_nodes_edge.push_back(m_edg->GetNode(1));
					// loop over the nodes of the edge to check if all of them are on the polyline
					for (int k = 0; k < int(vec_nodes_edge.size()); k++)
					{
						Point_on_Geo = false;
						for (std::size_t l = 0; l < nod_vector_at_geo.size(); l++)
						{
							if (vec_nodes_edge[k]->GetIndex() == nod_vector_at_geo[l])
							{
								Point_on_Geo = true;
								l = nod_vector_at_geo.size();
							}
						}
						// cout << "     Node: " << vec_nodes_edge[k]->GetIndex() << " " << Point_on_Geo << "\n";
						if (Point_on_Geo == false)
							Edge_on_Geo = false;
					}
					if (Edge_on_Geo == true)
						j = m_ele->GetEdgesNumber();
				}
				if ((m_edg->GetMark()) && (Use_Element == true) && (Edge_on_Geo == true))
				{
					vecConsideredEdges.push_back(m_edg->GetIndex()); // all edges that were already used are stored
					m_edg->SetNormalVector(m_ele->normal_vector, edg_normal_vector);
					edg_length = m_edg->getLength();
					// cout << "Element: " << m_ele->GetIndex() << " LÃ¤nge: " << edg_length << " Normalvektor: x=" <<
					// edg_normal_vector[0] << " y=" << edg_normal_vector[1] << " z=" << edg_normal_vector[2] << "\n";
					m_edg->GetEdgeVector(edge_vector);

					// cout << i << " " << m_ele->GetIndex() << "\n";
					// calculate the velocity vector perpendicular to the edge
					vn = MSkalarprodukt(v, edg_normal_vector, 3); // vn = MSkalarprodukt(v, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						vn_vec[j] = vn * edg_normal_vector[j];
					}

					//------------------------------------------------------------------------------
					// calculate advective mass flux = v_n * l^e * z^e * C^e
					for (size_t j = 0; j < 3; j++) // for (size_t j = 0; j < 3; j++)
					{
						// Darcy velocity times concentration times cross-section area
						J_adv[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea() * C_ele; // unit: [mol/s]
					}
					// cout << " element: " << m_ele->GetIndex() << " vx: " << vn_vec[0] << " vy: " << vn_vec[1] << "
					// vz: " << vn_vec[2] << " c1: " << ElementConcentration[0] << " c2: " << ElementConcentration[1] <<
					// " c3: " << ElementConcentration[2] << " c: " << C_ele << "\n";
					// cout << "      J_adv_x: " << J_adv[0] << " J_adv_y: " << J_adv[1] << " J_adv_z: " << J_adv[2] <<
					// "\n";

					//------------------------------------------------------------------------------
					// diffusive mass flux
					// calculate concentration gradient for the element
					Calc2DElementGradient(m_ele, ElementConcentration, ConcentrationGradient);
					// calculate diffusive mass flux
					for (size_t j = 0; j < 3; j++)
						j_diff[j] = Dm_eff * ConcentrationGradient[j]; // [mol/m2.s]
					// calculated the mass flux perpendicular to the edge
					temp_j = MSkalarprodukt(j_diff, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						j_diff[j] = temp_j * edg_normal_vector[j];
						J_diff[j] = j_diff[j] * edg_length * m_ele->GetFluxArea(); // [mol/s]
					}
					// cout << "      Dm: " << Dm_eff << " dc_x: " << ConcentrationGradient[0] << " dc_y: " <<
					// ConcentrationGradient[1] << " J_diff_x: " << J_diff[0] << " J_diff_y: " << J_diff[1] << "
					// J_diff_z: " << J_diff[2] << "\n";

					//------------------------------------------------------------------------------
					// dispersive mass flux
					// get dispersivities
					alpha_l = m_mat_mp->mass_dispersion_longitudinal;
					alpha_t = m_mat_mp->mass_dispersion_transverse;

					if (norm_v != 0)
					{
						Disp_xx = alpha_l * pow(v[0], 2) / norm_v + alpha_t * pow(v[1], 2) / norm_v; // ToDo:
						                                                                             // z-Dimension
						Disp_yy = alpha_l * pow(v[1], 2) / norm_v + alpha_t * pow(v[0], 2) / norm_v; // ToDo:
						                                                                             // z-Dimension
						Disp_xy = (alpha_l - alpha_t) * v[0] * v[1] / norm_v; // ToDo: z-Dimension
						// calculate dispersive mass flux
						j_disp[0] = porosity * (Disp_xx * ConcentrationGradient[0]
						                        + Disp_xy * ConcentrationGradient[1]); // [mol/m2.s]
						j_disp[1] = porosity * (Disp_yy * ConcentrationGradient[1]
						                        + Disp_xy * ConcentrationGradient[0]); // [mol/m2.s]
						j_disp[2] = 0; // [mol/m2.s]
					}
					else
						for (size_t j = 0; j < 3; j++)
							j_disp[j] = 0;
					// calculated the mass flux perpendicular to the edge
					temp_j = MSkalarprodukt(j_disp, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						j_disp[j] = temp_j * edg_normal_vector[j];
						J_disp[j] = j_disp[j] * edg_length * m_ele->GetFluxArea(); // [mol/s]
					}
				}
			}
		}
		if (m_ele->GetDimension() == 1)
		{
			// loop over the nodes of the element to check if the node at the polyline was already used, BG 11/2011
			Use_Element = true;
			for (size_t j = 0; j < static_cast<size_t>(m_ele->GetNodesNumber(false)); j++)
			{
				m_node = m_ele->GetNode(j);
				// check if edge was already used
				Node_already_used = false;
				for (int k = 0; k < int(vecConsideredNodes.size()); k++)
					if (m_node->GetIndex() == vecConsideredNodes[k])
						Node_already_used = true;
				if (Node_already_used == true)
				{
					Use_Element = false;
				}
			}

			if (Use_Element == true)
			{
				for (int l = 0; l < int(nod_vector_at_geo.size()); l++)
				{
					if (static_cast<size_t>(m_ele->GetNodeIndex(0)) == nod_vector_at_geo[l])
					{
						vecConsideredNodes.push_back(m_ele->GetNodeIndex(0));
						l = int(nod_vector_at_geo.size());
					}
					if (static_cast<size_t>(m_ele->GetNodeIndex(1)) == nod_vector_at_geo[l])
					{
						vecConsideredNodes.push_back(m_ele->GetNodeIndex(1));
						l = int(nod_vector_at_geo.size());
					}
				}
				// calculate advective mass flux = v_n * l^e * z^e * C^e
				// Darcy velocity times concentration times cross-section area
				J_adv[0] = norm_v * m_ele->GetFluxArea() * C_ele; // unit: [mol/s]
				// cout << " element: " << m_ele->GetIndex() << " vx: " << vn_vec[0] << " vy: " << vn_vec[1] << " vz: "
				// << vn_vec[2] << " c1: " << ElementConcentration[0] << " c2: " << ElementConcentration[1] << " c3: "
				// << ElementConcentration[2] << " c: " << C_ele << "\n";
				// cout << "      J_adv_x: " << J_adv[0] << " J_adv_y: " << J_adv[1] << " J_adv_z: " << J_adv[2] <<
				// "\n";

				//------------------------------------------------------------------------------
				// diffusive mass flux
				// calculate concentration gradient for the element
				Calc2DElementGradient(m_ele, ElementConcentration, ConcentrationGradient);
				// calculate diffusive mass flux
				j_diff[0] = Dm_eff * ConcentrationGradient[0]; // [mol/m2.s]
				J_diff[0] = j_diff[0] * m_ele->GetFluxArea(); // [mol/s]
				// cout << "      Dm: " << Dm_eff << " dc_x: " << ConcentrationGradient[0] << " dc_y: " <<
				// ConcentrationGradient[1] << " J_diff_x: " << J_diff[0] << " J_diff_y: " << J_diff[1] << " J_diff_z: "
				// << J_diff[2] << "\n";

				//------------------------------------------------------------------------------
				// dispersive mass flux
				// get dispersivity along the element direction
				alpha_l = m_mat_mp->mass_dispersion_longitudinal;

				// calculate dispersion coefficients based on velocity of the element
				if (norm_v != 0)
				{
					Disp_xx = alpha_l * norm_v;
					// calculate dispersive mass flux
					j_disp[0] = porosity * Disp_xx * ConcentrationGradient[0]; // [mol/m2.s]
				}
				else
					j_disp[0] = 0;
				// calculated the mass flux perpendicular to the edge
				J_disp[0] = j_disp[0] * m_ele->GetFluxArea(); // [mol/s]
			}
		}
		//--------------------------------------------------------------------
		// sum the mass flux over all considered elements
		J_adv_temp = MBtrgVec(J_adv, 3);
		J_disp_temp = MBtrgVec(J_disp, 3);
		J_diff_temp = MBtrgVec(J_diff, 3);
		J_temp = J_adv_temp + J_disp_temp + J_diff_temp;
		// mass check
		//		totalmassflux = 0.0; // 2012-08 TF not used
		//		totalmass = C_ele * porosity * m_ele->GetVolume(); // 2012-08 TF not used

		// 2012-08 TF, since totalmassflux is not used we can comment out the following two lines
		//		if (this->Tim->time_step_length > 0)
		//			totalmassflux = totalmass / this->Tim->time_step_length;
		// cout << "Element: " << i << " Masse: " << totalmass << " Abstrom_adv: " << this->Tim->this_stepsize *
		// MBtrgVec(J_adv, 3)  << " Abstrom_diff: " << this->Tim->this_stepsize * MBtrgVec(J_diff, 3) << " Abstrom_disp:
		// " << this->Tim->this_stepsize * MBtrgVec(J_disp, 3);
		// cout << "\n";
		// correct mass flux
		// if (J_temp > totalmassflux)
		//{
		//	J_adv_temp += J_adv_temp / J_temp * (totalmassflux - J_temp);
		//	J_disp_temp += J_disp_temp / J_temp * (totalmassflux - J_temp);
		//	J_diff_temp += J_diff_temp / J_temp * (totalmassflux - J_temp);
		//	J_temp = totalmassflux;
		//}
		J_adv_sum += J_adv_temp;
		J_disp_sum += J_disp_temp;
		J_diff_sum += J_diff_temp;
		J_sum += J_temp;
		if (this->Tim->time_step_length > 0)
			this->TotalMass[numberPolyline] += J_temp * this->Tim->time_step_length;
	}
	result[0] = J_adv_sum;
	result[1] = J_disp_sum;
	result[2] = J_diff_sum;
	result[3] = J_sum;
	result[4] = this->TotalMass[numberPolyline];
	ele_vector_at_geo.clear();

	delete[] ConcentrationGradient;
}

/**************************************************************************
GeoSys - Function: Calc2DElementGradient
Task: Calculate the Gradient for an 2D Element
Return: nothing
Programming: 05/2011 BG
Modification:
 **************************************************************************/
void CRFProcess::Calc2DElementGradient(MeshLib::CElem* m_ele, double ElementConcentration[4], double* grad)
{
	double coord_Point1[3];
	double coord_Point2[3];
	double coord_Point3[3];
	double coord_Point4[3];
	double* normal_vector;

	// get the nodes of the 2D Element
	vector<CNode*> vecElementNodes;
	m_ele->GetNodes(vecElementNodes);
	CPlaneEquation* PlaneEquation;
	PlaneEquation = new CPlaneEquation();

	if (m_ele->GetDimension() == 1)
	{
		if (m_ele->GetElementType() == MshElemType::LINE)
		{
			// calculate gradient based on the concentration difference and the length of the element
			grad[0] = (ElementConcentration[1] - ElementConcentration[0]) / (m_ele->GetVolume() / m_ele->GetFluxArea());
			grad[1] = 0;
			grad[2] = 0;
		}
		else
			cout << "This element option is not yet considered for calculating the concentration gradient!"
			     << "\n";
	}

	else if (m_ele->GetDimension() == 2)
	{
		if (m_ele->GetElementType()
		    == MshElemType::QUAD) // BG 04/2011: calculation of the normal vector of a quad element
		{
			// define the points used for the plane equation
			// order of points: 1, 2, 3, 4 against clock direction -> points are given in the order 1, 2, 4 to get a
			// positive normal vector
			double const* const p0(vecElementNodes[0]->getData());
			double const* const p1(vecElementNodes[1]->getData());
			double const* const p2(vecElementNodes[2]->getData());
			double const* const p3(vecElementNodes[3]->getData());

			coord_Point1[0] = p0[0];
			coord_Point1[1] = p0[1];
			coord_Point1[2] = ElementConcentration[0];
			coord_Point2[0] = p1[0];
			coord_Point2[1] = p1[1];
			coord_Point2[2] = ElementConcentration[1];
			coord_Point3[0] = p2[0];
			coord_Point3[1] = p2[1];
			coord_Point3[2] = ElementConcentration[2];
			coord_Point4[0] = p3[0];
			coord_Point4[1] = p3[1];
			coord_Point4[2] = ElementConcentration[3];

			// Calculate the plane equation
			PlaneEquation->CalculatePlaneEquationFrom3Points(coord_Point1, coord_Point2, coord_Point4);
			// check if 4. point lies within the plane
			if (PlaneEquation->CheckIfPointInPlane(coord_Point3) == false)
				return;
		}

		if (m_ele->GetElementType() == MshElemType::TRIANGLE)
		{
			// define the points used for the plane equation
			// order of points: 1, 2, 3 against clock direction -> points are given in the order 2, 3, 1 to get a
			// positive normal vector
			double const* const p0(vecElementNodes[0]->getData());
			double const* const p1(vecElementNodes[1]->getData());
			double const* const p2(vecElementNodes[2]->getData());

			coord_Point1[0] = p0[0];
			coord_Point1[1] = p0[1];
			coord_Point1[2] = ElementConcentration[0];
			coord_Point2[0] = p1[0];
			coord_Point2[1] = p1[1];
			coord_Point2[2] = ElementConcentration[1];
			coord_Point3[0] = p2[0];
			coord_Point3[1] = p2[1];
			coord_Point3[2] = ElementConcentration[2];

			// Calculate the plane equation
			PlaneEquation->CalculatePlaneEquationFrom3Points(coord_Point2, coord_Point3, coord_Point1);
		}
		else
			cout << "This element option is not yet considered for calculating the concentration gradient!"
			     << "\n";

		// calculate gradient
		normal_vector = PlaneEquation->GetNormalVector();
		grad[0] = normal_vector[0] / normal_vector[2];
		grad[1] = normal_vector[1] / normal_vector[2];
		grad[2] = 0; // grad[2] = normal_vector[2] / normal_vector[2];
	}

	else
		cout << "This element option is not yet considered for calculating the concentration gradient!"
		     << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2004 OK Implementation
   08/2006 OK new
**************************************************************************/
void CRFProcess::CalcELEVelocities(void)
{
	int eidx[3];

	// If not FLUID_MOMENTUM,
	//	if (_pcs_type_name.compare("RANDOM_WALK") != 0) {
	FiniteElement::ProcessType pcs_type(getProcessType()); // BG

	if (this->getProcessType() != FiniteElement::RANDOM_WALK)
	{
		int eidx[3];
		eidx[0] = GetElementValueIndex("VELOCITY1_X");
		eidx[1] = GetElementValueIndex("VELOCITY1_Y");
		eidx[2] = GetElementValueIndex("VELOCITY1_Z");
		for (size_t i = 0; i < 3; i++)
			if (eidx[i] < 0)
				cout << "Fatal error in CRFProcess::CalcELEVelocities - abort"
				     << "\n";
		// abort();	// PCH commented abort() out for FM.

		FiniteElement::ElementValue* gp_ele = NULL;
		double vx, vy, vz;
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
		{
			gp_ele = ele_gp_value[i];
			vx = gp_ele->Velocity(0, 0);
			SetElementValue(i, eidx[0], vx);
			SetElementValue(i, eidx[0] + 1, vx);
			vy = gp_ele->Velocity(1, 0);
			SetElementValue(i, eidx[1], vy);
			SetElementValue(i, eidx[1] + 1, vy);
			vz = gp_ele->Velocity(2, 0);
			SetElementValue(i, eidx[2], vz);
			SetElementValue(i, eidx[2] + 1, vz);
		}
	}
	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
	{
		eidx[0] = GetElementValueIndex("VELOCITY2_X");
		eidx[1] = GetElementValueIndex("VELOCITY2_Y");
		eidx[2] = GetElementValueIndex("VELOCITY2_Z");
		for (size_t i = 0; i < 3; i++)
			if (eidx[i] < 0)
				cout << "Fatal error in CRFProcess::CalcELEVelocities - abort"
				     << "\n";
		// abort();	// PCH commented abort() out for FM.

		FiniteElement::ElementValue* gp_ele = NULL;
		double vx, vy, vz;
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
		{
			gp_ele = ele_gp_value[i];
			vx = gp_ele->Velocity_g(0, 0);
			SetElementValue(i, eidx[0], vx);
			SetElementValue(i, eidx[0] + 1, vx);
			vy = gp_ele->Velocity_g(1, 0);
			SetElementValue(i, eidx[1], vy);
			SetElementValue(i, eidx[1] + 1, vy);
			vz = gp_ele->Velocity_g(2, 0);
			SetElementValue(i, eidx[2], vz);
			SetElementValue(i, eidx[2] + 1, vz);
		}
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   08/2006 OK Implementation
 **************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
//(vector<long>&ele_number_vector)
void CRFProcess::AssembleParabolicEquationRHSVector(CNode* m_nod)
{
	// cout << "CRFProcess::AssembleParabolicEquationRHSVector" << "\n";
	// int i;
	// WW long ldummy;
	// WW double ddummy;
	//----------------------------------------------------------------------
	// Init
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
		eqs->b[m_nod->getConnectedElementIDs()[i]] = 0.0;
	//----------------------------------------------------------------------
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	double edg_normal_vector[3];
	double edge_mid_point[3];
	vec<CEdge*> ele_edges_vector(15);
	int j;
	double aux_vector[3];
	double check_sign;
	//----------------------------------------------------------------------
	// Element velocity
	int v_eidx[3];
	// kg44  v_eidx[0] = GetElementValueIndex("VELOCITY1_X");
	// kg44  v_eidx[1] = GetElementValueIndex("VELOCITY1_Y");
	// kg44  v_eidx[2] = GetElementValueIndex("VELOCITY1_Z");
	CRFProcess* m_pcs_flow = NULL;
	//  if(_pcs_type_name.find("FLOW")!=string::npos) {
	if (isFlowProcess(this->getProcessType()))
		m_pcs_flow = this;
	else
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			cout << v_eidx[i] << i << " Warning in CRFProcess::AssembleParabolicEquationRHSVector - no PCS-VEL data"
			     << "\n";
			return;
		}
	double v[3];
	//======================================================================
	// Topology
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
	{
		m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
		m_ele->SetNormalVector(); // OK_BUGFIX
		v[0] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[0]);
		v[1] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[1]);
		v[2] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[2]);
		m_ele->SetMark(false);
		switch (m_ele->GetElementType())
		{
			//------------------------------------------------------------------
			// line elements
			case MshElemType::LINE:
			{
				v[1] = GetElementValue(m_ele->GetIndex(), v_eidx[0]);
				v[0] = GetElementValue(m_ele->GetIndex(), v_eidx[1]);
				if (m_nod->getConnectedElementIDs().size() == 1)
				{
					m_ele->SetMark(true);
					break;
				}
				double const* gravity_center(m_ele->GetGravityCenter());
				double const* const pnt(m_nod->getData());
				aux_vector[0] = gravity_center[0] - pnt[0];
				aux_vector[1] = gravity_center[1] - pnt[1];
				aux_vector[2] = gravity_center[2] - pnt[2];
				check_sign = MSkalarprodukt(v, aux_vector, 3);
				if (check_sign < 0.0)
					m_ele->SetMark(true);
				break;
			}
			//------------------------------------------------------------------
			// tri elements
			case MshElemType::TRIANGLE:
				m_ele->GetEdges(ele_edges_vector);
				for (j = 0; j < (int)m_ele->GetEdgesNumber(); j++)
				{
					m_edg = ele_edges_vector[j];
					if (m_edg->GetMark())
					{
						m_edg->SetNormalVector(m_ele->normal_vector, edg_normal_vector);
						break;
						/*
						            m_edg->GetEdgeMidPoint(edge_mid_point);
						           gravity_center = m_ele->GetGravityCenter();
						           aux_vector[0] = gravity_center[0] - edge_mid_point[0];
						           aux_vector[1] = gravity_center[1] - edge_mid_point[1];
						           aux_vector[2] = gravity_center[2] - edge_mid_point[2];
						           check_sign = MSkalarprodukt(edg_normal_vector,aux_vector,3);
						           if(check_sign<0.0) break;
						 */
					}
				}
			//----------------------------------------------------------------
			// ToDo
				//if (m_edg->GetMark())
				//	break;
                                break;
			default:
				cout << "Warning in CRFProcess::AssembleParabolicEquationRHSVector - not implemented for this element "
				        "type"
				     << "\n";
				break;
		} // switch
	}
	//======================================================================
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
	{
		m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
		switch (m_ele->GetElementType())
		{
			//------------------------------------------------------------------
			// line elements
			case MshElemType::LINE:
				if (m_ele->GetMark())
				{
					cout << m_ele->GetIndex() << "\n";
					// WW ldummy = m_nod->GetIndex();
					// WW ddummy = eqs->b[m_nod->GetIndex()];
					fem->ConfigElement(m_ele, false);
					fem->AssembleParabolicEquationRHSVector();
					// WW ddummy = eqs->b[m_nod->GetIndex()];
				}
				break;
			//------------------------------------------------------------------
			// tri elements
			case MshElemType::TRIANGLE:
			{
				m_edg->GetEdgeMidPoint(edge_mid_point);
				double const* gravity_center(m_ele->GetGravityCenter());
				aux_vector[0] = gravity_center[0] - edge_mid_point[0];
				aux_vector[1] = gravity_center[1] - edge_mid_point[1];
				aux_vector[2] = gravity_center[2] - edge_mid_point[2];
				check_sign = MSkalarprodukt(edg_normal_vector, aux_vector, 3);
				if (check_sign < 0.0)
					continue;
				{
					// cout << m_ele->GetIndex() << "\n";
					fem->ConfigElement(m_ele, false);
					fem->AssembleParabolicEquationRHSVector();
				}
				break;
			}
			//----------------------------------------------------------------
			// ToDo
			default:
				cout << "Warning in CRFProcess::AssembleParabolicEquationRHSVector - not implemented for this element "
				        "type"
				     << "\n";
				break;
		} // switch
	}
	//======================================================================
}
#endif

/**************************************************************************
PCSLib-Method:
08/2006 OK Implementation
compare with CMCDs PCSGetFluxProcess
03/2012 JT No more loops.
**************************************************************************/
CRFProcess* PCSGetFlow()
{
	if (pcs_number_flow >= 0)
		return pcs_vector[pcs_number_flow];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetHeat()
{
	if (pcs_number_heat >= 0)
		return pcs_vector[pcs_number_heat];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetDeformation()
{
	if (pcs_number_deformation >= 0)
		return pcs_vector[pcs_number_deformation];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetMass(size_t component_number)
{
	if (component_number < DOF_NUMBER_MAX)
	{ // don't exceed array dimensions
		if (pcs_number_mass[component_number] >= 0)
		{
			return pcs_vector[pcs_number_mass[component_number]];
		}
	}
	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   based on MMPCalcSecondaryVariables
   01/2007 OK Implementation
**************************************************************************/
void CRFProcess::SetBC()
{
	// WW CBoundaryCondition *m_bc = NULL;
	CBoundaryConditionNode* m_node = NULL;
	int nidx = GetNodeValueIndex(pcs_primary_function_name[0]);
	for (long i = 0; i < (long)bc_node_value.size(); i++)
	{
		m_node = bc_node_value[i];
		// WW m_bc = bc_node[i];
		// old time
		SetNodeValue(m_node->msh_node_number, nidx, m_node->node_value);
		// new time
		SetNodeValue(m_node->msh_node_number, nidx + 1, m_node->node_value);
	}
}

/**************************************************************************
   Task: Preprocessing function set flag pcs->flow_pcs_type for
      calculation of velocity and Saturation in mass transport element
      matrices

   Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void SetFlowProcessType()
{
	int i;
	int no_processes, flowtype;
	CRFProcess* m_pcs = NULL;

	m_pcs = PCSGetFlow();
	flowtype = m_pcs->type;
	no_processes = (int)pcs_vector.size();

	for (i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		// if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0)
		// m_pcs->twophaseflow=true;
		m_pcs->flow_pcs_type = flowtype;
	}
}

/**************************************************************************
   Task: Postprocessing function copies the new time step node values of
      secondary variables PRESSURE1 and SATURATION1 to the old time level

   Programing:
   13/2008   CB   Implementation                                          */
/**************************************************************************/
void CopyTimestepNODValuesSVTPhF()
{
	long i, j;
	int idx0, idx1;
	long nnodes;
	double val = 0;
	CRFProcess* m_pcs = NULL;

	nnodes = (long)fem_msh_vector[0]->nod_vector.size();

	for (j = 0; j < 2; j++) // pcs 1 and 2
	{
		m_pcs = pcs_vector[j];
		//		if (m_pcs->_pcs_type_name.compare("TWO_PHASE_FLOW") != 0)
		if (m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
			break;
		if (j == 0)
			// old timelevel
			idx0 = m_pcs->GetNodeValueIndex("PRESSURE2");
		else
			// old timelevel
			idx0 = m_pcs->GetNodeValueIndex("SATURATION1");
		idx1 = idx0 + 1;
		for (i = 0; i < nnodes; i++)
		{
			val = m_pcs->GetNodeValue(i, idx1);
			m_pcs->SetNodeValue(i, idx0, val);
		}
		// m_pcs->WriteAllVariables();
	}
}

/**************************************************************************
   PCSLib-Method:
   based on WriteSolution by WW
   01/2007 OK Implementation
**************************************************************************/
void CRFProcess::WriteAllVariables()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	string m_file_name = FileName + "_" + pcs_type_name + "_" + pcs_primary_function_name[0] + ".asc";
	ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}
	//
	int j;
	int idx[20];
	for (j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		os << pcs_primary_function_name[j] << " ";
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		os << pcs_primary_function_name[j] << " ";
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	if (this->type == 12) // output of old & new time step for secondary variables PRESSURE2 and SATURATION1
	{
		os << pcs_secondary_function_name[0] << " ";
		idx[2 * pcs_number_of_primary_nvals + 0] = GetNodeValueIndex(pcs_secondary_function_name[0]);
		os << pcs_secondary_function_name[0] << " ";
		idx[2 * pcs_number_of_primary_nvals + 1] = GetNodeValueIndex(pcs_secondary_function_name[0]) + 1;
		// other secondary variables
		for (j = 2; j < pcs_number_of_secondary_nvals; j++)
		{
			os << pcs_secondary_function_name[j] << " ";
			idx[2 * pcs_number_of_primary_nvals + j] = GetNodeValueIndex(pcs_secondary_function_name[j]);
		}
	}
	else
		for (j = 0; j < pcs_number_of_secondary_nvals; j++)
		{
			os << pcs_secondary_function_name[j] << " ";
			idx[2 * pcs_number_of_primary_nvals + j] = GetNodeValueIndex(pcs_secondary_function_name[j]);
		}
	os << "\n";
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			os << GetNodeValue(i, idx[j]) << "  ";
		for (j = 0; j < pcs_number_of_secondary_nvals; j++)
			os << GetNodeValue(i, idx[2 * pcs_number_of_primary_nvals + j]) << "  ";
		os << "\n";
	}
	os.close();
}

/**************************************************************************
   PCSLib-Method:
   based on MMPCalcSecondaryVariables
   01/2007 OK Implementation
   08/2008 CB NAPLdissolution
**************************************************************************/
void MMPCalcSecondaryVariablesNew(CRFProcess* m_pcs, bool NAPLdiss)
{
	long i;

	//----------------------------------------------------------------------
	int ndx_density_phase;
	int ndx_viscosity_phase;
	ndx_density_phase = -1; // WW
	ndx_viscosity_phase = -1; // WW
	CFEMesh* m_msh = m_pcs->m_msh; // PCH
	//----------------------------------------------------------------------
	m_pcs->SetBC();

	// For accessing the other process
	CRFProcess* cpl_pcs = NULL;
	if (m_pcs->pcs_type_number == 0)
		cpl_pcs = pcs_vector[m_pcs->pcs_number + 1];
	else if (m_pcs->pcs_type_number == 1)
		cpl_pcs = pcs_vector[m_pcs->pcs_number - 1];

	int ndx_pressure1, ndx_p_cap, ndx_pressure2, ndx_s_wetting, ndx_s_nonwetting;
	//======================================================================
	//----------------------------------------------------------------------
	// Capillary pressure - p_c (S) <- This is always the secondary variable
	// in both phase1 and phase2	// PCH
	CMediumProperties* mmp = NULL;

	//======================================================================
	switch (m_pcs->pcs_type_number)
	{
		case 0:
			//..................................................................
			//..................................................................
			// PCH
			// The primary variable is PRESSURE1
			// From PRESSURE1, we are assigning PRESSURE2 which is
			// the secondary variables of PRESSURE1.
			ndx_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1");
			ndx_pressure2 = m_pcs->GetNodeValueIndex("PRESSURE2");
			ndx_p_cap = m_pcs->GetNodeValueIndex("PRESSURE_CAP");
			double pressure1, pressure2, p_cap;
			for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
			{
				// New
				pressure1 = m_pcs->GetNodeValue(i, ndx_pressure1 + 1);
				// Old
				pressure2 = m_pcs->GetNodeValue(i, ndx_pressure1);

				// Let's get capillary pressure before updating pressure2
				// by accessing the primary variable of the saturation equation
				// not the secondary variable of it.
				int cpl_ndx_sat2 = cpl_pcs->GetNodeValueIndex("SATURATION2");
				double cpl_sat2 = cpl_pcs->GetNodeValue(i, cpl_ndx_sat2 + 1);

				if (mmp_vector.size() > 1)
				{
					double sum = 0.0;
					CNode* thisNode = m_msh->nod_vector[i];
					int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();
					// Harmonic mean
					for (int p = 0; p < NumOfNeighborElements; ++p)
					{
						// Mount neighboring elemenets and get the corresponding material group one by one.
						int eleIdx = thisNode->getConnectedElementIDs()[p];
						CElem* thisEle = m_msh->ele_vector[eleIdx];
						int matgrp = thisEle->GetPatchIndex();
						mmp = mmp_vector[matgrp];
						mmp->mode = 2;
						sum += 1.0 / cpl_sat2;
					}
					cpl_sat2 = (double)NumOfNeighborElements / sum;
				}
				// Assigning the secondary variable, Pc
				if (mmp_vector.size() > 1)
					p_cap = m_pcs->GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0 - cpl_sat2);
				else
					p_cap = mmp->CapillaryPressureFunction(1.0 - cpl_sat2);

				m_pcs->SetNodeValue(i, ndx_p_cap, p_cap);
				m_pcs->SetNodeValue(i, ndx_p_cap + 1, p_cap);

				pressure2 = pressure1 + p_cap;
				// Assigning the secondary variables
				// Previous
				m_pcs->SetNodeValue(i, ndx_pressure2, pressure2);
				// Now
				m_pcs->SetNodeValue(i, ndx_pressure2 + 1, pressure2);
			}
			//......................................................................
			ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY1");
			ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY1");
			printf("Pressure2 from the known Pressure1 is updated for Process 0\n");
			break;

		case 1:
			//..................................................................
			// PCH
			// Calc secondary variable saturation Snonwetting = 1-Swetting
			// Don't forget here the primary variable is SATURATION2
			// From SATURATION2, we are assigning SATURATION1 which is
			// the secondary variables of SATURATION2.
			ndx_s_wetting = m_pcs->GetNodeValueIndex("SATURATION1");
			ndx_s_nonwetting = m_pcs->GetNodeValueIndex("SATURATION2");
			ndx_p_cap = cpl_pcs->GetNodeValueIndex("PRESSURE_CAP");

			double s_wetting, s_nonwetting;
			for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
			{
				s_nonwetting = m_pcs->GetNodeValue(i, ndx_s_nonwetting + 1);
				// Due to the iterative solution scheme in solving Snw with no
				// explicit boundary condition for non-zero flux condition,
				// Snw may become negative particularly the density difference
				// between two fluids is big. To prevent negative Snw, the
				// saturation restriction added.
				if (mmp_vector.size() > 1)
				{
					double sum = 0.0;
					CNode* thisNode = m_msh->nod_vector[i];
					int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();
					// Harmonic mean
					for (int p = 0; p < NumOfNeighborElements; ++p)
					{
						// Mount neighboring elemenets and get the corresponding material group one by one.
						int eleIdx = thisNode->getConnectedElementIDs()[p];
						CElem* thisEle = m_msh->ele_vector[eleIdx];
						int matgrp = thisEle->GetPatchIndex();
						mmp = mmp_vector[matgrp];
						mmp->mode = 2;
						sum += 1.0 / s_nonwetting;
					}
					s_nonwetting = (double)NumOfNeighborElements / sum;
				}
				// Assigning the secondary variable, Pc
				if (mmp_vector.size() > 1)
					p_cap = m_pcs->GetCapillaryPressureOnNodeByNeighobringElementPatches(i, 2, 1.0 - s_nonwetting);
				else
					p_cap = mmp->CapillaryPressureFunction(1.0 - s_nonwetting);

				m_pcs->SetNodeValue(i, ndx_s_nonwetting, s_nonwetting);
				m_pcs->SetNodeValue(i, ndx_s_nonwetting + 1, s_nonwetting);
				s_wetting = 1.0 - s_nonwetting;

				// Assigning the secondary variables
				// Previous
				m_pcs->SetNodeValue(i, ndx_s_wetting, s_wetting);
				// Now
				m_pcs->SetNodeValue(i, ndx_s_wetting + 1, s_wetting);

				cpl_pcs->SetNodeValue(i, ndx_p_cap, p_cap);
				cpl_pcs->SetNodeValue(i, ndx_p_cap + 1, p_cap);
			}
			//......................................................................
			ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY2");
			ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY2");
			printf("Saturation1 from the known Saturation2 is updated for Process 1\n");
			break;
	}

	//----------------------------------------------------------------------
	// Fluid properties
	double density;
	double viscosity;
	CFluidProperties* m_mfp = NULL;
	m_mfp = mfp_vector[m_pcs->pcs_type_number];
	m_mfp->mode = 1;
	for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
	{
		// CB NAPL dissolution reqiuires update of Density based on new composition of NAPL phase
		// CB phase 2
		if (NAPLdiss == true && m_pcs->pcs_type_number == 1)
			density = CalcNAPLDens(i);
		else
			density = m_mfp->Density();
		m_pcs->SetNodeValue(i, ndx_density_phase, density);
		viscosity = m_mfp->Viscosity();
		m_pcs->SetNodeValue(i, ndx_viscosity_phase, viscosity);
	}
	m_mfp->mode = 0;
	//----------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::OBJRelations()
{
	bool succeed = true;
	std::cout << "OBJ->PCS relations" << '\n';

	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));

	// NUM
	std::cout << " - NUM->PCS" << '\n';
	m_num = NUMGet(pcs_type_name);
	if (!m_num)
	{
		std::cout << "Warning in CRFProcess::Create() - no NUM data - default"
		          << "\n";
		succeed = false;
	}
	else
	{
		num_type_name = "NEW"; // OK
	}

	// TIM
	cout << " - TIM->PCS" << '\n';
	Tim = TIMGet(pcs_type_name);
	if (Tim)
	{
		// Time unit factor //WW OK: -> TIM
		if (Tim->time_unit.find("MINUTE") != string::npos)
			time_unit_factor = 60.0;
		else if (Tim->time_unit.find("HOUR") != string::npos)
			time_unit_factor = 3600.0;
		else if (Tim->time_unit.find("DAY") != string::npos)
			time_unit_factor = 86400.0;
		else if (Tim->time_unit.find("MONTH") != string::npos)
			time_unit_factor = 2592000.0;
		else if (Tim->time_unit.find("YEAR") != string::npos)
			time_unit_factor = 31536000;
	}
	else
		cout << "Warning in CRFProcess::Create() - no TIM data - default"
		     << "\n";

	// OUT

	// MSH
	if (fem_msh_vector.size() == 1)
		m_msh = fem_msh_vector[0];
	else
		m_msh = MSHGet(pcs_type_name);
	if (!m_msh)
	{
		cout << "Warning in CRFProcess::Create() - no MSH data"
		     << "\n";
		succeed = false;
	}

	return succeed;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::NODRelations()
{
	bool succeed = true;
	size_t DOF = GetPrimaryVNumber(); // OK should be PCS member variable

	cout << "NOD->PCS relations" << '\n';

	// BC
	cout << " - BC->PCS" << '\n';
	if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos) // WW

		setBC_danymic_problems();
	else
	{
		// create BC groups for each process
		CBoundaryConditionsGroup* m_bc_group = NULL;
		std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
		for (size_t i = 0; i < DOF; i++)
		{
			BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
			m_bc_group = new CBoundaryConditionsGroup();
			// OK
			m_bc_group->setProcessTypeName(pcs_type_name);
			m_bc_group->setProcessPrimaryVariableName(pcs_primary_function_name[i]); // OK
			m_bc_group->Set(this, Shift[i]);
		}
	}

	// ST
	if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos) // WW

		setST_danymic_problems();

	CSourceTermGroup* m_st_group = NULL;
	if (WriteSourceNBC_RHS == 2) // Read from file
		ReadRHS_of_ST_NeumannBC();
	else // WW
	{ // Calculate directly
		std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
		for (size_t i = 0; i < DOF; i++)
		{
			m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
			if (!m_st_group)
			{
				m_st_group = new CSourceTermGroup();
				// OK
				m_st_group->pcs_type_name = pcs_type_name;
				// OK
				m_st_group->pcs_pv_name = pcs_primary_function_name[i];
				m_st_group->Set(this, Shift[i]);
			}
		}
		if (WriteSourceNBC_RHS == 1) // WW
			WriteRHS_of_ST_NeumannBC();
	}

	// NOD values
	cout << "->Config NOD values" << '\n';

	// Names
	nod_val_name_vector.clear();
	//
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		// new time
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
		// old time //need this MB!
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
	}
	for (int i = 0; i < pcs_number_of_secondary_nvals; i++)
		// new time
		nod_val_name_vector.push_back(pcs_secondary_function_name[i]);
	if ((int)nod_val_name_vector.size() != (2 * pcs_number_of_primary_nvals + pcs_number_of_secondary_nvals))
		succeed = false;

	// Values
	double* nod_values = NULL;
	const size_t nod_val_vector_size(nod_val_vector.size());
	for (size_t i = 0; i < nod_val_vector_size; i++)
	{
		delete[] nod_val_vector[i];
		nod_val_vector[i] = NULL;
	}
	nod_val_vector.clear();
	//
	// OK m_msh->NodesNumber_Quadratic;
	number_of_nvals = 2 * DOF + pcs_number_of_secondary_nvals;
	size_t nn = m_msh->nod_vector.size(); // 11.12.2012. WW
	for (long j = 0; j < number_of_nvals; j++)
	{
		nod_values = new double[nn];
		for (size_t i = 0; i < nn; i++)
			nod_values[i] = 0.0;
		nod_val_vector.push_back(nod_values);
	}
	if ((long)nod_val_vector.size() != (long)m_msh->nod_vector.size())
		succeed = false;

	// IC
	cout << "->Assign IC" << '\n';
	if (reload == 2 && type != 4 && type != 41)
		ReadSolution(); // WW
	SetIC();

	if (pcs_type_name_vector.size() && pcs_type_name_vector[0].find("DYNAMIC") != string::npos) // WW
		setIC_danymic_problems();
	return succeed;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::ELERelations()
{
	bool succeed = true;
	// OK->MB please shift to Config()
	//	if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0)
	if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		MSHDefineMobile(this);
	//
	if (type == 4 || type == 41)
		m_msh->SwitchOnQuadraticNodes(true);
	else
		m_msh->SwitchOnQuadraticNodes(false);
	CheckMarkedElement();

	// ELE - GP values
	AllocateMemGPoint();

	// ELE values
	double* ele_values = NULL; // PCH
	int number_of_evals = 2 * pcs_number_of_evals; // PCH, increase memory
	if (number_of_evals > 0) // WW added this "if" condition
	{
		// Names
		for (int i = 0; i < pcs_number_of_evals; i++)
		{
			// new time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
			// old time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
		}
		if (ele_val_name_vector.size() == 0)
			succeed = false;

		// Values
		size_t m_msh_ele_vector_size(m_msh->ele_vector.size());
		if (ele_val_vector.size() == 0)
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = new double[number_of_evals];
				size_eval += number_of_evals; // WW
				for (int i = 0; i < number_of_evals; i++)
					ele_values[i] = 0.0;
				ele_val_vector.push_back(ele_values);
			}
		else
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = ele_val_vector[j];
				/* //Comment by WW
				   #ifndef SX
				   #ifdef GCC
				   size = malloc_usable_size( ele_values )/sizeof(double);
				   #elif HORIZON
				   //KG44: malloc_usable_size and _msize are not available
				   #else
				   size= _msize( ele_values )/sizeof(double);
				   #endif
				   #endif
				 */
				ele_values = resize(ele_values, size_eval, size_eval + number_of_evals);
				size_eval += number_of_evals;
				ele_val_vector[j] = ele_values;
			}
		if (ele_val_vector.size() != m_msh->ele_vector.size())
			succeed = false;
	}

	// ELE matrices
	if (Memory_Type != 0)
	{
		AllocateLocalMatrixMemory();
		if ((long)Ele_Matrices.size() != (long)m_msh->ele_vector.size())
			succeed = false;
	}

	// Element matrix output. WW
	if (Write_Matrix)
	{
		cout << "->Write Matrix" << '\n';
		string m_file_name = FileName + "_" + convertProcessTypeToString(this->getProcessType())
		                     + "_element_matrix.txt";
		matrix_file = new fstream(m_file_name.c_str(), ios::trunc | ios::out);
		if (!matrix_file->good())
			cout << "Warning in GlobalAssembly: Matrix files are not found"
			     << "\n";
	}

	// FEM
	if (type == 4 || type == 41)
	{
		// Set initialization function
		CRFProcessDeformation* dm_pcs = (CRFProcessDeformation*)this;
		dm_pcs->Initialization();
		if (!dm_pcs->GetFEMAssembler())
			succeed = false;
	}
	else // Initialize FEM calculator
	{
		int Axisymm = 1; // ani-axisymmetry
		if (m_msh->isAxisymmetry())
			Axisymm = -1; // Axisymmetry is true
		// OK4801 needs NUM
		fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
		if (!fem)
			succeed = false;
	}

	return succeed;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
bool CRFProcess::CreateEQS()
{
	if (!m_num)
		return false; // OK46
	bool succeed = true;
	//----------------------------------------------------------------------------
	if (eqs)
		return false;
	//----------------------------------------------------------------------------
	int DOF = GetPrimaryVNumber(); // OK should be PCS member variable
	//----------------------------------------------------------------------------
	// EQS - create equation system
	cout << "->Create EQS" << '\n';
	//----------------------------------------------------------------------------
	if (type == 4)
	{
		eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF, DOF * m_msh->GetNodesNumber(true));
		InitializeLinearSolver(eqs, m_num);
		PCS_Solver.push_back(eqs); // WW
	}
	//----------------------------------------------------------------------------
	else if (type == 41)
	{
		if (num_type_name.find("EXCAVATION") != string::npos)
			eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF - 1, DOF * m_msh->GetNodesNumber(true));
		else
			eqs = CreateLinearSolverDim(m_num->ls_storage_method, DOF,
			                            (DOF - 1) * m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false));
		InitializeLinearSolver(eqs, m_num);
		PCS_Solver.push_back(eqs); // WW
	}
	//----------------------------------------------------------------------------
	else
	{
		/*
		    // If there is a solver exsiting. WW
		    CRFProcess* m_pcs = NULL;
		    for(int i=0; i<(int)pcs_vector.size(); i++)
		   {
		      m_pcs = pcs_vector[i];
		      if(m_pcs&&m_pcs->eqs)
		     {
		        if(m_pcs->_pcs_type_name.find("DEFORMATION")==string::npos)
		          break;
		     }
		   }
		   // If unique mesh
		   if(m_pcs&&m_pcs->eqs&&(fem_msh_vector.size()==1))
		   eqs = m_pcs->eqs;
		   //
		   else
		   {
		 */
		eqs = CreateLinearSolver(m_num->ls_storage_method, m_msh->GetNodesNumber(false));
		InitializeLinearSolver(eqs, m_num);
		PCS_Solver.push_back(eqs);
		//}
	}
	//----------------------------------------------------------------------------
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	strcpy(eqs->pcs_type_name, pcs_type_name.data());
	//----------------------------------------------------------------------------
	if ((int)PCS_Solver.size() == 0)
		succeed = false;
	//----------------------------------------------------------------------------
	return succeed;
}
#endif
/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 07.11.2008
void PCSCreateNew()
{
	int i;
	CRFProcess* m_pcs = NULL;
	//----------------------------------------------------------------------
	for (i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		m_pcs->CreateNew();
		//----------------------------------------------------------------------
	}
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
void CRFProcess::CreateNew()
{
	pcs_type_number = (int)pcs_vector.size();
	Config();
	m_bCheckOBJ = OBJRelations();
	m_bCheckEQS = CreateEQS();
	m_bCheckNOD = NODRelations();
	m_bCheckELE = ELERelations();
	MMP2PCSRelation(this);
	ConfigureCouplingForLocalAssemblier();
}
#endif
/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::Check()
{
	// MMP
	MSHTestMATGroups();
	return true;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool PCSCheck()
{
	if ((int)pcs_vector.size() == 0)
		return false;
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		// if(m_pcs->m_bCheck)
		if (!m_pcs->Check())
			return false;
#ifdef MFC
		FMRead(); // WW
#endif
	}

	return true;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC) // WW. 07.11.2008. 04.2012
void EQSDelete()
{
	LINEAR_SOLVER* eqs = NULL;
	//----------------------------------------------------------------------
	for (size_t i = 0; i < PCS_Solver.size(); i++)
	{
		eqs = PCS_Solver[i];
		if (eqs->unknown_vector_indeces)
			eqs->unknown_vector_indeces = (int*)Free(eqs->unknown_vector_indeces);
		if (eqs->unknown_node_numbers)
			eqs->unknown_node_numbers = (long*)Free(eqs->unknown_node_numbers);
		if (eqs->unknown_update_methods)
			eqs->unknown_update_methods = (int*)Free(eqs->unknown_update_methods);
		eqs = DestroyLinearSolver(eqs);
		PCS_Solver.erase((PCS_Solver.begin() + i));
	}
	// PCS_Solver.clear();
}
#endif
/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
void CRFProcess::NODRelationsDelete()
{
	std::string pcs_type_name(convertProcessTypeToString(getProcessType()));
	// BC
	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
	const size_t bc_node_value_size(bc_node_value.size());
	for (size_t i = 0; i < bc_node_value_size; i++)
	{
		delete bc_node_value[i];
		bc_node_value[i] = NULL;
	}
	bc_node_value.clear();

	// ST
	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		STGroupDelete(pcs_type_name, pcs_primary_function_name[i]);

	CNodeValue* nod_val = NULL;
	const size_t st_node_value_size(st_node_value.size());
	for (size_t i = 0; i < st_node_value_size; i++)
	{
		nod_val = st_node_value[i];
		// OK delete st_node_value[i];
		// OK st_node_value[i] = NULL;
		if (nod_val->check_me) // OK
		{
			nod_val->check_me = false;
			delete nod_val;
			nod_val = NULL;
		}
	}
	st_node_value.clear();

	// NOD values
	nod_val_name_vector.clear();
	const size_t nod_val_vector_size(nod_val_vector.size());
	for (size_t i = 0; i < nod_val_vector_size; i++)
	{
		delete nod_val_vector[i];
		nod_val_vector[i] = NULL;
	}
	nod_val_vector.clear();
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
void CRFProcess::ELERelationsDelete()
{
	long i;
	//----------------------------------------------------------------------
	// FEM element
	if (fem)
		delete fem; // WW
	fem = NULL;
	//----------------------------------------------------------------------
	// ELE matrices
	ElementMatrix* eleMatrix = NULL;
	ElementValue* gp_ele = NULL;
	if (Ele_Matrices.size() > 0)
	{
		for (i = 0; i < (long)Ele_Matrices.size(); i++)
		{
			eleMatrix = Ele_Matrices[i];
			delete eleMatrix;
			eleMatrix = NULL;
		}
		Ele_Matrices.clear();
	}
	//----------------------------------------------------------------------
	// ELE - GP values
	if (ele_gp_value.size() > 0)
	{
		for (i = 0; i < (long)ele_gp_value.size(); i++)
		{
			gp_ele = ele_gp_value[i];
			delete gp_ele;
			gp_ele = NULL;
		}
		ele_gp_value.clear();
	}
	//----------------------------------------------------------------------
	// ELE values
	ele_val_name_vector.clear();
	for (i = 0; i < (long)ele_val_vector.size(); i++)
	{
		delete ele_val_vector[i];
		// delete[] ele_val_vector[i];
		ele_val_vector[i] = NULL;
	}
	ele_val_vector.clear();
	//----------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
void CRFProcess::OBJRelationsDelete()
{
	//----------------------------------------------------------------------------
	cout << "OBJ->PCS relations delete" << '\n';
	//----------------------------------------------------------------------------
	m_num = NULL;
	Tim = NULL;
	m_msh = NULL;
	//----------------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
void CRFProcess::Delete()
{
	//----------------------------------------------------------------------------
	cout << "PCS  delete" << '\n';
	//----------------------------------------------------------------------------
	ELERelationsDelete();
	NODRelationsDelete();
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS                                //WW. 07.11.2008
	EQSDelete();
#endif
	OBJRelationsDelete();
	// MMP2PCSRelation(this);
	// ConfigureCouplingForLocalAssemblier();
	//----------------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW 07.11.2008
void CRFProcess::EQSDelete()
{
	std::string pcs_type_name(convertProcessTypeToString(this->getProcessType()));
	LINEAR_SOLVER* eqs = NULL;
	for (size_t i = 0; i < PCS_Solver.size(); i++)
	{
		eqs = PCS_Solver[i];
		if (pcs_type_name.compare(eqs->pcs_type_name) == 0)
		{
			if (eqs->unknown_vector_indeces)
				eqs->unknown_vector_indeces = (int*)Free(eqs->unknown_vector_indeces);
			if (eqs->unknown_node_numbers)
				eqs->unknown_node_numbers = (long*)Free(eqs->unknown_node_numbers);
			if (eqs->unknown_update_methods)
				eqs->unknown_update_methods = (int*)Free(eqs->unknown_update_methods);
			eqs = DestroyLinearSolver(eqs);
			eqs = NULL;
		}
	}

	for (size_t i = 0; i < PCS_Solver.size(); i++)
	{
		eqs = PCS_Solver[i];
		if (pcs_type_name.compare(eqs->pcs_type_name) == 0)
			PCS_Solver.erase((PCS_Solver.begin() + i));
	}
}
#endif

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   09/2007 WW Implementation
 **************************************************************************/
#if defined(USE_PETSC) // WW. 07.11.2008. 04.2012
// New solvers WW
#elif defined(NEW_EQS) // 1.09.2007 WW

void CreateEQS_LinearSolver()
{
	size_t i;
	// CB_merge_0513
	// int dof_DM = 1;
	int dof_DM = 0; // WW 02.2023. Pardiso
	CRFProcess* m_pcs = NULL;
	CFEMesh* a_msh = NULL;

#ifndef USE_MPI
	bool need_eqs = false; // WW 02.2013. Pardiso
	bool need_eqs_dof = false; // WW 02.2013. Pardiso
#endif
	int dof = 1;
	//
	// size_t dof_nonDM (1);     //WW 02.2013. Pardiso
	size_t dof_nonDM(0);

	int dof_max = 0;
	size_t dof_nonDM_max = 0;
	for (i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		if (m_pcs->type == 1212) // Important for parallel computing. 24.1.2011 WW
		{
			dof_nonDM = m_pcs->GetPrimaryVNumber();
			dof = dof_nonDM;
		}
		if (m_pcs->type == 4 || m_pcs->type / 10 == 4) // Deformation
		{
			dof_DM = m_pcs->GetPrimaryVNumber();
			dof = dof_DM;
			if (m_pcs->type == 42)
				dof = m_pcs->m_msh->GetMaxElementDim();
		}
		else // Monolithic scheme for the process with linear elements
		{
			// CB_merge_0513
			// if(dof_nonDM < m_pcs->GetPrimaryVNumber()) //WW 02.2023. Pardiso
			//{
			//   dof_nonDM = m_pcs->GetPrimaryVNumber();
			//   // PCH: DOF Handling for FLUID_MOMENTUM in case that the LIS and PARDISO solvers
			//   // are chosen.
			//   //				if(m_pcs->_pcs_type_name.compare("FLUID_MOMENTUM")==0)
			//   if(m_pcs->getProcessType() == FLUID_MOMENTUM)
			//      dof_nonDM = 1;
			//} //WW 02.2023. Pardiso

			// 02.2013. WW //WW 02.2023. Pardiso
			// Assume that the system with linear element only have one equation with DOF >1;
			if (m_pcs->GetPrimaryVNumber() > 1)
			{
				dof_nonDM = m_pcs->GetPrimaryVNumber();
				dof = dof_nonDM;
#ifndef USE_MPI
				need_eqs_dof = true;
#endif
			}
			else
			{
				dof = 1;
#ifndef USE_MPI
				need_eqs = true;
#endif
			} // WW 02.2023. Pardiso
		}
		dof_max = std::max(dof_max, dof);
		dof_nonDM_max = std::max(dof_nonDM_max, dof_nonDM);
	}
	// Check whether the JFNK method is employed for deformation problem 04.08.2010 WW
	CNumerics* num = NULL;
	for (i = 0; i < num_vector.size(); i++)
	{
		num = num_vector[i];
		if (num->nls_method == 2)
		{
			// Stiffness matrix of lower order grid may be used by more than one processes.
			// Therefore, memo allocation is performed for it.
			// Indicator: Not allocate memo for stiffness matrix for deformation
			dof_DM *= -1;
			break;
		}
	}

#ifndef USE_MPI
	SparseTable* sp = NULL;
	SparseTable* spH = NULL;
	Linear_EQS* eqs_dof = NULL; // WW 02.2023. Pardiso
	Linear_EQS* eqsH = NULL;
#endif
	//
	Linear_EQS* eqs = NULL;
	for (i = 0; i < fem_msh_vector.size(); i++)
	{
		a_msh = fem_msh_vector[i];
#if defined(USE_MPI)
		eqs = new Linear_EQS(a_msh->GetNodesNumber(true) * dof_max);
		EQS_Vector.push_back(eqs);
		EQS_Vector.push_back(eqs);
#else
		a_msh->CreateSparseTable();
		// sparse pattern with linear elements exists
		sp = a_msh->GetSparseTable();
		// sparse pattern with quadratic elements exists
		spH = a_msh->GetSparseTable(true);
		//
		eqs = NULL;
		eqsH = NULL;
		// CB_merge_0513
		eqs_dof = NULL; // WW 02.2023. Pardiso
		if (sp) // WW 02.2023. Pardiso
		{
			if (need_eqs) // 02.2013. WW
			{
				// eqs = new Linear_EQS(*sp, dof_nonDM);//WW 02.2023. Pardiso
				eqs = new Linear_EQS(*sp, 1);
			}
			if (need_eqs_dof)
			{
				eqs_dof = new Linear_EQS(*sp, dof_nonDM_max);
			}
		} // WW 02.2023. Pardiso
		if (spH)
			eqsH = new Linear_EQS(*spH, dof_DM);
		EQS_Vector.push_back(eqs);
		EQS_Vector.push_back(eqsH);
		EQS_Vector.push_back(eqs_dof); // WW 02.2023. Pardiso
#endif
	}
}

#else // NEW_EQS
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   09/2007 WW Implementation
 **************************************************************************/
#include <iomanip>
void CRFProcess::DumpEqs(string file_name)
{
	fstream eqs_out;
	eqs_out.open(file_name.c_str(), ios::out);
	eqs_out.setf(ios::scientific, ios::floatfield);
	setw(14);
	eqs_out.precision(14);
	//
	long nnode = eqs->dim / eqs->unknown_vector_dimension;
	for (long i = 0; i < eqs->dim; i++)
	{
		CNode const* const node(m_msh->nod_vector[i % nnode]);
		std::vector<size_t> const& connected_nodes(node->getConnectedNodes());
		const size_t n_connected_nodes(connected_nodes.size());
		for (int ii = 0; ii < eqs->unknown_vector_dimension; ii++)
			for (size_t j = 0; j < n_connected_nodes; j++)
			{
				const long k = ii * nnode + connected_nodes[j];
				if (k >= eqs->dim)
					continue;
				eqs_out << i << "  " << k << "  " << MXGet(i, k) << "\n";
			}
		eqs_out << i << "  " << eqs->b[i] << "  " << eqs->x[i] << "\n";
	}
	eqs_out.close();
}
#endif // ifdef NEW_QES

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   01/2008 WW Implementation
 **************************************************************************/
void CRFProcess::WriteBC()
{
	const size_t size_bc(bc_node_value.size());
	const size_t size_st(st_node_value.size());

	if (size_bc == 0 && size_st == 0)
		return;

#ifdef USE_PETSC
	std::string m_file_name = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_BC_ST_"
	                          + number2str(myrank) + ".asc";
#else
	std::string m_file_name = FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_BC_ST.asc";
#endif
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << "\n";
		abort();
	}
	os.setf(ios::scientific, ios::floatfield);
	os.precision(12);
	long nindex = 0;
	if (size_bc > 0)
	{
		os << "#Dirchilet BC  (from " << m_file_name << ".bc file) "
		   << "\n";
		os << "#Total BC nodes  " << size_bc << "\n";
		os << "#Node index, name, x, y, z,   value: "
		   << "\n";
		for (size_t i = 0; i < size_bc; i++)
		{
			nindex = bc_node_value[i]->geo_node_number;
#ifdef USE_PETSC
			std::size_t id_in_eqs = m_msh->nod_vector[nindex]->GetEquationIndex();
#else
			std::size_t id_in_eqs = nindex;
#endif
			//         anode = m_msh->nod_vector[nindex];
			//         os << nindex << "  " << bc_node_value[i]->pcs_pv_name << " "
			//            << std::setw(14) << anode->X() << " " << std::setw(14) << anode->Y()
			//            << " " << std::setw(14) << anode->Z() << " " << std::setw(14)
			//            << bc_node_value[i]->node_value << "\n";
			double const* const pnt(m_msh->nod_vector[nindex]->getData());
			os << id_in_eqs << "  " << bc_node_value[i]->pcs_pv_name << " " << std::setw(14) << pnt[0] << " "
			   << std::setw(14) << pnt[1] << " " << std::setw(14) << pnt[2] << " " << std::setw(14)
			   << bc_node_value[i]->node_value << "\n";
		}
	}
	if (size_st > 0)
	{
		os << "#Source term or Neumann BC  (from " << m_file_name << ".st file) "
		   << "\n";
		os << "#Total ST nodes  " << size_st << "\n";
		os << "#Node index, x, y, z, name    value: "
		   << "\n";
		for (size_t i = 0; i < size_st; i++)
		{
			nindex = st_node_value[i]->geo_node_number;
#ifdef USE_PETSC
			std::size_t id_in_eqs = m_msh->nod_vector[nindex]->GetEquationIndex();
#else
			std::size_t id_in_eqs = nindex;
#endif
			//         anode = m_msh->nod_vector[nindex];
			//         os << nindex << "  " << convertPrimaryVariableToString(
			//            st_node[i]->getProcessPrimaryVariable()) << " " << std::setw(14)
			//            << anode->X() << " " << std::setw(14) << anode->Y() << " "
			//            << std::setw(14) << anode->Z() << " " << std::setw(14)
			//            << st_node_value[i]->node_value << "\n";
			double const* const pnt(m_msh->nod_vector[nindex]->getData());
			os << id_in_eqs << "  " << convertPrimaryVariableToString(st_node[i]->getProcessPrimaryVariable()) << " "
			   << std::setw(14) << pnt[0] << " " << std::setw(14) << pnt[1] << " " << std::setw(14) << pnt[2] << " "
			   << std::setw(14) << st_node_value[i]->node_value << "\n";
		}
	}
	os << "#STOP"
	   << "\n";
	os.close();
}

/*************************************************************************
   GeoSys-Function:
   Task: PI time step contorl
   Programming:
   08/2008 WW Implementation
   10/2008 WW Node value criteria (test)
   03/2009 WW Euclidean norm
 **************************************************************************/
void CRFProcess::PI_TimeStepSize()
{
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	// Time step control
	double hmin;
	double hmax;
	double factor1; // 1/hmin
	double factor2; // 1/hmax
	double sfactor = 0.9;
	// WW double reject_factor;                          // BG

	double* u_n = _problem->GetBufferArray();

	double* eqs_x = NULL;
	if (m_num->nls_method == 1) // Newton-Raphson
	{
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
		eqs_x = eqs_new->GetGlobalSolution();
#else
#ifdef NEW_EQS
		eqs_x = eqs_new->x;
#else
		eqs_x = eqs->x;
#endif
#endif
	}

	//
	//
	hmax = Tim->GetMaximumTSizeRestric();
	hmin = Tim->GetMinimumTSizeRestric();
	//
	if (hmin < DBL_MIN)
		factor1 = 5.0;
	else
		factor1 = 1.0 / hmin;
	if (hmax < DBL_MIN)
		factor2 = 0.166666666666666666667e+00;
	else
		factor2 = 1.0 / hmax;
	if (factor1 < 1.0e0)
		factor1 = 5.0;
	if (factor2 > 1.0e0)
		factor2 = 0.166666666666666666667e+00;
	//
	hmax = Tim->max_time_step;
	if (hmax < DBL_MIN)
		hmax = fabs(Tim->time_end - aktuelle_zeit);
	//
	// esitmate the error
	double hnew;
	double err, fac;
	double factorGus;
	double hacc = Tim->GetHacc();
	double erracc = Tim->GetErracc();
//
#define aE_NORM
#ifdef E_NORM
	//
	long i;
	CElem* elem = NULL;
	bool Check2D3D;
	double norm_e, norm_en;
	double norm_e_rank, norm_en_rank;
	norm_e = norm_en = norm_e_rank = norm_en_rank = 0.;

	Check2D3D = false;
	if (type == 66) // Overland flow
		Check2D3D = true;
	//----------------------------------------------------------------------
	// DDC
	if (dom_vector.size() > 0)
	{
		cout << "      Domain Decomposition" << '\n';
		CPARDomain* m_dom = NULL;
		int j = 0;
//
#if defined(USE_MPI)
		j = myrank;
#else
			for (j = 0; j < (int)dom_vector.size(); j++)
			{
#endif
		m_dom = dom_vector[j];
		for (int ii = 0; ii < (int)continuum_vector.size(); ii++)
		{
			continuum = ii;
			//
			for (i = 0; i < (long)m_dom->elements.size(); i++)
			{
				elem = m_msh->ele_vector[m_dom->elements[i]];
				if (elem->GetMark())
				{
					elem->SetOrder(false);
					fem->SetElementNodesDomain(m_dom->element_nodes_dom[i]);
					fem->ConfigElement(elem, Check2D3D);
					fem->m_dom = m_dom;
					fem->CalcEnergyNorm(norm_e_rank, norm_en_rank);
					// _new
					if (ii == 1)
						fem->CalcEnergyNorm_Dual(norm_e_rank, norm_en_rank);
				}
			}
		}
#if defined(USE_MPI)
		MPI_Allreduce(&norm_e_rank, &norm_e, 1, MPI_DOUBLE, MPI_SUM, comm_DDC);
		MPI_Allreduce(&norm_en_rank, &norm_en, 1, MPI_DOUBLE, MPI_SUM, comm_DDC);
#else // USE_MPI
				norm_e += norm_e_rank;
				norm_en += norm_en_rank;
			}
//....................................................................
#endif
	}
	//----------------------------------------------------------------------
	// STD
	else
		for (int ii = 0; ii < (int)continuum_vector.size(); ii++)
		{
			continuum = ii;
			for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
			{
				elem = m_msh->ele_vector[i];
				if (elem->GetMark()) // Marked for use
				{
					elem->SetOrder(false);
					fem->ConfigElement(elem, Check2D3D);
					fem->CalcEnergyNorm(u_n, norm_e, norm_en);
					// _new
					if (ii == 1)
						fem->CalcEnergyNorm_Dual(u_n, norm_e, norm_en);
				}
			}
		}
	// compute energy norm as the error
	err = sqrt(fabs(norm_e / norm_en));
#else // ifdef E_NORM
	err = 0.0;
	//
	int ii, nidx1;
	long g_nnodes, j, k, l, size_x;
	double x0, x1;
	double Rtol = Tim->GetRTol();
	double Atol = Tim->GetATol();
	double* u_k = _problem->GetBufferArray(true);

	size_x = 0;
	for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
		nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
		g_nnodes = m_msh->getNumNodesLocal();
#else
		g_nnodes = m_msh->GetNodesNumber(false);
#endif
		size_x += g_nnodes;

		if (m_num->nls_method == 1) // Newton-Raphson
		{
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
				k = j;
				l = pcs_number_of_primary_nvals * j + ii;
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				l = j + ii * g_nnodes;
#endif
				x0 = u_n[l];
				x1 = GetNodeValue(k, nidx1);
				err += pow((eqs_x[l]) / (Atol + Rtol * max(fabs(x0), fabs(x1))), 2);
			}
		}
		else
		{
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
				k = j;
				l = pcs_number_of_primary_nvals * j + ii;
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				l = j + ii * g_nnodes;
#endif
				x0 = u_n[l];
				x1 = GetNodeValue(k, nidx1);
				err += pow((x1 - u_k[l]) / (Atol + Rtol * max(fabs(x0), fabs(x1))), 2);
			}
		}
	}

#if defined(USE_PETSC) // || defined(other parallel libs)//04.3012. WW
	double err_l = err;
	MPI_Allreduce(&err_l, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	long size_xloc = size_x;
	MPI_Allreduce(&size_xloc, &size_x, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
	err = sqrt(err / (double)size_x);
#endif

	//----------------------------------------------------------------------
	//
	// Set the reject factor for the first timestep BG
	// if (Tim->step_current == 1)
	//   Tim->reject_factor = 1;

	// compute hnew with the resriction: 0.2<=hnew/h<=6.0;
	fac = max(factor2, min(factor1, pow(err, 0.25) / sfactor));
	hnew = Tim->time_step_length
	       / fac; //*Tim->reject_factor;              // BG, use the reject factor to lower timestep after rejection BG

	// determine if the error is small enough
	if (err <= 1.0e0 && accepted) // step is accept (unless Newton diverged!)
	{
		accept_steps++;
		if (Tim->GetPITimeStepCrtlType() == 2) // Mod. predictive step size controller (Gustafsson)
		{
			if (accept_steps > 1)
			{
				factorGus = (hacc / Tim->time_step_length) * pow(err * err / erracc, 0.25) / sfactor;
				factorGus = max(factor2, min(factor1, factorGus));
				fac = max(fac, factorGus);
				hnew = Tim->time_step_length / fac;
			}
			hacc = Tim->time_step_length;
			erracc = max(1.0e-2, err);
			Tim->SetHacc(hacc);
			Tim->setErracc(erracc);
		}
		if (fabs(hnew) > hmax)
			hnew = hmax;
		if (!accepted)
			hnew = min(fabs(hnew), Tim->time_step_length);
		Tim->SetTimeStep(fabs(hnew));
		// store the used time steps for post-processing BG
		if (Tim->step_current == 1) // BG
			Tim->time_step_vector.push_back(Tim->time_step_length);
		Tim->time_step_vector.push_back(hnew);
		// WW reject_factor = 1;
		// end of time storage BG
	} // end if(err<=1.0e0)
	else
	{
		if (!accepted && err <= 1.0e0)
		{ // JT: Then error suggests success, but the iteration diverged.
			if (hnew / Tim->time_step_length > 0.99) // Shock the system to escape the stagnation.
				hnew = Tim->time_step_length * 0.8;
		}
		//
		// WW Tim->reject_factor = 1;                     //BG; if the time step is rejected the next timestep increase
		// is reduced by the reject factor (choose reject factor between 0.1 and 0.9); 1.0 means no change
		reject_steps++;
		accepted = false;
		// WW hnew = hnew / Tim->reject_factor;           //BG
		Tim->SetTimeStep(hnew);
		Tim->time_step_vector.push_back(hnew); // BG
		if (reject_steps > 100 && accept_steps == 0)
		{
			cout << "!!! More than 100 steps rejected and none of steps accepted. Quit the simulation now"
			     << "\n";
			exit(1);
		}
		// Recover solutions
	}
}

#ifdef NEW_EQS // WW
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  //For fluid momentum, WW
   Programming:
   01/2008 WW Implementation
 **************************************************************************/
void CRFProcess::EQSInitialize()
{
	eqs_new->Initialize();
}

#if defined(LIS) // 11.03.2008 WW
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  //For fluid momentum,
   Programming:
   02/2008 PCH Implementation
   03/2009 PCH option to tell if this is FLUID_MOMENTUM
 **************************************************************************/
void CRFProcess::EQSSolver(double* x)
{
	eqs_new->Solver(this->m_num); // NW

	// OK411
	for (int i = 0; i < (int)m_msh->nod_vector.size(); ++i)
		x[i] = eqs_new->X(i);
}
#endif
#endif

#ifdef GEM_REACT
void CRFProcess::IncorporateSourceTerms_GEMS(void)
{
	// Initialization
	long it; // iterator
	long N_Nodes; // Number of Nodes
	int nDC; // Number of mass transport components.
	int i = 0; // index of the component

	// Get a vector pointing to the REACT_GEM class
	if (m_vec_GEM)
	{
		// Get needed informations.----------------------------
		// Number of Nodes
		N_Nodes = (long)m_msh->GetNodesNumber(false);
		// Number of DC
		nDC = m_vec_GEM->nDC;
		// Identify which PCS it is and its sequence in mcp.
		i = this->pcs_component_number;
		// ----------------------------------------------------

		// Loop over all the nodes-----------------------------
		for (it = 0; it < N_Nodes /*Number of Nodes*/; it++)
		{
// Adding the rate of concentration change to the right hand side of the equation.
#ifdef NEW_EQS // 15.12.2008. WW
			eqs_new->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] / Tim->time_step_length;
#elif defined(USE_PETSC)
// eqs_new->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] / Tim->time_step_length;
#else
			eqs->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] / Tim->time_step_length;
#endif
		}
		// ----------------------------------------------------
	}
}
#endif
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   10/2008 //WW/CB Implementation
   07/2010 TF substituted GEOGetPLYByName
   03/2010  WW Read binary file of precipitation
   06/2010  WW Output top surface flux and head to a raster file
 **************************************************************************/
void CRFProcess::UpdateTransientBC()
{
	//--------------- 24.03.2010. WW
	long i;
	CSourceTerm* precip;
	CSourceTerm* a_st;
	precip = NULL;

	for (i = 0; i < (long)st_vector.size(); i++)
	{
		a_st = st_vector[i];
		if (a_st->getProcessDistributionType() == FiniteElement::PRECIPITATION
		    && a_st->getProcessType() == getProcessType())
		{
			precip = a_st;
			if (!m_msh->top_surface_checked) // 07.06--19.08.2010. WW
			{
				if (m_msh->GetCoordinateFlag() / 10 == 3)
					m_msh->MarkInterface_mHM_Hydro_3D();
				m_msh->top_surface_checked = true;
			}
			break;
		}
	}
	if (precip) // 08-07.06.2010.  WW
	{
		string ofile_name;
		ofile_name = precip->DirectAssign_Precipitation(aktuelle_zeit);

		if (m_msh->GetCoordinateFlag() / 10 == 3) // 19.08.2010. WW
		{
			/// Remove .bin from the file name
			i = (int)ofile_name.find_last_of(".");
			if (i > 0)
				ofile_name.erase(ofile_name.begin() + i, ofile_name.end());
			i = (int)ofile_name.find_last_of(".");
			if (i > 0)
				ofile_name.erase(ofile_name.begin() + i, ofile_name.end());
			string of_name = ofile_name + ".flx.asc";
			ofstream of_flux(of_name.c_str(), ios::trunc | ios::out);

			of_name = ofile_name + ".pri.asc";
			ofstream of_primary(of_name.c_str(), ios::trunc | ios::out);

			/// GIS_shape_head[0]:  ncols
			/// GIS_shape_head[1]:  nrows
			/// GIS_shape_head[2]:  xllcorner
			/// GIS_shape_head[3]:  yllcorner
			/// GIS_shape_head[4]:  cellsize
			/// GIS_shape_head[5]:  NONDATA_value
			double* g_para = precip->GIS_shape_head;
			long size = (long)(g_para[0] * g_para[1]);
			vector<double> cell_data_p(size);
			vector<double> cell_data_v(size);
			for (i = 0; i < size; i++)
			{
				cell_data_p[i] = g_para[5];
				cell_data_v[i] = g_para[5];
			}

			int j, k, nnodes;
			int node_xmin, node_xmax, node_ymin, node_ymax;
			long m, n, mm, nn, l;
			CElem* elem;
			double* cent;
			double vel_av[3], x1[3], x2[3], x3[3], sub_area[3], area, tol_a;

			double x_min, y_min, x_max, y_max;
			long row_min, col_min, row_max, col_max;

			int* vel_idx;
			int idx = fem->idx0 + 1;
			vel_idx = fem->idx_vel;

			long n_idx, irow, icol, nrow, ncol;
			nrow = (long)g_para[1];
			ncol = (long)g_para[0];

			node_xmin = node_xmax = node_ymin = node_ymax = 0;

			tol_a = 1.e-8;

			of_flux.setf(ios::fixed, ios::floatfield);
			of_primary.setf(ios::fixed, ios::floatfield);
			of_flux.precision(1);
			of_primary.precision(1);

			of_flux << "ncols" << setw(19) << ncol << "\n";
			of_flux << "nrows" << setw(19) << nrow << "\n";
			of_flux << "xllcorner" << setw(15) << g_para[2] << "\n";
			of_flux << "yllcorner" << setw(15) << g_para[3] << "\n";
			of_flux << "cellsize" << setw(16) << (long)g_para[4] << "\n";
			of_flux << "NODATA_value" << setw(11) << (long)g_para[5] << "\n";

			of_primary << "ncols" << setw(19) << ncol << "\n";
			of_primary << "nrows" << setw(19) << nrow << "\n";
			of_primary << "xllcorner" << setw(15) << g_para[2] << "\n";
			of_primary << "yllcorner" << setw(15) << g_para[3] << "\n";
			of_primary << "cellsize" << setw(16) << (long)g_para[4] << "\n";
			of_primary << "NODATA_value" << setw(11) << (long)g_para[5] << "\n";

			for (i = 0; i < (long)m_msh->face_vector.size(); i++)
			{
				elem = m_msh->face_vector[i];
				if (!elem->GetMark())
					continue;

				if (elem->GetElementType() != 4) /// If not triangle
					continue;

				//// In element
				nnodes = elem->GetNodesNumber(false);
				cent = elem->gravity_center;

				/// Find the range of this element
				x_min = y_min = 1.e+20;
				x_max = y_max = -1.e+20;
				for (k = 0; k < nnodes; k++)
				{
					double const* const pnt(elem->nodes[k]->getData());
					if (pnt[0] < x_min)
					{
						x_min = pnt[0];
						node_xmin = k;
					}
					if (pnt[0] > x_max)
					{
						x_max = pnt[0];
						node_xmax = k;
					}
					if (pnt[1] < y_min)
					{
						y_min = pnt[1];
						node_ymin = k;
					}
					if (pnt[1] > y_max)
					{
						y_max = pnt[1];
						node_ymax = k;
					}
				}

				/// Determine the cells that this element covers. 05.10. 2010
				col_min = (long)((x_min - g_para[2]) / g_para[4]);
				row_min = nrow - (long)((y_max - g_para[3]) / g_para[4]);
				col_max = (long)((x_max - g_para[2]) / g_para[4]);
				row_max = nrow - (long)((y_min - g_para[3]) / g_para[4]);

				double const* const pnt1(elem->nodes[0]->getData());
				x1[0] = pnt1[0];
				x1[1] = pnt1[1];
				double const* const pnt2(elem->nodes[1]->getData());
				x2[0] = pnt2[0];
				x2[1] = pnt2[1];
				double const* const pnt3(elem->nodes[2]->getData());
				x3[0] = pnt3[0];
				x3[1] = pnt3[1];

				x3[2] = x2[2] = x1[2] = 0.;
				cent[2] = 0.;

				for (m = col_min; m <= col_max; m++)
				{
					mm = m;
					if (m > ncol - 1)
						mm = ncol - 1;
					if (m < 0)
						mm = 0;
					cent[0] = g_para[2] + g_para[4] * (mm + 0.5);

					for (n = row_min; n <= row_max; n++)
					{
						nn = n;
						if (n > nrow - 1)
							nn = nrow - 1;
						if (nn < 0)
							nn = 0;
						cent[1] = g_para[3] + g_para[4] * (nrow - nn + 0.5);

						if (cent[0] < x_min)
						{
							double const* const pnt(elem->nodes[node_xmin]->getData());
							cent[0] = pnt[0];
							cent[1] = pnt[1];
						}
						if (cent[0] > x_max)
						{
							double const* const pnt(elem->nodes[node_xmax]->getData());
							cent[0] = pnt[0];
							cent[1] = pnt[1];
						}
						if (cent[1] < y_min)
						{
							double const* const pnt(elem->nodes[node_ymin]->getData());
							cent[0] = pnt[0];
							cent[1] = pnt[1];
						}
						if (cent[1] < y_max)
						{
							double const* const pnt(elem->nodes[node_ymax]->getData());
							cent[0] = pnt[0];
							cent[1] = pnt[1];
						}

						/// Check whether this point is in this element.
						sub_area[0] = ComputeDetTri(cent, x2, x3);
						sub_area[1] = ComputeDetTri(cent, x3, x1);
						sub_area[2] = ComputeDetTri(cent, x1, x2);
						area = ComputeDetTri(x1, x2, x3);

						/// This point locates within the element
						if (fabs(area - sub_area[0] - sub_area[1] - sub_area[2]) < tol_a)
						{
							/// Use sub_area[k] as shape function
							for (k = 0; k < 3; k++)
								sub_area[k] /= area;

							l = nn * ncol + mm;
							cell_data_p[l] = 0.0;

							for (k = 0; k < 3; k++)
								vel_av[k] = 0.;

							for (j = 0; j < nnodes; j++)
							{
								n_idx = elem->GetNodeIndex(j);

								cell_data_p[l] += sub_area[j] * GetNodeValue(n_idx, idx);

								for (k = 0; k < 3; k++)
									vel_av[k] += sub_area[j] * GetNodeValue(n_idx, vel_idx[k]);
							}
							cell_data_v[l] = 1000.0 * (vel_av[0] * (*elem->transform_tensor)(0, 2)
							                           + vel_av[1] * (*elem->transform_tensor)(1, 2)
							                           // 1000*:  m-->mm
							                           + vel_av[2] * (*elem->transform_tensor)(2, 2));
						}
					}
				}
			}

			of_flux.precision(4);
			// of_flux.setf(ios::scientific, ios::floatfield);
			of_primary.precision(2);
			for (irow = 0; irow < nrow; irow++)
			{
				for (icol = 0; icol < ncol; icol++)
				{
					m = irow * ncol + icol;
					of_flux << " " << setw(9) << cell_data_v[m];
					of_primary << setw(11) << cell_data_p[m];
				}
				of_flux << "\n";
				of_primary << "\n";
			}

			cell_data_p.clear();
			cell_data_v.clear();
			of_flux.close();
			of_primary.close();
		}
	}

	// transient boundary condition
	if (bc_transient_index.size() == 0)
		return;

	bool valid = false;
	long end_i = 0;
	double t_fac = 0.;
	std::vector<double> node_value;

	for (size_t i = 0; i < bc_transient_index.size(); i++)
	{
		std::vector<double> interpolation_points;
		std::vector<double> interpolation_values;

		CBoundaryCondition* bc = bc_node[bc_transient_index[i]];
		long start_i = bc_transient_index[i];
		if (i == bc_transient_index.size() - 1)
			end_i = (long)bc_node.size();
		else
			end_i = bc_transient_index[i + 1];

		// fetch points (representing mesh nodes) along polyline for interpolation
		std::vector<double> nodes_as_interpol_points;
		GEOLIB::Polyline const* ply(static_cast<GEOLIB::Polyline const*>(bc->getGeoObj()));
		m_msh->getPointsForInterpolationAlongPolyline(ply, nodes_as_interpol_points);

		valid = false;
		t_fac = 0.0;
		// Piecewise linear distributed.
		for (size_t i(0); i < bc->getDistribedBC().size(); i++)
		{
			for (size_t j = 0; j < ply->getNumberOfPoints(); j++)
				if (bc->getPointsWithDistribedBC()[i] == (int)ply->getPointID(j))
				{
					if (fabs(bc->getDistribedBC()[i]) < MKleinsteZahl)
						bc->getDistribedBC()[i] = 1.0e-20;
					interpolation_points.push_back(ply->getLength(j));
					interpolation_values.push_back(bc->getDistribedBC()[i]);

					CFunction* fct(FCTGet(bc->getPointsFCTNames()[i]));
					if (fct)
						t_fac = fct->GetValue(aktuelle_zeit, &valid);
					else
						std::cout << "Warning in CBoundaryConditionsGroup - no FCT data"
						          << "\n";

					if (valid)
						interpolation_values[interpolation_values.size() - 1] *= t_fac;

					break;
				}
		}
		std::vector<double> interpol_res;
		MathLib::PiecewiseLinearInterpolation(
		    interpolation_points, interpolation_values, nodes_as_interpol_points, interpol_res);

		for (long k = start_i; k < end_i; k++)
			bc_node_value[k]->node_value = interpol_res[k - start_i];
	}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::Add_GEMS_Water_ST
   Task: add one more record into the vector Water_ST_vec
   Programming:
   11/2008 //HS Implementation
 **************************************************************************/

void CRFProcess::Add_GEMS_Water_ST(long idx, double val)
{
	Water_ST_GEMS tmp_st;

	tmp_st.index_node = idx;
	tmp_st.water_st_value = val;

	Water_ST_vec.push_back(tmp_st);
}

void CRFProcess::Clean_Water_ST_vec()
{
	Water_ST_vec.clear();
}

void CRFProcess::configMaterialParameters()
{
	// Output material parameters
	// WW
	const size_t out_vector_size(out_vector.size());

	for (size_t i = 0; i < out_vector_size; i++)
	{
		COutput* out = out_vector[i];
		const size_t size(out->_nod_value_vector.size());
		for (size_t k = 0; k < size; k++)
			if (out->_nod_value_vector[k].find("PERMEABILITY_X1") != string::npos)
			{
				additioanl2ndvar_print = 1;
				break;
			}
		if (additioanl2ndvar_print == 1)
			break;
	}

	for (size_t i = 0; i < out_vector_size; i++)
	{
		COutput* out = out_vector[i];
		const size_t size(out->_nod_value_vector.size());
		for (size_t k = 0; k < size; k++)
		{
			if (out->_nod_value_vector[k].find("POROSITY") != string::npos)
			{
				if (additioanl2ndvar_print > 0)
					additioanl2ndvar_print = 2;
				else
					additioanl2ndvar_print = 3;
			}
			if (additioanl2ndvar_print > 1)
				break;
		}
		if (additioanl2ndvar_print > 1)
			break;
	}

	if (additioanl2ndvar_print > 0) // WW
	{
		if (additioanl2ndvar_print < 3)
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PERMEABILITY_X1";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "1/m^2";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PERMEABILITY_Y1";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "1/m^2";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
			if (max_dim == 2) // 3D
			{
				pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PERMEABILITY_Z1";
				pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "1/m^2";
				pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
				pcs_number_of_secondary_nvals++;
			}
		}
		if (additioanl2ndvar_print > 1) // WW
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "POROSITY";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
		}
	}
}

/*************************************************************************
   GeoSys - Function:
   06/2009 OK Implementation
 **************************************************************************/
bool PCSConfig()
{
	bool some_thing_done = false;
	for (int i = 0; i < (int)pcs_vector.size(); i++) // OK
	{
		pcs_vector[i]->Config();
		some_thing_done = true;
	}
	return some_thing_done;
}

/*************************************************************************
   GeoSys-Function: CalGPVelocitiesfromEclipse
   Task: Calculate gauss point velocities from Eclipse solution
      extrapolate velocities from nodes to gauss points
   Programming:
   09/2009 SB BG Implementation

 **************************************************************************/
void CRFProcess::CalGPVelocitiesfromECLIPSE(string path, int timestep, int phase_index, string phase)
{
	(void)path; // unused
	(void)timestep; // unused
	long i;
	MeshLib::CElem* elem = NULL;
	clock_t start, finish;
	double time;

	start = clock();

	cout << "        CalGPVelocitiesfromECLIPSE() ";

	// Test Output
	vector<string> vec_string;
	string tempstring;
	ostringstream temp;
	vec_string.push_back("Element; X; Y; Z; v_Geosys_x; v_Geosys_y; v_Geosys_z; v_Eclipse_x; v_Eclipse_y; v_Eclipse_z");

	// Loop over all elements
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i]; // get element
		if (elem->GetMark()) // Marked for use
		{
			// Test Output
			tempstring = "";
			temp.str("");
			temp.clear();
			temp << i;
			tempstring = temp.str();
			double const* gc(elem->GetGravityCenter());
			temp.str("");
			temp.clear();
			temp << gc[0];
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp << gc[1];
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp << gc[2];
			tempstring += "; " + temp.str();

			// Configure Element for interpolation of node velocities to GP velocities
			fem->ConfigElement(elem);
			// Interpolate from nodes to GP of actual element
			// cout << "Element: " << i << "\n";
			tempstring = fem->Cal_GP_Velocity_ECLIPSE(tempstring, true, phase_index, phase);

			// Test Output
			vec_string.push_back(tempstring);
		}
	} // end element loop

	// // Test Output
	// if (timestep == 1 || timestep % 10 == 0 ) {
	// int position = path.find_last_of("\\");
	// path = path.substr(0,position);
	// position = path.find_last_of("\\");
	// path = path.substr(0,position);
	// temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
	// string aus_file = path + "\\CheckVelocity_" + phase + "_" + tempstring + ".csv";
	// ofstream aus;
	// aus.open(aus_file.data(),ios::out);
	// for (i = 0; i < vec_string.size(); i++) {
	//	aus << vec_string[i] << "\n";
	//}
	// aus.close();
	//}

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "         Time: " << time << " seconds."
	     << "\n";
}

//-------------------------------------------------------------------------
// GeoSys - Function: CO2_H2O_NaCl_VLE_isobaric
// Task: Calculates the isobaric vapour liquid equilibrium for the CO2-H2O-NaCl System
// Return: nothing
// Programming: 11/2010 DL, BG
// Modification:
//-------------------------------------------------------------------------
void CRFProcess::CO2_H2O_NaCl_VLE_isobaric(
    double T, double P, Phase_Properties& vapor, Phase_Properties& liquid, Phase_Properties& solid, int f)
{
	// f ---  1 V-L-S, 2 V-L, 3 CO2-L-S, 4 CO2-L

	double tCO2, tH2O, tNaCl;
	double wH2O, mNaCl = 0, soluCO2 = 0, xCO2, yH2O, AW = 1.0; // ,mCO2; unused
	double a, b, tvCO2 = 0, tvH2O = 0, tH2Or = 0, er = 1.0e-6, err = 1.0e-6;
	int i, iter_max = 100;
	double Molweight_CO2, Molweight_H2O, Molweight_NaCl;

	Molweight_CO2 = 44.009;
	Molweight_H2O = 18.0148;
	Molweight_NaCl = 58.443; // ToDo: provide constants once in the whole project

	tCO2 = vapor.CO2 + liquid.CO2 + solid.CO2;
	tH2O = vapor.H2O + liquid.H2O + solid.H2O;
	tNaCl = vapor.NaCl + liquid.NaCl + solid.NaCl;
	yH2O = VLE::fraction_H2O(T, P, AW); // water fraction in vapor
	//	//cout << " yH2O " << yH2O << " xCO2 " << xCO2 << "\n";

	// cout << " flag " << f << "\n";
	// V-L-S
	if (f == 1 || f == 2)
	{
		// cout << " >>>>> phase equilibrium <<<<< " << "\n";
		a = 0.0;
		b = tCO2;
		for (i = 0; i < iter_max; i++)
		{
			tvCO2 = (a + b) / 2.0;
			tvH2O = tvCO2 * yH2O / (1 - yH2O);
			wH2O = (tH2O - tvH2O) * Molweight_H2O / 1000;

			mNaCl = tNaCl / wH2O;
			if (f == 1)
				if (mNaCl > 4.5)
					mNaCl = 4.5;
			// NaCl solubility calc

			soluCO2 = VLE::solubility_CO2(T, P, mNaCl);
			//???????????????????????????? (soluCO2 <= tCO2?)
			xCO2 = soluCO2 / (soluCO2 + 1000 / Molweight_H2O);
			tH2Or = tvH2O + (tCO2 - tvCO2) * (1 - xCO2) / xCO2;
			// cout << " tH2Or " << tH2Or << " tH2O " << tH2O << "\n";

			if (abs(tH2Or - tH2O) < er)
				break; // two phases coexisting
			if (tvCO2 < err) // single liquid phase
			{
				tvCO2 = 0;
				break;
			}
			if (tCO2 - tvCO2 < err) // single vapor phase
			{
				tvCO2 = tCO2;
				break;
			}
			if (tH2Or > tH2O)
				a = tvCO2;
			else
				b = tvCO2;
		}

		// single liquid phase
		if (tvCO2 == 0)
		{
			// cout << " >>>>> single liquid phase <<<<< " << "\n";
			vapor.CO2 = 0.0;
			vapor.H2O = 0.0;
			vapor.NaCl = 0.0;

			wH2O = tH2O * Molweight_H2O / 1000;
			mNaCl = tNaCl / wH2O;
			if (f == 1)
				if (mNaCl > 4.5)
					mNaCl = 4.5;
			// NaCl solubility calc

			liquid.CO2 = tCO2;
			liquid.H2O = tH2O;
			liquid.NaCl = mNaCl * wH2O;

			vapor.temperature = T;
			vapor.pressure = P;
			vapor.density = -1.0;
			vapor.viscosity = -1.0;
			vapor.volume = 0.0;
			vapor.mass = 0.0;

			liquid.temperature = T;
			liquid.pressure = P;
			liquid.density = density::CO2brine(T, P, mNaCl, liquid.CO2 / wH2O);
			liquid.viscosity = -1.0;
			liquid.mass = liquid.H2O * Molweight_H2O + liquid.CO2 * Molweight_CO2 + liquid.NaCl * Molweight_NaCl;
			liquid.volume = liquid.mass / liquid.density;
		}

		// single vapor phase
		else if (tvCO2 == tCO2)
		{
			// cout << " >>>>> single vapor phase <<<<< " << "\n";
			vapor.CO2 = tCO2;
			vapor.H2O = tH2O;
			vapor.NaCl = 0.0;

			liquid.CO2 = 0.0;
			liquid.H2O = 0.0;
			liquid.NaCl = 0.0;

			liquid.mass = 0.0;
			liquid.volume = 0.0;

			solid.NaCl = tNaCl;

			vapor.temperature = T;
			vapor.pressure = P;

			vapor.density = VLE::density_CO2(T, P); // to use mixture fluid EoS
			vapor.viscosity = -1.0;
			vapor.mass = vapor.CO2 * Molweight_CO2 + vapor.H2O * Molweight_H2O;
			vapor.volume = vapor.mass / vapor.density;
		}

		// two phases coexisting
		else
		{
			// cout << " >>>>> two phases coexisting <<<<< " << "\n";
			vapor.H2O = tvH2O;
			vapor.CO2 = tvCO2;
			vapor.NaCl = 0.0;

			liquid.H2O = tH2O - tvH2O;
			liquid.CO2 = tCO2 - tvCO2;
			liquid.NaCl = mNaCl * wH2O;

			vapor.temperature = T;
			vapor.density = VLE::density_CO2(T, P); // to use mixture fluid EoS
			vapor.viscosity = -1;
			vapor.mass = vapor.CO2 * Molweight_CO2 + vapor.H2O * Molweight_H2O;
			vapor.pressure = VLE::pressure_CO2(T, vapor.density);
			vapor.volume = vapor.mass / vapor.density;

			liquid.temperature = T;
			liquid.pressure = P;
			liquid.density = density::CO2brine(T, P, mNaCl, soluCO2);
			liquid.viscosity = -1.0;
			liquid.mass = liquid.CO2 * Molweight_CO2 + liquid.H2O * Molweight_H2O + liquid.NaCl * Molweight_NaCl;
			liquid.volume = liquid.mass / liquid.density;
		}
	}

	else if (f == 3 || f == 4)
	{
		wH2O = tH2O * Molweight_H2O / 1000;
		mNaCl = tNaCl / wH2O;
		if (f == 3)
			if (mNaCl > 4.5)
				mNaCl = 4.5;
		// NaCl solubility calc
		soluCO2 = VLE::solubility_CO2(T, P, mNaCl);

		// single liquid phase
		if (tCO2 < soluCO2 * wH2O)
		{
			// cout << " >>>>> single liquid phase <<<<< " << "\n";
			vapor.CO2 = 0.0;
			vapor.H2O = 0.0;
			vapor.NaCl = 0.0;

			liquid.CO2 = tCO2;
			liquid.H2O = tH2O;
			liquid.NaCl = mNaCl * wH2O;

			vapor.temperature = T;
			vapor.pressure = P;
			vapor.density = -1.0;
			vapor.viscosity = -1.0;
			vapor.volume = 0.0;
			vapor.mass = 0.0;

			liquid.temperature = T;
			liquid.pressure = P;
			liquid.density = density::CO2brine(T, P, mNaCl, liquid.CO2 / wH2O);
			liquid.viscosity = -1.0;
			liquid.mass = liquid.H2O * Molweight_H2O + liquid.CO2 * Molweight_CO2 + liquid.NaCl * Molweight_NaCl;
			liquid.volume = liquid.mass / liquid.density;
		}

		// two phases coexisting
		else
		{
			// cout << " >>>>> two phases coexisting <<<<< " << "\n";
			vapor.H2O = 0.0;
			vapor.CO2 = tCO2 - soluCO2 * wH2O;
			vapor.NaCl = 0.0;

			liquid.H2O = tH2O;
			liquid.CO2 = soluCO2 * wH2O;
			liquid.NaCl = mNaCl * wH2O;

			vapor.temperature = T;
			vapor.density = VLE::density_CO2(T, P); // to use mixture fluid EoS
			vapor.viscosity = -1;
			vapor.mass = vapor.CO2 * Molweight_CO2 + vapor.H2O * Molweight_H2O;
			vapor.pressure = VLE::pressure_CO2(T, vapor.density);
			vapor.volume = vapor.mass / vapor.density;

			liquid.temperature = T;
			liquid.pressure = P;
			liquid.density = density::CO2brine(T, P, mNaCl, soluCO2);
			liquid.viscosity = -1.0;
			liquid.mass = liquid.CO2 * Molweight_CO2 + liquid.H2O * Molweight_H2O + liquid.NaCl * Molweight_NaCl;
			liquid.volume = liquid.mass / liquid.density;
		}
	}
}

//-------------------------------------------------------------------------
// GeoSys - Function: CO2_H2O_NaCl_VLE_isobaric
// Task: Calculates the isobaric vapour liquid equilibrium for the CO2-H2O-NaCl System
// Return: nothing
// Programming: 11/2010 DL, BG
// Modification:
//-------------------------------------------------------------------------
void CRFProcess::CO2_H2O_NaCl_VLE_isochoric(Phase_Properties& vapor,
                                            Phase_Properties& liquid,
                                            Phase_Properties& solid,
                                            int f)
{
	// f ---  1 V-L-S, 2 V-L, 3 CO2-L-S, 4 CO2-L

	int i, i_max = 100;
	double V, V0, P1 = 0.0, P2 = 0.0;
	double err = 1.0e-4;
	double T, P;

	T = vapor.temperature;
	P = vapor.pressure;

	V0 = vapor.volume + liquid.volume;
	CO2_H2O_NaCl_VLE_isobaric(T, P, vapor, liquid, solid, f);
	V = vapor.volume + liquid.volume;

	// cout << " V input " << V << "\n";

	if (abs(V - V0) < err)
		return;
	else if (V > V0)
	{
		P1 = P;
		P2 = 5.0 * P; // P2=5.0*P;
	}
	else if (V < V0)
	{
		P1 = 0.2 * P; // P1=0.2*P;
		P2 = P;
	}

	for (i = 0; i < i_max; i++)
	{
		P = 0.5 * (P1 + P2);
		CO2_H2O_NaCl_VLE_isobaric(T, P, vapor, liquid, solid, f);
		V = vapor.volume + liquid.volume;
		// cout << P << " " << V << "\n";
		if (abs(V0 - V) < err)
			break;
		else if (V < V0)
			P2 = P;
		else if (V > V0)
			P1 = P;
	}
}

//-------------------------------------------------------------------------
// GeoSys - Function: CalculateFluidDensitiesAndViscositiesAtNodes
// Task: Calculates density and viscosity of a multiphase flow system
// Return: nothing
// Programming: 11/2010 BG, NB
// Modification:
//-------------------------------------------------------------------------
void CRFProcess::CalculateFluidDensitiesAndViscositiesAtNodes(CRFProcess* m_pcs)
{
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	double saturation_gas, saturation_liquid;
	double Molweight_CO2, Molweight_H2O, Molweight_NaCl;
	double Density_gas = 0, Density_liquid = 0, Density_pureCO2, Viscosity_liquid, Viscosity_gas;
	// ,c_H2OinLiquid; unused
	double c_CO2inLiquid, c_NaClinLiquid, c_H2OinGas;
	double b_CO2inPureWater, b_NaClinPureWater; // WW, b_H2OinPureCO2;
	double mass_gas;
	// ,Volume_eff; unused
	// WW double volume_liquid,
	double volume_gas, node_volume, porosity;
	// ,saturation_liquid_effective ,saturation_gas_effective; unused
	// WW  double saturation_gas_min, saturation_liquid_min;
	int variable_index, indexProcess;
	CMediumProperties* MediaProp;
	CFluidProperties* FluidProp;
	CRFProcess* pcs_MassTransport;
	double pressure, temperature;
	int MassTransportID[5] = {0, 0, 0, 0, 0};
	int MaterialGroup = 0;
	int TimeStepVariableIndex = 1;
	clock_t start, finish;
	double time;
	double a, b, temp_density; // unit: [kg/m³]
	double wCO2, wNaCl, wH2O; // unit: [kg/kg]
	double epsilon = 1e-2;
	int iter_max = 100;

	start = clock();

	// for initial conditions variable index is the first index, for other timesteps variable index is the second index
	if (m_pcs->Tim->step_current == 1)
		TimeStepVariableIndex = 0;
	else
		TimeStepVariableIndex = 1;

	Molweight_CO2 = 44.009; // [g/mol]
	Molweight_H2O = 18.0148; // [g/mol]
	Molweight_NaCl = 58.443; // [g/mol]				//ToDo: provide constants once in the whole project

	// WW saturation_gas_min = 0.00001;
	// WW saturation_liquid_min = 0.00001;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;

		// calculate Porevolume for the node based on connected elements
		for (size_t j = 0; j < m_node->getConnectedElementIDs().size(); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];
			MaterialGroup = m_ele->GetPatchIndex();
			MediaProp = mmp_vector[MaterialGroup];
			porosity = MediaProp->Porosity(m_ele->GetIndex(), m_pcs->m_num->ls_theta);
			// ToDo: Correct calculation of node volume
			node_volume = node_volume + m_ele->GetVolume() / m_ele->GetNodesNumber(false) * porosity;
		}

		// get pressure, temperature and saturations
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("PRESSURE2") + TimeStepVariableIndex;
		// unit: bar!!! Assuming that the non wetting phase pressure is the correct one
		pressure = m_pcs->GetNodeValue(i, variable_index) / 1e5;
		if (!T_Process)
		{
			// Get reference temperature if no heat transport is used
			FluidProp = MFPGet("LIQUID");
			temperature = FluidProp->getReferenceTemperature();
		}
		else
		{
			CRFProcess* HeatProcess = PCSGet("HEAT_TRANSPORT");
			//+1... new time level
			variable_index = HeatProcess->GetNodeValueIndex("TEMPERATURE1") + TimeStepVariableIndex;
			temperature = m_pcs->GetNodeValue(i, variable_index);
		}
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("SATURATION1") + TimeStepVariableIndex;
		saturation_liquid = m_pcs->GetNodeValue(i, variable_index);
		saturation_gas = 1 - saturation_liquid;

		// calculate liquid and gas volume from saturation
		// Volume_eff = node_volume * (1 - saturation_gas_min - saturation_liquid_min);		//[m³]
		//[m³]
		// WW volume_liquid = node_volume * saturation_liquid;
		volume_gas = node_volume * saturation_gas;

		// get ID's of mass transport processes
		for (int j = 0; j < int(pcs_vector.size()); j++)
		{
			pcs_MassTransport = pcs_vector[j];
			if (pcs_MassTransport->nod_val_name_vector[0] == "H2O")
				MassTransportID[0] = j;
			if ((pcs_MassTransport->nod_val_name_vector[0] == "C(4)")
			    || (pcs_MassTransport->nod_val_name_vector[0] == "CO2")) // "CO2_w" CB
				MassTransportID[1] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "NaCl")
				MassTransportID[2] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "CO2")
				MassTransportID[3] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "H2O_gas")
				MassTransportID[4] = j;
		}

		// Read CO2 concentration in water and calculate moles of CO2 in liquid phase
		indexProcess = MassTransportID[1];
		// +1: new timelevel
		variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
		                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
		//[mol/m³]
		c_CO2inLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);

		// Read NaCl concentration in water and calculate moles of NaCl in liquid phase
		indexProcess = MassTransportID[2];
		// +1: new timelevel
		variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
		                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
		//[mol/m³]
		c_NaClinLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);

		if (saturation_liquid > 0)
		{
			// Determine moles of water in liquid phase
			// calculate liquid density for initial conditions
			// estimate density -> calculate moles of components -> calculate real density -> compare with estimation ->
			// loop until similar values
			a = 800; // estimated range of valid density [kg/m³]
			b = 1800;
			for (int j = 0; j < iter_max; j++)
			{
				temp_density = (a + b) / 2.0;
				// mass fraction of kg CO2 per kg liquid (mol/m³ * m³/kg * g/mol / 1000) = kg/kg
				wCO2 = c_CO2inLiquid / temp_density * Molweight_CO2 / 1000;
				// mass fraction of kg NaCl per kg liquid (mol/m³ * m³/kg * g/mol / 1000) = kg/kg
				wNaCl = c_NaClinLiquid / temp_density * Molweight_NaCl / 1000;
				wH2O = 1 - wCO2 - wNaCl; // mass fraction of kg H2O per kg liquid
				// Molality = concentration in mol per kg of solution / mass fraction of solvent = mol per kg solvent
				// Molality of CO2 in mol per kg pure water (mol/m³ * m³/kg) = mol/kg
				b_CO2inPureWater = c_CO2inLiquid / temp_density * 1 / wH2O;
				// Molality of NaCl in mol per kg pure water (mol/m³ * m³/kg) = mol/kg
				b_NaClinPureWater = c_NaClinLiquid / temp_density * 1 / wH2O;
				// new estimate of density of brine with estimated concentration of CO2 and NaCl in mol per kg pure
				// water
				// function provides g/cm³ -> Density_liquid [kg/m³]
				Density_liquid = 1000 * density::CO2brine(temperature, pressure, b_NaClinPureWater, b_CO2inPureWater);
				// cout.precision(8);
				// cout << "Iteration: " << j << " estimated: " << temp_density << " calculated: " << Density_liquid <<
				// " a: " << a  << " b: " << b << "\n";
				// compare density with density estimation
				if (Density_liquid < 0)
					break;
				if (fabs(Density_liquid - temp_density) < epsilon)
					break;
				else if (Density_liquid < temp_density)
					b = temp_density;
				else if (Density_liquid > temp_density)
					a = temp_density;
			}
			// if (i == 230)
			//	cout << i << "\n";

			if (Density_liquid < 0)
			{
				cout << "Density calculation of water was not possible"
				     << "\n";
				//system("Pause");
				exit(0);
			}
			// calculate new moles of H2O
			// c_H2OinLiquid = (Density_liquid - c_CO2inLiquid * Molweight_CO2 * 1e-3 - c_NaClinLiquid * Molweight_NaCl
			// * 1e-3) / (Molweight_H2O * 1e-3);	//[mol]
			// indexProcess = MassTransportID[0];		//H2O in liquid
			// variable_index =
			// pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
			// // +1: new timelevel
			// pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_H2OinLiquid);

			// set new density to nodes
			variable_index = m_pcs->GetNodeValueIndex("DENSITY1");
			m_pcs->SetNodeValue(i, variable_index, Density_liquid);

			// set new viscosity to nodes
			variable_index = m_pcs->GetNodeValueIndex("VISCOSITY1");
			Viscosity_liquid = 5.1e-4;
			m_pcs->SetNodeValue(i, variable_index, Viscosity_liquid);
		}

		if (saturation_gas > 0)
		{
			// components in gas
			// Read H2O concentration in gas and calculate moles of H2O in gas phase
			indexProcess = MassTransportID[4];
			// +1: new timelevel
			variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
			                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
			//[mol/m³]
			c_H2OinGas = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);

			// calculate gas density for initial conditions
			a = 0.1; // estimated range of valid gas density [kg/m³]
			b = 1000;
			epsilon = 1e-3;
			iter_max = 100;
			for (int j = 0; j < iter_max; j++)
			{
				temp_density = (a + b) / 2.0;
				// mass fraction of kg H2O per kg gas  (mol/m³ * m³/kg * g/mol / 1000) = kg/kg
				wH2O = c_H2OinGas / temp_density * Molweight_H2O / 1000;
				wCO2 = 1 - wH2O; // mass fraction of kg CO2 per kg gas
				// Molality of H2O in pure CO2 (mol/m³ * m³/kg) = mol/kg)
				// WW b_H2OinPureCO2 = c_H2OinGas / temp_density * 1 / wCO2;
				// new estimate of density of brine with estimated concentration of CO2 and NaCl in mol per kg pure
				// water
				Density_pureCO2 = 1000 * VLE::density_CO2(temperature, pressure);
				Density_gas = Density_pureCO2;

				// Using Norberts EOS, uses mole fractions of components, temperature and pressure [Pa]
				double m_CO2, m_H2O; // mass [kg]
				double n_CO2, n_H2O; // [mol]
				double x_CO2; // WW, x_H2O;                  // mole fraction [mol/mol]
				mass_gas = volume_gas * Density_gas; // mass of gas phase [kg] = m³ * kg/m³
				m_CO2 = mass_gas * wCO2; // mas of CO2 [kg] = kg * kg/kg
				m_H2O = mass_gas * wH2O; // mass of H2O [kg] = kg * kg/kg
				n_CO2 = m_CO2 * Molweight_CO2 / 1000; // moles of CO2 [mol] = kg * g/mol / 1000 g/kg
				n_H2O = m_H2O * Molweight_H2O / 1000; // moles of H2O [mol] = kg * g/mol / 1000 g/kg
				x_CO2 = n_CO2 / (n_CO2 + n_H2O); // mole fraction of CO2 [mol/mol]
				// WW x_H2O = n_H2O / (n_CO2 + n_H2O);      // mole fraction of H2O [mol/mol]
				// ToDo: reorganisation of the code [kg/m³]
				Density_gas = DuansMixingRule(temperature, pressure * 1e5, x_CO2, 0, 1, 0);
				// cout.precision(8);
				// cout << "Iteration: " << j << " estimated: " << temp_density << " calculated: " << Density_gas << "
				// a: " << a  << " b: " << b << "\n";

				// compare density with density estimation
				if (Density_gas < 0)
					break;
				if (fabs(Density_gas - temp_density) < epsilon)
					break;
				if (Density_gas < temp_density)
					b = temp_density;
				if (Density_gas > temp_density)
					a = temp_density;
			}

			if (Density_gas < 0)
			{
				cout << "Density calculation of gas was not possible"
				     << "\n";
				//system("Pause");
				exit(0);
			}
			// set new density to nodes
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			m_pcs->SetNodeValue(i, variable_index, Density_gas);
			// set new viscosity to nodes
			variable_index = m_pcs->GetNodeValueIndex("VISCOSITY2");
			Viscosity_gas = 5.5e-5;
			m_pcs->SetNodeValue(i, variable_index, Viscosity_gas);
		}

		// cout << "   node: " << i << " Dichte_Wasser: " << Density_liquid << " Dichte_Gas: " << Density_gas << "\n";
	}
	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;

	cout << "        Calculating Density and Viscosity " << time << " seconds."
	     << "\n";
	cout << "\n";
}

//-------------------------------------------------------------------------
// GeoSys - Function: Phase_Transition_CO2
// Task: Calculations Phasetransition for CO2 between gas and liquid phase, calls a DLL
// Return: nothing
// Programming: 11/2010 BG, NB
// Modification:
//-------------------------------------------------------------------------
void CRFProcess::Phase_Transition_CO2(CRFProcess* m_pcs, int Step)
{
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	double saturation_gas, saturation_liquid;
	double Molweight_CO2, Molweight_H2O, Molweight_NaCl;
	// ,Density_liquid_old, Density_purewater, Density_pureCO2; unused
	double Volume_eff, Density_gas, Density_liquid, p_cap;
	// , c_CO2inLiquid_old; unused
	double c_H2OinLiquid, c_CO2inLiquid, c_NaClinLiquid, c_H2OinGas, c_CO2inGas;
	// double b_CO2inPureWater, b_NaClinPureWater, b_H2OinPureCO2;
	// double mass_gas, mass_liquid; unused
	double saturation_gas_min, saturation_liquid_min, saturation_gas_effective, saturation_liquid_effective;
	double porosity, node_volume;
	int variable_index, indexProcess;
	CMediumProperties* MediaProp;
	CFluidProperties* FluidProp;
	Phase_Properties gas;
	Phase_Properties liquid;
	Phase_Properties solid;
	CRFProcess* pcs_MassTransport;
	int MassTransportID[6] = {0, 0, 0, 0, 0, 0};
	int MaterialGroup = 0;
	int TimeStepVariableIndex = 1;
	// double c_oldCO2inLiquid; unused
	// WW double p_cap_1, pressure_1;

	//---------------------------------------------------------------------------------------------------------
	// Data preprocessing

	// for initial conditions variable index is the first index, for other timesteps variable index is the second index
	if (m_pcs->Tim->step_current == 1)
		TimeStepVariableIndex = 0;

	Molweight_CO2 = 44.009; // [g/mol]
	Molweight_H2O = 18.0148; // [g/mol]
	Molweight_NaCl = 58.443; // [g/mol]				//ToDo: provide constants once in the whole project

	saturation_gas_min = 0.00001;
	saturation_liquid_min = 0.00001;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		if (i == 230)
			cout << i << "\n";
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;
		// calculate Porevolume for the node based on connected elements
		for (size_t j = 0; j < m_node->getConnectedElementIDs().size(); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];
			MaterialGroup = m_ele->GetPatchIndex();
			MediaProp = mmp_vector[MaterialGroup];
			porosity = MediaProp->Porosity(m_ele->GetIndex(), m_pcs->m_num->ls_theta);
			// ToDo: Correct calculation of node volume
			node_volume = node_volume + m_ele->GetVolume() / m_ele->GetNodesNumber(false) * porosity;
		}

		// get pressure, temperature and saturations
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		p_cap = m_pcs->GetNodeValue(i, variable_index);
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("PRESSURE1");
		// WW p_cap_1 = m_pcs->GetNodeValue(i, variable_index);
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
		// unit: bar!!! Assuming that the non wetting phase pressure is the correct one
		gas.pressure = liquid.pressure = solid.pressure = m_pcs->GetNodeValue(i, variable_index) / 1e5;
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
		// WW pressure_1 = m_pcs->GetNodeValue(i, variable_index) / 1e5;
		if (!T_Process)
		{
			// Get reference temperature if no heat transport is used
			FluidProp = MFPGet("LIQUID");
			gas.temperature = liquid.temperature = solid.temperature = FluidProp->getReferenceTemperature();
		}
		else
		{
			CRFProcess* HeatProcess = PCSGet("HEAT_TRANSPORT");
			//+1... new time level
			variable_index = HeatProcess->GetNodeValueIndex("TEMPERATURE1") + 1;
			gas.temperature = liquid.temperature = solid.temperature = m_pcs->GetNodeValue(i, variable_index);
		}
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("SATURATION1");
		saturation_liquid = m_pcs->GetNodeValue(i, variable_index);
		saturation_gas = 1 - saturation_liquid;
		// calculate new effective saturation that sum up to 1
		saturation_liquid_effective = (saturation_liquid - saturation_liquid_min)
		                              / (1 - saturation_liquid_min - saturation_gas_min);
		saturation_gas_effective = (saturation_gas - saturation_gas_min)
		                           / (1 - saturation_liquid_min - saturation_gas_min);
		if (saturation_liquid_effective < 1e-10)
		{
			saturation_liquid_effective = 0;
			saturation_gas_effective = 1;
		}
		if (saturation_gas_effective < 1e-10)
		{
			saturation_gas_effective = 0;
			saturation_liquid_effective = 1;
		}

		// calculate liquid and gas volume from saturation
		//[m³]
		Volume_eff = node_volume * (1 - saturation_gas_min - saturation_liquid_min);
		//[m³]
		liquid.volume = Volume_eff * saturation_liquid_effective;
		//[m³]
		gas.volume = Volume_eff * saturation_gas_effective;

		// cout << "node: " << i << " p_cap: " << p_cap << " P2: " << gas.pressure << " sat_liquid: " <<
		// saturation_liquid << " sat_gas: " << saturation_gas << " V_liquid: " << liquid.volume << " V_gas: " <<
		// gas.volume << "\n";

		// get ID's of mass transport processes
		for (int j = 0; j < int(pcs_vector.size()); j++)
		{
			pcs_MassTransport = pcs_vector[j];
			if (pcs_MassTransport->nod_val_name_vector[0] == "H2O")
				MassTransportID[0] = j;
			if ((pcs_MassTransport->nod_val_name_vector[0] == "C(4)")
			    || (pcs_MassTransport->nod_val_name_vector[0] == "CO2")) // "CO2_w" CB
				MassTransportID[1] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "NaCl")
				MassTransportID[2] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "CO2")
				MassTransportID[3] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "H2O_gas")
				MassTransportID[4] = j;
			if (pcs_MassTransport->nod_val_name_vector[0] == "oldC(4)")
				MassTransportID[5] = j;
		}
		// get concentration of phase components and calculate moles of components per phase
		// Read CO2 concentration in water and calculate moles of CO2 in liquid phase
		indexProcess = MassTransportID[1];
		solid.CO2 = 0;
		// +1: new timelevel
		variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
		                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
		//[mol/m³]
		c_CO2inLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);
		liquid.CO2 = c_CO2inLiquid * liquid.volume; //[mol] = mol/m³ * m³

		// Read NaCl concentration in water and calculate moles of NaCl in liquid phase
		indexProcess = MassTransportID[2];
		gas.NaCl = 0;
		solid.NaCl = 0;
		// +1: new timelevel
		variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
		                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
		//[mol/m³]
		c_NaClinLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);
		//[mol] = mol/m³ * m³
		liquid.NaCl = c_NaClinLiquid * liquid.volume;

		// Determine moles of water in liquid phase
		// get liquid density from the last step
		solid.H2O = 0;
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("DENSITY1");
		// [kg/m³]
		Density_liquid = m_pcs->GetNodeValue(i, variable_index);

		// Test 1: Berechnung T_H2O aus alter Dichte und alter CO2-Konz,
		// variable_index = m_pcs->GetNodeValueIndex("DENSITY1_old"); //+1... new time level
		// if (m_pcs->Tim->step_current == 1)
		//	Density_liquid_old = Density_liquid;
		// else
		//	Density_liquid_old = m_pcs->GetNodeValue(i, variable_index);						// [kg/m³
		////set new Density as old value
		// m_pcs->SetNodeValue(i, variable_index, Density_liquid);
		//
		// indexProcess = MassTransportID[5];
		// variable_index =
		// pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) +
		// TimeStepVariableIndex; // +1: new timelevel
		// c_CO2inLiquid_old = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);		//[mol/m³]
		// if (i == 230)
		//	cout << i << "\n";
		// liquid.H2O = (Density_liquid_old * liquid.volume - c_CO2inLiquid_old * liquid.volume  * Molweight_CO2 * 1e-3
		// - liquid.NaCl * Molweight_NaCl * 1e-3 ) / (Molweight_H2O * 1e-3);	//[mol]

		if (i == 230)
			cout << i << "\n";

		// Test 2:
		if (Step == 1)
		{
			// phase transition after reactions -> use old liquid.H2O because no phase movement
			indexProcess = MassTransportID[0];
			// +1: new timelevel
			variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
			                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
			//[mol/m³]
			c_H2OinLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);
			//[mol] = mol/m³ * m³
			liquid.H2O = c_H2OinLiquid * liquid.volume;
		}
		else
			// calculate new liquid.H2O after multiphase flow calculation based on new pressure, new density because of
			// this new pressure, new phase saturations
			//[mol]
			liquid.H2O = (Density_liquid * liquid.volume - liquid.CO2 * Molweight_CO2 * 1e-3
			              - liquid.NaCl * Molweight_NaCl * 1e-3) / (Molweight_H2O * 1e-3);

		// if (m_pcs->Tim->step_current == 1)
		// liquid.H2O = (Density_liquid * liquid.volume - liquid.CO2 * Molweight_CO2 * 1e-3 - liquid.NaCl *
		// Molweight_NaCl * 1e-3 ) / (Molweight_H2O * 1e-3);	//[mol]
		// else {
		// indexProcess = MassTransportID[0];
		// variable_index =
		// pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1; //
		// +1: new timelevel
		// c_H2OinLiquid = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);		//[mol/m³]
		// liquid.H2O = c_H2OinLiquid * liquid.volume;										//[mol] = mol/m³ * m³
		//}

		// components in gas
		// Read H2O concentration in gas and calculate moles of H2O in gas phase
		indexProcess = MassTransportID[4];
		// +1: new timelevel
		variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
		                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + TimeStepVariableIndex;
		//[mol/m³]
		c_H2OinGas = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);
		gas.H2O = c_H2OinGas * gas.volume; //[mol] = mol/m³ * m³

		// Determine moles of CO2 in gas phase
		// get gas density from the last step
		//+1... new time level
		variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
		//[kg/m³]
		Density_gas = m_pcs->GetNodeValue(i, variable_index);

		// Test 2:
		if (Step == 1)
		{
			// phase transition after reactions -> use old liquid.H2O because no phase movement
			indexProcess = MassTransportID[3];
			// +1: new timelevel
			variable_index = pcs_vector[indexProcess]->GetNodeValueIndex(
			                     pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
			//[mol/m³]
			c_CO2inGas = pcs_vector[indexProcess]->GetNodeValue(i, variable_index);
			gas.CO2 = c_CO2inGas * gas.volume; //[mol] = mol/m³ * m³
		}
		else
			// calculate new liquid.H2O after multiphase flow calculation based on new pressure, new density because of
			// this new pressure, new phase saturations
			// [mol] = kg / (kg / mol)
			gas.CO2 = (Density_gas * gas.volume - gas.H2O * Molweight_H2O * 1e-3) / (Molweight_CO2 * 1e-3);

		//---------------------------------------------------------------------------------------------------------
		// Calculating phase transition with Dedong's method

		// recalculate units
		// kg/m³ -> g/cm³
		gas.density /= 1000;
		liquid.density /= 1000;
		solid.density /= 1000;
		// m³ -> cm³
		gas.volume *= 1e6;
		liquid.volume *= 1e6;
		solid.volume *= 1e6;
		// kg -> g
		gas.mass *= 1000;
		liquid.mass *= 1000;
		solid.mass *= 1000;

		if (m_pcs->Tim->step_current == 14)
			cout << i << "\n";

		if (i == 230)
			cout << i << " gas pressure: " << gas.pressure << " liquid pressure: " << liquid.pressure << "\n";

		// call Dedongs method to calculate CO2 phase transition
		CO2_H2O_NaCl_VLE_isochoric(gas, liquid, solid, 1);

		if (i == 230)
			cout << i << " gas pressure: " << gas.pressure << " liquid pressure: " << liquid.pressure << "\n";

		// g/cm³ -> kg/m³
		gas.density *= 1000;
		liquid.density *= 1000;
		solid.density *= 1000;
		// cm³ -> m³
		gas.volume /= 1e6;
		liquid.volume /= 1e6;
		solid.volume /= 1e6;
		// g -> kg
		gas.mass /= 1000;
		liquid.mass /= 1000;
		solid.mass /= 1000;
		// bar -> Pa
		gas.pressure *= 1e5;
		liquid.pressure *= 1e5;
		solid.pressure *= 1e5;

		//---------------------------------------------------------------------------------------------------------
		// Postprocessing of results to recalculate volume, pressure, saturation

		// calculate saturation of phases
		// double error = fabs((gas.volume + liquid.volume) - Volume_eff);
		// error = error / Volume_eff;
		if (fabs((gas.volume + liquid.volume) - Volume_eff) / Volume_eff > 1e-3)
		{
			cout << "The volume is not equal before and after the calculation of CO2 phase transition! "
			     << "\n";
			cout << "Before: " << Volume_eff << " After: " << gas.volume + liquid.volume << "\n";
			//system("Pause");
			exit(0);
		}

		// set new densities to nodes
		if (liquid.density >= 0)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY1");
			m_pcs->SetNodeValue(i, variable_index, liquid.density);
		}
		if (gas.density >= 0)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			m_pcs->SetNodeValue(i, variable_index, gas.density);
		}

		// set new viscosities to nodes
		if (liquid.viscosity >= 0)
		{
			variable_index = m_pcs->GetNodeValueIndex("VISCOSITY1");
			m_pcs->SetNodeValue(i, variable_index, liquid.viscosity);
		}
		if (gas.viscosity >= 0)
		{
			variable_index = m_pcs->GetNodeValueIndex("VISCOSITY2");
			m_pcs->SetNodeValue(i, variable_index, gas.viscosity);
		}

		// ToDo: check calculation
		saturation_gas = max(gas.volume / Volume_eff, saturation_gas_min);
		// saturation_liquid = max(liquid.volume / Volume_eff, saturation_liquid_min);
		saturation_liquid = max(1 - saturation_gas, saturation_liquid_min);

		// set new saturations to nodes
		if (Step == 1)
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("SATURATION1");
		else
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
		m_pcs->SetNodeValue(i, variable_index, saturation_liquid);

		// new calculation of capillary pressue

		// ToDo: works only for homogenouse properties
		MediaProp = mmp_vector[MaterialGroup];
		p_cap = MediaProp->CapillaryPressureFunction(saturation_liquid);

		if (i == 230)
			cout << i << " gas pressure: " << gas.pressure << " liquid pressure: " << liquid.pressure << "\n";

		// set new pressure to nodes
		if (Step == 1)
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE1");
		else
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		m_pcs->SetNodeValue(i, variable_index, p_cap);
		if (Step == 1)
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
		else
			//+1... new time level
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
		m_pcs->SetNodeValue(i, variable_index, gas.pressure);
		// cout << "node: " << i << " p_cap: " << p_cap << " P2: " << gas.pressure << " sat_liquid: " <<
		// saturation_liquid << " sat_gas: " << saturation_gas << " V_liquid: " << liquid.volume << " V_gas: " <<
		// gas.volume << "\n";

		// set new component concentration
		if (liquid.volume > 0)
		{
			//[mol/m³] = mol / m³
			c_H2OinLiquid = liquid.H2O / liquid.volume;
			//[mol/L] = mol / m³
			c_CO2inLiquid = liquid.CO2 / liquid.volume;
			//[mol/L] = mol / m³
			c_NaClinLiquid = liquid.NaCl / liquid.volume;
		}
		else
		{
			c_H2OinLiquid = 0;
			c_CO2inLiquid = 0;
			c_NaClinLiquid = 0;
		}

		indexProcess = MassTransportID[0]; // H2O in liquid
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_H2OinLiquid);

		indexProcess = MassTransportID[1]; // CO2 in liquid
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_CO2inLiquid);

		indexProcess = MassTransportID[2]; // NaCl in liquid
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_NaClinLiquid);

		if (gas.volume > 0)
		{
			c_CO2inGas = gas.CO2 / gas.volume; //[mol/L] = mol / m³
			c_H2OinGas = gas.H2O / gas.volume; //[mol/L] = mol / m³
		}
		else
		{
			c_CO2inGas = 0;
			c_H2OinGas = 0;
		}
		indexProcess = MassTransportID[3]; // CO2 in gas
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_CO2inGas);

		indexProcess = MassTransportID[4]; // H2O in gas
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_H2OinGas);

		indexProcess = MassTransportID[5]; // old CO2 in liquid
		// +1: new timelevel
		variable_index
		    = pcs_vector[indexProcess]->GetNodeValueIndex(pcs_vector[indexProcess]->pcs_primary_function_name[0]) + 1;
		pcs_vector[indexProcess]->SetNodeValue(i, variable_index, c_CO2inLiquid);
	}
}
//-------------------END------------------

// CB _ctx_

/***************************************************************************
    Calculate electric field from charges of ions present
    No external field included jet
    04/2010     CB    Implementation
****************************************************************************/

// CB _ctx_ CB_merge_0513 ??

/*
void Calc_ctx_Contribution(void){

cout << "      Calc_ctx_Contribution()" << "\n";

CRFProcess *m_pcs = NULL, *m_pcs2=NULL;
long nnodes, nelements, i, gp, k, n_gauss=0;
Mesh_Group::CElem *elem = NULL;
CFiniteElementStd* m_fem = NULL, *m_fem2=NULL;


double sum_C, sum_gradC[3], sum_vC[3];
double vel[3];
int gp_r=0, gp_s=0, gp_t=0;
double coef = 0.0, fkt=0.0;
int i_idx;
double Conc, valence, diff_coeff, gradConc[3], S_gradConc[3];
double *node_concs, conc, conc_help;
node_concs = new double [20];
int ele_dim, i_dim;

string _ctx_SpeciesName = "Substrate";
string BacteriaName = "Bacteria";
double X_s;

ElementValue* gp_ele = NULL; // Get gauss point values

//nnodes = (long) fem_msh_vector[0]->nod_vector.size();
nelements = (long) fem_msh_vector[0]->ele_vector.size();
// get number of gauss points


// get the bacterium process X
for(k=0;k<cp_vec.size();k++){
  if(cp_vec[k]->compname.compare(BacteriaName)==0)
    break;
}
m_pcs = PCSGet("MASS_TRANSPORT", cp_vec[k]->compname);
m_fem = m_pcs->GetAssember();
n_gauss = m_fem->GetNumGaussPoints();



X_s = cp_vec[k]->_ctx_Coefficient;

// get the substrate process x_s*grad(C)
for(k=0;k<cp_vec.size();k++){
  if(cp_vec[k]->compname.compare(_ctx_SpeciesName)==0)
    break;
}
m_pcs2 = PCSGet("MASS_TRANSPORT", _ctx_SpeciesName);

// element loop
for(i = 0; i < nelements; i++) {
  // get each element
  //if(i<5) cout << "\n" << " ele: " << i << "\n";
  elem = m_pcs->m_msh->ele_vector[i];
  ele_dim = elem->GetDimension();
  nnodes = elem->nodes_index.Size();

  if (elem->GetMark()){
    // for each element calculate grad(Substrate)
    // get number of gauss points for this element
    m_fem = m_pcs->GetAssember();
    m_fem->ConfigElement(elem);
    n_gauss = m_fem->GetNumGaussPoints();
    gp_ele = ele_gp_value[i];

    // gauss point loop
    for(gp=0;gp<n_gauss;gp++){
      // if(i<5) cout << "\n"<< "gauss point " <<  gp << "\n";
      sum_C = 0.0;
      for(i_dim=0;i_dim < ele_dim;i_dim++){
          sum_gradC[i_dim] = 0.0;
      }

      //gp_ele->getIPvalue_vec(gp,vel);
      //for(i_dim=0;i_dim < ele_dim;i_dim++) cout << " i_dim " << i_dim << ", vel:   " << vel[i_dim] << "\n";

      i_idx = 1; // as this is in the end of transport loop,  take c_new

      fkt = m_fem->GetGaussData(gp, gp_r, gp_s, gp_t);
      m_fem->ComputeShapefct(1);
      m_fem->ComputeGradShapefct(1);

        // Get C of substrate at Meshnodes of this ele
        for(int l=0; l<nnodes; l++){
          conc_help = m_pcs2->GetNodeValue(elem->GetNodeIndex(l), i_idx);
          node_concs[l] = fabs(conc_help);
          //if(i<5) cout << " l, node_conc[l]:  " << l << "   " << node_concs[l] << "\n";
      }

      // interpolate conc to GP
      //Conc = m_fem->interpolate(node_concs);
      //if(i< 5) cout << " Conc[gp]: " << Conc << "     " << "\n";

      // calc grad_C
      for(i_dim=0; i_dim < ele_dim; i_dim++)
        gradConc[i_dim] = 0.0;
      for(i_dim=0; i_dim < ele_dim; i_dim++){
        for(int l=0; l<nnodes; l++){
            gradConc[i_dim] += node_concs[l]* m_fem->Getdshapefct(i_dim*nnodes+l);
            //if(i< 5) cout << "i_dim, l, gradConc, node_concs, shapefunction : " << i_dim << ",  " << l ;
            //if(i< 5) cout << " ,  " << gradConc[i_dim] << ", " << node_concs[l] << " ,  " <<
m_fem->Getdshapefct(i_dim*m_fem->Getnnodes()+l) << "\n";
        }
      }


      // save result for element in gauss point
      for(i_dim=0; i_dim < ele_dim; i_dim++){
        S_gradConc[i_dim] = 0.0;
        //if(fabs(sum_C) > 1.0e-15)
        S_gradConc[i_dim] = -1.0 * X_s *  gradConc[i_dim] ;
        //if(i<5) cout << el_field[i_dim] << "  " << sum_gradC[i_dim] << "  " << sum_vC[i_dim] << "  "  << sum_C <<
"\n";

        m_fem->Set_ctx_(i,S_gradConc[i_dim], gp, i_dim);

        //gp_ele->_ctx_Gauss(i_dim,gp) = el_field[i_dim];
        //if(i<500) cout  << "   (" << gp << ", " << el_field[i_dim] << ") " << "\n";
        //cout << " Test: GetElectricField2: " << m_fem->GetElField(i,gp,i_dim) << "\n";
      }

    } // end gauss point loop
  }
} // end element loop
//cout << "\n" << "\n";
delete node_concs;

}

*/
