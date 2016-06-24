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
   Task:
   Programing:
   02/2003 OK Implementation
   11/2004 OK PCS2
**************************************************************************/
#ifndef rf_pcs_INC
#define rf_pcs_INC

#include <ctime>
#include <valarray>

#include "makros.h"

// MSHLib
#include "msh_lib.h"

// PCSLib
#include "ProcessInfo.h"
#include "rf_bc_new.h"
#include "rf_num_new.h"
#include "rf_tim_new.h"
#include "conversion_rate.h" // HS, 10.2011
#include "SparseMatrixDOK.h"

#include "Eigen/Eigen"

//#include "rf_st_new.h"//CMCD 02_06
// C++ STL
//#include <fstream>
//
// The follows are implicit declaration. WW
//---------------------------------------------------------------------------
#if defined(USE_PETSC) // || defined(using other parallel scheme)
namespace petsc_group
{
class PETScLinearSolver;
}
#endif

namespace FiniteElement
{
class CFiniteElementStd;
class CFiniteElementVec;
class ElementMatrix;
class ElementValue;
}

namespace MeshLib
{
class CFEMesh;
}

#ifdef NEW_EQS // WW
namespace Math_Group
{
class Linear_EQS;
}
using Math_Group::Linear_EQS;
#endif

// using Math_Group::SparseMatrixDOK;
//

/*
class etr_data
{
public:
    double x;
    double y;
    double val;
};
*/

class CSourceTermGroup;
class CSourceTerm;
class CNodeValue;
class Problem; // WW
class CECLIPSEData; // BG Coupling to Eclipse
class CDUMUXData; // BG Coupling to DuMux
class CPlaneEquation; // BG Calculating Plane Equation
using FiniteElement::CFiniteElementStd;
using FiniteElement::CFiniteElementVec;
using FiniteElement::ElementMatrix;
using FiniteElement::ElementValue;
using MeshLib::CFEMesh;
//---------------------------------------------------------------------------

#define PCS_FILE_EXTENSION ".pcs"

typedef struct /* Knotenwert-Informationen */
{
	char name[80]; /* Name der Knotengroesse */
	char einheit[10]; /* Einheit der Knotengroesse */
	int speichern; /* s.u., wird Wert gespeichert ? */
	int laden; /* s.u., wird Wert zu Beginn geladen ? */
	int restart; /* s.u., wird Wert bei Restart geladen ? */
	int adapt_interpol; /* Soll Wert bei Adaption auf Kinder interpoliert werden? */
	double vorgabe; /* Default-Wert fuer Initialisierung */
	int nval_index;
	int pcs_this;
	int timelevel;
} PCS_NVAL_DATA;

typedef struct /* element data info */
{
	char name[80]; /* Name der Elementgroesse */
	char einheit[10]; /* Einheit der Elementgroesse */
	int speichern; /* s.u., wird Wert gespeichert ? */
	int laden; /* s.u., wird Wert zu Beginn geladen ? */
	int restart; /* s.u., wird Wert bei Restart geladen ? */
	int adapt_interpol; /* Soll Wert bei Adaption auf Kinder vererbt werden? */
	double vorgabe; /* Default-Wert fuer Initialisierung */
	int eval_index;
	int index;
} PCS_EVAL_DATA;

typedef struct
{
	long index_node;
	double water_st_value;
} Water_ST_GEMS; // HS 11.2008

typedef struct // CB DL CO2 phase transition
{
	std::string name; // fluid name
	double temperature;
	double pressure;
	double density; // density g/cm^3
	double viscosity; // viscosity mPa.s
	double volume; // volume cm^3
	double mass; // weight g
	double CO2; // mole of CO2
	double H2O; // mole of H2O
	double NaCl; // mole of NaCl
	double C_GAS; // mole of CO2 in Gas
	double H2; // mole of H2
} Phase_Properties;

typedef struct
{
	double B;
	double C;
	double D;
	double E;
	double F;
	double b;
	double G;
	double Vc;
	double V2, V3, V5, V6;
	double Tc, Pc;
	double M;
	int id;
} VirialCoefficients;

#ifdef JFNK_H2M
/// Dirchlet BC at node. 09.2010. WW
typedef struct
{
	int var_idx;
	bool incremental;
	long bc_eqs_idx;
	long bc_node;
	double bc_value;
	double bc_value0;
} bc_JFNK;
#endif

// MB moved inside the Process object
// extern vector<double*>nod_val_vector; //OK
// extern vector<string>nod_val_name_vector; //OK
// extern void pcs->SetNodeValue(long,int,double); //OK
// extern double pcs->GetNodeValue(long,int); //OK
// extern int pcs->GetNodeValueIndex(string); //OK

// Element values for all process
// Moved inside Process object
// extern vector<double*>ele_val_vector; //PCH
// extern vector<string>eld_val_name_vector; //PCH
// extern void SetElementValue(long,int,double); //PCH
// extern double GetElementValue(long,int); //PCH
// extern int GetElementValueIndex(string); //PCH
/**************************************************************************
   FEMLib-Class:
   Task:
   Programing:
   03/2003 OK Implementation
   02/2005 WW Local element assembly (all protected members)
   12/2005 OK MSH_TYPE
   last modification:
   05/2010 TF inserted pointer to instance of class Problem inclusive access functions
**************************************************************************/

/**
 * class manages the physical processes
 */
class CRFProcess : public ProcessInfo
{
	//----------------------------------------------------------------------
	// Properties
private:
	/**
	 * _problem is a pointer to an instance of class Problem.
	 *  The pointer is used to get the geometric entities.
	 */
	Problem* _problem;

	void VariableStaticProblem();
	void VariableDynamics();
	bool compute_domain_face_normal; // WW
	int continuum;
	bool continuum_ic;

	bool isRSM; // WW
	double* eqs_x; //> Pointer to x array of eqs (added due to PETSC)

	std::vector<std::string> pcs_type_name_vector;
	bool _hasConstrainedBC;
	bool _hasConstrainedST;
	long _idxVx, _idxVy, _idxVz;

protected: // WW
	friend class FiniteElement::CFiniteElementStd;
	friend class FiniteElement::CFiniteElementVec;
	friend class FiniteElement::ElementValue;
	friend class ::CSourceTermGroup;
	friend class ::Problem;
	friend class CECLIPSEData; // BG
	friend class CDUMUXData; // SBG
	/// Number of nodes to a primary variable. 11.08.2010. WW
	int* p_var_index;
	long* num_nodes_p_var;

	/// Size of unknowns. 02.2011. WW
	long size_unknowns;
	// Assembler
	CFiniteElementStd* fem;
	// Time step control
	bool accepted; // 25.08.1008. WW
	int accept_steps; // 27.08.1008. WW
	int reject_steps; // 27.08.1008. WW
	//
	int dof; // WW
	long orig_size; // Size of source term nodes
	// ELE
	std::vector<FiniteElement::ElementMatrix*> Ele_Matrices;
	// Global matrix
	Math_Group::Vec* Gl_Vec; // NW
	Math_Group::Vec* Gl_Vec1; // NW
	Math_Group::Vec* Gl_ML; // NW
	Math_Group::SparseMatrixDOK* FCT_AFlux; // NW
#ifdef USE_PETSC
	Math_Group::SparseMatrixDOK* FCT_K;
	Math_Group::SparseMatrixDOK* FCT_d;
#endif
	/**
	 * Storage type for all element matrices and vectors
	 * Cases:
	 * 0. Do not keep them in the memory
	 * 1. Keep them to vector Ele_Matrices
	 */
	int Memory_Type;
	//....................................................................
	int additioanl2ndvar_print; // WW
	// TIM
	friend class CTimeDiscretization;
	CTimeDiscretization* Tim; // time
	bool femFCTmode; // NW
	void CopyU_n(); // 29.08.2008. WW
	// Time unit factor
	double time_unit_factor;
	int NumDeactivated_SubDomains;
	int* Deactivated_SubDomain;
// New equation and solver objects WW
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	petsc_group::PETScLinearSolver* eqs_new;
	int mysize;
	int myrank;
#elif defined(NEW_EQS)
#if defined(LIS) || defined(MKL)
public:
	Linear_EQS* eqs_new;
#else
	Linear_EQS* eqs_new;
#endif // LIS endif for Fluid Momentum	// PCH
	bool configured_in_nonlinearloop;
#else
	LINEAR_SOLVER* eqs;
#endif

//
#if defined(USE_MPI) || defined(USE_PETSC) // WW
	clock_t cpu_time_assembly;
#endif
	// Position of unknowns from different DOFs in the system equation
	//....................................................................
	// OUT
	// Write indices of the nodes with boundary conditons
	bool write_boundary_condition; // 15.01.2008. WW
	bool OutputMassOfGasInModel; // BG 05/2012
	// Element matrices output
	void Def_Variable_LiquidFlow();
	void Def_Variable_MultiPhaseFlow();
	bool Write_Matrix;
	std::fstream* matrix_file;
	// Write RHS from source or Neumann BC terms or BC to file
	// 0: Do nothing
	// 1: Write
	// 2: Read
	int WriteSourceNBC_RHS;
	int WriteProcessed_BC;
	// Write the current solutions/Read the previous solutions WW
	// -1, default. Do nothing
	// 1. Write
	// 2. Read
	// 3 read and write
	int reload;
	long nwrite_restart;
	void WriteRHS_of_ST_NeumannBC();
	void ReadRHS_of_ST_NeumannBC();
	void Write_Processed_BC(); // 05.08.2011. WW
	void Read_Processed_BC(); // 05.08.2011. WW

	friend bool PCSRead(std::string);
	//....................................................................
	// 1-GEO
	int ite_steps; /// Newton step index;
#ifdef JFNK_H2M

	/// Norn of F for Newton method
	double norm_F0;
	/// Buffer array for JFNK, 13.08.2010. WW
	bool JFNK_precond;
	double* norm_u_JFNK;
	double* array_u_JFNK;
	double* array_Fu_JFNK;
	std::vector<bc_JFNK> BC_JFNK;
#endif
public:
	// BG, DL Calculate phase transition of CO2
	void CO2_H2O_NaCl_VLE_isobaric(double T, double P, Phase_Properties& vapor, Phase_Properties& liquid,
	                               Phase_Properties& solid, int f);
	// BG, DL Calculate phase transition of CO2
	void CO2_H2O_NaCl_VLE_isochoric(Phase_Properties& vapor, Phase_Properties& liquid, Phase_Properties& solid, int f);
	// BG, NB Calculate phase transition of CO2
	void Phase_Transition_CO2(CRFProcess* m_pcs, int Step);
	int Phase_Transition_Model; // BG, NB flag of Phase_Transition_Model (1...CO2-H2O-NaCl)
	// BG 11/2010 Sets the initial conditions for multi phase flow if Phase_Transition_CO2 is used
	void CalculateFluidDensitiesAndViscositiesAtNodes(CRFProcess* m_pcs);
	/**
	 * Sets the value for pointer _problem.
	 * @param problem the value for _problem
	 */
	void setProblemObjectPointer(Problem* problem);

	/**
	 * get access to the instance of class Problem
	 * @return
	 */
	Problem* getProblemObjectPointer() const;
	std::string geo_type; // OK
	std::string geo_type_name; // OK
	//....................................................................
	// 2-MSH
	//....................................................................
	// 3-TIM
	//....................................................................
	// 4-IC
	//....................................................................
	// 5-BC
	// WW
	std::vector<CBoundaryConditionNode*> bc_node_value;
	std::vector<CBoundaryCondition*> bc_node; // WW
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	std::vector<long> bc_node_value_in_dom; // WW for domain decomposition
	std::vector<long> bc_local_index_in_dom; // WW for domain decomposition
	std::vector<long> rank_bc_node_value_in_dom; // WW
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	std::vector<long> bc_transient_index; // WW/CB
	std::vector<long> st_transient_index; // WW/CB...BG
	void UpdateTransientBC(); // WW/CB
	void UpdateTransientST(); // WW/CB...BG
	//....................................................................
	// 6-ST
	// Node values from sourse/sink or Neumann BC. WW
	std::vector<CNodeValue*> st_node_value; // WW
	std::vector<CSourceTerm*> st_node; // WW
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	std::vector<long> st_node_value_in_dom; // WW for domain decomposition
	std::vector<long> st_local_index_in_dom; // WW for domain decomposition
	std::vector<long> rank_st_node_value_in_dom; // WW
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	void RecordNodeVSize(const int Size) // WW
	{
		orig_size = Size;
	}
	int GetOrigNodeVSize() const // WW
	{
		return orig_size;
	}

	FiniteElement::CFiniteElementStd* getLinearFEMAssembler() { return fem; }

	//....................................................................
	// 7-MFP
	//....................................................................
	// 8-MSP
	//....................................................................
	// 9-MMP
	int GetContinnumType() const { return continuum; }
	// const int number_continuum=1;
	std::vector<double> continuum_vector;
	//....................................................................
	// 10-MCP
	//....................................................................
	// 11-OUT
	void WriteSolution(); // WW
	void ReadSolution(); // WW
	//....................................................................
	// 12-NUM
	//....................................................................
	// 13-ELE
	//----------------------------------------------------------------------
	// Methods
	//....................................................................
	// Construction / destruction
	CRFProcess(void);
	void Create(void);
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	void CreateNew(void);
	void CreateFDMProcess();
	void DestroyFDMProcess();
#endif
	virtual ~CRFProcess();
	//....................................................................
	// IO
	std::ios::pos_type Read(std::ifstream*);
	void PCSReadConfigurations();
	void Write(std::fstream*);
	//....................................................................
	// 1-GEO
	//....................................................................
	// 2-MSH
	//....................................................................
	// 3-TIM
	//....................................................................
	// 4-IC
	//....................................................................
	// 5-BC
	void CreateBCGroup();
	void SetBC(); // OK
	void WriteBC(); // 15.01.2008. WW
	//....................................................................
	// 6-ST
	void CreateSTGroup();

	void SetBoundaryConditionAndSourceTerm();

	//....................................................................
	// 7-MFP
	//....................................................................
	// 8-MSP
	//....................................................................
	// 9-MMP
	//....................................................................
	// 10-MCP
	//....................................................................
	// 11-OUT
	void WriteAllVariables(); // OK
	//....................................................................
	// 12-NUM
	//....................................................................
	// 13-ELE
	//....................................................................
	// 14-CPL
	void SetCPL(); // OK8 OK4310
	//....................................................................
	// 15-EQS
	// WW void AssembleParabolicEquationRHSVector(CNode*); //(vector<long>&); //OK
	//....................................................................
	int Shift[10];
	// 16-GEM  // HS 11.2008
	int flag_couple_GEMS; // 0-no couple; 1-coupled
	std::vector<Water_ST_GEMS> Water_ST_vec;
	std::vector<long> stgem_node_value_in_dom; // KG for domain decomposition
	std::vector<long> stgem_local_index_in_dom; // KG for domain decomposition
	// KG
	std::vector<long> rank_stgem_node_value_in_dom;

	void Clean_Water_ST_vec(void);
	void Add_GEMS_Water_ST(long idx, double val);
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	void SetSTWaterGemSubDomain(int myrank);
#endif
	// ECLIPSE interface:
	CDUMUXData* DuMuxData; // SBG
	CECLIPSEData* EclipseData; // BG
	void CalGPVelocitiesfromECLIPSE(std::string path, int timestep, int phase_index, std::string phase);
	std::string simulator; // which solver to use, i.e. GeoSys, ECLIPSE or DuMux
	std::string simulator_path; // path for executable of external simulator
	std::string simulator_model_path; // path to exclipse input data file (*.data), with extension
	bool PrecalculatedFiles; // defines if Eclipse or dumux is calculated or if precalculated files are used
	bool SaveEclipseDataFiles; // WTP: save Eclipse data input files for benchmarking on systems with no Eclipse
	                           // licenses; use in combination with Precalculated Files
	std::string simulator_well_path; // path to well schedule ( *.well), with extension
	// SB redo wtp
	std::string dissolved_co2_pcs_name; // Keyword DISSOLVED_CO2_PCS_NAME; Name of MASS_TRANSPORT Process which is used
	                                    // to store total dissolved CO2 from ECLIPSE
	std::string dissolved_co2_ingas_pcs_name;

	//....................................................................
	// Construction / destruction
	char pcs_name[MAX_ZEILE]; // string pcs_name;
	int pcs_number;
	int mobile_nodes_flag;

	int pcs_type_number;
	int type;
	int GetObjType() const { return type; }
	int pcs_component_number; // SB: counter for transport components
	int ML_Cap; // 23.01 2009 PCH
	int PartialPS; // 16.02 2009 PCH
	//
	// JT2012: Process type identifiers
	bool isPCSFlow;
	bool isPCSMultiFlow;
	bool isPCSHeat;
	bool isPCSMass;
	bool isPCSDeformation;

	int GetProcessComponentNumber() const // SB:namepatch
	{
		return pcs_component_number;
	}
	std::string file_name_base; // OK
	// Access to PCS
	// Configuration 1 - NOD
	PCS_NVAL_DATA* pcs_nval_data; /// OK
	int number_of_nvals;
	int pcs_number_of_primary_nvals;
	size_t GetPrimaryVNumber() const { return static_cast<size_t>(pcs_number_of_primary_nvals); }
	const char* pcs_primary_function_unit[7];
	const char* pcs_primary_function_name[7];
	const char* GetPrimaryVName(const int index) const { return pcs_primary_function_name[index]; }
	std::string primary_variable_name; // OK
	int pcs_number_of_secondary_nvals;
	size_t GetSecondaryVNumber() const { return static_cast<size_t>(pcs_number_of_secondary_nvals); }
	const char* pcs_secondary_function_name[PCS_NUMBER_MAX];
	const char* GetSecondaryVName(const int index) const { return pcs_secondary_function_name[index]; }
	const char* pcs_secondary_function_unit[PCS_NUMBER_MAX];
	int pcs_secondary_function_timelevel[PCS_NUMBER_MAX];
	int pcs_number_of_history_values;
	/*double pcs_secondary_function_time_history[PCS_NUMBER_MAX];//CMCD for analytical solution
	   double pcs_secondary_function_value_history[PCS_NUMBER_MAX];//CMCD for analytical solution
	   void Set_secondary_function_time_history(const int index, double value)
	   {pcs_secondary_function_time_history[index]=value;}//CMCD for analytical solution
	   void Set_secondary_function_value_history(const int index, double value)
	   {pcs_secondary_function_value_history[index]=value;}//CMCD for analytical solution
	   double Get_secondary_function_time_history(const int index){return
	   pcs_secondary_function_time_history[index];}//CMCD for analytical solution
	   double Get_secondary_function_value_history(const int index){return
	   pcs_secondary_function_value_history[index];}//CMCD for analytical solution*/
	// Configuration 2 - ELE
	PCS_EVAL_DATA* pcs_eval_data;
	int pcs_number_of_evals;
	const char* pcs_eval_name[PCS_NUMBER_MAX];
	const char* pcs_eval_unit[PCS_NUMBER_MAX];
	// Configuration 3 - ELE matrices
	// Execution
	// NUM
	std::string num_type_name;
	int rwpt_app;
	int rwpt_count; // YS 05.2013 Count the number of particles.
	int srand_seed;
	const char* pcs_num_name[2]; // For monolithic scheme
	TimType::type tim_type;
	const char* pcs_sol_name;
	std::string cpl_type_name;
	CNumerics* m_num;
	//
	bool selected; // OK
	bool saturation_switch; // JOD
	// MSH
	CFEMesh* m_msh; // OK
	std::string msh_type_name; // OK
	// MB-------------
	std::vector<double*> nod_val_vector; // OK
	// OK
	std::vector<std::string> nod_val_name_vector;
	void SetNodeValue(long, int, double); // OK
	double GetNodeValue(size_t, int); // OK
	double* getNodeValue_per_Variable(const int entry_id) const // WW
	{
		return nod_val_vector[entry_id];
	}
	int GetNodeValueIndex(const std::string&, bool reverse_order = false); // OK
	//-----------------------------

	std::vector<std::string> const& getElementValueNameVector() { return ele_val_name_vector; }

private:
	// PCH
	std::vector<std::string> ele_val_name_vector;

public:
	std::vector<double*> ele_val_vector; // PCH
	void SetElementValue(long, int, double); // PCH
	double GetElementValue(size_t, int); // PCH
	// PCH
	int GetElementValueIndex(const std::string&, bool reverse_order = false);
	// CB-----------------------------
	int flow_pcs_type;
	//----------------------------------------------------------------------
	// Methods
	// Access to PCS
	CRFProcess* GetProcessByFunctionName(char* name);
	CRFProcess* GetProcessByNumber(int);
	CFiniteElementStd* GetAssembler() { return fem; }
	int GetDOF() { return dof; }
	bool pcs_is_cpl_overlord;
	bool pcs_is_cpl_underling;
	CRFProcess* cpl_overlord;
	CRFProcess* cpl_underling;
	// CRFProcess *Get(string); // WW Removed
	// Configuration
	void Config();
	void ConfigGroundwaterFlow();
	void ConfigLiquidFlow();
	void ConfigNonIsothermalFlow();
	void ConfigNonIsothermalFlowRichards();
	void ConfigMassTransport();
	void ConfigHeatTransport();
	void ConfigDeformation();
	void ConfigMultiphaseFlow();
	void ConfigGasFlow();
	void ConfigUnsaturatedFlow(); // OK4104
	void ConfigFluidMomentum();
	void ConfigRandomWalk();
	void ConfigMultiPhaseFlow();
	void ConfigPS_Global(); // PCH
	void ConfigMULTI_COMPONENTIAL_FLOW(); // AKS/NB
	void ConfigTNEQ(); // HS,TN
	void ConfigTES(); // HS,TN
// Configuration 1 - NOD
#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
	void setSolver(petsc_group::PETScLinearSolver* petsc_solver);
	double CalcIterationNODError(int method); // OK // PETSC version in rf_pcs1.cpp WW
#endif

#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	void ConfigNODValues1(void);
	void ConfigNODValues2(void);
	void CreateNODValues(void);
	void SetNODValues(); // OK
	void CalcFluxesForCoupling(); // MB
	void SetNODFlux(); // OK
	//
	void AssembleParabolicEquationRHSVector(); // OK
	// 15-EQS
	//(vector<long>&); //OK
	void AssembleParabolicEquationRHSVector(MeshLib::CNode*);
#endif
	double CalcIterationNODError(FiniteElement::ErrorMethod method, bool nls_error, bool cpl_error = false); // OK
	// Add bool forward = true. WW
	void CopyTimestepNODValues(bool forward = true);
	void CopyTimestepELEValues(bool forward = true);
	// Coupling
	// WW double CalcCouplingNODError(); //MB
	void CopyCouplingNODValues();
	// WWvoid CalcFluxesForCoupling(); //MB
	// Configuration 2 - ELE
	void ConfigELEValues1(void);
	void ConfigELEValues2(void);
	void CreateELEValues(void);
	void CreateELEGPValues();
	void AllocateMemGPoint(); // NEW Gauss point values for CFEMSH WW
	void CalcELEVelocities(void);
	// void CalcELEMassFluxes();				//BG
	// WW   double GetELEValue(long index,double*gp,double theta,string nod_fct_name);
	void CheckMarkedElement(); // WW
	void CheckExcavedElement(); // WX
	// Configuration 3 - ELE matrices
	void CreateELEMatricesPointer(void);
	// Equation system
	//---WW
	CFiniteElementStd* GetAssember() { return fem; }
	void AllocateLocalMatrixMemory();
	virtual void GlobalAssembly(); // Make as a virtul function. //10.09.201l. WW
	/// For all PDEs excluding that for deformation. 24.11.2010l. WW
	void GlobalAssembly_std(const bool is_mixed_order, bool Check2D3D = false);
	/// Assemble EQS for deformation process.
	virtual void GlobalAssembly_DM(){};
#if defined(NEW_EQS) && defined(JFNK_H2M)
	/// Jacobian free methid to calculate J*v.
	// 11.08.2010.
	void Jacobian_Multi_Vector_JFNK(double* v = NULL, double* Jv = NULL);
	/// Recovery du from the temporary vector.
	void Recovery_du_JFNK(); // 02.11.2010.
	/// Line serach for Newton method.
	double LineSearch(); // 10.12.2010.
	/// Force term control for inexact Newton method. 01.2011.
	bool ForceTermCriterion(double* Jdx, const int num_iteration);
#endif
	void AddFCT_CorrectionVector(); // NW
	void ConfigureCouplingForLocalAssemblier();
	void CalIntegrationPointValue();
	bool cal_integration_point_value; // WW
	void CalGPVelocitiesfromFluidMomentum(); // SB 4900
	bool use_velocities_for_transport; // SB4900

	//---
	double Execute();
	double ExecuteNonLinear(int loop_process_number, bool print_pcs = true);
	void PrintStandardIterationInformation(bool write_std_errors = true);

	virtual void CalculateElementMatrices(void);
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	void DDCAssembleGlobalMatrix();
#endif
	virtual void AssembleSystemMatrixNew(void);

	// MW
	int getFirstNodeBelowGWL(size_t i);
	bool hasConstrainedBC() { return _hasConstrainedBC; }
	bool hasConstrainedST() { return _hasConstrainedST; }
	void hasConstrainedBC(const bool state) { _hasConstrainedBC = state; }
	void hasConstrainedST(const bool state) { _hasConstrainedST = state; }

	void setidxVx(int index) { _idxVx = index; }
	void setidxVy(int index) { _idxVy = index; }
	void setidxVz(int index) { _idxVz = index; }

	// ST
	void IncorporateSourceTerms(const int rank = -1);
// WW void CheckSTGroup(); //OK
#ifdef GEM_REACT
	void IncorporateSourceTerms_GEMS(void); // HS: dC/dt from GEMS chemical solver.
	int GetRestartFlag() const { return reload; }
#endif
	// BC
	void IncorporateBoundaryConditions(const int rank = -1);
	// PCH for FLUID_MOMENTUM
	void IncorporateBoundaryConditions(const int rank, const int axis);
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	void SetBoundaryConditionSubDomain(); // WW
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
// WW void CheckBCGroup(); //OK
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03~04.3012. WW
#ifdef NEW_EQS // WW
	void EQSInitialize();
	void EQSSolver(double* x); // PCH
#else
	int ExecuteLinearSolver(void);
	int ExecuteLinearSolver(LINEAR_SOLVER* eqs);
	LINEAR_SOLVER* getEQSPointer() const { return eqs; } // WW
#endif
#endif

	CTimeDiscretization* GetTimeStepping() const { return Tim; }
	double timebuffer; // YD
	//
	// NLS and CPL error and looping control
	int num_diverged;
	int num_notsatisfied;
	int iter_nlin;
	int iter_nlin_max;
	int iter_lin;
	int iter_lin_max;
	int iter_outer_cpl; // JT2012
	int iter_inner_cpl; // JT2012
	int pcs_num_dof_errors; // JT2012
	double pcs_relative_error[DOF_NUMBER_MAX]; // JT2012: for NLS, we store relative error for each DOF
	double pcs_absolute_error[DOF_NUMBER_MAX]; // JT2012: for NLS, we store error for each DOF
	double pcs_unknowns_norm;
	double cpl_max_relative_error; // JT2012: For CPL, we just store the maximum, not each dof value
	double cpl_max_relative_error_overall;
	double nls_max_relative_error;
	double cpl_absolute_error[DOF_NUMBER_MAX]; // JT2012:
	double temporary_absolute_error[DOF_NUMBER_MAX]; // JT2012:
	int temporary_num_dof_errors;
	int cpl_num_dof_errors; // JT2012
	bool first_coupling_iteration; // JT2012
	//
	// Specials
	void PCSMoveNOD();
	void PCSDumpModelNodeValues(void);
	// WW
	int GetNODValueIndex(const std::string& name, int timelevel);
	// BC for dynamic problems. WW
	inline void setBC_danymic_problems();
	inline void setST_danymic_problems();
	inline void setIC_danymic_problems();
	// Extropolate Gauss point values to node values. WW
	void Extropolation_GaussValue();
	void Extropolation_MatValue(); // WW
	// WW. 05.2009
	void Integration(std::vector<double>& node_velue);
	// Auto time step size control. WW
	void PI_TimeStepSize(); // WW
	bool TimeStepAccept() const { return accepted; }
	void SetDefaultTimeStepAccepted() { accepted = true; }
	// USER
	// ToDo
	double* TempArry; // MX
	void PCSOutputNODValues(void);
	void PCSSetTempArry(void); // MX
	double GetTempArryVal(int index) const // MX
	{
		return TempArry[index];
	}
	void LOPCopySwellingStrain(CRFProcess* m_pcs);
	VoidFuncInt PCSSetIC_USER;
	void SetIC();
	// Remove argument. WW
	void CalcSecondaryVariables(bool initial = false);
	void MMPCalcSecondaryVariablesRichards(int timelevel, bool update);
	// WW Reomve int timelevel, bool update
	// WW

	void CalcSecondaryVariablesTNEQ(); // HS
	void CalcSecondaryVariablesTES(); // HS
	void CalcSecondaryVariablesUnsaturatedFlow(bool initial = false);
	void CalcSecondaryVariablesPSGLOBAL(); // PCH
	void CalcSecondaryVariablesLiquidFlow(); // PCH
	double GetCapillaryPressureOnNodeByNeighobringElementPatches(int nodeIdx, int meanOption, double Sw);
	// JOD
	void CalcSaturationRichards(int timelevel, bool update);
	bool non_linear; // OK/CMCD
	// MX
	void InterpolateTempGP(CRFProcess*, std::string);
	// MX
	void ExtropolateTempGP(CRFProcess*, std::string);
	// Repeat Calculation,    JOD removed // HS reactivated
	void PrimaryVariableReload(); // YD
	void PrimaryVariableReloadRichards(); // YD
	void PrimaryVariableStorageRichards(); // YD
	bool adaption;
	void PrimaryVariableReloadTransport(); // kg44
	void PrimaryVariableStorageTransport(); // kg44
	// double GetNewTimeStepSizeTransport(double mchange); //kg44
	// FLX
	void CalcELEFluxes(const GEOLIB::Polyline* const ply, double* result);
	/**
	 * Necessary for the output of mass fluxes over polylines, BG 08/2011
	 * @param ply a pointer to a GEOLIB::Polyline
	 * @param NameofPolyline the name of the polyline
	 * @param result
	 */
	void CalcELEMassFluxes(const GEOLIB::Polyline* const ply, std::string const& NameofPolyline, double* result);
	double TotalMass[10]; // Necessary for the output of mass fluxes over polylines, BG 08/2011
	std::vector<std::string> PolylinesforOutput; // Necessary for the output of mass fluxes over polylines, BG 08/2011

	/**
	 * Necessary for the output of mass fluxes over polylines, BG 08/2011
	 * @param ele
	 * @param ElementConcentration
	 * @param grad
	 */
	void Calc2DElementGradient(MeshLib::CElem* ele, double ElementConcentration[4], double* grad);
	// NEW
	CRFProcess* CopyPCStoDM_PCS();
	bool OBJRelations(); // OK
	void OBJRelationsDelete(); // OK
	bool NODRelations(); // OK
	bool ELERelations(); // OK
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	bool CreateEQS(); // OK
	void EQSDelete(); // OK
	// Dumping matrix and RHS. WW
	void DumpEqs(std::string file_name);
#endif
	bool Check(); // OK
	void NODRelationsDelete(); // OK
	void ELERelationsDelete(); // OK
	bool m_bCheckOBJ; // OK
	bool m_bCheckNOD; // OK
	bool m_bCheckELE; // OK
	bool m_bCheckEQS; // OK
	void Delete(); // OK
	bool m_bCheck; // OK
	int ExcavMaterialGroup; // WX
	int ExcavDirection; // WX
	int ExcavCurve; // WX
	double ExcavBeginCoordinate; // WX
	int PCS_ExcavState; // WX
	int Neglect_H_ini; // WX
	int UpdateIniState; // WX
#if defined(USE_MPI) || defined(USE_PETSC) // WW
	void Print_CPU_time_byAssembly(std::ostream& os = std::cout) const
	{
		os << "\n***\nCPU time elapsed in the linear equation of " << convertProcessTypeToString(getProcessType())
		   << "\n";
		os << "--Global assembly: " << (double)cpu_time_assembly / CLOCKS_PER_SEC << "\n";
	}
#endif
	// HS 10.2011
	double m_rho_s_0;
	conversion_rate* m_conversion_rate;

#if defined(USE_PETSC) // 03.3012. WW
	/// Initialize the RHS array of the system of equations with the previous solution.
	void InitializeRHS_with_u0(const bool quad = false); // in rf_pcs1.cpp

	/// Initialize the unknows of equations with existing solutions within a index range.
	void initializeRHS_with_u0(const int min_id, const int max_id); // in rf_pcs1.cpp
#endif

private:
	/**
	 * Method configures the material parameters. For this purpose it searchs in all
	 * COutput objects (in the vector _nod_value_vector) for the values
	 * PERMEABILITY_X1 and POROSITY
	 */
	void configMaterialParameters();
	// method to check on constrained source terms
	bool checkConstrainedST(std::vector<CSourceTerm*>& st_vector, CSourceTerm const& st, CNodeValue const& st_node);
	// method to check on constrained boundary conditions
	bool checkConstrainedBC(CBoundaryCondition const& bc, CBoundaryConditionNode& bc_node, double& bc_value);
	std::valarray<double> getNodeVelocityVector(const long node_id);
	double calcPressureFromHead(CBoundaryCondition const& bc, std::size_t node_number, double pressure_value,
	                            double const& time_fac = 1, double const& fac = 1);
	double calcHeadFromPressure(CBoundaryCondition const& bc, std::size_t node_number, double pressure_value,
	                            double const& time_fac = 1, double const& fac = 1);

	double evaluteSwitchBC(CBoundaryCondition const& bc, CBoundaryConditionNode const& bc_node, double time_fac,
	                       double fac);

	int _pcs_constant_model;
	double _pcs_constant_value;
	int _pcs_constant_curve;
};

//========================================================================
// PCS
extern std::vector<CRFProcess*> pcs_vector;
extern std::vector<ElementValue*> ele_gp_value; // Gauss point value for velocity. WW
extern bool PCSRead(std::string);
extern void PCSWrite(std::string);
extern void RelocateDeformationProcess(CRFProcess* m_pcs);
extern void PCSDestroyAllProcesses(void);

extern CRFProcess* PCSGet(const std::string&);
/**
 * Function searchs in the global pcs_vector for a process with the process type pcs_type.
 * @param pcs_type process type
 * @return a pointer the the appropriate process or NULL (if not found)
 */
CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type); // TF

extern CRFProcess* PCSGetNew(const std::string&, const std::string&);
extern void PCSDelete();
extern void PCSDelete(const std::string&);
extern void PCSCreate();
// SB
extern int PCSGetPCSIndex(const std::string&, const std::string&);
// SB
extern CRFProcess* PCSGet(const std::string&, const std::string&);

/**
 * Function searchs in the global pcs_vector for a process
 * with the process type pcs_type and the primary function name
 * pv_name
 * @param pcs_type the process type
 * @param pv_name the name of the primary function
 * @return
 */
// TF
CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type, const std::string& pv_name);

// OK
extern CRFProcess* PCSGet(const std::string& variable_name, bool dummy);
extern CRFProcess* PCSGetFluxProcess(); // CMCD
extern CRFProcess* PCSGetFlow(); // OK//JT
extern CRFProcess* PCSGetHeat(); // JT
extern CRFProcess* PCSGetMass(size_t component_number); // JT
extern CRFProcess* PCSGetDeformation(); // JT
extern bool PCSConfig(); //
// NOD
extern int PCSGetNODValueIndex(const std::string&, int);
extern double PCSGetNODValue(long, char*, int);
extern void PCSSetNODValue(long, const std::string&, double, int);
// ELE
extern int PCSGetELEValueIndex(char*);
extern double PCSGetELEValue(long index, double* gp, double theta, const std::string& nod_fct_name);
// Specials
extern void PCSRestart();
extern std::string PCSProblemType();
// PCS global variables
extern int pcs_no_components;
extern bool pcs_monolithic_flow;
extern int pcs_deformation;

// ToDo
// SB
extern double PCSGetNODConcentration(long index, long component, long timelevel);
// SB
extern void PCSSetNODConcentration(long index, long component, long timelevel, double value);
extern char* GetCompNamehelp(char* name); // SB:namepatch - superseded by GetPFNamebyCPName
// SB4218
extern double PCSGetEleMeanNodeSecondary(long index,
                                         const std::string& pcs_name,
                                         const std::string& var_name,
                                         int timelevel);
// CB
extern double PCSGetEleMeanNodeSecondary_2(long index, int pcsT, const std::string& var_name, int timelevel);
extern std::string GetPFNamebyCPName(std::string line_string);

extern int memory_opt;

typedef struct /* Knotenwert-Informationen */
{
	char* name; /* Name der Knotengroesse */
	char* einheit; /* Einheit der Knotengroesse */
	int speichern; /* s.u., wird Wert gespeichert ? */
	int laden; /* s.u., wird Wert zu Beginn geladen ? */
	int restart; /* s.u., wird Wert bei Restart geladen ? */
	int adapt_interpol; /* Soll Wert bei Adaption auf Kinder interpoliert werden? */
	double vorgabe; /* Default-Wert fuer Initialisierung */
} NvalInfo;
extern int anz_nval; /* Anzahl der Knotenwerte */
extern int anz_nval0; /* Anzahl der Knotenwerte */
extern NvalInfo* nval_data;
extern int ModelsAddNodeValInfoStructure(char*, char*, int, int, int, int, double);

typedef struct /* Elementwert-Informationen */
{
	char* name; /* Name der Elementgroesse */
	char* einheit; /* Einheit der Elementgroesse */
	int speichern; /* s.u., wird Wert gespeichert ? */
	int laden; /* s.u., wird Wert zu Beginn geladen ? */
	int restart; /* s.u., wird Wert bei Restart geladen ? */
	int adapt_interpol; /* Soll Wert bei Adaption auf Kinder vererbt werden? */
	double vorgabe; /* Default-Wert fuer Initialisierung */
} EvalInfo;
extern int anz_eval; /* Anzahl der Elementwerte */
extern EvalInfo* eval_data;
extern int ModelsAddElementValInfoStructure(char*, char*, int, int, int, int, double);

extern int GetRFControlGridAdapt(void);
extern int GetRFControlModel(void);
extern int GetRFProcessChemicalModel(void);
extern int GetRFProcessFlowModel(void);
extern int GetRFProcessHeatReactModel(void);
extern int GetRFProcessNumPhases(void);
extern int GetRFProcessProcessing(char*);
extern int GetRFProcessProcessingAndActivation(const char*);
extern long GetRFProcessNumComponents(void);
extern int GetRFControlModex(void);
extern int GetRFProcessDensityFlow(void);
extern int GetRFProcessNumContinua(void);
extern int GetRFProcessNumElectricFields(void);
extern int GetRFProcessNumTemperatures(void);
extern int GetRFProcessSimulation(void);

extern void initializeConstrainedProcesses(std::vector<CRFProcess*>& pcs_vector);

// Coupling Flag. WW
extern bool T_Process; // Heat
extern bool H_Process; // Fluid
extern bool H2_Process; // Multi-phase
extern bool H3_Process; // 3-phase
extern bool M_Process; // Mechanical
extern bool RD_Process; // Richards
extern bool MH_Process; // MH monolithic scheme
extern bool MASS_TRANSPORT_Process; // Mass transport
extern bool FLUID_MOMENTUM_Process; // Momentum
extern bool RANDOM_WALK_Process; // RWPT
extern bool MULTI_COMPONENTIAL_FLOW_Process; // MCF
//
extern int pcs_number_deformation; // JT2012
extern int pcs_number_flow; // JT2012
extern int pcs_number_heat; // JT2012
extern std::vector<int> pcs_number_mass; // JT2012 (allow DOF_NUMBER_MAX components)
//
extern std::string project_title; // OK41
extern bool pcs_created;
extern std::vector<LINEAR_SOLVER*> PCS_Solver; // WW
// OK
extern void MMPCalcSecondaryVariablesNew(CRFProcess* m_pcs, bool NAPLdiss);
extern void SetFlowProcessType(); // CB 01/08
extern void CopyTimestepNODValuesSVTPhF(); // CB 13/08
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
//#ifndef NEW_EQS                                   //WW. 07.11.2008
extern void PCSCreateNew(); // OK
#endif
extern bool PCSCheck(); // OK
// New solvers WW
// Create sparse graph for each mesh    //1.11.2007 WW
#if defined(NEW_EQS) || defined(USE_PETSC) // || defined(other solver libs)//03.3012. WW
extern void CreateEQS_LinearSolver();
#endif

#ifdef GEM_REACT
class REACT_GEM;
extern REACT_GEM* m_vec_GEM;
#endif

#ifdef BRNS
class REACT_BRNS;
extern REACT_BRNS* m_vec_BRNS;
#endif

extern bool hasAnyProcessDeactivatedSubdomains; // NW
#endif
